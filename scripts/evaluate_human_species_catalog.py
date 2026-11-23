#!/usr/bin/env python3
"""Materialize #737 linear-TPM tumor observations and verify unified human catalog evaluation (#755)."""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import (  # noqa: E402
    human,
    human_screening_inputs,
    library_qc,
    samples,
    source_matrices,
    target_selection,
)
from oncoref.gene_matching import (  # noqa: E402
    harmonize_reference_expression,
    match_reference_genes,
)

COHORTS_737 = ("THYMCA", "BCC", "cSCC")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build_737_tumor_observations(out_dir: Path) -> dict:
    """Extract QC-passing, identity-harmonized per-sample linear-TPM rows for #737 cohorts."""
    ts = target_selection.results()
    candidates = ts[~ts["category"].eq("Unsafe")].copy()
    gene_to_entity = {
        g: r.entity_id
        for r in candidates.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }
    expected_members = {
        r.entity_id: len(json.loads(r.member_gene_ids)) for r in candidates.itertuples(index=False)
    }

    obs_rows: list[dict] = []
    cohort_receipts: list[dict] = []

    for code in COHORTS_737:
        ver = source_matrices.source_matrix_version(code)
        mpath = source_matrices.local_path(code, ver)
        if not mpath.is_file():
            raise FileNotFoundError(
                f"Missing local source matrix {mpath}; run `oncoref data fetch per-sample:{code}`"
            )
        msha = _sha256_file(mpath)
        s_manifest = samples.samples_for_cancer_code(code, included_only=True)
        source_cohort = str(s_manifest["source_cohort"].iloc[0])
        sample_ids = s_manifest["sample_id"].astype(str).tolist()

        raw = pq.read_table(mpath).to_pandas()
        signatures = library_qc.signatures(raw, sample_ids)
        signatures["reasons"] = signatures.apply(lambda r: ";".join(library_qc.reasons(r)), axis=1)
        qc_passing = signatures.index[signatures["reasons"].eq("")].tolist()

        id_rows = raw[["Ensembl_Gene_ID", "Symbol"]].copy()
        identity = match_reference_genes(
            id_rows,
            id_column="Ensembl_Gene_ID",
            symbol_column="Symbol",
            source=source_cohort,
            source_version=ver,
        )
        matched = identity[
            identity["identity_status"].eq("matched")
            & identity["canonical_gene_id"].isin(gene_to_entity)
        ]
        sub_raw = raw[raw["Ensembl_Gene_ID"].isin(matched["source_gene_id"])].copy()
        harmonized, _decisions = harmonize_reference_expression(
            sub_raw,
            sample_cols=qc_passing,
            source=source_cohort,
            source_version=ver,
            symbol_column="Symbol",
            duplicate_policy="quarantine",
        )
        harmonized["entity_id"] = harmonized["Ensembl_Gene_ID"].map(gene_to_entity)
        values = harmonized.groupby("entity_id")[qc_passing].sum(min_count=1)
        counts = harmonized.groupby("entity_id")[qc_passing].count()
        expected_series = pd.Series(expected_members).reindex(values.index)
        values = values.where(counts.eq(expected_series, axis=0))

        ident_df = samples.sample_identities(
            pd.DataFrame({"source_cohort": source_cohort, "sample_id": qc_passing})
        ).set_index("sample_id")

        kept_obs = 0
        for eid, s_vals in values.iterrows():
            valid_vals = s_vals.dropna()
            # Retain every sample >= 1.0 TPM (detection floor; 10x H0 requires >= 10 TPM)
            above_floor = valid_vals[valid_vals >= 1.0]
            for sid, val in above_floor.items():
                pid = str(ident_df.at[sid, "participant_id"] or "")
                obs_rows.append(
                    {
                        "entity_id": str(eid),
                        "cancer_code": code,
                        "source_cohort": source_cohort,
                        "sample_id": str(sid),
                        "participant_id": pid,
                        "value": round(float(val), 6),
                    }
                )
                kept_obs += 1

        cohort_receipts.append(
            {
                "cancer_code": code,
                "source_cohort": source_cohort,
                "source_matrix_version": ver,
                "matrix_sha256": msha,
                "included_samples": len(sample_ids),
                "qc_passing_samples": len(qc_passing),
                "retained_observations_ge_1_tpm": kept_obs,
            }
        )

    obs_df = (
        pd.DataFrame(obs_rows)
        .sort_values(["cancer_code", "entity_id", "sample_id"])
        .reset_index(drop=True)
    )
    csv_gz_path = out_dir / "human-737-tumor-observations.csv.gz"
    meta_path = out_dir / "human-737-tumor-observations.json"

    csv_bytes = obs_df.to_csv(index=False).encode("utf-8")
    with gzip.GzipFile(filename="", mode="wb", fileobj=csv_gz_path.open("wb"), mtime=0) as gz:
        gz.write(csv_bytes)

    digest = _sha256_file(csv_gz_path)
    receipt = {
        "schema": "oncoref.human-737-tumor-observations.v1",
        "baseline_policy_id": target_selection.POLICY_ID,
        "cohorts": cohort_receipts,
        "row_count": len(obs_df),
        "bytes": csv_gz_path.stat().st_size,
        "sha256": digest,
    }
    meta_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def verify_baseline_parity_and_write_preview(out_dir: Path) -> dict:
    """Assert 100% baseline parity on evaluate_species_catalog and write integrated preview."""
    baseline = target_selection.results()
    evaluated_baseline = human_screening_inputs.evaluate_human_catalog(
        include_home_panels=False,
        include_adjudications=False,
        include_737_cohorts=False,
    )
    eval_cat = evaluated_baseline.catalog

    mismatches = baseline["classification"].ne(eval_cat["classification"])
    if mismatches.any():
        bad = pd.DataFrame(
            {
                "entity_id": baseline.loc[mismatches, "entity_id"],
                "label": baseline.loc[mismatches, "protein_group_label"],
                "expected": baseline.loc[mismatches, "classification"],
                "got": eval_cat.loc[mismatches, "classification"],
            }
        )
        raise RuntimeError(
            f"Baseline classification parity failed on {int(mismatches.sum())} rows:\n"
            + bad.head(20).to_string(index=False)
        )

    expected_sex = human.sex_restricted_antigens().set_index("entity_id")["sex_restriction_code"]
    got_sex = evaluated_baseline.sex_restricted_antigens().set_index("entity_id")[
        "sex_restriction_code"
    ]
    if not expected_sex.equals(got_sex):
        diff = expected_sex.compare(got_sex)
        raise RuntimeError(
            f"Baseline sex-restriction parity failed on {len(diff)} rows:\n{diff.head(20)}"
        )

    changed, summary = human_screening_inputs.compute_next_policy_transitions()
    preview_csv = out_dir / "human-integrated-policy-preview.csv"
    preview_json = out_dir / "human-integrated-policy-preview.json"
    changed.to_csv(preview_csv, index=False)
    summary["sha256"] = _sha256_file(preview_csv)
    summary["bytes"] = preview_csv.stat().st_size
    preview_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    out_dir = ROOT / "oncoref" / "data"
    obs_receipt = build_737_tumor_observations(out_dir)
    print(f"Wrote #737 observations: {obs_receipt['row_count']} rows ({obs_receipt['bytes']} B)")
    from oncoref import load_dataset

    load_dataset._clear_cache()
    summary = verify_baseline_parity_and_write_preview(out_dir)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
