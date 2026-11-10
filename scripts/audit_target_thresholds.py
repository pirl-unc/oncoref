#!/usr/bin/env python3
"""Compare explicit cancer-expression rules without changing normal exclusions."""

from __future__ import annotations

import itertools
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_target_sample_ranks import CACHE, REBUILD, sha  # noqa: E402

from oncoref import reference_data, target_selection  # noqa: E402

OUT = REBUILD.parent / "thresholds"
CAPS = (None, 4, 8, 16, 32, 128)


def patient_values(values, samples):
    """One maximum per documented patient; unknown IDs never form a fake donor."""
    known = samples.participant_id.ne("")
    return values.loc[known].groupby(samples.loc[known, "participant_id"]).max()


def second_largest(values):
    a = np.where(np.isfinite(values), values, -np.inf)
    return np.partition(a, -2, axis=0)[-2] if len(a) >= 2 else np.full(a.shape[1], -np.inf)


def cohort_quantile(values, samples, q, minimum_patients=20):
    """Observed lower order statistic, with per-protein donor support."""
    summaries = []
    for _, meta in samples[samples.participant_id.ne("")].groupby("cancer_code"):
        donors = patient_values(values.loc[meta.index], meta)
        if len(donors) < minimum_patients:
            continue
        quantiles = donors.quantile(q, interpolation="lower")
        summaries.append(quantiles.where(donors.count().ge(minimum_patients)))
    if not summaries:
        return pd.Series(np.nan, index=values.columns)
    return pd.concat(summaries, axis=1).max(axis=1)


def passes(value, normal, fold, floor):
    """One enrichment denominator; equality to normal never qualifies."""
    return (value >= fold * np.maximum(normal, floor)) & (value > normal)


def disposable_ntpm(candidates):
    """Independent HPA nTPM cap observations; never mix with the TPM maximum."""
    source = Path(reference_data.local_path("hpa_rna_consensus", "v23"))
    identity_path = (
        ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz"
    )
    identities = pd.read_csv(identity_path).set_index("source_gene_id").canonical_gene_id
    lookup = {
        g: entity for entity, row in candidates.iterrows() for g in json.loads(row.member_gene_ids)
    }
    raw = pd.read_csv(source, sep="\t")
    raw = raw[raw.Tissue.isin(["appendix", "gallbladder", "tonsil"])].copy()
    raw["entity_id"] = raw.Gene.map(identities).map(lookup)
    raw = raw[raw.entity_id.notna()]
    rows = raw.groupby(["entity_id", "Tissue"]).nTPM.agg(["sum", "count"]).reset_index()
    expected = candidates.member_gene_ids.map(lambda v: len(json.loads(v)))
    rows["value"] = rows["sum"].where(rows["count"].eq(rows.entity_id.map(expected)))
    rows["unit"] = "nTPM"
    rows["source"] = "hpa_rna_consensus:v23"
    rows.to_csv(OUT / "disposable-bulk-ntpm.csv", index=False)
    maximum = rows.groupby("entity_id").value.max().reindex(candidates.index)
    assert maximum.notna().all()
    return maximum, {
        "source": "hpa_rna_consensus:v23",
        "path": str(source),
        "sha256": sha(source),
        "identity_path": str(identity_path.relative_to(ROOT)),
        "identity_sha256": sha(identity_path),
    }


def rule_specs():
    rules = []
    for floor, repeated_fold, single_fold in itertools.product(
        (1, 2, 5, 10), (1, 2, 3, 5), (5, 10, 20)
    ):
        rules.append(
            {
                "id": f"fold_L{floor}_R{repeated_fold}_S{single_fold}",
                "kind": "fold",
                "floor": floor,
                "repeated_fold": repeated_fold,
                "single_fold": single_fold,
                "repeated_rank": 0,
                "single_rank": 0,
            }
        )
    for repeated_fold in (1, 2):
        for repeated_rank, single_rank in ((0.8, 0.9), (0.8, 0.95), (0.9, 0.95), (0.95, 0.95)):
            rules.append(
                {
                    "id": f"rank_P{round(repeated_rank * 100)}_P{round(single_rank * 100)}_R{repeated_fold}",
                    "kind": "rank",
                    "floor": 1,
                    "repeated_fold": repeated_fold,
                    "single_fold": 5,
                    "repeated_rank": repeated_rank,
                    "single_rank": single_rank,
                }
            )
    for q, floor, fold, rare in itertools.product((0.8, 0.9, 0.95), (1, 5), (1, 2), (0, 5)):
        rules.append(
            {
                "id": f"cohort_Q{round(q * 100)}_L{floor}_R{fold}_S{rare}",
                "kind": "cohort",
                "floor": floor,
                "repeated_fold": fold,
                "single_fold": rare,
                "repeated_rank": 0,
                "single_rank": 0,
                "quantile": q,
            }
        )
    for fold in (1, 2, 5):
        rules.append(
            {
                "id": f"repeated_only_R{fold}",
                "kind": "fold",
                "floor": 1,
                "repeated_fold": fold,
                "single_fold": 0,
                "repeated_rank": 0,
                "single_rank": 0,
            }
        )
    for fold in (5, 10, 20):
        rules.append(
            {
                "id": f"single_only_S{fold}",
                "kind": "fold",
                "floor": 1,
                "repeated_fold": 0,
                "single_fold": fold,
                "repeated_rank": 0,
                "single_rank": 0,
            }
        )
    return rules


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    receipt = json.loads((CACHE / "receipt.json").read_text())
    for path, digest in receipt["implementation_sha256"].items():
        assert sha(ROOT / path) == digest
    for path, digest in receipt["outputs_sha256"].items():
        assert sha(CACHE / path) == digest
    matrix_path = ROOT / "tmp/target-catalog-cache/tumor-comparison.parquet"
    assert (
        sha(matrix_path)
        == receipt["candidate_matrix_sha256"]
        == target_selection.metadata()["comparison"]["tumor_matrix_cache_sha256"]
    )
    assert sha(ROOT / "oncoref/data/target-selection.csv.gz") == receipt["catalog_sha256"]
    assert sha(REBUILD / "tumor-samples.csv.gz") == receipt["samples_sha256"]
    catalog = target_selection.results().set_index("entity_id")
    x = catalog[catalog[list(target_selection.SHARED_CHECKS)].all(axis=1)].copy()
    m = pd.read_parquet(matrix_path).reindex(columns=x.index)
    ranks = pd.read_parquet(CACHE / "sample-ranks.parquet").reindex(index=m.index, columns=x.index)
    samples = (
        pd.read_csv(REBUILD / "tumor-samples.csv.gz", keep_default_na=False)
        .set_index("sample_key")
        .loc[m.index]
    )
    normal = x.normal_reference_tpm.to_numpy()
    normal_rows = pd.read_csv(REBUILD / "normal-comparison.csv.gz")
    assert set(normal_rows.unit) == {"TPM"}
    normal_rows = normal_rows[
        normal_rows.qc_pass & normal_rows.identity_resolved & normal_rows.value.notna()
    ]
    disposable = (
        normal_rows[normal_rows.scope.eq("disposable")]
        .groupby("entity_id")
        .value.max()
        .reindex(x.index)
    )
    adult_r = (
        normal_rows[normal_rows.scope.eq("reproductive")]
        .groupby("entity_id")
        .value.max()
        .reindex(x.index)
    )
    assert disposable.notna().all()
    disposable_n, ntpm_receipt = disposable_ntpm(x)
    peak = pd.concat([x.tumor_max_tpm, x.developmental_mean_tpm, adult_r], axis=1).max(axis=1)
    maximum_pass = peak.gt(np.maximum(x.normal_reference_tpm, disposable)) & peak.ge(1)
    x["disposable_max_tpm"] = disposable
    x["disposable_max_ntpm"] = disposable_n
    x["target_maximum_pass"] = maximum_pass
    x["tumor_max_within_sample_rank"] = ranks.max()
    top = np.nanargmax(m.to_numpy(), axis=0)
    x["rank_at_tumor_peak"] = ranks.to_numpy()[top, np.arange(len(x))]
    variants = {
        "full": m,
        "without_hl_merkel": m.loc[~samples.cancer_code.isin(["HL", "NEC_MERKEL"])],
    }
    without_peak = m.copy()
    a = without_peak.to_numpy(copy=True)
    a[top, np.arange(len(x))] = np.nan
    variants["without_largest_sample"] = pd.DataFrame(a, index=m.index, columns=m.columns)
    features = {}
    for variant, values in variants.items():
        meta = samples.loc[values.index]
        for rank in (0, 0.8, 0.9, 0.95):
            masked = values if rank == 0 else values.where(ranks.loc[values.index].ge(rank))
            features[(variant, rank, "second_patient")] = second_largest(
                patient_values(masked, meta).to_numpy()
            )
            features[(variant, rank, "maximum")] = masked.max().fillna(-np.inf).to_numpy()
        for q in (0.8, 0.9, 0.95):
            features[(variant, q, "cohort")] = (
                cohort_quantile(values, meta, q).fillna(-np.inf).to_numpy()
            )
    for q in (0.8, 0.9, 0.95):
        x[f"highest_cohort_p{round(q * 100)}_tpm"] = features[("full", q, "cohort")]
    rules = rule_specs()
    rows, decisions = [], {}
    nom = x.published_nomination.to_numpy()
    rbase = x.is_reproductive_antigen.to_numpy()
    old = x.is_target_candidate.to_numpy()
    families = (
        x.protein_group_label.str.match(r"^OR\d").to_numpy(),
        x.protein_group_label.str.startswith("KRTAP").to_numpy(),
    )
    for rule in rules:
        qualification = {}
        for variant in variants:
            repeated_value = (
                features[(variant, rule["quantile"], "cohort")]
                if rule["kind"] == "cohort"
                else features[(variant, rule["repeated_rank"], "second_patient")]
            )
            repeated = (
                passes(repeated_value, normal, rule["repeated_fold"], rule["floor"])
                if rule["repeated_fold"]
                else np.zeros(len(x), bool)
            )
            single = (
                passes(
                    features[(variant, rule["single_rank"], "maximum")],
                    normal,
                    rule["single_fold"],
                    rule["floor"],
                )
                if rule["single_fold"]
                else np.zeros(len(x), bool)
            )
            qualification[variant] = repeated | single
        decisions[rule["id"]] = qualification["full"]
        for cap in CAPS:
            eligible = maximum_pass.to_numpy() & (
                (disposable.le(cap) & disposable_n.le(cap)).to_numpy() if cap is not None else True
            )
            cancer = qualification["full"] & eligible
            reproductive = rbase & eligible
            cta = cancer & reproductive
            cancer_only = cancer & ~reproductive
            union = cancer | reproductive
            rows.append(
                {
                    **rule,
                    "disposable_cap_tpm": cap,
                    "disposable_cap_ntpm": cap,
                    "Cancer Antigens": int(cancer.sum()),
                    "CTA": int(cta.sum()),
                    "Cancer Antigen only": int(cancer_only.sum()),
                    "Reproductive Antigen only": int((reproductive & ~cancer).sum()),
                    "Target candidates": int(union.sum()),
                    "Other": int(
                        catalog.category.eq("Other").sum()
                        + (old & ~union).sum()
                        - (~old & union).sum()
                    ),
                    "Unsafe": int(catalog.category.eq("Unsafe").sum()),
                    "current_targets_to_other": int((old & ~union).sum()),
                    "cancer_membership_lost": int((x.is_cancer_antigen & ~cancer).sum()),
                    "published_cancer_nominations": int((cancer & nom).sum()),
                    "published_cta_nominations": int((cta & nom).sum()),
                    "olfactory_cancer_only": int((cancer_only & families[0]).sum()),
                    "krtap_cancer_only": int((cancer_only & families[1]).sum()),
                    "cancer_lost_without_hl_merkel": int(
                        (cancer & ~qualification["without_hl_merkel"]).sum()
                    ),
                    "cancer_lost_without_largest_sample": int(
                        (cancer & ~qualification["without_largest_sample"]).sum()
                    ),
                }
            )
    grid = pd.DataFrame(rows)
    grid.to_csv(OUT / "threshold-grid.csv", index=False)
    decisions = pd.DataFrame(decisions, index=x.index)
    decisions.to_csv(
        OUT / "cancer-rule-decisions.csv.gz", compression={"method": "gzip", "mtime": 0}
    )
    x.replace([-np.inf, np.inf], np.nan).to_csv(
        OUT / "protein-threshold-evidence.csv.gz", compression={"method": "gzip", "mtime": 0}
    )
    cohort_rows = []
    for code, meta in samples.groupby("cancer_code"):
        known = meta[meta.participant_id.ne("")]
        cohort_rows.append(
            {
                "cancer_code": code,
                "physical_samples": len(meta),
                "known_patients": known.participant_id.nunique(),
                "unknown_patient_samples": meta.participant_id.eq("").sum(),
            }
        )
    pd.DataFrame(cohort_rows).to_csv(OUT / "patient-coverage.csv", index=False)
    shutil.copyfile(CACHE / "sample-rank-cutoffs.csv.gz", OUT / "sample-rank-cutoffs.csv.gz")
    shutil.copyfile(CACHE / "receipt.json", OUT / "rank-build-receipt.json")
    baseline = grid[grid.id.eq("fold_L1_R1_S5") & grid.disposable_cap_tpm.isna()].iloc[0]
    assert baseline["Cancer Antigens"] == 882 and baseline.CTA == 632
    selected = [
        "fold_L1_R1_S5",
        "fold_L1_R2_S5",
        "fold_L1_R2_S10",
        "fold_L2_R2_S5",
        "fold_L5_R1_S5",
        "rank_P80_P95_R1",
        "rank_P80_P95_R2",
        "rank_P90_P95_R1",
        "rank_P95_P95_R1",
        "cohort_Q80_L1_R1_S0",
        "cohort_Q90_L1_R1_S0",
        "cohort_Q95_L1_R1_S0",
        "cohort_Q95_L1_R1_S5",
    ]
    grid[grid.id.isin(selected) & grid.disposable_cap_tpm.isna()].to_csv(
        OUT / "selected-rules.csv", index=False
    )
    summary = {
        "eligible_proteins": len(x),
        "rules": len(rules),
        "scenarios": len(grid),
        "physical_samples": len(samples),
        "known_patients": samples.loc[samples.participant_id.ne(""), "participant_id"].nunique(),
        "unknown_patient_samples": int(samples.participant_id.eq("").sum()),
        "cohorts": samples.cancer_code.nunique(),
        "cohorts_with_20_known_patients": sum(r["known_patients"] >= 20 for r in cohort_rows),
        "disposable_ntpm_source": ntpm_receipt,
        "input_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [
                matrix_path,
                REBUILD / "normal-comparison.csv.gz",
                ROOT / "oncoref/data/target-selection.csv.gz",
                CACHE / "sample-ranks.parquet",
            ]
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        grid[grid.id.isin(selected) & grid.disposable_cap_tpm.isna()][
            [
                "id",
                "Cancer Antigens",
                "CTA",
                "Cancer Antigen only",
                "Reproductive Antigen only",
                "Target candidates",
                "current_targets_to_other",
                "cancer_lost_without_hl_merkel",
                "cancer_lost_without_largest_sample",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    build()
