#!/usr/bin/env python3
"""Rebuild the human CTA source-intake and normal-tissue-scope comparison."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
from math import fsum
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import cta, cta_provenance, cta_regen, cta_sources, hpa, reference_data  # noqa: E402
from oncoref.cta_curation_plots import PRIMARY_SOURCES  # noqa: E402
from oncoref.cta_tissues import PROTEIN_DETECTED_LEVELS, adaptive_rna_threshold  # noqa: E402
from oncoref.load_dataset import get_data  # noqa: E402
from oncoref.normal_tissues import normal_tissue_taxonomy, resolve_normal_tissue  # noqa: E402
from oncoref.version import __version__  # noqa: E402

SCOPES = ("reproductive", "reproductive_or_disposable")
LABELS = {"reproductive": "Reproductive", "reproductive_or_disposable": "Reproductive + disposable"}
STRICT_IHC_LEVELS = PROTEIN_DETECTED_LEVELS | {"Ascending", "Descending"}
SOURCE_NAMES = {"bulk_rna": "hpa_rna_consensus", "ihc": "hpa_normal_tissue"}
PAPER_LABELS = {
    **{doi: name for name, doi in PRIMARY_SOURCES.items()},
    "10.1038/s41467-020-19141-w": "Bradley 2020",
}


def source_identity(row):
    """Collapse mapped aliases; preserve the native namespace of unmapped rows."""
    if row["Ensembl_Gene_ID"]:
        return row["Ensembl_Gene_ID"]
    native = row["source_gene_id"] or row["source_symbol"]
    namespace = row["source_id_type"] or ("ensembl" if native.startswith("ENSG") else "symbol")
    return f"unmapped:{namespace}:{native}"


def annotate_measurements(raw, modality, candidate_ids):
    """Preserve assay observations and resolve distinct source contexts once."""
    rows = raw.loc[raw.Gene.isin(candidate_ids)].copy()
    if modality == "ihc":
        rows = rows.loc[~rows[["Tissue", "Cell type", "Level"]].isna().all(axis=1)]
    rows = rows.rename(
        columns={
            "Gene": "Ensembl_Gene_ID",
            "Gene name": "Symbol",
            "Tissue": "tissue",
            "Cell type": "cell_type",
            "nTPM": "rna_ntpm",
            "Level": "ihc_level",
            "Reliability": "ihc_reliability",
        }
    )
    for col in ("tissue", "cell_type", "ihc_level", "ihc_reliability"):
        if col not in rows:
            rows[col] = ""
        rows[col] = rows[col].fillna("")
    if "rna_ntpm" not in rows:
        rows["rna_ntpm"] = np.nan
    if modality == "bulk_rna":
        values = pd.to_numeric(rows.rna_ntpm, errors="raise")
        if values.lt(0).any() or np.isinf(values).any():
            raise ValueError("Invalid normal-tissue RNA value")
        if rows.duplicated(["Ensembl_Gene_ID", "tissue"]).any():
            raise ValueError("Duplicate consensus RNA measurement")
        rows["valid_measurement"] = values.notna()
        rows["strict_positive"] = values.gt(0)
        rows["adaptive_positive"] = values.ge(1)
    else:
        rows["valid_measurement"] = rows.ihc_level.isin(STRICT_IHC_LEVELS | {"Not detected"})
        rows["strict_positive"] = rows.ihc_level.isin(STRICT_IHC_LEVELS)
        rows["adaptive_positive"] = rows.ihc_level.isin(PROTEIN_DETECTED_LEVELS)
    contexts = []
    for tissue, cell_type in (
        rows[["tissue", "cell_type"]].drop_duplicates().itertuples(index=False)
    ):
        annotation = resolve_normal_tissue(
            tissue or None,
            species="human",
            cell_type=cell_type or None,
            source_id=SOURCE_NAMES[modality],
            source_version="v23",
        )
        contexts.append(
            {
                "tissue": tissue,
                "cell_type": cell_type,
                "concept_id": annotation.concept_id or "",
                "mapping_status": annotation.mapping_status,
                "mapping_kind": annotation.mapping_kind,
                "mapping_references": ";".join(annotation.references),
                **{scope: scope in annotation.scopes for scope in SCOPES},
            }
        )
    mapped = pd.DataFrame(
        contexts,
        columns=[
            "tissue",
            "cell_type",
            "concept_id",
            "mapping_status",
            "mapping_kind",
            "mapping_references",
            *SCOPES,
        ],
    )
    rows = rows.merge(mapped, on=["tissue", "cell_type"], how="left", validate="many_to_one")
    # An empty assay table must retain Boolean masks when concatenated with RNA.
    # Object-dtype Python bools invert to integers (-1/-2), not logical negatives.
    for column in (*SCOPES, "valid_measurement", "strict_positive", "adaptive_positive"):
        rows[column] = rows[column].astype(bool)
    rows["modality"] = modality
    rows["source_name"] = SOURCE_NAMES[modality]
    rows["source_version"] = "v23"
    rows["context_label"] = rows.tissue + np.where(
        rows.cell_type.ne(""), " / " + rows.cell_type, ""
    )
    return rows


def evidence_status(rows, scope):
    """Compatibility of positive observations, with missingness kept distinct."""
    if rows.empty:
        return "unavailable"
    if not rows.valid_measurement.any():
        return "unresolved"
    positive = rows.loc[rows.strict_positive]
    if (positive.mapping_status.eq("mapped") & ~positive[scope]).any():
        return "outside_scope"
    if positive.mapping_status.ne("mapped").any() or (~rows.valid_measurement).any():
        return "unresolved"
    return "compatible" if len(positive) else "no_positive_observation"


def routine_contexts(raw, modality):
    """Source contexts measured for at least half of its measured genes."""
    columns = ["Tissue"] if modality == "bulk_rna" else ["Tissue", "Cell type"]
    counts = raw.dropna(subset=columns).groupby(columns).Gene.nunique()
    selected = counts[counts.ge(raw.Gene.nunique() / 2)].index
    return {label if isinstance(label, tuple) else (label, "") for label in selected}


def compare_candidates(table, observations, *, unreliable_ids, expected_contexts):
    """Two scope comparisons: existing adaptive policy and raw positive-evidence sensitivity.

    Adaptive decisions retain the extended panel's family, expression, nomination
    and manual review rules. Strict compatibility requires complete routine bulk
    RNA, positive evidence in scope, and no outside positive RNA/IHC. Unusable or
    reviewed unreliable IHC blocks strict compatibility. Missing IHC remains a
    separately labeled RNA-only result. Neither establishes adequate absence.
    """
    groups = {
        modality: dict(iter(frame.groupby("Ensembl_Gene_ID")))
        for modality, frame in observations.items()
    }
    empty = {modality: frame.iloc[:0] for modality, frame in observations.items()}
    records = []
    for scope in SCOPES:
        adaptive = table.copy().reset_index(drop=True)
        scope_records = []
        for i, row in adaptive.iterrows():
            gid = row.Ensembl_Gene_ID
            rna = groups["bulk_rna"].get(gid, empty["bulk_rna"])
            ihc = groups["ihc"].get(gid, empty["ihc"])
            rna_status = evidence_status(rna, scope)
            ihc_status = evidence_status(ihc, scope)
            seen_rna = set(
                zip(
                    rna.loc[rna.valid_measurement, "tissue"],
                    rna.loc[rna.valid_measurement, "cell_type"],
                )
            )
            seen_ihc = set(
                zip(
                    ihc.loc[ihc.valid_measurement, "tissue"],
                    ihc.loc[ihc.valid_measurement, "cell_type"],
                )
            )
            rna_missing = expected_contexts["bulk_rna"] - seen_rna
            ihc_missing = expected_contexts["ihc"] - seen_ihc
            rna_complete = bool(expected_contexts["bulk_rna"]) and not rna_missing
            ihc_unreliable = gid in unreliable_ids
            pos_ihc = ihc.loc[ihc.adaptive_positive]
            protein_present = len(pos_ihc) > 0 and not ihc_unreliable
            protein_ok = (
                bool(pos_ihc[scope].all() and pos_ihc.mapping_status.eq("mapped").all())
                if protein_present
                else True
            )
            ntpm = dict(zip(rna.tissue.str.lower(), rna.rna_ntpm))
            allowed = frozenset(rna.loc[rna[scope], "tissue"].str.lower())
            fraction = (
                round(cta_regen._fraction(ntpm, allowed, deflate=True), 4) if ntpm else np.nan
            )
            deflated_total = (
                fsum(rna.rna_ntpm.dropna().sub(1).clip(lower=0))
                if rna.valid_measurement.any()
                else np.nan
            )
            adaptive.at[i, "rna_deflated_reproductive_frac"] = fraction
            adaptive.at[i, "protein_reproductive"] = (
                str(protein_ok) if protein_present else "no data"
            )
            # All other adaptive fields stay bound to the extended snapshot.
            adaptive.at[i, "passes_filters"] = cta_regen._passes_filters_rule(adaptive.loc[i], 0.97)
            positives = pd.concat(
                [rna.loc[rna.strict_positive], ihc.loc[ihc.strict_positive]], ignore_index=True
            )
            outside = positives.loc[positives.mapping_status.eq("mapped") & ~positives[scope]]
            strict_unknown = (
                not rna_complete
                or rna_status in {"unavailable", "unresolved"}
                or ihc_status == "unresolved"
                or ihc_unreliable
            )
            if len(outside):
                strict_status = "outside_scope"
            elif strict_unknown:
                strict_status = "unresolved"
            elif positives.empty:
                strict_status = "no_positive_observation"
            else:
                strict_status = "compatible"
            outside_rna = rna.loc[
                rna.strict_positive & rna.mapping_status.eq("mapped") & ~rna[scope]
            ]
            outside_ihc = ihc.loc[
                ihc.strict_positive & ihc.mapping_status.eq("mapped") & ~ihc[scope]
            ]
            gate_failures = []
            if not protein_ok:
                gate_failures.append("protein_outside_selected_scope")
            if pd.isna(fraction):
                gate_failures.append("missing_rna")
            elif fraction < adaptive_rna_threshold(row.protein_reliability):
                gate_failures.append("rna_scope_fraction_below_threshold")
            strict_failures = []
            if len(outside_rna):
                strict_failures.append("positive_rna_outside_scope")
            if len(outside_ihc):
                strict_failures.append("positive_raw_ihc_outside_scope")
            if not rna_complete:
                strict_failures.append("incomplete_routine_rna")
            if rna_status in {"unavailable", "unresolved"}:
                strict_failures.append(f"rna_{rna_status}")
            if ihc_status == "unresolved":
                strict_failures.append("ihc_unresolved")
            if ihc_unreliable:
                strict_failures.append("reviewed_unreliable_ihc")
            if positives.empty:
                strict_failures.append("no_positive_observation")
            scope_records.append(
                {
                    "Ensembl_Gene_ID": gid,
                    "Symbol": row.Symbol,
                    "scope": scope,
                    "family_eligible": gid not in cta.NON_CTA_EXCLUDED_GENE_IDS,
                    "rna_deflated_scope_fraction": fraction,
                    "rna_max_ntpm": rna.rna_ntpm.max(),
                    "rna_deflated_total_ntpm": deflated_total,
                    "rna_deflated_denominator_zero": bool(deflated_total == 0),
                    "adaptive_rna_threshold": adaptive_rna_threshold(row.protein_reliability),
                    "adaptive_protein_pass": protein_ok,
                    "adaptive_gate_failures": ";".join(gate_failures),
                    "protein_reliability": row.protein_reliability,
                    "reviewed_unreliable_ihc": ihc_unreliable,
                    "rna_status": rna_status,
                    "ihc_status": ihc_status,
                    "strict_status": strict_status,
                    "strict_evidence_failures": ";".join(strict_failures),
                    "rna_observed_contexts": len(seen_rna),
                    "rna_routine_missing": len(rna_missing),
                    "ihc_observed_contexts": len(seen_ihc),
                    "ihc_routine_missing": len(ihc_missing),
                    "ihc_coverage": "unreliable"
                    if ihc_unreliable
                    else "unavailable"
                    if not seen_ihc
                    else "routine_complete"
                    if not ihc_missing
                    else "partial",
                    "outside_rna_max_ntpm": outside_rna.rna_ntpm.max()
                    if len(outside_rna)
                    else 0.0
                    if rna.valid_measurement.any()
                    else np.nan,
                    "outside_rna_tissues": ";".join(sorted(set(outside_rna.tissue))),
                    "outside_ihc_contexts": ";".join(sorted(set(outside_ihc.context_label))),
                    "disposable_positive_contexts": ";".join(
                        sorted(
                            set(
                                positives.loc[
                                    positives.reproductive_or_disposable & ~positives.reproductive,
                                    "context_label",
                                ]
                            )
                        )
                    ),
                    "absence_evidence": "not_established",
                }
            )
        reviewed = cta._with_specificity_columns(adaptive, tissue_scope="extended")
        # Holdbacks independent of the HPA gate also apply to the strict screen.
        gate_neutral = adaptive.copy()
        gate_neutral["passes_filters"] = True
        gate_neutral["never_expressed"] = False
        nomination_review = cta._with_specificity_columns(gate_neutral, tissue_scope="extended")
        for record, row, review in zip(
            scope_records,
            reviewed.itertuples(index=False),
            nomination_review.itertuples(index=False),
        ):
            record.update(
                {
                    "adaptive_hpa_pass": bool(row.passes_filters),
                    "adaptive_after_family": record["family_eligible"] and bool(row.passes_filters),
                    "adaptive_retained": record["family_eligible"]
                    and row.specificity_action == "include_default",
                    "adaptive_review_action": row.specificity_action,
                    "adaptive_review_status": row.specificity_status,
                    "adaptive_review_reason": row.specificity_rationale,
                    "nomination_review_eligible": review.specificity_action == "include_default",
                    "nomination_review_status": review.specificity_status,
                    "nomination_review_reason": review.specificity_rationale,
                    "strict_after_family": record["family_eligible"]
                    and record["strict_status"] == "compatible",
                    "strict_retained": record["family_eligible"]
                    and review.specificity_action == "include_default"
                    and record["strict_status"] == "compatible",
                }
            )
            records.append(record)
    result = pd.DataFrame(records)
    for flag in ("adaptive_retained", "strict_retained"):
        smaller = set(result.loc[result.scope.eq(SCOPES[0]) & result[flag], "Ensembl_Gene_ID"])
        larger = set(result.loc[result.scope.eq(SCOPES[1]) & result[flag], "Ensembl_Gene_ID"])
        if not smaller <= larger:
            raise ValueError(f"Scope union lost candidates: {flag}")
    return result


def summarize_sources(membership, candidates):
    rows = []
    for doi, group in membership.groupby("doi", sort=True):
        mapped = set(group.Ensembl_Gene_ID) - {""}
        coding = set(group.loc[group.protein_candidate_eligible, "Ensembl_Gene_ID"])
        for scope in SCOPES:
            subset = candidates.loc[
                candidates.scope.eq(scope) & candidates.Ensembl_Gene_ID.isin(coding)
            ]
            rows.append(
                {
                    "doi": doi,
                    "source": PAPER_LABELS.get(doi, group.citation.iloc[0]),
                    "scope": scope,
                    "source_identities": group.identity.nunique(),
                    "mapped": len(mapped),
                    "protein_coding": len(coding),
                    **{
                        field: int(subset[field].sum())
                        for field in (
                            "family_eligible",
                            "adaptive_hpa_pass",
                            "adaptive_after_family",
                            "adaptive_retained",
                            "strict_after_family",
                            "strict_retained",
                        )
                    },
                }
            )
    return pd.DataFrame(rows)


def write_csv(frame, path):
    compression = {"method": "gzip", "mtime": 0} if str(path).endswith(".gz") else None
    frame.to_csv(path, index=False, compression=compression, float_format="%.6g")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    membership = cta_provenance.paper_membership().fillna("")
    membership["identity"] = membership.apply(source_identity, axis=1)
    ids = set(membership.loc[membership.protein_candidate_eligible, "Ensembl_Gene_ID"])
    raw = get_data("cancer-testis-antigens-extended")
    if not ids <= set(raw.Ensembl_Gene_ID):
        raise ValueError("Imported protein-coding nominations lack candidate rows")
    table = (
        raw.loc[raw.Ensembl_Gene_ID.isin(ids)].sort_values("Ensembl_Gene_ID").reset_index(drop=True)
    )
    tables = {modality: hpa._read_hpa(name, "v23") for modality, name in SOURCE_NAMES.items()}
    receipts = {name: reference_data.provenance(name, "v23") for name in SOURCE_NAMES.values()}
    for receipt in receipts.values():
        if sha256(receipt["path"]) != receipt["sha256"]:
            raise ValueError("Cached HPA source does not match its receipt")
    observations = {
        modality: annotate_measurements(frame, modality, ids) for modality, frame in tables.items()
    }
    expected = {modality: routine_contexts(frame, modality) for modality, frame in tables.items()}
    candidates = compare_candidates(
        table,
        observations,
        unreliable_ids=cta_regen.cross_reactive_ihc(),
        expected_contexts=expected,
    )
    reproductive = candidates.loc[candidates.scope.eq("reproductive")]
    if (
        set(reproductive.loc[reproductive.adaptive_retained, "Ensembl_Gene_ID"])
        != cta.cta_extended_gene_ids()
    ):
        raise ValueError("Adaptive reproductive baseline differs from the published extended panel")
    source_counts = summarize_sources(membership, candidates)
    summary = []
    for scope, group in candidates.groupby("scope", sort=False):
        record = {
            "scope": scope,
            "source_identities": membership.identity.nunique(),
            "mapped": membership.loc[
                membership.Ensembl_Gene_ID.ne(""), "Ensembl_Gene_ID"
            ].nunique(),
            "protein_coding": len(group),
        }
        record.update(
            {
                field: int(group[field].sum())
                for field in (
                    "family_eligible",
                    "adaptive_hpa_pass",
                    "adaptive_after_family",
                    "adaptive_retained",
                    "strict_after_family",
                    "strict_retained",
                )
            }
        )
        record["strict_retained_without_ihc"] = int(
            (group.strict_retained & group.ihc_coverage.eq("unavailable")).sum()
        )
        record["strict_retained_partial_ihc"] = int(
            (group.strict_retained & group.ihc_coverage.eq("partial")).sum()
        )
        summary.append(record)
    summary = pd.DataFrame(summary)
    gains = []
    for policy in ("adaptive", "strict"):
        pivot = candidates.pivot(
            index="Ensembl_Gene_ID", columns="scope", values=f"{policy}_retained"
        )
        gained = set(pivot.index[pivot[SCOPES[1]] & ~pivot[SCOPES[0]]])
        extra = candidates.loc[
            candidates.scope.eq(SCOPES[1]) & candidates.Ensembl_Gene_ID.isin(gained)
        ].copy()
        extra.insert(0, "policy", policy)
        gains.append(extra)
    gains = pd.concat(gains, ignore_index=True)
    exports = {
        "candidate-decisions.csv": candidates,
        "source-counts.csv": source_counts,
        "scope-summary.csv": summary,
        "disposable-gains.csv": gains,
        "source-membership.csv.gz": membership,
        "publication-sources.csv": cta_sources.publication_sources(),
        "normal-observations.csv.gz": pd.concat(observations.values(), ignore_index=True),
        "historical-only-candidates.csv": raw.loc[~raw.Ensembl_Gene_ID.isin(ids)],
    }
    for filename, frame in exports.items():
        write_csv(frame, out / filename)
    (out / "tissue-taxonomy.json").write_text(json.dumps(normal_tissue_taxonomy(), indent=2) + "\n")
    from cta_scope_report_render import render

    render(out, candidates, source_counts, summary, gains)
    manifest = {
        "report_schema": "oncoref.cta-scope-report.v1",
        "species": "human",
        "taxon_id": 9606,
        "oncoref_version": __version__,
        "taxonomy_id": normal_tissue_taxonomy()["taxonomy_id"],
        "policy_version": normal_tissue_taxonomy()["policy_version"],
        "source_receipts": {
            name: {
                key: value
                for key, value in receipt.items()
                if key not in {"path", "checksum_matches", "verification_state"}
            }
            for name, receipt in receipts.items()
        },
        "routine_contexts": {modality: sorted(contexts) for modality, contexts in expected.items()},
        "candidate_universe": "All registered published nomination sets; source scopes and targeted selections are preserved. No minimum-cover preselection.",
        "strict_rna_rule": "Every positive estimate (>0 nTPM), complete routine bulk RNA, no positive outside scope; raw IHC Low/Medium/High/Ascending/Descending considered.",
        "unreliable_ihc_rule": "Existing adaptive policy ignores reviewed cross-reactive IHC; strict sensitivity retains raw outside detections and otherwise marks unreliable IHC unresolved.",
        "single_cell_rna": "Not part of these two historical RNA/IHC filtering policies; unsited cell labels cannot establish complete anatomical coverage.",
        "baseline_parity": "Exact adaptive reproductive gene-set equality to cta_extended_gene_ids verified",
        "python": platform.python_version(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ("pandas", "numpy", "matplotlib", "reportlab", "pypdf")
        },
        "implementation_sha256": {
            name: sha256(ROOT / name)
            for name in (
                "scripts/cta_scope_report.py",
                "scripts/cta_scope_report_render.py",
                "oncoref/cta.py",
                "oncoref/cta_regen.py",
                "oncoref/normal_tissues.py",
            )
        },
        "packaged_input_sha256": {
            name: sha256(ROOT / "oncoref/data" / name)
            for name in (
                "cancer-testis-antigens-extended.csv",
                "cta-specificity-audit.csv",
                "cta-ihc-unreliable.csv",
                "cta-publication-sources.csv",
                "cta-publication-membership.csv",
                "normal-tissue-taxonomy.yaml",
                "normal-tissue-mappings.csv",
            )
        },
        "outputs": {
            p.name: sha256(p)
            for p in sorted(out.iterdir())
            if p.is_file() and p.name != "run-manifest.json"
        },
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(f"Report: {out / 'cta-tissue-scope-report.pdf'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    build(args.out)


if __name__ == "__main__":
    main()
