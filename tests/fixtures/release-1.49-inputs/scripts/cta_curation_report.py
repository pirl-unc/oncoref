#!/usr/bin/env python3
"""Publish a readable current CTA report from genome-wide expression decisions."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_clinical_curation_report import clinical_disposition, read_verified  # noqa: E402
from cta_permitted_report import annotate_archive  # noqa: E402
from cta_report_pages import render  # noqa: E402
from cta_scope_report import annotate_measurements, sha256, write_csv  # noqa: E402

from oncoref import cta_regen, cta_sources, hpa, reference_data  # noqa: E402
from oncoref.cta_normal_screen import adaptive_normal_rna_limits  # noqa: E402
from oncoref.normal_tissues import normal_tissue_screening_policy  # noqa: E402

EVIDENCE = {
    "registered_intervention": "Trial registration",
    "approved_therapy": "Approved therapy",
    "shared_epitope_clinical": "Shared target in a clinical study",
    "family_label_clinical": "Antigen family named in a trial",
    "primary_clinical_report": "Published clinical study",
    "institutional_clinical_report": "Institutional treatment report",
    "clinical_safety_context": "Published safety findings",
}
RESULT = {
    "retained": "Included",
    "normal_protein_counterevidence": "Normal protein detected",
    "rna_threshold_review": "RNA above limit",
    "selectivity_reference_review": "Selectivity below limit",
    "coverage_or_curation_review": "Needs review",
    "external_safety_hold": "Safety hold",
}
REASON = {
    "rna_above_tissue_ceiling": "RNA exceeds a tissue's limit",
    "positive_ihc_outside_permitted": "Normal-tissue protein staining was detected",
    "reproductive_selectivity_below_threshold": "Reproductive expression is not sufficiently higher",
    "insufficient_positive_reproductive_evidence": "No measured reproductive expression",
    "incomplete_routine_rna": "Required RNA measurements are missing",
    "unusable_measurement": "An expression measurement is unusable",
    "unresolved_positive_anatomy": "A positive sample has unresolved tissue anatomy",
    "reviewed_unreliable_ihc": "The antibody evidence was reviewed as unreliable",
    "unresolved_reproductive_selectivity": "Reproductive RNA cannot establish selectivity",
    "family_exclusion": "Gene family is outside the curated CTA definition",
    "existing_nomination_or_manual_review": "Publication nomination or earlier review remains unresolved",
    "external_safety_hold": "Published normal-tissue safety evidence requires a hold",
}


def readable_reason(reasons):
    return "; ".join(
        REASON.get(reason, reason.replace("_", " ")) for reason in reasons.split(";") if reason
    )


def source_counts(current):
    """Count unique identities within each original source; overlaps are preserved."""
    membership = cta_sources.reviewed_publication_membership()
    mapped = set(current.Ensembl_Gene_ID)
    retained = set(current.loc[current.retained, "Ensembl_Gene_ID"])
    rows = []
    for source in cta_sources.publication_sources().itertuples():
        members = membership.loc[membership.source_tag.eq(source.source_tag)]
        ids = set(members.Ensembl_Gene_ID) - {""}
        rows.append(
            {
                "Publication / list": source.citation,
                "Source rows": source.published_rows,
                "Distinct identities": members.intake_identity.nunique(),
                "Mapped coding candidates": len(ids & mapped),
                "Included": len(ids & retained),
                "Source": source.source_url,
                "Exact table": source.source_table,
                "Evidence scope": source.evidence_scope,
            }
        )
    return pd.DataFrame(rows)


def threshold_comparison(out, current, categories, policy):
    """Recheck tissue-level RNA against raw pinned HPA and compare fixed/adaptive gates."""
    manifest = json.loads((categories / "run-manifest.json").read_text())
    expected = set(map(tuple, manifest["routine_contexts"]["ihc"]))
    parts, receipts = [], {}
    ids = set(current.Ensembl_Gene_ID)
    for modality, name in [("bulk_rna", "hpa_rna_consensus"), ("ihc", "hpa_normal_tissue")]:
        receipt = reference_data.provenance(name, "v23")
        assert sha256(receipt["path"]) == receipt["sha256"]
        assert receipt["sha256"] == manifest["sources"][name]["sha256"]
        receipts[name] = {key: receipt[key] for key in ["url", "sha256", "bytes", "version"]}
        parts.append(annotate_measurements(hpa._read_hpa(name, "v23"), modality, ids))
    observations = annotate_archive(pd.concat(parts, ignore_index=True))
    grouped = dict(iter(observations.groupby("Ensembl_Gene_ID", sort=False)))
    unreliable = cta_regen.cross_reactive_ihc()
    audits, comparisons = [], []
    tiers = ["vital", "protected", "manageable"]
    for row in current.itertuples():
        audit = adaptive_normal_rna_limits(
            grouped.get(row.Ensembl_Gene_ID, observations.iloc[:0]),
            routine_ihc_contexts=expected,
            reviewed_unreliable=row.Ensembl_Gene_ID in unreliable,
        )
        adaptive_rna = audit.rna_ntpm.gt(audit.rna_ceiling_ntpm)
        adaptive_selectivity = audit.reproductive_selectivity.lt(
            audit.minimum_reproductive_selectivity
        )
        failures = set(row.failure_reasons.split(";"))
        assert bool(adaptive_rna.any()) == ("rna_above_tissue_ceiling" in failures)
        assert bool(row.reproductive_positive and adaptive_selectivity.any()) == (
            "reproductive_selectivity_below_threshold" in failures
        )
        fixed_limit = audit.screening_tier.map(
            {tier: policy[tier]["maximum_rna"] for tier in tiers}
        )
        fixed_share = audit.screening_tier.map(
            {tier: policy[tier]["minimum_reproductive_selectivity"] for tier in tiers}
        )
        fixed_rna = audit.rna_ntpm.gt(fixed_limit)
        fixed_selectivity = audit.reproductive_selectivity.lt(fixed_share)
        zero_outside = audit.loc[audit.screening_tier.isin(tiers), "rna_ntpm"].gt(0).any()
        fixed_pass = row.retained and not fixed_rna.any() and not fixed_selectivity.any()
        comparisons.append(
            {
                "Gene": row.Symbol,
                "Gene ID": row.Ensembl_Gene_ID,
                "Fixed thresholds": bool(fixed_pass),
                "IHC-adaptive thresholds": row.retained,
                "Zero measured outside RNA": bool(row.retained and not zero_outside),
                "Fixed RNA findings": "; ".join(audit.loc[fixed_rna, "tissue"]),
                "Fixed selectivity findings": "; ".join(audit.loc[fixed_selectivity, "tissue"]),
            }
        )
        audit["Symbol"] = row.Symbol
        audits.append(audit)
    comparison = pd.DataFrame(comparisons)
    write_csv(comparison, out / "threshold-comparison.csv")
    write_csv(pd.concat(audits, ignore_index=True), out / "tissue-evidence.csv.gz")
    summary = pd.DataFrame(
        [
            {"Screen": label, "Included": int(comparison[label].sum())}
            for label in [
                "Fixed thresholds",
                "IHC-adaptive thresholds",
                "Zero measured outside RNA",
            ]
        ]
    )
    write_csv(summary, out / "threshold-comparison-summary.csv")
    return summary, receipts


def build(out, categories, reference_checks):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    previous_dir = ROOT / "docs/audits/cta-curation-20261021"
    old = read_verified(previous_dir, "candidate-decisions.csv.gz")
    intake = read_verified(previous_dir, "intake-decisions.csv.gz")
    genes = read_verified(categories, "gene-categories.csv.gz")
    metadata = old[
        [
            "Ensembl_Gene_ID",
            "family_eligible",
            "nomination_review_eligible",
            "nomination_review_reason",
            "source_dois",
            "retained",
        ]
    ].rename(columns={"retained": "previous_retained"})
    current = metadata.merge(genes, on="Ensembl_Gene_ID", how="left", validate="one_to_one")
    assert len(current) == len(old) and current.category.notna().all()
    current["evidence_pass"] = current.category.eq("Reproductive")
    current["pre_safety_retained"] = (
        current.evidence_pass & current.family_eligible & current.nomination_review_eligible
    )
    safety = pd.read_csv(
        ROOT / "oncoref/data/cta-external-safety-reviews.csv", keep_default_na=False
    )
    current["external_safety_hold"] = current.Ensembl_Gene_ID.isin(
        safety.loc[safety.action.eq("hold"), "Ensembl_Gene_ID"]
    )
    current["retained"] = current.pre_safety_retained & ~current.external_safety_hold
    reasons = []
    for row in current.itertuples():
        parts = [row.failure_reasons, row.quality_flags]
        if not row.reproductive_positive:
            parts.append("insufficient_positive_reproductive_evidence")
        if not row.family_eligible:
            parts.append("family_exclusion")
        if not row.nomination_review_eligible:
            parts.append("existing_nomination_or_manual_review")
        if row.external_safety_hold:
            parts.append("external_safety_hold")
        reasons.append(";".join(filter(None, parts)))
    current["failure_reasons"] = reasons
    current["Result"] = [RESULT[clinical_disposition(row)] for row in current.itertuples()]
    current["Reason"] = current.failure_reasons.map(readable_reason)
    current["Change"] = [
        "Newly included"
        if r.retained and not r.previous_retained
        else "No longer included"
        if r.previous_retained and not r.retained
        else "Unchanged"
        for r in current.itertuples()
    ]
    clinical = pd.read_csv(ROOT / "oncoref/data/cta-trial-evidence.csv", keep_default_na=False)
    clinical["Clinical evidence"] = clinical.evidence_kind.map(EVIDENCE)
    assert clinical["Clinical evidence"].notna().all()
    clinical = clinical.merge(current, on=["Symbol", "Ensembl_Gene_ID"], validate="one_to_one")
    assert len(clinical) == 30
    current["Clinical evidence"] = current.Ensembl_Gene_ID.map(
        clinical.set_index("Ensembl_Gene_ID")["Clinical evidence"]
    ).fillna("")
    decisions = current.set_index("Ensembl_Gene_ID")
    for name in ["retained", "failure_reasons"]:
        intake[name] = [
            decisions.loc[gid, name] if gid in decisions.index else old_value
            for gid, old_value in zip(intake.Ensembl_Gene_ID, intake[name])
        ]
    stages = [
        ("Source nominations", len(intake)),
        ("Mapped protein-coding candidates", len(current)),
    ]
    mask = pd.Series(True, index=current.index)
    failed = current.failure_reasons.str.split(";").map(set)
    for name, excluded in [
        ("Reproductive expression measured", {"insufficient_positive_reproductive_evidence"}),
        (
            "Usable evidence and RNA coverage",
            {
                "incomplete_routine_rna",
                "unusable_measurement",
                "unresolved_positive_anatomy",
                "reviewed_unreliable_ihc",
                "unresolved_reproductive_selectivity",
            },
        ),
        ("No disallowed normal protein staining", {"positive_ihc_outside_permitted"}),
        ("RNA within tissue limits", {"rna_above_tissue_ceiling"}),
        ("Reproductive selectivity passes", {"reproductive_selectivity_below_threshold"}),
        (
            "Family and nomination review passes",
            {"family_exclusion", "existing_nomination_or_manual_review"},
        ),
        ("Included after safety review", {"external_safety_hold"}),
    ]:
        mask &= ~failed.map(lambda value, excluded=excluded: bool(value & excluded))
        stages.append((name, int(mask.sum())))
    assert mask.equals(current.retained)
    funnel = pd.DataFrame(stages, columns=["Step", "Genes"])
    human = current[
        [
            "Symbol",
            "Ensembl_Gene_ID",
            "category",
            "Result",
            "Reason",
            "rna_exclusions",
            "ihc_exclusions",
            "selectivity_exclusions",
            "Clinical evidence",
            "Change",
        ]
    ].rename(
        columns={
            "Symbol": "Gene",
            "Ensembl_Gene_ID": "Gene ID",
            "category": "Expression group",
            "rna_exclusions": "RNA findings",
            "ihc_exclusions": "Protein findings",
            "selectivity_exclusions": "Selectivity findings",
        }
    )
    for name, frame in {
        "candidate-decisions.csv.gz": current,
        "gene-results.csv.gz": human,
        "clinical-target-audit.csv": clinical,
        "intake-decisions.csv.gz": intake,
        "funnel.csv": funnel,
        "changes.csv": human.loc[human.Change.ne("Unchanged")],
    }.items():
        write_csv(frame, out / name)
    checks = {p.stem: json.loads(p.read_text()) for p in sorted(reference_checks.glob("*.json"))}
    (out / "ensembl-identity-checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    policy = normal_tissue_screening_policy("ihc_adaptive_cta")
    (out / "screening-policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    sources = source_counts(current)
    write_csv(sources, out / "publication-results.csv")
    comparison, receipts = threshold_comparison(out, current, categories, policy)
    render(
        out, current, intake, clinical, genes, funnel, EVIDENCE, RESULT, policy, sources, comparison
    )
    inputs = [
        previous_dir / "candidate-decisions.csv.gz",
        previous_dir / "intake-decisions.csv.gz",
        categories / "gene-categories.csv.gz",
        categories / "run-manifest.json",
        ROOT / "oncoref/data/cta-trial-evidence.csv",
        ROOT / "oncoref/data/cta-external-safety-reviews.csv",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/cta-publication-membership.csv",
        ROOT / "oncoref/data/cta-publication-sources.csv",
        ROOT / "oncoref/data/cta-intake-identity-reviews.csv",
    ]
    manifest = {
        "schema": "oncoref.cta-curation-report.v2",
        "policy_id": policy["policy_id"],
        "sources": receipts,
        "inputs": {str(p.resolve().relative_to(ROOT)): sha256(p) for p in inputs},
        "implementation_sha256": {
            name: sha256(ROOT / name)
            for name in [
                "scripts/cta_curation_report.py",
                "scripts/cta_report_pages.py",
                "scripts/cta_scope_report.py",
                "scripts/cta_permitted_report.py",
                "scripts/cta_clinical_curation_report.py",
                "oncoref/cta_normal_screen.py",
                "oncoref/cta_sources.py",
            ]
        },
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(funnel.to_string(index=False))
    print(
        current.loc[
            current.Symbol.isin(["XAGE1A", "XAGE1B", "LY6K"]), ["Symbol", "Result", "Reason"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["out", "categories", "reference-checks"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    build(args.out, args.categories, args.reference_checks)
