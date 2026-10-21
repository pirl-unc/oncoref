#!/usr/bin/env python3
"""Reproduce the published-candidate screen with vital/protected RNA ceilings."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_scope_report import annotate_measurements, sha256, write_csv  # noqa: E402

from oncoref.normal_tissues import (  # noqa: E402
    normal_tissue_screening_details,
    normal_tissue_screening_policy,
    normal_tissue_taxonomy,
    resolve_normal_tissue,
)

ARCHIVE = ROOT / "docs/audits/cta-tissue-scopes-20261019"
EXPANDED = "reproductive_or_permitted_normal"
COMBINED = "reproductive_or_disposable"
SCOPES = ("reproductive", COMBINED, EXPANDED)
SCENARIOS = {
    "previous_1": (COMBINED, 1.0, 1.0),
    "expanded_1": (EXPANDED, 1.0, 1.0),
    "tiered": (EXPANDED, 1.0, 5.0),
}
TARGETS = {
    "PRAME": "PRAME",
    "CT83": "KK-LC-1",
    "CTAG1B": "NY-ESO-1",
    "CTAG2": "LAGE-1",
    "MAGEA1": "MAGE-A1",
    "MAGEA3": "MAGE-A3",
    "MAGEA4": "MAGE-A4",
    "MAGEA10": "MAGE-A10",
    "MAGEA11": "MAGE-A11",
    "MAGEC1": "MAGE-C1 / CT7",
    "MAGEC2": "MAGE-C2",
    "XAGE1B": "XAGE1B",
    "XAGE1A": "XAGE1A",
    "SSX1": "SSX1",
    "SSX2": "SSX2",
    "ROPN1": "ROPN1",
    "ROPN1B": "ROPN1B",
    "LDHC": "LDHC",
    "HORMAD1": "HORMAD1",
    "HORMAD2": "HORMAD2",
    "PAGE4": "PAGE4",
    "PAGE2B": "PAGE2B",
    "FATE1": "FATE1",
    "LIN28B": "LIN28B",
    "PIWIL1": "PIWIL1",
    "TDRD1": "TDRD1",
    "BRDT": "BRDT",
    "CABYR": "CABYR",
    "EGFL6": "EGFL6",
    "GAGE12F": "GAGE12F",
    "GAGE12G": "GAGE12G",
    "GAGE12J": "GAGE12J",
    "GAGE13": "GAGE13",
    "MAGEB6": "MAGEB6",
    "CT45A8": "CT45A8",
    "CT45A10": "CT45A10",
    "DDX43": "HAGE",
    "CTCFL": "BORIS",
    "AKAP4": "AKAP4",
    "PLAC1": "PLAC1",
    "PASD1": "PASD1",
    "LY6K": "LY6K",
}


def annotate_archive(observations):
    """Re-resolve original human HPA labels using the current cited taxonomy."""
    frames = []
    ids = set(observations.Ensembl_Gene_ID)
    for modality, group in observations.groupby("modality", sort=False):
        raw = group.rename(
            columns={
                "Ensembl_Gene_ID": "Gene",
                "Symbol": "Gene name",
                "tissue": "Tissue",
                "cell_type": "Cell type",
                "rna_ntpm": "nTPM",
                "ihc_level": "Level",
                "ihc_reliability": "Reliability",
            }
        )
        columns = ["Gene", "Gene name", "Tissue"]
        columns += ["nTPM"] if modality == "bulk_rna" else ["Cell type", "Level", "Reliability"]
        if modality == "bulk_rna":
            raw["nTPM"] = pd.to_numeric(raw.nTPM, errors="coerce")
        frames.append(annotate_measurements(raw[columns], modality, ids, scopes=SCOPES))
    result = pd.concat(frames, ignore_index=True)
    contexts = []
    for row in (
        result[["source_name", "source_version", "tissue", "cell_type"]]
        .drop_duplicates()
        .itertuples(index=False)
    ):
        annotation = resolve_normal_tissue(
            row.tissue or None,
            species="human",
            cell_type=row.cell_type or None,
            source_id=row.source_name,
            source_version=row.source_version,
        )
        details = normal_tissue_screening_details(annotation)
        contexts.append(
            {
                **row._asdict(),
                "tissue_category": details["tissue_category"],
                "tissue_category_label": details["tissue_category_label"],
                "screening_tier": details["screening_tier"],
                "screening_tier_label": details["screening_tier_label"],
                "tissue_review_flags": ";".join(details["review_flags"]),
                "tissue_review_notes": " ".join(
                    review["rationale"] for review in details["context_reviews"].values()
                ),
            }
        )
    return result.merge(pd.DataFrame(contexts), validate="many_to_one")


def evaluate_gene(
    rows,
    *,
    scope,
    vital_limit,
    protected_limit,
    expected,
    unreliable=False,
    rna_limits=None,
    manageable_limit=None,
    selectivity_audit=None,
):
    """Separate RNA ceilings and a strict outside-IHC gate, with explicit missingness."""
    config = normal_tissue_screening_policy()
    rna = rows.loc[rows.modality.eq("bulk_rna")]
    ihc = rows.loc[rows.modality.eq("ihc")]
    mapped = rows.mapping_status.eq("mapped")
    allowed = mapped & rows[scope]
    reproductive = allowed & rows[config["required_positive_scope"]] & rows.valid_measurement
    tiers = rows.screening_tier.copy()
    # Scope-only comparison: held-out requested exceptions receive the protected
    # ceiling; scenarios with the narrower scope use the same 1 nTPM ceiling in both tiers.
    tiers.loc[tiers.eq("permitted") & ~allowed] = "protected"
    tiers.loc[allowed] = "permitted"
    limits = tiers.map(
        {
            "vital": vital_limit,
            "protected": protected_limit,
            "manageable": config["manageable"]["maximum_rna"]
            if manageable_limit is None
            else manageable_limit,
        }
    ).astype(float)
    if rna_limits is not None:
        if not rna_limits.index.equals(rna.index):
            raise ValueError("RNA ceilings must match every RNA row in order")
        if not rna.empty:
            limits.loc[rna.index] = rna_limits
    forbidden_ihc = rows.loc[rows.modality.eq("ihc") & rows.strict_positive & mapped & ~allowed]
    above = rows.loc[
        rows.modality.eq("bulk_rna")
        & mapped
        & tiers.isin(["vital", "protected", "manageable"])
        & rows.rna_ntpm.gt(limits)
    ]
    positive_inside = bool(
        (
            reproductive
            & rows.modality.eq("bulk_rna")
            & rows.rna_ntpm.gt(config["minimum_reproductive_rna"])
        ).any()
        or (reproductive & rows.modality.eq("ihc") & rows.strict_positive).any()
    )
    seen = {
        name: set(
            zip(
                frame.loc[frame.valid_measurement, "tissue"],
                frame.loc[frame.valid_measurement, "cell_type"],
            )
        )
        for name, frame in (("bulk_rna", rna), ("ihc", ihc))
    }
    missing = {name: set(expected[name]) - seen[name] for name in expected}
    failures, unknown = [], []
    selectivity_failures = []
    selectivity_minima = {}
    if selectivity_audit is not None:
        if not selectivity_audit.index.equals(rna.index):
            raise ValueError("Selectivity audit must match every RNA row in order")
        if not selectivity_audit.Ensembl_Gene_ID.eq(rna.Ensembl_Gene_ID).fillna(
            False
        ).all() or not (selectivity_audit.tissue.eq(rna.tissue).fillna(False).all()):
            raise ValueError("Selectivity audit must match the gene and tissue observations")
        selected = selectivity_audit.loc[
            tiers.loc[rna.index].isin(["vital", "protected", "manageable"])
        ]
        if selected.empty or selected.reproductive_selectivity.isna().any():
            unknown.append("unresolved_reproductive_selectivity")
        if selected.minimum_reproductive_selectivity.isna().any():
            raise ValueError("Screened tissues require selectivity thresholds")
        below = selected.loc[
            selected.reproductive_selectivity.lt(selected.minimum_reproductive_selectivity)
        ]
        if len(below):
            failures.append("reproductive_selectivity_below_threshold")
        selectivity_failures = [
            f"{r.tissue}: {r.reproductive_selectivity:.4%} (minimum "
            f"{r.minimum_reproductive_selectivity:.0%}, {r.ihc_reliability_used})"
            for r in below.itertuples(index=False)
        ]
        selectivity_minima = {
            f"{tier}_worst_selectivity_margin": group.reproductive_selectivity.sub(
                group.minimum_reproductive_selectivity
            ).min()
            for tier, group in selected.groupby("screening_tier")
        }
    if len(above):
        failures.append("rna_above_tissue_ceiling")
    if len(forbidden_ihc):
        failures.append("positive_ihc_outside_permitted")
    if not positive_inside:
        failures.append("insufficient_positive_reproductive_evidence")
    if not expected["bulk_rna"] or missing["bulk_rna"]:
        unknown.append("incomplete_routine_rna")
    if (~rows.valid_measurement).any():
        unknown.append("unusable_measurement")
    if (rows.strict_positive & (~mapped | tiers.eq("unresolved"))).any():
        unknown.append("unresolved_positive_anatomy")
    if unreliable:
        unknown.append("reviewed_unreliable_ihc")
    ihc_coverage = (
        "unreliable"
        if unreliable
        else "unavailable"
        if not seen["ihc"]
        else "partial"
        if missing["ihc"]
        else "routine_complete"
    )
    review_rows = rows.loc[rows.strict_positive & rows.tissue_review_flags.ne("")]
    return {
        "evidence_pass": not failures and not unknown,
        "evidence_status": "excluded" if failures else "unresolved" if unknown else "compatible",
        "failure_reasons": ";".join(failures + unknown),
        "positive_evidence_inside": positive_inside,
        "positive_evidence_reproductive": positive_inside,
        "selectivity_gate_applied": selectivity_audit is not None,
        "selectivity_exclusions": "; ".join(selectivity_failures),
        **selectivity_minima,
        "rna_routine_missing": len(missing["bulk_rna"]),
        "ihc_routine_missing": len(missing["ihc"]),
        "ihc_coverage": ihc_coverage,
        "vital_rna_max_ntpm": rows.loc[
            rows.modality.eq("bulk_rna") & tiers.eq("vital"), "rna_ntpm"
        ].max(),
        "protected_rna_max_ntpm": rows.loc[
            rows.modality.eq("bulk_rna") & tiers.eq("protected"), "rna_ntpm"
        ].max(),
        "manageable_rna_max_ntpm": rows.loc[
            rows.modality.eq("bulk_rna") & tiers.eq("manageable"), "rna_ntpm"
        ].max(),
        "rna_exclusions": "; ".join(
            f"{r.tissue}: {r.rna_ntpm:g} nTPM ({tiers.loc[i]}, limit {limits.loc[i]:g})"
            for i, r in above.sort_values("rna_ntpm", ascending=False).iterrows()
        ),
        "ihc_exclusions": "; ".join(
            f"{r.context_label}: {r.ihc_level} ({r.ihc_reliability or 'reliability unavailable'})"
            for r in forbidden_ihc.itertuples(index=False)
        ),
        "ihc_forbidden_contexts_observed": len(
            set(
                zip(
                    rows.loc[
                        rows.modality.eq("ihc") & rows.valid_measurement & mapped & ~allowed,
                        "tissue",
                    ],
                    rows.loc[
                        rows.modality.eq("ihc") & rows.valid_measurement & mapped & ~allowed,
                        "cell_type",
                    ],
                )
            )
        ),
        "absence_evidence": "not_established",
        "tissue_review_status": "required" if len(review_rows) else "not_triggered",
        "tissue_review_flags": ";".join(
            sorted({flag for value in review_rows.tissue_review_flags for flag in value.split(";")})
        ),
        "tissue_review_contexts": "; ".join(sorted(set(review_rows.context_label))),
    }


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    archive_manifest = json.loads((ARCHIVE / "run-manifest.json").read_text())
    inputs = (
        "candidate-decisions.csv",
        "normal-observations.csv.gz",
        "source-membership.csv.gz",
        "source-counts.csv",
    )
    for name in inputs:
        if sha256(ARCHIVE / name) != archive_manifest["outputs"][name]:
            raise ValueError(f"Archived input checksum mismatch: {name}")
    base = pd.read_csv(ARCHIVE / inputs[0], keep_default_na=False)
    base = base.loc[base.scope.eq(COMBINED)].set_index("Ensembl_Gene_ID")
    observations = annotate_archive(
        pd.read_csv(ARCHIVE / inputs[1], keep_default_na=False, low_memory=False)
    )
    expected = {
        key: {tuple(context) for context in contexts}
        for key, contexts in archive_manifest["routine_contexts"].items()
    }
    groups = dict(iter(observations.groupby("Ensembl_Gene_ID")))
    records = []
    for gid, row in base.iterrows():
        for scenario, (scope, vital, protected) in SCENARIOS.items():
            record = evaluate_gene(
                groups.get(gid, observations.iloc[:0]),
                scope=scope,
                vital_limit=vital,
                protected_limit=protected,
                expected=expected,
                unreliable=row.reviewed_unreliable_ihc,
            )
            record.update(
                {
                    "Ensembl_Gene_ID": gid,
                    "Symbol": row.Symbol,
                    "scenario": scenario,
                    "scope": scope,
                    "vital_rna_ceiling_ntpm": vital,
                    "protected_rna_ceiling_ntpm": protected,
                    "family_eligible": row.family_eligible,
                    "nomination_review_eligible": row.nomination_review_eligible,
                    "nomination_review_reason": row.nomination_review_reason,
                    "adaptive_combined_retained": row.adaptive_retained,
                    "retained": bool(
                        record["evidence_pass"]
                        and row.family_eligible
                        and row.nomination_review_eligible
                    ),
                }
            )
            if not row.family_eligible:
                record["failure_reasons"] += ";family_exclusion"
            if not row.nomination_review_eligible:
                record["failure_reasons"] += ";existing_nomination_or_manual_review"
            record["failure_reasons"] = record["failure_reasons"].strip(";")
            records.append(record)
    candidates = pd.DataFrame(records)
    sets = {
        scenario: set(group.loc[group.retained, "Ensembl_Gene_ID"])
        for scenario, group in candidates.groupby("scenario")
    }
    if not sets["previous_1"] <= sets["expanded_1"] <= sets["tiered"]:
        raise ValueError("Expansion or ceiling relaxation lost candidates")
    # The previous sensitivity is an independently computed, documented anchor.
    if len(sets["previous_1"]) != 547:
        raise ValueError("Previous 1 nTPM combined-scope sensitivity no longer reproduces 547")
    summary = []
    for scenario in SCENARIOS:
        group = candidates.loc[candidates.scenario.eq(scenario)]
        retained = group.loc[group.retained]
        summary.append(
            {
                "scenario": scenario,
                "protein_coding_candidates": len(group),
                "evidence_pass": int(group.evidence_pass.sum()),
                "retained": len(retained),
                "without_ihc": int(retained.ihc_coverage.eq("unavailable").sum()),
                "partial_ihc": int(retained.ihc_coverage.eq("partial").sum()),
                "complete_routine_ihc": int(retained.ihc_coverage.eq("routine_complete").sum()),
                "no_forbidden_ihc_contexts_observed": int(
                    retained.ihc_forbidden_contexts_observed.eq(0).sum()
                ),
                "retained_from_adaptive_combined": int(retained.adaptive_combined_retained.sum()),
            }
        )
    summary = pd.DataFrame(summary)
    membership = pd.read_csv(ARCHIVE / "source-membership.csv.gz", keep_default_na=False)
    sources = pd.read_csv(ARCHIVE / "source-counts.csv")
    sources = sources.loc[sources.scope.eq(COMBINED)].copy()
    doi_labels = sources.set_index("doi").source.to_dict()
    citation_groups = (
        membership.loc[membership.Ensembl_Gene_ID.ne("")]
        .groupby("Ensembl_Gene_ID")
        .doi.agg(lambda values: sorted(set(values)))
    )
    candidates["source_dois"] = candidates.Ensembl_Gene_ID.map(citation_groups).map(
        lambda values: ";".join(values)
    )
    candidates["source_publications"] = candidates.Ensembl_Gene_ID.map(citation_groups).map(
        lambda values: ";".join(doi_labels.get(doi, doi) for doi in values)
    )
    for scenario, ids in sets.items():
        counts = (
            membership.loc[membership.Ensembl_Gene_ID.isin(ids)]
            .groupby("doi")
            .Ensembl_Gene_ID.nunique()
        )
        sources[scenario] = sources.doi.map(counts).fillna(0).astype(int)
    targets = candidates.loc[
        candidates.scenario.eq("tiered") & candidates.Symbol.isin(TARGETS)
    ].copy()
    targets["name"] = targets.Symbol.map(TARGETS)
    for scenario, ids in sets.items():
        targets[f"retained_{scenario}"] = targets.Ensembl_Gene_ID.isin(ids)
    targets["reason"] = targets.apply(
        lambda r: (
            "Included; IHC coverage: " + r.ihc_coverage
            if r.retained
            else "; ".join(
                filter(
                    None,
                    [
                        r.rna_exclusions,
                        r.ihc_exclusions,
                        r.failure_reasons,
                        r.nomination_review_reason if not r.nomination_review_eligible else "",
                    ],
                )
            )
        ),
        axis=1,
    )
    targets = targets.set_index("Symbol").reindex(TARGETS).reset_index()
    if targets.Ensembl_Gene_ID.isna().any():
        raise ValueError("A named published target is absent from the source universe")
    exports = {
        "candidate-decisions.csv": candidates,
        "scope-summary.csv": summary,
        "source-counts.csv": sources,
        "notable-targets.csv": targets,
        "normal-observations.csv.gz": observations,
        "tissue-contexts.csv": observations[
            [
                "source_name",
                "source_version",
                "tissue",
                "cell_type",
                "concept_id",
                "mapping_status",
                "mapping_references",
                "screening_tier",
                *SCOPES,
            ]
        ].drop_duplicates(),
    }
    for name, frame in exports.items():
        write_csv(frame, out / name)
    (out / "tissue-policy.json").write_text(json.dumps(normal_tissue_taxonomy(), indent=2) + "\n")
    from cta_permitted_report_render import render

    render(out, candidates, sources, summary, targets)
    (out / "run-manifest.json").write_text(
        json.dumps(
            {
                "report_schema": "oncoref.cta-permitted-report.v1",
                "taxonomy_id": normal_tissue_taxonomy()["taxonomy_id"],
                "policy": normal_tissue_screening_policy(),
                "input_report": str(ARCHIVE.relative_to(ROOT)),
                "input_manifest_sha256": sha256(ARCHIVE / "run-manifest.json"),
                "inputs": {name: sha256(ARCHIVE / name) for name in inputs},
                "scenarios": SCENARIOS,
                "routine_contexts": archive_manifest["routine_contexts"],
                "manual_reviews": "Existing family and nomination/manual-review decisions preserved and exposed separately from the new expression gates.",
                "checks": [
                    "547-candidate previous-sensitivity parity",
                    "nested scope expansion and ceiling relaxation",
                ],
                "implementation_sha256": {
                    name: sha256(ROOT / name)
                    for name in (
                        "scripts/cta_permitted_report.py",
                        "scripts/cta_permitted_report_render.py",
                        "scripts/cta_scope_report.py",
                        "oncoref/normal_tissues.py",
                        "oncoref/data/normal-tissue-taxonomy.yaml",
                        "oncoref/data/normal-tissue-mappings.csv",
                    )
                },
                "outputs": {
                    path.name: sha256(path)
                    for path in sorted(out.iterdir())
                    if path.is_file() and path.name != "run-manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(summary.to_string(index=False))
    print(targets[["Symbol", "retained", "reason"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
