#!/usr/bin/env python3
"""Compare fixed and IHC-conditioned absolute RNA ceilings on published CTAs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_permitted_report import (  # noqa: E402
    EXPANDED,
    TARGETS,
    annotate_archive,
    evaluate_gene,
    sha256,
    write_csv,
)

from oncoref.cta_normal_screen import adaptive_normal_rna_limits  # noqa: E402
from oncoref.normal_tissues import normal_tissue_screening_policy  # noqa: E402

ARCHIVE = ROOT / "docs/audits/cta-permitted-tissues-20261019"
ORIGIN = ROOT / "docs/audits/cta-tissue-scopes-20261019"


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ARCHIVE / "run-manifest.json").read_text())
    input_names = ("candidate-decisions.csv", "normal-observations.csv.gz", "source-counts.csv")
    for name in input_names:
        if sha256(ARCHIVE / name) != manifest["outputs"][name]:
            raise ValueError(f"Archived input hash mismatch: {name}")
    original = json.loads((ORIGIN / "run-manifest.json").read_text())
    if (
        sha256(ORIGIN / "source-membership.csv.gz")
        != original["outputs"]["source-membership.csv.gz"]
    ):
        raise ValueError("Publication membership hash mismatch")
    base = pd.read_csv(ARCHIVE / input_names[0], keep_default_na=False)
    base = base.loc[base.scenario.eq("tiered")].set_index("Ensembl_Gene_ID")
    observations = annotate_archive(
        pd.read_csv(ARCHIVE / input_names[1], keep_default_na=False, low_memory=False)
    )
    expected = {
        key: {tuple(pair) for pair in contexts}
        for key, contexts in manifest["routine_contexts"].items()
    }
    groups = dict(iter(observations.groupby("Ensembl_Gene_ID")))
    decisions, audits = [], []
    for number, (gid, row) in enumerate(base.iterrows(), start=1):
        if number == 1 or number % 250 == 0:
            print(f"Evaluating candidate {number}/{len(base)}", flush=True)
        local = groups.get(gid, observations.iloc[:0])
        audit = adaptive_normal_rna_limits(
            local,
            routine_ihc_contexts=expected["ihc"],
            reviewed_unreliable=row.ihc_coverage == "unreliable",
        )
        audit["Symbol"] = row.Symbol
        audits.append(audit)
        for scenario in ("fixed", "adaptive"):
            result = evaluate_gene(
                local,
                scope=EXPANDED,
                vital_limit=1,
                protected_limit=5,
                expected=expected,
                unreliable=row.ihc_coverage == "unreliable",
                rna_limits=audit.rna_ceiling_ntpm if scenario == "adaptive" else None,
            )
            if scenario == "fixed" and result["evidence_pass"] != row.evidence_pass:
                raise ValueError(f"Fixed expression gate changed for {gid}")
            result.update(
                {
                    "Ensembl_Gene_ID": gid,
                    "Symbol": row.Symbol,
                    "scenario": scenario,
                    "retained": bool(
                        result["evidence_pass"]
                        and row.family_eligible
                        and row.nomination_review_eligible
                    ),
                    "family_eligible": row.family_eligible,
                    "nomination_review_eligible": row.nomination_review_eligible,
                    "nomination_review_reason": row.nomination_review_reason,
                    "source_dois": row.source_dois,
                    "source_publications": row.source_publications,
                    "relaxed_tissue_count": int(audit.relaxed.sum())
                    if scenario == "adaptive"
                    else 0,
                    "rescued_rna_contexts": "; ".join(
                        f"{r.tissue}: {r.rna_ntpm:g} <= {r.rna_ceiling_ntpm:g} nTPM; "
                        f"{r.ihc_evidence_status}"
                        for r in audit.loc[
                            audit.rna_ntpm.gt(audit.baseline_ceiling_ntpm)
                            & audit.rna_ntpm.le(audit.rna_ceiling_ntpm)
                        ].itertuples()
                    )
                    if scenario == "adaptive"
                    else "",
                }
            )
            if not row.family_eligible:
                result["failure_reasons"] += ";family_exclusion"
            if not row.nomination_review_eligible:
                result["failure_reasons"] += ";existing_nomination_or_manual_review"
            result["failure_reasons"] = result["failure_reasons"].strip(";")
            if scenario == "fixed" and result["retained"] != row.retained:
                raise ValueError(f"Fixed final decision changed for {gid}")
            decisions.append(result)
    candidates = pd.DataFrame(decisions)
    sets = {
        key: set(group.loc[group.retained, "Ensembl_Gene_ID"])
        for key, group in candidates.groupby("scenario")
    }
    if len(sets["fixed"]) != 629 or not sets["fixed"] <= sets["adaptive"]:
        raise ValueError("Fixed 629-candidate parity or nested relaxation failed")
    audits = pd.concat(audits)
    summary = []
    for key in ("fixed", "adaptive"):
        rows = candidates.loc[candidates.scenario.eq(key)]
        retained = rows.loc[rows.retained]
        summary.append(
            {
                "scenario": key,
                "candidates": len(rows),
                "evidence_pass": int(rows.evidence_pass.sum()),
                "retained": len(retained),
                "without_ihc": int(retained.ihc_coverage.eq("unavailable").sum()),
                "partial_ihc": int(retained.ihc_coverage.eq("partial").sum()),
                "complete_routine_ihc": int(retained.ihc_coverage.eq("routine_complete").sum()),
                "recovered": len(sets[key] - sets["fixed"]),
            }
        )
    summary = pd.DataFrame(summary)
    sources = pd.read_csv(ARCHIVE / "source-counts.csv")
    sources["fixed"] = sources.tiered
    membership = pd.read_csv(ORIGIN / "source-membership.csv.gz", keep_default_na=False)
    counts = (
        membership.loc[membership.Ensembl_Gene_ID.isin(sets["adaptive"])]
        .groupby("doi")
        .Ensembl_Gene_ID.nunique()
    )
    sources["adaptive"] = sources.doi.map(counts).fillna(0).astype(int)
    targets = candidates.loc[
        candidates.scenario.eq("adaptive") & candidates.Symbol.isin(TARGETS)
    ].copy()
    targets["name"] = targets.Symbol.map(TARGETS)
    targets["fixed_retained"] = targets.Ensembl_Gene_ID.isin(sets["fixed"])
    targets = targets.set_index("Symbol").reindex(TARGETS).reset_index()
    if targets.Ensembl_Gene_ID.isna().any():
        raise ValueError("Named target missing from candidate universe")
    for name, frame in {
        "candidate-decisions.csv": candidates,
        "rna-threshold-audit.csv.gz": audits.reset_index(drop=True),
        "scope-summary.csv": summary,
        "source-counts.csv": sources,
        "notable-targets.csv": targets,
        "recovered-candidates.csv": candidates.loc[
            candidates.scenario.eq("adaptive")
            & candidates.Ensembl_Gene_ID.isin(sets["adaptive"] - sets["fixed"])
        ],
    }.items():
        write_csv(frame, out / name)
    policy = normal_tissue_screening_policy("ihc_adaptive_cta")
    (out / "screening-policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    from cta_adaptive_report_render import render

    render(out, summary, sources, targets)
    implementation = (
        "scripts/cta_adaptive_report.py",
        "scripts/cta_adaptive_report_render.py",
        "scripts/cta_permitted_report.py",
        "scripts/cta_scope_report.py",
        "oncoref/cta_normal_screen.py",
        "oncoref/normal_tissues.py",
        "oncoref/data/normal-tissue-taxonomy.yaml",
        "oncoref/data/normal-tissue-mappings.csv",
    )
    (out / "run-manifest.json").write_text(
        json.dumps(
            {
                "report_schema": "oncoref.cta-ihc-adaptive-report.v1",
                "policy": policy,
                "input_report": str(ARCHIVE.relative_to(ROOT)),
                "input_manifest_sha256": sha256(ARCHIVE / "run-manifest.json"),
                "inputs": {name: sha256(ARCHIVE / name) for name in input_names},
                "publication_membership_sha256": sha256(ORIGIN / "source-membership.csv.gz"),
                "routine_contexts": manifest["routine_contexts"],
                "checks": [
                    "All fixed expression and final decisions match",
                    "629 fixed candidates",
                    "Fixed retained set is a subset of adaptive",
                ],
                "implementation_sha256": {name: sha256(ROOT / name) for name in implementation},
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
    print(
        targets[
            ["Symbol", "fixed_retained", "retained", "rna_exclusions", "rescued_rna_contexts"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
