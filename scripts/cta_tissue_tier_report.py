#!/usr/bin/env python3
"""Audit tissue reclassification against the frozen preceding adaptive CTA screen."""

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

from oncoref import figure_style  # noqa: E402
from oncoref.cta_normal_screen import adaptive_normal_rna_limits  # noqa: E402
from oncoref.normal_tissues import (  # noqa: E402
    normal_tissue_screening_policy,
    normal_tissue_taxonomy,
)

PREVIOUS = ROOT / "docs/audits/cta-tissue-tiers-20261020"
OBSERVATIONS = ROOT / "docs/audits/cta-permitted-tissues-20261019"
TRACKED = set(TARGETS) | {
    "TEX15",
    "CCDC172",
    "CDRT15",
    "CIB4",
    "CSAG3",
    "DCST1",
    "DNAH17",
    "FSIP2",
    "KLK15",
    "MS4A6E",
    "NLRP4",
    "NOXRED1",
    "RPL10L",
    "TKTL1",
    "UBE2U",
    "ZPBP",
}


def checked_input(directory, name):
    manifest = json.loads((directory / "run-manifest.json").read_text())
    path = directory / name
    if sha256(path) != manifest["outputs"][name]:
        raise ValueError(f"Archived input hash mismatch: {path}")
    return path


def build(out):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory; historical reports are immutable")
    out.mkdir(parents=True, exist_ok=True)
    previous_path = checked_input(PREVIOUS, "candidate-decisions.csv")
    observation_path = checked_input(OBSERVATIONS, "normal-observations.csv.gz")
    prior = json.loads((PREVIOUS / "run-manifest.json").read_text())
    policy = normal_tissue_screening_policy("ihc_adaptive_cta")
    if policy["negative_ihc_rna_limits"] != prior["policy"]["negative_ihc_rna_limits"]:
        raise ValueError("This comparison isolates anatomy changes; RNA ladders must match")
    old = pd.read_csv(previous_path, keep_default_na=False)
    old = old.loc[old.scenario.eq("current_adaptive")].set_index("Ensembl_Gene_ID")
    observations = annotate_archive(
        pd.read_csv(observation_path, keep_default_na=False, low_memory=False)
    )
    expected = {key: set(map(tuple, pairs)) for key, pairs in prior["routine_contexts"].items()}
    groups = dict(iter(observations.groupby("Ensembl_Gene_ID")))
    decisions, audits = [], []
    for number, (gid, row) in enumerate(old.iterrows(), start=1):
        if number == 1 or number % 250 == 0:
            print(f"Evaluating {number}/{len(old)}", flush=True)
        local = groups.get(gid, observations.iloc[:0])
        unreliable = row.ihc_coverage == "unreliable"
        audit = adaptive_normal_rna_limits(
            local, routine_ihc_contexts=expected["ihc"], reviewed_unreliable=unreliable
        )
        audit["Symbol"] = row.Symbol
        audits.append(audit)
        for scenario in ("current_fixed", "current_adaptive"):
            result = evaluate_gene(
                local,
                scope=EXPANDED,
                vital_limit=policy["vital"]["maximum_rna"],
                protected_limit=policy["protected"]["maximum_rna"],
                expected=expected,
                unreliable=unreliable,
                rna_limits=audit.rna_ceiling_ntpm if scenario == "current_adaptive" else None,
            )
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
                    **{
                        key: row[key]
                        for key in (
                            "family_eligible",
                            "nomination_review_eligible",
                            "nomination_review_reason",
                            "source_dois",
                            "source_publications",
                        )
                    },
                }
            )
            for key, reason in (
                ("family_eligible", "family_exclusion"),
                ("nomination_review_eligible", "existing_nomination_or_manual_review"),
            ):
                if not row[key]:
                    result["failure_reasons"] = ";".join(
                        filter(None, [result["failure_reasons"], reason])
                    )
            decisions.append(result)
    current = pd.DataFrame(decisions)
    fixed = current.loc[current.scenario.eq("current_fixed")]
    adaptive = current.loc[current.scenario.eq("current_adaptive")].set_index("Ensembl_Gene_ID")
    if not set(fixed.loc[fixed.retained, "Ensembl_Gene_ID"]) <= set(
        adaptive.index[adaptive.retained]
    ):
        raise ValueError("Fixed retained set must be a subset of the current adaptive set")
    comparison = adaptive.copy()
    for column in (
        "retained",
        "evidence_pass",
        "rna_exclusions",
        "ihc_exclusions",
        "failure_reasons",
    ):
        comparison["previous_" + column] = old[column]
    comparison["change"] = [
        "gained" if now and not before else "lost" if before and not now else "unchanged"
        for now, before in zip(comparison.retained, comparison.previous_retained)
    ]
    changes = comparison.loc[comparison.change.ne("unchanged")].sort_values(["change", "Symbol"])
    summary = []
    for name, frame in (
        ("previous_adaptive", old),
        ("current_fixed", fixed),
        ("current_adaptive", adaptive),
    ):
        kept = frame.loc[frame.retained]
        summary.append(
            {
                "scenario": name,
                "candidates": len(frame),
                "evidence_pass": int(frame.evidence_pass.sum()),
                "retained": len(kept),
                "without_ihc": int(kept.ihc_coverage.eq("unavailable").sum()),
                "partial_ihc": int(kept.ihc_coverage.eq("partial").sum()),
                "complete_routine_ihc": int(kept.ihc_coverage.eq("routine_complete").sum()),
            }
        )
    summary = pd.DataFrame(summary)
    context_columns = [
        "source_name",
        "source_version",
        "tissue",
        "cell_type",
        "concept_id",
        "mapping_status",
        "mapping_kind",
        "mapping_references",
        "screening_tier",
        "tissue_category",
        "screening_tier_label",
        "tissue_review_flags",
        "tissue_review_notes",
    ]
    for name, frame in {
        "candidate-decisions.csv": current,
        "candidate-comparison.csv": comparison.reset_index(),
        "membership-changes.csv": changes.reset_index(),
        "notable-targets.csv": comparison.loc[comparison.Symbol.isin(TRACKED)].reset_index(),
        "rna-threshold-audit.csv.gz": pd.concat(audits).reset_index(drop=True),
        "tissue-contexts.csv": observations[context_columns].drop_duplicates(),
        "scope-summary.csv": summary,
    }.items():
        write_csv(frame, out / name)
    taxonomy = normal_tissue_taxonomy()
    (out / "tissue-taxonomy.json").write_text(json.dumps(taxonomy, indent=2) + "\n")
    render(out, summary, changes, adaptive, policy)
    write_manifest(out, policy, taxonomy, prior, previous_path, observation_path)
    print(summary.to_string(index=False))
    print(
        changes[["Symbol", "change", "rna_exclusions", "previous_rna_exclusions"]].to_string(
            index=False
        )
    )


def write_manifest(out, policy, taxonomy, prior, previous_path, observation_path):
    implementation = [
        "scripts/cta_tissue_tier_report.py",
        "scripts/cta_permitted_report.py",
        "scripts/cta_scope_report.py",
        "oncoref/cta_normal_screen.py",
        "oncoref/normal_tissues.py",
        "oncoref/data/normal-tissue-taxonomy.yaml",
        "oncoref/data/normal-tissue-mappings.csv",
    ]
    (out / "run-manifest.json").write_text(
        json.dumps(
            {
                "report_schema": "oncoref.cta-tissue-tier-report.v1",
                "policy": policy,
                "taxonomy_id": taxonomy["taxonomy_id"],
                "routine_contexts": prior["routine_contexts"],
                "inputs": {
                    str(p.relative_to(ROOT)): sha256(p)
                    for p in (
                        previous_path,
                        observation_path,
                        PREVIOUS / "run-manifest.json",
                        OBSERVATIONS / "run-manifest.json",
                    )
                },
                "checks": [
                    "Archived input hashes verified",
                    "RNA ladder identical to previous adaptive report",
                    "Every candidate reevaluated from original measurements",
                    "Current fixed retained set is a subset of current adaptive",
                ],
                "implementation_sha256": {name: sha256(ROOT / name) for name in implementation},
                "outputs": {
                    p.name: sha256(p)
                    for p in sorted(out.iterdir())
                    if p.is_file() and p.name != "run-manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )


def render(out, summary, changes, adaptive, policy):
    figure_style.apply()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 3.5))
    labels = [
        "Previous anatomy\nadaptive",
        "Current anatomy\nfixed 1 / 5",
        "Current anatomy\nadaptive",
    ]
    left = [0] * len(summary)
    for column, label, color in (
        ("without_ihc", "IHC unavailable", "#888888"),
        ("partial_ihc", "Partial routine IHC", "#377eb8"),
        ("complete_routine_ihc", "Complete routine IHC", "#237a57"),
    ):
        values = summary[column].tolist()
        ax.barh(labels, values, left=left, label=label, color=color, height=0.55)
        left = [a + b for a, b in zip(left, values)]
    for i, total in enumerate(left):
        ax.text(total + 5, i, str(total), va="center")
    ax.invert_yaxis()
    ax.set_xlim(0, max(left) * 1.12)
    ax.set_xlabel("Retained published candidates")
    ax.set_title("Tissue assignments change; adaptive RNA ceilings stay fixed")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=8)
    fig.subplots_adjust(left=0.23, right=0.97, top=0.88, bottom=0.28)
    for suffix in ("png", "svg"):
        figure_style.save(fig, out / f"policy-comparison.{suffix}", keep_titles=True)
    plt.close(fig)
    gained = ", ".join(changes.loc[changes.change.eq("gained"), "Symbol"]) or "None"
    lost = ", ".join(changes.loc[changes.change.eq("lost"), "Symbol"]) or "None"
    flagged = int((adaptive.retained & adaptive.tissue_review_status.eq("required")).sum())
    summary_table = "\n".join(
        "| " + " | ".join(map(str, row)) + " |"
        for row in [
            summary.columns,
            ["---"] * len(summary.columns),
            *summary.itertuples(index=False, name=None),
        ]
    )
    (out / "index.md").write_text(
        "# CTA tissue contexts and removable organs — 20 October 2026\n\n"
        f"Policy `{policy['policy_id']}`. All 2,537 published candidates are reevaluated "
        "from the same pinned human HPA v23 observations.\n\n"
        "## Tissue changes\n\n"
        "- Appendix and gallbladder join the conditional **Disposable** category.\n"
        "- Native observations expose **Reproductive, Disposable, Protected, Vital** separately.\n"
        "- Skin and salivary allowances gain explicit extent/control review annotations.\n"
        "- Earlier organ allowances remain project choices, pending the manageability-policy review.\n"
        "- [Patient-context gene filters](../../normal-tissue-gene-filters.md) require separate "
        "sex/stage-aware negative coverage; this adult HPA report cannot certify those filters.\n\n"
        "The adaptive vital/protected ceilings remain **1/5, 2/5, 3/7.5 and 5/10 nTPM** "
        "for missing/uncertain, Approved, Supported and Enhanced negative IHC. "
        "The proposed 1/2/3/6 and 2/4/6/12 ladder is not adopted in this comparison.\n\n"
        "## Results\n\n" + summary_table + "\n\n"
        "![Retained candidates and IHC coverage](policy-comparison.png)\n\n"
        f"**Gained:** {gained}.\n\n**Lost:** {lost}.\n\n"
        f"**{flagged} retained candidates have positive observations requiring tissue review.** "
        "This is a research-screen result, not completed therapeutic acceptance. Review flags "
        "record RNA >0 nTPM or positive IHC in contexts with review conditions; "
        "they do not independently change expression-pass or final-retention decisions.\n\n"
        "## Evidence and interpretation\n\n"
        "The full [tissue definitions and citations](../../normal-tissue-scopes.md) explain "
        "conditional organ allowances and their limits. Native tissue/cell labels and "
        "citations remain in [tissue-contexts.csv](tissue-contexts.csv). An organ allowance "
        "does not allow expression in the same immune-cell population elsewhere. Unsited "
        "immune cells remain unresolved. Mouse has no mapped tonsil; canine and mouse "
        "anatomical assignments do not transfer human thresholds or clinical tolerability.\n\n"
        "Any positive IHC outside acceptable tissues, including Low, still excludes. "
        "RNA relaxation requires valid negative staining across every routine cell context "
        "of the same anatomical tissue. Missing/partial IHC is explicit. Complete routine "
        "coverage means 50 RNA contexts and 80 IHC tissue/cell pairs, not exhaustive anatomy. "
        "Whole eye and inner ear have no routine RNA measurements in this input; assigning "
        "their tier does not establish absence of expression in those organs. "
        "Existing family and manual curation gates remain separate. CTAG2's historical "
        "manual holdback is retained. Negative observations do not prove antigen absence.\n\n"
        "## Downloads and reproduction\n\n"
        "- [All candidate comparisons and reasons](candidate-comparison.csv)\n"
        "- [Gained/lost candidates](membership-changes.csv)\n"
        "- [Named targets, including the preceding sensitivity-analysis genes](notable-targets.csv)\n"
        "- [Current fixed and adaptive decisions](candidate-decisions.csv)\n"
        "- [Every RNA value, ceiling and supporting IHC](rna-threshold-audit.csv.gz)\n"
        "- [Taxonomy, review conditions and bibliography](tissue-taxonomy.json)\n"
        "- [Run manifest and checksums](run-manifest.json)\n\n"
        "Run `python scripts/cta_tissue_tier_report.py --out <new-empty-directory>` "
        "at this report's implementation revision. The preceding "
        "[tissue-tier report](../cta-tissue-tiers-20261020/index.md) remains an immutable baseline.\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
