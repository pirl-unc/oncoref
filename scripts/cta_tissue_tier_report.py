#!/usr/bin/env python3
"""Audit tissue reclassification and the category-by-IHC grid against a frozen baseline."""

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

PREVIOUS = ROOT / "docs/audits/cta-tissue-contexts-20261020"
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
        # Apply the archived ladder to the newly classified anatomy. Caution
        # receives the former Protected ladder in this anatomy-only comparator.
        audit["previous_grid_ceiling_ntpm"] = [
            prior["policy"]["negative_ihc_rna_limits"][r.ihc_reliability_used][
                "protected" if r.screening_tier == "caution" else r.screening_tier
            ]
            if r.screening_tier in {"vital", "protected", "caution"}
            else float("nan")
            for r in audit.itertuples(index=False)
        ]
        audits.append(audit)
        for scenario in ("reclassified_previous_limits", "current_fixed", "current_adaptive"):
            limits = {
                "reclassified_previous_limits": audit.previous_grid_ceiling_ntpm,
                "current_fixed": None,
                "current_adaptive": audit.rna_ceiling_ntpm,
            }[scenario]
            result = evaluate_gene(
                local,
                scope=EXPANDED,
                vital_limit=policy["vital"]["maximum_rna"],
                protected_limit=policy["protected"]["maximum_rna"],
                caution_limit=policy["caution"]["maximum_rna"],
                expected=expected,
                unreliable=unreliable,
                rna_limits=limits,
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
    reclassified = current.loc[current.scenario.eq("reclassified_previous_limits")].set_index(
        "Ensembl_Gene_ID"
    )
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
        comparison["reclassified_" + column] = reclassified[column]
    comparison["change"] = [
        "gained" if now and not before else "lost" if before and not now else "unchanged"
        for now, before in zip(comparison.retained, comparison.previous_retained)
    ]
    changes = comparison.loc[comparison.change.ne("unchanged")].sort_values(["change", "Symbol"])
    summary = []
    for name, frame in (
        ("previous_adaptive", old),
        ("reclassified_previous_limits", reclassified),
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
        "rna-threshold-audit.csv.gz": pd.concat(audits).reset_index(),
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
                    "Anatomy-only comparator uses archived Vital/Protected ladders; Caution uses archived Protected",
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

    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    labels = [
        "Previous anatomy\nadaptive",
        "New anatomy\nprevious RNA grid",
        "New anatomy\nfixed 1 / 2 / 5",
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
    ax.set_title("Tissue reclassification and category-specific IHC/RNA thresholds")
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
    grid = "\n".join(
        "| " + " | ".join(map(str, row)) + " |"
        for row in [
            ["IHC evidence", "Vital", "Protected", "Caution"],
            ["---", "---:", "---:", "---:"],
            *[
                [
                    key,
                    *[
                        policy["negative_ihc_rna_limits"][key][tier]
                        for tier in ("vital", "protected", "caution")
                    ],
                ]
                for key in ("Missing", "Uncertain", "Approved", "Supported", "Enhanced")
            ],
        ]
    )
    (out / "index.md").write_text(
        "# CTA tissue priorities and IHC-dependent thresholds — 20 October 2026\n\n"
        f"Policy `{policy['policy_id']}`. All {len(adaptive):,} published candidates are "
        "reevaluated from the same pinned human HPA v23 observations.\n\n"
        "## Tissue categories\n\n"
        "- **Reproductive:** existing reproductive/developmental unions.\n"
        "- **Disposable:** cutaneous melanocytes, follicles, tonsil, appendix and gallbladder.\n"
        "- **Caution:** whole skin, spleen, thymus, thyroid and adipose.\n"
        "- **Protected:** includes salivary glands, cartilage, esophagus, other protected "
        "gastrointestinal tissues and lymph nodes.\n"
        "- **Vital / critical function:** existing highest-protection organs.\n\n"
        "Caution is screened: it is not a permitted normal-expression exception. "
        "Specific permitted skin compartments take precedence over the whole-skin category. "
        "A bulk skin measurement cannot establish compartment restriction.\n\n"
        "## Maximum RNA by tissue category and IHC evidence (nTPM)\n\n" + grid + "\n\n"
        "Approved/Supported/Enhanced relaxation requires valid **Not detected** staining "
        "in every routine IHC context of the same organ; the weakest reliability wins. "
        "Missing, incomplete, uncertain or unusable IHC earns no relaxation. "
        "Any positive IHC, including Low, in Caution, Protected or Vital excludes. "
        "Reliability labels are not calibrated probabilities of true protein absence. "
        "These are project parameters, not clinical safety limits or GTEx-to-HPA conversions.\n\n"
        "## Results and attribution\n\n" + summary_table + "\n\n"
        "The **reclassified_previous_limits** scenario changes tissue assignments while "
        "retaining the archived RNA grid (new Caution receives the previous Protected grid). "
        "Comparing it with **current_adaptive** isolates the RNA-grid change. "
        "The fixed screen uses baseline ceilings without IHC-based relaxation.\n\n"
        "![Retained candidates and IHC coverage](policy-comparison.png)\n\n"
        f"**Gained against the previous release:** {gained}.\n\n"
        f"**Lost against the previous release:** {lost}.\n\n"
        f"**{flagged} retained candidates have positive observations triggering tissue review.** "
        "Flags record RNA >0 nTPM or positive IHC in flagged contexts; they do not "
        "establish completed therapeutic acceptance or independently veto retention.\n\n"
        "## Evidence and interpretation\n\n"
        "[Tissue definitions](../../normal-tissue-scopes.md) and the "
        "[literature review](../../normal-tissue-manageability.md) explain why published "
        "Caution/Hazard categories support relative priorities, while exact memberships "
        "and thresholds remain project definitions. [Native contexts](tissue-contexts.csv) "
        "retain source labels, mappings and citations.\n\n"
        "An organ allowance does not allow the same immune-cell population elsewhere. "
        "Unknown anatomy remains unresolved. Human numeric thresholds and tolerability "
        "do not transfer to dog or mouse. Routine coverage is 50 RNA contexts and 80 IHC "
        "tissue/cell pairs, not exhaustive anatomy. Whole eye, inner ear and cartilage lack "
        "routine RNA measurements here; classification does not establish their negativity. "
        "The [patient-context filters](../../normal-tissue-gene-filters.md) require separate "
        "maternal/fetal and sex-specific coverage absent from this adult HPA report. "
        "Existing family and manual-curation gates, including CTAG2's holdback, are retained.\n\n"
        "## Downloads and reproduction\n\n"
        "- [All candidate comparisons and exact exclusion reasons](candidate-comparison.csv)\n"
        "- [Gained/lost candidates](membership-changes.csv)\n"
        "- [Named published targets](notable-targets.csv)\n"
        "- [All three current-anatomy scenarios](candidate-decisions.csv)\n"
        "- [Every RNA ceiling and supporting IHC](rna-threshold-audit.csv.gz)\n"
        "- [Taxonomy, review conditions and bibliography](tissue-taxonomy.json)\n"
        "- [Run manifest and checksums](run-manifest.json)\n\n"
        "Run `python scripts/cta_tissue_tier_report.py --out <new-empty-directory>` "
        "at this report's implementation revision. The "
        "[previous tissue-context report](../cta-tissue-contexts-20261020/index.md) is the "
        "immutable baseline.\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
