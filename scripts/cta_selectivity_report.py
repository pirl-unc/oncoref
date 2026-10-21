#!/usr/bin/env python3
"""Audit reproductive admission and tissue-by-tissue absolute/selectivity gates."""

from __future__ import annotations

import argparse
import json
import sys
from itertools import pairwise
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_permitted_report import (  # noqa: E402
    EXPANDED,
    annotate_archive,
    evaluate_gene,
    sha256,
    write_csv,
)
from cta_tissue_tier_report import TRACKED, checked_input  # noqa: E402

from oncoref import figure_style  # noqa: E402
from oncoref.cta_normal_screen import adaptive_normal_rna_limits  # noqa: E402
from oncoref.normal_tissues import (  # noqa: E402
    normal_tissue_screening_policy,
    normal_tissue_taxonomy,
)

PREVIOUS = ROOT / "docs/audits/cta-tissue-priorities-20261020"
OBSERVATIONS = ROOT / "docs/audits/cta-permitted-tissues-20261019"
TIERS = ("vital", "protected", "manageable")


def build(out):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory; historical reports are immutable")
    out.mkdir(parents=True, exist_ok=True)
    previous_path = checked_input(PREVIOUS, "candidate-decisions.csv")
    observation_path = checked_input(OBSERVATIONS, "normal-observations.csv.gz")
    prior = json.loads((PREVIOUS / "run-manifest.json").read_text())
    policy = normal_tissue_screening_policy("ihc_adaptive_cta")
    taxonomy = normal_tissue_taxonomy()
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
        fixed = audit.copy()
        fixed["rna_ceiling_ntpm"] = fixed.baseline_ceiling_ntpm
        fixed["minimum_reproductive_selectivity"] = fixed.screening_tier.map(
            {tier: policy[tier]["minimum_reproductive_selectivity"] for tier in TIERS}
        )
        fixed["selectivity_margin"] = (
            fixed.reproductive_selectivity - fixed.minimum_reproductive_selectivity
        )
        fixed["ihc_reliability_used"] = "Baseline (adaptation disabled)"
        audit["baseline_minimum_reproductive_selectivity"] = fixed.minimum_reproductive_selectivity
        audits.append(audit)
        for scenario, active, combined in (
            ("reproductive_required", audit, False),
            ("current_fixed", fixed, True),
            ("current_adaptive", audit, True),
        ):
            result = evaluate_gene(
                local,
                scope=EXPANDED,
                vital_limit=policy["vital"]["maximum_rna"],
                protected_limit=policy["protected"]["maximum_rna"],
                manageable_limit=policy["manageable"]["maximum_rna"],
                expected=expected,
                unreliable=unreliable,
                rna_limits=active.rna_ceiling_ntpm,
                selectivity_audit=active if combined else None,
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
    scenarios = {
        name: frame.set_index("Ensembl_Gene_ID") for name, frame in current.groupby("scenario")
    }
    adaptive = scenarios["current_adaptive"]
    chain = [scenarios["current_fixed"], adaptive, scenarios["reproductive_required"], old]
    for smaller, larger in pairwise(chain):
        assert set(smaller.index[smaller.retained]) <= set(larger.index[larger.retained])
    comparison = adaptive.copy()
    for column in (
        "retained",
        "evidence_pass",
        "rna_exclusions",
        "ihc_exclusions",
        "failure_reasons",
    ):
        comparison["previous_" + column] = old[column]
        comparison["reproductive_required_" + column] = scenarios["reproductive_required"][column]
    comparison["change"] = [
        "lost" if before and not now else "gained" if now and not before else "unchanged"
        for now, before in zip(comparison.retained, comparison.previous_retained)
    ]
    comparison["new_exclusion_stage"] = [
        "reproductive_evidence"
        if before and not positive
        else "selectivity"
        if positive and not now
        else ""
        for before, positive, now in zip(
            comparison.previous_retained,
            comparison.reproductive_required_retained,
            comparison.retained,
        )
    ]
    changes = comparison.loc[comparison.change.ne("unchanged")].sort_values(
        ["new_exclusion_stage", "Symbol"]
    )
    summary = []
    for name, frame in [
        ("previous_adaptive", old),
        *[
            (name, scenarios[name])
            for name in ("reproductive_required", "current_fixed", "current_adaptive")
        ],
    ]:
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
        "tissue_category_label",
        "screening_tier_label",
        "tissue_review_flags",
        "tissue_review_notes",
    ]
    for name, frame in {
        "candidate-decisions.csv": current,
        "candidate-comparison.csv": comparison.reset_index(),
        "membership-changes.csv": changes.reset_index(),
        "notable-targets.csv": comparison.loc[
            comparison.Symbol.isin(TRACKED | set(changes.Symbol))
        ].reset_index(),
        "rna-selectivity-audit.csv.gz": pd.concat(audits).reset_index(),
        "tissue-contexts.csv": observations[context_columns].drop_duplicates(),
        "scope-summary.csv": summary,
    }.items():
        write_csv(frame, out / name)
    (out / "tissue-taxonomy.json").write_text(json.dumps(taxonomy, indent=2) + "\n")
    render(out, summary, changes, policy)
    implementation = [
        "scripts/cta_selectivity_report.py",
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
                "report_schema": "oncoref.cta-selectivity-report.v1",
                "policy": policy,
                "taxonomy_id": taxonomy["taxonomy_id"],
                "routine_contexts": prior["routine_contexts"],
                "inputs": {
                    str(p.relative_to(ROOT)): sha256(p)
                    for p in [
                        previous_path,
                        observation_path,
                        PREVIOUS / "run-manifest.json",
                        OBSERVATIONS / "run-manifest.json",
                    ]
                },
                "checks": [
                    "Archived input hashes verified",
                    "Every candidate reevaluated from original measurements",
                    "Fixed combined is a subset of adaptive combined, reproductive-required, and previous retained sets",
                    "Selectivity and absolute ceilings use the same local IHC evidence; positive IHC remains an independent veto",
                ],
                "implementation_sha256": {name: sha256(ROOT / name) for name in implementation},
                "outputs": {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
            },
            indent=2,
        )
        + "\n"
    )
    print(summary.to_string(index=False))
    print(changes[["Symbol", "new_exclusion_stage", "selectivity_exclusions"]].to_string())


def render(out, summary, changes, policy):
    figure_style.apply()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    labels = [
        "v1.10.0\nabsolute screen",
        "Require reproductive\npositive evidence",
        "Combined gates\nfixed baseline",
        "Combined gates\nIHC adaptive",
    ]
    left = [0] * len(summary)
    for column, label, color in [
        ("without_ihc", "IHC unavailable", "#888888"),
        ("partial_ihc", "Partial routine IHC", "#377eb8"),
        ("complete_routine_ihc", "Complete routine IHC", "#237a57"),
    ]:
        values = summary[column].tolist()
        ax.barh(labels, values, left=left, label=label, color=color, height=0.55)
        left = [a + b for a, b in zip(left, values)]
    for i, total in enumerate(left):
        ax.text(total + 5, i, str(total), va="center")
    ax.invert_yaxis()
    ax.set_xlim(0, max(left) * 1.12)
    ax.set_xlabel("Retained published candidates")
    ax.set_title("Reproductive evidence, absolute RNA and tissue selectivity")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=8)
    fig.subplots_adjust(left=0.25, right=0.97, top=0.88, bottom=0.28)
    for suffix in ("png", "svg"):
        figure_style.save(fig, out / f"policy-comparison.{suffix}", keep_titles=True)
    plt.close(fig)

    def table(rows):
        return "\n".join("| " + " | ".join(map(str, row)) + " |" for row in rows)

    grid = table(
        [
            ["IHC confidence", "Vital", "Protected", "Manageable"],
            ["---", "---", "---", "---"],
            *[
                [
                    confidence,
                    *[
                        f"{policy['negative_ihc_rna_limits'][confidence][tier]:g} nTPM / {policy['negative_ihc_minimum_reproductive_selectivity'][confidence][tier]:.0%}"
                        for tier in TIERS
                    ],
                ]
                for confidence in ("Missing", "Uncertain", "Approved", "Supported", "Enhanced")
            ],
        ]
    )
    summary_table = table(
        [
            summary.columns,
            ["---"] * len(summary.columns),
            *summary.itertuples(index=False, name=None),
        ]
    )
    losses = (
        "\n".join(
            f"- **{stage}:** " + ", ".join(sorted(group.Symbol)) + "."
            for stage, group in changes.groupby("new_exclusion_stage")
        )
        or "- None."
    )
    (out / "index.md").write_text(
        "# CTA reproductive evidence and tissue selectivity — 20 October 2026\n\n"
        f"Policy `{policy['policy_id']}`. Same 2,537 published candidates and pinned human HPA v23 measurements as the v1.10.0 comparison.\n\n"
        "## Policy\n\n"
        "- **Reproductive:** positive evidence required: RNA ≥1 nTPM or positive IHC.\n"
        "- **Disposable:** expression is exempt from exclusion; it cannot supply reproductive evidence or inflate the selectivity reference.\n"
        "- **Manageable, Protected, Vital:** every measured tissue must pass both gates below. Any positive IHC, including Low, independently excludes.\n\n"
        "`Selectivity = Rmax / (Rmax + tissue RNA)`, where Rmax is the highest valid raw RNA in a mapped Reproductive tissue for that gene. "
        "No tissue sums, deflation or rounding before gating. This is pairwise expression selectivity, not the percentage of total-body RNA or protein. "
        "Zero or missing reproductive RNA cannot establish selectivity, even with positive reproductive IHC. Missing routine RNA and unusable data continue to block retention.\n\n"
        "## Maximum nTPM / minimum selectivity\n\n" + grid + "\n\n"
        "Both limits are inclusive. The same local IHC evidence selects both thresholds: complete routine Not detected staining in the RNA organ, using the weakest reliability. "
        "Other organs and narrower regions cannot supply this evidence. Missing, uncertain, incomplete or unusable IHC uses baseline thresholds. "
        "These are user-approved research parameters, not validated therapeutic safety limits.\n\n"
        "## Results\n\n" + summary_table + "\n\n"
        "The reproductive-required scenario isolates the admission-rule change. The combined adaptive scenario then adds selectivity; "
        "the fixed scenario uses both baseline gates without IHC relaxation. Current retained sets are nested as checked in the manifest.\n\n"
        "![Retained candidates and IHC coverage](policy-comparison.png)\n\n"
        "### Newly excluded against v1.10.0\n\n" + losses + "\n\n"
        "## Coverage and interpretation\n\n"
        "The panel contains 50 routine RNA contexts and 80 routine IHC tissue/cell pairs. This is not exhaustive anatomy: cartilage, whole eye and inner ear lack routine RNA measurements here. "
        "Missing IHC is explicit and earns no relaxation. Adult HPA data does not certify maternal/fetal profiles. Human numerical limits do not transfer to dog or mouse. "
        "Family and nomination/manual-review gates, including CTAG2’s holdback, are unchanged. The frozen public CTA panel and its older global deflated fraction remain separate.\n\n"
        "[Definitions and source mappings](../../normal-tissue-scopes.md), [literature and limitations](../../normal-tissue-manageability.md), "
        "[v1.10.0 baseline](../cta-tissue-priorities-20261020/index.md). HPA also uses expression enrichment against other tissues in its "
        "[specificity definitions](https://v23.proteinatlas.org/humanproteome/tissue/tissue%2Bspecific); our reproductive reference and numerical grid are project definitions.\n\n"
        "## Downloads\n\n"
        "- [All candidate decisions and exclusion causes](candidate-comparison.csv)\n"
        "- [Newly excluded candidates](membership-changes.csv)\n"
        "- [Named targets and changed candidates](notable-targets.csv)\n"
        "- [All current scenarios](candidate-decisions.csv)\n"
        "- [Every RNA value, reproductive reference, selectivity and supporting IHC](rna-selectivity-audit.csv.gz)\n"
        "- [Native tissue mappings](tissue-contexts.csv)\n"
        "- [Policy, anatomy and bibliography](tissue-taxonomy.json)\n"
        "- [Checksums and provenance](run-manifest.json)\n\n"
        "Reproduce with `python scripts/cta_selectivity_report.py --out <empty-directory>` at this report’s implementation revision.\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
