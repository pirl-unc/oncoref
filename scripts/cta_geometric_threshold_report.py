#!/usr/bin/env python3
"""Reassess the complete pinned normal-RNA evidence with a two-axis doubling grid."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import sha256 as digest  # noqa: E402
from cta_report_common import write_csv as write  # noqa: E402

from oncoref import cell_rna_policy as policy  # noqa: E402

MODES = ["legacy", "ihc_doubling", "two_axis_doubling"]
ASSAYS = {
    "hpa_single_cell_tissue": "HPA tissue clusters",
    "hpa_immune_cell": "HPA immune",
    "hpa_immune_monaco": "Monaco immune",
    "dice": "DICE immune",
    "cornea_collin_adult_counts": "Adult cornea",
}
NAMED = [
    "PRAME",
    "MAGEA4",
    "CTAG1A",
    "CTAG1B",
    "XAGE1A",
    "XAGE1B",
    "CTAG2",
    "SSX2",
    "CT83",
    "TBXT",
    "SMC1B",
    "PAGE1",
]


def cell_assessment(
    cells: pd.DataFrame, mode: str, *, max_ihc_grade: str = "Enhanced"
) -> pd.DataFrame:
    d = cells.copy()
    levels = {"Missing/uncertain": 0, **policy.GRADES}
    ranks = d.negative_ihc_grade.map(levels)
    if ranks.isna().any() or max_ihc_grade not in levels:
        raise ValueError("Unknown cell IHC grade")
    d["negative_ihc_grade"] = ranks.clip(upper=levels[max_ihc_grade]).map(dict(enumerate(levels)))
    d = policy.assess(
        d,
        policy.linked_thresholds(mode=mode),
        ihc_factors=dict.fromkeys(levels, 1),
        inclusive_ceiling=True,
    )
    d["local_exclusion"] = d.cta_excluded | ~d.mapping_status.eq("mapped")
    d.loc[d.category.eq(policy.NEUTRAL_REPRODUCTIVE_IMMUNE), "local_exclusion"] = False
    return d


def normal_passes(genes: pd.DataFrame, bulk: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    """Return one aligned candidate-pass mask per normal RNA assay."""
    result = pd.DataFrame(index=genes.index)
    result["Normal bulk"] = genes.gene_id.map(
        bulk.groupby("Ensembl_Gene_ID").local_exclusion.any()
    ).eq(False)
    for assay, label in ASSAYS.items():
        result[label] = genes.gene_id.map(
            cells[cells.assay.eq(assay)].groupby("gene_id").local_exclusion.any()
        ).eq(False)
    return result


def eligibility(genes: pd.DataFrame, tumor: str = "maximum10") -> pd.Series:
    tumor_pass = {
        "recurrence": genes.tumor_positive,
        "maximum10": genes.tumor_maximum_max_ptpm.ge(10),
        "maximum30": genes.tumor_maximum_max_ptpm.ge(30),
    }[tumor]
    return genes.positive_context & tumor_pass & ~genes.ihc_detected & genes.rna_coverage_complete


def funnel(
    genes: pd.DataFrame, gates: pd.DataFrame, *, mode: str
) -> tuple[pd.DataFrame, pd.Series]:
    """Count sequential removals and retain each gene's first failed gate."""
    steps = [
        ("Reference genes", pd.Series(True, index=genes.index)),
        ("Reproductive or embryonic association", genes.positive_context),
        ("Tumor sample ≥10 pTPM", genes.tumor_maximum_max_ptpm.ge(10)),
        ("No detected nonpermitted IHC", ~genes.ihc_detected),
        ("Complete required RNA panels", genes.rna_coverage_complete),
        *[(col, gates[col]) for col in gates],
    ]
    keep = pd.Series(True, index=genes.index)
    first_failure = pd.Series("Selected", index=genes.index)
    rows = []
    for number, (label, mask) in enumerate(steps):
        before = int(keep.sum())
        first_failure.loc[keep & ~mask] = label
        keep &= mask
        rows.append(
            {
                "mode": mode,
                "step": number,
                "filter": label,
                "before": before,
                "removed": before - int(keep.sum()),
                "remaining": int(keep.sum()),
            }
        )
    return pd.DataFrame(rows), first_failure


def changed_contexts(old: pd.DataFrame, new: pd.DataFrame, *, kind: str) -> pd.DataFrame:
    id_col = "Ensembl_Gene_ID" if kind == "bulk" else "gene_id"
    value = "rna_ntpm" if kind == "bulk" else "value"
    cutoff = "effective_ceiling_ntpm" if kind == "bulk" else "cta_cutoff"
    columns = [id_col, "tissue", value]
    if kind == "cell":
        columns += ["cell_type", "assay", "unit", "category", "negative_ihc_grade"]
    else:
        columns += ["screening_tier", "negative_ihc_grade"]
    # Both assessments preserve the same row ordering from the pinned input.
    assert old[id_col].equals(new[id_col])
    changed = old.local_exclusion.ne(new.local_exclusion)
    r = new.loc[changed, columns].rename(columns={id_col: "gene_id", value: "value"}).copy()
    r["kind"] = kind
    r["previous_ceiling"] = old.loc[changed, cutoff]
    r["new_ceiling"] = new.loc[changed, cutoff]
    r["previous_excluded"] = old.loc[changed, "local_exclusion"]
    r["new_excluded"] = new.loc[changed, "local_exclusion"]
    return r


def build(source, out: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    manifest = json.loads((source / "run-manifest.json").read_text())
    for name in ["gene-comparison.csv.gz", "bulk-contexts.csv.gz", "cell-contexts.csv.gz"]:
        if digest(source / name) != manifest["outputs"][name]:
            raise ValueError(f"Pinned evidence hash mismatch: {name}")
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(source / "gene-comparison.csv.gz", low_memory=False)
    bulk = pd.read_csv(source / "bulk-contexts.csv.gz", low_memory=False)
    cells = pd.read_csv(
        source / "cell-contexts.csv.gz",
        usecols=[
            "gene_id",
            "source",
            "version",
            "taxonomy_id",
            "tissue",
            "cell_type",
            "state",
            "stage",
            "unit",
            "observations",
            "supported_observations",
            "raw_maximum",
            "value",
            "mapping_status",
            "concept_id",
            "category",
            "neutral_reproductive_immune",
            "ihc_exact_match",
            "ihc_positive",
            "negative_ihc_grade",
        ],
        low_memory=False,
    )
    union = eligibility(genes, "recurrence") | eligibility(genes, "maximum10")
    if not union.equals(genes.ihc_contexts_evaluated):
        raise ValueError("Evidence pool must cover every potentially selected gene")
    if set(cells.gene_id) != set(genes.loc[union, "gene_id"]) or set(bulk.Ensembl_Gene_ID) != set(
        cells.gene_id
    ):
        raise ValueError("Incomplete local evidence universe")
    scenarios, cascades, mode_gates, assessments = [], [], {}, {}
    for mode in MODES:
        b = policy.apply_bulk_ladder(bulk, mode=mode)
        c = cell_assessment(cells, mode)
        gates = normal_passes(genes, b, c)
        mode_gates[mode] = gates
        for tumor in ["recurrence", "maximum10", "maximum30"]:
            chosen = eligibility(genes, tumor) & gates.all(axis=1)
            genes[f"{mode}_{tumor}_selected"] = chosen
            scenarios.append({"mode": mode, "tumor": tumor, "selected": int(chosen.sum())})
        cascade, failures = funnel(genes, gates, mode=mode)
        genes[f"{mode}_first_failure"] = failures
        cascades.append(cascade)
        if mode != "ihc_doubling":
            assessments[mode] = (b, c)
        print(mode, scenarios[-3:], flush=True)
    for tumor, old_col in [
        ("recurrence", "linked_recurrence_selected"),
        ("maximum10", "linked_maximum_tumor_selected"),
        ("maximum30", "linked_tumor30_selected"),
    ]:
        if not genes[f"legacy_{tumor}_selected"].equals(genes[old_col]):
            raise ValueError("Independent reapplication did not replay the earlier screen")
    selected = genes.two_axis_doubling_maximum10_selected
    old_selected = genes.legacy_maximum10_selected
    genes["membership_change"] = np.select(
        [selected & ~old_selected, ~selected & old_selected, selected & old_selected],
        ["Gained", "Lost", "Retained"],
        default="Excluded in both",
    )
    new_bulk, new_cells = assessments["two_axis_doubling"]
    old_bulk, old_cells = assessments["legacy"]
    changes = pd.concat(
        [
            changed_contexts(old_bulk, new_bulk, kind="bulk"),
            changed_contexts(old_cells, new_cells, kind="cell"),
        ],
        ignore_index=True,
    )
    changes["symbol"] = changes.gene_id.map(genes.set_index("gene_id").symbol)
    changes["membership_change"] = changes.gene_id.map(genes.set_index("gene_id").membership_change)
    write(changes, out / "changed-contexts.csv.gz")
    write(
        changes[changes.membership_change.isin(["Gained", "Lost"])],
        out / "membership-change-reasons.csv",
    )
    write(pd.concat(cascades, ignore_index=True), out / "filter-cascade.csv")
    write(pd.DataFrame(scenarios), out / "scenario-counts.csv")
    write(policy.linked_thresholds(mode="two_axis_doubling"), out / "cell-thresholds.csv")
    write(new_bulk, out / "bulk-contexts.csv.gz")
    write(new_cells, out / "cell-contexts.csv.gz")
    # Failure counts in a common denominator, plus leave-one-filter-out effects.
    pre = eligibility(genes)
    gates = mode_gates["two_axis_doubling"]
    independent = []
    for label in gates:
        independent.append(
            {
                "filter": label,
                "eligible_before_normal_rna": int(pre.sum()),
                "independent_failures": int((pre & ~gates[label]).sum()),
                "uniquely_excluded": int(
                    (pre & ~gates[label] & gates.drop(columns=label).all(axis=1)).sum()
                ),
            }
        )
    write(pd.DataFrame(independent), out / "filter-independent-effects.csv")
    for label in gates:
        genes[f"pass_{label}"] = gates[label]
    # All numeric IHC steps are evaluated with identical tissue and source facts.
    ihc_rows = []
    for grade in ["Missing", *policy.GRADES]:
        b = policy.apply_bulk_ladder(bulk, mode="two_axis_doubling", max_ihc_grade=grade)
        c = cell_assessment(
            cells,
            "two_axis_doubling",
            max_ihc_grade="Missing/uncertain" if grade == "Missing" else grade,
        )
        chosen = pre & normal_passes(genes, b, c).all(axis=1)
        genes[f"ihc_cap_{grade}_selected"] = chosen
        ihc_rows.append({"maximum_ihc_allowance": grade, "selected": int(chosen.sum())})
    write(pd.DataFrame(ihc_rows), out / "ihc-choice-counts.csv")
    write(genes, out / "gene-comparison.csv.gz")
    write(genes[selected], out / "selected-genes.csv")
    write(genes[genes.membership_change.isin(["Gained", "Lost"])], out / "membership-changes.csv")
    write(genes[genes.symbol.isin(NAMED)], out / "named-gene-comparison.csv")
    # Normalize only by each observation's own policy ceiling; keep source units.
    margin_rows = []
    for label, frame, id_col, numerator, denominator in [
        ("Normal bulk", new_bulk, "Ensembl_Gene_ID", "rna_ntpm", "effective_ceiling_ntpm"),
        *[
            (
                label,
                new_cells[new_cells.assay.eq(assay) & new_cells.category.isin(policy.NORMAL)],
                "gene_id",
                "value",
                "cta_cutoff",
            )
            for assay, label in ASSAYS.items()
        ],
    ]:
        f = frame.assign(ratio=frame[numerator] / frame[denominator])
        best = f[f.ratio.notna()].sort_values("ratio", ascending=False).drop_duplicates(id_col)
        for row in best.to_dict("records"):
            margin_rows.append(
                {
                    "gene_id": row[id_col],
                    "filter": label,
                    "ratio": row["ratio"],
                    "tissue": row["tissue"],
                    "cell_type": row.get("cell_type", ""),
                    "value": row[numerator],
                    "ceiling": row[denominator],
                    "unit": row.get("unit", "nTPM"),
                }
            )
    margins = pd.DataFrame(margin_rows)
    margins["symbol"] = margins.gene_id.map(genes.set_index("gene_id").symbol)
    write(margins, out / "normal-rna-margins.csv.gz")
    write(margins[margins.symbol.isin(NAMED)], out / "named-normal-rna-margins.csv")
    summary = {
        "policy": policy.GEOMETRIC_POLICY,
        "candidate_formula": "2 * 2**tissue_rank * 2**ihc_rank",
        "reference_review_multiplier": policy.REFERENCE_REVIEW_MULTIPLIER,
        "universe": len(genes),
        "evaluated_evidence_pool": int(union.sum()),
        "selected": int(selected.sum()),
        "gained": int((selected & ~old_selected).sum()),
        "lost": int((~selected & old_selected).sum()),
        "retained": int((selected & old_selected).sum()),
        "historical_labels": genes.loc[selected, "Expression category"].value_counts().to_dict(),
        "named_selected": genes.loc[selected & genes.symbol.isin(NAMED), "symbol"].tolist(),
        "scenarios": scenarios,
        "ihc_choices": ihc_rows,
        "scope": "All candidate associations and source measurements are unchanged; the entire prior eligible union is reassessed, not just previously selected genes.",
        "cascade_order": "Association, tumor, IHC detection, RNA coverage, bulk, HPA tissue clusters, HPA immune, Monaco, DICE, adult cornea. Sequential removals depend on order; independent failures overlap.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    inputs = {
        str((source / name).relative_to(ROOT)): digest(source / name)
        for name in [
            "run-manifest.json",
            "gene-comparison.csv.gz",
            "bulk-contexts.csv.gz",
            "cell-contexts.csv.gz",
        ]
    }
    for name in [
        "scripts/cta_geometric_threshold_report.py",
        "oncoref/cell_rna_policy.py",
        "oncoref/data/normal-tissue-taxonomy.yaml",
    ]:
        inputs[name] = digest(ROOT / name)
    outputs = {p.name: digest(p) for p in out.iterdir() if p.is_file()}
    (out / "run-manifest.json").write_text(
        json.dumps({"inputs": inputs, "sources": manifest["sources"], "outputs": outputs}, indent=2)
        + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=ROOT / "docs/audits/cta-linked-thresholds-20261027"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.source.resolve(), args.output)
