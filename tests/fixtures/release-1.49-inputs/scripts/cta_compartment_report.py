#!/usr/bin/env python3
"""Replay the versioned human CTA research policy across the full reference."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_cell_policy_report import summarize  # noqa: E402
from cta_expanded_funnel import markdown_table  # noqa: E402
from cta_expanded_normal import bulk_assessment  # noqa: E402
from cta_harmonized_report import decision_layers, figures  # noqa: E402
from cta_quantitative_threshold_report import positive_profile_controls  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import (  # noqa: E402
    cell_rna_policy,
    hpa,
    normal_cell_data,
    normal_cell_screen,
    reference_data,
)
from oncoref.compartment_rna_policy import (  # noqa: E402
    IHC_RANKS,
    POLICY_ID,
    TIER_BASES,
    assess_contexts,
    assess_imported_protein_reviews,
    collapse_hpa_ihc,
)
from oncoref.cta_funnel import classify_layers  # noqa: E402

BASE = ROOT / "docs/audits/cta-harmonized-20261027"
EXPANDED = ROOT / "docs/audits/cta-expanded-funnel-20261027"
PRIOR = ROOT / "docs/audits/cta-compartment-policy-20261027"
OUTPUT = ROOT / "docs/audits/cta-release-1.21.0"
PROBLEMATIC = {
    "CEACAM5",
    "CA9",
    "ERBB2",
    "MAGEA3",
    "MAGEA12",
    "TTN",
    "FAP",
    "IL5",
    "IL13",
    "ABCB5",
    "CTAG2",
}


def aggregate_contexts(frame: pd.DataFrame) -> pd.DataFrame:
    """Retain independently contributing holds and strong evidence per gene."""
    return frame.groupby("gene_id")[["candidate_exclusion", "unsafe_evidence"]].any()


def label_baseline_locus_decisions(frame: pd.DataFrame) -> pd.DataFrame:
    """Prevent older locus masks from masquerading as the current group policy."""
    names = {c: "baseline_" + c for c in frame if c.startswith(("hold:", "unsafe:"))}
    names.update(
        {
            "gene_category": "baseline_gene_category",
            "gene_exclusion_reasons": "baseline_gene_exclusion_reasons",
            "Unsafe_reasons": "baseline_gene_unsafe_reasons",
            "policy": "baseline_policy",
            "Expression category": "historical_expression_category",
            "Expression detail": "historical_expression_detail",
            "Expression forbidden": "historical_expression_forbidden",
        }
    )
    return frame.rename(columns=names)


def apply_hpa_protein_levels(facts: pd.DataFrame, output: Path) -> set[str]:
    """Treat Low as an observed flag without erasing stronger HPA protein findings."""
    d = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    d = d.rename(columns={"Gene": "gene_id", "Tissue": "tissue", "Cell type": "cell_type"})
    d["source"], d["version"], d["taxonomy_id"] = "hpa_normal_tissue", "v23", 9606
    d = normal_cell_screen.contexts(d)
    normal = d.category.isin(["vital", "protected", "manageable"])
    d["low_protein"] = normal & d.Level.eq("Low")
    d["reliable_low_protein"] = d.low_protein & d.Reliability.isin(cell_rna_policy.GRADES)
    d["protein_exclusion"] = normal & d.Level.isin(cell_rna_policy.POSITIVE - {"Low"})
    d["reliable_protein_exclusion"] = d.protein_exclusion & d.Reliability.isin(
        cell_rna_policy.GRADES
    )
    summary = d.groupby("gene_id")[
        ["low_protein", "reliable_low_protein", "protein_exclusion", "reliable_protein_exclusion"]
    ].any()
    write_csv(summary.reset_index(), output / "hpa-protein-level-summary.csv.gz")
    members = facts.member_gene_ids.map(json.loads)
    low_ids = set(summary.index[summary.low_protein])
    other_ids = set(summary.index[summary.protein_exclusion])
    only_low = members.map(
        lambda group: bool(set(group) & low_ids) and not bool(set(group) & other_ids)
    )
    reliable_low_ids = set(summary.index[summary.reliable_low_protein])
    reliable_other_ids = set(summary.index[summary.reliable_protein_exclusion])
    only_reliable_low = members.map(
        lambda group: (
            bool(set(group) & reliable_low_ids) and not bool(set(group) & reliable_other_ids)
        )
    )
    # Relax only a baseline HPA finding whose observed positives are all Low.
    # Unexplained prior findings cannot disappear merely because this extraction
    # lacks a matching row or uses a different source annotation.
    facts["hpa_low_only_relaxation"] = facts.ihc_detected & only_low
    facts["ihc_detected"] &= ~only_low
    facts["hpa_reliable_low_relaxation"] = facts.ihc_reliable_detected & only_reliable_low
    facts["ihc_reliable_detected"] &= ~only_reliable_low
    return low_ids


def group_ihc_evidence(
    facts: pd.DataFrame, output: Path
) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Rebuild group IHC levels and require positive controls for every member."""
    groups = facts[facts.member_count.gt(1)]
    raw = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    combined = collapse_hpa_ihc(raw, groups)
    controls = positive_profile_controls(raw)
    group_controls = {
        row.representative_gene_id
        for row in groups.itertuples()
        if set(json.loads(row.member_gene_ids)) <= controls
    }
    conflicts = set(groups.loc[groups["Antibody concern"], "representative_gene_id"])
    write_csv(combined, output / "group-ihc.csv.gz")
    return combined, group_controls, conflicts


def apply_assay(
    facts: pd.DataFrame, assessment: pd.DataFrame, assay: str, *, subset: pd.Series | None = None
) -> None:
    """Assign one recomputed assay, preserving absence of required source rows."""
    selected = pd.Series(True, index=facts.index) if subset is None else subset
    for prefix, column in [("hold:", "candidate_exclusion"), ("unsafe:", "unsafe_evidence")]:
        observed = facts.representative_gene_id.map(assessment[column])
        use = selected & observed.notna()
        facts.loc[use, prefix + assay] = observed[use].astype(bool)


def cluster_assessments(
    facts: pd.DataFrame,
    output: Path,
    grouped_ihc: tuple[pd.DataFrame, set[str], set[str]],
    *,
    reuse: bool = False,
) -> tuple[pd.DataFrame, list[dict]]:
    """Stream every singleton and use already summed multi-locus cluster values."""
    receipts = [
        reference_data.provenance(source, verify_content=True)
        for source in ["hpa_single_cell_tissue", "hpa_normal_tissue"]
    ]
    sources = [{k: r[k] for k in ["name", "version", "url", "sha256"]} for r in receipts]
    cache_inputs = input_hashes(
        [
            ROOT / "oncoref/compartment_rna_policy.py",
            ROOT / "oncoref/normal_cell_screen.py",
            ROOT / "oncoref/cell_rna_policy.py",
            ROOT / "oncoref/normal_tissues.py",
            ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
            ROOT / "oncoref/data/normal-tissue-mappings.csv",
            ROOT / "scripts/cta_quantitative_threshold_report.py",
            EXPANDED / "gene-classifications.csv.gz",
            BASE / "protein-group-classifications.csv.gz",
            BASE / "group-normal/cluster-contexts.csv.gz",
            output / "group-ihc.csv.gz",
        ]
    )
    signature = hashlib.sha256(
        (
            inspect.getsource(cluster_assessments)
            + inspect.getsource(aggregate_contexts)
            + inspect.getsource(group_ihc_evidence)
        ).encode()
    ).hexdigest()
    cache_path = output / "cluster-replay.json"
    if reuse and cache_path.exists():
        cache = json.loads(cache_path.read_text())
        if (
            cache["inputs"] != cache_inputs
            or cache["sources"] != sources
            or cache["function_sha256"] != signature
        ):
            raise ValueError("Cluster cache inputs changed; run without --reuse-clusters")
        for name, digest in cache["outputs"].items():
            if sha256(output / name) != digest:
                raise ValueError("Cluster cache output changed")
        return pd.read_csv(output / "cluster-assessments.csv.gz").set_index("gene_id"), sources
    singleton_ids = set(facts.loc[facts.member_count.eq(1), "representative_gene_id"])
    old = pd.read_csv(EXPANDED / "gene-classifications.csv.gz", low_memory=False)
    conflicts = set(old.loc[old["Antibody concern"], "gene_id"])
    ihc = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    controls = positive_profile_controls(ihc)
    receipt = receipts[0]
    summaries, changed, peaks = [], [], []

    def collect(d: pd.DataFrame) -> None:
        result = assess_contexts(d)
        summaries.append(aggregate_contexts(result))
        measured = result[result.candidate_ceiling.notna() & result.value.notna()].copy()
        measured["ratio"] = measured.value / measured.candidate_ceiling
        peaks.append(measured.sort_values("ratio", ascending=False).drop_duplicates("gene_id"))
        altered = result.compartment_basis.ne("anatomical_default") & result.value.notna()
        # Preserve every measured context whose compartment rule changed.
        changed.append(
            result.loc[
                altered,
                [
                    "gene_id",
                    "source",
                    "tissue",
                    "cell_type",
                    "unit",
                    "value",
                    "category",
                    "selection_tier",
                    "compartment_basis",
                    "stem_progenitor_caution",
                    "negative_ihc_grade",
                    "candidate_ceiling",
                    "strong_boundary",
                    "candidate_exclusion",
                    "unsafe_evidence",
                ],
            ]
        )

    for number, chunk in enumerate(pd.read_csv(receipt["path"], sep="\t", chunksize=250000)):
        d = (
            chunk[chunk.Gene.isin(singleton_ids)]
            .rename(
                columns={
                    "Gene": "gene_id",
                    "Tissue": "tissue",
                    "Cell type": "cell_type",
                    "nTPM": "value",
                }
            )
            .copy()
        )
        if d.empty:
            continue
        d["source"], d["version"], d["taxonomy_id"] = "hpa_single_cell_tissue", "v23", 9606
        d["state"], d["stage"], d["unit"] = "unspecified", "adult", "nTPM"
        d["value"] = d.value.where(d["Read count"].ge(20))
        d = normal_cell_screen.contexts(d)
        d = cell_rna_policy.reproductive_immune_contexts(d)
        d = cell_rna_policy.protected_immune_contexts(d)
        d = cell_rna_policy.match_negative_ihc(
            d,
            ihc[ihc.Gene.isin(d.gene_id.unique())],
            positive_controls=controls,
            conflicting_genes=conflicts,
        )
        collect(d)
        if number % 10 == 0:
            print(f"Cluster source chunk {number}", flush=True)
    grouped = pd.read_csv(BASE / "group-normal/cluster-contexts.csv.gz", low_memory=False)
    grouped = grouped.drop(
        columns=["Level", "Reliability", "ihc_positive", "ihc_exact_match", "negative_ihc_grade"],
        errors="ignore",
    )
    combined, group_controls, group_conflicts = grouped_ihc
    grouped = cell_rna_policy.match_negative_ihc(
        grouped, combined, positive_controls=group_controls, conflicting_genes=group_conflicts
    )
    collect(grouped)
    write_csv(pd.concat(changed, ignore_index=True), output / "changed-compartment-contexts.csv.gz")
    aggregated = pd.concat(summaries).groupby(level=0).any()
    write_csv(aggregated.reset_index(), output / "cluster-assessments.csv.gz")
    write_csv(
        pd.concat(peaks, ignore_index=True)
        .sort_values("ratio", ascending=False)
        .drop_duplicates("gene_id"),
        output / "cluster-peak-evidence.csv.gz",
    )
    cache = {
        "inputs": cache_inputs,
        "sources": sources,
        "function_sha256": signature,
        "outputs": {
            name: sha256(output / name)
            for name in [
                "cluster-assessments.csv.gz",
                "changed-compartment-contexts.csv.gz",
                "cluster-peak-evidence.csv.gz",
            ]
        },
    }
    cache_path.write_text(json.dumps(cache, indent=2) + "\n")
    return aggregated, sources


def bulk_assessments(
    facts: pd.DataFrame, output: Path, grouped_ihc: tuple[pd.DataFrame, set[str], set[str]]
) -> pd.DataFrame:
    """Reassess every matched bulk context, preserving missing-measurement holds."""
    single = pd.read_csv(EXPANDED / "normal/bulk-contexts.csv.gz", low_memory=False)
    grouped = pd.read_csv(BASE / "group-normal/bulk-contexts.csv.gz", low_memory=False)
    combined, controls, conflicts = grouped_ihc
    grouped = bulk_assessment(
        grouped[["Ensembl_Gene_ID", "tissue", "rna_ntpm"]].rename(
            columns={"Ensembl_Gene_ID": "Gene", "tissue": "Tissue", "rna_ntpm": "nTPM"}
        ),
        combined,
        controls,
        conflicts,
    )
    d = pd.concat(
        [single[~single.Ensembl_Gene_ID.isin(grouped.Ensembl_Gene_ID)], grouped], ignore_index=True
    )
    d = d.rename(
        columns={"Ensembl_Gene_ID": "gene_id", "rna_ntpm": "value", "screening_tier": "category"}
    )
    d["source"], d["taxonomy_id"], d["stage"], d["unit"], d["cell_type"] = (
        "hpa_rna_consensus",
        9606,
        "adult",
        "nTPM",
        "bulk",
    )
    d["mapping_status"] = np.where(d.category.eq("unresolved"), "unmapped", "mapped")
    r = assess_contexts(d)
    r["candidate_exclusion"] |= r.category.isin(TIER_BASES) & r.value.isna()
    write_csv(
        r[
            [
                "gene_id",
                "tissue",
                "category",
                "value",
                "unit",
                "negative_ihc_grade",
                "candidate_ceiling",
                "strong_boundary",
                "candidate_exclusion",
                "unsafe_evidence",
            ]
        ],
        output / "bulk-contexts.csv.gz",
    )
    return aggregate_contexts(r)


def other_assay_assessments(facts: pd.DataFrame, output: Path) -> list[dict]:
    """Replay complete immune/corneal contexts instead of rescaling a peak.

    A peak RNA ratio cannot encode a second, unmapped context or an independent
    protein finding. Reassessing all contexts preserves those concurrent holds.
    Multi-locus DICE and cornea use their dedicated uncertainty bounds afterward.
    """
    singleton_ids = set(facts.loc[facts.member_count.eq(1), "representative_gene_id"])
    records: list[pd.DataFrame] = []
    receipts = []
    ambiguous = []

    def receipt(source: str) -> None:
        item = reference_data.provenance(source, verify_content=True)
        receipts.append({k: item[k] for k in ["name", "version", "url", "sha256"]})

    def collect(data: pd.DataFrame, assay: str) -> None:
        if data.empty:
            return
        d = normal_cell_screen.contexts(
            data.drop(
                columns=[
                    "category",
                    "mapping_status",
                    "concept_id",
                    "prior_context_category",
                    "prior_immune_category",
                ],
                errors="ignore",
            )
        )
        d = cell_rna_policy.reproductive_immune_contexts(d)
        d = cell_rna_policy.protected_immune_contexts(d)
        # No reviewed cell/IHC equivalences exist for these source assays.
        d["negative_ihc_grade"] = "Missing/uncertain"
        d["ihc_positive"] = False
        d = assess_contexts(d)
        d["assay"] = assay
        records.append(d)

    for source, label in [
        ("hpa_immune_cell", "HPA immune"),
        ("hpa_immune_monaco", "Monaco immune"),
    ]:
        receipt(source)
        d = normal_cell_screen.contexts(normal_cell_data.observations(source, genes=singleton_ids))
        collect(summarize(d, support=d.value.notna()), label)
        collect(pd.read_csv(BASE / f"group-normal/{source}-contexts.csv.gz"), label)
    for source in normal_cell_data.sources().source:
        if not source.startswith("dice_"):
            continue
        receipt(source)
        d = normal_cell_data.dice_observations(source, genes=singleton_ids)
        duplicated = set(d.loc[d.duplicated(["gene_id", "replicate_slot"], keep=False), "gene_id"])
        ambiguous.extend(
            {
                "gene_id": gid,
                "source": source,
                "reason": "ambiguous replicate identity; no invented donor pairing",
            }
            for gid in sorted(duplicated)
        )
        d = normal_cell_screen.contexts(d[~d.gene_id.isin(duplicated)])
        collect(
            summarize(d, support=d.value.notna(), rank=3, unique="replicate_slot"), "DICE immune"
        )
    for source in ["cornea_collin_adult_counts", "cornea_collin_adultcornea_cells"]:
        receipt(source)
    d = normal_cell_screen.contexts(normal_cell_data.cornea_summary(genes=singleton_ids))
    collect(
        summarize(d, support=d.positive_cells.ge(5) & d.read_count.ge(20), rank=2, unique="sample"),
        "Adult cornea",
    )
    contexts = pd.concat(records, ignore_index=True)
    write_csv(contexts, output / "immune-cornea-contexts.csv.gz")
    write_csv(
        pd.DataFrame(ambiguous, columns=["gene_id", "source", "reason"]),
        output / "optional-ambiguous-observations.csv",
    )
    for assay, part in contexts.groupby("assay"):
        apply_assay(facts, aggregate_contexts(part), assay)
    return receipts


def uncertain_group_assays(facts: pd.DataFrame, output: Path) -> None:
    """Recompute unpaired DICE and corneal detection uncertainty at new limits."""
    multi = facts.member_count.gt(1)
    dice = pd.read_csv(BASE / "group-normal/dice-rank-bounds.csv.gz")
    d = assess_contexts(dice)
    d["candidate_exclusion"] |= d.measured_members.eq(d.member_count) & d.value_upper.gt(
        d.candidate_ceiling
    )
    apply_assay(facts, aggregate_contexts(d), "DICE immune", subset=multi)
    write_csv(d, output / "dice-group-contexts.csv.gz")
    lower = pd.read_csv(BASE / "group-normal/cornea-contexts.csv.gz")
    d = assess_contexts(lower)
    raw = pd.read_csv(BASE / "group-normal/cornea-detection-bounds.csv.gz", low_memory=False)
    upper = summarize(
        normal_cell_screen.contexts(raw),
        support=raw.positive_cells_upper.ge(5) & raw.read_count.ge(20),
        rank=2,
        unique="sample",
    )
    keys = ["gene_id", "source", "version", "taxonomy_id", "tissue", "cell_type", "state", "stage"]
    d = d.merge(
        upper[[*keys, "value"]].rename(columns={"value": "upper_supported_value"}),
        on=keys,
        how="left",
        validate="one_to_one",
    )
    d["candidate_exclusion"] |= d.upper_supported_value.gt(d.candidate_ceiling)
    apply_assay(facts, aggregate_contexts(d), "Adult cornea", subset=multi)
    write_csv(d, output / "cornea-group-contexts.csv.gz")


def threshold_figure(output: Path) -> None:
    """Replace the reused report's old threshold figure with this policy's grid."""
    import matplotlib.pyplot as plt
    import plot_cta_geometric as style

    grid = np.array([[v * 2**j for j in range(4)] for v in TIER_BASES.values()])
    fig, ax = plt.subplots(figsize=(12, 7.5))
    fig.subplots_adjust(left=0.23, right=0.96, top=0.78, bottom=0.24)
    style.title(
        fig,
        "Provisional compartment-aware RNA ceilings",
        "Research prioritization · values are policy choices, not validated toxicity thresholds",
    )
    ax.imshow(np.log2(grid), cmap="YlGnBu", aspect="auto")
    ax.set_xticks(
        range(4),
        ["Missing / unusable", "Approved negative", "Supported negative", "Enhanced negative"],
    )
    ax.set_yticks(
        range(4), ["Vital", "Protected", "Manageable / fibroblasts", "Specified bowel epithelium"]
    )
    for i, j in np.ndindex(grid.shape):
        ax.text(
            j,
            i,
            f"{grid[i, j]:g}\nStrong ≥{2 * grid[i, j]:g}",
            ha="center",
            va="center",
            color="white" if grid[i, j] >= 32 else style.INK,
            fontsize=14,
        )
    style.footer(
        fig,
        "Same tier and IHC grade use linked native TPM/nTPM numbers; CP10K is nominal /100.\nOnly named adult HPA cell compartments change tiers. Bulk keeps its organ tier.\nUndifferentiated bowel cells retain a stem/progenitor caution; Medium/High HPA protein and clinical exclusions remain active.",
    )
    style.save(fig, output / "figures", "threshold-grid")


def component_effects(facts: pd.DataFrame, contexts: pd.DataFrame) -> pd.DataFrame:
    """Count cumulative policy changes, restoring stricter masks for ablations.

    Counts are order-dependent contributions, not independent causal effects.
    Restoring the anatomical threshold can only add exclusions to the fully
    replayed cluster result; all unchanged cell contexts remain represented.
    """
    baseline = pd.read_csv(BASE / "protein-group-classifications.csv.gz", low_memory=False)
    rows = []
    steps = [
        ("Previous policy", None, False, False),
        ("Double linked RNA grid", set(), False, False),
        ("Also classify fibroblasts as Manageable", {"identified_fibroblast"}, False, False),
        (
            "Also allow specified large-bowel epithelium",
            {"identified_fibroblast", "large_bowel_epithelial_population"},
            False,
            False,
        ),
        ("Also flag Low HPA protein", None, True, False),
        ("Also flag explicitly Low imported protein", None, True, True),
    ]
    for label, enabled, low_hpa, low_imported in steps:
        if label == "Previous policy":
            result = baseline
        else:
            variant = facts.copy()
            if not low_hpa:
                for field in ["ihc_detected", "ihc_reliable_detected"]:
                    variant[field] = baseline[field]
            if not low_imported:
                variant["Accepted normal protein"] = baseline["Accepted normal protein"]
            if enabled is not None:
                restore = contexts[~contexts.compartment_basis.isin(enabled)].copy()
                ceiling = (
                    restore.category.map(TIER_BASES)
                    * 2.0 ** restore.negative_ihc_grade.map(IHC_RANKS)
                    * restore.unit.map({"TPM": 1.0, "nTPM": 1.0, "CP10K": 0.01})
                )
                restore["candidate_exclusion"] = restore.value.gt(ceiling)
                restore["unsafe_evidence"] = restore.value.ge(2 * ceiling)
                extra = aggregate_contexts(restore)
                for prefix, column in [
                    ("hold:", "candidate_exclusion"),
                    ("unsafe:", "unsafe_evidence"),
                ]:
                    variant[prefix + "HPA tissue clusters"] |= (
                        variant.representative_gene_id.map(extra[column]).fillna(False).astype(bool)
                    )
            result = classify_layers(
                variant, decision_layers(variant), entity_column="protein_equivalence_group_id"
            ).classifications
        counts = {
            name: int(result.category.eq(name).sum()) for name in ["CTA", "Unsafe", "Neither"]
        }
        rows.append({"cumulative_policy": label, **counts})
    frame = pd.DataFrame(rows)
    frame["additional_CTA"] = frame.CTA.diff().fillna(0).astype(int)
    return frame


def review_changes(result: pd.DataFrame, output: Path) -> dict:
    """Compare this review with the actual prior proposal, including layer masks."""
    prior = pd.read_csv(PRIOR / "protein-group-classifications.csv.gz", low_memory=False).set_index(
        "protein_equivalence_group_id"
    )
    current = result.set_index("protein_equivalence_group_id").reindex(prior.index)
    changes = []
    for field in [
        "category",
        *[c for c in current if c.startswith(("hold:", "unsafe:"))],
        "stem_progenitor_expression_flag",
    ]:
        mask = current[field].ne(prior[field])
        changes.extend(
            {
                "protein_equivalence_group_id": gid,
                "label": current.loc[gid, "protein_group_label"],
                "field": field,
                "previous": prior.loc[gid, field],
                "reviewed": current.loc[gid, field],
            }
            for gid in current.index[mask]
        )
    write_csv(
        pd.DataFrame(
            changes,
            columns=["protein_equivalence_group_id", "label", "field", "previous", "reviewed"],
        ),
        output / "review-decision-changes.csv",
    )
    old = pd.read_csv(BASE / "group-normal/group-ihc.csv.gz")
    new = pd.read_csv(output / "group-ihc.csv.gz")
    ihc = old.merge(
        new,
        on=["Gene", "Tissue", "Cell type"],
        suffixes=("_previous", "_reviewed"),
        validate="one_to_one",
    )
    changed = ihc.Level_previous.ne(ihc.Level_reviewed)
    write_csv(ihc[changed], output / "review-ihc-level-changes.csv.gz")
    return {
        "changed_final_categories": int(current.category.ne(prior.category).sum()),
        "changed_layer_decisions": len(changes),
        "corrected_group_ihc_levels": int(changed.sum()),
    }


def named_normal_evidence(named: pd.DataFrame, output: Path) -> pd.DataFrame:
    """Report the maximum measured fraction of each assay's own policy ceiling."""
    records = []
    for filename, assay in [
        ("bulk-contexts.csv.gz", "Normal bulk"),
        ("cluster-peak-evidence.csv.gz", "HPA tissue clusters"),
        ("immune-cornea-contexts.csv.gz", None),
        ("dice-group-contexts.csv.gz", "DICE immune"),
        ("cornea-group-contexts.csv.gz", "Adult cornea"),
    ]:
        frame = pd.read_csv(output / filename, low_memory=False)
        frame = frame[
            frame.gene_id.isin(named.representative_gene_id)
            & frame.value.notna()
            & frame.candidate_ceiling.notna()
        ].copy()
        if assay is not None:
            frame["assay"] = assay
        frame["ratio"] = frame.value / frame.candidate_ceiling
        frame["cell_type"] = frame.get("cell_type", "bulk")
        records.append(frame)
    combined = (
        pd.concat(records, ignore_index=True)
        .sort_values("ratio", ascending=False)
        .drop_duplicates(["gene_id", "assay"])
    )
    combined = combined.drop(columns=["protein_group_label"], errors="ignore").merge(
        named[["representative_gene_id", "protein_group_label", "category"]],
        left_on="gene_id",
        right_on="representative_gene_id",
        suffixes=("_tissue", "_decision"),
        validate="many_to_one",
    )
    columns = [
        "protein_group_label",
        "category_decision",
        "assay",
        "tissue",
        "cell_type",
        "value",
        "unit",
        "negative_ihc_grade",
        "candidate_ceiling",
        "strong_boundary",
        "ratio",
    ]
    write_csv(combined[columns], output / "named-normal-evidence.csv")
    return combined[columns]


def build(output: Path = OUTPUT, *, reuse_clusters: bool = False) -> None:
    """Publish complete decisions, admission changes, known targets and figures."""
    output.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(BASE / "protein-group-classifications.csv.gz", low_memory=False)
    facts["baseline_category"] = facts.category
    facts["baseline_exclusion_reasons"] = facts.exclusion_reasons
    facts = facts.rename(columns={"member_gene_categories": "baseline_member_gene_categories"})
    facts["baseline_accepted_normal_protein"] = facts["Accepted normal protein"]
    grouped_ihc = group_ihc_evidence(facts, output)
    apply_assay(facts, bulk_assessments(facts, output, grouped_ihc), "Normal bulk")
    clusters, receipts = cluster_assessments(facts, output, grouped_ihc, reuse=reuse_clusters)
    apply_assay(facts, clusters, "HPA tissue clusters")
    receipts.extend(other_assay_assessments(facts, output))
    uncertain_group_assays(facts, output)
    low_hpa = apply_hpa_protein_levels(facts, output)
    protein = assess_imported_protein_reviews(
        pd.read_csv(ROOT / "oncoref/data/cta-policy-findings.csv"),
        pd.read_csv(output / "protein-intensity-input.csv"),
    )
    write_csv(protein, output / "imported-protein-review-decisions.csv")
    retained = set(protein.loc[protein.retained_protein_exclusion, "Gene ID"])
    limited = set(protein.loc[protein.limited_review_flag, "Gene ID"])
    members = facts.member_gene_ids.map(json.loads)
    facts["Accepted normal protein"] = members.map(lambda ids: bool(set(ids) & retained))
    facts["limited_protein_review_flag"] = members.map(lambda ids: bool(set(ids) & limited))
    low_ids = low_hpa | set(protein.loc[protein.low_normal_protein_flag, "Gene ID"])
    facts["low_normal_protein_flag"] = members.map(lambda ids: bool(set(ids) & low_ids))
    result, funnel, independent = classify_layers(
        facts, decision_layers(facts), entity_column="protein_equivalence_group_id"
    )
    reverse = classify_layers(
        facts, list(reversed(decision_layers(facts))), entity_column="protein_equivalence_group_id"
    ).classifications
    if not result.category.equals(reverse.category):
        raise ValueError("Classification changed with layer order")
    result["policy"] = POLICY_ID
    cautions = pd.read_csv(output / "changed-compartment-contexts.csv.gz")
    components = component_effects(facts, cautions)
    write_csv(components, output / "policy-component-effects.csv")
    stem_ids = set(cautions.loc[cautions.stem_progenitor_caution, "gene_id"])
    result["stem_progenitor_expression_flag"] = result.representative_gene_id.isin(stem_ids)
    reviewed = review_changes(result, output)
    write_csv(result, output / "protein-group-classifications.csv.gz")
    for category in ["CTA", "Unsafe", "Neither"]:
        write_csv(
            result[result.category.eq(category)],
            output / f"{category.lower()}-protein-groups.csv.gz",
        )
    write_csv(funnel, output / "filtering-funnel.csv")
    write_csv(independent, output / "independent-layer-effects.csv")
    changed = result[result.category.ne(result.baseline_category)]
    added = changed[changed.category.eq("CTA")]
    write_csv(changed, output / "classification-changes.csv")
    write_csv(added, output / "newly-selected-groups.csv")
    named = result[result.member_named_controls.fillna("").ne("")]
    write_csv(named, output / "named-antigen-decisions.csv")
    named_normal_evidence(named, output)
    write_csv(
        result[result.protein_group_label.isin(PROBLEMATIC)],
        output / "problematic-target-checks.csv",
    )
    known = pd.read_csv(BASE / "known-cta-decisions.csv.gz", low_memory=False)
    known = label_baseline_locus_decisions(known)
    lookup = result.set_index("protein_equivalence_group_id")
    known["baseline_protein_group_category"] = known.protein_group_category
    known["protein_group_category"] = known.protein_equivalence_group_id.map(lookup.category)
    known["policy"] = POLICY_ID
    known["protein_group_exclusion_reasons"] = known.protein_equivalence_group_id.map(
        lookup.exclusion_reasons
    )
    write_csv(known, output / "known-cta-decisions.csv.gz")
    crosswalk = pd.read_csv(BASE / "gene-to-protein-group-decisions.csv.gz", low_memory=False)
    crosswalk = label_baseline_locus_decisions(crosswalk)
    crosswalk["baseline_protein_group_category"] = crosswalk.protein_group_category
    crosswalk["protein_group_category"] = crosswalk.protein_equivalence_group_id.map(
        lookup.category
    )
    crosswalk["protein_group_exclusion_reasons"] = crosswalk.protein_equivalence_group_id.map(
        lookup.exclusion_reasons
    )
    crosswalk["policy"] = POLICY_ID
    write_csv(crosswalk, output / "gene-to-protein-group-decisions.csv.gz")
    lists = known.dropna(subset=["membership_sources"]).copy()
    lists["source_tag"] = lists.membership_sources.str.split(";")
    lists = lists.explode("source_tag").drop_duplicates(
        ["source_tag", "protein_equivalence_group_id"]
    )
    source_rows, sets = [], {}
    for tag, part in lists.groupby("source_tag"):
        source_rows.append(
            {
                "source_tag": tag,
                "mapped_groups": len(part),
                **{
                    c: int(part.protein_group_category.eq(c).sum())
                    for c in ["CTA", "Unsafe", "Neither"]
                },
            }
        )
        sets[tag] = set(
            part.loc[part.protein_group_category.eq("CTA"), "protein_equivalence_group_id"]
        )
    sources = pd.DataFrame(source_rows)
    overlaps = pd.DataFrame(
        [
            {"source_a": a, "source_b": b, "selected_overlap": len(av & bv)}
            for a, av in sets.items()
            for b, bv in sets.items()
        ]
    )
    write_csv(sources, output / "source-category-counts.csv")
    write_csv(overlaps, output / "source-overlaps.csv")
    figures(output, result, funnel, independent, named, sources, overlaps)
    threshold_figure(output)
    summary = {
        "policy": POLICY_ID,
        "categories": result.category.value_counts().to_dict(),
        "selected_member_loci": int(result.loc[result.category.eq("CTA"), "member_count"].sum()),
        "new_selected_groups": len(added),
        "lost_selected_groups": int(
            (result.baseline_category.eq("CTA") & result.category.ne("CTA")).sum()
        ),
        "named_categories": named.set_index("protein_group_label").category.to_dict(),
        "order_independent": True,
        "review": reviewed,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "index.md").write_text(f"""# OncoRef 1.21.0 CTA/oncofetal research selection

**{summary["categories"].get("CTA", 0):,} selected groups ({summary["selected_member_loci"]:,} loci), {summary["categories"].get("Unsafe", 0):,} Unsafe and {summary["categories"].get("Neither", 0):,} Neither** across the complete 19,863-group reference. There are **{len(added)} new admissions** and **{summary["lost_selected_groups"]} lost selections** relative to the previous 343-group funnel.

This is a deliberately more permissive research prioritization policy, not a clinically validated safety model. PRAME informed the parameter choice and is therefore a calibration example, not independent validation. No gene name or clinical product bypasses the rules.

Use the [shared human, dog and mouse API](../../cta-selection.md). The [literature and implementation review](literature-review.md) records the corrected assumptions and code. [Slides](slides/cta-selection-1.21.0.pptx) and the [PDF report](cta-selection-1.21.0.pdf) use these same data.

CTA passes every gate. Unsafe indicates stronger evidence under this screening policy, not demonstrated toxicity. Neither includes intermediate RNA, missing required coverage and unresolved evidence. Every layer evaluates every group, even after an earlier hold.

The release review corrected **{reviewed["corrected_group_ihc_levels"]:,} grouped IHC levels** and **{reviewed["changed_layer_decisions"]} intermediate decisions**, with **{reviewed["changed_final_categories"]} final category changes** relative to the prior 463-group proposal. [Decision changes](review-decision-changes.csv) and [IHC corrections](review-ihc-level-changes.csv.gz) retain the audit.

## Uniform rules

Reproductive/embryonic association, tumor ≥2 TPM/pTPM in one sample, required RNA coverage, Medium/High HPA protein findings, validated external protein findings and clinical-safety exclusions retain their existing rules. DICE and corneal coverage remain optional; measured findings still count.

Missing/unusable-IHC RNA ceilings are **Vital 4, Protected 8, Manageable 16, specified renewing large-bowel epithelium 32**. Each negative-IHC confidence step doubles its row. The stronger RNA boundary is twice the candidate ceiling. Assays retain native units and support guards; CP10K uses the existing nominal /100 policy, not a TPM conversion.

Mapped adult HPA **fibroblasts** use Manageable. Colon/rectum distal enterocytes, intestinal goblet cells, Paneth cells and undifferentiated cells use the renewing-mucosa row. Other cells keep their anatomical tier. Undifferentiated cells retain a stem/progenitor caution: their admission is an explicit tolerance assumption, not evidence that the regenerative compartment is disposable. Cardiomyocytes, endothelial cells and airway basal cells receive no such relaxation. Bulk retains its organ tier.

Explicitly **Low** normal HPA staining or an Accepted literature finding with extracted **Low** intensity remains an observed protein flag when RNA passes. Medium/High HPA staining retains its hold/exclusion, and Accepted literature positives without an extracted intensity retain their exclusion. Low staining never supplies negative-IHC credit. Review access (abstract/full text) is recorded separately and never erases a known positive. Every imported finding is audited in [the protein-review table](imported-protein-review-decisions.csv), alongside [HPA protein-level decisions](hpa-protein-level-summary.csv.gz).

Author-described low expression and HPA Low retain different reporting scales. This policy tolerates both categories without claiming equal intensity, a known stained fraction or a calibrated H-score. Group aggregation retains the strongest observed level and its own confidence; negative-IHC credit requires complete member coverage and a positive profile control for every member.

## Contribution of each policy change

These are cumulative changes in the displayed order; admissions are not independent effects.

{markdown_table(components)}

## Named targets

{markdown_table(named[["protein_group_label", "category", "exclusion_reasons"]].fillna("All gates pass"))}

[All new admissions](newly-selected-groups.csv), [every changed classification](classification-changes.csv), [known CTA/source-list decisions](known-cta-decisions.csv.gz), [problematic target checks](problematic-target-checks.csv), [all group decisions](protein-group-classifications.csv.gz), and [all loci with identity and group decisions](gene-to-protein-group-decisions.csv.gz).

[Named-target normal evidence](named-normal-evidence.csv) shows the largest measured fraction of each assay's own tissue/IHC ceiling. This ratio compares policy margins and does not equate native assay units or predict toxicity.

## Filtering funnel

![Three-way funnel](figures/classification-funnel.png)

{markdown_table(funnel)}

![Threshold grid](figures/threshold-grid.png)

![Named antigens](figures/named-antigen-distribution.png)

![Source-list overlaps](figures/published-list-overlaps.png)

## Evidence and limitations

PRAME's 12-nTPM cardiac fibroblast measurement and 20.3-nTPM colon undifferentiated-cell measurement are retained; the new ceilings are 16 and 32. Its adrenal bulk 3.6 nTPM is below the new Vital ceiling 4. Its adrenal protein finding is explicitly **Low** in the primary abstract and remains a positive protein flag; it is not discarded because of limited source access. IMA203 clinical evidence remains separate and does not participate in admission. See the [target evidence assessment](../cta-target-evidence-20261027/index.md) and [primary adrenal-protein report](https://pubmed.ncbi.nlm.nih.gov/33890816/).

Controlled fibroblast depletion has been tolerated in mouse studies; that motivates separating stroma from cardiomyocytes but does not validate these human cutoffs. [Primary study](https://elifesciences.org/articles/69854). Direct epithelial targeting can produce severe colitis, including a dose-limiting signal in a small CEA-directed TCR trial; allowing greater mucosal expression is a research risk preference. [Primary trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC3048186/).

No numerical cutoff here predicts peptide presentation or injury severity. The permissive mucosal assumption and incomplete older protein reviews require particular scrutiny. These outputs do not overwrite the [prior 343-group policy](../cta-harmonized-20261027/index.md), released classifications or forbidden sequence references.

## Reproduction

Run `python scripts/cta_compartment_report.py`, then `python scripts/build_cta_selection_release.py`. [Changed cell compartments](changed-compartment-contexts.csv.gz), [bulk contexts](bulk-contexts.csv.gz), [immune/corneal contexts](immune-cornea-contexts.csv.gz), [per-layer effects](independent-layer-effects.csv), [counts](summary.json), and [provenance](manifest.json) retain the decisions. Multi-locus values are summed before assessment; unpaired DICE and corneal detection uncertainty remain explicit. Figures also have SVG versions. The earlier 463-group audit preserves its original implementation bytes under `source-code/` and remains reproducible from commit `d79f4d4`.
""")
    inputs = [
        Path(__file__),
        ROOT / "oncoref/compartment_rna_policy.py",
        ROOT / "scripts/cta_harmonized_report.py",
        ROOT / "scripts/cta_cell_policy_report.py",
        ROOT / "scripts/cta_quantitative_threshold_report.py",
        ROOT / "scripts/cta_expanded_normal.py",
        ROOT / "scripts/cta_report_common.py",
        ROOT / "oncoref/normal_cell_data.py",
        ROOT / "oncoref/normal_cell_screen.py",
        ROOT / "oncoref/normal_tissues.py",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/normal-tissue-mappings.csv",
        ROOT / "oncoref/cell_rna_policy.py",
        ROOT / "oncoref/cta_funnel.py",
        ROOT / "oncoref/data/cta-policy-findings.csv",
        output / "protein-intensity-input.csv",
        BASE / "manifest.json",
        BASE / "protein-group-classifications.csv.gz",
        BASE / "known-cta-decisions.csv.gz",
        BASE / "gene-to-protein-group-decisions.csv.gz",
        BASE / "group-normal/manifest.json",
        EXPANDED / "normal/manifest.json",
    ]
    # Presentation/package outputs have their own final receipt. Excluding them
    # here prevents a report -> package -> report provenance cycle on reruns.
    ancillary = {
        "manifest.json",
        "release-manifest.json",
        "species-source-inventory.csv",
        "literature-review.md",
        "release-notes.md",
        "validation.json",
        "presentation-data.json",
        "cta-selection-1.21.0.pdf",
    }
    manifest = {
        "inputs": input_hashes(inputs),
        "sources": receipts,
        "outputs": {
            str(p.relative_to(output)): sha256(p)
            for p in sorted(output.rglob("*"))
            if p.is_file()
            and p.name not in ancillary
            and "slides" not in p.relative_to(output).parts
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reuse-clusters",
        action="store_true",
        help="Reuse cluster replay only after verifying all input and output hashes",
    )
    build(reuse_clusters=parser.parse_args().reuse_clusters)
