#!/usr/bin/env python3
"""Replay a provisional compartment-aware CTA policy across the full reference."""

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
from cta_harmonized_report import ASSAYS, decision_layers, figures  # noqa: E402
from cta_quantitative_threshold_report import positive_profile_controls  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import cell_rna_policy, hpa, normal_cell_screen, reference_data  # noqa: E402
from oncoref.compartment_rna_policy import (  # noqa: E402
    IHC_RANKS,
    POLICY_ID,
    TIER_BASES,
    assess_contexts,
    assess_imported_protein_reviews,
)
from oncoref.cta_funnel import classify_layers  # noqa: E402

BASE = ROOT / "docs/audits/cta-harmonized-20261027"
EXPANDED = ROOT / "docs/audits/cta-expanded-funnel-20261027"
OUTPUT = ROOT / "docs/audits/cta-compartment-policy-20261027"
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
    d["protein_exclusion"] = normal & d.Level.isin(cell_rna_policy.POSITIVE - {"Low"})
    d["reliable_protein_exclusion"] = d.protein_exclusion & d.Reliability.isin(
        cell_rna_policy.GRADES
    )
    summary = d.groupby("gene_id")[
        ["low_protein", "protein_exclusion", "reliable_protein_exclusion"]
    ].any()
    write_csv(summary.reset_index(), output / "hpa-protein-level-summary.csv.gz")
    members = facts.member_gene_ids.map(json.loads)
    low_ids = set(summary.index[summary.low_protein])
    other_ids = set(summary.index[summary.protein_exclusion])
    only_low = members.map(
        lambda group: bool(set(group) & low_ids) and not bool(set(group) & other_ids)
    )
    # Relax only a baseline HPA finding whose observed positives are all Low.
    # Unexplained prior findings cannot disappear merely because this extraction
    # lacks a matching row or uses a different source annotation.
    facts["hpa_low_only_relaxation"] = facts.ihc_detected & only_low
    facts["ihc_detected"] &= ~only_low
    facts["ihc_reliable_detected"] &= ~only_low
    return low_ids


def linked_assays(facts: pd.DataFrame) -> None:
    """Replay uniformly doubled RNA limits, retaining unresolved nonnumeric holds.

    The maximum normalized ratio is sufficient when every ceiling scales equally.
    Protein/unknown-context residuals are retained separately. Group DICE and cornea
    are replaced below using their uncertainty bounds, rather than assuming paired
    donors or cell detections. Bulk and clusters are also replaced by full contexts.
    """
    single = pd.read_csv(EXPANDED / "normal/normal-peak-evidence.csv.gz")
    grouped = pd.read_csv(BASE / "group-normal/normal-peak-evidence.csv.gz")
    peaks = pd.concat([single[~single.gene_id.isin(grouped.gene_id)], grouped])
    for assay in ASSAYS:
        ratio = facts.representative_gene_id.map(
            peaks[peaks.assay.eq(assay)].set_index("gene_id").ratio
        )
        old_hold, old_strong = facts["hold:" + assay], facts["unsafe:" + assay]
        residual_hold = old_hold & ~ratio.gt(1)
        residual_strong = old_strong & ~ratio.ge(2)
        facts["hold:" + assay] = ratio.gt(2) | residual_hold
        facts["unsafe:" + assay] = ratio.ge(4) | residual_strong


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
    facts: pd.DataFrame, output: Path, *, reuse: bool = False
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
        ]
    )
    signature = hashlib.sha256(
        (inspect.getsource(cluster_assessments) + inspect.getsource(aggregate_contexts)).encode()
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
    summaries, changed = [], []

    def collect(d: pd.DataFrame) -> None:
        result = assess_contexts(d)
        summaries.append(aggregate_contexts(result))
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
    collect(pd.read_csv(BASE / "group-normal/cluster-contexts.csv.gz", low_memory=False))
    write_csv(pd.concat(changed, ignore_index=True), output / "changed-compartment-contexts.csv.gz")
    aggregated = pd.concat(summaries).groupby(level=0).any()
    write_csv(aggregated.reset_index(), output / "cluster-assessments.csv.gz")
    cache = {
        "inputs": cache_inputs,
        "sources": sources,
        "function_sha256": signature,
        "outputs": {
            name: sha256(output / name)
            for name in ["cluster-assessments.csv.gz", "changed-compartment-contexts.csv.gz"]
        },
    }
    cache_path.write_text(json.dumps(cache, indent=2) + "\n")
    return aggregated, sources


def bulk_assessments(facts: pd.DataFrame, output: Path) -> pd.DataFrame:
    """Reassess every matched bulk context, preserving missing-measurement holds."""
    single = pd.read_csv(EXPANDED / "normal/bulk-contexts.csv.gz", low_memory=False)
    grouped = pd.read_csv(BASE / "group-normal/bulk-contexts.csv.gz", low_memory=False)
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


def uncertain_group_assays(facts: pd.DataFrame) -> None:
    """Recompute unpaired DICE and corneal detection uncertainty at new limits."""
    multi = facts.member_count.gt(1)
    dice = pd.read_csv(BASE / "group-normal/dice-rank-bounds.csv.gz")
    d = assess_contexts(dice)
    d["candidate_exclusion"] |= d.measured_members.eq(d.member_count) & d.value_upper.gt(
        d.candidate_ceiling
    )
    apply_assay(facts, aggregate_contexts(d), "DICE immune", subset=multi)
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


def build(output: Path = OUTPUT, *, reuse_clusters: bool = False) -> None:
    """Publish complete decisions, admission changes, known targets and figures."""
    output.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(BASE / "protein-group-classifications.csv.gz", low_memory=False)
    facts["baseline_category"] = facts.category
    facts["baseline_exclusion_reasons"] = facts.exclusion_reasons
    facts = facts.rename(columns={"member_gene_categories": "baseline_member_gene_categories"})
    facts["baseline_accepted_normal_protein"] = facts["Accepted normal protein"]
    linked_assays(facts)
    apply_assay(facts, bulk_assessments(facts, output), "Normal bulk")
    clusters, receipts = cluster_assessments(facts, output, reuse=reuse_clusters)
    apply_assay(facts, clusters, "HPA tissue clusters")
    uncertain_group_assays(facts)
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
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "index.md").write_text(f"""# Provisional compartment-aware CTA candidate set

**{summary["categories"].get("CTA", 0):,} selected groups ({summary["selected_member_loci"]:,} loci), {summary["categories"].get("Unsafe", 0):,} Unsafe and {summary["categories"].get("Neither", 0):,} Neither** across the complete 19,863-group reference. There are **{len(added)} new admissions** and **{summary["lost_selected_groups"]} lost selections** relative to the previous 343-group funnel.

This is a deliberately more permissive research prioritization policy, not a clinically validated safety model. PRAME informed the parameter choice and is therefore a calibration example, not independent validation. No gene name or clinical product bypasses the rules.

## Uniform rules

Reproductive/embryonic association, tumor ≥2 TPM/pTPM in one sample, required RNA coverage, Medium/High HPA protein findings, validated external protein findings and clinical-safety exclusions retain their existing rules. DICE and corneal coverage remain optional; measured findings still count.

Missing/unusable-IHC RNA ceilings are **Vital 4, Protected 8, Manageable 16, specified renewing large-bowel epithelium 32**. Each negative-IHC confidence step doubles its row. The stronger RNA boundary is twice the candidate ceiling. Assays retain native units and support guards; CP10K uses the existing nominal /100 policy, not a TPM conversion.

Mapped adult HPA **fibroblasts** use Manageable. Colon/rectum distal enterocytes, intestinal goblet cells, Paneth cells and undifferentiated cells use the renewing-mucosa row. Other cells keep their anatomical tier. Undifferentiated cells retain a stem/progenitor caution: their admission is an explicit tolerance assumption, not evidence that the regenerative compartment is disposable. Cardiomyocytes, endothelial cells and airway basal cells receive no such relaxation. Bulk retains its organ tier.

Explicitly **Low** normal HPA staining or an Accepted literature finding with extracted **Low** intensity remains an observed protein flag when RNA passes. Medium/High HPA staining retains its hold/exclusion, and Accepted literature positives without an extracted intensity retain their exclusion. Low staining never supplies negative-IHC credit. Review access (abstract/full text) is recorded separately and never erases a known positive. Every imported finding is audited in [the protein-review table](imported-protein-review-decisions.csv), alongside [HPA protein-level decisions](hpa-protein-level-summary.csv.gz).

## Contribution of each policy change

These are cumulative changes in the displayed order; admissions are not independent effects.

{markdown_table(components)}

## Named targets

{markdown_table(named[["protein_group_label", "category", "exclusion_reasons"]].fillna("All gates pass"))}

[All new admissions](newly-selected-groups.csv), [every changed classification](classification-changes.csv), [known CTA/source-list decisions](known-cta-decisions.csv.gz), [problematic target checks](problematic-target-checks.csv), [all group decisions](protein-group-classifications.csv.gz), and [all loci with identity and group decisions](gene-to-protein-group-decisions.csv.gz).

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

Run `python scripts/cta_compartment_report.py`. [Changed cell compartments](changed-compartment-contexts.csv.gz), [bulk contexts](bulk-contexts.csv.gz), [per-layer effects](independent-layer-effects.csv), [counts](summary.json), and [provenance](manifest.json) retain the decisions. Multi-locus values are summed before assessment; unpaired DICE and corneal detection uncertainty remain explicit. Figures also have SVG versions.
""")
    inputs = [
        Path(__file__),
        ROOT / "oncoref/compartment_rna_policy.py",
        ROOT / "scripts/cta_harmonized_report.py",
        ROOT / "scripts/cta_cell_policy_report.py",
        ROOT / "scripts/cta_quantitative_threshold_report.py",
        ROOT / "oncoref/normal_cell_screen.py",
        ROOT / "oncoref/normal_tissues.py",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/normal-tissue-mappings.csv",
        ROOT / "oncoref/cell_rna_policy.py",
        ROOT / "oncoref/cta_funnel.py",
        ROOT / "oncoref/data/cta-policy-findings.csv",
        output / "protein-intensity-input.csv",
        BASE / "manifest.json",
        BASE / "group-normal/manifest.json",
        EXPANDED / "normal/manifest.json",
    ]
    manifest = {
        "inputs": input_hashes(inputs),
        "sources": receipts,
        "outputs": {
            str(p.relative_to(output)): sha256(p)
            for p in sorted(output.rglob("*"))
            if p.is_file() and p.name != "manifest.json"
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
