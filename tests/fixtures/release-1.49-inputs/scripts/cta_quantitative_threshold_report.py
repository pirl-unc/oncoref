#!/usr/bin/env python3
"""Compare tumor admission and local IHC-conditioned or linked normal-RNA gates.

The archived 447-gene proposal stays reproducible. This audit applies the existing
experimental bulk and cell IHC rules as an additional exclusion layer. It never
changes the released expression categories or the forbidden-protein reference.
The linked-tiers mode re-evaluates all eligible genes without the old RNA caps.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from cta_cell_policy_report import scenario_thresholds, summarize  # noqa: E402
from cta_report_common import sha256 as digest  # noqa: E402
from cta_report_common import write_csv as write  # noqa: E402

from oncoref import (  # noqa: E402
    cell_rna_policy,
    cta_quantitative,
    cta_tissues,
    hpa,
    hpa_cancer,
    normal_cell_data,
    normal_cell_screen,
    reference_data,
)
from oncoref.cta_normal_screen import adaptive_normal_rna_limits  # noqa: E402

POLICY = "oncoref.cta-quantitative-ihc-comparison.2026-10-27.1"
NEUTRAL_POLICY = "oncoref.cta-quantitative-tissue-context.2026-10-27.1"
FACTORS = {"Missing/uncertain": 1, "Approved": 1.25, "Supported": 1.5, "Enhanced": 2}
CONTROLS = [
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
]


def tumor_comparison(facts: pd.DataFrame, tumor: pd.DataFrame, *, floor=30) -> pd.DataFrame:
    """Separate a sample maximum from cohort prevalence, preserving native pTPM."""
    if not np.isfinite(floor) or floor <= 0:
        raise ValueError("Positive finite native-pTPM floor required")
    if tumor.duplicated(["gene_id", "cancer_code", "cohort"]).any():
        raise ValueError("Duplicate gene/cancer/cohort summaries")
    observed = tumor[tumor.samples.gt(0)].copy()
    if observed.max_ptpm.isna().any() or not np.isfinite(observed.max_ptpm).all():
        raise ValueError("Measured cohorts require finite maxima")
    best = (
        observed.sort_values(["max_ptpm", "cancer_code", "cohort"], ascending=[False, True, True])
        .drop_duplicates("gene_id")
        .set_index("gene_id")
    )
    result = facts.copy()
    for col in ["max_ptpm", "cancer_code", "cohort"]:
        result["tumor_maximum_" + col] = result.gene_id.map(best[col])
    result["tumor_maximum_floor_ptpm"] = floor
    result["tumor_maximum_positive"] = result.tumor_maximum_max_ptpm.ge(floor)
    other = result.positive_context & result.rna_coverage_complete & ~result.ihc_detected
    for col in [
        "adult_bulk_exclusion",
        "adult_cell_exclusion",
        "immune_exclusion",
        "cornea_exclusion",
    ]:
        other &= ~result[col]
    if not (other & result.tumor_positive).equals(result.selected):
        raise ValueError("Frozen baseline does not replay")
    result["maximum_tumor_selected"] = other & result.tumor_maximum_positive
    return result


def positive_profile_controls(ihc: pd.DataFrame) -> set[str]:
    """Measured reproductive-positive HPA profiles, never a named-gene allowlist."""
    d = ihc.rename(columns={"Tissue": "tissue", "Cell type": "cell_type"}).copy()
    d["source"], d["version"], d["taxonomy_id"] = "hpa_normal_tissue", "v23", 9606
    d = normal_cell_screen.contexts(d)
    return set(
        d.loc[
            d.category.eq("reproductive")
            & d.Level.isin(cell_rna_policy.POSITIVE)
            & d.Reliability.isin(cell_rna_policy.GRADES),
            "Gene",
        ]
    )


def bulk_limits(
    bulk: pd.DataFrame, ihc: pd.DataFrame, *, positive_controls, conflicts, legacy_cap=10
) -> pd.DataFrame:
    """Local HPA negative-IHC ladder, with an optional historical uniform cap.

    Complete routine tissue IHC is required; one cell type cannot clear a whole
    tissue. The prior selectivity ladder is deliberately not an admission gate.
    """
    routine = tuple(
        ihc[["Tissue", "Cell type"]].dropna().drop_duplicates().itertuples(index=False, name=None)
    )
    groups = dict(iter(ihc.groupby("Gene")))
    output = []
    for gid, rna in bulk.groupby("Gene", sort=True):
        local = groups.get(gid, ihc.iloc[:0])
        rna = rna.rename(
            columns={"Gene": "Ensembl_Gene_ID", "Tissue": "tissue", "nTPM": "rna_ntpm"}
        ).copy()
        rna["source_name"], rna["source_version"], rna["modality"] = (
            "hpa_rna_consensus",
            "v23",
            "bulk_rna",
        )
        rna["cell_type"], rna["ihc_level"], rna["ihc_reliability"] = "", "", ""
        rna["valid_measurement"] = rna.rna_ntpm.notna()
        protein = local.rename(
            columns={
                "Gene": "Ensembl_Gene_ID",
                "Tissue": "tissue",
                "Cell type": "cell_type",
                "Level": "ihc_level",
                "Reliability": "ihc_reliability",
            }
        ).copy()
        protein["source_name"], protein["source_version"], protein["modality"] = (
            "hpa_normal_tissue",
            "v23",
            "ihc",
        )
        protein["rna_ntpm"], protein["valid_measurement"] = np.nan, True
        audit = adaptive_normal_rna_limits(
            pd.concat([rna, protein], ignore_index=True),
            routine_ihc_contexts=routine,
            reviewed_unreliable=gid in conflicts,
        ).reset_index(drop=True)
        audit = audit[~audit.tissue.isin(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES)].copy()
        audit["profile_positive_control"] = gid in positive_controls
        if gid not in positive_controls:
            audit["rna_ceiling_ntpm"] = audit.baseline_ceiling_ntpm
            audit["relaxed"] = False
        audit["effective_ceiling_ntpm"] = (
            audit.rna_ceiling_ntpm.clip(upper=legacy_cap)
            if legacy_cap is not None
            else audit.rna_ceiling_ntpm
        )
        # Already-permitted compartments retain the original bulk cap, never a
        # new IHC waiver. Unresolved anatomy is a hold, not an infinite ceiling.
        if legacy_cap is not None:
            audit.loc[audit.screening_tier.eq("permitted"), "effective_ceiling_ntpm"] = legacy_cap
        audit["local_exclusion"] = (
            audit.effective_ceiling_ntpm.isna()
            | audit.rna_ntpm.isna()
            | audit.rna_ntpm.gt(audit.effective_ceiling_ntpm)
        )
        if legacy_cap is None:
            audit.loc[audit.screening_tier.eq("permitted"), "local_exclusion"] = False
            audit["reference_review_ceiling_ntpm"] = (
                audit.effective_ceiling_ntpm * cell_rna_policy.REFERENCE_REVIEW_MULTIPLIER
            )
            audit["reference_review"] = audit.rna_ntpm.ge(audit.reference_review_ceiling_ntpm)
        output.append(audit)
    return pd.concat(output, ignore_index=True)


def rebase_reproductive_immune(facts: pd.DataFrame) -> pd.DataFrame:
    """Reprocess all genes before selecting a pool; neutral contexts may rescue genes.

    The original 426 observed clusters remain the coverage contract. Neutral
    observations count as observed coverage but never as a normal exclusion or
    as positive reproductive admission. Counts of screened contexts are separate.
    """
    ids = set(facts.gene_id)
    source = reference_data.provenance("hpa_single_cell_tissue", verify_content=True)
    parts = []
    for chunk in pd.read_csv(source["path"], sep="\t", chunksize=250000):
        d = (
            chunk[chunk.Gene.isin(ids)]
            .rename(
                columns={
                    "Gene": "gene_id",
                    "Tissue": "tissue",
                    "Cell type": "cell_type",
                    "nTPM": "value",
                    "Read count": "read_count",
                }
            )
            .copy()
        )
        if d.empty:
            continue
        d["source"], d["version"], d["taxonomy_id"], d["stage"] = (
            "hpa_single_cell_tissue",
            "v23",
            9606,
            "adult",
        )
        d = normal_cell_screen.contexts(d)
        d = d[~d.category.isin(cell_rna_policy.PERMITTED)].copy()
        d = cell_rna_policy.reproductive_immune_contexts(d)
        normal = ~d.neutral_reproductive_immune
        d["screened"] = normal
        d["supported"] = d.value.where(normal & d.read_count.ge(20))
        d["raw_screened"] = d.value.where(normal)
        d["raw_neutral"] = d.value.where(~normal)
        parts.append(
            d.groupby("gene_id").agg(
                cell_records=("value", "size"),
                cell_screened_records=("screened", "sum"),
                cell_neutral_records=("neutral_reproductive_immune", "sum"),
                cell_supported_records=("supported", "count"),
                cell_max_ntpm=("supported", "max"),
                cell_raw_max_ntpm=("raw_screened", "max"),
                cell_neutral_max_ntpm=("raw_neutral", "max"),
            )
        )
    combined = pd.concat(parts)
    counts = [c for c in combined if c.endswith("records")]
    combined = combined.groupby(level=0).agg({c: "sum" if c in counts else "max" for c in combined})
    d = facts.copy()
    d["historical_uniform_selected"] = d.selected
    d["historical_cell_max_ntpm"] = d.cell_max_ntpm
    for col in combined:
        d[col] = d.gene_id.map(combined[col])
    # The pinned IHC export has no additional blood-lineage contexts sampled in
    # reproductive tissues. Assert this rather than silently retaining a veto
    # from a neutral source if the source inventory changes in a future build.
    h = (
        hpa._read_hpa("hpa_normal_tissue", "v23")[["Tissue", "Cell type"]]
        .dropna()
        .drop_duplicates()
    )
    h = h.rename(columns={"Tissue": "tissue", "Cell type": "cell_type"})
    h["source"], h["version"], h["taxonomy_id"], h["stage"] = (
        "hpa_normal_tissue",
        "v23",
        9606,
        "adult",
    )
    h = normal_cell_screen.contexts(h)
    h = h[~h.category.isin(cell_rna_policy.PERMITTED)].copy()
    if cell_rna_policy.reproductive_immune_contexts(h).neutral_reproductive_immune.any():
        raise ValueError("Recompute IHC exclusions for new neutral IHC contexts")
    result = cta_quantitative.classify(d)
    result["policy"] = NEUTRAL_POLICY
    return result


def build(
    out: Path,
    *,
    neutral_reproductive_immune: bool = False,
    linked_tiers: bool = False,
    tumor_floor=30,
) -> None:
    """Build the versioned audit tables, report and source provenance."""
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    facts = cta_quantitative.results()
    if not cta_quantitative.classify(facts).selected.equals(facts.selected):
        raise ValueError("Quantitative baseline has changed")
    if neutral_reproductive_immune:
        print("Reprocessing every reference gene with local immune-context neutrality", flush=True)
        facts = rebase_reproductive_immune(facts)
        write(facts, out / "context-reprocessed-facts.csv.gz")
    tumor = hpa_cancer.hpa_cancer_rna_prevalence()
    genes = tumor_comparison(facts, tumor, floor=tumor_floor)
    # This is an exclusion-only layer. Every gene outside this union already
    # fails an unchanged gate, so it cannot be rescued by local negative IHC.
    pool = genes.selected | genes.maximum_tumor_selected
    if linked_tiers:
        # Re-evaluate every otherwise eligible gene. Old 10-bulk/30-cluster/
        # 5-immune/1-cornea caps must not silently override the shared ladder.
        pool = (
            genes.positive_context
            & genes.rna_coverage_complete
            & ~genes.ihc_detected
            & (
                genes.tumor_positive
                | genes.tumor_maximum_positive
                | genes.tumor_maximum_max_ptpm.ge(30)
            )
        )
    ids = set(genes.loc[pool, "gene_id"])
    print(f"Assessing the {len(ids)} genes eligible under either tumor rule", flush=True)
    inputs = {
        "oncoref/data/cta-quantitative-gene-sets.csv.gz": digest(
            ROOT / "oncoref/data/cta-quantitative-gene-sets.csv.gz"
        )
    }
    receipts = []

    def receipt(source):
        r = reference_data.provenance(source, verify_content=True)
        receipts.append({k: r[k] for k in ["name", "version", "url", "sha256"]})
        return r

    for source in ["hpa_rna_consensus", "hpa_normal_tissue", "hpa_cancer_rna"]:
        receipt(source)
    ihc = hpa._read_hpa("hpa_normal_tissue", "v23").dropna(subset=["Tissue", "Cell type", "Level"])
    controls = positive_profile_controls(ihc)
    conflicts = set(facts.loc[facts["Antibody concern"], "gene_id"])
    bulk = hpa._read_hpa("hpa_rna_consensus", "v23")
    bulk = bulk[bulk.Gene.isin(ids)]
    audit = bulk_limits(
        bulk,
        ihc,
        positive_controls=controls,
        conflicts=conflicts,
        legacy_cap=None if linked_tiers else 10,
    )
    write(audit, out / "bulk-contexts.csv.gz")
    print("Local bulk IHC comparison complete", flush=True)
    frames, raw_clusters = [], []
    r = receipt("hpa_single_cell_tissue")
    for chunk in pd.read_csv(r["path"], sep="\t", chunksize=250000):
        d = chunk[chunk.Gene.isin(ids)].copy()
        if d.empty:
            continue
        d = d.rename(
            columns={
                "Gene": "gene_id",
                "Tissue": "tissue",
                "Cell type": "cell_type",
                "nTPM": "value",
                "Read count": "read_count",
            }
        )
        d["source"], d["version"], d["taxonomy_id"] = "hpa_single_cell_tissue", "v23", 9606
        d["state"], d["stage"], d["unit"] = "unspecified", "adult", "nTPM"
        raw_clusters.append(d)
    raw = normal_cell_screen.contexts(pd.concat(raw_clusters, ignore_index=True))
    raw = raw[~raw.category.isin(cell_rna_policy.PERMITTED)].copy()
    frames.append(summarize(raw, support=raw.read_count.ge(20)))
    write(
        (cell_rna_policy.reproductive_immune_contexts(raw) if neutral_reproductive_immune else raw)[
            raw.gene_id.isin(set(facts.loc[facts.symbol.isin(CONTROLS), "gene_id"]))
        ],
        out / "named-cell-measurements.csv",
    )
    for source, unit in [("hpa_immune_cell", "nTPM"), ("hpa_immune_monaco", "TPM")]:
        r = receipt(source)
        d = normal_cell_data.observations(source, genes=ids)
        if not d.unit.eq(unit).all():
            raise ValueError("Immune source units changed")
        d = normal_cell_screen.contexts(d)
        frames.append(summarize(d, support=d.value.notna()))
    for source in normal_cell_data.sources().source:
        if source.startswith("dice_"):
            receipt(source)
            d = normal_cell_data.dice_observations(source, genes=ids)
            # PAR loci are not present in this eligible union. Never pool their
            # source records as if they were additional independent samples.
            if d.duplicated(["gene_id", "replicate_slot"]).any():
                raise ValueError("Ambiguous DICE identity requires separate-record handling")
            d = normal_cell_screen.contexts(d)
            frames.append(summarize(d, support=d.value.notna(), rank=3, unique="replicate_slot"))
    for source in ["cornea_collin_adult_counts", "cornea_collin_adultcornea_cells"]:
        receipt(source)
    d = normal_cell_screen.contexts(normal_cell_data.cornea_summary(genes=ids))
    frames.append(
        summarize(d, support=d.positive_cells.ge(5) & d.read_count.ge(20), rank=2, unique="sample")
    )
    contexts = pd.concat(frames, ignore_index=True)
    if neutral_reproductive_immune or linked_tiers:
        contexts = normal_cell_screen.contexts(
            contexts.drop(columns=["category", "mapping_status"])
        )
    if neutral_reproductive_immune:
        contexts = cell_rna_policy.reproductive_immune_contexts(contexts)
    if linked_tiers:
        contexts = cell_rna_policy.protected_immune_contexts(contexts)
    contexts = cell_rna_policy.match_negative_ihc(
        contexts,
        ihc,
        positive_controls=controls,
        conflicting_genes=conflicts,
    )
    thresholds = cell_rna_policy.linked_thresholds() if linked_tiers else scenario_thresholds()
    factors = dict.fromkeys(FACTORS, 1) if linked_tiers else FACTORS
    assessed = cell_rna_policy.assess(
        contexts, thresholds, ihc_factors=factors, inclusive_ceiling=linked_tiers
    )
    assessed["unresolved_context"] = ~assessed.mapping_status.eq("mapped")
    assessed["local_exclusion"] = assessed.cta_excluded | assessed.unresolved_context
    if neutral_reproductive_immune:
        assessed.loc[assessed.neutral_reproductive_immune, "local_exclusion"] = False
    write(assessed, out / "cell-contexts.csv.gz")
    write(assessed[assessed.local_exclusion], out / "cell-exclusions.csv.gz")
    write(thresholds, out / "cell-thresholds.csv")
    genes["ihc_contexts_evaluated"] = pool
    genes["ihc_bulk_exclusion"] = genes.gene_id.map(
        audit.groupby("Ensembl_Gene_ID").local_exclusion.any()
    )
    genes["ihc_cell_exclusion"] = genes.gene_id.map(
        assessed.groupby("gene_id").local_exclusion.any()
    )
    # Every potentially admitted gene must have been assessed in both layers.
    if genes.loc[pool, ["ihc_bulk_exclusion", "ihc_cell_exclusion"]].isna().any().any():
        raise ValueError("Missing conditional-IHC comparison")
    local_pass = genes.ihc_bulk_exclusion.eq(False) & genes.ihc_cell_exclusion.eq(False)
    genes["ihc_conditioned_selected"] = genes.selected & local_pass
    genes["maximum_tumor_ihc_conditioned_selected"] = genes.maximum_tumor_selected & local_pass
    if linked_tiers:
        genes["linked_recurrence_selected"] = pool & genes.tumor_positive & local_pass
        genes["linked_maximum_tumor_selected"] = pool & genes.tumor_maximum_positive & local_pass
        genes["linked_tumor30_selected"] = pool & genes.tumor_maximum_max_ptpm.ge(30) & local_pass
        # Quantify the separate effect of linking the corneal count scale.
        # This is explicitly a sensitivity, not an unreported assay override.
        legacy_cornea = assessed.copy()
        cornea = legacy_cornea.assay.eq("cornea_collin_adult_counts")
        legacy_cornea.loc[cornea, "local_exclusion"] = (
            legacy_cornea.loc[cornea, "value"].ge(1)
            | legacy_cornea.loc[cornea, "unresolved_context"]
        )
        genes["legacy_cornea_cell_exclusion"] = genes.gene_id.map(
            legacy_cornea.groupby("gene_id").local_exclusion.any()
        )
        sensitivity_pass = (
            pool & genes.ihc_bulk_exclusion.eq(False) & genes.legacy_cornea_cell_exclusion.eq(False)
        )
        genes["linked_recurrence_legacy_cornea_selected"] = sensitivity_pass & genes.tumor_positive
        genes["linked_maximum_tumor_legacy_cornea_selected"] = (
            sensitivity_pass & genes.tumor_maximum_positive
        )
        genes["linked_tumor30_legacy_cornea_selected"] = (
            sensitivity_pass & genes.tumor_maximum_max_ptpm.ge(30)
        )
    if neutral_reproductive_immune and not linked_tiers:
        # Explicit sensitivity only: allow supported protected-tier cell RNA
        # below 30 nTPM, while keeping vital, immune, bulk and protein gates.
        alternative = assessed.copy()
        protected = alternative.assay.eq("hpa_single_cell_tissue") & alternative.category.eq(
            "protected"
        )
        alternative.loc[protected, "local_exclusion"] = (
            alternative.loc[protected, "value"].ge(30)
            | alternative.loc[protected, "ihc_positive"]
            | alternative.loc[protected, "unresolved_context"]
        )
        genes["protected30_cell_exclusion"] = genes.gene_id.map(
            alternative.groupby("gene_id").local_exclusion.any()
        )
        genes["protected30_selected"] = (
            genes.selected
            & genes.ihc_bulk_exclusion.eq(False)
            & genes.protected30_cell_exclusion.eq(False)
        )
        genes["maximum_tumor_protected30_selected"] = (
            genes.maximum_tumor_selected
            & genes.ihc_bulk_exclusion.eq(False)
            & genes.protected30_cell_exclusion.eq(False)
        )
    policy = (
        cell_rna_policy.LINKED_POLICY
        if linked_tiers
        else NEUTRAL_POLICY
        if neutral_reproductive_immune
        else POLICY
    )
    genes["comparison_policy"] = policy
    write(genes, out / "gene-comparison.csv.gz")
    if linked_tiers:
        write(genes[genes.linked_maximum_tumor_selected], out / "selected-genes.csv")
    write(genes[genes.symbol.isin(CONTROLS)], out / "named-gene-comparison.csv")
    write(hpa_cancer.hpa_cancer_rna_cohorts(), out / "tumor-cohorts.csv")
    named_ids = set(facts.loc[facts.symbol.isin(CONTROLS), "gene_id"])
    broad_path = ROOT / "docs/audits/cta-neither-review-20261024/rna-context-maxima.csv.gz"
    inputs[str(broad_path.relative_to(ROOT))] = digest(broad_path)
    broad = pd.concat(
        [
            d[d.gene_id.isin(named_ids) & ~d.permitted]
            for d in pd.read_csv(broad_path, chunksize=100000, low_memory=False)
        ],
        ignore_index=True,
    )
    broad["symbol"] = broad.gene_id.map(facts.set_index("gene_id").symbol)
    if neutral_reproductive_immune:
        broad["taxonomy_id"], broad["stage"] = 9606, "adult"
        broad = normal_cell_screen.contexts(broad.drop(columns=["category", "mapping_status"]))
        broad = cell_rna_policy.reproductive_immune_contexts(broad)
    write(broad, out / "named-normal-contexts.csv.gz")
    peak_pool = broad[~broad.neutral_reproductive_immune] if neutral_reproductive_immune else broad
    peaks = peak_pool.sort_values("maximum", ascending=False).drop_duplicates(["gene_id", "unit"])
    write(peaks, out / "named-normal-maxima.csv")
    write(bulk[bulk.Gene.isin(named_ids)], out / "named-bulk-measurements.csv")
    rows = []
    scenarios = [
        "selected",
        "maximum_tumor_selected",
        "ihc_conditioned_selected",
        "maximum_tumor_ihc_conditioned_selected",
    ]
    if linked_tiers:
        scenarios = [
            "linked_recurrence_selected",
            "linked_maximum_tumor_selected",
            "linked_tumor30_selected",
            "linked_recurrence_legacy_cornea_selected",
            "linked_maximum_tumor_legacy_cornea_selected",
            "linked_tumor30_legacy_cornea_selected",
        ]
    elif neutral_reproductive_immune:
        scenarios += ["protected30_selected", "maximum_tumor_protected30_selected"]
    for col in scenarios:
        s = genes[genes[col]]
        rows.append(
            {
                "scenario": col,
                "genes": len(s),
                "controls": int(s.symbol.isin(CONTROLS).sum()),
                **s["Expression category"].value_counts().to_dict(),
            }
        )
    summary = {
        "policy": policy,
        "scenarios": rows,
        "ihc_evaluated_genes": len(ids),
        "maximum_tumor_lost": int((genes.selected & ~genes.maximum_tumor_selected).sum()),
        "maximum_tumor_gained": int((~genes.selected & genes.maximum_tumor_selected).sum()),
        "exact_cell_ihc_matches": int(assessed.ihc_exact_match.sum()),
        "cell_ihc_relaxations": int(assessed.negative_ihc_grade.ne("Missing/uncertain").sum()),
        "assessed_cell_contexts": len(assessed),
        "ihc_factors": factors,
        "maximum_tumor_floor_ptpm": tumor_floor,
        "scope": "Exclusion-only IHC comparison for the union of genes eligible under either tumor rule; other genes retain their unchanged failing gates.",
        "tumor_source": "HPA 25.1 native pTPM; 21 TCGA and 10 validation cohorts; other OncoRef cohorts not queried.",
        "bulk_rule": "Existing local HPA negative-IHC absolute ladder, capped at 10 nTPM; no legacy selectivity gate.",
        "cell_rule": "Existing illustrative native-unit dual-cutoff matrix and exact-cell negative-IHC multipliers; original 30 nTPM adult-cell cap also remains.",
        "protein_rule": "Original detected-normal-IHC exclusion retained, including uncertain staining; uncertain signal is a hold, not verified protein.",
        "age_scope": "Normal panels used by the adult research screen; no pediatric normal panel or age-stratified validation.",
    }
    if neutral_reproductive_immune:
        summary["context_rule"] = (
            "Immune cells sampled in an explicit reproductive or developmental context are neutral for normal exclusion; blood and other tissues remain screened. Neutral observations cannot provide positive reproductive admission."
        )
        summary["historical_uniform_genes"] = int(facts.historical_uniform_selected.sum())
        summary["context_neutral_uniform_gained"] = int(
            (facts.selected & ~facts.historical_uniform_selected).sum()
        )
        summary["screened_cluster_panel"] = sorted(
            set(facts.loc[facts.cell_records.eq(426), "cell_screened_records"].astype(int))
        )
        summary["neutral_cluster_panel"] = sorted(
            set(facts.loc[facts.cell_records.eq(426), "cell_neutral_records"].astype(int))
        )
        if not linked_tiers:
            summary["protected30_sensitivity"] = (
                "Supported protected-tier HPA cluster RNA <30 nTPM instead of the IHC-dependent 5-to-10 ladder; every other gate is unchanged. This does not replace the main tiered scenario or alter reference-review cutoffs."
            )
    if linked_tiers:
        summary.update(
            scope="Shared normal-RNA ladder applied to every gene passing association, protein and coverage gates under either tumor rule; old uniform RNA caps do not restrict the pool.",
            bulk_rule="Uncapped existing tissue/IHC absolute ladder; inclusive ceilings; no legacy selectivity gate.",
            cell_rule="Same tissue/IHC ladder and inclusive boundary as bulk; identified immune cells default to Protected, with no Manageable subset exceptions.",
            cornea_rule="Native CP10K ceiling is the shared per-million ceiling divided by 100. This is nominal policy linkage, not calibrated equivalence to TPM. Original 1 CP10K is reported only as a sensitivity.",
            reference_rule="Review begins at 2x the shared candidate ceiling. This provisional multiplier is not a validated toxicity boundary or an automatic sequence-reference change.",
            assay_deviations=[],
            historical_uniform_maximum_tumor_lost=summary["maximum_tumor_lost"],
            historical_uniform_maximum_tumor_gained=summary["maximum_tumor_gained"],
            maximum_tumor_lost=int(
                (genes.linked_recurrence_selected & ~genes.linked_maximum_tumor_selected).sum()
            ),
            maximum_tumor_gained=int(
                (~genes.linked_recurrence_selected & genes.linked_maximum_tumor_selected).sum()
            ),
        )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    outputs = {p.name: digest(p) for p in out.iterdir() if p.is_file()}
    inputs[str(Path(__file__).relative_to(ROOT))] = digest(__file__)
    for name in [
        "scripts/cta_cell_policy_report.py",
        "oncoref/cell_rna_policy.py",
        "oncoref/cta_normal_screen.py",
        "oncoref/cta_quantitative.py",
        "oncoref/normal_cell_screen.py",
        "oncoref/normal_tissues.py",
        "oncoref/cta_tissues.py",
        "oncoref/data/normal-tissue-taxonomy.yaml",
        "oncoref/data/normal-tissue-mappings.csv",
    ]:
        inputs[name] = digest(ROOT / name)
    (out / "run-manifest.json").write_text(
        json.dumps({"inputs": inputs, "sources": receipts, "outputs": outputs}, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--neutral-reproductive-immune", action="store_true")
    parser.add_argument("--linked-tiers", action="store_true")
    parser.add_argument("--tumor-floor", type=float, default=30)
    args = parser.parse_args()
    build(
        args.output,
        neutral_reproductive_immune=args.neutral_reproductive_immune,
        linked_tiers=args.linked_tiers,
        tumor_floor=args.tumor_floor,
    )
