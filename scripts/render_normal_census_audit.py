#!/usr/bin/env python3
"""Render the normal-atlas inventory, QC requirements and threshold sensitivity."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--siletti", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    inventory = json.loads((args.inventory / "manifest.json").read_text())
    run = json.loads((args.counts / "manifest.json").read_text())
    populations = pd.read_parquet(args.counts / "populations.parquet")
    ages = populations.groupby("age_scope").n_cells.sum()
    tissues = (
        populations.groupby(["tissue_general", "age_scope"]).n_cells.sum().unstack(fill_value=0)
    )
    tissues.to_csv(args.output / "tissue-inventory.csv")
    sources = (
        populations.groupby(
            ["collection_id", "collection_name", "collection_doi", "age_scope"], dropna=False
        )
        .n_cells.sum()
        .unstack(fill_value=0)
        .reset_index()
    )
    sources["source_url"] = "https://cellxgene.cziscience.com/collections/" + sources.collection_id
    sources["paper_url"] = sources.collection_doi.fillna("").map(
        lambda x: "https://doi.org/" + x if x else ""
    )
    for field in [
        "donor_health_review",
        "ambient_RNA_review",
        "doublet_review",
        "alignment_review",
        "feature_filter_review",
    ]:
        sources[field] = "pending"
    sources.to_csv(args.output / "source-registry.csv", index=False)
    pd.read_csv(args.inventory / "datasets.csv").query(
        "dataset_id in @populations.dataset_id"
    ).to_csv(args.output / "dataset-versions.csv", index=False)
    populations.groupby(
        ["development_stage", "development_stage_ontology_term_id", "age_scope"]
    ).n_cells.sum().to_csv(args.output / "age-mapping.csv")

    base = ROOT / "docs/audits/plac1-brain-20261110"
    catalog = pd.read_csv(base / "cta-comparison/adult-granule-current-ctas.csv")
    support = pd.read_csv(base / "count-qc/count-support.csv")
    members = catalog[["entity_id", "protein", "member_gene_ids"]].copy()
    members["gene_id"] = members.member_gene_ids.map(json.loads)
    joined = members.explode("gene_id").merge(support, on="gene_id", how="inner")
    # Apply any-member expression separately within each donor, never sum paralog rates.
    joined = joined.sort_values("pct_positive", ascending=False).drop_duplicates(
        ["entity_id", "scope"]
    )
    sensitivity = []
    for threshold in [1, 2, 5]:
        positive = joined[joined.pct_positive.ge(threshold) & joined.n_positive.ge(10)]
        donors = positive[positive.scope.ne("all")]
        sensitivity.append(
            {
                "threshold_percent": threshold,
                "pooled": positive[positive.scope.eq("all")].entity_id.nunique(),
                "at_least_one_donor": donors.entity_id.nunique(),
                "both_donors": int(donors.groupby("entity_id").scope.nunique().ge(2).sum()),
            }
        )
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity.to_csv(args.output / "cerebellum-threshold-sensitivity.csv", index=False)
    joined.to_csv(
        args.output / "cerebellum-donor-support.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    top = tissues.adult.nlargest(10).sort_values()
    axes[0].barh(top.index, top / 1e6, color="#236b8e")
    axes[0].set_xlabel("Adult normal-labelled observations (millions)")
    axes[0].set_title("Large inventory; source QC remains separate")
    x = np.arange(len(sensitivity))
    for offset, field, label, color in [
        (-0.25, "pooled", "Donors pooled", "#c4813d"),
        (0, "at_least_one_donor", "Either donor", "#236b8e"),
        (0.25, "both_donors", "Both donors", "#59814e"),
    ]:
        bars = axes[1].bar(x + offset, sensitivity[field], width=0.24, label=label, color=color)
        axes[1].bar_label(bars, padding=2, fontsize=9)
    axes[1].set_xticks(x, [f"{n}%" for n in sensitivity.threshold_percent])
    axes[1].set_xlabel("RNA detection cutoff (any positive count)")
    axes[1].set_ylabel("Current CTAs above cutoff")
    axes[1].set_title("Cerebellum: calls depend on the counting rule")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(args.output / "inventory-and-thresholds.png", dpi=180)
    fig.savefig(args.output / "inventory-and-thresholds.svg")
    plt.close(fig)
    blocks = sorted(args.counts.glob("counts-*.parquet"))
    n_features = len(pd.read_csv(args.counts / "selected-features.csv"))
    expected = (n_features + run["block_size"] - 1) // run["block_size"]
    stats = {
        "inventory_observations": inventory["n_cells"],
        "adult_observations": int(ages.adult),
        "adult_brain_observations": int(tissues.loc["brain", "adult"]),
        "source_collections": len(sources),
        "source_datasets": populations.dataset_id.nunique(),
        "expression_blocks_complete": len(blocks),
        "expression_blocks_expected": expected,
        "expression_proteins_selected": run["n_proteins"],
        "expression_genes_available": n_features,
        "classification_changed": False,
    }
    (args.output / "summary.json").write_text(json.dumps(stats, indent=2) + "\n")
    rows = "\n".join(
        f"| {r.threshold_percent}% | {r.pooled} | {r.at_least_one_donor} | {r.both_donors} |"
        for r in sensitivity.itertuples()
    )
    text = f"""# Normal-cell atlas reassessment

Audit: 10 November 2026. **The normal-tissue inventory contains {inventory["n_cells"]:,} primary cells/nuclei from {len(sources)} collections and {populations.dataset_id.nunique()} datasets. These are observations retained by metadata filters, not cells independently certified as healthy or free of technical noise.**

The expression build selects **{run["n_proteins"]:,} proteins**: every current non-Unsafe record, every protein excluded only by the old rare-cell rule, and the familiar-CTA comparison panel. This includes PLAC1 and all 638 current CTAs. {n_features:,} exact Ensembl gene IDs are available; {len(run["unavailable_genes"])} requested IDs are absent from this Census release. Other Unsafe proteins retain their independent exclusions; they are not newly reassessed by this expression extraction.

Expression blocks completed when this page was rendered: **{len(blocks)}/{expected}**. The inventory and Sepp sensitivity analysis are complete. The broad expression build and source-level QC are separate work. **No catalog classifications have changed in this audit.**

## Inventory

| Age scope | Observations |
|---|---:|
| Adult | {int(ages.adult):,} |
| Childhood | {int(ages.child):,} |
| Prenatal | {int(ages.prenatal):,} |
| Age unresolved | {int(ages.unresolved):,} |

Adult brain alone contributes {int(tissues.loc["brain", "adult"]):,} observations. [Tissues](tissue-inventory.csv), [all source collections and papers](source-registry.csv), [dataset versions](dataset-versions.csv), [age mapping](age-mapping.csv).

The build pins CELLxGENE Census **{run["census_release"]}**, the stable release returned by the service on the audit date. It retains primary observations with the exact disease label `normal` and `tissue_type=tissue`; this excludes labelled organoids and cell cultures. Compound or uncertain disease labels are not silently accepted. Adult status comes from donor-stage metadata, never a cell-type name. Childhood, prenatal and ambiguous stages stay separate. A normal label does not exclude unrecorded illness, treatment, tissue pathology or postmortem effects; these require source review.

For this audit, adult means a stated age or age range entirely **18 years or older**. HsapDv's broad adult, prime-adult and young-adult terms start at **15**, so those labels remain unresolved until the original donor ages are recovered. They are not admitted to the confirmed-adult count. [Ontology definitions](adult-age-ontology.json), [HsapDv ontology](https://www.ebi.ac.uk/ols4/ontologies/hsapdv).

`is_primary_data` removes duplicates identified by Census curation. Dataset partitions and repeated libraries are not independent biological replicates. Donors are scoped to their collection; unknown or pooled donor identities cannot count as independent donors. Cross-collection donor overlap still needs review. Census and HPA summaries of the same study count as one source, and the current Sepp H5AD is a version comparison, not an independent replication of Sepp in Census.

Sources: [Census releases](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html), [Census schema](https://github.com/chanzuckerberg/cellxgene-census/blob/main/docs/cellxgene_census_schema.md), [incremental raw-count API](https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_compute_over_X.html).

## What changes at 5% per donor?

Among current CTAs measured in the two adult Sepp granule-cell donors:

| Detection cutoff | Donors pooled | At least one donor | Both donors |
|---|---:|---:|---:|
{rows}

These are sensitivity counts, not revised classifications. They require at least ten positive nuclei in the relevant comparison; the smaller donor has 2,731 granule nuclei and therefore does not independently satisfy the old 3,000-cell support rule. Missing members are not negative. The 5% single-donor hit is ENSG00000279983. [Threshold table](cerebellum-threshold-sensitivity.csv), [protein and donor results](cerebellum-donor-support.csv.gz).

PLAC1 has 514/16,833 positive adult granule nuclei; 506 contain one count, eight contain two. Its detection rates are 3.25% and 2.01% in the two donors. In library SN233 from the 42-year-old donor, detection is 2.08% in the lowest-depth quarter and 5.65% in the highest-depth quarter. This association does not prove contamination: genuine weak RNA also becomes easier to detect at greater depth. It does show why a raw percentage is not a portable safety threshold. [Count support](../plac1-brain-20261110/count-qc/count-support.csv) · [Depth comparison](../plac1-brain-20261110/count-qc/depth-sensitivity.csv.gz).

The denominator is the sampled named cell population, not the entire organ. Five percent of a specialized population is not five percent organ loss. The fraction of nuclei with captured RNA also does not identify the fraction with recognizable protein or peptide–HLA.

![Inventory and threshold sensitivity](inventory-and-thresholds.png)

## Evidence required before a rare-cell finding changes a category

| Check | Required interpretation |
|---|---|
| Source identity | Pin study, dataset revision, reference annotation and count layer. Do not combine inconsistent source versions. |
| Adult normal tissue | Verify donor age, health, tissue origin and treatment; retain other ages and states separately. |
| Cell QC | Record source cell-calling, damaged-cell and doublet filters; inspect depth, detected genes and appropriate mitochondrial/nuclear metrics within each library and cell type. A universal QC cutoff can remove real rare populations. |
| Count support | Report positive cells, total molecules, one-count and multi-count support per donor. Keep UMI protocols, full-length read counts and fractional corrected values distinct. |
| Background | Compare uncorrected and corrected counts where available; use empty droplets or validated negative controls. A larger percentage alone does not establish expression above background. |
| Alignment and identity | Inspect unique mapping, strandedness, exonic versus intronic support, overlapping genes and paralogs for decisive sparse findings. A gene-count matrix cannot provide these checks. |
| Replication | Look for a localized signal across donors and independent studies. One strong, verified donor is not cancelled by many negative donors. One weak library is not enough for an automatic exclusion. |
| Missingness | Preserve genes missing or censored by upstream filtering. Census presence=False does not distinguish unmeasured genes from genes with no positive counts; verify the original feature list before claiming a measured negative. |
| Protein evidence | Keep IHC, selective MS, spatial RNA and receptor-relevant protein/peptide evidence linked to the exact gene and compartment. RNA replication is not protein confirmation. |

Use **5% as a reported sensitivity cutoff**, alongside 1% and 2%, until calibration against source-specific background and orthogonal evidence supports an operational rule. Do not turn “below 5%” into a safety allowance. Nor should “above 1%” automatically establish Unsafe when the evidence consists of unresolved single-molecule signals.

The intended decision flow stays within the existing categories: credible normal evidence is assessed under its tissue tier; unresolved evidence that could change admission belongs in **Other:Review**; a demonstrated artifact can be dismissed with its evidence; and a passing normal-expression screen permits the existing cancer/reproductive admission rules. This audit supplies the data and QC record before any catalog rebuild adopts that change.

Methods supporting these checks: [SoupX](https://doi.org/10.1093/gigascience/giaa151), [sequencing depth and detection](https://doi.org/10.1038/s41467-020-14482-y), [biological replication in single-cell analysis](https://doi.org/10.1038/s41467-021-25960-2). None of these papers establishes a 1%, 2% or 5% antigen-toxicity threshold.

## Cerebellum count provenance

The authors' [workflow](https://gitlab.com/kaessmannlab/mammalian-cerebellum/-/blob/main/README.org) explicitly counts UMIs with Cell Ranger against mature-mRNA and pre-mRNA references. One count therefore represents a deduplicated molecule assignment, not necessarily one sequencing read. Their [available preprocessing script](https://gitlab.com/kaessmannlab/mammalian-cerebellum/-/blob/main/scripts/R/001_batch_preprocessing.R) selects nuclei by intronic fraction, uses Scrublet, and removes genes detected in fewer than 25 cells within a library. Its merged matrices can therefore contain censored zeros. The current H5AD does not retain the full per-library gene filtering mask or separate exonic/intronic layers. The exact export and background-correction provenance still needs reconciliation with the paper's final workflow.

Those issues warrant checking; they do not establish that the CTA signals are alignment errors. Some reproductive-gene expression in brain is real—the Sepp study independently examined examples such as ZP2—but that does not validate the specific sparse hits in this audit. [Sepp et al.](https://doi.org/10.1038/s41586-023-06884-x).

## Reproduction

Use an isolated environment with `cellxgene-census`, `pyarrow`, `h5py`, `scipy`, `pandas` and `matplotlib`. The shared Python environment was not modified. Large matrices, observation metadata and count blocks are cached under `tmp/normal-atlas`; compact reports are retained here.

```sh
python scripts/inventory_normal_census.py --release {run["census_release"]} --output tmp/normal-atlas/census-{run["census_release"]}
python scripts/audit_normal_census.py --inventory tmp/normal-atlas/census-{run["census_release"]} --output tmp/normal-atlas/counts-{run["census_release"]}
python scripts/render_normal_census_audit.py --inventory tmp/normal-atlas/census-{run["census_release"]} --counts tmp/normal-atlas/counts-{run["census_release"]} --output docs/audits/normal-cell-atlas-20261110
```
"""
    if args.siletti:
        proteins = pd.read_csv(args.siletti / "granule-protein-detection.csv.gz")
        genes = pd.read_csv(args.siletti / "granule-gene-detection.csv.gz")
        sepp_genes = support[support.scope.eq("all")]
        all_genes = genes[genes.scope.eq("all")]
        quoted = [
            "CBY3",
            "PLAC1",
            "DMRT1",
            "CAGE1",
            "TMEM269",
            "NUTM1",
            "SYCP1",
            "MAGEB10",
            "PRAME",
        ]
        comparison = sepp_genes[sepp_genes.symbol.isin(quoted)].merge(
            all_genes,
            on=["gene_id", "symbol"],
            suffixes=("_sepp", "_siletti"),
            validate="one_to_one",
        )
        comparison = comparison.set_index("symbol").loc[quoted].reset_index()
        comparison.to_csv(args.output / "independent-study-comparison.csv", index=False)
        old_hits = pd.read_csv(base / "cta-comparison/current-cta-granule-hits.csv")
        independent = proteins[proteins.scope.eq("all")]
        hit_comparison = old_hits[["entity_id", "protein", "pct_positive"]].merge(
            independent.drop(columns="protein"),
            on="entity_id",
            suffixes=("_sepp", "_siletti"),
            validate="one_to_one",
        )
        hit_comparison.to_csv(args.output / "original-46-independent-study.csv", index=False)
        ctas = proteins[proteins.category.eq("CTA")]
        ctas.to_csv(
            args.output / "independent-study-all-ctas.csv.gz",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
        for filename in [
            "manifest.json",
            "population-qc.csv",
            "granule-gene-detection.csv.gz",
            "granule-protein-detection.csv.gz",
        ]:
            shutil.copyfile(args.siletti / filename, args.output / ("siletti-" + filename))
        measured = hit_comparison.coverage.ne("unmeasured")
        rows = "\n".join(
            f"| {r.symbol} | {r.pct_positive_sepp:.3f}% | {r.pct_positive_siletti:.3f}% | {r.n_positive_siletti:,} | {r.n_ge2_siletti:,} |"
            for r in comparison.itertuples()
        )
        independent_text = f"""## Independent adult-cerebellum replication

The three disjoint cerebellar dissections from [Siletti et al.](https://doi.org/10.1126/science.add7046) contain **151,947 unique nuclei**, including **88,302 granule nuclei** from three adults aged 29, 42 and 50. Cell IDs were checked across files; all retained observations are primary and labelled normal. The audit evaluates **all 19,893 catalog protein records**, using exact gene IDs and reporting unavailable members separately. The PLAC1 totals were independently reproduced with AnnData slicing.

| Protein | Sepp: 16,833 granule nuclei | Siletti: 88,302 granule nuclei | Siletti positive nuclei | Siletti nuclei with ≥2 counts |
|---|---:|---:|---:|---:|
{rows}

Among the original 46 CTA hits, **{int(measured.sum())} have an available member feature**, **{int(hit_comparison.n_positive.gt(0).sum())} have positive counts**, and **{int(hit_comparison.pct_positive_siletti.ge(1).sum())} reaches 1%** in the independent pooled granule population: TMEM269. Five records have no measured member. These results corroborate sparse RNA assignments, not the original detection frequencies or protein production. Shared technical artifacts can recur across studies.

Across all 638 current CTAs, the independent atlas provides complete member features for 573, partial coverage for 21 and no member feature for 44. Three CTAs reach 1% pooled detection: MANSC4 (1.80%), TMEM269 (1.04%) and ENSG00000258691 (1.02%). None reaches 5% pooled detection, and none reaches 1% if a positive nucleus must contain at least two counts. These alternate cutoffs remain descriptive comparisons.

PLAC1 is detected in each independent donor: 99/16,112 (0.61%), 173/40,112 (0.43%) and 95/32,078 (0.30%). MAGEA4, CT83 and XAGE1A/B have available features but zero granule counts in these files. That is a measured-matrix result, not proof of biological absence: the original reference excludes some overlapping genes/transcripts, and upstream processing remains relevant. [Siletti reference and counting workflow](https://github.com/linnarsson-lab/adult-human-brain).

[All original 46 comparisons](original-46-independent-study.csv) · [All current CTAs by donor](independent-study-all-ctas.csv.gz) · [All catalog proteins by donor](siletti-granule-protein-detection.csv.gz) · [Gene counts](siletti-granule-gene-detection.csv.gz) · [Population QC](siletti-population-qc.csv) · [Source versions and hashes](siletti-manifest.json).

"""
        text = text.replace(
            "## Evidence required before", independent_text + "## Evidence required before"
        )
        stats["independent_cerebellum_unique_nuclei"] = 151947
        stats["independent_granule_nuclei"] = 88302
        stats["independent_catalog_proteins_audited"] = 19893
        (args.output / "summary.json").write_text(json.dumps(stats, indent=2) + "\n")
    (args.output / "index.md").write_text(text)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
