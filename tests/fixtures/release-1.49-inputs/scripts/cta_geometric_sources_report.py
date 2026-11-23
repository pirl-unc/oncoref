#!/usr/bin/env python3
"""Explain the geometric funnel's sources, tumor evidence gaps and RNA coverage."""

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

import matplotlib.pyplot as plt  # noqa: E402
import plot_cta_geometric as style  # noqa: E402
from cta_report_common import sha256 as digest  # noqa: E402
from cta_report_common import write_csv as write  # noqa: E402
from cta_selection_figures import LABELS  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib_venn import venn2  # noqa: E402

from oncoref import cta_provenance  # noqa: E402
from oncoref.cta_sources import publication_sources  # noqa: E402


def plot_routes(routes: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.8))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.77, bottom=0.18, wspace=0.32)
    style.title(
        fig,
        "Two quantitative nomination routes",
        "HPA reproductive RNA and Tyser embryonic RNA are alternative entry routes; every later filter still applies.",
    )
    for ax, scope, heading in zip(
        axes,
        ["association", "selected"],
        ["Before tumor and normal filters", "After the complete funnel"],
        strict=True,
    ):
        r = routes[routes.scope.eq(scope)].set_index("route").genes
        v = venn2(
            subsets=(int(r["HPA only"]), int(r["Tyser only"]), int(r["Both"])),
            set_labels=("HPA reproductive", "Tyser embryonic"),
            set_colors=(style.TEAL, "#8374AE"),
            alpha=0.65,
            ax=ax,
        )
        for label in v.subset_labels:
            if label is not None:
                label.set_fontsize(16)
                label.set_fontweight("bold")
                label.set_text(f"{int(label.get_text()):,}")
        for label in v.set_labels:
            label.set_fontsize(11)
        total = int(r.sum())
        ax.set_title(f"{heading}\n{total:,} genes in the union", fontsize=13, weight="bold", pad=12)
    style.footer(
        fig,
        "HPA v23: testis/ovary/placenta ≥5 nTPM and greater than nonreproductive bulk. Tyser E-MTAB-9388: supported lineage mean ≥5 TPM.\nTyser is one early embryo, not a curated oncofetal-antigen list. Circle areas are comparable only within each panel.",
    )
    style.save(fig, out, "quantitative-source-overlap")


def plot_papers(sources: pd.DataFrame, overlaps: pd.DataFrame, out: Path) -> None:
    order = sources.sort_values(["selected_genes", "mapped_genes"], ascending=False).doi.tolist()
    labels = [LABELS.get(doi, doi) for doi in order]
    counts = overlaps.pivot(index="source_a", columns="source_b", values="selected_overlap").loc[
        order, order
    ]
    jac = overlaps.pivot(index="source_a", columns="source_b", values="selected_jaccard").loc[
        order, order
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 9), gridspec_kw={"width_ratios": [0.90, 1.4]})
    fig.subplots_adjust(left=0.12, right=0.96, top=0.78, bottom=0.25, wspace=0.55)
    style.title(
        fig,
        "How the selected genes overlap published lists",
        "306 of 319 candidates occur in at least one imported publication list; 13 have no mapped membership in this registry.",
    )
    s = sources.set_index("doi").loc[order]
    y = np.arange(len(order))
    axes[0].barh(y, s.mapped_genes, color=style.PALE, height=0.66, label="All mapped genes")
    axes[0].barh(y, s.selected_genes, color=style.TEAL, height=0.38, label="In selected 319")
    axes[0].set_yticks(y, labels, fontsize=9)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 1340)
    axes[0].set_xlabel("Genes per paper · overlap allowed")
    axes[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2, fontsize=9)
    axes[0].set_title("Selected / all mapped", loc="left", weight="bold", fontsize=12, pad=15)
    for i, r in enumerate(s.itertuples()):
        axes[0].text(
            1330, i, f"{r.selected_genes} / {r.mapped_genes:,}", va="center", ha="right", fontsize=9
        )
    axes[1].imshow(jac.to_numpy(), cmap="YlGnBu", norm=Normalize(0, 1), aspect="equal")
    axes[1].set_xticks(y, labels, fontsize=8, rotation=55, ha="right")
    axes[1].set_yticks(y, labels, fontsize=8)
    axes[1].tick_params(length=0)
    axes[1].set_title(
        "Shared selected genes · diagonal = paper total",
        loc="left",
        weight="bold",
        fontsize=11,
        pad=15,
    )
    for i, j in np.ndindex(counts.shape):
        axes[1].text(
            j,
            i,
            str(counts.iloc[i, j]),
            ha="center",
            va="center",
            fontsize=7.5,
            color="white" if jac.iloc[i, j] > 0.55 else style.INK,
        )
    style.footer(
        fig,
        "Paper memberships annotate provenance; they are not admission requirements for the quantitative funnel. Tables from the same paper are deduplicated.\nMatrix color = Jaccard overlap; numbers = gene counts. Targeted reports are included. Boel/BAGE has no resolved eligible locus in this registry.\nUnmapped, noncoding and ambiguous entries do not silently become protein-coding candidates. Registry coverage is not all published literature.",
    )
    style.save(fig, out, "publication-source-overlap")


def plot_gaps(tumor: pd.DataFrame, coverage: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.8))
    fig.subplots_adjust(left=0.19, right=0.96, top=0.74, bottom=0.20, wspace=0.90)
    style.title(
        fig,
        "Low tumor expression is different from missing data",
        "The two funnel stages remove 112 and 122 genes for different reasons; neither stage proves biological absence.",
    )
    specs = [
        (
            axes[0],
            ["Measured, max <10 pTPM", "No usable tumor measurement"],
            [
                int(tumor.tumor_reason.eq("Measured maximum below 10 pTPM").sum()),
                int(tumor.tumor_reason.eq("No usable mapped tumor measurement").sum()),
            ],
            "112 fail tumor admission",
        ),
        (
            axes[1],
            coverage.missing_panels.tolist(),
            coverage.genes.tolist(),
            "122 fail normal-RNA coverage",
        ),
    ]
    for ax, labels, values, heading in specs:
        ax.barh(
            range(len(labels)),
            values,
            color=[style.TEAL] + [style.MUTED] * (len(labels) - 1),
            height=0.55,
        )
        ax.set_yticks(range(len(labels)), labels, fontsize=10)
        ax.invert_yaxis()
        ax.set_xlim(0, max(values) * 1.25)
        ax.set_title(heading, loc="left", fontsize=12, weight="bold", pad=18)
        ax.set_xlabel("Genes")
        ax.tick_params(axis="y", length=0)
        for i, n in enumerate(values):
            ax.text(n + max(values) * 0.03, i, str(n), va="center", fontsize=13, weight="bold")
    style.footer(
        fig,
        "All 48 measured tumor failures have positive maxima (0.5288–9.9904 pTPM); none has a measured zero maximum.\nCoverage groups are mutually exclusive. Independently, 114 lack corneal coverage, 31 lack bulk coverage and 18 lack DICE coverage.",
    )
    style.save(fig, out, "tumor-and-coverage-gaps")


def build(audit: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    out = audit / "source-review"
    out.mkdir(exist_ok=True)
    g = pd.read_csv(audit / "gene-comparison.csv.gz", low_memory=False)
    manifest = json.loads((audit / "run-manifest.json").read_text())
    assert digest(audit / "gene-comparison.csv.gz") == manifest["outputs"]["gene-comparison.csv.gz"]
    chosen = g.two_axis_doubling_maximum10_selected
    routes = []
    for scope, keep in [("association", g.positive_context), ("selected", chosen)]:
        for label, mask in [
            ("HPA only", g.reproductive_route & ~g.developmental_route),
            ("Tyser only", ~g.reproductive_route & g.developmental_route),
            ("Both", g.reproductive_route & g.developmental_route),
        ]:
            routes.append({"scope": scope, "route": label, "genes": int((keep & mask).sum())})
    routes = pd.DataFrame(routes)
    write(routes, out / "quantitative-source-overlaps.csv")
    memberships = cta_provenance.paper_membership()
    sources = publication_sources()
    coding = memberships[
        memberships.protein_candidate_eligible & memberships.Ensembl_Gene_ID.isin(g.gene_id)
    ].drop_duplicates(["doi", "Ensembl_Gene_ID"])
    selected_ids = set(g.loc[chosen, "gene_id"])
    sets = {
        doi: set(coding.loc[coding.doi.eq(doi), "Ensembl_Gene_ID"]) for doi in sources.doi.unique()
    }
    rows = []
    for doi, ids in sets.items():
        rows.append(
            {
                "doi": doi,
                "paper": LABELS.get(doi, doi),
                "mapped_genes": len(ids),
                "selected_genes": len(ids & selected_ids),
                "source_tags": ";".join(sorted(sources.loc[sources.doi.eq(doi), "source_tag"])),
                "evidence_scope": " | ".join(
                    sources.loc[sources.doi.eq(doi), "evidence_scope"].drop_duplicates()
                ),
            }
        )
    paper_counts = pd.DataFrame(rows)
    overlaps = []
    for a, ids_a in sets.items():
        for b, ids_b in sets.items():
            sa, sb = ids_a & selected_ids, ids_b & selected_ids
            overlaps.append(
                {
                    "source_a": a,
                    "source_b": b,
                    "mapped_overlap": len(ids_a & ids_b),
                    "selected_overlap": len(sa & sb),
                    "selected_jaccard": len(sa & sb) / len(sa | sb) if sa | sb else 0,
                }
            )
    overlaps = pd.DataFrame(overlaps)
    assert len(set().union(*sets.values()) & selected_ids) == 306
    write(paper_counts, out / "publication-source-counts.csv")
    write(overlaps, out / "publication-source-overlaps.csv")
    write(
        coding[coding.Ensembl_Gene_ID.isin(selected_ids)],
        out / "selected-publication-membership.csv.gz",
    )
    write(
        g[chosen & ~g.gene_id.isin(set().union(*sets.values()))],
        out / "selected-without-mapped-publication.csv",
    )
    failed = g[g.positive_context & ~g.tumor_maximum_max_ptpm.ge(10)].copy()
    failed["tumor_reason"] = np.where(
        failed.tumor_maximum_max_ptpm.isna(),
        "No usable mapped tumor measurement",
        "Measured maximum below 10 pTPM",
    )
    failed["display_gene"] = failed.symbol.fillna(failed.gene_id)
    tumor_cols = [
        "gene_id",
        "display_gene",
        "tumor_reason",
        "tumor_maximum_max_ptpm",
        "tumor_maximum_cancer_code",
        "tumor_maximum_cohort",
        "tumor_cohorts_measured",
        "reproductive_route",
        "developmental_route",
    ]
    write(failed[tumor_cols], out / "tumor-admission-failures.csv")
    write(
        failed[failed.tumor_maximum_max_ptpm.notna()][tumor_cols].sort_values(
            "tumor_maximum_max_ptpm"
        ),
        out / "tumor-measured-below-10.csv",
    )
    write(failed[failed.tumor_maximum_max_ptpm.isna()][tumor_cols], out / "tumor-unmeasured.csv")
    assert len(failed) == 112 and failed.tumor_maximum_max_ptpm.eq(0).sum() == 0
    coverage = g[
        g.positive_context
        & g.tumor_maximum_max_ptpm.ge(10)
        & ~g.ihc_detected
        & ~g.rna_coverage_complete
    ].copy()
    missing = {
        "Bulk": coverage.normal_bulk_tissues.ne(50)
        | coverage.reproductive_ntpm.isna()
        | coverage.adult_bulk_max_ntpm.isna(),
        "DICE": coverage.dice_sources.ne(15) | coverage.dice_third_max_tpm.isna(),
        "Cornea": coverage.cornea_records.ne(84),
    }
    for key, values in missing.items():
        coverage[f"missing_{key}"] = values
    coverage["missing_panels"] = [
        " + ".join(k for k in missing if coverage.loc[i, f"missing_{k}"]) for i in coverage.index
    ]
    assert coverage.missing_panels.ne("").all() and len(coverage) == 122
    coverage_groups = (
        coverage.groupby("missing_panels")
        .size()
        .rename("genes")
        .reset_index()
        .sort_values("genes", ascending=False)
    )
    write(coverage, out / "normal-rna-coverage-failures.csv")
    write(coverage_groups, out / "coverage-failure-combinations.csv")
    write(
        g[chosen]
        .groupby(["Expression category", "Expression detail"])
        .size()
        .rename("genes")
        .reset_index(),
        out / "selected-historical-labels.csv",
    )
    write(g[chosen & g["Expression category"].eq("Unsafe")], out / "selected-historical-unsafe.csv")
    membership = g[
        [
            "gene_id",
            "symbol",
            "reproductive_route",
            "developmental_route",
            "positive_context",
            "two_axis_doubling_maximum10_selected",
        ]
    ].copy()
    membership["mapped_publication_count"] = (
        membership.gene_id.map(coding.groupby("Ensembl_Gene_ID").doi.nunique())
        .fillna(0)
        .astype(int)
    )
    write(membership, out / "source-route-membership.csv.gz")
    style.setup()
    plot_routes(routes, out)
    plot_papers(paper_counts, overlaps, out)
    plot_gaps(failed, coverage_groups, out)
    inputs = {
        str((audit / "gene-comparison.csv.gz").relative_to(ROOT)): digest(
            audit / "gene-comparison.csv.gz"
        )
    }
    for name in [
        "scripts/cta_geometric_sources_report.py",
        "scripts/plot_cta_geometric.py",
        "scripts/cta_selection_figures.py",
        "oncoref/cta_provenance.py",
        "oncoref/cta_sources.py",
        "oncoref/data/cta-publication-membership.csv",
        "oncoref/data/cta-publication-sources.csv",
    ]:
        inputs[name] = digest(ROOT / name)
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": inputs,
                "outputs": {
                    p.name: digest(p)
                    for p in out.iterdir()
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(
        paper_counts[["paper", "mapped_genes", "selected_genes"]]
        .sort_values("selected_genes", ascending=False)
        .to_string(index=False)
    )
    print(routes.to_string(index=False))
    print(coverage_groups.to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("audit", type=Path)
    build(p.parse_args().audit.resolve())
