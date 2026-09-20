#!/usr/bin/env python3
"""Normal-tissue RNA and IHC figures for the current compact-list protein identities."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from cta_report_common import seal_stage, verify_analysis

from oncoref import reference_data
from oncoref.cta_tissues import PERMISSIVE_REPRODUCTIVE_TISSUES
from oncoref.hpa import _read_hpa, hpa_normal_tissue

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

ROOT = Path(__file__).resolve().parents[1]
KEY = "proteoform_key"
ID = "Ensembl_Gene_ID"
LEVELS = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
REPRODUCTIVE = PERMISSIVE_REPRODUCTIVE_TISSUES | {"endometrium 1", "endometrium 2"}
MEMBER_COLORS = ["#127C80", "#C8772B", "#725C9E"]


def ordered_tissues(tissues):
    names = set(tissues)
    priority = [
        "testis",
        "ovary",
        "placenta",
        "epididymis",
        "seminal vesicle",
        "prostate",
        "fallopian tube",
        "endometrium",
        "endometrium 1",
        "endometrium 2",
        "cervix",
        "vagina",
        "breast",
        "lactating breast",
    ]
    reproductive = [t for t in priority if t in names]
    other = sorted(names - set(reproductive))
    assert set(reproductive) == names & REPRODUCTIVE
    return reproductive, other


def ihc_matrix(observations, members, tissues):
    scored = observations[observations.Level.isin(LEVELS)].copy()
    scored["score"] = scored.Level.map(LEVELS)
    return (
        scored.groupby(["Gene", "Tissue"])
        .score.max()
        .unstack()
        .reindex(index=members, columns=tissues)
    )


def sum_member_rna(matrix):
    """Only return a tissue sum when every registered member has a source value."""
    return matrix.sum(axis=0, min_count=len(matrix))


def routine_ihc_tissues(observations):
    scored = observations[observations.Level.isin(LEVELS)]
    genes = scored.Gene.nunique()
    per_tissue = scored.groupby("Tissue").Gene.nunique()
    return set(per_tissue[per_tissue >= 0.5 * genes].index)


def tissue_maximum(values, n_members):
    """Keep rounded-zero status, its bound, and all ties separate from absence."""
    maximum = values.max()
    return {
        "value": maximum,
        "tissues": ";".join(sorted(values.index[values.eq(maximum)])),
        "status": "unavailable"
        if pd.isna(maximum)
        else "reported_zero"
        if maximum == 0
        else "positive_estimate",
        "zero_upper_bound": 0.05 * n_members if maximum == 0 else np.nan,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_proteoform_report_20260917")
    out = parser.parse_args().out.resolve()
    dest = out / "hpa_tissues"
    dest.mkdir(exist_ok=True)
    selected = pd.read_csv(
        out / "selections/min20/prevalence_gt50_transcriptome_p90/proteoforms_ranked.csv"
    )
    strict = pd.read_csv(
        out / "selections/min20/prevalence_gt70_transcriptome_p90/proteoforms_ranked.csv"
    )
    members = pd.read_csv(out / "gene_to_proteoform_mapping.csv")
    members = members[members[KEY].isin(selected[KEY])].copy()
    assert members[ID].is_unique
    verify_analysis(out)
    rna_all, ihc_all = _read_hpa("hpa_rna_consensus", "v23"), hpa_normal_tissue("v23")
    routine = routine_ihc_tissues(ihc_all)
    rna = rna_all[rna_all.Gene.isin(members[ID])].copy()
    ihc = ihc_all[ihc_all.Gene.isin(members[ID])].copy()
    assert not rna.duplicated(["Gene", "Tissue"]).any()
    assert rna.nTPM.ge(0).all()
    rna[KEY] = rna.Gene.map(members.set_index(ID)[KEY])
    ihc[KEY] = ihc.Gene.map(members.set_index(ID)[KEY])
    rna["tissue_group"] = np.where(
        rna.Tissue.isin(REPRODUCTIVE), "Reproductive (including breast)", "Other"
    )
    ihc["tissue_group"] = np.where(
        ihc.Tissue.isin(REPRODUCTIVE), "Reproductive (including breast)", "Other"
    )
    rna.to_csv(dest / "member_rna_ntpm.csv", index=False)
    ihc.to_csv(dest / "member_ihc_observations.csv", index=False)
    links = members.copy()
    links["hpa_tissue_url"] = [
        f"https://v23.proteinatlas.org/{r.Ensembl_Gene_ID}-{r.Symbol}/tissue"
        for r in links.itertuples()
    ]
    links.to_csv(dest / "hpa_gene_links.csv", index=False)
    rep, other = ordered_tissues(rna_all.Tissue.dropna().unique())
    tissues = rep + other
    ihc_rep, ihc_other = ordered_tissues(ihc_all.Tissue.dropna().unique())
    ihc_tissues = ihc_rep + ihc_other
    rna_wide = rna.pivot(index="Gene", columns="Tissue", values="nTPM").reindex(columns=tissues)
    ihc_wide = ihc_matrix(ihc, members[ID], ihc_tissues)
    protein_rna, protein_ihc, summary_rows = [], [], []
    for r in selected.itertuples():
        key = getattr(r, KEY)
        ids = members.loc[members[KEY].eq(key), ID].tolist()
        summed = sum_member_rna(rna_wide.reindex(ids))
        strongest = ihc_wide.reindex(ids).max(axis=0)
        protein_rna.append(summed.rename(key))
        protein_ihc.append(strongest.rename(key))
        sub = ihc[ihc.Gene.isin(ids) & ihc.Level.isin(LEVELS)]
        other_max = tissue_maximum(summed[other], len(ids))
        summary_rows.append(
            {
                KEY: key,
                "Symbol": r.Symbol,
                "protein_length_aa": r.protein_length_aa,
                "in_stricter_list": key in set(strict[KEY]),
                "member_symbols": r.member_symbols,
                "n_rna_tissues_complete": int(summed.notna().sum()),
                "max_reproductive_summed_ntpm": summed[rep].max(),
                "max_reproductive_tissue": ";".join(
                    sorted(summed[rep].index[summed[rep].eq(summed[rep].max())])
                ),
                "max_other_summed_ntpm": other_max["value"],
                "max_other_tissue": other_max["tissues"],
                "max_other_status": other_max["status"],
                "max_other_zero_upper_bound_ntpm": other_max["zero_upper_bound"],
                "n_ihc_tissues_with_scored_data": int(strongest.notna().sum()),
                "n_ihc_routine_tissues_expected": len(routine),
                "n_ihc_routine_tissues_scored": int(
                    strongest.reindex(sorted(routine)).notna().sum()
                ),
                "n_ihc_special_study_tissues_scored": int(
                    strongest.drop(index=sorted(routine)).notna().sum()
                ),
                "n_ihc_tissues_any_detected": int(strongest.gt(0).sum()),
                "ihc_reliability": ";".join(sorted(sub.Reliability.dropna().unique())),
            }
        )
    pr = pd.DataFrame(protein_rna).rename_axis(KEY)
    pi = pd.DataFrame(protein_ihc).rename_axis(KEY)
    pr.to_csv(dest / "protein_summed_rna_ntpm.csv")
    pi.to_csv(dest / "protein_max_ihc_score.csv")
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(dest / "protein_tissue_summary.csv", index=False)
    pd.DataFrame(
        {
            "Tissue": tissues,
            "group": ["Reproductive (including breast)"] * len(rep) + ["Other"] * len(other),
        }
    ).to_csv(dest / "rna_tissue_display_groups.csv", index=False)
    assert pr.shape == (len(selected), len(tissues))
    # Unmeasured members remain unavailable in sums; never fill them with zero.

    np.random.seed(20260920)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.hashsalt": "oncoref-cta-report",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    plot_rows = []

    def save(fig, name, title):
        fig.savefig(
            out / "plots" / f"{name}.png", dpi=165, bbox_inches="tight", metadata={"Date": None}
        )
        fig.savefig(out / "plots" / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
        plt.close(fig)
        plot_rows.append({"name": name, "title": title, "category": "hpa"})

    names = [
        f"{r.Symbol}{'*' if r.in_stricter_list else ''} ({r.protein_length_aa} aa)"
        for r in summary.itertuples()
    ]
    fig, ax = plt.subplots(figsize=(23, 9), layout="constrained")
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#D4D4D4")
    cmap.set_under("white")
    im = ax.imshow(np.log10(pr.to_numpy() + 1), aspect="auto", cmap=cmap, vmin=1e-12)
    ax.axvline(len(rep) - 0.5, color=MEMBER_COLORS[1], lw=3)
    for i in range(len(pr)):
        for j, v in enumerate(pr.iloc[i]):
            txt = "NA" if pd.isna(v) else f"{v:.0f}" if v >= 10 else (f"{v:.1f}" if v > 0 else "0*")
            ax.text(
                j,
                i,
                txt,
                ha="center",
                va="center",
                fontsize=6,
                color="white" if np.log10(v + 1) > im.norm.vmax * 0.58 else "#203746",
            )
    ax.set(
        xticks=range(len(tissues)),
        xticklabels=tissues,
        yticks=range(len(names)),
        yticklabels=names,
        title="Normal-tissue RNA for the compact CTA lists: reproductive tissues on the left, other tissues on the right",
    )
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    colorbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    ticks = [n for n in [0, 1, 10, 100, 1000, 10000] if np.log10(n + 1) <= im.norm.vmax]
    colorbar.set_ticks(np.log10(np.array(ticks) + 1), labels=[str(n) for n in ticks])
    colorbar.set_label("Summed member RNA (nTPM; logarithmic color)")
    fig.supxlabel(
        "Source: HPA v23 RNA consensus via oncoref. Values are sums of member-gene tissue reference nTPM, not measured protein abundance.\n"
        "Left group includes breast; thymus remains among other tissues. * = member of the stricter shortlist. All tissues and rounded estimates are shown.",
        fontsize=10,
    )
    save(
        fig,
        "hpa_rna_tissue_overview",
        "HPA RNA tissue expression: reproductive versus other tissues",
    )

    ihc_cmap = ListedColormap(["#FFFFFF", "#B9DAD9", "#55A5A6", "#146468"])
    ihc_cmap.set_bad("#D4D4D4")
    ihc_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], ihc_cmap.N)
    fig, ax = plt.subplots(figsize=(23, 8.5), layout="constrained")
    im = ax.imshow(pi.to_numpy(), aspect="auto", cmap=ihc_cmap, norm=ihc_norm)
    ax.axvline(len(ihc_rep) - 0.5, color=MEMBER_COLORS[1], lw=3)
    ax.set(
        xticks=range(len(ihc_tissues)),
        xticklabels=ihc_tissues,
        yticks=range(len(names)),
        yticklabels=names,
        title="HPA normal-tissue protein staining: strongest reported observation per protein entry and tissue",
    )
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    colorbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, ticks=[0, 1, 2, 3])
    colorbar.ax.set_yticklabels(["Not detected", "Low", "Medium", "High"])
    fig.supxlabel(
        "HPA v23 IHC via oncoref. Maximum ordinal level across measured member genes and cell types; levels are not added.\n"
        "Gray = no scored observation; white = measured 'Not detected'. Reliability is retained in the observation table. Left of the divider: reproductive tissues, including breast.\n"
        "RNA and IHC have different source tissue vocabularies. * = member of the stricter shortlist.",
        fontsize=10,
    )
    save(fig, "hpa_ihc_tissue_overview", "HPA protein staining: reproductive versus other tissues")

    for r in selected.itertuples():
        key = getattr(r, KEY)
        mm = members[members[KEY].eq(key)]
        ids, syms = mm[ID].tolist(), mm.Symbol.tolist()
        fig = plt.figure(figsize=(20, 9), layout="constrained")
        gs = fig.add_gridspec(2, 2, width_ratios=[len(rep), len(other)], height_ratios=[3.4, 1.2])
        ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        maximum = float(pr.loc[key].max())
        for ax, group, title in [
            (ax1, rep, "Reproductive tissues (including breast)"),
            (ax2, other, "Other tissues"),
        ]:
            bottom = np.zeros(len(group))
            for i, (gene, sym) in enumerate(zip(ids, syms)):
                values = rna_wide.reindex([gene])[group].iloc[0].to_numpy()
                ax.bar(
                    range(len(group)),
                    values,
                    bottom=bottom,
                    color=MEMBER_COLORS[i % len(MEMBER_COLORS)],
                    label=sym,
                    width=0.8,
                )
                bottom += values
            ax.set(
                xticks=range(len(group)),
                xticklabels=group,
                title=title,
                ylim=(0, max(1, maximum) * 2),
            )
            ax.set_yscale("symlog", linthresh=1)
            ax.tick_params(axis="x", rotation=90, labelsize=8)
            ax.grid(axis="y", alpha=0.15)
        ax1.set_ylabel("RNA nTPM (stacked member genes; symlog scale)")
        ax2.legend(title="HPA gene measurement", fontsize=9, title_fontsize=9, loc="upper right")
        ax3 = fig.add_subplot(gs[1, :])
        obs = ihc_wide.reindex(ids)
        if obs.notna().any().any():
            im = ax3.imshow(obs.to_numpy(), aspect="auto", cmap=ihc_cmap, norm=ihc_norm)
            ax3.axvline(len(ihc_rep) - 0.5, color=MEMBER_COLORS[1], lw=3)
            ax3.set(
                xticks=range(len(ihc_tissues)),
                xticklabels=ihc_tissues,
                yticks=range(len(ids)),
                yticklabels=syms,
                title="Protein staining (IHC): strongest scored cell type per tissue",
            )
            ax3.tick_params(axis="x", rotation=90, labelsize=6.5)
            colorbar = fig.colorbar(im, ax=ax3, fraction=0.015, pad=0.01, ticks=[0, 1, 2, 3])
            colorbar.ax.set_yticklabels(["Not detected", "Low", "Medium", "High"], fontsize=7)
        else:
            ax3.axis("off")
            ax3.text(
                0.5,
                0.55,
                "No scored HPA v23 protein-staining observations for these member genes.\nRNA measurements above do not establish protein abundance.",
                ha="center",
                va="center",
                transform=ax3.transAxes,
                fontsize=13,
            )
        row = summary[summary[KEY].eq(key)].iloc[0]
        fig.suptitle(
            f"{r.Symbol} ({r.protein_length_aa} aa): HPA normal-tissue expression",
            fontsize=18,
            weight="bold",
        )
        maximum = (
            f"< {row.max_other_zero_upper_bound_ntpm:g}"
            if row.max_other_status == "reported_zero"
            else f"{row.max_other_summed_ntpm:g}"
        )
        max_tissues = row.max_other_tissue.split(";") if row.max_other_tissue else []
        where = max_tissues[0] if len(max_tissues) == 1 else f"{len(max_tissues)} tied tissues"
        fig.supxlabel(
            f"Source: HPA v23, accessed through oncoref. RNA: {len(tissues)} tissues; identical-protein member genes are stacked. Same RNA axis on both sides.\n"
            f"Largest other-tissue RNA sum: {maximum} nTPM in {where}; zeros are rounded estimates. IHC gray = unavailable, white = measured not detected; IHC values are not summed.",
            fontsize=9,
        )
        save(
            fig,
            f"hpa_protein_{r.Symbol.replace('/', '_')}",
            f"{r.Symbol}: HPA RNA and normal-tissue protein evidence",
        )

    index = pd.read_csv(out / "plot_index.csv")
    pd.concat([index[index.category.ne("hpa")], pd.DataFrame(plot_rows)], ignore_index=True).to_csv(
        out / "plot_index.csv", index=False
    )
    provenance = [
        reference_data.provenance(name, "v23", verify_content=True)
        for name in ["hpa_rna_consensus", "hpa_normal_tissue"]
    ]
    assert all(p["checksum_matches"] for p in provenance)
    manifest = {
        "sources": provenance,
        "n_protein_entries": len(selected),
        "n_member_genes": len(members),
        "n_rna_tissues": len(tissues),
        "n_ihc_source_tissues": len(ihc_tissues),
        "n_figures": len(plot_rows),
        "n_proteins_with_scored_ihc": int(summary.n_ihc_tissues_with_scored_data.gt(0).sum()),
        "rna_aggregation": "sum member-gene consensus nTPM within tissue; NA if any member missing",
        "ihc_aggregation": "maximum scored level across member genes and cell types; no addition; missing remains unavailable",
        "reproductive_definition": "oncoref PERMISSIVE_REPRODUCTIVE_TISSUES, including breast; thymus classified Other",
        "ihc_reliability": sorted(ihc.loc[ihc.Level.isin(LEVELS), "Reliability"].dropna().unique()),
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = [
        "# HPA normal-tissue expression for the compact CTA protein lists",
        "",
        f"{len(selected)} protein identities; {len(members)} member genes; {len(tissues)} RNA tissues.",
        "Reported zero is a rounded estimate below 0.05 nTPM per member, not biological absence. Equal maxima retain all tied tissues.",
        "Figures are redrawn from oncoref's cached HPA v23 source tables, rather than screenshots of the HPA website.",
        "RNA bars stack the member-gene nTPM measurements; their total is a tissue-reference RNA sum, not protein abundance or patient-specific coexpression. The overview uses logarithmic color; individual bars use a symlog axis with a linear region from 0 to 1 nTPM.",
        "Reproductive tissues appear first and include breast. Thymus is shown among other tissues; none of the 50 RNA tissues is hidden. This display grouping differs from the older restriction annotation, which excludes thymus from its somatic maximum.",
        "IHC protein staining is shown separately and retains source tissue labels. For an overview cell, take the strongest scored observation across member genes/cell types. Individual pages retain the gene rows. Levels are not added. Missing observations are gray; measured Not detected is white. The summary records scored IHC counts and reliability; absence of a score is not evidence of absence.",
        "The earlier shortlist normal-tissue annotation is a maximum single-member somatic nTPM; this new RNA view sums members and also displays thymus and breast explicitly. They answer different summary questions.",
        "",
        "## HPA source pages",
        "",
    ]
    for r in selected.itertuples():
        ll = links[links[KEY].eq(getattr(r, KEY))]
        lines.append(
            f"- {r.Symbol}: "
            + ", ".join(f"[{q.Symbol}]({q.hpa_tissue_url})" for q in ll.itertuples())
        )
    lines += [
        "",
        "## Data",
        "",
        "- protein_tissue_summary.csv: lengths, maxima and IHC availability.",
        "- protein_summed_rna_ntpm.csv: the protein-by-tissue RNA matrix.",
        "- member_rna_ntpm.csv: all original member-gene tissue values.",
        "- member_ihc_observations.csv: tissue/cell type, level and reliability.",
        "- protein_max_ihc_score.csv: 0 Not detected, 1 Low, 2 Medium, 3 High; missing blank.",
        "- hpa_gene_links.csv: direct versioned HPA pages.",
        "- manifest.json: source URLs, version and verified checksums.",
    ]
    (dest / "report.md").write_text("\n".join(lines) + "\n")
    seal_stage(
        out,
        "hpa",
        [Path(__file__), out / "analysis_receipt.json"],
        [
            *sorted(dest.glob("*.csv")),
            dest / "manifest.json",
            *[
                out / "plots" / f"{row['name']}.{ext}"
                for row in plot_rows
                for ext in ("png", "svg")
            ],
        ],
    )
    print(summary.to_string(index=False))
    print(f"Created {len(plot_rows)} HPA figures")


if __name__ == "__main__":
    main()
