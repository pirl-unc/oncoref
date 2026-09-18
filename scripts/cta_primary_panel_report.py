#!/usr/bin/env python3
"""A summary-first presentation of the existing 13-protein CTA shortlist."""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/cta_proteoform_report_20260917"
KEY = "proteoform_key"
SETTINGS = [(50, 70), (70, 70), (50, 90), (70, 90)]
COHORTS = [
    "LUAD",
    "LUSC",
    "COAD",
    "READ",
    "LIHC",
    "BRCA",
    "STAD",
    "PAAD",
    "ESCA",
    "HNSC",
    "PRAD",
    "CESC",
]
NAMES = [
    "Lung adeno.",
    "Lung squamous",
    "Colon",
    "Rectum",
    "Liver",
    "Breast",
    "Stomach",
    "Pancreas",
    "Esophagus",
    "Head/neck",
    "Prostate",
    "Cervix",
]


def prepare(out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from plot_cta_proteoform_report import label_coverage_points, protein_gene_legend

    dest = out / "primary_panel"
    dest.mkdir(exist_ok=True)
    lists = {
        (f, p): pd.read_csv(
            out / f"selections/min20/prevalence_gt{f}_transcriptome_p{p}/proteoforms_ranked.csv"
        )
        for f, p in SETTINGS
    }
    primary = lists[50, 90]
    primary.to_csv(dest / "primary_proteoforms.csv", index=False)
    (dest / "primary_proteoforms.txt").write_text("\n".join(primary.Symbol) + "\n")
    metrics = pd.read_csv(out / "all_proteoform_cohort_metrics.csv.gz")
    records = []
    for code in COHORTS:
        vals = pd.read_parquet(out / f"checkpoints/{code}_patients.parquet")
        cuts = (
            pd.read_csv(out / f"checkpoints/{code}_cutoffs.csv")
            .set_index("patient_id")
            .loc[vals.columns]
        )
        for (f, p), selected in lists.items():
            keys = selected[KEY]
            mm = (
                metrics[metrics.cancer_code.eq(code) & metrics.transcriptome_percentile.eq(90)]
                .set_index(KEY)
                .loc[keys]
            )
            assert mm.complete_measurement.all(), code
            values = vals.loc[keys].to_numpy()
            assert np.isfinite(values).all()
            for evaluation in (70, 90):
                positive = (values > cuts[f"p{evaluation}_cutoff"].to_numpy()).any(axis=0)
                records.append(
                    {
                        "cancer_code": code,
                        "selection_prevalence": f,
                        "selection_percentile": p,
                        "n_proteins": len(keys),
                        "evaluation_percentile": evaluation,
                        "n_patients": len(positive),
                        "n_expressing_any": int(positive.sum()),
                        "fraction_expressing_any": float(positive.mean()),
                    }
                )
    coverage = pd.DataFrame(records)
    coverage.to_csv(dest / "same_threshold_panel_comparison.csv", index=False)
    selected_metrics = metrics[
        metrics.cancer_code.isin(COHORTS)
        & metrics.transcriptome_percentile.eq(90)
        & metrics[KEY].isin(primary[KEY])
    ]
    selected_metrics.to_csv(dest / "primary_common_cancer_heatmap.csv", index=False)
    cohort_sizes = coverage.drop_duplicates("cancer_code").set_index("cancer_code").n_patients
    rows90 = (
        coverage[coverage.evaluation_percentile.eq(90)]
        .pivot(index="cancer_code", columns="n_proteins", values="fraction_expressing_any")
        .loc[COHORTS]
    )
    wins = int(rows90[13].gt(rows90[38]).sum())
    assert wins == 10
    manifest = {
        "primary_selection": ">50% prevalence above each patient's p90 in at least one eligible cohort",
        "n_primary_proteins": 13,
        "n_benchmark_cohorts": len(COHORTS),
        "n_cohorts_13_beats_38_at_p90": wins,
        "benchmark_cohorts": COHORTS,
        "all_panel_members_measured_in_benchmark": True,
        "comparison_rule": "Patient-level OR at identical p70 or p90 evaluation thresholds; no pooled cohort or subtype maxima",
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_under("white")
    cmap.set_bad("#D4D4D4")

    def save(fig, name):
        fig.savefig(dest / f"{name}.png", dpi=165, bbox_inches="tight")
        fig.savefig(dest / f"{name}.svg", bbox_inches="tight")
        plt.close(fig)

    labels = [f"{name}\n{code} (n={cohort_sizes[code]})" for name, code in zip(NAMES, COHORTS)]
    data = selected_metrics.pivot(
        index=KEY, columns="cancer_code", values="fraction_expressing"
    ).reindex(index=primary[KEY], columns=COHORTS)
    matrix = np.vstack([rows90[13].to_numpy(), data.to_numpy()])
    fig, ax = plt.subplots(figsize=(16, 10), layout="constrained")
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=np.nextafter(0.0, 1.0), vmax=1)
    for i in range(len(matrix)):
        for j, value in enumerate(matrix[i]):
            label = "<1%" if 0 < value < 0.005 else f"{value:.0%}"
            ax.text(
                j,
                i,
                label,
                ha="center",
                va="center",
                fontsize=10,
                color="white" if value > 0.55 else "#203746",
            )
    ax.axhline(0.5, color="#C8772B", lw=3)
    ax.set(
        xticks=range(12),
        xticklabels=labels,
        yticks=range(14),
        yticklabels=[
            "ANY of the 13 proteins",
            *[f"{r.Symbol} ({r.protein_length_aa} aa)" for r in primary.itertuples()],
        ],
        title="Primary 13-protein panel: expression above each patient's p90",
    )
    ax.tick_params(axis="x", rotation=40)
    plt.setp(ax.get_xticklabels(), ha="right")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Fraction of patients expressing")
    fig.supxlabel(
        "Broad cohorts for common cancers, without selecting a favorable molecular subtype. All listed proteins have complete member-gene data here.\nANY counts each patient once if one or more listed proteins exceed p90. White = exactly zero; low nonzero frequencies remain colored.",
        fontsize=10,
    )
    save(fig, "primary_common_cancer_heatmap")

    fig, axes = plt.subplots(2, 1, figsize=(16, 9), layout="constrained")
    for ax, eval_p in zip(axes, [70, 90]):
        arr = []
        for f, p in SETTINGS:
            sub = coverage[
                coverage.selection_prevalence.eq(f)
                & coverage.selection_percentile.eq(p)
                & coverage.evaluation_percentile.eq(eval_p)
            ].set_index("cancer_code")
            arr.append(sub.loc[COHORTS].fraction_expressing_any.to_numpy())
        arr = np.array(arr)
        im = ax.imshow(arr, aspect="auto", cmap=cmap, vmin=np.nextafter(0.0, 1.0), vmax=1)
        for i in range(4):
            for j, value in enumerate(arr[i]):
                ax.text(
                    j,
                    i,
                    "<1%" if 0 < value < 0.005 else f"{value:.0%}",
                    ha="center",
                    va="center",
                    color="white" if value > 0.55 else "#203746",
                )
        ax.add_patch(Rectangle((-0.49, 1.51), 11.98, 0.98, fill=False, ec="#C8772B", lw=2.5))
        ax.set(
            xticks=range(12),
            xticklabels=NAMES,
            yticks=range(4),
            yticklabels=[f"{len(lists[f, p])} proteins: >{f}% / p{p}" for f, p in SETTINGS],
            title=f"Same evaluation threshold for every list: above patient p{eval_p}",
        )
        ax.tick_params(axis="x", rotation=25)
    fig.colorbar(
        im, ax=axes, fraction=0.018, pad=0.02, label="Patients expressing ANY listed protein"
    )
    fig.suptitle(
        "Why lead with 13 proteins? Compare list size at the same expression threshold", fontsize=17
    )
    fig.supxlabel(
        "Orange box = recommended primary panel. The 38-protein list offers wider moderate-expression coverage at p70;\nthe 13-protein list has higher p90 coverage in 10 of these 12 broad cohorts. The two lists are not nested. White = zero.",
        fontsize=10,
    )
    save(fig, "list_size_coverage_comparison")

    fig, ax = plt.subplots(figsize=(12, 7), layout="constrained")
    x = "normal_rna_max_member_somatic_ntpm"
    y = "p90_fraction_cancer_type_groups_gt10"
    ax.scatter(
        primary[x],
        primary[y] * 100,
        c=np.where(primary.n_member_genes.gt(1), "#C8772B", "#127C80"),
        s=65,
    )
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set(
        xlim=(-0.015, primary[x].max() * 1.6),
        ylim=(-3, 52),
        xlabel="Maximum member-gene normal somatic RNA (HPA nTPM)",
        ylabel="Cancer-type groups with >10% p90-positive patients (%)",
        title="Primary panel: cancer coverage and normal-tissue RNA evidence",
    )
    protein_gene_legend(ax)
    label_coverage_points(ax, primary, y, top_n=13, x_column=x)
    fig.supxlabel(
        "All 13 proteins are labeled. Coverage uses the fixed 69 cancer-type groups across the complete analysis.\nThe x-axis is the existing single-member HPA annotation; it excludes reproductive tissues and thymus. RNA abundance is not protein abundance.",
        fontsize=10,
    )
    save(fig, "primary_normal_tissue_scatter")
    print(f"Prepared primary panel; p90 comparison favors 13 over 38 in {wins}/12 cohorts")


def render(out):
    from PIL import Image
    from pypdf import PdfReader, PdfWriter
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A3, landscape, letter
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.pdfgen import canvas
    from reportlab.platypus import Paragraph, Table, TableStyle

    dest = out / "primary_panel"
    front = io.BytesIO()
    c = canvas.Canvas(front, pagesize=letter)
    teal = colors.HexColor("#127C80")
    ink = colors.HexColor("#203746")
    style = ParagraphStyle("body", fontName="Helvetica", fontSize=11, leading=15, textColor=ink)
    c.setFillColor(teal)
    c.rect(0, 750, 612, 42, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica", 10)
    c.drawString(36, 768, "CTA ANALYSIS | PRIMARY PANEL AND SUPPORTING EVIDENCE")
    c.setFillColor(ink)
    c.setFont("Helvetica-Bold", 23)
    c.drawString(36, 712, "Lead with the 13-protein panel")
    y = 687

    def para(text):
        nonlocal y
        p = Paragraph(text, style)
        _, h = p.wrap(540, 1000)
        p.drawOn(c, 36, y - h)
        y -= h + 14

    para(
        "<b>Recommendation:</b> present the existing <b>>50% prevalence / p90</b> shortlist as the primary panel. It balances a manageable number of protein identities with strong-expression coverage in several common cancers. Keep the 38-protein p70 list as the broader, moderate-expression alternative."
    )
    para(
        "<b>Selection:</b> each protein exceeds the patient's own transcriptome p90 in >50% of patients in at least one cohort with >=20 patient groups. Identical-protein gene loci are combined. Selection in one cancer does not imply high prevalence in every cancer."
    )
    para(
        "<b>Why this list:</b> using the same p90 evaluation threshold, 13 proteins cover more patients than the 38-protein list in <b>10 of 12 broad common-cancer cohorts</b>. At p70, the 38-protein panel often covers more patients. The two panels are not nested; this is an expression-threshold tradeoff, not simply a list-size effect."
    )
    coverage = pd.read_csv(dest / "same_threshold_panel_comparison.csv")
    vals = coverage[coverage.n_proteins.eq(13) & coverage.evaluation_percentile.eq(90)].set_index(
        "cancer_code"
    )
    rows = [["Common-cancer cohort", "Patients expressing >=1 panel protein at p90"]]
    for code, name in [
        ("LUAD", "Lung adenocarcinoma"),
        ("LUSC", "Lung squamous carcinoma"),
        ("HNSC", "Head/neck squamous carcinoma"),
        ("BRCA", "Breast invasive carcinoma"),
    ]:
        row = vals.loc[code]
        rows.append([f"{name} (n={int(row.n_patients)})", f"{row.fraction_expressing_any:.1%}"])
    table = Table(rows, colWidths=[285, 255])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), teal),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 9),
                ("TOPPADDING", (0, 0), (-1, -1), 9),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#EDF5F5"), colors.white]),
            ]
        )
    )
    _, h = table.wrap(540, 1000)
    table.drawOn(c, 36, y - h)
    y -= h + 16
    para(
        "<b>Coverage gaps:</b> p90 coverage is 0.6% in pancreatic adenocarcinoma, 3.4% in prostate, 6.2% in colon and 12.0% in rectal cancer. These are observed cohort fractions, not population-weighted or global estimates. No favorable subtype maxima are used in this opening comparison."
    )
    para(
        "<b>Read first:</b> page 2 shows all 13 proteins across common cancers; page 3 compares all four lists at matched thresholds; pages 4-6 show labeled normal-tissue evidence; page 7 gives the complete primary shortlist. The full analysis follows."
    )
    assert y > 36, y
    c.showPage()
    for filename in [
        "primary_common_cancer_heatmap",
        "list_size_coverage_comparison",
        "primary_normal_tissue_scatter",
    ]:
        path = dest / f"{filename}.png"
        with Image.open(path) as im:
            iw, ih = im.size
        size = landscape(A3)
        w, h = size
        c.setPageSize(size)
        scale = min((w - 40) / iw, (h - 54) / ih)
        c.drawImage(
            str(path),
            (w - iw * scale) / 2,
            (h - ih * scale) / 2 + 8,
            width=iw * scale,
            height=ih * scale,
        )
        c.showPage()
    c.save()
    front.seek(0)
    original = PdfReader(out / "cta-analysis-all-figures.pdf")
    writer = PdfWriter()
    titles = [
        "Recommendation: primary 13-protein panel",
        "Primary panel across common cancers",
        "List size versus matched-threshold coverage",
        "Labeled normal-tissue evidence",
    ]
    for page, title in zip(PdfReader(front).pages, titles):
        writer.add_page(page)
        writer.add_outline_item(title, len(writer.pages) - 1)
    index = pd.read_csv(out / "pdf_page_index.csv")
    first = [
        int(index.loc[index.figure_name.eq(name), "pdf_page"].iloc[0]) - 1
        for name in ["hpa_rna_tissue_overview", "hpa_ihc_tissue_overview"]
    ]
    first.append(
        int(index.loc[index.title.eq(">50% / p90: shortlist rows 1-13"), "pdf_page"].iloc[0]) - 1
    )
    for old in first:
        writer.add_page(original.pages[old])
        writer.add_outline_item(index.iloc[old].title, len(writer.pages) - 1)
    parent = writer.add_outline_item(
        "Supporting analysis: methods, alternatives and all figures", len(writer.pages)
    )
    page_map = []
    for old, page in enumerate(original.pages):
        if old in first:
            continue
        writer.add_page(page)
        writer.add_outline_item(index.iloc[old].title, len(writer.pages) - 1, parent=parent)
        page_map.append(
            {
                "original_page": old + 1,
                "new_page": len(writer.pages),
                "title": index.iloc[old].title,
            }
        )
    for number, page in enumerate(writer.pages, 1):
        w, h = float(page.mediabox.width), float(page.mediabox.height)
        overlay = io.BytesIO()
        cc = canvas.Canvas(overlay, pagesize=(w, h))
        cc.setFillColor(colors.white)
        cc.rect(w - 76, 0, 76, 32, fill=1, stroke=0)
        cc.setFillColor(ink)
        cc.setFont("Helvetica", 8)
        cc.drawRightString(w - 36, 18, str(number))
        cc.save()
        overlay.seek(0)
        page.merge_page(PdfReader(overlay).pages[0])
    writer.add_metadata(
        {"/Title": "CTA analysis: primary 13-protein panel and supporting evidence"}
    )
    output = out / "cta-analysis-primary-panel.pdf"
    with output.open("wb") as f:
        writer.write(f)
    pd.DataFrame(page_map).to_csv(dest / "supporting_page_map.csv", index=False)
    (dest / "README.md").write_text(
        "# Primary CTA panel\n\nRecommended existing list: 13 proteins, >50% prevalence above each patient's p90 in at least one eligible cohort.\n\n- primary_proteoforms.csv and .txt: chosen identities, lengths and existing evaluation fields.\n- same_threshold_panel_comparison.csv: all four lists evaluated at both p70 and p90 in the same 12 broad cohorts.\n- primary_common_cancer_heatmap.csv: all individual-protein prevalence values.\n- supporting_page_map.csv: original-to-new supporting page mapping.\n\nThe 13-protein and 38-protein lists are not nested. Coverage is patient-level OR within cohorts, with no pooled patients or favorable subtype selection. All member-gene data are available in these benchmark cohorts.\n"
    )
    print(f"Created {output}: {len(writer.pages)} pages")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["prepare", "render"])
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    (prepare if args.phase == "prepare" else render)(args.out)
