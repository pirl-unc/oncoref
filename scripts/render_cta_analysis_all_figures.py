#!/usr/bin/env python3
"""Assemble the complete CTA lists and figure atlas using reportlab."""

from __future__ import annotations

import argparse
import json
import shutil
import zipfile
from html import escape
from pathlib import Path

import pandas as pd
from cta_report_protein_lengths import LENGTH_METHOD, annotate_report
from PIL import Image
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape, letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--out", type=Path, default=root / "outputs/cta_threshold_report_20260917")
    parser.add_argument("--base", type=Path, default=root / "outputs/cta_threshold_report_20260916")
    args = parser.parse_args()
    out, base = args.out.resolve(), args.base.resolve()
    proteins = pd.read_csv(out / "cta_protein_lengths.csv")
    annotate_report(out, proteins)
    summary = pd.read_csv(out / "selection_summary.csv")
    sizes = pd.read_csv(out / "cohort_sizes.csv")
    availability = pd.read_csv(out / "cohort_size_availability.csv")
    strict = pd.read_csv(out / "strictest_genes_ranked.csv")
    strict_pairs = pd.read_csv(out / "strictest_passing_cancer_types.csv")
    focused_strict_pairs = pd.read_csv(
        out / "selections/large_linear/prevalence_gt70_transcriptome_p90/passing_gene_cohorts.csv"
    )
    plots = pd.read_csv(out / "plot_index.csv")
    groups = pd.read_csv(out / "cohort_overlap_groups.csv")
    manifest = json.loads((out / "run_manifest.json").read_text())
    min_size, large_size = manifest["min_cohort_size"], manifest["large_min_cohort_size"]
    n_small = int(sizes.passes_default_size.sum())
    n_large = int(sizes.passes_large_linear_filter.sum())
    n_group = groups.overlap_group.nunique()
    base_manifest = json.loads((base / "run_manifest.json").read_text())
    validation = json.loads((out / "validation.json").read_text())
    pdf_path = out / "cta-analysis-all-figures.pdf"
    pdf = canvas.Canvas(str(pdf_path), pagesize=letter)
    pdf.setTitle("CTA analysis: complete gene lists and all figures")
    pdf.setAuthor("oncoref analysis")
    teal, ink, pale = (
        colors.HexColor("#127C80"),
        colors.HexColor("#203746"),
        colors.HexColor("#EDF5F5"),
    )
    body = ParagraphStyle("body", fontName="Helvetica", fontSize=10, leading=14, textColor=ink)
    small = ParagraphStyle("small", parent=body, fontSize=8, leading=10)
    method_style = ParagraphStyle("methods", parent=body, fontSize=9, leading=12)
    cell_style = ParagraphStyle("cell", parent=small, fontSize=7.5, leading=9)
    white_cell = ParagraphStyle("white", parent=cell_style, textColor=colors.white)
    page = 0
    page_index = []

    def para(text, y, width=540, style=body, x=36):
        obj = Paragraph(text, style)
        _, h = obj.wrap(width, 2000)
        obj.drawOn(pdf, x, y - h)
        return y - h - 10

    def header(title, subtitle=None):
        pdf.setPageSize(letter)
        pdf.setFillColor(teal)
        pdf.rect(0, 758, 612, 34, fill=1, stroke=0)
        pdf.setFillColor(colors.white)
        pdf.setFont("Helvetica-Bold", 11)
        pdf.drawString(36, 770, "CTA ANALYSIS | SEPTEMBER 17, 2026")
        pdf.setFillColor(ink)
        pdf.setFont("Helvetica-Bold", 20)
        pdf.drawString(36, 724, title)
        return para(subtitle, 699) if subtitle else 696

    def finish(title, category, name=None, size=letter):
        nonlocal page
        page += 1
        key = f"page_{page}"
        pdf.bookmarkPage(key)
        pdf.addOutlineEntry(title, key, level=0)
        w, _ = size
        pdf.setFillColor(colors.HexColor("#667A86"))
        pdf.setFont("Helvetica", 7)
        pdf.drawString(
            36, 18, "oncoref | Within-patient transcriptome percentiles | Patient groups after QC"
        )
        pdf.drawRightString(w - 36, 18, str(page))
        page_index.append(
            {"pdf_page": page, "category": category, "title": title, "figure_name": name or ""}
        )
        pdf.showPage()

    def table_at(data, widths, y, x=36, font_size=8):
        rows = []
        for i, row in enumerate(data):
            style = white_cell if i == 0 else cell_style
            rows.append([Paragraph(escape(str(v)), style) for v in row])
        t = Table(rows, colWidths=widths, repeatRows=1)
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), teal),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [pale, colors.white]),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                    ("LEFTPADDING", (0, 0), (-1, -1), 5),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ]
            )
        )
        _, h = t.wrap(sum(widths), 2000)
        if y - h < 38:
            raise ValueError(f"Table does not fit on page {page + 1}: bottom={y - h}")
        t.drawOn(pdf, x, y - h)
        return y - h - 14

    y = header(
        "CTA analysis - all figures",
        "Complete gene lists, cohort-size support, overlap-adjusted breadth and patient-level expression.",
    )
    y = para(
        "<b>Focused settings:</b> prevalence strictly above 50% or 70% of a cohort, with expression strictly above the patient's transcriptome p70 or p90. The original six-gene result at >75% / p90 remains a separate case study.",
        y,
    )
    y = para(
        f"<b>Three cohort policies:</b> all passing cohorts; a configurable default minimum of {min_size} patient groups; and the requested larger-cohort analysis using at least {large_size} patient groups and linear-TPM-comparable sources only.",
        y,
    )
    rows = [
        [
            "Prevalence",
            "Transcriptome",
            "All cohorts",
            f"n >= {min_size}",
            f"n >= {large_size}, linear TPM",
        ]
    ]
    for _, s in summary.iterrows():
        rows.append(
            [
                f">{int(s.prevalence_gt_pct)}%",
                f"p{int(s.transcriptome_percentile)}",
                int(s.n_genes_all),
                int(s.n_genes_min_size),
                int(s.n_genes_large_linear),
            ]
        )
    y = table_at(rows, [95, 100, 100, 100, 145], y)
    y = para(
        f"<b>Cohort sizes:</b> {len(sizes)} cohort views span {int(sizes.n_patients.min())}-{int(sizes.n_patients.max()):,} patient groups; median {sizes.n_patients.median():g}. The default n >= {min_size} retains {n_small} cohorts; n >= {large_size} plus linear TPM retains {n_large}. Known patient overlap links the 142 views into {n_group} support groups.",
        y,
    )
    y = para(
        "<b>What largest passing cohort means:</b> the maximum cohort size among cohorts that actually satisfy the gene's expression and prevalence filters. A large cohort in which the gene fails contributes nothing to this measure.",
        y,
    )
    y = para(
        "<b>Six-gene case study:</b> PRAME's largest passing cohort is melanoma (422/466 expressing); DPPA5's is testicular germ-cell tumors (112/148). CTAG1A and CTAG1B pass in a 28-patient myxoid liposarcoma proxy-scale cohort. EGFL6 and POTEF have maximum passing sizes of 1 and 9, respectively.",
        y,
    )
    y = para(
        f"<b>Navigation:</b> the PDF bookmarks contain every list and figure. There are {len(plots)} figures, including all-patient distributions for the six genes in all 142 cohorts. The same plots are supplied as PNG and SVG; CSV and GMT exports contain the full numerical results and gene sets.",
        y,
    )
    finish("Overview and focused gene-set counts", "overview")

    methods = [
        (
            "Input and expression",
            f"The analytical snapshot contains 293 canonical default CTA genes and 142 registered cancer cohorts, using oncoref {base_manifest['oncoref_version']}. Expression is gene-level biological clean TPM, with technical/censored genes excluded from the transcriptome denominator. All finite biological rows, including zeros, define each patient's linear-interpolated quantile. This is a within-patient expression percentile, not a percentile across patients.",
        ),
        (
            "Patient grouping and QC",
            "Known repeat specimens were averaged in linear clean-expression space before quantile calculation. Existing QC failures and two NUTM samples containing negative TPM-labeled values were excluded. The remaining valid NUTM sample has n=1. In 27 cohorts, some source samples lack verified donor mappings and stand in for patients. Thus n means patient groups, not a universally verified number of independent people.",
        ),
        (
            "Positive patients, means and boundaries",
            "Expression must be strictly greater than the patient's quantile; ties do not count. Prevalence must be strictly greater than 50% or 70%. Minimum cohort size is inclusive: n=10 passes a minimum of 10. Positive-subset expression means are arithmetic means in linear units. No positives gives a missing mean, not zero. Missing gene measurements do not qualify.",
        ),
        (
            "Largest passing cohort",
            "For each selected gene and threshold combination, take the maximum n only among passing cohorts. All tied largest cohorts are recorded in CSV. The size-only default and the n>=20 linear-TPM filter are separate policies: proxy measurements can pass the size-only filter while failing the linear-TPM filter.",
        ),
        (
            "Overlap-adjusted breadth",
            "Cohort views sharing a known patient identity are connected. Each connected component contributes at most one vote to a gene's breadth score. Components are fixed using all 142 cohort views before filtering, so molecular subtypes connected through a parent remain one group even if the parent itself does not pass. Separate disjoint histologies remain separate unless a shared identity connects them.",
        ),
        (
            "Identity and ranking limitations",
            "TCGA and TARGET case IDs have global namespaces; recognized Treehouse donor IDs are matched across Treehouse views; otherwise IDs are scoped to physical source cohorts. This removes observed duplicate support but cannot discover unrecorded overlap across sources. Counts describe conservative support groups, not independent studies or causal evidence. Equal breadth scores share a rank; largest passing n orders ties for display.",
        ),
        (
            "Proxy scale and varying gene universes",
            "Ten cohorts use microarray, CPM or single-cell proxy scales. Their positive fractions remain available, but their absolute expression means are labeled source proxy and are not interchangeable with TPM. The linear-only policy excludes them. Transcriptome percentiles depend on each assay's measured biological gene universe.",
        ),
        (
            "Patient-level plots",
            "Every displayed point is one retained patient group. Teal indicates expression above that group's p90, gray indicates below or tied; black diamonds mark the mean among positive groups. Expression/p90 ratio >1 is positive. Absolute-expression panels retain labeled source units. Symmetric-log axes include zero while showing a wide dynamic range.",
        ),
        ("Protein length (aa)", LENGTH_METHOD),
        (
            "Validation",
            f"Seven focused numerical tests cover quantile ties, exact prevalence and minimum-size boundaries, missingness, repeat-donor averaging, largest-support selection and transitive overlap. Full-data checks validate all 12 focused policy-specific gene sets and nested thresholds. Reconstructed patient-level data reproduce positive counts and means in {validation['patient_expression_gene_cohort_checks']} fully measured gene/cohort pairs.",
        ),
    ]
    y = header("Methods and definitions")
    for title, text in methods:
        block = f"<b>{escape(title)}.</b> {escape(text)}"
        probe = Paragraph(block, method_style)
        _, h = probe.wrap(540, 2000)
        if y - h < 55:
            finish("Methods and definitions", "methods")
            y = header("Methods and definitions (continued)")
        y = para(block, y, style=method_style)
    finish("Methods and definitions", "methods")

    y = header("Available cohort sizes")
    y = para(
        "Counts below use QC-eligible patient groups after known-repeat merging. These are cohort views, so totals include parent/subtype overlap.",
        y,
    )
    y = table_at(
        [
            ["Minimum size", "All scales: cohorts", "Linear TPM: cohorts"],
            *[
                [int(r.min_cohort_size), int(r.n_cohorts), int(r.n_linear_cohorts)]
                for _, r in availability.iterrows()
            ],
        ],
        [160, 190, 190],
        y,
    )
    y = para("Largest available cohorts", y)
    y = table_at(
        [
            ["Cohort", "Cancer type", "Patient groups"],
            *[
                [r.cancer_code, r.cancer_name, int(r.n_patients)]
                for _, r in sizes.nlargest(8, "n_patients").iterrows()
            ],
        ],
        [110, 330, 100],
        y,
    )
    finish("Available cohort sizes and largest cohorts", "sizes")

    # Every gene in every focused list is represented, with explicit policy membership.
    for _, s in summary.iterrows():
        all_genes = pd.read_csv(out / "selections/all" / s.scenario / "genes_ranked.csv")
        large = pd.read_csv(out / "selections/large_linear" / s.scenario / "genes_ranked.csv")
        rank_map = large.set_index("Symbol").breadth_rank.to_dict()
        group_map = large.set_index("Symbol").n_overlap_adjusted_groups.to_dict()
        all_genes["large_rank"] = all_genes.Symbol.map(rank_map)
        all_genes = all_genes.sort_values(
            ["large_rank", "n_overlap_adjusted_groups", "largest_passing_cohort_size", "Symbol"],
            ascending=[True, False, False, True],
            na_position="last",
        )
        for start in range(0, len(all_genes), 27):
            chunk = all_genes.iloc[start : start + 27]
            title = (
                f"Gene lists: >{int(s.prevalence_gt_pct)}% above p{int(s.transcriptome_percentile)}"
            )
            y = header(
                title,
                f"Complete list, rows {start + 1}-{start + len(chunk)} of {len(all_genes)}. All: {int(s.n_genes_all)} genes; n >= {min_size}: {int(s.n_genes_min_size)}; n >= {large_size}, linear TPM: {int(s.n_genes_large_linear)}.",
            )
            y = para(
                f"Every row belongs to the all-cohort set. Protein aa is the longest Ensembl 112 protein. Yes defines each filtered list. Rank and groups use only n >= {large_size}, linear-TPM passing cohorts. Largest n uses all passing cohorts.",
                y,
                style=small,
            )
            rows = [
                [
                    "Gene",
                    "Protein aa",
                    "Largest n",
                    "Largest passing cohort(s)",
                    "All overlap groups",
                    f"n >= {min_size}",
                    f"n >= {large_size} linear",
                    "Linear rank / groups",
                ]
            ]
            for _, g in chunk.iterrows():
                rows.append(
                    [
                        g.Symbol,
                        int(g.protein_length_aa),
                        int(g.largest_passing_cohort_size),
                        g.largest_passing_cohorts,
                        int(g.n_overlap_adjusted_groups),
                        "Yes" if g.passes_min_cohort_size else "No",
                        "Yes" if g.passes_large_linear_filter else "No",
                        f"{rank_map[g.Symbol]} / {group_map[g.Symbol]}"
                        if g.Symbol in rank_map
                        else "--",
                    ]
                )
            table_at(rows, [75, 42, 44, 115, 64, 48, 65, 87], y)
            finish(f"{title} - rows {start + 1}-{start + len(chunk)}", "gene_lists")

    y = header("Six-gene case study: >75% / p90")
    y = para(
        "This retains the original strictest six-gene result as a reference. It is separate from the new >70% focused setting. Size alone and linear-TPM scale are distinct requirements.",
        y,
    )
    rows = [
        [
            "Gene",
            "Protein aa",
            "Largest passing n",
            "Largest passing cohort",
            f"n >= {min_size}",
            f"n >= {large_size}, linear TPM",
        ]
    ]
    for _, g in strict.iterrows():
        rows.append(
            [
                g.Symbol,
                int(g.protein_length_aa),
                int(g.largest_passing_cohort_size),
                g.largest_passing_cohorts,
                "Yes" if g.passes_min_cohort_size else "No",
                "Yes" if g.passes_large_linear_filter else "No",
            ]
        )
    y = table_at(rows, [62, 45, 75, 148, 85, 125], y)
    y = para(
        "The two genes retained by the larger linear-TPM policy are <b>DPPA5 and PRAME</b>. EGFL6's only passing cohort has one patient group; POTEF's largest has nine. CTAG1A/B pass in myxoid liposarcoma, but that cohort uses a microarray-derived proxy scale.",
        y,
    )
    y = para(
        "The following pages list every passing cancer type, positive count, prevalence and positive-subset mean for all six. Proxy means are marked explicitly. Parent/subtype overlap is retained in the detailed rows and removed only when computing breadth scores.",
        y,
    )
    finish("Six-gene case study and largest support", "strictest")
    for start in range(0, len(strict_pairs), 17):
        chunk = strict_pairs.iloc[start : start + 17]
        y = header(
            "Passing cancer types for the six genes",
            f">75% prevalence above p90 | rows {start + 1}-{start + len(chunk)} of {len(strict_pairs)}",
        )
        rows = [
            [
                "Gene",
                "Protein aa",
                "Cancer type (cohort)",
                "Positive / n",
                "Positive %",
                "Mean positive",
                "Scale",
                f"n >= {large_size} linear",
            ]
        ]
        for _, r in chunk.iterrows():
            rows.append(
                [
                    r.Symbol,
                    int(r.protein_length_aa),
                    f"{r.cancer_name} ({r.cancer_code})",
                    f"{int(r.n_expressing)}/{int(r.n_patients)}",
                    f"{100 * r.fraction_expressing:.1f}%",
                    f"{r.mean_expression_expressing:.2f}",
                    "TPM" if r.linear_tpm_comparable else "Proxy",
                    "Yes" if r.passes_large_linear_filter else "No",
                ]
            )
        table_at(rows, [49, 39, 188, 62, 48, 60, 34, 60], y)
        finish(f"Six genes: passing cancer types {start + 1}-{start + len(chunk)}", "strictest")

    for start in range(0, len(focused_strict_pairs), 17):
        chunk = focused_strict_pairs.iloc[start : start + 17]
        y = header(
            "Stringent focus: passing cancer types",
            f">70% prevalence above p90 | n >= {large_size}, linear TPM only | "
            f"rows {start + 1}-{start + len(chunk)} of {len(focused_strict_pairs)}",
        )
        y = para(
            "<b>DPPA5 and PRAME</b> pass this policy. Compared with >75%, PRAME adds "
            "angiosarcoma (15/20) and MSI endometrial cancer (30/41). All pairs on the "
            "preceding >75% pages also pass >70% before applying the cohort-size and scale requirements.",
            y,
        )
        rows = [
            [
                "Gene",
                "Protein aa",
                "Cancer type (cohort)",
                "Positive / n",
                "Positive %",
                "Mean positive TPM",
            ]
        ]
        for _, r in chunk.iterrows():
            rows.append(
                [
                    r.Symbol,
                    int(r.protein_length_aa),
                    f"{r.cancer_name} ({r.cancer_code})",
                    f"{int(r.n_expressing)}/{int(r.n_patients)}",
                    f"{100 * r.fraction_expressing:.1f}%",
                    f"{r.mean_expression_expressing:.2f}",
                ]
            )
        table_at(rows, [52, 42, 226, 70, 60, 90], y)
        finish(
            f">70% / p90: larger linear-TPM passing cancer types {start + 1}-{start + len(chunk)}",
            "focused_strictest",
        )

    # Natural figure dimensions are preserved on A3 pages for legible dense axes.
    for category in ["main", "appendix"]:
        subset = plots[plots.category.eq(category)]
        for _, r in subset.iterrows():
            path = out / "plots" / f"{r['name']}.png"
            with Image.open(path) as im:
                iw, ih = im.size
            page_size = landscape(A3) if iw / ih > 1.15 else A3
            w, h = page_size
            pdf.setPageSize(page_size)
            scale = min((w - 42) / iw, (h - 62) / ih)
            rw, rh = iw * scale, ih * scale
            pdf.drawImage(str(path), (w - rw) / 2, (h - rh) / 2 + 8, width=rw, height=rh)
            finish(str(r.title), category, str(r["name"]), page_size)
    pdf.save()
    pd.DataFrame(page_index).to_csv(out / "pdf_page_index.csv", index=False)

    # Carry the numerical source snapshot and code with the final atlas.
    for name in [
        "cohort_audit.csv",
        "cta_input_universe.csv",
        "sample_patient_audit.csv.gz",
        "patient_transcriptome_cutoffs.csv.gz",
    ]:
        shutil.copy2(base / name, out / name)
    master = pd.read_csv(base / "all_gene_cohort_metrics.csv.gz")
    master[master.transcriptome_percentile.isin([70, 90])].to_csv(
        out / "all_gene_cohort_metrics_p70_p90.csv.gz", index=False
    )
    shutil.copy2(base / "run_manifest.json", out / "base_run_manifest.json")
    repro = out / "reproducibility"
    repro.mkdir(exist_ok=True)
    for name in [
        "cta_threshold_report.py",
        "cta_threshold_report_details.py",
        "render_cta_analysis_all_figures.py",
        "cta_report_protein_lengths.py",
    ]:
        shutil.copy2(root / "scripts" / name, repro / name)
    for name in ["test_cta_threshold_report.py", "test_cta_threshold_report_details.py"]:
        shutil.copy2(root / "tests" / name, repro / name)
    methods_md = "\n\n".join(f"**{title}.** {text}" for title, text in methods)
    (out / "methodology.md").write_text("# CTA analysis methodology\n\n" + methods_md + "\n")
    md = [
        "# CTA analysis - complete lists and all figures",
        "",
        "Updated September 17, 2026.",
        "",
        f"[Complete PDF atlas: {page} pages](cta-analysis-all-figures.pdf) | [Download bundle](cta-analysis-bundle.zip)",
        "",
        f"Focus: **>50% and >70% cohort prevalence x transcriptome p70 and p90**. Default minimum cohort size: **{min_size}**. Requested larger-cohort policy: **n >= {large_size}, linear TPM only**.",
        "",
        "| Prevalence | Transcriptome | All cohorts | Minimum n=10 | n>=20, linear TPM |",
        "|---|---|---:|---:|---:|",
    ]
    for _, s in summary.iterrows():
        md.append(
            f"| >{int(s.prevalence_gt_pct)}% | p{int(s.transcriptome_percentile)} | {int(s.n_genes_all)} | {int(s.n_genes_min_size)} | {int(s.n_genes_large_linear)} |"
        )
    md += [
        "",
        f"Cohort sizes range from {int(sizes.n_patients.min())} to {int(sizes.n_patients.max()):,}, median {sizes.n_patients.median():g}. {n_small} meet the default size floor; {n_large} meet the larger linear-TPM policy. Known-overlap adjustment yields {n_group} support groups from 142 cohort views.",
        "",
        "## Largest passing cohorts: original six-gene list",
        "",
        "| Gene | Protein length (aa) | Largest passing n | Cancer type(s) | Pass n>=10 | Pass n>=20, linear |",
        "|---|---:|---:|---|---|---|",
    ]
    for _, g in strict.iterrows():
        md.append(
            f"| {g.Symbol} | {int(g.protein_length_aa)} | {int(g.largest_passing_cohort_size)} | {g.largest_passing_cancer_types} | {'Yes' if g.passes_min_cohort_size else 'No'} | {'Yes' if g.passes_large_linear_filter else 'No'} |"
        )
    md += [
        "",
        "[All passing cancer types, counts and positive means](strictest_passing_cancer_types.csv)",
        "[Focused >70% / p90 passing cancer types, n>=20 linear TPM](selections/large_linear/prevalence_gt70_transcriptome_p90/passing_gene_cohorts.csv)",
        "",
        "## Gene sets and rankings",
        "",
        "Protein lengths use the longest annotated Ensembl 112 protein per gene. [Lengths and chosen protein IDs](cta_protein_lengths.csv) | [Every annotated isoform](cta_protein_isoforms.csv) | [Protein source and method](protein_length_manifest.json)",
        "",
        "[All-cohort GMT](gene_sets_all.gmt) | [Minimum-size GMT](gene_sets_min_size.gmt) | [Larger linear-TPM GMT](gene_sets_large_linear.gmt)",
        "",
    ]
    for _, s in summary.iterrows():
        md.append(
            f"- >{int(s.prevalence_gt_pct)}% / p{int(s.transcriptome_percentile)}: "
            + " | ".join(
                f"[{label}](selections/{policy}/{s.scenario}/genes_ranked.csv)"
                for policy, label in [
                    ("all", "All passing cohorts"),
                    ("min_size", f"Minimum n={min_size}"),
                    ("large_linear", f"n>={large_size}, linear TPM"),
                ]
            )
        )
    md += [
        "",
        "## Figures",
        "",
        f"All {len(plots)} figures appear in the PDF, with bookmarks and a [page index](pdf_page_index.csv). PNG and SVG versions are in `plots/`.",
        "",
    ]
    for _, r in plots.iterrows():
        md.append(f"- [{r.title}](plots/{r['name']}.png)")
    md += [
        "",
        "## Reproduce",
        "",
        "From this oncoref checkout:",
        "",
        "```sh",
        "MPLCONFIGDIR=$PWD/tmp/cta_threshold_report/matplotlib CANCERDATA_PER_SAMPLE_CACHE=1 .venv/bin/python scripts/cta_threshold_report_details.py --min-cohort-size 10 --large-min-cohort-size 20 --reuse-patient-values",
        ".venv/bin/python scripts/cta_report_protein_lengths.py",
        "```",
        "",
        "Use `--reuse-patient-values` only when the base snapshot and source matrices are unchanged. The PDF renderer needs reportlab, pandas and Pillow. The original numerical inputs and reconstruction audits remain in the September 16 report directory.",
        "",
        "[Methods and limitations](methodology.md) | [Validation](validation.json) | [Cohort-size availability](cohort_size_availability.csv) | [Overlap group map](cohort_overlap_groups.csv)",
        "",
        "Patient counts use known donor grouping; 27 cohorts have incomplete donor identities. Proxy units and unknown cross-study overlap remain explicit limitations.",
    ]
    (out / "report.md").write_text("\n".join(md) + "\n")
    zip_path = out / "cta-analysis-bundle.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.rglob("*")):
            if (
                p.is_file()
                and p != zip_path
                and p.suffix not in [".log", ".parquet"]
                and "qa" not in p.parts
            ):
                z.write(p, p.relative_to(out))
    print(
        f"Created {page}-page PDF with {len(plots)} figures and complete lists; bundle {zip_path.stat().st_size / 1e6:.1f} MB"
    )


if __name__ == "__main__":
    main()
