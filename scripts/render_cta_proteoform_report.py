#!/usr/bin/env python3
"""Assemble the focused proteoform shortlists and every figure into one PDF."""

from __future__ import annotations

import argparse
import json
import shutil
import zipfile
from html import escape
from pathlib import Path

import pandas as pd
from PIL import Image
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape, letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
KEY = "proteoform_key"
FOCUSED = [(50, 70), (70, 70), (50, 90), (70, 90)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/cta_proteoform_report_20260917")
    out = parser.parse_args().out.resolve()
    universe = pd.read_csv(out / "proteoform_universe.csv")
    cohorts = pd.read_csv(out / "cohort_audit.csv")
    eligible = cohorts[cohorts.n_patients.ge(20)]
    summary = pd.read_csv(out / "selection_summary.csv")
    plots = pd.read_csv(out / "plot_index.csv")
    n_types = int(summary.n_eligible_cancer_type_groups.iloc[0])
    n_views = len(eligible)
    lists = {
        (f, p): pd.read_csv(
            out / f"selections/min20/prevalence_gt{f}_transcriptome_p{p}/proteoforms_ranked.csv"
        )
        for f, p in FOCUSED
    }
    text_reports = []
    for (f, p), selected in lists.items():
        folder = out / f"selections/min20/prevalence_gt{f}_transcriptome_p{p}"
        (folder / "proteoform_symbols.txt").write_text("\n".join(selected.Symbol) + "\n")
        columns = [
            "Symbol",
            "protein_length_aa",
            "member_symbols",
            "p70_fraction_cancer_type_groups_gt10",
            "p90_fraction_cancer_type_groups_gt10",
            "largest_passing_cohorts",
            "largest_passing_cohort_size",
            "largest_cohort_prevalence",
            "largest_cohort_positive_mean_tpm",
        ]
        text = (
            f"CTA PROTEOFORM SHORTLIST: >{f}% prevalence above patient p{p}\n"
            f"{len(selected)} protein identities. At least 20 patient groups per cohort.\n"
            "Each patient's own transcriptome percentile determines positivity.\n"
            "Coverage columns use >10% positive patients across 69 cancer-type groups.\n"
            "Fractions are on a 0-1 scale. Mean TPM uses positive patients in the largest passing cohort; NA = unavailable/non-TPM.\n\n"
            + selected[columns].to_string(
                index=False, na_rep="NA", float_format=lambda value: f"{value:.4g}"
            )
            + "\n"
        )
        (folder / "shortlist_report.txt").write_text(text)
        text_reports.append(text)
    (out / "all_cutoff_reports.txt").write_text("\n\n".join(text_reports))
    stringent = lists[70, 90]
    # Compact cross-list evaluation export with explicit fixed denominators.
    base = lists[50, 70]
    chosen_columns = [
        KEY,
        "Symbol",
        "member_symbols",
        "member_gene_ids",
        "n_member_genes",
        "additional_member_symbols",
        "protein_length_aa",
        "protein_ids",
        "normal_rna_max_member_somatic_ntpm",
        "normal_rna_max_member_and_tissue",
        "normal_rna_members_with_evidence",
        "member_restriction_confidence",
        "member_annotation_flags",
        "max_member_isoform_count",
    ]
    chosen_columns += [
        f"p{p}_{name}"
        for p in (70, 90)
        for name in [
            "n_cancer_type_groups_gt10",
            "n_eligible_cancer_type_groups",
            "fraction_cancer_type_groups_gt10",
            "n_cohort_views_gt10",
            "n_eligible_cohort_views",
            "fraction_cohort_views_gt10",
            "n_fully_measured_cohort_views",
            "n_cohort_views_gt10_allow_partial",
        ]
    ]
    evaluation = base[chosen_columns].copy()
    for (f, p), data in lists.items():
        key = f"gt{f}_p{p}"
        evaluation[f"selected_{key}"] = evaluation[KEY].isin(data[KEY])
        for column in [
            "largest_passing_cohort_size",
            "largest_passing_cohorts",
            "largest_cohort_prevalence",
            "largest_cohort_positive_mean_tpm",
            "n_passing_cancer_type_groups",
        ]:
            evaluation[f"{key}_{column}"] = evaluation[KEY].map(data.set_index(KEY)[column])
    evaluation = evaluation.sort_values(
        ["p90_fraction_cancer_type_groups_gt10", "p70_fraction_cancer_type_groups_gt10", "Symbol"],
        ascending=[False, False, True],
    )
    evaluation.to_csv(out / "shortlist_evaluation.csv", index=False)

    pdf = canvas.Canvas(str(out / "cta-analysis-all-figures.pdf"), pagesize=letter)
    pdf.setTitle("CTA proteoform analysis: focused shortlists and all figures")
    pdf.setAuthor("oncoref analysis")
    teal, ink, pale = (
        colors.HexColor("#127C80"),
        colors.HexColor("#203746"),
        colors.HexColor("#EDF5F5"),
    )
    body = ParagraphStyle("body", fontName="Helvetica", fontSize=10, leading=14, textColor=ink)
    small = ParagraphStyle("small", parent=body, fontSize=9, leading=12)
    cell = ParagraphStyle("cell", parent=body, fontSize=8, leading=10)
    white_cell = ParagraphStyle("white", parent=cell, textColor=colors.white)
    page, page_index = 0, []

    def paragraph(text, y, width=540, style=body, x=36):
        item = Paragraph(text, style)
        _, h = item.wrap(width, 2000)
        item.drawOn(pdf, x, y - h)
        return y - h - 10

    def header(title, subtitle=None, size=letter):
        pdf.setPageSize(size)
        w, h = size
        pdf.setFillColor(teal)
        pdf.rect(0, h - 34, w, 34, fill=1, stroke=0)
        pdf.setFillColor(colors.white)
        pdf.setFont("Helvetica-Bold", 11)
        pdf.drawString(36, h - 22, "CTA PROTEOFORM ANALYSIS | SEPTEMBER 17, 2026")
        pdf.setFillColor(ink)
        pdf.setFont("Helvetica-Bold", 20)
        pdf.drawString(36, h - 68, title)
        return paragraph(subtitle, h - 94, w - 72) if subtitle else h - 94

    def finish(title, category, size=letter, figure=""):
        nonlocal page
        page += 1
        key = f"page_{page}"
        pdf.bookmarkPage(key)
        pdf.addOutlineEntry(title, key, level=0)
        pdf.setFillColor(colors.HexColor("#657985"))
        pdf.setFont("Helvetica", 7)
        pdf.drawString(
            36, 18, "oncoref | Identical-protein expression sums | Cohort floor: 20 patient groups"
        )
        pdf.drawRightString(size[0] - 36, 18, str(page))
        pdf.showPage()
        page_index.append(
            {"pdf_page": page, "title": title, "category": category, "figure_name": figure}
        )

    def table(rows, widths, y):
        cooked = [
            [Paragraph(escape(str(v)), white_cell if i == 0 else cell) for v in row]
            for i, row in enumerate(rows)
        ]
        item = Table(cooked, colWidths=widths, repeatRows=1)
        item.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), teal),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [pale, colors.white]),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                ]
            )
        )
        _, h = item.wrap(sum(widths), 2000)
        if y - h < 40:
            raise ValueError(f"Table overflow, page {page + 1}, bottom {y - h}")
        item.drawOn(pdf, 36, y - h)
        return y - h - 14

    def frac(row, p, types=True):
        noun = "cancer_type_groups" if types else "cohort_views"
        n, d = int(row[f"p{p}_n_{noun}_gt10"]), int(row[f"p{p}_n_eligible_{noun}"])
        return f"{n}/{d} ({n / d:.1%})"

    y = header("CTA analysis at the proteoform level")
    y = paragraph(
        "<b>293 CTA genes map to 268 protein identities.</b> Eighteen identities combine multiple identical-protein loci, including CTAG1A/CTAG1B as <b>NY-ESO-1</b>. Seven additional loci encoding those same proteins are included, for 300 member genes.",
        y,
    )
    y = paragraph(
        f"<b>Focused settings:</b> >50% or >70% positive patients, at transcriptome p70 or p90. All {n_views} cohorts with at least 20 patient groups are eligible, including non-TPM sources. These represent {n_types} cancer-type groups after subtype and overlap adjustment.",
        y,
    )
    rows = [["Required cohort prevalence", "Above p70", "Above p90"]]
    for f in (50, 70):
        rows.append([f">{f}%", len(lists[f, 70]), len(lists[f, 90])])
    y = table(rows, [270, 135, 135], y)
    y = paragraph(
        "<b>Stringent shortlist: >70% / p90.</b> Broader coverage below is measured separately at >10% patient prevalence. It includes cancer types with lower frequencies than the selection threshold.",
        y,
    )
    rows = [["Proteoform", "Length aa", ">10% at p70: cancer types", ">10% at p90: cancer types"]]
    for _, r in stringent.iterrows():
        rows.append([r.Symbol, int(r.protein_length_aa), frac(r, 70), frac(r, 90)])
    y = table(rows, [120, 70, 175, 175], y)
    y = paragraph(
        "<b>Heatmaps preserve lower frequencies.</b> Columns enter when any shortlisted protein is positive in >10% of that cohort. Every cell in an included column is displayed, including low nonzero frequencies and measured zeros. Orange outlines mark the cells that pass the stringent selection rule.",
        y,
    )
    y = paragraph(
        f"<b>Evaluation tables:</b> complete shortlists with length, both coverage denominators, largest passing cohort, positive-patient mean TPM and available normal-tissue RNA annotations. All {len(plots)} figures are included with bookmarks. CSV, protein-identity GMT, PNG and SVG exports accompany the report.",
        y,
    )
    finish("Overview and stringent shortlist", "overview")

    methods = [
        (
            "Protein identity",
            "The genome-wide oncoref registry groups genes with a byte-identical longest Ensembl 112 protein sequence. Every selected group was checked against the cached GRCh38 peptide FASTA. Tied longest isoforms use a shared identical sequence. This operational protein identity does not resolve tumor-specific isoforms or post-translational modifications.",
        ),
        (
            "Recalculated expression and thresholds",
            "Within each retained patient group, member-gene clean expression is summed. Known repeat specimens retain the original equal-weight linear averaging. Identical-protein genes across the biological transcriptome are collapsed before recalculating each patient's p70 and p90; noncoding and ungrouped biological rows remain background entries. Zeros remain in the quantile calculation. Ties at the quantile do not count as positive.",
        ),
        (
            "Patient and cohort scope",
            "The original 142 source-matrix hashes, sample QC and donor assignments are retained. The report uses only the 100 cohorts with at least 20 patient groups. Existing QC failures and negative-TPM NUTM samples remain excluded. Some cohorts have unverified donor identities, so patient groups are not universally confirmed independent people.",
        ),
        (
            "Missing loci and assay scale",
            "An entry has all member-gene data available when every registered member gene has an expression value in every patient group in that cohort. Recorded zeros count as available data; missing values do not. Only entries with all these values can establish primary selection and coverage. Partial member sums are retained and flagged in the heatmaps and patient plots. Non-TPM cohorts are eligible for prevalence. TPM means and medians are populated only for linear-TPM sources; proxy-scale means, when exported, retain their own labeled units.",
        ),
        (
            "Selection versus coverage",
            "Selection requires strictly >50% or >70% prevalence above p70 or p90 in at least one eligible cohort. The separate coverage question uses strictly >10% patient prevalence at each of p70 and p90. Equality is excluded in both cases. Mean expression is computed only among positive patients for the corresponding percentile.",
        ),
        (
            "Coverage denominator",
            f"Raw coverage is the number of qualifying cohort views divided by all {n_views} eligible views. Adjusted coverage counts each of the {n_types} cancer-type groups once if any eligible member cohort qualifies. Molecular/risk subtypes and evidence-scope views map to their registry parent; cohorts sharing known patients also merge transitively. Distinct histologic types remain separate unless linked by those rules. A subtype hit supports its group, not every patient in the parent cancer type.",
        ),
        (
            "Unavailable and overlapping data",
            "The fixed denominator retains unmeasured entries, so coverage is observed support, not evidence of absence in unevaluable types. CSVs also give measured-only coverage and partial-member sensitivity counts. The group map documents every merge. This adjustment removes specified duplicate support; unknown cross-source patient overlap may remain.",
        ),
        (
            "Protein length and normal-tissue evidence",
            "Length is the full longest annotated protein in amino acids, including signal peptides and propeptides and excluding terminal stop markers. Normal somatic RNA is the largest single-member HPA nTPM annotation, with gene and tissue named; it is not the summed proteoform value. Member restriction confidence and incomplete evidence are exported. RNA abundance and reference restriction are not direct peptide-presentation or clinical-safety measurements.",
        ),
        (
            "Validation",
            "All source hashes, 142 cohort patient counts and per-patient expression mass conservation were checked. Full-data validation checks complete member coverage, exact fractions, nested selections and largest-cohort support. Twelve numerical tests include NY-ESO-1 summation and recalculated cutoffs, fixed coverage denominators, subtype grouping and retention of low heatmap frequencies.",
        ),
    ]
    y = header("Definitions and interpretation")
    for title, text in methods:
        block = f"<b>{title}.</b> {escape(text)}"
        obj = Paragraph(block, small)
        _, h = obj.wrap(540, 2000)
        if y - h < 48:
            finish("Definitions and interpretation", "methods")
            y = header("Definitions and interpretation (continued)")
        y = paragraph(block, y, style=small)
    finish("Definitions and interpretation", "methods")

    # Mapping provenance makes the distinction from a gene-level shortlist explicit.
    multi = universe[universe.n_member_genes.gt(1)]
    for start in range(0, len(multi), 12):
        part = multi.iloc[start : start + 12]
        y = header(
            "Identical-protein CTA groups",
            "Only groups touching the original CTA universe are included. All listed member genes contribute when measured.",
        )
        rows = [["Proteoform", "Identical-protein member genes", "aa", "Additional loci"]]
        for _, r in part.iterrows():
            rows.append(
                [
                    r.Symbol,
                    r.member_symbols.replace(";", ", "),
                    int(r.protein_length_aa),
                    r.additional_member_symbols if pd.notna(r.additional_member_symbols) else "--",
                ]
            )
        table(rows, [105, 285, 40, 110], y)
        finish(f"Identical-protein groups {start + 1}-{start + len(part)}", "identity")

    wide = landscape(A3)
    widths = [112, 40, 45, 102, 102, 115, 143, 67, 106, 144, 80]
    # Table fits within the wide page's 1118pt content area.
    for f, p in FOCUSED:
        data = lists[f, p]
        for start in range(0, len(data), 19):
            part = data.iloc[start : start + 19]
            y = header(
                f"Proteoform shortlist: >{f}% prevalence above p{p}",
                f"Complete list: {len(data)} entries; rows {start + 1}-{start + len(part)}. Coverage columns use >10% positive patients, separately from the selection threshold.",
                size=wide,
            )
            y = paragraph(
                f"Type coverage uses {n_types} cancer-type groups; cohort coverage uses {n_views} cohort views. Mean TPM refers to the largest passing cohort and is blank for a proxy-scale source. HPA nTPM is the largest single-member normal somatic value, with its gene/tissue, and is not a proteoform sum.",
                y,
                width=wide[0] - 72,
                style=small,
            )
            rows = [
                [
                    "Proteoform",
                    "aa",
                    "Member loci",
                    ">10% at p70: cancer types",
                    ">10% at p90: cancer types",
                    ">10% cohort views: p70 / p90",
                    "Largest passing cohort(s) / n",
                    "Positive %",
                    "Positive mean TPM",
                    "Normal somatic RNA: max member nTPM / tissue",
                    "HPA confidence",
                ]
            ]
            for _, r in part.iterrows():
                mean = (
                    f"{r.largest_cohort_positive_mean_tpm:.2f}"
                    if pd.notna(r.largest_cohort_positive_mean_tpm)
                    else "-- (proxy)"
                )
                normal = (
                    f"{r.normal_rna_max_member_somatic_ntpm:g}; {r.normal_rna_max_member_and_tissue}"
                    if pd.notna(r.normal_rna_max_member_somatic_ntpm)
                    else "Unavailable"
                )
                rows.append(
                    [
                        r.Symbol,
                        int(r.protein_length_aa),
                        int(r.n_member_genes),
                        frac(r, 70),
                        frac(r, 90),
                        f"{int(r.p70_n_cohort_views_gt10)}/{n_views} / {int(r.p90_n_cohort_views_gt10)}/{n_views}",
                        f"{r.largest_passing_cohorts}; n={int(r.largest_passing_cohort_size)}",
                        f"{r.largest_cohort_prevalence:.1%}",
                        mean,
                        normal,
                        r.member_restriction_confidence
                        if pd.notna(r.member_restriction_confidence)
                        else "--",
                    ]
                )
            table(rows, widths, y)
            finish(
                f">{f}% / p{p}: shortlist rows {start + 1}-{start + len(part)}",
                "shortlists",
                size=wide,
            )

    strict_pairs = pd.read_csv(
        out / "selections/min20/prevalence_gt70_transcriptome_p90/passing_proteoform_cohorts.csv"
    )
    strict_pairs = strict_pairs.sort_values(["Symbol", "n_patients"], ascending=[True, False])
    for start in range(0, len(strict_pairs), 18):
        part = strict_pairs.iloc[start : start + 18]
        y = header(
            "Stringent shortlist: passing cancer cohorts",
            f">70% / p90; cohort size >=20. Rows {start + 1}-{start + len(part)} of {len(strict_pairs)}. Broader lower-frequency coverage appears in the heatmaps.",
        )
        rows = [
            ["Protein", "aa", "Cancer cohort", "Positive / n", "Positive %", "Mean positive TPM"]
        ]
        for _, r in part.iterrows():
            rows.append(
                [
                    r.Symbol,
                    int(r.protein_length_aa),
                    f"{r.cancer_name} ({r.cancer_code})",
                    f"{r.n_expressing}/{r.n_patients}",
                    f"{r.fraction_expressing:.1%}",
                    f"{r.mean_expressing_tpm:.2f}"
                    if pd.notna(r.mean_expressing_tpm)
                    else "-- (proxy)",
                ]
            )
        table(rows, [60, 36, 226, 68, 60, 90], y)
        finish(f"Stringent passing cancer cohorts {start + 1}-{start + len(part)}", "stringent")

    mortality_dir = out / "mortality"
    if (mortality_dir / "compact_list_coverage_summary.csv").exists():
        mortality = pd.read_csv(mortality_dir / "compact_list_coverage_summary.csv")
        top = pd.read_csv(mortality_dir / "top10_world_mortality.csv")
        y = header("Compact lists across leading mortality categories")
        y = paragraph(
            "Mortality ranking comes directly from <b>oncoref.cancer_burden_df()</b>: curated GLOBOCAN 2022 global mortality shares. Cohorts use <b>oncoref.burden_category()</b>. The residual other/unknown bucket is excluded from this ranking of named cancers. Internal absolute death-count fields are unavailable.",
            y,
        )
        y = paragraph(
            "Coverage below is the <b>best observed cohort</b> in each category, including molecular subtypes. It counts patients expressing <b>at least one</b> listed protein above their own p90. Supporting cohort and patient count are shown. Overlapping cohorts are never pooled; these percentages are not worldwide patient prevalence. A + marks a lower bound when some panel proteins lack complete data.",
            y,
        )
        rows = [
            [
                "Rank",
                "Cancer category",
                "Global deaths share",
                "13-protein list: best cohort ANY coverage",
                "3-protein list: best cohort ANY coverage",
            ]
        ]
        for r in top.itertuples():
            cells = []
            for f in (50, 70):
                q = mortality[
                    mortality.prevalence_gt_pct.eq(f)
                    & mortality.burden_category.eq(r.burden_category)
                ].iloc[0]
                cells.append(
                    f"{q.best_panel_fraction:.1%}{'+' if q.best_panel_is_lower_bound else ''}; {q.best_panel_cohort}; n={int(q.best_panel_cohort_n)}"
                )
            rows.append([r.mortality_rank, r.cancer_label, f"{r.world_mortality_pct:g}%", *cells])
        y = table(rows, [28, 85, 57, 185, 185], y)
        paragraph(
            "Head/neck includes salivary cohorts under the current oncoref mapping. Its best result is ADCC (adenoid cystic carcinoma), not the broad HNSC cohort; exact salivary inclusion in the mortality aggregate needs review.",
            y,
            style=small,
        )
        finish("Top mortality categories: compact-list patient coverage", "mortality_summary")

        y = header("Mortality reference compatibility check")
        y = paragraph(
            "The already retrieved external IARC figures describe <b>2024</b> (published 8 July 2026), while oncoref cites <b>2022</b>. The top-five order agrees. The table below compares like broad categories where possible; differences between years alone do not prove an error. External values are an audit comparison only and do not replace the internal ranking.",
            y,
        )
        comparison = pd.read_csv(mortality_dir / "mortality_reference_comparison.csv")
        rows = [
            [
                "Category",
                "oncoref 2022 share",
                "External 2024 share",
                "Difference (percentage points)",
            ]
        ]
        for r in comparison.itertuples():
            rows.append(
                [
                    r.burden_category.replace("_", " "),
                    f"{r.internal_world_mortality_pct:.1f}%",
                    f"{r.external_world_mortality_pct:.2f}%",
                    f"{r.difference_percentage_points:+.2f}",
                ]
            )
        y = table(rows, [145, 125, 125, 145], y)
        y = paragraph(
            "<b>Category definitions:</b> oncoref groups oral cavity, pharynx, larynx and nasopharynx as head/neck and splits AML from other leukemia. GLOBOCAN's site table separates the head/neck sites and combines leukemia. The external head/neck comparison sums lip/oral cavity, oropharynx, hypopharynx, nasopharynx and larynx.",
            y - 10,
            style=small,
        )
        y = paragraph(
            "<b>Mapping to review:</b> oncoref's liver burden includes intrahepatic bile duct, but CHOL maps to gallbladder_biliary. This report follows the current API mapping, so CHOL is outside the liver coverage summary. The issue is recorded for reference-data review.",
            y,
            style=small,
        )
        y = paragraph(
            'Head/neck mortality source composition is also unresolved for salivary cancers. The external comparison sums five sites and excludes salivary gland; including it adds 0.21 percentage points. Follow-ups: <link href="https://github.com/pirl-unc/oncoref/issues/542" color="#127C80">#542: mortality refresh</link> and <link href="https://github.com/pirl-unc/oncoref/issues/543" color="#127C80">#543: site mapping</link>.',
            y,
            style=small,
        )
        paragraph(
            'External audit source: Sung et al., Global cancer statistics 2024, Table 1. <link href="https://doi.org/10.3322/caac.70090" color="#127C80">DOI:10.3322/caac.70090</link>. Internal source: oncoref cancer-incidence-mortality.csv, citing GLOBOCAN 2022 (DOI:10.3322/caac.21834). Full provenance and differences are in mortality/mortality_reference_comparison.csv.',
            y,
            style=small,
        )
        finish("Internal and external mortality reference comparison", "mortality_audit")

    hpa_dir = out / "hpa_tissues"
    if (hpa_dir / "protein_tissue_summary.csv").exists():
        hpa_summary = pd.read_csv(hpa_dir / "protein_tissue_summary.csv")
        hpa_links = pd.read_csv(hpa_dir / "hpa_gene_links.csv")
        y = header("HPA normal-tissue expression for the compact lists")
        y = paragraph(
            "The 13-protein p90 shortlist contains the strictest three entries and maps to 16 genes. All have HPA v23 RNA consensus measurements for 50 tissues, loaded through oncoref with verified source checksums. The following figures are redrawn from those HPA tables; original gene pages are linked below.",
            y,
        )
        y = paragraph(
            "<b>RNA:</b> identical-protein member genes are stacked, giving summed tissue-reference nTPM. This is an RNA proxy, not measured protein abundance. Reproductive tissues appear on the left, including breast; other tissues, including thymus, appear on the right. Each protein page uses the same RNA axis on both sides.",
            y,
        )
        y = paragraph(
            "<b>Protein staining (IHC):</b> five of the 13 entries have scored tissue observations, all with Enhanced reliability. The overview takes the strongest reported member-gene/cell-type level per tissue; individual pages retain member genes. Levels are never summed. Gray means unavailable; white means a measured Not detected result. IHC and RNA retain their different tissue vocabularies.",
            y,
        )
        y = paragraph(
            "The older shortlist normal-tissue column is the maximum single-member somatic nTPM. These new plots show member sums and all tissues, including breast and thymus, so their maxima can differ.",
            y,
            style=small,
        )
        for r in hpa_summary.itertuples():
            links = hpa_links[hpa_links[KEY].eq(getattr(r, KEY))]
            rendered = ", ".join(
                f'<link href="{q.hpa_tissue_url}" color="#127C80">{q.Symbol}</link>'
                for q in links.itertuples()
            )
            evidence = (
                f"IHC: {r.n_ihc_tissues_with_scored_data} tissues"
                if r.n_ihc_tissues_with_scored_data
                else "IHC unavailable"
            )
            y = paragraph(
                f"<b>{r.Symbol}</b> ({r.protein_length_aa} aa): {rendered}. {evidence}.",
                y,
                style=small,
            )
        finish("HPA tissue evidence: methods and original source pages", "hpa_summary")

    for category in [
        "overview",
        "mortality",
        "hpa",
        "selection",
        "heatmaps",
        "patients",
        "appendix",
    ]:
        for _, r in plots[plots.category.eq(category)].iterrows():
            path = out / "plots" / f"{r['name']}.png"
            with Image.open(path) as im:
                iw, ih = im.size
            size = landscape(A3) if iw / ih > 1.12 else A3
            w, h = size
            pdf.setPageSize(size)
            scale = min((w - 42) / iw, (h - 62) / ih)
            rw, rh = iw * scale, ih * scale
            pdf.drawImage(str(path), (w - rw) / 2, (h - rh) / 2 + 8, width=rw, height=rh)
            finish(r.title, category, size=size, figure=r["name"])
    pdf.save()
    pd.DataFrame(page_index).to_csv(out / "pdf_page_index.csv", index=False)
    (out / "methodology.md").write_text(
        "# Proteoform analysis methods\n\n"
        + "\n\n".join(f"**{t}.** {b}" for t, b in methods)
        + "\n"
    )
    readme = [
        "# CTA proteoform shortlists",
        "",
        f"[Complete PDF: {page} pages, {len(plots)} figures](cta-analysis-all-figures.pdf)",
        "",
        "Four settings only: >50% or >70% prevalence, p70 or p90; minimum cohort size 20. All expression scales are eligible; TPM means use only linear-TPM cohorts.",
        "",
        "[Combined shortlist evaluation table](shortlist_evaluation.csv) | [All protein identities and lengths](proteoform_universe.csv) | [Gene membership and normal-tissue annotations](proteoform_member_annotations.csv)",
        "",
        "| Selection | Proteoforms | Full shortlist |",
        "|---|---:|---|",
    ]
    for f, p in FOCUSED:
        readme.append(
            f"| >{f}% / p{p} | {len(lists[f, p])} | [Evaluation table](selections/min20/prevalence_gt{f}_transcriptome_p{p}/proteoforms_ranked.csv) |"
        )
    readme += [
        "",
        f"Coverage uses >10% positive patients at each percentile, across all {n_views} eligible cohort views and {n_types} subtype/overlap-adjusted cancer-type groups. Fixed denominators retain unavailable entries. Measured-only fractions and partial-member counts are available in CSV.",
        "",
        "| Stringent protein | aa | >10% p70 type coverage | >10% p90 type coverage |",
        "|---|---:|---:|---:|",
    ]
    for _, r in stringent.iterrows():
        readme.append(
            f"| {r.Symbol} | {int(r.protein_length_aa)} | {frac(r, 70)} | {frac(r, 90)} |"
        )
    readme += [
        "",
        "Heatmap columns enter when any selected proteoform exceeds 10% prevalence. All values in included columns are retained; low nonzero values are displayed, zeros remain zeros, incomplete sums are marked, and unavailable cells are gray. Asterisks and orange borders mark selection-passing cells.",
        "",
        "## Data and provenance",
        "",
        "- [Coverage for all proteoforms](proteoform_coverage_gt10.csv)",
        "- [Cancer-type and overlap group mapping](cohort_overlap_groups.csv)",
        "- [Every displayed heatmap cell](shortlist_heatmap_all_displayed_cells.csv.gz)",
        "- [All cohort-level measurements](all_proteoform_cohort_metrics.csv.gz)",
        "- [Gene-to-proteoform changes](gene_vs_proteoform_selection.csv)",
        "- [Patient-level case studies](case_study_patient_expression.csv.gz)",
        "- [Protein sequence and HPA provenance](annotation_manifest.json)",
        "- [Definitions and limitations](methodology.md)",
        "- [Validation](validation.json)",
        "",
        "## Reproduce from this checkout",
        "",
        "```sh",
        "CANCERDATA_PER_SAMPLE_CACHE=1 .venv/bin/python scripts/cta_proteoform_report.py",
        "MPLCONFIGDIR=$PWD/tmp/cta_threshold_report/matplotlib .venv/bin/python scripts/plot_cta_proteoform_report.py",
        "MPLCONFIGDIR=$PWD/tmp/cta_threshold_report/matplotlib .venv/bin/python scripts/cta_mortality_coverage_report.py",
        "MPLCONFIGDIR=$PWD/tmp/cta_threshold_report/matplotlib .venv/bin/python scripts/cta_hpa_tissue_report.py",
        "python scripts/render_cta_proteoform_report.py",
        "```",
        "",
        f"The renderer requires reportlab, pandas and Pillow; plotting also requires matplotlib and adjustText. Matrix caches and the original QC/donor audit must be available; use --resume only when those inputs and the numerical calculation are unchanged. The final PDF retains all {len(plots)} figures; checkpoints and raw matrices are excluded from the download bundle.",
        "",
        "The identity registry uses exact longest-protein sequence identity, not isoform-resolved proteomics. Normal HPA values are member-level reference annotations, not a direct protein-level safety assessment.",
    ]
    readme += [
        "",
        "## Readable text and mortality coverage",
        "",
        "[All four cutoff reports in plain text](all_cutoff_reports.txt). Each selection folder also contains shortlist_report.txt, proteoform_symbols.txt and proteoform_keys.txt alongside proteoforms_ranked.csv.",
        "[Top mortality report](mortality/report.md) | [Mortality coverage CSV](mortality/compact_list_coverage_summary.csv) | [Reference compatibility audit](mortality/mortality_reference_comparison.csv)",
        "[HPA tissue report](hpa_tissues/report.md) | [HPA tissue summary CSV](hpa_tissues/protein_tissue_summary.csv) | [Original HPA gene pages](hpa_tissues/hpa_gene_links.csv)",
    ]
    (out / "report.md").write_text("\n".join(readme) + "\n")
    repro = out / "reproducibility"
    repro.mkdir(exist_ok=True)
    for filename in [
        "cta_proteoform_report.py",
        "cta_mortality_coverage_report.py",
        "cta_hpa_tissue_report.py",
        "plot_cta_proteoform_report.py",
        "render_cta_proteoform_report.py",
        "cta_report_protein_lengths.py",
        "cta_threshold_report_details.py",
        "cta_threshold_report.py",
    ]:
        shutil.copy2(ROOT / "scripts" / filename, repro / filename)
    for filename in [
        "test_cta_proteoform_report.py",
        "test_cta_mortality_coverage_report.py",
        "test_cta_hpa_tissue_report.py",
        "test_cta_threshold_report.py",
        "test_cta_threshold_report_details.py",
    ]:
        shutil.copy2(ROOT / "tests" / filename, repro / filename)
    heatmaps = pd.read_csv(out / "shortlist_heatmap_all_displayed_cells.csv.gz")
    validation = json.loads((out / "validation.json").read_text())
    validation["n_subthreshold_nonzero_heatmap_cells"] = int(
        (
            heatmaps.available_sum_complete_measurement
            & heatmaps.n_expressing.gt(0)
            & ~heatmaps.passes_selection
        ).sum()
    )
    for (f, p), selected in lists.items():
        scenario = f"prevalence_gt{f}_transcriptome_p{p}"
        shown = heatmaps[heatmaps.scenario.eq(scenario)]
        assert len(shown) == len(selected) * shown.cancer_code.nunique()
        assert set(shown[KEY]) == set(selected[KEY])
        assert set(shown.cancer_code) == set(
            shown.loc[
                shown.available_sum_complete_measurement
                & (100 * shown.n_expressing > 10 * shown.n_patients),
                "cancer_code",
            ]
        )
    validation["heatmap_full_grid_and_column_rule_verified"] = True
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    archive = out / "cta-proteoform-analysis-bundle.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(out.rglob("*")):
            if (
                path.is_file()
                and path != archive
                and not {"checkpoints", "qa"} & set(path.parts)
                and path.suffix != ".log"
            ):
                z.write(path, path.relative_to(out))
    print(
        f"Created {page}-page PDF, {len(plots)} figures, {archive.stat().st_size / 1e6:.1f} MB bundle"
    )


if __name__ == "__main__":
    main()
