"""Figures and readable documentation for the reproducible CTA scope audit."""

from __future__ import annotations

import textwrap
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd

from oncoref import figure_style
from oncoref.normal_tissues import normal_tissue_taxonomy

SCOPES = ("reproductive", "reproductive_or_disposable")
LABELS = ("Reproductive", "Reproductive + disposable")
FIGURES = {
    "source-intake": "Published nomination intake",
    "filtering-stages": "Source intake through tissue filtering",
    "source-retention": "Retained candidates by source",
    "evidence-coverage": "Observed compatibility and missing evidence",
}
CAPTIONS = {
    "source-intake": "Each paper contributes its registered tables or targeted selection. Aliases collapse to canonical gene IDs; unmapped identities remain in the intake. Nested tables are deduplicated within each paper. Sources overlap and cannot be summed. Bai 2016 is a targeted EGFL6 nomination; Bradley 2020 contributes the ten candidates in Figure 3.",
    "filtering-stages": "Both scopes use the same source universe. The adaptive policy retains the current reliability-dependent RNA thresholds, expression rule, nomination holdbacks and reviews. The strict sensitivity screen requires complete routine bulk RNA, positive evidence in scope, and no positive outside RNA/IHC observation. Strict compatibility remains conditional on measured coverage.",
    "source-retention": "Green shows candidates retained with reproductive tissues alone. Blue shows additional candidates admitted by the disposable compartments. Each candidate may occur in several sources. Both methods preserve family exclusions and nomination/review holdbacks. The strict screen also retains separate missing-IHC coverage flags in the candidate table.",
    "evidence-coverage": "Counts cover all protein-coding nominations before family or review exclusions. RNA positivity means any reported estimate above zero. Raw IHC positivity includes Low, Medium, High, Ascending and Descending. The joint strict screen requires complete routine bulk RNA; unresolved mapping, unusable measurements or reviewed unreliable IHC blocks compatibility. Missing IHC is recorded separately and never establishes protein absence.",
}


def save_figure(fig, out, stem, caption):
    import matplotlib.pyplot as plt

    fig.text(
        0.02, 0.02, textwrap.fill(caption, width=133), fontsize=8, va="bottom", linespacing=1.35
    )
    for suffix in ("png", "pdf"):
        figure_style.save(fig, out / f"{stem}.{suffix}", keep_titles=True)
    plt.close(fig)


def draw_figures(out, candidates, sources, summary):
    figure_style.apply()
    import matplotlib.pyplot as plt

    source = sources.loc[sources.scope.eq(SCOPES[0])].sort_values(
        "source_identities", ascending=False
    )
    y = np.arange(len(source))
    fig, ax = plt.subplots(figsize=(11, 8.5))
    for offset, column, label, color in (
        (-0.25, "source_identities", "Published identities", "#9ca3af"),
        (0, "mapped", "Mapped gene IDs", "#377eb8"),
        (0.25, "protein_coding", "Protein-coding", "#237a57"),
    ):
        bars = ax.barh(y + offset, source[column], height=0.23, color=color, label=label)
        ax.bar_label(bars, padding=3, fontsize=7.8)
    ax.set_yticks(y, source.source)
    ax.invert_yaxis()
    ax.set_xlabel("Unique identities within each paper")
    ax.set_xlim(0, source.source_identities.max() * 1.16)
    ax.set_title("Human CTA nomination sources | registered source selections")
    ax.legend(loc="lower right")
    fig.subplots_adjust(left=0.2, right=0.98, top=0.93, bottom=0.19)
    save_figure(fig, out, "source-intake", CAPTIONS["source-intake"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.7))
    for ax, policy, title in zip(
        axes,
        ("adaptive", "strict"),
        ("Current adaptive policy", "Strict positive-evidence sensitivity"),
    ):
        columns = [
            "source_identities",
            "mapped",
            "protein_coding",
            "family_eligible",
            f"{policy}_after_family",
            f"{policy}_retained",
        ]
        labels = [
            "Published identities",
            "Mapped genes",
            "Protein-coding",
            "Family eligible",
            "Tissue-scope screen",
            "Final review rules",
        ]
        for j, (scope, name, color) in enumerate(zip(SCOPES, LABELS, ("#237a57", "#377eb8"))):
            values = summary.set_index("scope").loc[scope, columns].astype(int)
            bars = ax.barh(
                np.arange(len(columns)) + (j - 0.5) * 0.36,
                values,
                height=0.34,
                color=color,
                label=name,
            )
            ax.bar_label(bars, padding=3, fontsize=8)
        ax.set_yticks(np.arange(len(columns)), labels)
        ax.invert_yaxis()
        ax.set_xlabel("Unique identities / genes")
        ax.set_xlim(0, summary.source_identities.max() * 1.19)
        ax.set_title(title)
    axes[1].legend(loc="lower right", fontsize=8)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.91, bottom=0.23, wspace=0.66)
    save_figure(fig, out, "filtering-stages", CAPTIONS["filtering-stages"])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 8))
    order = source.doi.tolist()
    for ax, policy, title in zip(
        axes,
        ("adaptive", "strict"),
        ("Current adaptive policy", "Strict positive-evidence sensitivity"),
    ):
        values = sources.pivot(index="doi", columns="scope", values=f"{policy}_retained").loc[order]
        base = values[SCOPES[0]]
        gain = values[SCOPES[1]] - base
        ax.barh(y, base, height=0.64, color="#237a57", label="Reproductive")
        ax.barh(y, gain, left=base, height=0.64, color="#377eb8", label="Additional disposable")
        extent = max(values[SCOPES[1]].max(), 1)
        for k, (n, extra) in enumerate(zip(base, gain)):
            ax.text(n + extra + extent * 0.015, k, f"{n:,} + {extra:,}", va="center", fontsize=8)
        ax.set_yticks(y, source.source)
        ax.invert_yaxis()
        ax.set_xlim(0, extent * 1.3)
        ax.set_title(title)
        ax.set_xlabel("Retained genes (reproductive + additional disposable)")
    axes[1].legend(loc="lower right", fontsize=8)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.93, bottom=0.2, wspace=0.6)
    save_figure(fig, out, "source-retention", CAPTIONS["source-retention"])

    statuses = [
        ("compatible", "Compatible observed pattern", "#237a57"),
        ("outside_scope", "Outside-scope positive", "#c05a52"),
        ("no_positive_observation", "No positive observation", "#a1a1a1"),
        ("unresolved", "Unresolved", "#dfaf36"),
        ("unavailable", "Unavailable", "#e5e7eb"),
    ]
    labels, counts = [], []
    for modality, label in (
        ("rna", "Bulk RNA"),
        ("ihc", "Raw IHC"),
        ("strict", "Joint strict screen"),
    ):
        for scope, scope_label in zip(SCOPES, LABELS):
            subset = candidates.loc[candidates.scope.eq(scope), f"{modality}_status"]
            labels.append(f"{label}\n{scope_label}")
            counts.append(subset.value_counts().to_dict())
    fig, ax = plt.subplots(figsize=(12, 6.8))
    left = np.zeros(len(labels))
    for key, label, color in statuses:
        values = np.array([row.get(key, 0) for row in counts])
        bars = ax.barh(
            np.arange(len(labels)), values, left=left, height=0.6, color=color, label=label
        )
        for bar, n in zip(bars, values):
            if n >= 60:
                ax.text(
                    bar.get_x() + n / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{n:,}",
                    ha="center",
                    va="center",
                    color="white" if key in {"compatible", "outside_scope"} else "#222222",
                    fontsize=8,
                )
        left += values
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, left.max())
    ax.set_xlabel("Protein-coding candidates (exact counts in evidence-status-counts.csv)")
    ax.set_title("Human HPA v23 | compatibility of the observations actually available")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.15), ncol=3, fontsize=8)
    fig.subplots_adjust(left=0.23, right=0.98, top=0.91, bottom=0.31)
    save_figure(fig, out, "evidence-coverage", CAPTIONS["evidence-coverage"])
    pd.DataFrame(
        [
            {"view": label, **{key: row.get(key, 0) for key, _, _ in statuses}}
            for label, row in zip(labels, counts)
        ]
    ).to_csv(out / "evidence-status-counts.csv", index=False)


def report_text(summary, candidates, sources, gains):
    indexed = summary.set_index("scope")
    r, d = indexed.loc[SCOPES[0]], indexed.loc[SCOPES[1]]
    n_papers = sources.doi.nunique()
    ihc_missing = int(r.strict_retained_without_ihc)
    added = gains.loc[gains.policy.eq("adaptive"), "Symbol"].sort_values().tolist()
    gain_note = f" The adaptive additions are {', '.join(added)}." if added else ""
    zero_denominator = gains.loc[gains.policy.eq("adaptive") & gains.rna_deflated_denominator_zero]
    floor_note = " ".join(
        f"{row.Symbol} has maximum RNA {row.rna_max_ntpm:g} nTPM and a zero deflated RNA denominator; its adaptive admission relies on the existing zero-denominator convention and {row.protein_reliability} in-scope IHC."
        for row in zero_denominator.itertuples(index=False)
    )
    return {
        "overview": f"This human analysis starts with {n_papers} registered published nomination sources: {int(r.source_identities):,} distinct source identities and {int(r.protein_coding):,} mapped protein-coding candidates. It evaluates reproductive tissues and their union with cutaneous melanocytes and hair-follicle compartments using the same candidate universe.",
        "results": f"The adaptive reproductive baseline reproduces all {int(r.adaptive_retained):,} genes in the published extended panel. Adding disposable compartments retains {int(d.adaptive_retained):,} genes ({int(d.adaptive_retained - r.adaptive_retained):,} additional). The strict sensitivity screen retains {int(r.strict_retained):,} and {int(d.strict_retained):,}, respectively ({int(d.strict_retained - r.strict_retained):,} additional)."
        + gain_note
        + (" " + floor_note if floor_note else ""),
        "coverage": f"Strict means no outside positive observation under the stated assay rules. Adequate evidence of absence is not established. Of the {int(r.strict_retained):,} reproductive strict-screen candidates, {ihc_missing:,} lack usable IHC observations; {int(r.strict_retained_partial_ihc):,} have partial routine IHC coverage. A routine context is represented for at least half of the genes in that source; its exact tissue/cell list is recorded in the manifest. These coverage distinctions remain explicit in every candidate record.",
        "adaptive": "The adaptive analysis reproduces the existing extended CTA policy. RNA uses the sum of max(nTPM - 1, 0) inside scope divided by that sum across all measured tissues, rounded to four decimals. A zero denominator receives fraction 1 by the existing policy; this is not positive RNA support, and the expression/IHC rules still apply. Thresholds are 0.80 (Enhanced IHC), 0.90 (Supported), 0.95 (Approved), and 0.97 (Uncertain or missing). Low/Medium/High IHC detections must all map inside scope. Existing cross-reactive-antibody exclusions, expression/rescue rules, family exclusions, normal-only nomination holdbacks and extended-scope manual reviews remain in force. The disposable comparison changes only the tissue allowance; it does not replace the public CTA panel.",
        "strict": "The strict sensitivity screen uses every positive bulk RNA estimate (>0 nTPM) and raw IHC Low/Medium/High/Ascending/Descending. Any mapped outside positive blocks compatibility. Missing routine bulk RNA, unresolved positive mapping, unusable measurements, or reviewed unreliable IHC prevents a compatible result. A positive in-scope observation is required. Missing IHC remains explicitly RNA-only evidence; partial IHC remains partial. Family and nomination/review exclusions also apply. This screen does not use the adaptive RNA fraction or low-expression rescue, so the two methods need not form nested gene sets.",
        "anatomy": "Reproductive includes male and female reproductive anatomy (including mammary), placenta and other extraembryonic compartments, and explicitly embryonic/fetal contexts. The combined scope additionally admits cutaneous melanocytes and hair-follicle compartments. Whole skin, unsited melanocytes, ocular pigment and retinal pigment epithelium do not qualify. HPA bulk RNA has no separate disposable-compartment measurement, so its allowed tissue set is identical for the two scopes. Only resolved IHC tissue/cell pairs can add the disposable allowance in this report.",
        "sources": "Source rows preserve DOI, URL, file/sheet/row location, original identifier, mapping method and assay resolution. Nested tables are deduplicated within a paper; papers overlap. The intake uses every registered published set, without selecting papers based on the final retained panel. Targeted source selections remain labeled. Unmapped and noncoding nominations are retained in the intake audit and do not enter the protein-coding screen. Historical-only nominations remain in a separate export.",
        "species": "The numerical filtering comparison is human-only and uses pinned HPA v23 bulk RNA and normal-tissue IHC. The taxonomy also supports dog and mouse anatomy. Canine native expression and publication evidence requires its own assay-specific analysis; human thresholds and coverage do not transfer. No native mouse expression/CTA atlas is bundled. Developmental contexts absent from HPA are unassessed here. HPA single-cell RNA is not included in these RNA/IHC policy comparisons; unsited cell types cannot supply full anatomical coverage.",
    }


def write_pdf(out, text, summary, sources):
    from pypdf import PdfReader, PdfWriter
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        KeepTogether,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="ReportBody", fontName="Helvetica", fontSize=9.5, leading=13, spaceAfter=9
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReportTitle",
            fontName="Helvetica-Bold",
            fontSize=21,
            leading=25,
            textColor=colors.HexColor("#183b4e"),
            spaceAfter=16,
        )
    )
    styles.add(
        ParagraphStyle(
            name="ReportSmall", fontName="Helvetica", fontSize=8.5, leading=11, spaceAfter=6
        )
    )
    story = [
        Paragraph("CTA sources and normal-tissue scopes", styles["ReportTitle"]),
        Paragraph(
            f"Human HPA v23 | Taxonomy policy {normal_tissue_taxonomy()['policy_version']}",
            styles["ReportSmall"],
        ),
    ]
    for key in ("overview", "results", "coverage"):
        story.append(Paragraph(escape(text[key]), styles["ReportBody"]))
    data = [["Screen", "Reproductive", "+ disposable", "Additional"]]
    indexed = summary.set_index("scope")
    for field, label in (
        ("adaptive_retained", "Current adaptive policy"),
        ("strict_retained", "Strict sensitivity"),
    ):
        a, b = (int(indexed.loc[scope, field]) for scope in SCOPES)
        data.append([label, str(a), str(b), str(b - a)])
    table = Table(data, colWidths=[2.5 * inch, 0.95 * inch, 0.95 * inch, 0.8 * inch], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eff3")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("LINEBELOW", (0, -1), (-1, -1), 0.5, colors.lightgrey),
            ]
        )
    )
    story.extend([Spacer(1, 6), table, Spacer(1, 15)])
    for heading, key in (
        ("Tissue definitions and assay resolution", "anatomy"),
        ("Current adaptive policy", "adaptive"),
        ("Strict sensitivity screen", "strict"),
        ("Source intake and provenance", "sources"),
        ("Species and coverage boundaries", "species"),
    ):
        story.append(
            KeepTogether(
                [
                    Paragraph(heading, styles["Heading2"]),
                    Paragraph(escape(text[key]), styles["ReportBody"]),
                ]
            )
        )
    story.append(Paragraph("Sources and reproducibility", styles["Heading2"]))
    story.append(
        Paragraph(
            'HPA v23: <link href="https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip" color="blue">bulk RNA</link> and <link href="https://v23.proteinatlas.org/download/normal_tissue.tsv.zip" color="blue">normal-tissue IHC</link>. The companion CSVs, taxonomy bibliography and run manifest preserve source hashes, code hashes, policy definitions and candidate decisions.',
            styles["ReportBody"],
        )
    )
    for row in sources.drop_duplicates("doi").sort_values("source").itertuples(index=False):
        story.append(
            Paragraph(
                f'{escape(row.source)}: <link href="https://doi.org/{escape(row.doi)}" color="blue">doi:{escape(row.doi)}</link>',
                styles["ReportSmall"],
            )
        )
    story.append(
        Paragraph(
            'Tissue rationale: Blagosklonny 2003, doi:10.1038/sj.bjc.6601256; Johnson et al. 2009, <link href="https://pmc.ncbi.nlm.nih.gov/articles/PMC2929689/" color="blue">PMC2929689</link>; Purba et al. 2019, <link href="https://pubmed.ncbi.nlm.nih.gov/31512803/" color="blue">PMID31512803</link>; <link href="https://www.cancer.gov/about-cancer/treatment/side-effects/hair-loss" color="blue">NCI alopecia</link>. Complete anatomy and species references accompany tissue-taxonomy.json.',
            styles["ReportSmall"],
        )
    )
    intro = out / "report-methods.pdf"

    def footer(canvas, doc):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(42, 28, "oncoref | CTA source and tissue-scope audit")
        canvas.drawRightString(A4[0] - 42, 28, str(doc.page))

    SimpleDocTemplate(
        str(intro), pagesize=A4, rightMargin=42, leftMargin=42, topMargin=40, bottomMargin=43
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    writer = PdfWriter()
    writer.append(str(intro), outline_item="Results, methods and references")
    for stem, title in FIGURES.items():
        writer.append(str(out / f"{stem}.pdf"), outline_item=title)
    writer.add_metadata(
        {
            "/Title": "CTA sources: reproductive and reproductive + disposable",
            "/Author": "oncoref",
            "/Subject": "Human HPA v23 source intake, filtering and coverage",
        }
    )
    target = out / "cta-tissue-scope-report.pdf"
    writer.write(str(target))
    writer.close()
    if len(PdfReader(target).pages) != len(PdfReader(intro).pages) + len(FIGURES):
        raise ValueError("Report is missing a figure")
    intro.unlink()


def write_index(out, text, summary, sources, gains):
    lines = [
        "# CTA sources: reproductive and reproductive + disposable",
        "",
        f"**Human HPA v23 · taxonomy policy {normal_tissue_taxonomy()['policy_version']}**",
        "",
        "[Download the complete report (PDF)](cta-tissue-scope-report.pdf) · [Exact scope definitions](../../normal-tissue-scopes.md)",
        "",
        text["overview"],
        "",
        text["results"],
        "",
        "| Filtering method | Reproductive | Reproductive + disposable | Added by disposable |",
        "|---|---:|---:|---:|",
    ]
    indexed = summary.set_index("scope")
    for field, label in (
        ("adaptive_retained", "Current adaptive policy"),
        ("strict_retained", "Strict positive-evidence sensitivity"),
    ):
        a, b = (int(indexed.loc[scope, field]) for scope in SCOPES)
        lines.append(f"| {label} | {a:,} | {b:,} | {b - a:,} |")
    lines.extend(
        [
            "",
            text["coverage"],
            "",
            "## What changes when disposable tissues are included?",
            "",
            text["anatomy"],
            "",
        ]
    )
    if gains.empty:
        lines.append(
            "No additional genes pass either complete filtering policy with the disposable allowance in this source snapshot. This does not imply that the categories are equivalent: the allowance changes resolved tissue membership, while other measured expression can still block admission."
        )
    else:
        for policy, group in gains.groupby("policy"):
            lines.append(
                f"**{policy.title()} additions:** {', '.join(group.Symbol.sort_values())}."
            )
    lines.extend(["", "[Gene-level gains and supporting contexts](disposable-gains.csv)", ""])
    for stem, title in FIGURES.items():
        lines.extend(
            [
                f"## {title}",
                "",
                CAPTIONS[stem],
                "",
                f"[Vector PDF]({stem}.pdf)",
                "",
                f"![{title}]({stem}.png)",
                "",
            ]
        )
    lines.extend(["## Methods", ""])
    for heading, key in (
        ("Adaptive filtering", "adaptive"),
        ("Strict sensitivity", "strict"),
        ("Source intake", "sources"),
        ("Species and coverage", "species"),
    ):
        lines.extend([f"### {heading}", "", text[key], ""])
    lines.extend(
        [
            "## Audit files",
            "",
            "| File | Contents |",
            "|---|---|",
            "| [Candidate decisions](candidate-decisions.csv) | Both scopes and policies; RNA/IHC status, outside-tissue counterevidence, coverage, review actions, disposable contexts |",
            "| [Source counts](source-counts.csv) | Deduplicated intake and filtering counts for every paper |",
            "| [Scope summary](scope-summary.csv) | Exact overall counts and missing-IHC counts |",
            "| [Evidence status counts](evidence-status-counts.csv) | Exact segment counts behind the coverage figure |",
            "| [Source membership](source-membership.csv.gz) | Original identities, DOI, source table/sheet/rows, mapping method and assay resolution |",
            "| [Normal observations](normal-observations.csv.gz) | All candidate bulk RNA/IHC rows with original labels, anatomical mappings and citations |",
            "| [Publication registry](publication-sources.csv) | Registered source selections, URLs and original artifact hashes |",
            "| [Historical-only candidates](historical-only-candidates.csv) | Archived nominations outside the published-source intake |",
            "| [Tissue taxonomy](tissue-taxonomy.json) | Definitions, unions, species applicability and bibliography |",
            "| [Run manifest](run-manifest.json) | Source, code and output hashes; versions; exact routine coverage contexts |",
            "",
            "`adaptive_retained` is the existing adaptive policy evaluated under the named scope; only the reproductive baseline is the public extended panel. `strict_retained` is a sensitivity-screen result. Neither field claims adequate negative coverage or clinical safety. `outside_rna_max_ntpm` is zero only when RNA was measured without an outside positive; unavailable RNA remains blank.",
            "",
            "## Reproduce",
            "",
            "```sh",
            "python scripts/cta_scope_report.py --out outputs/cta-tissue-scopes-current",
            "```",
            "",
            "Requires the report dependencies and pinned HPA v23 bulk RNA and IHC files. The report checks source hashes, equality of the reproductive adaptive baseline with the public extended panel, and monotonicity when the tissue union expands. Historical audits remain archived under their original dates.",
            "",
            "## Citations",
            "",
            "Measurements: [HPA v23 bulk RNA](https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip) and [normal-tissue IHC](https://v23.proteinatlas.org/download/normal_tissue.tsv.zip). Anatomical and disposable-tissue rationales are cited in the [scope definitions](../../normal-tissue-scopes.md). Source nominations:",
            "",
        ]
    )
    for row in sources.drop_duplicates("doi").sort_values("source").itertuples(index=False):
        lines.append(
            f"- [{row.source}](https://doi.org/{row.doi}) — the exact imported selection is recorded in the publication registry."
        )
    (out / "index.md").write_text("\n".join(lines) + "\n")


def render(out, candidates, sources, summary, gains):
    out = Path(out)
    draw_figures(out, candidates, sources, summary)
    text = report_text(summary, candidates, sources, gains)
    write_pdf(out, text, summary, sources)
    write_index(out, text, summary, sources, gains)
