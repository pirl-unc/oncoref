#!/usr/bin/env python3
"""Build a telescoping CTA handoff guide from the verified, cached analysis."""

from __future__ import annotations

import io
import shutil
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd
from PIL import Image
from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/cta_proteoform_report_20260917"
DEST = OUT / "primary_panel"
W, H = landscape(letter)
TEAL = colors.HexColor("#127C80")
INK = colors.HexColor("#203746")
PALE = colors.HexColor("#EDF5F5")
MUTED = colors.HexColor("#506473")


def main():
    primary = pd.read_csv(DEST / "primary_proteoforms.csv")
    burden = pd.read_csv(DEST / "primary_burden_coverage.csv")
    primary = primary.merge(
        burden.drop(columns="Symbol"), on="proteoform_key", validate="one_to_one"
    )
    ranked = primary.sort_values("p90_n_cancer_type_groups_gt10", ascending=False, kind="stable")
    members = pd.read_csv(OUT / "proteoform_member_annotations.csv").set_index("Symbol")
    genes = sorted({g for s in primary.member_symbols for g in s.split(";")})
    ids = sorted({g for s in primary.member_gene_ids for g in s.split(";")})
    assert len(primary) == 13 and len(genes) == len(ids) == 16
    (DEST / "primary_gene_symbols.txt").write_text("\n".join(genes) + "\n")
    (DEST / "primary_ensembl_gene_ids.txt").write_text("\n".join(ids) + "\n")
    handoff = []
    for r in ranked.itertuples():
        for symbol in r.member_symbols.split(";"):
            gene_id = members.loc[symbol, "Ensembl_Gene_ID"]
            assert gene_id in r.member_gene_ids.split(";")
            assert members.loc[symbol, "protein_sequence_sha256"] == r.protein_sequence_sha256
            handoff.append(
                {
                    "protein_group": r.Symbol,
                    "gene_symbol": symbol,
                    "ensembl_gene_id": gene_id,
                    "protein_length_aa": r.protein_length_aa,
                    "proteoform_key": r.proteoform_key,
                    "protein_sequence_sha256": r.protein_sequence_sha256,
                    "cancer_type_groups_gt10_p90": r.p90_n_cancer_type_groups_gt10,
                    "cancer_type_groups_denominator": 69,
                    "world_incidence_pct_represented": r.world_incidence_pct_represented,
                    "world_mortality_pct_represented": r.world_mortality_pct_represented,
                }
            )
    pd.DataFrame(handoff).to_csv(DEST / "primary_gene_handoff.csv", index=False)
    stream = io.BytesIO()
    c = canvas.Canvas(stream, pagesize=(W, H))
    titles = []
    body = ParagraphStyle("body", fontName="Helvetica", fontSize=11, leading=15, textColor=INK)
    small = ParagraphStyle("small", parent=body, fontSize=9, leading=12)
    cell = ParagraphStyle("cell", parent=body, fontSize=9, leading=12)
    y = 0

    def text(s, x=36, width=720, style=body, gap=13):
        nonlocal y
        p = Paragraph(s, style)
        _, h = p.wrap(width, 1000)
        p.drawOn(c, x, y - h)
        y -= h + gap

    def start(title, takeaway, level):
        nonlocal y
        titles.append(title)
        c.setFillColor(TEAL)
        c.rect(0, H - 31, W, 31, fill=1, stroke=0)
        c.setFillColor(colors.white)
        c.setFont("Helvetica", 9)
        c.drawString(36, H - 20, f"CTA WORKING LIST  |  {level}")
        c.setFillColor(INK)
        c.setFont("Helvetica-Bold", 23)
        c.drawString(36, H - 67, title)
        y = H - 86
        text(takeaway, gap=17)

    def end():
        assert y >= 40, (titles[-1], y)
        c.setFillColor(MUTED)
        c.setFont("Helvetica", 8)
        c.drawString(
            36, 20, "Oncoref cached analysis | 17 September 2026 | RNA-based candidate selection"
        )
        c.drawRightString(W - 36, 20, str(len(titles)))
        c.showPage()

    def table(rows, widths, font=9):
        nonlocal y
        cs = ParagraphStyle("tablecell", parent=cell, fontSize=font, leading=font + 3)
        hs = ParagraphStyle(
            "tablehead", parent=cs, textColor=colors.white, fontName="Helvetica-Bold"
        )
        data = [
            [Paragraph(str(v), hs if i == 0 else cs) for v in row] for i, row in enumerate(rows)
        ]
        t = Table(data, colWidths=widths, hAlign="LEFT")
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), TEAL),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [PALE, colors.white]),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("TOPPADDING", (0, 0), (-1, -1), 6),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                ]
            )
        )
        _, h = t.wrap(720, 1000)
        t.drawOn(c, 36, y - h)
        y -= h + 14

    def figure(title, takeaway, path, footer, level="2 / CHECK THE EVIDENCE"):
        nonlocal y
        start(title, takeaway, level)
        bottom = 82
        with Image.open(path) as im:
            iw, ih = im.size
        scale = min(720 / iw, (y - bottom) / ih)
        c.drawImage(
            str(path),
            (W - iw * scale) / 2,
            bottom + (y - bottom - ih * scale) / 2,
            width=iw * scale,
            height=ih * scale,
        )
        y = 70
        text(footer, style=small, gap=0)
        end()

    start(
        "Start with these 13 protein groups",
        "A practical starting list for downstream work: compact enough to inspect individually, with strong-expression coverage in several common cancers.",
        "1 / CHOOSE AND TAKE THE LIST",
    )
    text(
        "<b>Recommended list: more than 50% of patients above their own transcriptome p90 in at least one cohort of 20 or more patient groups.</b> p90 means expression exceeds the 90th percentile of that patient's transcriptome. It does not mean 90% of patients express the gene."
    )
    text(
        "<b>The list</b><br/>PRAME, XAGE1A/B, MAGEA3, MAGEA6, NY-ESO-1, EGFL6, MAGEA2/B, MAGEC1, ZFP42, GFY, DPPA5, TRIML2, UTF1."
    )
    table(
        [
            ["Why carry this forward?", "What to keep in mind"],
            [
                "At the same p90 cutoff, this list covers more patients than the 38-protein alternative in 10 of 12 broad common-cancer cohorts.",
                "The 13 entries are an existing threshold-selected list, not an optimized minimal panel. Some qualify through a specific cancer and contribute little to common-cancer coverage.",
            ],
            [
                "Patients expressing at least one entry: lung adenocarcinoma 71.9%; lung squamous 59.8%; head/neck 36.1%; breast 20.1%.",
                "Coverage is low in pancreas (0.6%), prostate (3.4%) and colon (6.2%). These are cohort results, not worldwide patient estimates.",
            ],
        ],
        [360, 360],
    )
    text(
        "<b>Use the right identifiers:</b> 13 protein groups correspond to <b>16 gene symbols</b>. Page 2 gives the mapping. The import files use real member-gene symbols; grouped labels such as NY-ESO-1 are display names.",
        gap=12,
    )
    text(
        "<b>Read only as far as you need:</b> pages 1-2 choose and export; pages 3-8 check cancer coverage, normal tissues and alternatives; pages 9-10 explain the rules and files. The full report then opens into a bookmarked evidence library.",
        style=small,
    )
    end()

    start(
        "The list you can take into other work",
        "Sorted by breadth at p90. Keep the group mapping when your downstream analysis operates on genes rather than protein identities.",
        "1 / CHOOSE AND TAKE THE LIST",
    )
    rows = [
        [
            "Protein group",
            "Gene symbols to use",
            "Length (aa)",
            "Cancer-type groups >10%*",
            "World incidence represented**",
            "World mortality represented**",
        ]
    ]
    for r in ranked.itertuples():
        rows.append(
            [
                escape(r.Symbol),
                escape(r.member_symbols.replace(";", ", ")),
                str(r.protein_length_aa),
                f"{r.p90_n_cancer_type_groups_gt10}/69 ({r.p90_fraction_cancer_type_groups_gt10:.1%})",
                f"{r.world_incidence_pct_represented:.1f}%",
                f"{r.world_mortality_pct_represented:.1f}%",
            ]
        )
    table(rows, [110, 180, 65, 125, 120, 120], font=9)
    text(
        "*142 cohort views -> 100 with >=20 patients -> 69 groups after combining overlapping and molecular subtype views. A group counts when any eligible cohort has >10% p90-positive patients.",
        style=small,
        gap=8,
    )
    text(
        "**Sum of distinct burden categories hit by each protein, using oncoref world 2022 shares. Each category counts once. A subtype hit represents the whole category: these are <b>burden represented, not percentages of patients expressing the protein</b>. Details and gene imports: page 10.",
        style=small,
    )
    end()

    figure(
        "Where the list covers common cancers",
        "<b>Read the top row first:</b> the fraction of patients expressing any entry. Then inspect individual rows for the proteins relevant to your cancer.",
        DEST / "primary_common_cancer_heatmap.png",
        "Broad cohorts are used consistently. White = exactly zero; low nonzero values remain colored. Selection in one cancer does not imply high prevalence elsewhere.",
    )

    start(
        "The entries serve different purposes",
        "Use the full 13 as the starting candidate pool; prioritize within it according to the cancer types relevant to your next analysis.",
        "2 / CHECK THE EVIDENCE",
    )
    table(
        [
            ["What the existing data suggest", "How to use that information"],
            [
                "<b>Broadest observed p90 breadth:</b> PRAME (29/69 types), MAGEA3 (19/69), MAGEA6 (17/69), XAGE1A/B (13/69).",
                "Begin inspection here when breadth matters. XAGE1A/B contributes strongly in lung adenocarcinoma. Breadth is not the same as unique contribution to a combined panel.",
            ],
            [
                "<b>More cancer-specific entries:</b> DPPA5, TRIML2 and UTF1 each reach >10% in 1/69 types; ZFP42 and GFY in 2/69.",
                "Retain them when those cancers matter. The largest selection-passing cohort for DPPA5, TRIML2, UTF1 and ZFP42 is testicular germ cell tumors (n=148); GFY qualifies in rhabdoid tumors (n=63).",
            ],
            [
                "<b>Other useful context:</b> NY-ESO-1 reaches >10% in 8/69 types; MAGEA2/B and EGFL6 in 7/69; MAGEC1 in 3/69.",
                "Check the full cancer heatmaps before discarding a narrower entry. Largest qualifying cohorts include synovial sarcoma for NY-ESO-1 (n=47), meningioma for EGFL6 (n=379), and multiple myeloma for MAGEC1 (n=81).",
            ],
        ],
        [360, 360],
        font=10,
    )
    text(
        "<b>Do not interpret this ordering as a new filter.</b> Closely related proteins may cover many of the same patients. This report has not optimized unique patient coverage per added gene or applied a normal-tissue exclusion threshold."
    )
    text(
        "<b>A good working list depends on the next task.</b> For a general expression screen, carry all 13 groups / 16 genes. For a cancer-specific screen, start with the relevant column on page 3 and use the detailed cohort heatmaps in Appendix A.",
        style=small,
    )
    end()

    figure(
        "Normal tissues: inspect before narrowing",
        "<b>EGFL6 stands out for higher normal somatic RNA.</b> All 13 entries are labeled; use this plot to identify follow-up questions, not as a pass/fail safety screen.",
        DEST / "primary_normal_tissue_scatter.png",
        "This scatter uses the maximum single-member HPA annotation and excludes reproductive tissues and thymus. The next tissue maps include thymus among other tissues.",
    )
    figure(
        "Which normal tissues express these genes?",
        "<b>Look across the other-tissue block as well as the reproductive block.</b> The map shows HPA RNA, summed across genes that encode the same protein identity.",
        OUT / "plots/hpa_rna_tissue_overview.png",
        "HPA v23 consensus nTPM. Reproductive grouping includes breast; thymus is in other tissues. Summed RNA is not a protein measurement or evidence of expression in the same cells.",
    )
    figure(
        "Protein evidence is available for only some",
        "<b>Gray means no scored data, not no protein.</b> Five of the 13 entries have scored normal-tissue staining; use the individual HPA pages for the underlying evidence.",
        OUT / "plots/hpa_ihc_tissue_overview.png",
        "HPA immunohistochemistry: strongest score across member genes and cell types. White = measured Not detected. RNA and staining are distinct evidence types.",
    )
    figure(
        "When should you choose a different list?",
        "<b>Use 13 for a compact high-expression starting point.</b> Consider 38 if broader moderate expression is useful; the 3-entry list favors stringent prevalence in at least one cancer.",
        DEST / "list_size_coverage_comparison.png",
        "Compare rows within the same panel. The 13- and 38-entry lists are not nested. A stricter selection threshold does not guarantee broader coverage or lower normal-tissue expression.",
    )

    start(
        "The rules, in plain language",
        "These definitions explain what a selected entry means and let another analyst reproduce the list without guessing.",
        "3 / UNDERSTAND AND REUSE",
    )
    table(
        [
            ["Step", "What was done"],
            [
                "1. Define a protein entry",
                "Genes with identical longest Ensembl release 112 protein sequences are grouped. This is sequence-identity grouping, not an isoform-resolved proteomics result. MAGEA3 and MAGEA6 remain distinct.",
            ],
            [
                "2. Define a positive patient",
                "Sum expression for genes in the same protein group. Collapse the transcriptome in the same way, then compare to that patient's p70 or p90, including zeros. Expression must be strictly above the cutoff.",
            ],
            [
                "3. Select a candidate",
                "For the primary list, strictly more than 50% of patients must be positive in at least one cohort with at least 20 patient groups. All member genes must have measurements; missing is not zero.",
            ],
            [
                "4. Describe coverage",
                "Within a cohort, ANY counts each patient once if at least one listed entry is positive. Of 142 cohort views, 100 have >=20 patients. They form 73 patient-overlap groups and 69 groups after also merging molecular/evidence subtypes. These 69 are the breadth denominator; they are not the full registry.",
            ],
            [
                "5. Describe expression",
                "Mean TPM is calculated among positive patients only and only where linear TPM is available. Other supported expression scales can contribute to percentile-based selection.",
            ],
        ],
        [120, 600],
        font=10,
    )
    text(
        "<b>Burden columns:</b> qualifying cohort codes map to distinct oncoref burden categories; their global incidence or mortality shares are summed once per protein. No prevalence multiplier or renormalization is applied. The rounded reference totals are 99.9% incidence and 100.2% mortality. Subtype hits and broad categories limit interpretation; category mappings are retained for audit.",
        style=small,
    )
    text(
        "<b>What this establishes:</b> a reproducible RNA-expression candidate pool. Protein abundance, antigen presentation, targetability and suitability for a particular downstream assay are not established by these filters.",
        style=small,
    )
    end()

    start(
        "Use the files; open detail only as needed",
        "All files below are in the download bundle. Start with the gene handoff, then join the evidence you need using the stable protein-group key.",
        "3 / UNDERSTAND AND REUSE",
    )
    table(
        [
            ["Your next task", "File or section"],
            [
                "Load genes into another tool",
                "primary_panel/primary_gene_symbols.txt<br/>primary_panel/primary_ensembl_gene_ids.txt",
            ],
            [
                "Keep genes mapped to identical proteins",
                "primary_panel/primary_gene_handoff.csv<br/>Gene mapping, lengths and burden sums. primary_burden_coverage_details.csv lists every category included; coverage_denominator_audit.csv explains the 69-group denominator.",
            ],
            [
                "Inspect all annotations for the 13 groups",
                "primary_panel/primary_proteoforms.csv<br/>Includes qualifying cohorts, largest cohort size, positive-patient mean TPM where available, breadth and normal-tissue annotations.",
            ],
            [
                "Replot common-cancer coverage",
                "primary_panel/primary_common_cancer_heatmap.csv<br/>primary_panel/same_threshold_panel_comparison.csv",
            ],
            [
                "Choose another cutoff or inspect HPA",
                "selections/min20/ contains CSV and text files for all four lists.<br/>hpa_tissues/hpa_gene_links.csv links to original HPA pages.",
            ],
        ],
        [240, 480],
        font=10,
    )
    text(
        "<b>Evidence library in the full report:</b> A (p11) - primary-list cancer details; B (p18) - normal tissues and individual protein pages; C (p34) - alternatives; D (p68) - mortality context; E (p77) - patient distributions; F (p100) - methods and identity audit. Expand the PDF bookmarks to jump directly to a figure.",
        gap=12,
    )
    text(
        "<b>Keep the distinction:</b> primary_proteoforms.txt has 13 <i>display labels</i>; primary_gene_symbols.txt has the 16 <i>gene symbols</i> for gene-based tools. NY-ESO-1 = CTAG1A + CTAG1B; XAGE1A/B = XAGE1A + XAGE1B; MAGEA2/B = MAGEA2 + MAGEA2B.",
        style=small,
    )
    end()
    c.save()
    stream.seek(0)
    front = PdfReader(stream)
    assert len(front.pages) == 10
    guide = PdfWriter()
    for i, page in enumerate(front.pages):
        guide.add_page(page)
        guide.add_outline_item(titles[i], i)
    guide.add_metadata({"/Title": "CTA working list: 10-page decision guide"})
    with (OUT / "cta-working-list-guide.pdf").open("wb") as f:
        guide.write(f)

    original = PdfReader(OUT / "cta-analysis-all-figures.pdf")
    index = pd.read_csv(OUT / "pdf_page_index.csv").fillna("")
    writer = PdfWriter()
    root = writer.add_outline_item("START HERE: 10-page decision guide", 0)
    for i, page in enumerate(front.pages):
        writer.add_page(page)
        writer.add_outline_item(titles[i], i, parent=root)
    groups = [
        ("A. Primary list: cancer coverage and selection", [11, 46, 68, 69, 70, 19, 20]),
        ("B. Normal tissues: inspect each protein", [16, *range(29, 44)]),
        (
            "C. Alternative lists and selection thresholds",
            [*range(5, 11), 12, 13, 17, 44, 45, 47, *range(48, 68), 71, 72],
        ),
        ("D. Mortality context: category and cohort detail", [14, 22, 23, 24, 25, 26, 27, 28, 15]),
        ("E. Patient-level distributions (selected case studies)", list(range(73, 96))),
        ("F. Methods, identity mapping and cohort audit", [2, 3, 4, 18, 21, 96, 97, 98, 1]),
    ]
    included = [n for _, nums in groups for n in nums]
    assert sorted(included) == list(range(1, 99))
    mapping = []
    for group, nums in groups:
        parent = writer.add_outline_item(group, len(writer.pages))
        for old in nums:
            page = original.pages[old - 1]
            writer.add_page(page)
            new = len(writer.pages)
            title = index.iloc[old - 1].title
            writer.add_outline_item(title, new - 1, parent=parent)
            mapping.append(
                {"original_page": old, "new_page": new, "section": group, "title": title}
            )
    for i, page in enumerate(writer.pages[10:], 11):
        w, h = float(page.mediabox.width), float(page.mediabox.height)
        overlay = io.BytesIO()
        cc = canvas.Canvas(overlay, pagesize=(w, h))
        cc.setFillColor(colors.white)
        cc.rect(w - 84, 0, 84, 32, fill=1, stroke=0)
        cc.setFillColor(INK)
        cc.setFont("Helvetica", 8)
        cc.drawRightString(w - 36, 18, str(i))
        cc.save()
        overlay.seek(0)
        page.merge_page(PdfReader(overlay).pages[0])
    writer.add_metadata(
        {"/Title": "CTA working list: decision guide and complete evidence library"}
    )
    with (OUT / "cta-analysis-primary-panel.pdf").open("wb") as f:
        writer.write(f)
    pd.DataFrame(mapping).to_csv(DEST / "supporting_page_map.csv", index=False)
    (DEST / "decision_guide_review.md").write_text("""# Editorial reassessment

The previous opening repeated the selection rationale before providing a usable gene list, mixed figure galleries with decisions, and left readers to distinguish display labels from importable genes. The revised guide moves the list and identity mapping to page 2, separates breadth from cancer-specific usefulness, explains normal-tissue evidence before comparison options, and defers the complete 98-page evidence collection to six bookmarked sections. All 98 original pages are retained once. Original selections and numerical estimates are unchanged.

The 13-entry recommendation is a starting candidate pool, not an optimized minimal panel. Gene exports contain 16 member genes, preserving the 13 identity groups. No new safety filter or cancer-specific selection was applied.
""")
    readme = DEST / "README.md"
    readme.write_text(
        readme.read_text().split("\n## Start here")[0]
        + "\n## Start here\n\nUse ../cta-working-list-guide.pdf for the 10-page decision guide; ../cta-analysis-primary-panel.pdf adds the full bookmarked evidence library. primary_gene_symbols.txt and primary_ensembl_gene_ids.txt contain the 16 member genes. primary_gene_handoff.csv maps those genes to 13 protein groups. The original primary_proteoforms.txt remains a display-label list.\n"
    )
    report = OUT / "report.md"
    content = report.read_text()
    if "[10-page decision guide]" not in content:
        content = content.replace(
            "# CTA proteoform shortlists\n",
            "# CTA proteoform shortlists\n\n[10-page decision guide](cta-working-list-guide.pdf) | [Gene handoff CSV](primary_panel/primary_gene_handoff.csv) | [Importable gene symbols](primary_panel/primary_gene_symbols.txt)\n",
            1,
        )
    report.write_text(content)
    shutil.copy2(Path(__file__), OUT / "reproducibility" / Path(__file__).name)
    print(
        f"Built guide: {len(guide.pages)} pages; full report: {len(writer.pages)} pages; {len(genes)} genes / {len(primary)} protein groups"
    )


if __name__ == "__main__":
    main()
