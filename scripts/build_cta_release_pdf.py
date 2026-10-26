#!/usr/bin/env python3
"""Render the 1.20.0 CTA decision brief from its pinned policy and candidate table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd
from pypdf import PdfReader
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    Flowable,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
BLUE, TEAL, RED = "#152F47", "#087E8B", "#A5413A"
BASE = "https://pirl-unc.github.io/oncoref/"


class Funnel(Flowable):
    """Two proportional stacked bars; same denominator, mutually exclusive outcomes."""

    def __init__(self, old, new):
        super().__init__()
        self.old, self.new, self.width, self.height = old, new, 480, 115

    def draw(self):
        c = self.canv
        palette = [TEAL, "#C99236", RED]
        for y, label, counts in [
            (78, "Previous release", self.old),
            (22, "Reviewed policy", self.new),
        ]:
            c.setFillColor(colors.HexColor(BLUE))
            c.setFont("Helvetica-Bold", 10)
            c.drawString(0, y + 20, label)
            x = 0
            for value, color in zip(counts, palette):
                w = 480 * value / sum(counts)
                c.setFillColor(colors.HexColor(color))
                c.rect(x, y, w, 15, fill=1, stroke=0)
                x += w
            c.setFillColor(colors.HexColor(BLUE))
            c.setFont("Helvetica", 9)
            c.drawString(
                0,
                y - 12,
                f"CTA {counts[0]:,}    Neither {counts[1]:,}    Unsafe {counts[2]:,}    Total {sum(counts):,}",
            )


def build(output):
    report = ROOT / "docs/audits/cta-policy-20261025"
    m = json.loads((report / "policy.json").read_text())
    assert m["policy"] == "oncoref.target-sets.2026-10-25.2"
    counts = m["counts"]
    rows = pd.read_csv(report / "candidate-results.csv.gz", low_memory=False).set_index("Gene")
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            "TitleCTA",
            fontName="Helvetica-Bold",
            fontSize=25,
            leading=29,
            textColor=colors.HexColor(BLUE),
            spaceAfter=16,
        )
    )
    styles.add(
        ParagraphStyle(
            "HCTA",
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=17,
            textColor=colors.HexColor(TEAL),
            spaceBefore=15,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(
            "BodyCTA",
            fontName="Helvetica",
            fontSize=10.5,
            leading=15,
            textColor=colors.HexColor(BLUE),
            spaceAfter=10,
        )
    )
    styles.add(
        ParagraphStyle(
            "CellCTA", fontName="Helvetica", fontSize=9, leading=12, textColor=colors.HexColor(BLUE)
        )
    )
    styles.add(
        ParagraphStyle(
            "SmallCTA",
            fontName="Helvetica",
            fontSize=9,
            leading=12,
            textColor=colors.HexColor(BLUE),
            spaceAfter=8,
        )
    )
    story = []

    def p(text, style="BodyCTA"):
        return Paragraph(text, styles[style])

    def add(text, style="BodyCTA"):
        story.append(p(text, style))

    def table(data, widths):
        t = Table(
            [[p(str(c), "CellCTA") for c in row] for row in data],
            colWidths=widths,
            repeatRows=1,
            hAlign="LEFT",
        )
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#DCEBF0")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 8),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                    ("TOPPADDING", (0, 0), (-1, -1), 8),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
                    ("LINEBELOW", (0, 0), (-1, -1), 0.4, colors.HexColor("#D6DEE4")),
                ]
            )
        )
        story.append(t)

    def link(label, url):
        return f'<link href="{escape(url)}" color="{TEAL}"><u>{label}</u></link>'

    add("CTA decisions<br/>Evidence first, explicit boundaries", "TitleCTA")
    add("OncoRef 1.20.0 | 25 October 2026", "HCTA")
    add("213 CTA genes / 359 distinct complete target proteins", "HCTA")
    add(
        "The research target set now separates credible normal expression from unresolved evidence. A completed review can end in a hold: it does not require pretending an ambiguous stain has been experimentally resolved."
    )
    story.append(
        Funnel([399, 561, 1580], [counts["nominated"][g] for g in ["CTA", "Neither", "Unsafe"]])
    )
    story.append(Spacer(1, 18))
    table(
        [
            ["Outcome", "Meaning"],
            [
                "CTA",
                "Affirmative reproductive admission and all applicable normal gates pass; source-backed nomination required.",
            ],
            [
                "Neither",
                "Material uncertainty, inadequate evidence/coverage, or insufficient selectivity. Neither target nor default forbidden background.",
            ],
            [
                "Unsafe",
                "Accepted normal protein, excessive normal RNA, or an applicable safety veto. Included in the forbidden background.",
            ],
        ],
        [85, 395],
    )
    add("What changed", "HCTA")
    add(
        "From the released 399-gene list: <b>6 move to Unsafe and 180 to Neither</b>. From the unpublished 391-gene draft: all 178 antibody cases receive dispositions (177 Neither, 1 Unsafe). No member of that discrepancy queue remains retained."
    )
    add(
        "Same 2,540 nominated genes in both bars. Numerical RNA/IHC settings and anatomy did not change. CTA does not mean clinically proven safe; Unsafe does not mean a therapy has demonstrated injury.",
        "SmallCTA",
    )
    add(
        link(
            "Full report, every candidate and downloadable tables",
            BASE + "audits/cta-policy-20261025/",
        ),
        "SmallCTA",
    )

    story.append(PageBreak())
    add("One decision hierarchy", "TitleCTA")
    for title, text in [
        (
            "1. Identify",
            "Keep species, gene/reference, normal cell state, assay, native units and source version. Missing evidence is unknown, never zero.",
        ),
        (
            "2. Exclude",
            "Accepted non-permitted normal protein or curated positive IHC (including Low), a failing absolute RNA ceiling, or an applicable safety veto gives Unsafe. Any audited non-permitted/unsited value above 30 TPM or nTPM gives Unsafe.",
        ),
        (
            "3. Hold",
            "Unresolved normal-expression or specificity conflicts, uncertain antibody/cell-RNA evidence, or qualified donor variability give Neither. A reviewed Variable finding can turn its corresponding ordinary RNA failure into Neither, but never overrides &gt;30, accepted normal protein, positive IHC or a safety veto.",
        ),
        (
            "4. Admit and qualify",
            "Require reproductive RNA at least 5 nTPM, reliable reproductive protein, or accepted reproductive plus tumor/recognition evidence. Then require nomination, measured coverage and restriction. Protein support bypasses the admission floor only.",
        ),
    ]:
        add(title, "HCTA")
        add(text)
    add("Bulk RNA ceilings (nTPM) / selectivity settings", "HCTA")
    table(
        [
            ["Tissue", "Missing/uncertain", "Approved -", "Supported -", "Enhanced -"],
            ["Vital", "2 / 80%", "4 / 75%", "6 / 70%", "12 / 65%"],
            ["Protected", "6 / 75%", "8 / 70%", "12 / 65%", "24 / 60%"],
            ["Manageable", "10 / 70%", "20 / 65%", "30 / 60%", "60 / 55%*"],
        ],
        [82, 110, 96, 96, 96],
    )
    add(
        "R = maximum reproductive nTPM; N = normal tissue nTPM; C = ceiling; s = setting.<br/><b>Require N &lt;= C, N &lt;= R x (1 - s) / s + 0.1, and N &lt; R.</b><br/>*The independent &gt;30 rule caps the effective Manageable ceiling at 30.",
        "SmallCTA",
    )
    add(
        "Relaxation requires complete routine negative local IHC at the weakest usable confidence and reliable reproductive-positive profile support. RNA/IHC donors are not presumed paired. Reproductive and Disposable are permitted; the other three classes are screened. The numerical policy is human-only.",
        "SmallCTA",
    )

    story.append(PageBreak())
    add("Important target changes", "TitleCTA")
    reasons = [
        ("FTMT", "Normal human neuronal protein; small postmortem series."),
        ("ABCC11", "Reported normal neural protein; full reagent methods remain an access limit."),
        ("YBX2", "Normal vascular smooth-muscle/adrenal protein; not cardiomyocyte evidence."),
        ("CCDC158", "Reported renal tubular protein/function; isoform-specific rescue unresolved."),
        ("ACTRT1", "Endogenous protein in differentiating normal keratinocytes."),
        ("MAGEC2", "Purkinje-cell staining: gene/family specificity unresolved."),
        ("LGALS13", "Extra-placental reactivity: galectin identity and anatomy unresolved."),
        (
            "PLCZ1",
            "Corrected to qualitative pancreatic RNA. The original staining was CFTR, not PLCZ1.",
        ),
        (
            "PRAME",
            "Independent normal adrenal-protein finding accepted; another negative panel does not erase it.",
        ),
        (
            "MAGEA4",
            "Raw neuronal/epidermal stains unresolved. Independent OTI1F9 validation is reassuring but does not resolve the older reagents.",
        ),
    ]
    table(
        [["Gene", "Outcome", "Reason and limitation"]]
        + [[g, rows.loc[g, "Category"], r] for g, r in reasons],
        [75, 70, 335],
    )
    add("Preserve clinical evidence at its actual scope", "HCTA")
    add(
        "PRAME/IMA203 and MAGE-A4/afami-cel retain product-specific clinical evidence. Peptide, HLA, receptor and treatment context matter. This gene-wide policy is not a judgement that those products are clinically unsafe."
    )
    add(
        "<b>Still CTA:</b> CTAG1A/B, XAGE1A/B and CT55. <b>Still Neither:</b> CTAG2/Variable and FOXR2/Antibody. Familiarity with a target is not an exemption from normal-evidence rules.",
        "SmallCTA",
    )
    add(
        link(
            "Citations, full access-depth notes and all 178 individual adjudications",
            BASE + "audits/cta-adjudication-20261025/",
        ),
        "SmallCTA",
    )

    story.append(PageBreak())
    add("Coverage, caveats and adoption", "TitleCTA")
    add("What is complete", "HCTA")
    add(
        "All eight initial cases and the 178 remaining antibody cases have explicit decisions, source-linked findings, limitations and reconsideration requirements. The antibody work reviews pinned annotations, normal-cell contexts, reagents, profile notes and RNA context; it is not independent rereading of every histology image or experimental resolution of every stain."
    )
    add("What remains uncertain", "HCTA")
    add(
        "The initial literature search covered all 2,540 nominations; only 32 received new individual primary-paper assessment in that initial pass, with 46 specific notes including prior work. Additional targeted papers informed this adjudication. Comprehensive gene-by-gene literature review, continuous normal-expression breadth and calibrated RNA uncertainty remain unfinished. Bulk RNA can miss rare or activated normal cells."
    )
    add("Cell evidence stays in native units", "HCTA")
    add(
        "HPA single-cell: at least 5 nTPM and 20 reads, with reproductive comparison/coverage checks. DICE: at least 5 TPM in three distribution values; these are not stable cross-gene donor IDs. HPA/Monaco immune summaries: at least 5 in native units. Cornea: mean CP10K at least 1, five positive cells and 20 counts in two groups. These give review holds unless a stronger exclusion applies; CP10K is never converted to TPM.",
        "SmallCTA",
    )
    add("Use matching policies and proteins", "HCTA")
    table(
        [
            ["Scope", "CTA", "Neither", "Unsafe"],
            ["Nominated genes", "213", "741", "1,586"],
            ["Canonical genes", "213", "2,179", "17,705"],
            ["Complete proteins", "359", "3,230", "126,396"],
        ],
        [165, 105, 105, 105],
    )
    add(
        "A shared exact protein takes Unsafe, then Neither, then CTA precedence. Forbidden = unsafe.fasta.gz + additional-forbidden.fasta.gz. Neither is not forbidden. Full proteins can still share k-mers; subtract the forbidden union without crossing record boundaries or ambiguous residues.",
        "SmallCTA",
    )
    add(
        "<b>Explicit migration required:</b> target_sets.results() supplies this policy. cta_policy.manifest() supplies active version, scope, thresholds, counts and hashes. cta_policy.findings() and adjudications() expose the evidence trail. The historical 624-gene CTA API is unchanged; upgrading the package alone does not select this list.",
        "SmallCTA",
    )
    add(
        link(
            "Release and matching sequence assets",
            "https://github.com/pirl-unc/oncoref/releases/tag/v1.20.0",
        )
        + " | "
        + link("Rules and source references", BASE + "target-sets/")
        + " | "
        + link(
            "Reproducible PDF script",
            "https://github.com/pirl-unc/oncoref/blob/v1.20.0/scripts/build_cta_release_pdf.py",
        ),
        "SmallCTA",
    )

    output.parent.mkdir(parents=True, exist_ok=True)

    def footer(canvas, doc):
        canvas.setFillColor(colors.HexColor(BLUE))
        canvas.setFont("Helvetica", 8)
        canvas.drawString(58, 25, "OncoRef 1.20.0 | " + m["policy"])
        canvas.drawRightString(538, 25, str(doc.page))

    SimpleDocTemplate(
        str(output),
        pagesize=(596, 842),
        leftMargin=58,
        rightMargin=58,
        topMargin=46,
        bottomMargin=45,
        title="OncoRef 1.20.0 CTA decisions",
        author="OncoRef",
        subject=m["policy"],
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    pdf = PdfReader(output)
    assert len(pdf.pages) == 4, len(pdf.pages)
    for page in pdf.pages:
        assert page.extract_text().strip()
    print(f"{output}: 4 pages; policy {m['policy']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "output/pdf/cta-policy-1.20.0.pdf")
    build(parser.parse_args().output)
