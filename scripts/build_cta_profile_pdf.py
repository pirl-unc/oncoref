#!/usr/bin/env python3
"""Render the auditable ZNF730 calculation tables with live source hyperlinks."""

import xml.etree.ElementTree as ET
from pathlib import Path

import markdown
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/audits/cta-unassigned-20261029/znf730.md"
OUTPUT = ROOT / "output/pdf/znf730-diffuseness-tables.pdf"
WEB = "https://pirl-unc.github.io/oncoref/audits/cta-unassigned-20261029/"


def build(source=SOURCE, output=OUTPUT, *, title="ZNF730 diffuseness: every tissue calculation"):
    """Keep every native context, full reference table and interpretable weight."""
    soup = ET.fromstring(
        "<root>"
        + markdown.markdown(
            source.read_text().replace("—", "-").replace("–", "-").replace("→", "to"),
            extensions=["tables"],
        )
        + "</root>"
    )
    for anchor in soup.findall(".//a"):
        if not anchor.attrib["href"].startswith(("https://", "http://")):
            anchor.attrib["href"] = WEB + anchor.attrib["href"]
        anchor.attrib["color"] = "#087E8B"
    for tag in soup.findall(".//code"):
        tag.tag = "font"
        tag.attrib["name"] = "Courier"
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            "BodyCTA",
            fontName="Helvetica",
            fontSize=9,
            leading=13,
            textColor=colors.HexColor("#173B55"),
            spaceAfter=7,
        )
    )
    styles.add(
        ParagraphStyle("CellCTA", parent=styles["BodyCTA"], fontSize=8.5, leading=10, spaceAfter=0)
    )
    styles.add(
        ParagraphStyle(
            "HeadCTA", parent=styles["CellCTA"], textColor=colors.white, fontName="Helvetica-Bold"
        )
    )
    for name in ["Heading1", "Heading2", "Heading3"]:
        styles[name].textColor = colors.HexColor("#173B55")
        styles[name].alignment = TA_LEFT
        styles[name].keepWithNext = 1
    story = []
    width = A4[0] - 32 * mm

    def content(node):
        return (node.text or "") + "".join(ET.tostring(child, encoding="unicode") for child in node)

    for node in soup:
        if node.tag in {"h1", "h2", "h3"}:
            story.append(Paragraph(content(node), styles["Heading" + node.tag[1]]))
        elif node.tag == "p":
            story.append(Paragraph(content(node), styles["BodyCTA"]))
        elif node.tag == "table":
            rows = [
                [
                    Paragraph(content(cell), styles["HeadCTA" if cell.tag == "th" else "CellCTA"])
                    for cell in list(row)
                ]
                for row in node.findall(".//tr")
            ]
            n = len(rows[0])
            weights = (
                [2.9, 1, 1, 1.35, 1.35]
                if n == 5
                else [2, 1]
                if n == 2
                else [1.1, 1, 1, 1, 0.8, 1.3, 1.3]
                if n == 7
                else [1] * n
            )
            table = Table(rows, colWidths=[width * w / sum(weights) for w in weights], repeatRows=1)
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#173B55")),
                        (
                            "ROWBACKGROUNDS",
                            (0, 1),
                            (-1, -1),
                            [colors.white, colors.HexColor("#F3F6F8")],
                        ),
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ("TOPPADDING", (0, 0), (-1, -1), 4),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                        ("LINEBELOW", (0, 0), (-1, -1), 0.3, colors.HexColor("#D9E2E8")),
                    ]
                )
            )
            story.extend([table, Spacer(1, 9)])
    output.parent.mkdir(parents=True, exist_ok=True)

    def footer(canvas, _doc):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#536675"))
        canvas.drawString(16 * mm, 10 * mm, "OncoRef | 29 October 2026 | RNA profile identity")
        canvas.drawRightString(A4[0] - 16 * mm, 10 * mm, str(canvas.getPageNumber()))

    SimpleDocTemplate(
        str(output),
        pagesize=A4,
        leftMargin=16 * mm,
        rightMargin=16 * mm,
        topMargin=16 * mm,
        bottomMargin=18 * mm,
        title=title,
        author="OncoRef",
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    print(output)


if __name__ == "__main__":
    build()
    build(
        SOURCE.with_name("diffuseness-impact.md"),
        OUTPUT.with_name("cta-diffuseness-impact.pdf"),
        title="CTA diffuseness: chosen method and measured impact",
    )
