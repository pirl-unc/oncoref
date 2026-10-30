#!/usr/bin/env python3
"""Render the three-question CTA definition with independent safety and admission.

This explanatory figure reads the documented rule contract; the selection replay supplies decisions.
The source script, PDF and PNG preserve the distinction between safety evidence
and the independent candidate admission requirements.
"""

import shutil
import subprocess
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
INK = colors.HexColor("#173B55")
TEAL = colors.HexColor("#087E8B")


def build():
    """Write one readable landscape page and the matching documentation image."""
    destination = ROOT / "output/pdf/cta-evidence-axes.pdf"
    destination.parent.mkdir(parents=True, exist_ok=True)
    images = ROOT / "docs/images"
    images.mkdir(exist_ok=True)
    page = canvas.Canvas(str(destination), pagesize=(792, 612), invariant=True)
    page.setTitle("CTA evidence: safety and identity")
    page.setAuthor("OncoRef")
    page.setFillColor(TEAL)
    page.setFont("Helvetica-Bold", 12)
    page.drawString(36, 575, "ONCOREF  /  SAFETY AND CANDIDATE IDENTITY")
    page.setFillColor(INK)
    page.setFont("Helvetica-Bold", 26)
    page.drawString(36, 537, "Safety, cancer and Reproductive evidence")
    paragraph = ParagraphStyle(
        "body", fontName="Helvetica", fontSize=12, leading=16, textColor=INK, alignment=TA_LEFT
    )
    heading = ParagraphStyle("heading", parent=paragraph, textColor=colors.white)

    def p(text, header=False):
        return Paragraph(text, heading if header else paragraph)

    entries = [
        ["Assessment", "Required evidence", "Current implementation"],
        [
            "<b>Safety</b>",
            "Acceptable normal-tissue RNA / protein evidence across Vital, Protected and Manageable compartments.",
            "Active. Only this assessment can supply Unsafe evidence. Manageable excess blocks CTA under the current rules.",
        ],
        [
            "<b>Cancer expression</b>",
            "Expression in at least one measured cancer sample.",
            "Active. At least 2 in native TPM or HPA pTPM. One sample qualifies. No tumor-versus-normal ratio.",
        ],
        [
            "<b>Reproductive evidence</b>",
            "Affirmative expression in reproductive, placental, embryonic or fetal contexts.",
            "Active quantified routes. Reproductive or supported developmental expression can qualify.",
        ],
    ]
    table = Table(
        [[p(value, row == 0) for value in values] for row, values in enumerate(entries)],
        colWidths=[150, 275, 295],
        rowHeights=[31, 84, 65, 77],
    )
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), INK),
                ("BACKGROUND", (0, 3), (-1, 3), colors.HexColor("#F3F6F8")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 9),
                ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                ("TOPPADDING", (0, 0), (-1, -1), 7),
                ("LINEBELOW", (0, 0), (-1, -1), 0.5, colors.HexColor("#D9E2E8")),
            ]
        )
    )
    table.wrapOn(page, 720, 400)
    table.drawOn(page, 36, 235)
    statement = p(
        "<b>CTA:</b> safety AND cancer expression AND Reproductive evidence, with complete "
        'required coverage, including the 52-context <link href="https://gtexportal.org/home/downloads/adult-gtex">GTEx v11</link> core.<br/>'
        "<b>Reproductive</b> includes reproductive organs, placenta, extraembryonic, embryonic and fetal contexts.<br/>"
        "<b>Diffuseness:</b> HPA average N / (N + R) >10% blocks a complete adult/placental route. "
        "Disposable counts; Permissive thymus is excluded. Embryonic concentration is unassessed.<br/>"
        "<b>Only safety establishes Unsafe.</b> Missing core data remain unassigned. "
        "The former maximum-versus-maximum restriction rule is removed."
    )
    _, height = statement.wrap(720, 120)
    statement.drawOn(page, 36, 205 - height)
    page.setFont("Helvetica", 9)
    page.drawString(
        36,
        25,
        "Local policy: oncoref.cta-selection.2026-10-29.4. Clinical RNA credit retains tested epitope limits.",
    )
    page.save()
    shutil.copyfile(destination, images / destination.name)
    subprocess.run(
        [
            "pdftoppm",
            "-singlefile",
            "-scale-to",
            "1800",
            "-png",
            str(destination),
            str(images / "cta-evidence-axes"),
        ],
        check=True,
    )
    print(destination)


if __name__ == "__main__":
    build()
