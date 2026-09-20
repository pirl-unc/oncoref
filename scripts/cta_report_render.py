"""Data-driven, deterministic rendering shared by the CTA report entry points."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

import pandas as pd
from cta_report_common import seal_stage, verify_analysis
from PIL import Image
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle


def figure_rows(out):
    path = Path(out) / "plot_index.csv"
    if not path.exists():
        return pd.DataFrame(columns=["name", "title", "category"])
    rows = pd.read_csv(path).fillna({"category": "other"})
    if rows.name.duplicated().any():
        raise ValueError("Duplicate figure names")
    for name in rows.name:
        if not (Path(out) / "plots" / f"{name}.png").is_file():
            raise ValueError(f"Missing indexed figure: {name}")
    return rows.sort_values(["category", "name"], kind="stable")


class Report:
    def __init__(self, path):
        self.size = landscape(A3)
        self.width, self.height = self.size
        self.pdf = canvas.Canvas(str(path), pagesize=self.size, invariant=1)
        self.pdf.setTitle("CTA expression analysis")
        self.index = []
        self.style = ParagraphStyle("body", fontName="Helvetica", fontSize=13, leading=18)
        self.cell = ParagraphStyle("cell", fontName="Helvetica", fontSize=10, leading=13)

    def start(self, title):
        self.pdf.setFillColor(colors.HexColor("#127C80"))
        self.pdf.rect(0, self.height - 42, self.width, 42, fill=1, stroke=0)
        self.pdf.setFillColor(colors.white)
        self.pdf.setFont("Helvetica-Bold", 19)
        self.pdf.drawString(36, self.height - 28, title)
        self.pdf.setFillColor(colors.HexColor("#203746"))
        self.y = self.height - 70

    def paragraph(self, text):
        paragraph = Paragraph(escape(text), self.style)
        _, height = paragraph.wrap(self.width - 72, 2000)
        if self.y - height < 44:
            raise ValueError("Report paragraph exceeds page; split the content")
        paragraph.drawOn(self.pdf, 36, self.y - height)
        self.y -= height + 18

    def finish(self, title, category, figure=""):
        number = len(self.index) + 1
        self.pdf.setFont("Helvetica", 9)
        self.pdf.drawString(
            36, 20, "RNA expression candidates; atlas measurements do not establish clinical safety"
        )
        self.pdf.drawRightString(self.width - 36, 20, str(number))
        self.pdf.bookmarkPage(str(number))
        self.pdf.addOutlineEntry(title, str(number), 0)
        self.pdf.showPage()
        self.index.append(
            {"pdf_page": number, "title": title, "category": category, "figure_name": figure}
        )

    def table(self, frame, title, category="selections"):
        # Fixed small chunks keep arbitrary result counts readable and paginated.
        for start in range(0, max(1, len(frame)), 12):
            part = frame.iloc[start : start + 12].fillna("unavailable")
            page_title = (
                f"{title}: rows {start + 1}-{start + len(part)}"
                if len(part)
                else f"{title}: no selections"
            )
            self.start(page_title)
            data = [[Paragraph(escape(str(v)), self.cell) for v in frame.columns]]
            data += [
                [Paragraph(escape(str(v)), self.cell) for v in row]
                for row in part.itertuples(index=False, name=None)
            ]
            table = Table(
                data,
                colWidths=[(self.width - 72) / len(frame.columns)] * len(frame.columns),
                repeatRows=1,
            )
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EDED")),
                        (
                            "ROWBACKGROUNDS",
                            (0, 1),
                            (-1, -1),
                            [colors.white, colors.HexColor("#F0F6F6")],
                        ),
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ("TOPPADDING", (0, 0), (-1, -1), 8),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
                    ]
                )
            )
            _, height = table.wrap(self.width - 72, self.height)
            if height > self.y - 44:
                raise ValueError(f"Table exceeds page: {title}")
            table.drawOn(self.pdf, 36, self.y - height)
            self.finish(page_title, category)

    def figure(self, path, title, category, name):
        self.start(title)
        with Image.open(path) as im:
            width, height = im.size
        scale = min((self.width - 40) / width, (self.y - 45) / height)
        self.pdf.drawImage(
            str(path),
            (self.width - width * scale) / 2,
            44 + (self.y - 45 - height * scale) / 2,
            width=width * scale,
            height=height * scale,
        )
        self.finish(title, category, name)

    def save(self):
        self.pdf.save()
        return pd.DataFrame(self.index)


def render_report(out):
    out = Path(out)
    verify_analysis(out)
    manifest = json.loads((out / "run_manifest.json").read_text())
    plots = figure_rows(out)
    report = Report(out / "cta-analysis-all-figures.pdf")
    report.start("CTA expression analysis")
    report.paragraph(
        "Candidates are selected by strict within-patient expression and cohort prevalence thresholds. All cohorts use the same declared biological gene background; proteoform analysis collapses that background consistently before computing quantiles. Additional source rows cannot change the cutoff."
    )
    report.paragraph(
        f"This run contains {len(plots)} indexed figures. Every figure category is included. Tables below are generated from the current selections; gene identifiers and identical-protein identities remain separate."
    )
    report.paragraph(
        "Missing measurements cannot qualify. Cohort views can share patients or represent molecular subtypes; breadth uses the recorded overlap and cancer-type groups. Small cohorts and p30/p50 gene screens are exploratory. Low-percentile cutoffs can equal zero and should not be interpreted as strong expression."
    )
    report.paragraph(
        "Normal RNA sums the same member loci as tumor RNA. HPA zeros are rounded estimates below 0.05 nTPM per locus. RNA and IHC do not measure peptide presentation. Burden-category shares describe represented categories, not worldwide expression prevalence; mortality cohort estimates use broad histologies."
    )
    report.finish("Analysis definitions and limitations", "methods")
    selected_fields = [
        k
        for k in (
            "common_background_genes",
            "minimum_patients_for_selection",
            "cohort_policies",
            "min_cohort_size",
            "large_min_cohort_size",
            "cohorts_analyzed",
            "percentiles",
            "qc_policy",
            "repeat_sample_policy",
        )
        if k in manifest
    ]
    report.table(
        pd.DataFrame(
            {"Run setting": selected_fields, "Value": [str(manifest[k]) for k in selected_fields]}
        ),
        "Run settings",
        "methods",
    )
    paths = sorted((out / "selections").glob("*/*/*ranked.csv"))
    if not paths:
        paths = sorted((out / "gene_sets").glob("*/genes.csv"))
    for path in paths:
        frame = pd.read_csv(path)
        columns = [
            c
            for c in [
                "Symbol",
                "member_symbols",
                "protein_length_aa",
                "n_passing_cancer_type_groups",
                "n_overlap_adjusted_groups",
                "n_qualifying_cohorts",
                "largest_passing_cohort_size",
                "best_cohort_n_patients",
            ]
            if c in frame
        ]
        report.table(frame[columns], " / ".join(path.relative_to(out).parts[1:-1]))
    for row in plots.itertuples():
        report.figure(out / "plots" / f"{row.name}.png", row.title, row.category, row.name)
    index = report.save()
    index.to_csv(out / "pdf_page_index.csv", index=False)
    if set(index.loc[index.figure_name.ne(""), "figure_name"]) != set(plots.name):
        raise ValueError("PDF does not contain exactly the indexed figures")
    lines = [
        "# CTA expression analysis",
        "",
        "[Complete tables and figures](cta-analysis-all-figures.pdf)",
        "",
        f"{len(plots)} figures across all categories; {len(index)} PDF pages.",
        "",
        "See run_manifest.json, validation.json and analysis_receipt.json for definitions and verified hashes.",
    ]
    lines += [f"- [{r.title}](plots/{r.name}.png)" for r in plots.itertuples()]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    seal_stage(
        out,
        "render",
        [
            Path(__file__),
            out / "analysis_receipt.json",
            *[out / "plots" / f"{name}.png" for name in plots.name],
        ],
        [out / "cta-analysis-all-figures.pdf", out / "pdf_page_index.csv"],
    )
    return index
