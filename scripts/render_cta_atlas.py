#!/usr/bin/env python3
"""Create a navigable overview and complete CTA/Neither/Unsafe PDF catalogs."""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/audits/cta-atlas-20261024"
OUT = ROOT / "output/pdf"
TMP = ROOT / "tmp/pdfs"
BASE_URL = "https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/"
INK = colors.HexColor("#163047")
MUTED = colors.HexColor("#526673")
TEAL = colors.HexColor("#087F79")
PALE = colors.HexColor("#F0F5F7")
STATUS = {"CTA": "#087F79", "Neither": "#A76A12", "Unsafe": "#A33F45"}
W, H = A4
WIDTH = W - 86


def clean(value):
    text = str(value) if pd.notna(value) else ""
    return (
        text.replace("→", " -> ")
        .replace("≥", ">=")
        .replace("≤", "<=")
        .replace("—", " - ")
        .replace("–", "-")
        .replace("‑", "-")
    )


def esc(value):
    return html.escape(clean(value))


def short(value, limit=350):
    text = clean(value).strip()
    if len(text) <= limit:
        return text
    return text[:limit].rsplit(" ", 1)[0] + " ... [full text in table]"


def link(url, label):
    return f'<a href="{html.escape(url, quote=True)}" color="#087F79">{esc(label)}</a>'


def paragraph(text, style="Body", raw=False):
    content = text if raw else esc(text)
    if style == "Code":
        content = content.replace("\n", "<br/>")
    return Paragraph(content, STYLES[style])


def section(title, key=None):
    p = paragraph(title, "Section")
    p.bookmark = key or re.sub(r"\W+", "-", title.lower())
    return p


def table(headers, rows, widths=None, font="Table"):
    cells = [[paragraph(v, "TableHead") for v in headers]]
    cells += [[paragraph(v, font) for v in row] for row in rows]
    t = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.spaceAfter = 8
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), INK),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 7),
                ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                ("TOPPADDING", (0, 0), (-1, -1), 7),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
                ("LINEBELOW", (0, 0), (-1, 0), 0.7, INK),
                ("LINEBELOW", (0, 1), (-1, -1), 0.25, colors.HexColor("#D7E0E5")),
            ]
        )
    )
    return t


class Report(BaseDocTemplate):
    def __init__(self, path, title, catalog=False):
        super().__init__(
            str(path),
            pagesize=A4,
            leftMargin=43,
            rightMargin=43,
            topMargin=54,
            bottomMargin=47,
            title=title,
            author="oncoref",
            pageCompression=1,
        )
        self.report_title = title
        self.gene_pages = {}
        self.catalog = catalog
        self.addPageTemplates(
            PageTemplate(
                id="main",
                frames=[
                    Frame(
                        43,
                        47,
                        WIDTH,
                        H - 101,
                        leftPadding=0,
                        rightPadding=0,
                        topPadding=0,
                        bottomPadding=0,
                    )
                ],
                onPage=self.decorate,
            )
        )

    def decorate(self, canvas, doc):
        canvas.saveState()
        canvas.setFillColor(INK)
        canvas.setFont("AtlasBold", 9)
        canvas.drawString(43, H - 27, "ONCOREF / CTA EVIDENCE ATLAS")
        canvas.setStrokeColor(colors.HexColor("#CCD9DE"))
        canvas.line(43, H - 36, W - 43, H - 36)
        canvas.setFont("Atlas", 8)
        canvas.setFillColor(MUTED)
        canvas.drawString(43, 29, "24 OCT 2026 | Policy 2026-10-24.1 | Research classifications")
        canvas.drawRightString(W - 43, 29, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if hasattr(flowable, "bookmark"):
            self.canv.bookmarkPage(flowable.bookmark)
            self.canv.addOutlineEntry(flowable.getPlainText(), flowable.bookmark, 0, False)
        if hasattr(flowable, "gene_id"):
            self.gene_pages[flowable.gene_id] = self.page
            self.canv.bookmarkPage(flowable.gene_id)
            self.canv.addOutlineEntry(flowable.getPlainText(), flowable.gene_id, 0, False)


def cover(title, subtitle):
    return [Spacer(1, 18), paragraph(title, "Title"), paragraph(subtitle, "Deck"), Spacer(1, 10)]


def figure_heatmap(genes):
    rna = pd.read_csv(DATA / "cancer-rna.csv.gz")
    symbols = [
        "PRAME",
        "MAGEA4",
        "MAGEA3",
        "MAGEA10",
        "CTAG1A",
        "CTAG1B",
        "CTAG2",
        "XAGE1A",
        "XAGE1B",
        "CT83",
        "SSX2",
        "ROPN1",
        "CT55",
        "FOXR2",
        "HORMAD1",
        "TDRD1",
        "BRDT",
        "DAZL",
        "TEX101",
        "ACRBP",
    ]
    q = genes.loc[genes.Gene.isin(symbols), ["Gene", "Gene ID", "Category"]]
    d = rna.loc[rna.cohort.eq("TCGA") & rna.samples.ge(20)].merge(
        q, left_on="gene_id", right_on="Gene ID"
    )
    matrix = d.pivot(index="Gene", columns="cancer_code", values="prevalence_ptpm_ge_5").reindex(
        symbols
    )
    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#DDE2E5")
    image = ax.imshow(matrix.values * 100, aspect="auto", vmin=0, vmax=100, cmap=cmap)
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(matrix.index)), matrix.index, fontsize=10)
    category = q.set_index("Gene").Category
    for label in ax.get_yticklabels():
        label.set_color(STATUS[category[label.get_text()]])
    ax.set_title(
        "Cancer RNA occurrence: fraction of measured samples >=5 pTPM",
        loc="left",
        fontsize=12,
        pad=15,
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025, label="Samples (%)")
    fig.tight_layout()
    path = TMP / "cancer-rna-heatmap.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def reason_text(row):
    category, detail = row["Category"], row["Detail"]
    if category == "CTA":
        return "Passes the frozen research screen. Clinical efficacy and safety are not implied."
    detail_text = {
        "Protein": "Normal IHC is positive outside permitted compartments.",
        "RNA": "Normal bulk RNA exceeds a category/IHC-dependent ceiling.",
        "High RNA": "At least one non-permitted or unsited RNA observation is >30 TPM or nTPM.",
        "Normal": "Published protein production in non-permitted normal cells.",
        "Safety": "A separately reviewed safety concern overrides expression nomination.",
        "Low": "Insufficient affirmative reproductive evidence under the 5 nTPM/protein-support admission rule.",
        "Antibody": "Primary/curated antibody evidence needs resolution.",
        "Cellular": "Cell-resolved or immune RNA conflicts require review.",
        "Variable": "Vital-organ donor RNA varies from reported zero to a failing value; negative IHC is from other specimens.",
        "Selectivity": "Insufficient reproductive enrichment; selectivity alone is not proof of unsafe targeting.",
        "Unassessed": "Identity or measurement coverage is insufficient.",
        "Screening": "Other nomination, coverage or screening requirements remain unmet.",
    }
    return detail_text.get(detail, row["Reason"])


def normal_evidence(row, limit=300):
    if row["Highest normal nTPM context"] or row["Highest normal TPM context"]:
        return short(
            "; ".join(
                v
                for v in [row["Highest normal nTPM context"], row["Highest normal TPM context"]]
                if v
            ),
            limit,
        )
    return short(
        row["IHC evidence"]
        or row["RNA evidence"]
        or "Expanded maxima were not audited for this inherited exclusion; see original evidence in the table.",
        limit,
    )


def catalog(genes, category, filename, label):
    subset = genes.loc[genes.Category.eq(category)]
    doc = Report(OUT / filename, f"CTA atlas: {label}", catalog=True)
    story = cover(
        label,
        f"{len(subset):,} mapped genes | Alphabetical catalog with a PDF bookmark for every gene",
    )
    story += [
        paragraph(
            "Every entry carries its current classification, cancer occurrence, normal-expression evidence and curated biology. The complete tables retain unabridged annotations, all tumor cohorts, every normal context and source links."
        ),
        paragraph(
            "Cancer RNA summaries rank cohorts by the fraction >=5 pTPM (at least 20 measured samples). Cancer IHC summaries rank Medium/High patient counts (at least five scored patients). These are different cohorts and native measures, not assay sensitivity estimates. No values are imputed for missing genes."
        ),
        paragraph(
            "Biology excerpts come from reviewed UniProt entries (CC BY 4.0). They preserve qualifications but are not an independent review of every underlying paper. Clinical outcomes are reviewed only for the explicitly marked targets; a trial link alone establishes neither efficacy nor safety."
        ),
        paragraph(
            "Read alongside the overview for the threshold grid, priority reconsideration list, known therapeutic risks, assay limitations and unmapped BAGE/GAGE2B records."
        ),
        paragraph(
            link(BASE_URL + "explore.html", "Open the searchable atlas and full tables"), raw=True
        ),
        PageBreak(),
    ]
    for r in subset.to_dict("records"):
        header = paragraph(
            f"{r['Gene'] or r['Gene ID']}  |  {r['Category']} / {r['Detail']}", "Gene"
        )
        header.gene_id = r["Gene ID"]
        identity = paragraph(r["Gene ID"], "Small")
        identity.keepWithNext = True
        lines = [paragraph(f"<b>Decision.</b> {esc(reason_text(r))}", "Catalog", raw=True)]
        if r["Literature conflict"]:
            lines.append(
                paragraph(
                    f"<b>{esc(r['Literature conflict priority'])}.</b> {esc(r['Literature conflict'])} "
                    + link(r["Literature conflict source"], "Evidence")
                    + " Frozen classification unchanged; adjudicate before target use.",
                    "Catalog",
                    raw=True,
                )
            )
        if r["Boundary review"]:
            lines.append(
                paragraph(f"<b>Boundary.</b> {esc(r['Boundary review'])}", "Catalog", raw=True)
            )
        lines.append(
            paragraph(f"<b>Normal RNA / IHC.</b> {esc(normal_evidence(r))}", "Catalog", raw=True)
        )
        lines.append(
            paragraph(
                f"<b>Cancer RNA.</b> {esc(short(r['Cancer RNA at >=5 pTPM'], 250))}",
                "Catalog",
                raw=True,
            )
        )
        lines.append(
            paragraph(
                f"<b>Cancer IHC.</b> {esc(short(r['Cancer IHC Medium/High'], 190))}",
                "Catalog",
                raw=True,
            )
        )
        biology = (
            r["Biological function"]
            or r["Curated tissue biology"]
            or "No functional description is linked by an exact Ensembl cross-reference in this reviewed UniProt snapshot."
        )
        lines.append(paragraph(f"<b>Biology.</b> {esc(short(biology, 440))}", "Catalog", raw=True))
        if r["Curated tissue biology"] and r["Biological function"]:
            lines.append(
                paragraph(
                    f"<b>Curated tissue note.</b> {esc(short(r['Curated tissue biology'], 210))}",
                    "Catalog",
                    raw=True,
                )
            )
        if r["Review question"]:
            lines.append(
                paragraph(
                    f"<b>Reconsideration.</b> {esc(short(r['Reconsideration evidence'], 250))}",
                    "Catalog",
                    raw=True,
                )
            )
        evidence = r["Reviewed evidence types"] or r["Trial evidence"]
        lines.append(
            paragraph(f"<b>Therapeutic evidence.</b> {esc(evidence)}.", "Catalog", raw=True)
        )
        sources = [
            link(u.strip(), "UniProt " + a.strip())
            for u, a in zip(r["Biology sources"].split(" | "), r["UniProt entries"].split(" | "))
            if u.strip()
        ]
        sources += [
            link(u.strip(), f"Reviewed source {i + 1}")
            for i, u in enumerate(r["Reviewed literature sources"].split(" | "))
            if u.strip()
        ]
        if not r["Reviewed literature sources"]:
            sources += [
                link(u.strip(), "Trial/source record")
                for u in r["Trial sources"].split(" | ")
                if u.strip()
            ]
        sources += [
            link(
                f"https://www.proteinatlas.org/{r['Gene ID']}",
                "HPA portal (live; not the pinned values)",
            )
        ]
        lines.append(paragraph(" | ".join(sources), "Small", raw=True))
        if r["Shared UniProt entry"]:
            lines.append(
                paragraph(
                    "Shared UniProt entry: annotation is not independently proven for each encoding locus.",
                    "Small",
                )
            )
        lines.append(Spacer(1, 12))
        story += [KeepTogether([header, identity, *lines])]
    doc.build(story)
    return [{"Gene ID": gid, "PDF": filename, "Page": page} for gid, page in doc.gene_pages.items()]


def overview(genes, summary, reviews):
    story = cover("The new CTA set", "What survives, what fails, and what deserves another look")
    story += [
        table(
            ["Retained CTA candidates", "Neither", "Excluded / Unsafe"],
            [["399", "561", "1,580"]],
            [WIDTH / 3] * 3,
        ),
        Spacer(1, 12),
        paragraph(
            "The report covers all 2,540 mapped CTA nominations. CTA means a research-screen pass, Unsafe means exclusion from targeting and inclusion in the forbidden background, and Neither means neither selected nor added to that background. No category constitutes clinical clearance."
        ),
        paragraph(
            "The latest change excludes any observed RNA >30 TPM or >30 nTPM outside explicitly mapped Reproductive/Disposable contexts. It removes 37 previously retained CTAs and 253 previous Neither genes across the full canonical universe. Existing stricter RNA and protein rules remain."
        ),
        paragraph(
            "MAGEA3 and ROPN1 are the two new clinical/requested benchmark losses. CT55 stays retained. CTAG2 stays Variable; FOXR2 stays under antibody review. The review questions and exact reasons are preserved rather than waived for familiar targets."
        ),
        paragraph(
            f"New literature review flags {summary['coverage']['Retained genes with literature conflicts']} retained genes for normal-tissue conflicts. ABCC11 and MAGEC2 warrant urgent protein-evidence review before target use. These flags are shown beside the frozen decisions; they are not silently treated as clean passes."
        ),
        section("How to use this report"),
        table(
            ["Question", "Where to look"],
            [
                [
                    "Which genes made it?",
                    "Retained catalog: all 399; retained.csv for full evidence.",
                ],
                [
                    "Which could be reconsidered?",
                    "Priority follow-up and boundary sections; Neither catalog (561) and excluded catalog (1,580).",
                ],
                [
                    "How often are they expressed in cancer?",
                    "Cancer occurrence heatmap and tables; complete RNA/IHC cohort exports.",
                ],
                [
                    "Which have benefit or toxicity evidence?",
                    "Primary clinical literature and 30 clinical-target entries; registration is not an outcome.",
                ],
                [
                    "What is known about their biology?",
                    "Every gene has a UniProt annotation or an explicit mapping/function gap; source-linked full text is in all-genes.csv.gz.",
                ],
            ],
            [WIDTH * 0.4, WIDTH * 0.6],
        ),
        Spacer(1, 10),
        paragraph(
            link(BASE_URL + "explore.html", "Search all genes, evidence and sources online"),
            raw=True,
        ),
        PageBreak(),
    ]
    story += [
        section("The current rules"),
        paragraph(
            "Admission and normal counterevidence are separate. A nomination needs reproductive RNA >=5 nTPM, reliable positive reproductive protein, or cited reproductive plus tumor/recognition support. Good positive protein evidence can satisfy admission below 5 nTPM. None of these bypass normal-expression gates."
        ),
        table(
            ["Tissue category", "RNA ceilings, nTPM", "Selectivity minima"],
            [
                ["Vital", "2 / 4 / 6 / 12", "80 / 75 / 70 / 65%"],
                ["Protected", "6 / 8 / 12 / 24", "75 / 70 / 65 / 60%"],
                ["Manageable", "10 / 20 / 30 / 60", "70 / 65 / 60 / 55%"],
            ],
            [120, 200, WIDTH - 320],
        ),
        paragraph(
            "Order: missing/uncertain, Approved-negative, Supported-negative, Enhanced-negative. Relaxation requires complete local negative IHC and a reliable reproductive positive control. Any accepted positive normal IHC, including Low, excludes. Uncertain/cross-reactive antibodies are handled separately."
        ),
        paragraph(
            "For each bulk tissue, require N <= C, N <= R(1-s)/s + 0.1, and N < R. R is the maximum reproductive RNA, N the tissue RNA, C its ceiling and s its selectivity minimum. The fixed 0.1 nTPM allowance is not a fitted uncertainty interval. The new >30 rule independently caps every audited individual/cluster observation, regardless of enrichment or IHC."
        ),
        paragraph(
            "Normal-cell and immune data also trigger review at lower values: source-specific rules inspect clusters, DICE distributions, immune summaries, corneal CP10K and primary antibody conflicts. These can produce Neither even when the bulk RNA/IHC grid passes. CP10K is never compared numerically with TPM/nTPM."
        ),
        paragraph(
            "Reproductive and Disposable contexts are permitted by this research policy. Immune cells sampled inside those organs retain their immune identity. Unresolved anatomy does not establish permission. Human thresholds have not been validated for dogs or mice."
        ),
        section("Genes and proteins are different scopes"),
        table(
            ["Scope", "CTA", "Unsafe", "Neither"],
            [
                ["Mapped CTA nominations", "399", "1,580", "561"],
                ["Canonical human genes", "399", "17,699", "1,999"],
                ["Complete unique proteins", "723", "126,384", "2,878"],
            ],
            [WIDTH - 240, 80, 80, 80],
        ),
        paragraph(
            "Exact shared sequences use Unsafe > Neither > CTA. An additional 391 noncanonical reference loci have high-RNA exclusions. Gene expression does not establish which isoform is expressed; full-protein disjointness does not make their k-mers disjoint."
        ),
        PageBreak(),
    ]
    detail_names = {
        "Protein": "Positive normal IHC",
        "RNA": "Bulk RNA above applicable ceiling",
        "High RNA": "New high normal RNA rule",
        "Normal": "Published normal-cell protein",
        "Safety": "Reviewed safety hold",
        "Antibody": "Antibody evidence unresolved",
        "Cellular": "Cell / immune RNA conflict",
        "Low": "Weak affirmative reproductive evidence",
        "Screening": "Other screening requirements",
        "Selectivity": "Insufficient reproductive enrichment",
        "Unassessed": "Coverage / identity gap",
        "Variable": "Vital donor zero/high pattern",
    }
    story += [
        section("Where candidates leave the funnel"),
        paragraph(
            "Counts below are mutually exclusive final reason labels, not the number of every failing measurement. A gene can have several counterevidence findings even though it has one final label."
        ),
        table(
            ["Result", "Principal reason", "Genes"],
            [
                [r["Category"], detail_names.get(r["Detail"], r["Detail"]), str(r["Genes"])]
                for r in summary["details"]
                if r["Category"] != "CTA"
            ],
            [85, WIDTH - 145, 60],
        ),
        Spacer(1, 12),
        paragraph(
            "Neither is not a safe reserve list. Low means insufficient target evidence; Cellular and Antibody mean unresolved conflicts. The maximum audited non-permitted values remaining anywhere in canonical Neither are 29.9398 TPM and 30 nTPM. Missing measurements remain missing."
        ),
        PageBreak(),
    ]
    story += [
        section("Cancer occurrence: the measured landscape"),
        Image(str(figure_heatmap(genes)), width=WIDTH, height=WIDTH * 0.8),
        paragraph(
            "HPA v25.1-1 TCGA cohorts; RNA >=5 native pTPM. Gray cells are missing exact-ID coverage, not negatives. Label colors show CTA (teal), Neither (amber) and Unsafe (red). Full tables include both TCGA and independent validation cohorts and thresholds 0.1, 1, 2 and 5 pTPM."
        ),
        paragraph(
            "This does not estimate prevalence in all patients with cancer: cohorts differ in histology, stage, sampling and selection. A bulk tumor signal can arise in tumor-associated normal cells. Closely related genes may share sequence and antibody recognition. Gene-level expression cannot prove target-peptide presentation."
        ),
        PageBreak(),
    ]
    important = [
        "MAGEA4",
        "PRAME",
        "CTAG1A",
        "CTAG1B",
        "XAGE1A",
        "XAGE1B",
        "CT83",
        "CT55",
        "MAGEA3",
        "ROPN1",
        "CTAG2",
        "FOXR2",
    ]
    story += [
        section("Cancer occurrence for key targets"),
        table(
            ["Gene / result", "Highest ranked RNA cohorts (>=5 pTPM)", "Cancer IHC Medium/High"],
            [
                [
                    f"{r['Gene']}\n{r['Category']}",
                    short(r["Cancer RNA at >=5 pTPM"], 200),
                    short(r["Cancer IHC Medium/High"], 150),
                ]
                for r in genes.to_dict("records")
                if r["Gene"] in important
            ],
            [80, 265, WIDTH - 345],
        ),
        paragraph(
            "RNA and IHC are unpaired. Small cancer IHC groups (often about 8-12 patients) do not calibrate antibody sensitivity. The aggregate download lacks antibody identity and reliability; discordance should prompt assay review rather than numerical averaging. Synovial sarcoma trial enrichment is not represented by a dedicated HPA RNA cohort here."
        ),
        PageBreak(),
    ]
    story += [
        section("Retained genes with notable normal RNA"),
        paragraph(
            "These are the 15 retained genes with the largest audited non-permitted nTPM values. Native TPM is shown separately in the downloadable table. Single-cell nTPM is not directly comparable to bulk-tissue nTPM; the source and cell identity matter. A retained label does not mean no normal expression."
        ),
        table(
            ["Gene", "Highest normal observation", "Reproductive bulk nTPM"],
            [
                [
                    r["Gene"],
                    r["Highest normal nTPM context"],
                    f"{float(r['Reproductive RNA (nTPM)']):.3g}",
                ]
                for r in genes.loc[genes.Category.eq("CTA")]
                .sort_values("Normal RNA maximum (nTPM)", ascending=False)
                .head(15)
                .to_dict("records")
            ],
            [75, WIDTH - 160, 85],
        ),
        paragraph(
            "Examples near 30 are TEX14 (breast macrophages), MORC1 (PBMC platelets), and KIAA1210. Immune cells inside breast are not granted reproductive permission. These findings warrant source/identity/protein follow-up; the present report does not change their frozen decisions."
        ),
        PageBreak(),
    ]
    conflicts = yaml.safe_load((DATA / "retained-literature-conflicts.yml").read_text())["records"]
    story += [
        section("Retained does not resolve these literature conflicts"),
        paragraph(
            "The new UniProt review was screened manually for positive normal-tissue assertions, excluding sentences that merely mention cancer or negative expression. These are follow-up flags, not automatic rule changes. ABCC11 and MAGEC2 have primary normal neural protein reports and deserve urgent adjudication. Other entries are RNA-, isoform- or context-dependent and need their original assays checked."
        ),
        table(
            ["Gene(s)", "Priority / reported conflict"],
            [[", ".join(r["genes"]), r["priority"] + ": " + r["finding"]] for r in conflicts],
            [120, WIDTH - 120],
        ),
        paragraph(
            "Full citations and evidence required for each review are in retained-literature-conflicts.csv and in the relevant gene entries. No change to the 399-gene frozen set or forbidden sequences is implied by this reporting layer."
        ),
        PageBreak(),
    ]
    watch = genes.loc[genes["Review question"].ne("") & genes.Category.ne("CTA")]
    story += [
        section("Priority reconsideration queue"),
        paragraph(
            "Priority reflects prior biological/clinical interest and a specific resolvable question, not a numerical claim of safety. A positive result on one follow-up does not waive the remaining gates."
        ),
        table(
            ["Gene / result", "Question and evidence needed"],
            [
                [
                    f"{r['Gene']}\n{r['Category']} / {r['Detail']}",
                    r["Review question"] + " " + r["Reconsideration evidence"],
                ]
                for r in watch.to_dict("records")
            ],
            [95, WIDTH - 95],
        ),
        paragraph(
            "For CTAG2, pinned heart RNA donors are 13.5, 0, 0 and 7.1 nTPM; negative heart IHC comes from different specimens. Both nonzero donors fail the existing local rules. The expanded 24.9 nTPM maximum is from another normal context and does not replace the heart evidence. FOXR2's normal bulk RNA is low, but uncertain antibody staining still requires adjudication."
        ),
        PageBreak(),
    ]
    boundaries = genes.loc[genes["Boundary review"].ne("") & genes.Category.ne("CTA")]
    story += [
        section("Which exclusions sit near an arbitrary cutoff?"),
        paragraph(
            "These transparent review bands are descriptive: within 10% above the new 30 ceiling, or all failed bulk numeric comparisons within 20% of their limits. They are not new acceptance thresholds. RNA-only does not mean clinically safe, and many genes still have additional holds."
        ),
        table(
            ["Gene", "Result", "Boundary flag"],
            [
                [r["Gene"], r["Category"] + " / " + r["Detail"], r["Boundary review"]]
                for r in boundaries.to_dict("records")
            ],
            [80, 100, WIDTH - 180],
        ),
        PageBreak(),
    ]
    story += [
        section("What the clinical literature actually supports"),
        paragraph(
            "This is a targeted primary-source review, not a systematic review of all 2,540 genes. Outcomes below retain the study/product context. Trials with multiple antigens cannot attribute benefit to every component. Reported CRS, conditioning-related cytopenias and receptor cross-reactivity are distinct from demonstrated on-target injury to normal tissue."
        ),
    ]
    for review in reviews:
        if review["evidence"] in {
            "Normal function",
            "Cancer mechanism",
            "Human tumor recognition",
            "Tumor expression and immunity",
        }:
            continue
        story += [
            section(" / ".join(review["genes"]) + " | " + review["evidence"], review["id"]),
            paragraph(review["finding"]),
            paragraph(link(review["source"], review["title"]), "Small", raw=True),
        ]
    story += [PageBreak(), section("Biological promise is a separate question")]
    for review in reviews:
        if review["evidence"] not in {
            "Normal function",
            "Cancer mechanism",
            "Human tumor recognition",
            "Tumor expression and immunity",
        }:
            continue
        story += [
            section(" / ".join(review["genes"]) + " | " + review["evidence"], review["id"]),
            paragraph(review["finding"]),
            paragraph(link(review["source"], review["title"]), "Small", raw=True),
        ]
    story += [
        PageBreak(),
        section("All clinical-target intake records"),
        paragraph(
            "The 30 curated trial/therapy-linked genes are shown whether retained or excluded. A registry entry is a traceable nomination, not a demonstrated response, approved use, current recruitment status or absence of toxicity."
        ),
        table(
            ["Gene", "Current decision", "Evidence recorded"],
            [
                [
                    r["Gene"],
                    r["Category"] + " / " + r["Detail"],
                    r["Reviewed evidence types"] or r["Trial evidence"],
                ]
                for r in genes.loc[genes["Trial sources"].ne("")].to_dict("records")
            ],
            [85, 120, WIDTH - 205],
        ),
        PageBreak(),
    ]
    story += [
        section("Coverage, gaps and interpretation"),
        table(
            ["Evidence layer", "Genes covered / 2,540"],
            [[k, f"{v:,}"] for k, v in summary["coverage"].items()],
            [WIDTH - 120, 120],
        ),
        paragraph(
            "UniProt annotation is available through an exact Ensembl cross-reference for 2,505 genes, but only 1,873 have a function annotation. No function annotation is not evidence of no function. Shared UniProt records remain flagged. Cancer mechanisms and normal physiological roles are reproduced as curated assertions, preserving 'may' and 'by similarity' qualifiers."
        ),
        paragraph(
            "BAGE and GAGE2B remain in intake but outside the 2,540 mapped-gene denominator. Their historical antigen identities and absent matched reference/HPA rows need explicit locus/sequence resolution. No paralog or zero expression is substituted. SPAG9 is mapped and excluded for clear normal RNA/protein expression."
        ),
        paragraph(
            "Earlier Unsafe genes were not all re-audited against every expanded normal source. A blank normal maximum for these genes is an audit-coverage gap, not a negative result. Source cohorts overlap (HPA summaries can contain the same DICE/immune data); repetitions are not independent biological validation."
        ),
        paragraph(
            "The current tissue permission policy is a research selection assumption, not a clinical statement that reproductive or disposable-tissue injury is harmless. Pregnancy- and fertility-specific subsets require their separate anatomy definitions. Neither absence of toxicity in a small trial nor negative IHC proves universal safety."
        ),
        section("Reproduce and inspect"),
        paragraph("python scripts/build_cta_atlas.py\npython scripts/render_cta_atlas.py", "Code"),
        paragraph(
            "from oncoref import protein_annotations, reference_data\nprotein_annotations.gene_annotations(['ENSG00000147381'])\nreference_data.provenance('uniprot_reviewed_human')",
            "Code",
        ),
        paragraph(
            "Shared sources are checksum pinned: HPA normal v23; HPA cancer v25.1-1; UniProt 2026_03. The report adds evidence and navigation without changing policy 2026-10-24.1. Human data and thresholds are not automatically transferred across species."
        ),
        paragraph(
            link(BASE_URL + "provenance.json", "Source provenance and checksums")
            + " | "
            + link("https://github.com/pirl-unc/oncoref/pull/619", "Merged policy review PR #619")
            + " | "
            + link(
                "https://github.com/pirl-unc/oncoref/releases/tag/v1.18.0", "Pinned protein assets"
            ),
            "Small",
            raw=True,
        ),
    ]
    Report(OUT / "cta-overview.pdf", "The new CTA set").build(story)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(DATA / "all-genes.csv.gz", low_memory=False).fillna("")
    summary = json.loads((DATA / "summary.json").read_text())
    reviews = yaml.safe_load((DATA / "literature.yml").read_text())["records"]
    overview(genes, summary, reviews)
    pages = []
    for category, filename, title in [
        ("CTA", "cta-retained-catalog.pdf", "Retained CTA candidates"),
        ("Neither", "cta-neither-catalog.pdf", "Neither: unresolved and below-threshold genes"),
        ("Unsafe", "cta-excluded-catalog.pdf", "Excluded: Unsafe research classifications"),
    ]:
        pages.extend(catalog(genes, category, filename, title))
    pd.DataFrame(pages).to_csv(DATA / "pdf-gene-index.csv", index=False)
    print(json.dumps({p.name: p.stat().st_size for p in OUT.glob("cta-*.pdf")}, indent=2))


font_dir = Path(matplotlib.get_data_path()) / "fonts/ttf"
pdfmetrics.registerFont(TTFont("Atlas", str(font_dir / "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("AtlasBold", str(font_dir / "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFontFamily(
    "Atlas", normal="Atlas", bold="AtlasBold", italic="Atlas", boldItalic="AtlasBold"
)
STYLES = getSampleStyleSheet()
for name, size, leading, color, extra in [
    ("Body", 9.4, 14, INK, {"spaceAfter": 9}),
    ("Deck", 12, 18, MUTED, {"spaceAfter": 12}),
    (
        "Section",
        14,
        19,
        INK,
        {"spaceBefore": 12, "spaceAfter": 9, "fontName": "AtlasBold", "keepWithNext": True},
    ),
    (
        "Gene",
        13,
        18,
        INK,
        {"spaceBefore": 8, "spaceAfter": 2, "fontName": "AtlasBold", "keepWithNext": True},
    ),
    ("Catalog", 8.7, 12.3, INK, {"spaceAfter": 4}),
    ("Small", 7.7, 10.5, MUTED, {"spaceAfter": 5, "keepWithNext": False}),
    ("Table", 8.4, 11.5, INK, {}),
    ("TableHead", 8.5, 12, colors.white, {"fontName": "AtlasBold"}),
]:
    STYLES.add(
        ParagraphStyle(
            name,
            fontName="Atlas",
            fontSize=size,
            leading=leading,
            textColor=color,
            alignment=TA_LEFT,
            **extra,
        )
    ) if "fontName" not in extra else STYLES.add(
        ParagraphStyle(
            name, fontSize=size, leading=leading, textColor=color, alignment=TA_LEFT, **extra
        )
    )
STYLES["Title"].fontName = "AtlasBold"
STYLES["Title"].fontSize = 29
STYLES["Title"].leading = 35
STYLES["Title"].alignment = TA_LEFT
STYLES["Title"].textColor = INK
STYLES["Title"].spaceAfter = 14
STYLES["Code"].fontName = "Atlas"
STYLES["Code"].fontSize = 8
STYLES["Code"].leading = 12
STYLES["Code"].spaceAfter = 10

if __name__ == "__main__":
    main()
