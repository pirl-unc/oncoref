"""Figures, PDF and web report for the tiered permitted-tissue CTA screen."""

from __future__ import annotations

from xml.sax.saxutils import escape

import numpy as np

from oncoref import figure_style
from oncoref.normal_tissues import normal_tissue_screening_policy

NAMES = {
    "previous_1": "Prior scope | all outside RNA <=1",
    "expanded_1": "Expanded scope | all outside RNA <=1",
    "tiered": "Expanded scope | vital <=1, protected <=5",
}
REFERENCES = {
    "HPA v23 bulk RNA": "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
    "HPA v23 normal IHC": "https://v23.proteinatlas.org/download/normal_tissue.tsv.zip",
    "NCI organ inflammation": "https://www.cancer.gov/about-cancer/treatment/side-effects/organ-inflammation",
    "NIDDK adrenal insufficiency": "https://www.niddk.nih.gov/health-information/endocrine-diseases/adrenal-insufficiency-addisons-disease/definition-facts",
    "PRAME clinical TCR study": "https://doi.org/10.1038/s41591-025-03650-6",
    "KK-LC-1 clinical trial": "https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2022-07699",
    "ROPN1 preclinical TCR study": "https://doi.org/10.1158/2159-8290.CD-24-0168",
}
METHOD = (
    "All 2,537 protein-coding candidates from the 11 registered publication sources are screened. "
    "Each valid RNA measurement is compared with its tissue-specific ceiling using raw HPA v23 nTPM: "
    "vital <=1; protected <=5. There is no averaging across tissues or subtraction of a noise floor. "
    "Any Low, Medium, High, Ascending or Descending IHC observation outside the permitted set blocks retention. "
    "At least one permitted RNA observation >=1 nTPM or positive permitted IHC is required. "
    "All 50 routine bulk RNA contexts must be measured. Unmapped or tier-unresolved positive anatomy, "
    "unusable measurements and previously reviewed unreliable IHC prevent retention. "
    "Family exclusions and existing nomination/manual reviews also apply. "
    "Missing IHC is explicitly RNA-only observed compatibility; partial IHC remains partial. "
    "Neither establishes protein absence. Complete routine IHC means all 80 specified tissue/cell contexts. "
    "Routine contexts are those represented for at least half the genes in their source; the manifest lists them."
)
SCOPE = (
    "Permitted: all reproductive organs, mammary, placenta/extraembryonic and explicit embryonic/fetal contexts; "
    "cutaneous melanocytes and hair follicles; plus whole skin, salivary glands and thymus. "
    "Whole-skin permission includes its mapped cellular compartments. Generic unsited cell types, eye, retina "
    "and inner ear do not become skin. This is an explicit project screening allowance."
)
COVERAGE = (
    "The tier assignment is a conservative project policy, not a clinically validated list of expendable organs. "
    "The clinical references motivate review of organ injury; they do not establish the 1/5 nTPM thresholds. "
    "The numerical analysis is human HPA v23 bulk RNA and normal IHC only. Anatomy transfers to dog and mouse; "
    "these human RNA thresholds and tolerability assumptions do not. Single-cell RNA, tumor expression and "
    "peptide-HLA presentation are not evaluated here. Developmental tissues absent from HPA remain unassessed."
)


def short_reason(row):
    if row.retained:
        return f"Included; IHC {row.ihc_coverage.replace('_', ' ')}"
    pieces = []
    for text in (row.rna_exclusions, row.ihc_exclusions):
        if text:
            values = text.split("; ")
            pieces.append(
                "; ".join(values[:2]) + (f"; +{len(values) - 2} more" if len(values) > 2 else "")
            )
    for flag in row.failure_reasons.split(";"):
        if flag not in {"rna_above_tissue_ceiling", "positive_ihc_outside_permitted"}:
            pieces.append(flag.replace("_", " "))
    return "; ".join(pieces)


def figures(out, sources, summary, targets):
    figure_style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    def save(fig, stem):
        for suffix in ("png", "svg"):
            figure_style.save(fig, out / f"{stem}.{suffix}", keep_titles=True)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4.6))
    left = np.zeros(len(summary))
    for field, label, color in (
        ("without_ihc", "IHC unavailable", "#a7adb5"),
        ("partial_ihc", "Partial routine IHC", "#377eb8"),
        ("complete_routine_ihc", "Complete routine IHC", "#237a57"),
    ):
        values = summary[field].to_numpy()
        ax.barh(np.arange(3), values, left=left, color=color, height=0.55, label=label)
        for i, value in enumerate(values):
            if value >= 25:
                ax.text(left[i] + value / 2, i, str(value), ha="center", va="center", color="white")
        left += values
    for i, value in enumerate(left):
        ax.text(value + 8, i, str(int(value)), va="center", fontweight="bold")
    ax.set_yticks(np.arange(3), [NAMES[s] for s in summary.scenario])
    ax.invert_yaxis()
    ax.set_xlim(0, left.max() * 1.1)
    ax.set_xlabel("Retained published candidates")
    ax.set_title("Scope expansion and RNA ceilings | IHC coverage stays explicit")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, -0.5), ncol=3, fontsize=8)
    fig.subplots_adjust(left=0.43, right=0.97, top=0.86, bottom=0.28)
    save(fig, "policy-comparison")

    source = sources.sort_values("source_identities", ascending=False)
    fig, ax = plt.subplots(figsize=(10.5, 7.5))
    y = np.arange(len(source))
    for j, (key, color) in enumerate(zip(NAMES, ("#a7adb5", "#377eb8", "#237a57"))):
        bars = ax.barh(y + (j - 1) * 0.24, source[key], height=0.22, color=color, label=NAMES[key])
        ax.bar_label(bars, padding=2, fontsize=7)
    ax.set_yticks(y, source.source)
    ax.invert_yaxis()
    ax.set_xlim(0, source.tiered.max() * 1.15)
    ax.set_xlabel("Retained genes within each published source (sources overlap)")
    ax.set_title("Published nomination sources | retention under each scenario")
    ax.legend(loc="lower right", fontsize=8)
    fig.subplots_adjust(left=0.22, right=0.97, top=0.92, bottom=0.1)
    save(fig, "source-retention")

    focus = targets.head(25)
    values = np.array(
        [[int(row[f"retained_{key}"]) for key in NAMES] for _, row in focus.iterrows()]
    )
    fig, ax = plt.subplots(figsize=(9, 10))
    ax.imshow(values, cmap=ListedColormap(["#c66a60", "#237a57"]), vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(np.arange(len(focus)), focus.name)
    ax.set_xticks(
        np.arange(3),
        ["Prior scope\n1 nTPM", "Expanded scope\n1 nTPM", "Expanded scope\n1 / 5 nTPM"],
    )
    ax.xaxis.tick_top()
    for i in range(len(focus)):
        for j in range(3):
            ax.text(
                j,
                i,
                "Keep" if values[i, j] else "Out",
                color="white",
                ha="center",
                va="center",
                fontsize=9,
            )
    ax.set_title("Selected published targets | all evidence and curation gates", pad=45)
    ax.legend(
        handles=[
            Patch(color="#237a57", label="Retained"),
            Patch(color="#c66a60", label="Excluded or unresolved"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.09),
        ncol=2,
    )
    fig.subplots_adjust(left=0.28, right=0.95, top=0.86, bottom=0.1)
    save(fig, "published-targets")


def render(out, candidates, sources, summary, targets):
    figures(out, sources, summary, targets)
    config = normal_tissue_screening_policy()
    counts = summary.set_index("scenario")
    final = counts.loc["tiered"]
    targets = targets.copy()
    targets["short_reason"] = targets.apply(short_reason, axis=1)
    markdown_names = {key: name.replace("|", "&#124;") for key, name in NAMES.items()}
    table = "\n".join(
        f"| {markdown_names[r.scenario]} | {r.retained} | {r.without_ihc} | {r.partial_ihc} | {r.complete_routine_ihc} |"
        for r in summary.itertuples(index=False)
    )
    target_table = "\n".join(
        f"| {r.name} (`{r.Symbol}`) | {'Included' if r.retained else 'Excluded / unresolved'} | {r.short_reason} |"
        for r in targets.itertuples(index=False)
    )
    vital = ", ".join(name.replace("_", " ") for name in config["vital"]["roots"])
    protected = ", ".join(name.replace("_", " ") for name in config["protected"]["roots"])
    refs = "\n".join(f"- [{name}]({url})" for name, url in REFERENCES.items())
    text = f"""# CTA screen: expanded permitted tissues and protected organs

**Human HPA v23 · policy `{config["policy_id"]}`**

[Complete PDF report](cta-permitted-tissue-report.pdf) · [Tissue definitions](../../normal-tissue-scopes.md) · [Original source-intake report](../cta-tissue-scopes-20261019/index.md)

## Current rule

{SCOPE}

| Tissue tier | RNA rule | IHC rule |
|---|---|---|
| Permitted | No upper ceiling | Positive staining allowed |
| Vital | Every measured value ≤1 nTPM | Any positive staining excludes, including Low |
| Protected normal | Every measured value ≤5 nTPM | Any positive staining excludes, including Low |
| Unresolved anatomy | Positive evidence prevents retention | Positive evidence prevents retention |

**Vital roots:** {vital}. Brain includes all mapped brain regions; spinal cord is explicit. Organ subcompartments inherit their organ's tier.

**Protected roots:** {protected}. Retina includes retinal pigment epithelium. Head, foot and generic soft tissue remain unresolved because these broad labels can span tiers. No unlisted tissue silently receives the 5 nTPM allowance.

## Results and IHC coverage

| Scenario | Retained | No usable IHC | Partial routine IHC | Complete routine IHC |
|---|---:|---:|---:|---:|
{table}

The selected policy retains **{int(final.retained)}** candidates. **{int(final.no_forbidden_ihc_contexts_observed)}** have no valid IHC observations in forbidden contexts; these are not IHC-negative findings. Requiring all routine IHC contexts would retain **{int(final.complete_routine_ihc)}**. Adequate evidence of absence remains unestablished.

The prior-scope comparison is reproductive + disposable with a uniform 1 nTPM ceiling and the same positive-in-scope requirement. Expansion adds thymus, whole skin and salivary gland. The final step changes only the protected-tissue RNA ceiling from 1 to 5. The historical adaptive panel and zero-positive strict screen remain separately documented in the original report.

![Policy counts and coverage](policy-comparison.png)

[Vector figure](policy-comparison.svg)

## Published source retention

All 11 registered published sources enter the same 2,537-gene protein-coding universe. Source nominations overlap; counts cannot be summed. Original DOI/table/row mappings and targeted source selections remain in the [source-intake audit](../cta-tissue-scopes-20261019/index.md#published-nomination-intake).

![Source retention](source-retention.png)

[Vector figure](source-retention.svg)

## Known and published targets

This selected list mixes clinically studied and research targets; inclusion does not imply equal clinical maturity. PRAME and KK-LC-1 have clinical TCR studies, whereas the cited ROPN1 work is preclinical. Every published candidate, not just these examples, has a decision in the download below.

![Published target comparison](published-targets.png)

[Vector figure](published-targets.svg)

| Target | Result | Evidence / reason |
|---|---|---|
{target_table}

MAGE-A11's existing manual holdback remains separate from the expression gates. Reviewed-unreliable IHC remains an evidence-quality exclusion, with raw staining preserved. Tissue maxima are gene-level source observations; closely related paralogs and shared assays retain the original source limitations.

## Methods and limits

{METHOD}

{COVERAGE}

## Downloads and reproduction

- [All candidate decisions](candidate-decisions.csv)
- [Selected published targets and full exclusion reasons](notable-targets.csv)
- [Source counts](source-counts.csv) and [overall counts](scope-summary.csv)
- [Original observations with current mappings and tiers](normal-observations.csv.gz)
- [Exact measured tissue/cell mappings](tissue-contexts.csv)
- [Complete taxonomy and policy](tissue-policy.json)
- [Input, code and output checksums; exact routine contexts](run-manifest.json)

```sh
python scripts/cta_permitted_report.py --out outputs/cta-permitted-current
```

The rebuild verifies archived input hashes, reproduces the previous 547-gene sensitivity, and checks that scope expansion and RNA-ceiling relaxation cannot lose candidates. It uses the committed source audit and does not require live downloads.

## References

{refs}
"""
    (out / "index.md").write_text(text)
    pdf(out, summary, targets, vital, protected)


def pdf(out, summary, targets, vital, protected):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.platypus import (
        Image,
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="SmallBody", fontName="Helvetica", fontSize=8.5, leading=11, spaceAfter=6
        )
    )
    styles.add(ParagraphStyle(name="CellBody", fontName="Helvetica", fontSize=8, leading=10))
    styles["Heading2"].keepWithNext = True
    story = []

    def paragraph(text, style="SmallBody"):
        return Paragraph(escape(str(text)), styles[style])

    def table(rows, widths):
        result = Table(
            [[paragraph(cell, "CellBody") for cell in row] for row in rows],
            colWidths=widths,
            repeatRows=1,
            hAlign="LEFT",
        )
        result.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8edf2")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#9ca3af")),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 6),
                    (
                        "ROWBACKGROUNDS",
                        (0, 1),
                        (-1, -1),
                        [colors.white, colors.HexColor("#f6f8fa")],
                    ),
                ]
            )
        )
        return result

    story += [
        paragraph("CTA screen: expanded permitted tissues", "Title"),
        paragraph("Human HPA v23 | vital RNA <=1 nTPM; protected RNA <=5 nTPM", "Heading2"),
        paragraph(SCOPE),
    ]
    story += [
        paragraph("Exact rules", "Heading2"),
        paragraph(METHOD),
        paragraph(
            "Vital roots: "
            + vital
            + ". All mapped brain regions and organ descendants inherit these roots."
        ),
        paragraph(
            "Protected roots: "
            + protected
            + ". Head, foot and generic soft tissue remain unresolved."
        ),
        paragraph(COVERAGE),
    ]
    story += [
        paragraph("Results", "Heading2"),
        table(
            [["Scenario", "Retained", "No IHC", "Partial IHC", "Complete IHC"]]
            + [
                [
                    NAMES[r.scenario],
                    r.retained,
                    r.without_ihc,
                    r.partial_ihc,
                    r.complete_routine_ihc,
                ]
                for r in summary.itertuples(index=False)
            ],
            [270, 60, 60, 65, 65],
        ),
    ]
    story += [
        PageBreak(),
        paragraph("Policy comparison", "Heading2"),
        Image(str(out / "policy-comparison.png"), width=520, height=218),
        Spacer(1, 12),
        paragraph("Published source retention", "Heading2"),
        paragraph(
            "All 11 registered sources enter the same 2,537-gene protein-coding universe. Sources overlap and cannot be summed."
        ),
        Image(str(out / "source-retention.png"), width=520, height=371),
    ]
    story += [
        PageBreak(),
        paragraph("Selected published targets", "Heading2"),
        paragraph(
            "Green means retained under the stated evidence and curation gates. Missing IHC remains unassessed. All exact exclusion reasons follow in the table."
        ),
        Image(str(out / "published-targets.png"), width=500, height=556),
    ]
    story += [
        PageBreak(),
        paragraph("Known and published targets: decisions", "Heading2"),
        paragraph(
            "Selected clinical and research targets. Exclusions summarize up to two observations per assay; the CSV contains every failing context. Existing manual reviews remain separate from the new expression gates."
        ),
    ]
    story.append(
        table(
            [["Target", "Result", "Evidence / reason"]]
            + [
                [r.name, "Included" if r.retained else "Out / unresolved", r.short_reason]
                for r in targets.itertuples(index=False)
            ],
            [100, 80, 340],
        )
    )
    story.append(paragraph("Evidence and citations", "Heading2"))
    for name, url in REFERENCES.items():
        story.append(
            Paragraph(
                f'<link href="{escape(url)}" color="#205c8f">{escape(name)}</link>',
                styles["SmallBody"],
            )
        )
    story.append(
        Paragraph(
            '<link href="https://pirl-unc.github.io/oncoref/audits/cta-tissue-scopes-20261019/" color="#205c8f">Original source-intake audit and publication citations</link>',
            styles["SmallBody"],
        )
    )

    def footer(canvas, doc):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#64748b"))
        canvas.drawString(44, 24, "OncoRef | 2026-10-19 | Human HPA v23 | Project screening policy")
        canvas.drawRightString(568, 24, str(doc.page))

    SimpleDocTemplate(
        str(out / "cta-permitted-tissue-report.pdf"),
        pagesize=letter,
        leftMargin=44,
        rightMargin=44,
        topMargin=38,
        bottomMargin=40,
        title="CTA expanded permitted-tissue screening report",
    ).build(story, onFirstPage=footer, onLaterPages=footer)
