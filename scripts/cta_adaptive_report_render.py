"""Web, figures and PDF for the IHC-conditioned absolute RNA screen."""

from __future__ import annotations

from xml.sax.saxutils import escape

import numpy as np

from oncoref import figure_style
from oncoref.normal_tissues import normal_tissue_screening_policy

NAMES = {"fixed": "Fixed RNA ceilings (1 / 5)", "adaptive": "Tissue-specific negative IHC"}
REFERENCES = {
    "HPA v23 consensus RNA input": "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
    "HPA v23 normal IHC input": "https://v23.proteinatlas.org/download/normal_tissue.tsv.zip",
    "HPA normalization and maximum-source consensus": "https://www.proteinatlas.org/humanproteome/brain/method/transcriptomics",
    "CTAG2 current HPA and GTEx heart data": "https://www.proteinatlas.org/ENSG00000126890-CTAG2/tissue/heart%2Bmuscle",
    "PRAME current adrenal data": "https://www.proteinatlas.org/ENSG00000185686-PRAME/tissue/adrenal%2Bgland",
    "IMA203 phase 1 primary report": "https://doi.org/10.1038/s41591-025-03650-6",
}
METHOD = (
    "The same 2,537 published candidates and permitted tissue union are evaluated in both "
    "screens. The fixed screen uses vital <=1 nTPM and protected <=5 nTPM. The adaptive "
    "screen increases those ceilings only when all routine IHC tissue/cell contexts for "
    "the same anatomical tissue have valid Not detected results. The routine inventory "
    "is frozen from contexts measured for at least half of the source's measured genes. "
    "The weakest reliability "
    "among all IHC observations in that tissue determines the ceiling. Extra positive or "
    "invalid observations prevent relaxation. Missing contexts, unknown reliability, "
    "Uncertain reliability and reviewed-unreliable IHC cannot increase a ceiling. "
    "A narrow brain region cannot clear a broader RNA region. Negative evidence in another "
    "organ cannot supply support. These are tissue-specific observation rules, not claims "
    "of complete biological sampling."
)
GATES = (
    "Any positive IHC outside the permitted set, including Low, remains an independent "
    "exclusion. Complete routine bulk RNA (50 contexts), positive evidence in the permitted "
    "set (RNA >=1 nTPM or positive IHC), resolved positive anatomy and usable observations "
    "remain required. Existing family exclusions and manual nomination decisions remain "
    "separate gates. Complete routine IHC means 80 tissue/cell contexts; missing and partial "
    "coverage remain explicit. RNA is pinned human HPA v23 consensus nTPM with inclusive "
    "ceilings, no rounding to pass a cutoff, and no subtraction or averaging across organs."
)
SOURCE_NOTE = (
    "HPA is both a study and a portal that displays other studies. HPA consensus selects "
    "the maximum tissue estimate from HPA and GTEx; it is not evidence that both datasets "
    "agree. The pinned v23 CTAG2 heart value used here is 5.1 nTPM. The previously reviewed "
    "v23 HPA donors were 13.5, 7.1, 0 and 0 nTPM. The current HPA page, checked 2026-10-20, "
    "reports its own cohort at 5.0, versus HPA-normalized GTEx at 0.4 (atrial appendage) and "
    "0.9 (left ventricle). These are different datasets and release-specific estimates. "
    "Current page values do not replace the pinned input. CTAG2 heart IHC is Not detected "
    "with Enhanced reliability, so the new ceiling is 5; pinned 5.1 still fails. "
    "PRAME adrenal RNA is 3.6 in the pinned consensus and its IHC is Not detected, Enhanced. "
    "The current adrenal page separates HPA 0.5 from GTEx 3.6."
)
LIMITS = (
    "The 1/5, 2/5, 3/7.5 and 5/10 nTPM values are user-approved starting project parameters, "
    "not clinically validated safety limits. Enhanced is an HPA validation label, not a "
    "guaranteed detection limit for a negative stain. Routine negative IHC does not establish "
    "absence of protein or peptide-HLA presentation. Anatomical definitions apply across "
    "human, dog and mouse; these numerical rules apply only to human HPA v23. The existing "
    "core/extended CTA API panels retain their separately documented selection rules."
)


def reason(row):
    if row.retained:
        if not row.fixed_retained:
            values = row.rescued_rna_contexts.split("; ")
            # Preserve the tissue value and its IHC support together.
            return (
                "Recovered: "
                + "; ".join(values[:4])
                + ("; further contexts in CSV" if len(values) > 4 else "")
            )
        return "Retained; IHC " + row.ihc_coverage.replace("_", " ")
    pieces = []
    for text in (row.rna_exclusions, row.ihc_exclusions):
        if text:
            values = text.split("; ")
            pieces.append(
                "; ".join(values[:2]) + (f"; +{len(values) - 2} more" if len(values) > 2 else "")
            )
    for flag in row.failure_reasons.split(";"):
        if flag == "existing_nomination_or_manual_review":
            pieces.append("Existing review: " + row.nomination_review_reason)
        elif flag not in {"rna_above_tissue_ceiling", "positive_ihc_outside_permitted"}:
            pieces.append(flag.replace("_", " "))
    return "; ".join(pieces)


def figures(out, summary, sources, targets):
    figure_style.apply()
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    def save(fig, name):
        for suffix in ("png", "svg"):
            figure_style.save(fig, out / f"{name}.{suffix}", keep_titles=True)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    left = np.zeros(len(summary))
    for field, label, color in (
        ("without_ihc", "IHC unavailable", "#a7adb5"),
        ("partial_ihc", "Partial routine IHC", "#377eb8"),
        ("complete_routine_ihc", "Complete routine IHC", "#237a57"),
    ):
        values = summary[field].to_numpy()
        ax.barh(np.arange(len(summary)), values, left=left, label=label, color=color, height=0.5)
        for i, value in enumerate(values):
            if value >= 25:
                ax.text(left[i] + value / 2, i, str(value), color="white", ha="center", va="center")
        left += values
    for i, value in enumerate(left):
        ax.text(value + 8, i, str(int(value)), va="center", fontweight="bold")
    ax.set_yticks(np.arange(len(summary)), [NAMES[k] for k in summary.scenario])
    ax.set_xlim(0, left.max() * 1.12)
    ax.invert_yaxis()
    ax.set_xlabel("Retained published candidates")
    ax.set_title("Absolute RNA ceilings conditioned on local negative IHC")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, -0.55), ncol=3, fontsize=8)
    fig.subplots_adjust(left=0.32, right=0.97, top=0.86, bottom=0.3)
    save(fig, "policy-comparison")

    source = sources.sort_values("source_identities", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 7))
    y = np.arange(len(source))
    for i, (key, color) in enumerate(zip(NAMES, ("#a7adb5", "#237a57"))):
        bars = ax.barh(y + (i - 0.5) * 0.3, source[key], height=0.28, color=color, label=NAMES[key])
        ax.bar_label(bars, padding=2, fontsize=7)
    ax.set_yticks(y, source.source)
    ax.invert_yaxis()
    ax.set_xlim(0, source.adaptive.max() * 1.15)
    ax.set_xlabel("Retained genes per published source; overlapping nominations")
    ax.set_title("Published source retention")
    ax.legend(loc="lower right", fontsize=8)
    fig.subplots_adjust(left=0.23, right=0.97, top=0.91, bottom=0.1)
    save(fig, "source-retention")

    values = np.array([[int(r.fixed_retained), int(r.retained)] for r in targets.itertuples()])
    fig, ax = plt.subplots(figsize=(8, 11))
    ax.imshow(values, cmap=ListedColormap(["#c66a60", "#237a57"]), vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(np.arange(len(targets)), targets.name, fontsize=8)
    ax.set_xticks([0, 1], ["Fixed 1 / 5 nTPM", "IHC-adaptive ceilings"])
    ax.xaxis.tick_top()
    for i in range(len(targets)):
        for j in range(2):
            ax.text(
                j,
                i,
                "Keep" if values[i, j] else "Out",
                color="white",
                ha="center",
                va="center",
                fontsize=7,
            )
    ax.set_title("42 published targets | all evidence and curation gates", pad=35)
    fig.subplots_adjust(left=0.3, right=0.96, top=0.91, bottom=0.03)
    save(fig, "published-targets")


def render(out, summary, sources, targets):
    targets = targets.copy()
    targets["reason"] = targets.apply(reason, axis=1)
    figures(out, summary, sources, targets)
    config = normal_tissue_screening_policy("ihc_adaptive_cta")
    final = summary.set_index("scenario").loc["adaptive"]
    ladder = [["Local negative IHC reliability", "Vital maximum nTPM", "Protected maximum nTPM"]]
    for key in ("Missing", "Uncertain", "Approved", "Supported", "Enhanced"):
        limits = config["negative_ihc_rna_limits"][key]
        ladder.append([key, f"{limits['vital']:g}", f"{limits['protected']:g}"])
    counts = [["Scenario", "Retained", "No IHC", "Partial IHC", "Complete IHC"]]
    for row in summary.itertuples():
        counts.append(
            [
                NAMES[row.scenario],
                row.retained,
                row.without_ihc,
                row.partial_ihc,
                row.complete_routine_ihc,
            ]
        )
    target_rows = [["Target", "Fixed", "Adaptive", "Evidence / reason"]]
    for row in targets.itertuples():
        target_rows.append(
            [
                f"{row.name} ({row.Symbol})",
                "Keep" if row.fixed_retained else "Out",
                "Keep" if row.retained else "Out",
                row.reason,
            ]
        )

    def md_table(rows):
        rendered = [
            "| " + " | ".join(str(cell).replace("|", "&#124;") for cell in row) + " |"
            for row in rows
        ]
        rendered.insert(1, "|" + "---|" * len(rows[0]))
        return "\n".join(rendered)

    refs = "\n".join(f"- [{name}]({url})" for name, url in REFERENCES.items())
    text = f"""# CTA screen: IHC-adaptive absolute RNA ceilings

**Human HPA v23 · policy `{config["policy_id"]}` · 2026-10-20**

[PDF report](cta-ihc-adaptive-report.pdf) · [Tissue definitions](../../normal-tissue-scopes.md) · [Fixed-screen report](../cta-permitted-tissues-20261019/index.md)

## Approved starting thresholds

{md_table(ladder)}

{METHOD}

{GATES}

{LIMITS}

## Results

{md_table(counts)}

The adaptive policy retains **{final.retained}** candidates, including **{final.recovered}** recovered from the fixed screen. **{final.without_ihc}** retained candidates have no usable IHC and receive no IHC-based relaxation. Coverage is a separate result; routine completeness does not establish biological absence.

![Counts and IHC coverage](policy-comparison.png)

[Vector figure](policy-comparison.svg)

## Published source retention

![Overlapping publication sources](source-retention.png)

[Vector figure](source-retention.svg)

## CTAG2: HPA and GTEx are different datasets

{SOURCE_NOTE}

See the [earlier CTAG2 donor, isoform and clinical review](../ctag2-heart-evidence.md). Bulk RNA and IHC donors differ. Neither the low GTEx central estimate nor negative IHC resolves the expressing cell type or proves absence of peptide presentation.

## Known and published targets

![42 targets under both policies](published-targets.png)

[Vector figure](published-targets.svg)

{md_table(target_rows)}

The selected list contains both clinical and research targets. Manual exclusions remain independent: MAGE-A11 can pass expression gates while retaining its earlier nomination holdback. CTAG2 also retains its prior review decision. The complete CSV preserves every failing tissue and every rescued RNA context.

## Audit and reproduction

- [All candidate decisions](candidate-decisions.csv)
- [Every RNA measurement, chosen ceiling and local IHC support](rna-threshold-audit.csv.gz)
- [Recovered candidates](recovered-candidates.csv)
- [All 42 selected targets](notable-targets.csv)
- [Overall counts](scope-summary.csv) and [publication counts](source-counts.csv)
- [Exact policy](screening-policy.json) and [checksums / routine contexts](run-manifest.json)
- [Pinned original observations](../cta-permitted-tissues-20261019/normal-observations.csv.gz)

```sh
python scripts/cta_adaptive_report.py --out outputs/cta-adaptive-current
```

The rebuild verifies input hashes, reproduces every fixed expression and final decision (629 retained), and requires the fixed set to be a subset of the adaptive set. Only the absolute RNA ceilings change; input measurements and other gates are held constant.

## References

{refs}
"""
    (out / "index.md").write_text(text)
    pdf(out, ladder, counts, target_rows)


def pdf(out, ladder, counts, target_rows):
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
        ParagraphStyle(name="Small", fontName="Helvetica", fontSize=8.5, leading=11, spaceAfter=6)
    )
    styles.add(ParagraphStyle(name="Cell", fontName="Helvetica", fontSize=7.5, leading=9.5))

    def p(text, style="Small"):
        return Paragraph(escape(str(text)), styles[style])

    def table(rows, widths):
        obj = Table(
            [[p(cell, "Cell") for cell in row] for row in rows],
            colWidths=widths,
            repeatRows=1,
            hAlign="LEFT",
        )
        obj.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8edf2")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#9ca3af")),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    (
                        "ROWBACKGROUNDS",
                        (0, 1),
                        (-1, -1),
                        [colors.white, colors.HexColor("#f6f8fa")],
                    ),
                ]
            )
        )
        return obj

    story = [
        p("CTA screen: IHC-adaptive RNA ceilings", "Title"),
        p("Human HPA v23 | approved project policy | 2026-10-20", "Heading2"),
        table(ladder, [260, 130, 130]),
        Spacer(1, 10),
        p(METHOD),
        p(GATES),
        p(LIMITS),
        p("Results and coverage", "Heading2"),
        table(counts, [260, 65, 65, 65, 65]),
        PageBreak(),
        p("Policy comparison", "Heading2"),
        Image(str(out / "policy-comparison.png"), width=520, height=208),
        p("Published source retention", "Heading2"),
        Image(str(out / "source-retention.png"), width=520, height=364),
        p(
            "Sources overlap and cannot be summed. The complete audit includes all 2,537 published candidates."
        ),
        PageBreak(),
        p("Known and published targets", "Heading2"),
        Image(str(out / "published-targets.png"), width=460, height=632.5),
        PageBreak(),
        p("Decisions and supporting evidence", "Heading2"),
        p(
            "Selected clinical and research targets. The CSV preserves all failing and rescued contexts; summaries below are abbreviated."
        ),
        table(target_rows, [105, 35, 45, 335]),
        PageBreak(),
        p("Source interpretation: CTAG2 and PRAME", "Heading2"),
        p(SOURCE_NOTE),
        p(
            "RNA and IHC donors differ. Neither a low central RNA estimate nor negative IHC establishes absence of protein or peptide-HLA. Existing core/extended API membership is separate from these comparative screens."
        ),
        p("References", "Heading2"),
    ]
    for name, url in REFERENCES.items():
        story.append(
            Paragraph(
                f'<link href="{escape(url)}" color="#205c8f">{escape(name)}</link>', styles["Small"]
            )
        )
    story.append(
        Paragraph(
            '<link href="https://pirl-unc.github.io/oncoref/audits/cta-ihc-adaptive-20261020/" color="#205c8f">Full web report, source citations and downloadable audit</link>',
            styles["Small"],
        )
    )

    def footer(canvas, doc):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#64748b"))
        canvas.drawString(44, 24, "OncoRef | 2026-10-20 | Human HPA v23 | Project screening policy")
        canvas.drawRightString(568, 24, str(doc.page))

    SimpleDocTemplate(
        str(out / "cta-ihc-adaptive-report.pdf"),
        pagesize=letter,
        leftMargin=44,
        rightMargin=44,
        topMargin=38,
        bottomMargin=40,
        title="CTA IHC-adaptive absolute RNA screening report",
    ).build(story, onFirstPage=footer, onLaterPages=footer)
