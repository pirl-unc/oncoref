#!/usr/bin/env python3
"""Build the linked TCR-T SAE review, tables and PDF from the curated ledger."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import cta_selection, tcr_trials  # noqa: E402
from oncoref.cta_selection_policy import POLICY_ID  # noqa: E402
from oncoref.version import __version__  # noqa: E402

OUTPUT = ROOT / "docs/audits/cta-tcr-trials-20261029"
REPLAY = ROOT / "docs/audits/cta-selection-20261028"


def summaries() -> pd.DataFrame:
    """One product per row, with all tested/cross-recognized protein origins."""
    ledger = tcr_trials.evidence()
    decisions = cta_selection.results()
    by_gene = {g: r for r in decisions.to_dict("records") for g in r["gene_ids"]}
    rows = []
    for p in tcr_trials.products().to_dict("records"):
        cohorts = [c for c in ledger["cohorts"] if c["product_id"] == p["product_id"]]
        events = [e for e in ledger["events"] if e["product_id"] == p["product_id"]]
        proteins = []
        for t in p["targets"]:
            r = by_gene.get(t["gene_id"], {})
            category = r.get("category")
            proteins.append(f"{t['symbol']}: {category if pd.notna(category) else 'Unassigned'}")
        rows.append(
            {
                "product_id": p["product_id"],
                "Product / receptor": p["name"],
                "Proteins and current categories": "; ".join(proteins),
                "Peptides": "; ".join(sorted({t["peptide"] for t in p["targets"]})),
                "HLA": "; ".join(p["hla"]),
                "Treated cohorts": "; ".join(str(c["patients_treated"]) for c in cohorts),
                "Serious / severe events": " ".join(e["event"] + "." for e in events)
                or "No sentinel event entered; see cohort safety summary.",
                "Clinical interpretation": " ".join(c["finding"] for c in cohorts),
                "RNA credit": "One step, capped" if p["rna_credit"] else "None",
                "Sources": "; ".join(
                    dict.fromkeys(ledger["sources"][c["source"]]["url"] for c in cohorts)
                ),
            }
        )
    return pd.DataFrame(rows)


def summary_markdown() -> str:
    """Compact section embedded directly in the main curation report."""
    d = summaries()
    lines = [
        "## Serious adverse events from CTA-directed TCR-T",
        "",
        "**Treatment toxicity and antigen-specific injury are different findings.** "
        "SAE seriousness (for example hospitalization or death) is distinct from CTCAE severity. "
        "This targeted review includes important severe events and ungraded serious cases; "
        "it is not an exhaustive SAE database or a pooled incidence estimate.",
        "",
        "| Protein origins (current screen) | Receptor / product | Important events | RNA credit |",
        "|---|---|---|---|",
    ]
    for row in d.to_dict("records"):
        lines.append(
            "| "
            + " | ".join(
                [
                    row["Proteins and current categories"],
                    f"[{row['Product / receptor']}]({row['Sources'].split('; ')[0]})",
                    row["Serious / severe events"],
                    row["RNA credit"],
                ]
            )
            + " |"
        )
    lines += [
        "",
        "The MAGE-A3 **a3a** cardiac cases implicated **titin (TTN)**, an unintended target; "
        "the separate NCI MAGE-A3/A9 receptor neurological cases implicated the cross-recognized "
        "**MAGE-A12** peptide. These findings do not establish that MAGE-A3 itself was present in heart or brain. "
        "TBI-1301's grade 3 lung injury retains an unresolved-specificity hold for that product. "
        "c259's late cardiac death and Guillain-Barre cases remain caveats without an established peptide-specific mechanism.",
        "",
        "[Full attribution, epitope/HLA mappings, cohort denominators and limitations]"
        "(../cta-tcr-trials-20261029/index.md). "
        "A favorable reviewed product supplies one capped RNA step, never protein or clinical-injury clearance. "
        "The exact tested peptide matters: NY-ESO-1 and **CTAG2/LAGE-1A** share SLLMWITQC; "
        "this does not clear LAGE-1B or other CTAG2 peptides.",
    ]
    return "\n".join(lines)


def build() -> None:
    """Rebuild all formats from the same ledger and released selection snapshot."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ledger = tcr_trials.evidence()
    summary = summaries()
    baseline_meta = json.loads((OUTPUT / "baseline-metadata.json").read_text())
    baseline_file = OUTPUT / "baseline-human.csv.gz"
    if sha256(baseline_file) != baseline_meta["data_sha256"]["cta-selection-human.csv.gz"]:
        raise ValueError("Changed baseline selection snapshot")
    before = pd.read_csv(baseline_file, keep_default_na=False).set_index("entity_id")
    after = cta_selection.results().set_index("entity_id")
    if set(before.index) != set(after.index):
        raise ValueError("Clinical-credit comparison requires the same protein universe")
    comparison = pd.DataFrame(
        {"label": after.label, "before": before.category, "after": after.category.fillna("")}
    )
    changes = comparison[comparison.before.ne(comparison.after)].reset_index()
    write_csv(changes, OUTPUT / "membership-changes.csv")
    write_csv(summary, OUTPUT / "tcr-trial-summary.csv")
    write_csv(tcr_trials.cohorts(), OUTPUT / "cohorts.csv")
    write_csv(tcr_trials.events(), OUTPUT / "sentinel-events.csv")
    (OUTPUT / "evidence.json").write_text(json.dumps(ledger, indent=2) + "\n")
    contexts = pd.read_csv(REPLAY / "widely-studied-contexts.csv.gz", low_memory=False)
    selected = contexts[
        contexts.clinical_rna_credit & contexts.value.notna() & contexts.candidate_ceiling.notna()
    ]
    write_csv(selected, OUTPUT / "clinical-credit-contexts.csv.gz")
    # The complete replay records every context for each clinically reviewed protein.
    grouped = (
        selected.groupby(["gene_id", "assay"], sort=True)
        .agg(
            measured_contexts=("value", "size"),
            raised_ceilings=("clinical_credit_applied", "sum"),
            baseline_exclusions=("baseline_exclusion", "sum"),
            adjusted_exclusions=("candidate_exclusion", "sum"),
        )
        .reset_index()
    )
    write_csv(grouped, OUTPUT / "clinical-credit-impact.csv")
    intro = """# TCR-T clinical evidence, serious events and RNA screening credit

This is a targeted, source-linked review of published CTA-directed TCR-T experience.
**A favorable trial is evidence for a particular peptide, HLA, receptor and regimen.**
The approved policy grants a modest protein-profile RNA credit, with those limitations
attached; it does not declare the whole protein clinically safe.

## How the evidence changes the screen

1. Require reviewed results from an actually treated cohort, with its denominator.
2. Identify the peptide, its protein origins, the HLA and the receptor/product.
3. Review whether any grade 4-5 event was suspected or confirmed to involve recognition
   of healthy tissue. Such a signal blocks credit for that receptor. Unassessed
   specificity and ungraded suspected specificity injury also block credit.
4. Otherwise grant **one RNA ladder step**, with the Enhanced-negative ceiling as the cap.
   Trials do not stack. The IHC grade is not changed or invented.
5. Keep independent protein and clinical-injury exclusions, required coverage,
   cancer expression, Reproductive evidence and diffuseness rules.

| Tier | Baseline: missing / Approved / Supported / Enhanced IHC | With qualifying clinical credit |
|---|---|---|
| Vital | 4 / 8 / 16 / 32 | 8 / 16 / 32 / 32 |
| Protected | 8 / 16 / 32 / 64 | 16 / 32 / 64 / 64 |
| Manageable | 16 / 32 / 64 / 128 | 32 / 64 / 128 / 128 |

Values are nominal TPM/nTPM screening parameters, not clinically measured toxicity limits.
Corneal CP10K retains its separate nominal /100 scale. An unfavorable product does not
automatically make every peptide of its intended protein Unsafe. For example, titin
cross-recognition implicates a receptor and a different protein. Accepted gene-level
normal-protein/clinical counterevidence is still handled by the existing safety layer.

CTAG2 was already CTA before this credit. Its HPA v23 heart value is 5.1 nTPM with
Enhanced-negative IHC and a 32 nTPM ceiling; the cap leaves that ceiling at 32.
The shared clinical peptide is in **LAGE-1A**, not every CTAG2 isoform. Reference
isoform IDs, peptide positions and sequence checksums are in the evidence ledger.
Gene-level RNA cannot identify which isoform supplied a normal-tissue measurement.

## Reading the safety findings

“No identified grade 4-5 specificity signal” is a reviewed attribution conclusion,
not an author-reported zero incidence or proof that risk is zero. Small cohorts,
HLA selection, follow-up, conditioning and incomplete mechanistic investigations
limit the inference. This initial review is not systematic or exhaustive. Unlisted
products remain unreviewed; registration, preclinical work and missing reporting do
not qualify. No minimum response rate is imposed on the independent safety credit.

The sentinel table includes deaths, significant organ injury and severe inflammatory
events. It is not a substitute for each paper's complete AE/SAE tables. Unknown CTCAE
grades remain null. Case reports do not supply new treated-cohort denominators.
"""
    sections = [
        intro,
        summary_markdown().replace("../cta-tcr-trials-20261029/index.md", "index.md"),
    ]
    sections.append(
        f"\n## Measured impact\n\n**{len(changes)} protein category changes** versus "
        f"the pinned {baseline_meta['species']['human']['policy_id']} snapshot. "
        "[Membership comparison](membership-changes.csv). RNA-ceiling changes and "
        "individual exclusion effects are in the linked context and assay tables; "
        "a changed ceiling need not change the final category.\n"
    )
    for p in tcr_trials.products().to_dict("records"):
        sections += [
            f"\n## {p['name']}\n",
            f"**RNA credit: {'one capped step' if p['rna_credit'] else 'none'}.** {p['review_scope']}\n",
            "HLA: " + ", ".join(p["hla"]) + ".\n",
        ]
        for t in p["targets"]:
            sections.append(
                f"- **{t['symbol']}** (`{t['gene_id']}`): `{t['peptide']}`. {t['isoform_scope']}.\n"
            )
        for c in [c for c in ledger["cohorts"] if c["product_id"] == p["product_id"]]:
            source = ledger["sources"][c["source"]]
            sections.append(
                f"\n**{c['patients_treated']} treated**, {c['population']}. {c['finding']} "
                f"{c['qualification']} {c['follow_up']} "
                f"[{source['title']}]({source['url']}).\n"
            )
        for e in [e for e in ledger["events"] if e["product_id"] == p["product_id"]]:
            source = ledger["sources"][e["source"]]
            sections.append(
                f"\n**{e['event']}** ({e['patients']} patient(s); grade {e['grade'] or 'unreported'}). "
                f"Treatment attribution: {e['treatment_attribution'].replace('_', ' ')}. "
                f"Specificity: {e['specificity_attribution'].replace('_', ' ')}. "
                f"{e['details']} [{source['title']}]({source['url']}).\n"
            )
    sections += [
        """
## API, source storage and reproduction

```python
from oncoref import tcr_trials, reference_data

tcr_trials.products()       # product-level eligibility, exact peptide origins and HLA
tcr_trials.cohorts()        # treated denominators and explicit attribution review
tcr_trials.events()         # sentinel events, nullable severity and attribution
tcr_trials.protein_credit(["ENSG00000126890"])  # CTAG2; not isoform-wide clearance
reference_data.ensure("tcr_trial_lete_mrcls_2025", "snapshot-2026-10-29")
```

`oncoref tcr-trials --table events` exports JSON. `--table all` includes the complete
ledger. Original papers/protocols use the same checksum-pinned source API and
`<cache>/sources/<source>/<version>/<filename>` layout as other references;
`CANCERDATA_DATA_DIR` selects the cache root. Paper URLs can change upstream:
checksum mismatch fails explicitly instead of silently updating the reviewed snapshot.

Reproduce with `python scripts/build_tcr_trial_report.py` after the selection replay.
The curation lives in `oncoref/data/tcr-trial-evidence.json`; both the ledger and policy
code participate in the replay receipts. The report does not re-adjudicate literature
on each build. New evidence requires an explicit reviewed ledger update.

[Summary table](tcr-trial-summary.csv) · [Cohorts](cohorts.csv) ·
[Sentinel events](sentinel-events.csv) · [Complete evidence](evidence.json) ·
[Context-level RNA effects](clinical-credit-contexts.csv.gz) · [RNA effects by assay](clinical-credit-impact.csv)
"""
    ]
    name = f"cta-tcr-safety-{__version__}.pdf"
    sections.append(f"\n[PDF report]({name})\n")
    (OUTPUT / "index.md").write_text("\n".join(sections))
    build_pdf(summary, ledger, name)
    (OUTPUT / "manifest.json").write_text(
        json.dumps(
            {
                "policy_id": POLICY_ID,
                "inputs": input_hashes(
                    [
                        Path(__file__),
                        ROOT / "oncoref/tcr_trials.py",
                        ROOT / "oncoref/cta_selection_policy.py",
                        ROOT / "oncoref/data/tcr-trial-evidence.json",
                        ROOT / "oncoref/data/reference-releases.json",
                        ROOT / "oncoref/data/cta-selection.json",
                        ROOT / "scripts/verify_tcr_trial_sequences.py",
                        OUTPUT / "baseline-metadata.json",
                        baseline_file,
                        REPLAY / "widely-studied-contexts.csv.gz",
                    ]
                ),
                "outputs": {
                    p.name: sha256(p)
                    for p in sorted(OUTPUT.iterdir())
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )


def build_pdf(summary, ledger, name):
    """Use short product pages; keep event attribution and source links readable."""
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import (
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    styles = getSampleStyleSheet()
    for key in ["BodyText", "Normal"]:
        styles[key].fontName = "Helvetica"
        styles[key].fontSize = 10
        styles[key].leading = 14
    styles["BodyText"].spaceAfter = 8
    styles["Title"].alignment = TA_LEFT
    styles["Title"].textColor = colors.HexColor("#173B55")
    styles["Heading2"].textColor = colors.HexColor("#087E8B")
    story = []

    def para(text, style="BodyText"):
        return Paragraph(text, styles[style])

    story += [
        para("CTA-directed TCR-T: serious events and clinical evidence", "Title"),
        para(f"OncoRef {__version__} | Reviewed 29 October 2026 | {POLICY_ID}"),
        para(
            "Clinical experience can modestly relax the RNA screen. It cannot clear every peptide, isoform or receptor from the parent protein.",
            "Heading2",
        ),
        para(
            "Grant one existing RNA ladder step for a reviewed treated TCR-T cohort without an identified grade 4-5 specificity-related signal. Cap at Enhanced-negative. Preserve normal-protein and clinical-injury exclusions. Trial registration and unknown reporting do not qualify."
        ),
        para(
            "SAE seriousness is different from CTCAE severity. The following pages include fatal, serious and important severe events; this is a targeted review, not an exhaustive SAE list or pooled incidence estimate."
        ),
        para("Clinical-adjusted RNA ceilings", "Heading2"),
    ]
    table = Table(
        [
            [para(x) for x in r]
            for r in [
                ["Tier", "Baseline IHC ladder", "With clinical credit"],
                ["Vital", "4 / 8 / 16 / 32", "8 / 16 / 32 / 32"],
                ["Protected", "8 / 16 / 32 / 64", "16 / 32 / 64 / 64"],
                ["Manageable", "16 / 32 / 64 / 128", "32 / 64 / 128 / 128"],
            ]
        ],
        colWidths=[105, 185, 210],
    )
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E7F0F2")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 9),
            ]
        )
    )
    story += [
        table,
        Spacer(1, 12),
        para(
            "Columns follow missing/unusable, Approved-negative, Supported-negative and Enhanced-negative IHC. Native TPM/nTPM retain nominal policy values; corneal CP10K uses its existing /100 policy. These are not measured clinical toxicity limits."
        ),
        para("CTAG2 / LAGE-1A", "Heading2"),
        para(
            "CTAG2 already passes the current screen. LAGE-1A shares SLLMWITQC with NY-ESO-1, recognized by c259 in lete-cel. HPA v23 heart RNA is 5.1 nTPM; Enhanced-negative IHC already gives the capped 32 nTPM ceiling. Other CTAG2 isoforms and peptides are not clinically cleared."
        ),
    ]
    for row in summary.to_dict("records"):
        pid = row["product_id"]
        story += [
            PageBreak(),
            para(escape(row["Product / receptor"]), "Title"),
            para(escape(row["Proteins and current categories"])),
            para("Peptide(s): " + escape(row["Peptides"]) + "<br/>" + escape(row["HLA"])),
            para("RNA credit: " + row["RNA credit"], "Heading2"),
        ]
        for c in [c for c in ledger["cohorts"] if c["product_id"] == pid]:
            s = ledger["sources"][c["source"]]
            story += [
                para(f"{c['patients_treated']} treated: " + escape(c["population"]), "Heading2"),
                para(escape(c["finding"])),
                para(escape(c["qualification"])),
                para(
                    f'<link href="{escape(s["url"])}" color="#087E8B">{escape(s["title"])}</link>'
                ),
            ]
        for e in [e for e in ledger["events"] if e["product_id"] == pid]:
            s = ledger["sources"][e["source"]]
            story += [
                para(escape(e["event"]), "Heading2"),
                para(
                    f"{e['patients']} patient(s); grade {e['grade'] or 'unreported'}. "
                    f"Treatment: {e['treatment_attribution'].replace('_', ' ')}. "
                    f"Specificity: {e['specificity_attribution'].replace('_', ' ')}."
                ),
                para(escape(e["details"])),
                para(
                    f'<link href="{escape(s["url"])}" color="#087E8B">Primary event source</link>'
                ),
            ]

    def footer(canvas, doc):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#536675"))
        canvas.drawString(42, 25, "OncoRef | Research screening, not clinical safety clearance")
        canvas.drawRightString(570, 25, str(doc.page))

    destination = ROOT / "output/pdf" / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(
        str(destination),
        pagesize=(612, 792),
        rightMargin=42,
        leftMargin=42,
        topMargin=42,
        bottomMargin=44,
        title="CTA TCR-T safety evidence",
        author="OncoRef",
        invariant=1,
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    (OUTPUT / name).write_bytes(destination.read_bytes())


if __name__ == "__main__":
    build()
