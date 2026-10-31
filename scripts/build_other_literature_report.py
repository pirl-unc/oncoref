#!/usr/bin/env python3
"""Rebuild a complete per-protein Other audit and readable PDF from pinned evidence."""

from __future__ import annotations

import hashlib
import json
from html import escape
from pathlib import Path

import pandas as pd

from oncoref import cta_selection_policy, gene_knowledge

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-other-literature-20261030"
CURRENT = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
DESCRIPTIONS = {
    "Low": "Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target.",
    "Weak": "No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps.",
    "Normal": "Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify.",
    "Disposable": "Appreciable expression is confined to permitted Disposable contexts under the measured screen; no qualifying CTA identity.",
    "Diffuse": "The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe.",
    "Tumor-supported": "Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target.",
    "Reproductive": "Reproductive evidence qualifies but cancer presence is unestablished.",
    "Inducible": "Normal activation evidence requires safety adjudication.",
}


def markdown(frame):
    def cell(value):
        return str(value).replace("|", " / ").replace("\n", " ")

    return "\n".join(
        [
            "| " + " | ".join(frame.columns) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |",
            *[
                "| " + " | ".join(map(cell, row)) + " |"
                for row in frame.itertuples(index=False, name=None)
            ],
        ]
    )


def write_csv(frame, path):
    frame.to_csv(
        path,
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


def build():
    baseline = pd.read_csv(OUT / "baseline-proteins.csv.gz", low_memory=False)
    current = pd.read_csv(CURRENT, low_memory=False).set_index("protein_equivalence_group_id")
    knowledge = gene_knowledge.findings()
    genes = pd.read_csv(OUT / "genes.csv")
    searches = pd.read_csv(OUT / "searches.csv")
    annotations = pd.read_csv(OUT / "uniprot-annotations.csv").fillna("")
    observations = pd.read_csv(OUT / "normal-contexts.csv.gz")
    rows = []
    for old in baseline.sort_values("protein_group_label").itertuples():
        new = current.loc[old.protein_equivalence_group_id]
        member_ids = set(json.loads(old.member_gene_ids))
        findings = knowledge[
            knowledge.gene_ids.map(
                lambda values, member_ids=member_ids: bool(member_ids & set(values))
            )
        ]
        decisive = findings[findings.effect.isin(["exclude", "review", "reproductive"])]
        lookup = genes[genes.protein_id.eq(old.protein_equivalence_group_id)]
        receipts = searches[searches.gene_id.isin(member_ids)]
        obs = observations[observations.gene_id.eq(old.representative_gene_id)]
        supported = obs[
            obs.supported_at_least_floor & obs.category.isin(["vital", "protected", "manageable"])
        ]
        normals = "; ".join(sorted({f"{r.tissue} / {r.cell_type}" for r in supported.itertuples()}))
        if pd.isna(new.category):
            reason = "Normal evidence is unresolved; review before direct targeting or shared-sequence clearance."
        elif new.category == "Unsafe":
            reason = "Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury."
        elif new.category == "CTA":
            reason = "Independent cancer presence, reproductive presence, concentration, normal safety and required coverage pass."
        else:
            reason = DESCRIPTIONS[new.subtype]
        if not decisive.empty:
            reason += " " + " ".join(dict.fromkeys(decisive.finding))
        if new.category == "CTA" and old.protein_group_label in ["CDY1/B", "OR2T7"]:
            reason += " Positive cancer RNA below 2 TPM now establishes presence; it does not establish therapeutic abundance."
        local_ann = annotations[annotations.gene_id.isin(member_ids)]
        rows.append(
            {
                "protein_id": old.protein_equivalence_group_id,
                "protein": old.protein_group_label,
                "genes": "; ".join(lookup.gene_symbol),
                "gene_ids": json.dumps(sorted(member_ids)),
                "baseline": old.classification,
                "classification": new.classification,
                "direct_target_candidate": new.category == "CTA"
                if pd.notna(new.category)
                else False,
                "shared_sequence": "review"
                if pd.isna(new.category)
                else "exclude"
                if new.category == "Unsafe"
                else "allow",
                "decision": reason,
                "review_scope": "Individual accepted finding / explicit unresolved finding"
                if len(decisive)
                else "Measured-profile review with search/annotation triage; no accepted individual override",
                "findings": json.dumps(findings.id.tolist()),
                "source_links": "; ".join(findings.source_url.drop_duplicates()),
                "searches": len(receipts),
                "search_complete": bool(receipts.complete.all()),
                "annotation_available": bool(len(local_ann)),
                "supported_normal_contexts": normals,
                "maximum_normal_ntpm": new.healthy_maximum_ntpm,
                "maximum_normal_tpm": new.healthy_maximum_tpm,
                "maximum_normal_cp10k": new.healthy_maximum_cp10k,
                "reproductive_ntpm": new.reproductive_ntpm,
                "tumor_tpm": new.oncoref_maximum_tpm,
                "hpa_tumor_ptpm": new.hpa_tumor_group_lower_bound_ptpm,
            }
        )
    decisions = pd.DataFrame(rows)
    assert len(decisions) == 600 and decisions.protein_id.is_unique
    write_csv(decisions, OUT / "protein-decisions.csv")
    export_findings = knowledge.copy()
    for column in ["gene_ids", "gene_symbols"]:
        export_findings[column] = export_findings[column].map(json.dumps)
    write_csv(export_findings, OUT / "reviewed-findings.csv")
    counts = current.classification.value_counts().to_dict()
    transitions = (
        decisions.groupby(["baseline", "classification"]).size().rename("proteins").reset_index()
    )
    write_csv(transitions, OUT / "transitions.csv")
    changes = decisions[decisions.baseline.ne(decisions.classification)]
    write_csv(changes, OUT / "changed-proteins.csv")
    # The >30 raw observations are retained in a review ledger, never hidden as zero.
    high = observations[
        observations.category.isin(["vital", "protected", "manageable"])
        & observations.unit.isin(["TPM", "nTPM"])
        & observations.raw_value.gt(30)
    ].copy()
    high["classification"] = high.protein_label.map(decisions.set_index("protein").classification)
    if high.classification.str.startswith("Other:").any():
        raise ValueError("Unresolved >30 normal observations remain in Other")
    write_csv(high, OUT / "high-raw-normal-review.csv")
    sections = [
        "# Other proteins: literature and category audit",
        "",
        "This audit reviews the **600 Other proteins / 611 member genes in 1.25.0**. Every member was searched. Search retrieval and annotation triage are not equivalent to reading and adjudicating every paper. The per-protein table states the depth of review; absence of a search hit never proves normal absence.",
        "",
        "[Full decisions](protein-decisions.csv) · [Changed proteins](changed-proteins.csv) · [Reviewed findings and sources](reviewed-findings.csv) · [PDF report](other-protein-review.pdf) · [Current complete CTA report](../cta-source-walkthrough-20261028/index.md)",
        "",
        "## Current complete universe",
        "",
        markdown(pd.DataFrame(counts.items(), columns=["Classification", "Proteins"])),
        "",
        "## What changed among the 600",
        "",
        markdown(transitions),
        "",
        "## Rules and consequences",
        "",
        "1. **Normal safety comes first.** Existing tissue/IHC RNA ceilings remain unchanged. Reviewed adult human normal protein, or localized RNA plus orthogonal normal-cell function, can establish exclusion in Vital/Protected. Rare cells and normal activation states count. Literature findings retain assay, species, stage, exact identity and limitations.",
        "2. **Unresolved normal evidence stays out of Other.** Sparse sensory-cell RNA, unresolved species/isoform attribution, and high raw normal observations without sufficient support receive a review hold. This does not label them clinically toxic.",
        "3. **Cancer presence means >0 measured TPM or HPA pTPM in one compatible sample, or reliable tumor protein.** Abundance is reported separately. No cancer/normal ratio is required. Microarray proxies, pseudobulk and invalid legacy scales cannot be compared numerically with TPM.",
        "4. **Reproductive presence accepts RNA or reliable protein.** The RNA-only floor remains 5 nTPM, with the existing read/cell support rules. Reviewed reproductive protein is an independent route and must cover the selected identical-sequence members. Concentration, safety and coverage still apply.",
        "5. **CTA / Unsafe / Other / Unassessed have distinct sequence roles.** CTA is a candidate direct target. Unsafe supplies an exclusion background. Other supplies no exclusion under this policy but is not a direct-target category. Unassessed requires review and cannot silently be treated as allowed. These are research roles, not epitope/HLA or all-isoform clearance.",
        "",
        "## Not all Other:Normal proteins are skin proteins",
        "",
        "Among the 170 baseline Normal proteins with supported nonreproductive RNA at or above the descriptive floor, **41 have skin measurements; 33 reach the floor only in skin**. Overlapping counts include **81 Vital**, **59 Protected**, **50 Manageable**, and **4 Disposable**. These counts overlap and must not be added. Most observations remain below the existing safety ceiling. The full table retains cell names and native scales.",
        "",
        "[Every qualifying tissue/cell measurement](normal-category-tissues.csv) · [Tissue summary](normal-tissue-summary.csv) · [Raw and supported observations](normal-contexts.csv.gz)",
        "",
        "## Reconsideration and limitations",
        "",
        "- The human olfactory supplement identifies 68 baseline Other receptor genes in 1–20 positive neurons. Sparse detection is a review hold, not proof of protein abundance. Unlisted receptor genes are not thereby negative.",
        "- Species, stage and isoforms matter: rat kidney IHC cannot become human IHC; fetal cochlear cloning does not prove adult protein; DUX4 full-length and shorter healthy-muscle isoforms need sequence-specific treatment.",
        "- Three raw normal observations exceed 30 TPM/nTPM despite failing recurrence/read support. They are retained as review holds with the raw values visible, not relabeled zero.",
        "- Published normal production is not a demonstrated therapy-related adverse event. Functional receptor inhibition also differs from cytotoxic targeting of the expressing cell.",
        "- The baseline universe is representative selected protein sequences, not every splice isoform. Normal atlases still miss some rare cells, activation states and donor variability.",
        "",
        "## Cancer measurements for CDY1/B and OR2T7",
        "",
        "[All 143 registered cohort views, with native-scale and missingness labels](requested-cancer-sources.csv). Compatible source maxima remain approximately 1.597 TPM for CDY1/B and 1.184 TPM for OR2T7. HPA independently reports 1.208 pTPM for CDY1/B. Positive low values now attest to cancer presence; they are not claims of high tumor abundance.",
        "",
        "## Reusable API and pinned sources",
        "",
        '```python\nfrom oncoref import cta_selection, gene_knowledge, reference_data\n\nfindings = gene_knowledge.findings("GFRAL")\nroles = cta_selection.sequence_dispositions()\nsource = reference_data.ensure("olfactory_durante_supplement", "2020-02-18")\nreceipt = reference_data.provenance("olfactory_durante_supplement", "2020-02-18")\n```',
        "",
        "The knowledge reference ships in the wheel. The olfactory workbook uses the same versioned, checksummed reference-data cache as HPA/DICE/cornea. Native counts/fractions stay in their own units. Human findings never transfer to dog or mouse by orthology. Legacy exact-protein backgrounds remain separately versioned; downstream code must explicitly adopt these roles.",
        "",
        "## Search accounting",
        "",
        f"{len(searches):,} search queries; {searches.gene_id.nunique()} genes; {int(searches.complete.sum()):,} fully paginated queries. Four broad IFNA1/NPS queries reached the retrieval cap; aliases also create unrelated hits. The bibliography is a lead list, not accepted evidence. Accepted findings are listed separately with their limitations.",
        "",
        "## Every baseline Other protein",
        "",
        markdown(
            decisions[["protein", "baseline", "classification", "decision"]].rename(
                columns={
                    "protein": "Protein",
                    "baseline": "Before",
                    "classification": "Now",
                    "decision": "Reason",
                }
            )
        ),
        "",
        "## Reviewed source findings",
        "",
    ]
    for r in knowledge.itertuples():
        sections += [
            f"### {', '.join(r.gene_symbols)} — {r.effect}",
            "",
            r.finding,
            "",
            f"**Context:** {r.tissue}; {r.cell_type}; {r.state}. **Method:** {r.assay}.",
            "",
            f"**Consequence:** {r.decision_reason}",
            "",
            f"**Limitations:** {r.limitations}",
            "",
            f"[{r.source_title}]({r.source_url}) — {r.source_location}.",
            "",
        ]
    (OUT / "index.md").write_text("\n".join(sections))
    make_pdf(decisions, transitions, counts)
    manifest = json.loads((OUT / "search-manifest.json").read_text())
    manifest["files"] = {
        name: hashlib.sha256((OUT / name).read_bytes()).hexdigest() for name in manifest["files"]
    }
    (OUT / "search-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    receipt = {
        "policy_id": cta_selection_policy.POLICY_ID,
        "baseline_proteins": 600,
        "member_genes": 611,
        "findings": len(knowledge),
        "inputs": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [CURRENT, ROOT / "oncoref/data/gene-knowledge.json", Path(__file__)]
        },
        "outputs": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(OUT.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(transitions.to_string(index=False))


def make_pdf(decisions, transitions, counts):
    from reportlab.lib import colors
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.platypus import (
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", fontName="Helvetica", fontSize=8, leading=10))

    def p(s, style="BodyText"):
        return Paragraph(escape(str(s)).replace("\n", "<br/>"), styles[style])

    def table(frame, widths):
        cells = [[p(c, "Small") for c in frame.columns]] + [
            [p(v, "Small") for v in row]
            for row in frame.fillna("").itertuples(index=False, name=None)
        ]
        t = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E5EFF2")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 6),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                    ("LINEBELOW", (0, 0), (-1, -1), 0.25, colors.HexColor("#D1DEE4")),
                ]
            )
        )
        return t

    story = [
        p("OncoRef: Other protein review", "Title"),
        p("600 baseline proteins / 611 member genes", "Heading2"),
        p(
            "This report separates confirmed normal production, unresolved normal findings, and proteins that remain Other under the measured screen. Every member was searched; retrieved papers are leads until individually adjudicated. Research categories are not clinical safety guarantees."
        ),
        Spacer(1, 14),
        p("Complete current universe", "Heading2"),
        table(pd.DataFrame(counts.items(), columns=["Classification", "Proteins"]), [390, 110]),
        Spacer(1, 14),
        p("Independent gates", "Heading2"),
        p(
            "Normal safety establishes Unsafe. Cancer presence (>0 TPM/pTPM or reliable tumor protein), reproductive presence (supported RNA or protein), concentration, and required coverage independently establish CTA eligibility. Missing or conflicting normal evidence is Unassessed."
        ),
        p("Sequence roles", "Heading2"),
        p(
            "CTA: candidate direct target. Unsafe: exclude shared sequence. Other: no exclusion under this research policy, but not a direct target. Unassessed: review; never silently allow. Selected reference proteins do not cover every isoform or establish peptide presentation."
        ),
        PageBreak(),
        p("Transitions among the 600", "Title"),
        table(
            transitions.rename(
                columns={"baseline": "Before", "classification": "After", "proteins": "Proteins"}
            ),
            [190, 240, 70],
        ),
        Spacer(1, 12),
        p("Tissue audit", "Heading2"),
        p(
            "Of 170 baseline Other:Normal proteins with appreciable nonreproductive RNA, 41 have skin measurements and 33 reach the floor only in skin. Overlapping counts include 81 Vital and 59 Protected. Most remain below the existing safety ceilings. Rare normal cells and activation states need separate evidence."
        ),
        p("Evidence limitations", "Heading2"),
        p(
            "Sparse olfactory-neuron counts are review evidence, not TPM or protein. Human/animal, adult/fetal and gene/isoform evidence are not interchangeable. Three raw normal observations above 30 remain explicit review holds despite inadequate read or recurrence support."
        ),
        p("Linked complete evidence", "Heading2"),
        Paragraph(
            '<link href="https://pirl-unc.github.io/oncoref/audits/cta-other-literature-20261030/">Online audit, all source citations, per-protein CSV, native tissue measurements, and 143-cohort cancer audit</link>',
            styles["BodyText"],
        ),
        PageBreak(),
        p("Every baseline Other protein", "Title"),
        p(
            "The table includes every member of the baseline protein set. Full native measurements, reviewed sources and search coverage are provided in the linked CSV and online reference."
        ),
        Spacer(1, 12),
        table(
            decisions[["protein", "classification", "decision"]].rename(
                columns={"protein": "Protein", "classification": "Now", "decision": "Reason"}
            ),
            [85, 115, 300],
        ),
    ]
    dest = ROOT / "output/pdf/other-protein-review.pdf"
    dest.parent.mkdir(parents=True, exist_ok=True)

    def footer(c, doc):
        c.setFont("Helvetica", 8)
        c.drawString(36, 20, "OncoRef | " + cta_selection_policy.POLICY_ID)
        c.drawRightString(559, 20, str(doc.page))

    SimpleDocTemplate(
        str(dest),
        pagesize=(595, 842),
        rightMargin=45,
        leftMargin=45,
        topMargin=40,
        bottomMargin=36,
        title="OncoRef Other protein review",
        author="OncoRef",
        invariant=True,
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    (OUT / "other-protein-review.pdf").write_bytes(dest.read_bytes())


if __name__ == "__main__":
    build()
