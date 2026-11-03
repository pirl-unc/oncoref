#!/usr/bin/env python3
"""Generate the central tissue-hazard ontology report, figures and inventories."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import tissue_hazards as hazards  # noqa: E402
from oncoref.normal_tissues import normal_tissue_mappings  # noqa: E402

COLORS = {
    "vital": "#963D4B",
    "protected": "#B27324",
    "manageable": "#307B89",
    "disposable": "#6A7197",
    "unresolved": "#80858C",
}
DEFAULT = ROOT / "docs/audits/tissue-hazards-1.29.0"
INTRO = """This ontology assigns adult human compartments to four protection priorities.
It separates immediate functional loss, persistent injury, regenerative capacity
and available compensation. The assignments are policy judgments supported by
different kinds of evidence, not measured probabilities of fatal toxicity.

The initiating target, inflammatory amplifiers and the tissue causing final organ
failure can be different. A preserved stem-cell pool does not ensure recovery
under persistent recognition. Conversely, the label 'stem cell' alone does not
prove irreplaceability: alternative progenitors sometimes compensate.

Reproductive/developmental membership remains an independent anatomical axis.
RNA detection, accessible surface antigen, peptide-HLA presentation, functional
recognition and demonstrated tissue injury are distinct evidence states. This
registry classifies the consequence of injury; it does not manufacture missing
recognition evidence. Adult human calls are not validated animal or pediatric
tolerability claims.
"""
POLICY = """Bulk measurements use organ-level rules. A narrower exception requires a
resolved tissue and cell identity from a cellular assay. Unknown cell labels
inherit organ concern; mixed labels cannot earn a downgrade. Unknown anatomy,
conflicting rules and contexts outside adult human scope remain unresolved.

The RNA and IHC cascade reads these same assignments. Existing numeric RNA/IHC
ceilings are research parameters, not newly calibrated by the literature review.
In the existing brain/spinal-cord/heart/lung rare-cell panel, qualifying Vital
findings retain the Unsafe screen, while Protected findings create a candidate
review hold. The 1% detection line is a screening parameter, not a risk estimate.
Other organs have ontology entries without an implied new rare-cell data panel.

Disposable is a narrow conditional allowance for specified compartments. A
removable organ does not make its widely shared cell types disposable elsewhere.
The existing reproductive permission remains separate from the functional hazard
reported here. Historical reports retain their original versioned policy.
"""


def write_csv(frame, path):
    d = frame.copy()
    for c in d:
        if d[c].map(lambda v: isinstance(v, (list, dict))).any():
            d[c] = d[c].map(
                lambda v: json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else v
            )
    d.to_csv(path, index=False, lineterminator="\n")


def figures(out, doc):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.hashsalt": "oncoref-hazards-1.29.0",
        }
    )
    paths = []

    def save(fig, name):
        fig.savefig(out / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
        fig.savefig(out / f"{name}.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        paths.append(name)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.set(xlim=(0, 11), ylim=(0, 7))
    ax.axis("off")
    ax.text(0, 6.7, "One ontology, separate questions", size=23, weight="bold", color="#172B3A")
    boxes = [
        (0, 4.6, "IDENTITY", "Organ + cell subset\nAdult human scope\nBulk stays broad"),
        (3.8, 4.6, "CONSEQUENCE", "Immediate function\nPersistent injury\nRepair and compensation"),
        (
            7.6,
            4.6,
            "EVIDENCE",
            "Clinical injury / recognition\nAnimal depletion / repair\nFunctional analogy / inference",
        ),
    ]
    for x, y, title, body in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y - 1.25),
                3.4,
                1.8,
                boxstyle="round,pad=0.12",
                facecolor="#EDF2F4",
                edgecolor="#D1DADD",
            )
        )
        ax.text(x + 0.15, y + 0.2, title, weight="bold", color="#172B3A")
        ax.text(x + 0.15, y - 0.2, body, va="top", linespacing=1.5, size=11)
    for i, t in enumerate(hazards.TIERS):
        x = i * 2.8
        ax.add_patch(
            FancyBboxPatch(
                (x, 1.7),
                2.5,
                0.65,
                boxstyle="round,pad=0.07",
                facecolor=COLORS[t],
                edgecolor="none",
            )
        )
        ax.text(
            x + 1.25,
            2.02,
            t.title(),
            color="white",
            ha="center",
            va="center",
            weight="bold",
            size=13,
        )
    ax.annotate(
        "",
        xy=(5.5, 2.6),
        xytext=(5.5, 3.2),
        arrowprops={"arrowstyle": "->", "color": "#172B3A", "lw": 2},
    )
    ax.text(
        5.5,
        0.85,
        "Protection priority is a policy judgment, not a fatality probability.",
        ha="center",
        size=13,
        weight="bold",
    )
    ax.text(
        5.5,
        0.25,
        "Reproductive membership and actual receptor recognition remain independent axes.",
        ha="center",
        size=11,
    )
    save(fig, "ontology-logic")

    a = pd.DataFrame(doc["assessments"])
    counts = pd.crosstab(a.system, a.tier).reindex(columns=hazards.TIERS, fill_value=0).sort_index()
    fig, ax = plt.subplots(figsize=(10, 7))
    left = np.zeros(len(counts))
    for tier in hazards.TIERS:
        ax.barh(
            counts.index,
            counts[tier],
            left=left,
            color=COLORS[tier],
            label=tier.title(),
            height=0.7,
        )
        left += counts[tier]
    ax.invert_yaxis()
    ax.set_xlabel("Number of curated rules (not cells, samples or toxicity events)")
    ax.set_title("Coverage across organ systems", loc="left", weight="bold", pad=44)
    ax.legend(ncols=4, loc="lower left", bbox_to_anchor=(0, 1.01), frameon=False)
    fig.tight_layout()
    save(fig, "organ-coverage")

    ids = [
        "lung",
        "alveolar_at2",
        "alveolar_at1",
        "airway_club",
        "airway_basal",
        "airway_ciliated",
        "airway_goblet",
        "airway_secretory",
        "airway_ionocyte",
    ]
    rows = [next(r for r in doc["assessments"] if r["id"] == i) for i in ids]
    designs = [
        ("cell_therapy_clinical", "Human cell therapy"),
        ("selective_ablation", "Animal depletion"),
        ("antigen_specific_t_cells", "Animal antigen-specific T cells"),
        ("functional_disease_cohort", "Human functional analogy"),
        ("functional_gene_deletion", "Gene deletion"),
        ("disease_tissue_and_mechanistic", "Human disease/mechanism"),
    ]

    def design_presence(record, design):
        cited = {doc["references"][k]["design"] for k in record["references"]}
        if design == "selective_ablation":
            present = any(
                d.startswith("selective_ablation") or d == "infection_and_ablation" for d in cited
            )
        else:
            present = design in cited
        if not present:
            return 0
        # Distinguish studies of the named lineage from references borrowed for
        # a neighboring compartment or shared defense function. In particular,
        # pneumocyte pathology cannot prove AT1/AT2 clinical causality, and
        # club-cell or mucin experiments cannot prove selective goblet killing.
        focused = {
            "lung": {d for d, _ in designs},
            "alveolar_at2": {"selective_ablation", "antigen_specific_t_cells"},
            "alveolar_at1": {"selective_ablation"},
            "airway_club": {"selective_ablation"},
            "airway_basal": {"selective_ablation", "disease_tissue_and_mechanistic"},
            "airway_ciliated": {"selective_ablation", "functional_disease_cohort"},
        }
        return 2 if design in focused.get(record["id"], set()) else 1

    matrix = np.array([[design_presence(r, d) for d, _ in designs] for r in rows])
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.imshow(
        matrix,
        cmap=matplotlib.colors.ListedColormap(["#F0F3F5", "#B9CCD8", "#526E80"]),
        vmin=0,
        vmax=2,
        aspect="auto",
    )
    ax.set_yticks(range(len(rows)), [r["label"] + "  /  " + r["tier"].title() for r in rows])
    ax.set_xticks(range(len(designs)), [s for _, s in designs], rotation=25, ha="right")
    for n, r in enumerate(rows):
        ax.get_yticklabels()[n].set_color(COLORS[r["tier"]])
    ax.set_title("Lung: evidence designs cited for each call", loc="left", weight="bold", pad=15)
    fig.text(
        0.5,
        0.01,
        "Dark = cited study focused on the named lineage/compartment. Pale = adjacent function or unresolved lineage.\nBlank = design not cited here. Neither color establishes tolerability or quantifies clinical risk.",
        ha="center",
        size=9,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, "lung-evidence")

    examples = [
        ("Lung", "lung", ["airway_ciliated", "airway_club", "alveolar_at2"]),
        ("Pancreas", "pancreas", ["beta", "pancreas_acinar"]),
        ("Blood / marrow", "blood", ["mature_b", "hsc", "marrow_stroma"]),
        ("Skin", "skin", ["keratinocyte", "skin_melanocyte", "epidermal_basal"]),
        ("Bowel", "large_bowel", ["enterocyte", "intestinal_crypt", "intestinal_paneth"]),
    ]
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set(xlim=(0, 11), ylim=(0, 6))
    ax.axis("off")
    ax.text(
        0,
        5.7,
        "Bulk tissue does not prove cell-restricted exposure",
        size=19,
        weight="bold",
        color="#172B3A",
    )
    for n, (label, bulk, subs) in enumerate(examples):
        y = 4.8 - n
        br = next(r for r in doc["assessments"] if r["id"] == bulk)
        ax.text(
            0,
            y,
            label + "\n" + br["tier"].title(),
            color=COLORS[br["tier"]],
            weight="bold",
            va="center",
        )
        for j, s in enumerate(subs):
            r = next(r for r in doc["assessments"] if r["id"] == s)
            x = 2.4 + j * 2.9
            ax.add_patch(
                FancyBboxPatch(
                    (x, y - 0.32),
                    2.6,
                    0.65,
                    boxstyle="round,pad=.08",
                    facecolor=COLORS[r["tier"]],
                    edgecolor="none",
                )
            )
            short = {
                "airway_ciliated": "Ciliated",
                "airway_club": "Club",
                "alveolar_at2": "AT2",
                "beta": "Beta cells",
                "pancreas_acinar": "Acinar",
                "mature_b": "Mature B cells",
                "hsc": "Hematopoietic progenitors",
                "marrow_stroma": "Marrow stroma",
                "keratinocyte": "Mature keratinocytes",
                "skin_melanocyte": "Cutaneous melanocytes",
                "epidermal_basal": "Epidermal basal",
                "enterocyte": "Mature enterocytes",
                "intestinal_crypt": "Crypt progenitors",
                "intestinal_paneth": "Paneth",
            }[s]
            ax.text(
                x + 1.3,
                y,
                short + "\n" + r["tier"].title(),
                ha="center",
                va="center",
                color="white",
                size=10,
            )
    ax.text(
        0,
        0.02,
        "Examples are qualitative policy distinctions; all require the stated anatomical and injury-context conditions.",
        size=10,
    )
    save(fig, "bulk-cell-distinctions")
    return paths


def native_inventory():
    rows = []
    for r in normal_tissue_mappings(species="human"):
        tissue = r["concept_id"] or r["source_tissue"]
        cellular = bool(r["source_cell_type"])
        finding = hazards.assess(
            tissue, r["source_cell_type"], resolution="cell" if cellular else "bulk"
        )
        rows.append(
            {
                **r,
                **{
                    f"hazard_{k}": v
                    for k, v in finding.items()
                    if k in ["tier", "assessment_id", "reason", "references"]
                },
            }
        )
    return pd.DataFrame(rows)


def delta(out):
    before_path = out / "baseline-classifications.csv.gz"
    current = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
    if not before_path.exists() or not current.exists():
        return pd.DataFrame()
    before = pd.read_csv(before_path, low_memory=False)
    after = pd.read_csv(current, low_memory=False)
    columns = ["protein_equivalence_group_id", "category", "classification"]
    joined = (
        before[[*columns, "exclusion_reasons"]]
        .rename(columns={"exclusion_reasons": "prior_exclusion_reasons"})
        .merge(
            after[[*columns, "protein_group_label", "exclusion_reasons"]],
            on=columns[0],
            suffixes=("_before", "_after"),
            validate="one_to_one",
            how="outer",
        )
    )
    changed = joined[joined.classification_before.ne(joined.classification_after)].fillna("")
    write_csv(changed, out / "classification-changes.csv")
    peaks_path = current.parent / "normal-peak-evidence.csv.gz"
    if peaks_path.exists():
        peaks = pd.read_csv(peaks_path, low_memory=False)
        labels = after.set_index("representative_gene_id").protein_group_label
        peaks["label"] = peaks.gene_id.map(labels)
        peaks = peaks[peaks.label.isin(changed.protein_group_label)]
        write_csv(
            peaks[
                [
                    "label",
                    "gene_id",
                    "assay",
                    "source",
                    "tissue",
                    "cell_type",
                    "selection_tier",
                    "value",
                    "unit",
                    "candidate_ceiling",
                    "negative_ihc_grade",
                    "unsafe_evidence",
                ]
            ],
            out / "changed-protein-rna-contexts.csv",
        )
    return changed


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    doc = hazards.registry()
    figs = figures(out, doc)
    records = doc["assessments"]
    refs = doc["references"]
    anatomy = hazards.anatomy_inventory()
    native = native_inventory()
    write_csv(hazards.assessments(), out / "assessments.csv")
    write_csv(hazards.references(), out / "references.csv")
    write_csv(anatomy, out / "canonical-anatomy.csv")
    write_csv(native, out / "native-mappings.csv")
    changes = delta(out)
    (out / "registry.json").write_text(json.dumps(doc, indent=2) + "\n")
    counts = f"{len(records)} curated rules; {len(refs)} references; {len(anatomy)} human anatomy concepts; {len(native)} native mapping rows."
    lines = [
        "# Normal-tissue hazard ontology",
        "",
        "Version 1.29.0 | reviewed " + doc["review_date"],
        "",
        counts,
        "",
        INTRO,
        "",
        "## How the policy is applied",
        "",
        POLICY,
        "",
    ]
    for tier, definition in doc["tier_definitions"].items():
        lines += [f"**{tier.title()}:** {definition}", ""]
    for f in figs:
        lines += [f"![{f.replace('-', ' ').title()}]({f}.png)", ""]
    lines += [
        "## Inventories and reproducibility",
        "",
        "[PDF report](tissue-hazard-ontology.pdf) · [Interactive report](report.html) · [Assessments](assessments.csv) · [References](references.csv) · [Every human anatomy concept](canonical-anatomy.csv) · [Native mappings](native-mappings.csv) · [Registry JSON](registry.json) · [Classification changes](classification-changes.csv)",
        "",
        f"{len(changes)} protein-group classifications differ from the 1.28.0 baseline. See the complete change table; counts are policy outcomes, not clinical events.",
        "",
        "The change table includes prior and current exclusion reasons. [Current RNA contexts for changed proteins](changed-protein-rna-contexts.csv) retain native units, tiers and ceilings. These are source-specific peak summaries, not every observed cell. The rare-cell screen supplies separate detection-frequency evidence.",
        "",
        "Reproduce: `python scripts/build_tissue_hazard_report.py`. Figures show rules and evidence designs, never calibrated toxicity scores. Native mappings record coverage, not absence of expression.",
        "",
    ]
    for system in sorted({r["system"] for r in records}):
        lines += ["## " + system, ""]
        for r in [r for r in records if r["system"] == system]:
            links = ", ".join(f"[{k}]({refs[k]['url']})" for k in r["references"])
            lines += [
                f"### {r['label']} — {r['tier'].title()}",
                "",
                f"`{r['id']}` | {r['level']} rule | adult human | {r['inference'].replace('_', ' ')}",
                "",
                r["rationale"],
                "",
                "**Acute loss:** " + r["acute_loss"],
                "",
                "**Persistent injury:** " + r["persistent_loss"],
                "",
                "**Repair:** " + r["repair"],
                "",
                "**Compensation:** " + r["management"],
                "",
                "**Recognition:** " + r["recognition_evidence"],
                "",
                "**Evidence scope:** " + r["evidence_scope"],
                "",
                "**Evidence:** " + links,
                "",
            ]
    lines += ["## Evidence ledger", ""]
    for k, r in refs.items():
        lines += [
            f"### {k}",
            "",
            f"[{r['title']}]({r['url']})",
            "",
            f"Model: {r['model']}; design: {r['design'].replace('_', ' ')}.",
            "",
            r["finding"],
            "",
            "**Limit:** " + r["limitation"],
            "",
        ]
    (out / "index.md").write_text("\n".join(lines) + "\n")
    render_html(out, doc, counts, figs, changes)
    render_pdf(out, doc, counts, figs, changes)
    source_files = [
        ROOT / "oncoref/data/normal-tissue-hazards.yaml",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/normal-tissue-mappings.csv",
        ROOT / "oncoref/tissue_hazards.py",
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        ROOT / "docs/audits/cta-selection-20261028/normal-peak-evidence.csv.gz",
        Path(__file__),
    ]

    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()

    manifest = {
        "schema": "oncoref.tissue-hazard-report.v1",
        "registry_id": doc["registry_id"],
        "review_date": doc["review_date"],
        "baseline": {
            "version": "1.28.0",
            "revision": "24cca2808de612940036df74a12243af43b5ca1f",
            "file": "baseline-classifications.csv.gz",
        },
        "counts": {
            "assessments": len(records),
            "references": len(refs),
            "anatomy": len(anatomy),
            "native_mappings": len(native),
            "classification_changes": len(changes),
        },
        "inputs": {str(p.relative_to(ROOT)): sha(p) for p in source_files},
        "outputs": {
            p.name: sha(p)
            for p in sorted(out.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(counts, "Classification changes:", len(changes), flush=True)


def render_html(out, doc, counts, figs, changes):
    e = html.escape
    cards = []
    for r in doc["assessments"]:
        links = " · ".join(
            f'<a href="{e(doc["references"][k]["url"])}">{e(k)}</a>' for k in r["references"]
        )
        content = "".join(
            f"<p><b>{label}</b> {e(r[field])}</p>"
            for label, field in [
                ("Immediate loss:", "acute_loss"),
                ("Persistent injury:", "persistent_loss"),
                ("Repair:", "repair"),
                ("Compensation:", "management"),
                ("Evidence scope:", "evidence_scope"),
            ]
        )
        cards.append(
            f'<article data-tier="{r["tier"]}" data-search="{e((r["label"] + " " + r["system"] + " " + r["rationale"]).lower())}"><span class="tag {r["tier"]}">{r["tier"].title()}</span><small>{e(r["system"])} / {r["level"]}</small><h3>{e(r["label"])}</h3><p>{e(r["rationale"])}</p><details><summary>Consequences, repair and evidence</summary>{content}<p>{e(r["recognition_evidence"])}</p><p>{links}</p><code>{e(r["id"])}</code></details></article>'
        )
    refcards = "".join(
        f'<p><b>{e(k)}</b> · {e(r["model"])} · {e(r["design"].replace("_", " "))}<br><a href="{e(r["url"])}">{e(r["title"])}</a><br>{e(r["finding"])}<br><em>Limit: {e(r["limitation"])}</em></p>'
        for k, r in doc["references"].items()
    )
    images = "".join(
        f'<figure><img src="{name}.svg" alt="{e(name.replace("-", " "))}"></figure>'
        for name in figs
    )

    def paragraphs(s):
        return "".join("<p>" + e(p.replace("\n", " ")) + "</p>" for p in s.strip().split("\n\n"))

    style = """body{font:16px/1.6 system-ui,sans-serif;margin:0;color:#172b3a;background:#f4f5f3}main{max-width:1150px;margin:auto;padding:48px 24px}h1{font-size:42px;line-height:1.15}h2{margin-top:48px}a{color:#245f75}header{max-width:860px}small{display:block;color:#58656b}.filters{position:sticky;top:0;background:#f4f5f3;padding:16px 0;z-index:1;display:flex;gap:12px}input,select{padding:12px;font:inherit;border:1px solid #aab5b7;border-radius:6px}input{flex:1}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:18px}article{padding:22px;background:white;border:1px solid #dbe0df;border-radius:10px}h3{margin:8px 0}article p{font-size:14px}.tag{float:right;color:white;padding:2px 9px;border-radius:4px;font-size:12px}.vital{background:#963d4b}.protected{background:#b27324}.manageable{background:#307b89}.disposable{background:#6a7197}summary{cursor:pointer;font-weight:600;font-size:14px}figure{margin:24px 0;background:white;padding:12px;border:1px solid #dbe0df}img{width:100%;height:auto}code{font-size:12px}.ledger p{padding:14px 0;border-bottom:1px solid #dbe0df}article[hidden]{display:none}@media print{.filters{position:static}details{display:block}}"""
    script = """const q=document.querySelector('#q'),t=document.querySelector('#tier');function filter(){let n=0;document.querySelectorAll('article').forEach(a=>{a.hidden=!(a.dataset.search.includes(q.value.toLowerCase())&&(!t.value||a.dataset.tier===t.value));if(!a.hidden)n++});document.querySelector('#n').textContent=n+' assessments shown'}q.addEventListener('input',filter);t.addEventListener('change',filter);filter();"""
    (out / "report.html").write_text(
        f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>OncoRef tissue-hazard ontology 1.29.0</title><style>{style}</style><main><header><small>ONCOREF / EVIDENCE &amp; POLICY / 1.29.0</small><h1>What normal-tissue injury can cost</h1><p>{e(counts)}</p>{paragraphs(INTRO)}<p><a href="tissue-hazard-ontology.pdf">Download PDF</a> · <a href="assessments.csv">Assessment data</a> · <a href="classification-changes.csv">Classification changes ({len(changes)})</a></p></header><h2>How the ontology works</h2>{paragraphs(POLICY)}{images}<h2>Browse every assessment</h2><div class="filters"><input id="q" aria-label="Search assessments" placeholder="Search an organ, cell or rationale"><select id="tier" aria-label="Filter tier"><option value="">All tiers</option>'
        + "".join(f'<option value="{t}">{t.title()}</option>' for t in hazards.TIERS)
        + f'</select></div><p id="n"></p><div class="cards">{"".join(cards)}</div><h2>Evidence ledger</h2><div class="ledger">{refcards}</div><p>Policy fingerprint: <code>{e(doc["registry_id"])}</code></p></main><script>{script}</script></html>'
    )


def render_pdf(out, doc, counts, figs, changes):
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="BodySmall",
            fontName="Helvetica",
            fontSize=9,
            leading=12,
            spaceAfter=6,
            textColor=colors.HexColor("#172B3A"),
        )
    )
    styles.add(
        ParagraphStyle(
            name="Entry",
            fontName="Helvetica-Bold",
            fontSize=12,
            leading=15,
            spaceBefore=12,
            spaceAfter=6,
            keepWithNext=True,
        )
    )
    styles["Title"].textColor = colors.HexColor("#172B3A")
    styles["Heading1"].textColor = colors.HexColor("#172B3A")
    styles["Heading1"].spaceBefore = 16
    e = html.escape

    def p(text):
        return Paragraph(text, styles["BodySmall"])

    story = [
        Paragraph("Normal-tissue hazard ontology", styles["Title"]),
        p("OncoRef 1.29.0 | Reviewed " + doc["review_date"]),
        p(counts),
        Spacer(1, 0.12 * inch),
    ]
    story += [p(e(s.replace("\n", " "))) for s in INTRO.strip().split("\n\n")]
    story.append(Image(str(out / "ontology-logic.png"), width=7 * inch, height=4.1 * inch))
    story.append(PageBreak())
    story.append(Paragraph("Interpreting the four tiers", styles["Heading1"]))
    for t, d in doc["tier_definitions"].items():
        story.append(p(f'<font color="{COLORS[t]}"><b>{t.title()}</b></font> - {e(d)}'))
    story += [p(e(s.replace("\n", " "))) for s in POLICY.strip().split("\n\n")]
    story.append(
        p(
            f"<b>Release comparison:</b> {len(changes)} protein-group classifications differ from the recorded 1.28.0 baseline. The complete CSV identifies each group and its current reasons. These are screening outcomes, not clinical event rates."
        )
    )
    for f in figs[1:]:
        story.append(PageBreak())
        story.append(Paragraph(f.replace("-", " ").title(), styles["Heading1"]))
        from PIL import Image as PILImage

        im = PILImage.open(out / f"{f}.png")
        w, h = im.size
        story.append(Image(str(out / f"{f}.png"), width=7 * inch, height=7 * inch * h / w))
    for system in sorted({r["system"] for r in doc["assessments"]}):
        story += [PageBreak(), Paragraph(system, styles["Heading1"])]
        for r in [r for r in doc["assessments"] if r["system"] == system]:
            heading = Paragraph(
                e(r["label"]) + f' <font color="{COLORS[r["tier"]]}"> / {r["tier"].title()}</font>',
                styles["Entry"],
            )
            body = [
                heading,
                p(e(r["rationale"])),
                p(
                    "<b>Acute:</b> "
                    + e(r["acute_loss"])
                    + " <b>Persistent:</b> "
                    + e(r["persistent_loss"])
                ),
                p(
                    "<b>Repair:</b> "
                    + e(r["repair"])
                    + " <b>Compensation:</b> "
                    + e(r["management"])
                ),
                p("<b>Evidence scope:</b> " + e(r["evidence_scope"])),
                p(
                    "<b>Evidence:</b> "
                    + ", ".join(
                        f'<link href="{e(doc["references"][k]["url"])}" color="#245f75">{e(k)}</link>'
                        for k in r["references"]
                    )
                ),
                p(
                    f'<font color="#64727a">{r["id"]} | {r["level"]} | adult human | policy inference; recognition evaluated separately</font>'
                ),
            ]
            story.append(KeepTogether(body))
    story += [PageBreak(), Paragraph("Evidence ledger and limitations", styles["Heading1"])]
    for k, r in doc["references"].items():
        story.append(
            KeepTogether(
                [
                    Paragraph(e(k), styles["Entry"]),
                    p(f'<link href="{e(r["url"])}" color="#245f75">{e(r["title"])}</link>'),
                    p(e(r["model"]) + " / " + e(r["design"].replace("_", " "))),
                    p(e(r["finding"])),
                    p("<b>Limit:</b> " + e(r["limitation"])),
                ]
            )
        )

    def footer(canvas, _doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#64727A"))
        canvas.drawString(
            0.65 * inch,
            0.35 * inch,
            "OncoRef 1.29.0 | Qualitative policy, not a toxicity probability",
        )
        canvas.drawRightString(7.85 * inch, 0.35 * inch, str(_doc.page))
        canvas.restoreState()

    pdf = SimpleDocTemplate(
        str(out / "tissue-hazard-ontology.pdf"),
        pagesize=(8.5 * inch, 11 * inch),
        rightMargin=0.65 * inch,
        leftMargin=0.65 * inch,
        topMargin=0.6 * inch,
        bottomMargin=0.65 * inch,
        title="OncoRef normal-tissue hazard ontology 1.29.0",
        author="OncoRef",
        invariant=1,
    )
    pdf.build(story, onFirstPage=footer, onLaterPages=footer)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT)
    build(parser.parse_args().out)
