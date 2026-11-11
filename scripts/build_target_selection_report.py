#!/usr/bin/env python3
"""Render the target-set contract and link its measured catalog rebuild."""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyBboxPatch

from oncoref import target_selection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/target-selection-20261110"


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    path = OUT / f"{name}.svg"
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=160)
    plt.close(fig)


def figures():
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 11, "svg.hashsalt": "oncoref-target-sets-v1"}
    )
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set(xlim=(0, 11), ylim=(0, 7))
    ax.axis("off")
    ax.text(
        5.5, 6.7, "One shared screen, two antigen sets", ha="center", fontsize=19, weight="bold"
    )
    ax.text(
        5.5,
        6.25,
        "All proteins → normal screen + coverage + identity/translation + bounds + sequence",
        ha="center",
    )
    ax.add_patch(
        FancyBboxPatch(
            (0.2, 0.7),
            10.6,
            4.9,
            boxstyle="round,pad=0.02",
            facecolor="#fafafa",
            edgecolor="#8993a0",
        )
    )
    ax.text(5.5, 5.25, "Eligible proteins", ha="center", weight="bold", fontsize=13)
    ax.add_patch(Ellipse((4.1, 3.2), 5.3, 3.3, color="#3e8bba", alpha=0.26))
    ax.add_patch(Ellipse((6.9, 3.2), 5.3, 3.3, color="#d49b42", alpha=0.3))
    ax.text(3.2, 3.7, "Cancer Antigens", ha="center", fontsize=14, weight="bold")
    ax.text(3.05, 2.9, "Qualifying tumor\nevidence", ha="center", linespacing=1.5)
    ax.text(7.8, 3.7, "Reproductive Antigens", ha="center", fontsize=14, weight="bold")
    ax.text(7.95, 2.9, "Reproductive or\ndevelopmental evidence", ha="center", linespacing=1.5)
    ax.text(5.5, 3.2, "CTA", ha="center", fontsize=20, weight="bold", color="#1d3537")
    ax.text(5.5, 2.7, "intersection", ha="center", fontsize=10)
    ax.text(
        5.5,
        1.3,
        "Target candidates = the union of both sets",
        ha="center",
        fontsize=12,
        weight="bold",
    )
    ax.text(5.5, 0.95, "Neither association qualifies → Other", ha="center", fontsize=10)
    ax.text(
        5.5,
        5.85,
        "Accepted exclusion → Unsafe     |     Failed shared admission check → Other",
        ha="center",
        color="#8b3039",
    )
    ax.text(
        5.5,
        0.2,
        "Published target support annotates every category. Areas do not represent counts.",
        ha="center",
        fontsize=10,
    )
    save(fig, "sets")

    fig, ax = plt.subplots(figsize=(10, 5.4))
    n = np.linspace(0, 10, 501)
    ax.plot(
        n,
        np.maximum(n, target_selection.DETECTION_FLOOR),
        color="#247aa5",
        lw=2.5,
        label="Two independent patients: each ≥1 and >normal maximum",
    )
    ax.plot(
        n,
        target_selection.SINGLE_TUMOR_FOLD * np.maximum(n, target_selection.DETECTION_FLOOR),
        color="#a36a20",
        lw=2.5,
        label=f"One sample: ≥{target_selection.SINGLE_TUMOR_FOLD:g} × max(normal maximum, 1)",
    )
    ax.set(
        xlabel="Highest eligible normal value N (same comparison scale)",
        ylabel="Tumor expression required",
        xlim=(0, 10),
        ylim=(0, 105),
    )
    ax.set_title("Cancer RNA evidence: repeated elevation or one strong tumor", loc="left", pad=15)
    ax.spines[["right", "top"]].set_visible(False)
    ax.grid(alpha=0.15)
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    fig.text(
        0.13,
        -0.01,
        "Blue line: equality to N does not pass. Both routes still require the shared normal screen.",
        fontsize=10,
    )
    save(fig, "tumor-rule")


def svg(name):
    text = (OUT / f"{name}.svg").read_text()
    # Standalone figures use Matplotlib's repeated figure/font IDs. Namespace
    # them when embedding both SVGs in one HTML document.
    text = re.sub(r' id="([^"]+)"', rf' id="{name}-\1"', text)
    text = text.replace('xlink:href="#', f'xlink:href="#{name}-')
    text = text.replace("url(#", f"url(#{name}-")
    return "<svg" + text.split("<svg", 1)[1]


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    definition = target_selection.definition()
    (OUT / "definition.json").write_text(json.dumps(definition, indent=2) + "\n")
    figures()
    rows = "".join(
        f"<tr><td>{n:g}</td><td>{'≥1' if n < 1 else f'&gt;{n:g}'} each</td>"
        f"<td>≥{target_selection.SINGLE_TUMOR_FOLD * max(n, target_selection.DETECTION_FLOOR):g}</td></tr>"
        for n in (0, 0.4, 1, 3, 8)
    )
    rare = definition["rare_cell_rna"]
    fold, cap = target_selection.SINGLE_TUMOR_FOLD, target_selection.DISPOSABLE_CAP
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cancer and Reproductive Antigens — OncoRef</title>
<style>body{{font:17px/1.6 system-ui,sans-serif;color:#23313a;background:#fff;margin:0}}
main{{max-width:1020px;padding:35px 24px;margin:auto}}h1{{font-size:36px;line-height:1.2}}
h2{{margin-top:2em}}.status{{border-left:4px solid #a36a20;background:#faf6ec;padding:15px}}
svg{{width:100%;height:auto}}figure{{margin:25px 0}}table{{width:100%;border-collapse:collapse}}
td,th{{padding:10px;text-align:left;border-bottom:1px solid #d9dfe3}}th{{background:#eff3f5}}
code{{font-size:.9em}}a{{color:#176689}}.note{{color:#53636c}}</style><main>
<h1>Cancer and Reproductive Antigens</h1>
<p><strong>CTA = Cancer Antigens ∩ Reproductive Antigens.</strong><br>
Target candidates = Cancer Antigens ∪ Reproductive Antigens.</p>
<p class="status">Rebuilt catalog: <code>{html.escape(definition["policy_id"])}</code>.
{", ".join(f"{target_selection.metadata()['summary']['categories'][c]:,} {c}" for c in target_selection.CATEGORIES)}.
<a href="rebuild/report.html">Filterable protein table, measurements and category changes</a>.
The released 1.36.0 catalog remains a separate snapshot under its earlier policy.</p>
<figure>{svg("sets")}</figure>
<h2>One shared eligibility screen</h2>
<p>Accepted normal-tissue or applicable clinical exclusion → <strong>Unsafe</strong>.
Required coverage, identity, translation or unresolved evidence → <strong>Other:Review</strong>.
Nested or Unsafe-contained sequence → <strong>Other:Sequence</strong>.
Target maximum or Disposable cap failure → <strong>Other:Expression</strong>.
The remaining proteins enter the two independent association checks.</p>
<h2>Rare normal-cell RNA</h2>
<p>Each adult donor is assessed separately in each named Vital or Protected cell population.
A count matrix cannot show where a read aligned: one mis-assigned read per cell from a paralog,
pseudogene or overlapping gene looks like sparse expression. Cells with a single molecule
therefore never change a category.</p>
<table><thead><tr><th>Question</th><th>Requirement</th></tr></thead><tbody>
<tr><td>Is there a signal?</td><td>{html.escape(rare["signal"])}</td></tr>
<tr><td>Is it replicated?</td><td>{html.escape(rare["replication"])} {html.escape(rare["single_study"])}</td></tr>
<tr><td>Could the reads belong elsewhere?</td><td>{html.escape(rare["mapping_risk"])}</td></tr>
<tr><td>Does another assay agree?</td><td>{html.escape(rare["confirmation"])}</td></tr></tbody></table>
<table><thead><tr><th>Grade</th><th>Meaning</th><th>Vital</th><th>Protected</th></tr></thead>
<tbody><tr><td>Established</td><td>{html.escape(rare["established"])}</td><td>Unsafe</td><td>Other:Review</td></tr>
<tr><td>Unresolved</td><td>{html.escape(rare["unresolved"])}</td><td>Other:Review</td><td>Recorded</td></tr>
<tr><td>Below the screen</td><td>Neither</td><td colspan="2">No finding. This is not evidence of absence or safety.</td></tr></tbody></table>
<p>1%, 2% and 5% detection remain reported comparison metrics. {html.escape(rare["release"])}
<a href="../normal-cell-atlas-20261110/rare-cell/">Donor-level evidence</a>.</p>
<h2>Cancer evidence</h2>
<p>Let N be the highest qualified nonexempt normal measurement in a complete, compatible
comparison panel. Two independent patients each need T ≥1 and T &gt;N; alternatively,
one tumor needs T ≥{fold:g} × max(N, 1). Sample quality and exact protein identity must pass.
The highest tumor, reproductive or developmental value must also exceed every measured
nonreproductive normal and Disposable value, and Disposable expression must stay at or
below {cap:g} TPM and {cap:g} nTPM, each on its own scale.</p>
<table><thead><tr><th>Normal maximum N</th><th>Two independent tumors</th><th>One tumor</th></tr></thead>
<tbody>{rows}</tbody></table>
<figure>{svg("tumor-rule")}</figure>
<p>Do not divide tumor TPM by an unrelated normal nTPM or CP10K maximum. Record the
comparison methods, normal contexts, source statistics and sample identities.
Missing required normal data leave the comparison unassessed. All other normal panels
continue to enforce their exclusions.
<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5903355/">Wang et al., 2018</a>.</p>
<p>Reviewed reliable tumor protein is a separate positive evidence route. Qualitative
staining does not supply an RNA fold change. A high tumor value never cancels a normal exclusion.</p>
<h2>Reproductive and developmental evidence</h2>
<p>Adult reproductive, placental, extraembryonic, embryonic and fetal associations count.
Use ≥1 native TPM/nTPM with the existing cell/read support, or reliable mapped protein
evidence. Record the actual branch and stage. The old adult reproductive concentration
score becomes descriptive; it is not an admission gate.</p>
<p>The 1-unit detection convention is distinct from protein production or antigen
presentation. <a href="https://v23.proteinatlas.org/about/assays%2Bannotation">HPA methods</a>.</p>
<h2>Published target support</h2>
<p>A curated target nomination adds supporting evidence. Preserve the exact protein–paper
link and what it establishes: nomination, biology, tumor expression, protein identity,
immune recognition, or normal-tissue counterevidence. Source-wide roles do not turn
every gene in a paper into a target.</p>
<p>Core papers cover the catalog efficiently; supplementary papers retain independent
experiments, identity limits and counterevidence. Both remain available through
<code>cta_literature.attestations()</code>, including links attached to Other and Unsafe proteins.</p>
<h2>Read the decision precisely</h2>
<p>“Cancer Antigen only” and “Reproductive Antigen only” describe membership in one qualified
set. They do not establish absence everywhere else. Unassessed associations stay visible.
Target candidates is the union's name; it does not assert clinical safety.</p>
<p><a href="../../target-selection/">Full definitions and API</a> ·
<a href="definition.json">Machine-readable contract</a></p>
</main></html>"""
    (OUT / "report.html").write_text(page)


if __name__ == "__main__":
    build()
