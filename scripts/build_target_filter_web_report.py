"""Render the web summary of categories, funnels, curation and known-CTA placements.

Writes docs/audits/target-filters-20261112/report.html and known-ctas.csv.
Usage: PYTHONPATH=. python scripts/build_target_filter_web_report.py
"""

import html
import json
from pathlib import Path

import pandas as pd

from oncoref import target_curation, target_selection
from oncoref import target_filters as tf
from oncoref.version import __version__

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/target-filters-20261112"
FUNNEL = ROOT / "docs/audits/target-selection-20261110/rebuild/funnel.csv"
KNOWN = ROOT / "oncoref/data/cancer-testis-antigens.csv"
CATEGORIES = ("CTA", "Cancer Antigen only", "Reproductive Antigen only", "Other", "Unsafe")
COLORS = {
    "CTA": "#2f7f6f",
    "Cancer Antigen only": "#3e8bba",
    "Reproductive Antigen only": "#d49b42",
    "Other": "#9aa5ad",
    "Unsafe": "#b5535c",
}
MODES = ("loose", "standard", "strict")
LAYER_NAMES = {
    "off_target": "Off target",
    "warnings": "Warnings",
    "curation": "Curation",
    "sets": "Sets",
    "papers": "Papers",
    "exclusions": "Exclusions",
}
STYLE = """body{font:17px/1.6 system-ui,sans-serif;color:#23313a;background:#fff;margin:0}
main{max-width:1060px;padding:35px 24px;margin:auto}h1{font-size:36px;line-height:1.2}
h2{margin-top:2.2em}.status{border-left:4px solid #a36a20;background:#faf6ec;padding:15px}
table{width:100%;border-collapse:collapse;margin:14px 0}
td,th{padding:8px 10px;text-align:left;border-bottom:1px solid #d9dfe3;vertical-align:top}
th{background:#eff3f5}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
code{font-size:.9em}a{color:#176689}.note{color:#53636c}
.bar{height:22px;border-radius:3px;min-width:2px}.barcell{width:42%}
.pill{display:inline-block;padding:1px 9px;border-radius:11px;color:#fff;font-size:.82em;white-space:nowrap}
.ok{color:#2f7f6f;font-weight:600}.out{color:#8b3039}
.scroll{max-height:560px;overflow:auto;border:1px solid #d9dfe3}
.scroll table{margin:0}.scroll th{position:sticky;top:0}
input{font:inherit;padding:6px 10px;width:320px;max-width:100%;margin:8px 0}"""


def esc(value):
    return html.escape(str(value))


def table(headers, rows, numeric=()):
    head = "".join(
        f'<th class="{"n" if i in numeric else ""}">{h}</th>' for i, h in enumerate(headers)
    )
    body = "".join(
        "<tr>"
        + "".join(f'<td class="{"n" if i in numeric else ""}">{c}</td>' for i, c in enumerate(row))
        + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def bar(value, total, color):
    return f'<div class="bar" style="width:{100 * value / total:.2f}%;background:{color}"></div>'


def pill(category):
    return f'<span class="pill" style="background:{COLORS[category]}">{esc(category)}</span>'


def outcome(text):
    if text == "✓":
        return '<span class="ok">kept</span>'
    if text == "—":
        return '<span class="note">—</span>'
    return f'<span class="out">{esc(LAYER_NAMES.get(text, text)).lower()}</span>'


def known_ctas(measured, applied):
    """CTpedia genes mapped to catalog proteins, with each mode's outcome."""
    known = pd.read_csv(KNOWN, low_memory=False)
    known = known[known.source_databases.fillna("").str.contains("CTpedia")]
    by_gene = {}
    for index, members in zip(measured.index, measured.member_gene_ids):
        for gene_id in json.loads(members) if isinstance(members, str) else list(members):
            by_gene[gene_id] = index
    known = known.assign(row=known.Ensembl_Gene_ID.map(by_gene))
    symbols = known.groupby("row").Symbol.agg(lambda s: ", ".join(sorted(s)))
    out = measured.loc[symbols.index.astype(int)].copy()
    out["ctpedia_symbols"] = symbols.to_numpy()
    for mode in MODES:
        result = applied[mode].loc[out.index]
        verdict = result.filter_failed.replace("", "✓")
        out[mode] = verdict.where(out.is_target_candidate, "—")
    return len(known), out.sort_values("protein_group_label")


def main():
    measured = tf.measures(target_selection.results())
    candidates = measured[measured.is_target_candidate]
    applied = {mode: tf.apply(f, measured) for mode, f in tf.FILTER_MODES.items()}
    counts = measured.category.value_counts()
    total = len(measured)

    category_rows = [
        (pill(c), f"{int(counts[c]):,}", bar(int(counts[c]), total, COLORS[c])) for c in CATEGORIES
    ]

    steps = pd.read_csv(FUNNEL)
    funnel_rows = [("All proteins assessed", "", "", f"{total:,}", bar(total, total, "#5b6b75"))]
    for r in steps.itertuples():
        funnel_rows.append(
            (
                esc(r.step),
                esc(r.result),
                f"−{r.removed:,}",
                f"{r.remaining:,}",
                bar(r.remaining, total, "#5b6b75"),
            )
        )

    mode_rows = []
    for mode in tf.FILTER_MODES:
        kept = applied[mode][applied[mode].keep]
        by = kept.filtered_category.value_counts()
        mode_rows.append(
            (
                f"<code>{mode}</code>",
                f"{len(kept):,}",
                *(f"{int(by.get(c, 0)):,}" for c in CATEGORIES[:3]),
                bar(len(kept), len(candidates), COLORS["CTA"]),
            )
        )
    layer_rows = []
    for mode in MODES:
        steps_ = tf.funnel(mode, measured)
        layer_rows.append(
            (f"<code>{mode}</code>", *(f"{int(r.remaining):,}" for r in steps_.itertuples()))
        )
    rule_blocks = "".join(
        f"<h3><code>{mode}</code></h3><ul>"
        + "".join(f"<li>{esc(line)}</li>" for line in tf.describe(mode).splitlines()[1:])
        + "</ul>"
        for mode in MODES
    )

    comps = target_curation.compartments().set_index("compartment")
    removed_by = {
        "vital": "<code>standard</code>, <code>strict</code>",
        "protected": "<code>strict</code>",
        "manageable": "no mode",
        "disposable": "never",
    }
    tier_rows = []
    for tier in target_curation.TIERS:
        rows = candidates[candidates.curated_tier.eq(tier)]
        keys = sorted(
            {
                k
                for text in rows.curated_compartments
                for k in text.split("; ")
                if comps.tier[k] == tier
            }
        )
        tier_rows.append(
            (tier.capitalize(), len(rows), esc(", ".join(comps.label[keys])), removed_by[tier])
        )
    proposed = comps[comps.tier_source == "proposed"]
    proposed_rows = [
        (esc(r.label), r.tier.capitalize(), r.confidence, esc(r.tier_references))
        for r in proposed.itertuples()
    ]
    families = target_curation.families()
    family_rows = [
        (
            esc(r.label),
            int(candidates.protein_group_label.str.match(r.pattern).sum()),
            esc(comps.label[r.compartment]),
            r.tier.capitalize() if r.acts else f"{r.tier.capitalize()} (recorded only)",
            esc(r.sources),
        )
        for r in families.itertuples()
    ]
    acting = target_curation.findings()
    acting = acting[acting.acts.astype(str).str.lower().eq("true")]
    acting = acting[
        acting.tier.isin(["vital", "protected"]) & ~acting.gene.str.match(r"^(OR\d|TAS2R)")
    ]
    finding_rows = [
        (
            esc(r.gene),
            r.tier.capitalize(),
            esc(r.tissue),
            esc(r.claim[:260]),
            esc(r.best_source[:120]),
        )
        for r in acting.drop_duplicates("gene")
        .sort_values(["tier", "gene"], ascending=[False, True])
        .itertuples()
    ]
    loci = target_curation.readthrough()

    genes, known = known_ctas(measured, applied)
    known_counts = known.category.value_counts()
    placement_rows = [
        (
            pill(c),
            int(known_counts.get(c, 0)),
            *(
                int((known.loc[known.category.eq(c), mode] == "✓").sum())
                if c in CATEGORIES[:3]
                else "—"
                for mode in MODES
            ),
        )
        for c in CATEGORIES
    ]
    placement_rows.append(
        ("<strong>All</strong>", len(known), *(int((known[mode] == "✓").sum()) for mode in MODES))
    )
    reasons = (
        known[known.category.eq("Unsafe")].Unsafe_reasons.str.split(" | ", regex=False).explode()
    )
    reason_rows = [(esc(k), int(v)) for k, v in reasons.value_counts().items()]
    known_rows = [
        (
            esc(r.protein_group_label),
            pill(r.category),
            f"{r.off_target_tpm:.1f}" if pd.notna(r.off_target_tpm) else "",
            ""
            if pd.isna(r.independent_patients_above_normal)
            else f"{int(r.independent_patients_above_normal):,}",
            int(r.papers_cta),
            esc(r.curated_tier),
            outcome(r.loose),
            outcome(r.standard),
            outcome(r.strict),
            esc(r.Unsafe_reasons if isinstance(r.Unsafe_reasons, str) else ""),
        )
        for r in known.itertuples()
    ]
    known[
        [
            "protein_group_label",
            "ctpedia_symbols",
            "category",
            "Unsafe_reasons",
            "off_target_tpm",
            "off_target_ntpm",
            "independent_patients_above_normal",
            "papers_cta",
            "curated_tier",
            "curated_compartments",
            *MODES,
        ]
    ].to_csv(OUT / "known-ctas.csv", index=False)

    page = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Target candidates: categories, funnel and filters — OncoRef {__version__}</title>
<style>{STYLE}</style><main>
<h1>Target candidates: categories, funnel and filters</h1>
<p><strong>Of {total:,} proteins assessed, {len(candidates)} are Target candidates.
Filter modes narrow them to {mode_rows[1][1]}, {mode_rows[2][1]} and {mode_rows[3][1]}.</strong></p>
<p class="status">OncoRef {__version__} · catalog <code>{target_selection.definition()["policy_id"]}</code> ·
filters <code>{tf.FILTER_POLICY_ID}</code> · curation <code>{target_curation.registry()["policy"]}</code>.<br>
These are outputs of screening rules. They are not validated tumor antigens.</p>
<p><a href="index.md">Full filter report</a> · <a href="normal-function-sweep.md">Literature sweep</a> ·
<a href="../../target-filters.md">Rules reference</a> · <a href="../target-selection-20261110/report.html">Catalog rules</a> ·
<a href="known-ctas.csv">Known CTA table (CSV)</a></p>

<h2>1. Categories</h2>
<p>CTA = Cancer Antigens ∩ Reproductive Antigens. Target candidates = Cancer Antigens ∪ Reproductive Antigens:
the first three rows.</p>
{table(["Category", "Proteins", ""], category_rows, numeric=(1,))}

<h2>2. Catalog funnel</h2>
<p>How {total:,} proteins become {len(candidates)} candidates. Each step removes proteins into the category shown.</p>
{table(["Step", "Removed into", "Removed", "Remaining", ""], funnel_rows, numeric=(2, 3))}

<h2>3. Filter modes</h2>
<p>Each mode keeps a subset of the one before. Every mode keeps {", ".join(tf.POSITIVE_CONTROLS)}.
Categories are the ones after the evidence knobs have tightened each set.</p>
{table(["Filter mode", "Kept", *CATEGORIES[:3], ""], mode_rows, numeric=(1, 2, 3, 4))}
<p>Candidates remaining after each layer, in order:</p>
{table(["Filter mode", "Start", *LAYER_NAMES.values()], layer_rows, numeric=tuple(range(1, 8)))}
{rule_blocks}

<h2>4. Curation layer</h2>
<p>Cited literature places some candidates in normal adult tissue that the expression panels do not
sample. Each tissue carries a tier from the hazard registry, or a cited proposal where the registry had none.</p>
{table(["Curated tier", "Candidates", "Tissues", "Removed in"], tier_rows, numeric=(1,))}
<h3>Family rules</h3>
{table(["Family", "Candidates", "Normal site", "Tier", "Sources"], family_rows, numeric=(1,))}
<h3>Tiers proposed for tissues the registry did not cover</h3>
{table(["Compartment", "Proposed tier", "Confidence", "References"], proposed_rows)}
<h3>Single proteins with a confirmed finding in a Vital or Protected tissue</h3>
{table(["Protein", "Tier", "Where", "Finding", "Source"], finding_rows)}
<h3>Read-through loci</h3>
<p>{len(loci)} candidates share exons with two or more other protein-coding genes and are removed in every
mode but <code>all</code>: {esc(", ".join(loci.protein_group_label))}.</p>

<h2>5. Where the known CTAs end up</h2>
<p>CTpedia lists {genes} genes. They map to {len(known)} catalog proteins, because identical copies are
grouped into one protein.</p>
{table(["Catalog category", "Proteins", *(f"Kept by <code>{m}</code>" for m in MODES)], placement_rows, numeric=(1, 2, 3, 4))}
<h3>Why {int(known_counts.get("Unsafe", 0))} are Unsafe</h3>
<p>A protein can have more than one reason.</p>
{table(["Reason", "Proteins"], reason_rows, numeric=(1,))}
<h3>Every known CTA</h3>
<p class="note">A mode column shows “kept” or the first filter layer that removes the protein.</p>
<input id="q" placeholder="Filter by name, category or reason" oninput="f(this.value)">
<div class="scroll" id="known">
{table(["Protein", "Category", "Off-target TPM", "Tumors above normal", "CTA papers", "Curated tier", "loose", "standard", "strict", "Unsafe reasons"], known_rows, numeric=(2, 3, 4))}
</div>
<script>function f(q){{q=q.toLowerCase();document.querySelectorAll('#known tbody tr').forEach(function(r){{
r.style.display=r.textContent.toLowerCase().indexOf(q)<0?'none':''}})}}</script>
</main></html>
"""
    (OUT / "report.html").write_text(page)
    print(f"{genes} CTpedia genes, {len(known)} proteins")
    print(known.category.value_counts().to_dict())
    for mode in MODES:
        failed = known.loc[known.is_target_candidate, mode].value_counts().to_dict()
        print(mode, failed)


if __name__ == "__main__":
    main()
