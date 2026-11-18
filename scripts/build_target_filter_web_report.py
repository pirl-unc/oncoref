"""Render the web summary of categories, funnels, curation and known-CTA placements.

Writes docs/audits/target-filters-20261112/report.html and known-ctas.csv.
Usage: PYTHONPATH=. python scripts/build_target_filter_web_report.py
"""

import html
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import target_curation, target_selection, tissue_tiers  # noqa: E402
from oncoref import target_filters as tf  # noqa: E402
from oncoref.version import __version__  # noqa: E402

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
    # Reported for reading; no knob acts on it.
    measured["off_target_ntpm"] = measured[["healthy_maximum_ntpm", "disposable_max_ntpm"]].max(
        axis=1, skipna=False
    )
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
    profiles = tf.profile_summary(measured)
    profile_rows = [
        (
            f"<code>{esc(r.profile)}</code>",
            esc(r.protected_compartments),
            f"<strong>{r.CTA:,}</strong>",
            f"{r.CA_only:,}",
            f"{r.RA_only:,}",
            f"{r.Union:,}",
            bar(r.CTA, int(counts["CTA"]), COLORS["CTA"]),
        )
        for r in profiles.itertuples()
    ]
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

    removed_by = {
        "vital": "<code>standard</code>, <code>strict</code>",
        "protected": "<code>strict</code>",
        "manageable": "no mode",
        "disposable": "never",
    }
    tier_rows = []
    for tier in target_curation.TIERS:
        rows = candidates[candidates.curated_tier.eq(tier)]
        tissues = sorted(
            {t.split(": ")[0] for text in rows.curated_compartments for t in text.split("; ")}
        )
        tier_rows.append((tier.capitalize(), len(rows), esc(", ".join(tissues)), removed_by[tier]))
    added = target_curation.compartments()
    references = tissue_tiers.references().set_index("reference_id").url
    proposed_rows = [
        (
            esc(r.label),
            r.tier.capitalize(),
            r.confidence,
            ", ".join(f'<a href="{references[k]}">{esc(k)}</a>' for k in r.references),
        )
        for r in added.itertuples()
    ]
    families = target_curation.families()
    family_counts = candidates.literature_family.value_counts()
    family_rows = [
        (
            esc(r.label),
            int(family_counts.get(r.id, 0)),
            esc(r.tissue),
            r.tier.capitalize() if r.filter_acts else f"{r.tier.capitalize()} (recorded only)",
            esc(r.sources),
        )
        for r in families.itertuples()
    ]
    entries = target_curation.entries()
    excluded = entries[entries.effect.eq("exclude") & entries.review.eq("candidate_sweep")]
    finding_rows = [
        (
            esc(r.protein),
            r.tissue_tier.capitalize(),
            esc(r.decision_tissue),
            esc(r.decision_reason),
            " ".join(
                f'<a href="{esc(u)}">source {i}</a>'
                for i, u in enumerate(r.sources.split(" | ")[:2], 1)
            ),
        )
        for r in excluded.sort_values(
            ["tissue_tier", "protein"], ascending=[False, True]
        ).itertuples()
    ]
    recorded = entries[
        entries.effect.eq("none") & entries.recorded_tier.isin(["vital", "protected"])
    ]
    recorded_rows = [
        (
            esc(r.protein),
            r.recorded_tier.capitalize(),
            esc(r.recorded_tissue),
            esc(r.decision_reason),
        )
        for r in recorded.sort_values(
            ["recorded_tier", "protein"], ascending=[False, True]
        ).itertuples()
    ]
    loci = target_curation.readthrough().drop_duplicates("entity_id")
    loci = loci.merge(measured[["entity_id", "classification"]], on="entity_id")
    moved = loci[loci.classification.eq("Other:Sequence")]

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
            "" if pd.isna(r.tumor_fold) else f"{r.tumor_fold:,.0f}",
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
            "tumor_fold",
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
filters <code>{tf.FILTER_POLICY_ID}</code> · tissue-tier updates <code>{target_curation.compartments().updates_id.iloc[0][-12:]}</code>.<br>
These are outputs of screening rules. They are not validated tumor antigens.</p>
<p><a href="index.md">Full filter report</a> · <a href="normal-function-sweep.md">Literature sweep</a> ·
<a href="../../target-filters.md">Rules reference</a> · <a href="../target-selection-20261110/report.html">Catalog rules</a> ·
<a href="known-ctas.csv">Known CTA table (CSV)</a></p>

<h2>1. Categories</h2>
<p>CTA = Cancer Antigens ∩ Reproductive Antigens. Target candidates = Cancer Antigens ∪ Reproductive Antigens:
the first three rows.</p>
{table(["Category", "Proteins", ""], category_rows, numeric=(1,))}

<h2>2. Catalog funnel</h2>
<p>How {total:,} proteins become {len(candidates)} candidates. Each step removes proteins into the category shown.
A Cancer Antigen needs one QC-passing tumor sample, from any comparable TPM cohort, at 10 × every non-reproductive normal value.</p>
{table(["Step", "Removed into", "Removed", "Remaining", ""], funnel_rows, numeric=(2, 3))}

<h2>3. Filter modes and host toxicity profiles</h2>
<p>Each mode keeps a subset of the one before. Every mode keeps {", ".join(tf.POSITIVE_CONTROLS)}.
Categories are the ones after the fold knob has tightened Cancer Antigens.</p>
{table(["Filter mode", "Kept", *CATEGORIES[:3], ""], mode_rows, numeric=(1, 2, 3, 4))}
<h3>Host anatomy and tissue-tolerance profiles</h3>
<p>Host selection is separate from biological membership. No reproductive tissue is protected by default. Explicit protection, accepted loss and absence are distinct choices. Missing coverage and unresolved protein findings fail the host screen; native TPM and nTPM have separate limits. Milk-program labels are inferred. <a href="../../host-selection/">Contract and limitations</a>.</p>
{table(["Host profile", "Protected compartments", "CTA", "CA only", "RA only", "CA ∪ RA", ""], profile_rows, numeric=(2, 3, 4, 5))}
<p>Candidates remaining after each layer, in order:</p>
{table(["Filter mode", "Start", *LAYER_NAMES.values()], layer_rows, numeric=tuple(range(1, 8)))}
{rule_blocks}

<h2>4. Curation: one set of records for catalog and filters</h2>
<p>The catalog and the filters read the same literature review and the same tissue-tier registry.
A gene-specific report in normal adult Vital or Protected tissue, with a quote retrieved from the source,
makes a protein Unsafe in the catalog. A confirmed finding that falls short of that, a family rule by
name or the nose and skin reference is recorded with its tier. When that home is a tissue no panel
resolves and no reproductive evidence qualifies, the catalog makes the protein Other:Tissue; otherwise
only the filters act on the tier.</p>
<h3>Made Unsafe by the candidate literature sweep</h3>
{table(["Protein", "Tier", "Tissue", "Basis", "Sources"], finding_rows)}
<h3>Recorded in a Vital or Protected tissue, still a candidate</h3>
{table(["Protein", "Tier", "Tissue", "Why it does not meet the rubric"], recorded_rows)}
<h3>What the curation filter layer removes</h3>
{table(["Recorded tier", "Candidates", "Tissues", "Removed in"], tier_rows, numeric=(1,))}
<h3>Family rules</h3>
<p>A family member without reproductive evidence is Other:Tissue in the catalog when the family's home is a tissue no panel resolves. The counts are the members that remain candidates.</p>
{table(["Family", "Candidates", "Tissue", "Tier", "Sources"], family_rows, numeric=(1,))}
<h3>Tissues added to the tissue-tier registry</h3>
<p>Written from abstracts on 2026-11-12. Olfactory epithelium moved to Manageable on 2026-11-15 on the
oncology precedent: anosmia is graded at the lowest severity and accepted in standard cancer treatment
(<a href="../../tissue-tiers.md">tissue tiers</a>).</p>
{table(["Tissue", "Tier", "Confidence", "References"], proposed_rows)}
<h3>Read-through loci</h3>
<p>{len(loci)} catalog loci have every coding transcript joining two separate genes. The catalog's sequence
check places {len(moved)} of them in Other:Sequence; the rest are Unsafe or Other for another reason.
Other:Sequence: {esc(", ".join(sorted(moved.protein_group_label)))}.</p>

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
{table(["Protein", "Category", "Off-target TPM", "Tumor fold", "CTA papers", "Curated tier", "loose", "standard", "strict", "Unsafe reasons"], known_rows, numeric=(2, 3, 4))}
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
