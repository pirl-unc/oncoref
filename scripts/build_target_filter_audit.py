"""Rebuild docs/audits/target-filters-20261112: the filter report and its tables.

Every number in ``index.md`` is computed here from the packaged catalog.
Usage: PYTHONPATH=. python scripts/build_target_filter_audit.py
"""

from dataclasses import replace
from pathlib import Path

import pandas as pd

from oncoref import target_filters as tf
from oncoref import target_selection
from oncoref.cta_literature import sources

OUT = Path(__file__).resolve().parents[1] / "docs/audits/target-filters-20261112"
SETS = ("CTA", "Cancer Antigen only", "Reproductive Antigen only")
VIEWS = ("loose", "standard", "strict")

#: Familiar published CTAs, by group. Names are catalog labels or their first member.
FAMILIAR = {
    "MAGE-A": "MAGEA1 MAGEA2 MAGEA3 MAGEA4 MAGEA6 MAGEA10 MAGEA12",
    "MAGE-B/C": "MAGEB1 MAGEB2 MAGEB4 MAGEC1 MAGEC2",
    "NY-ESO / LAGE / PRAME": "NY-ESO-1 CTAG2 PRAME LAGE3",
    "SSX / GAGE / PAGE / XAGE": "SSX1 SSX2 SSX4 GAGE1 GAGE2A PAGE1 PAGE4 XAGE1A/B",
    "Testis / meiosis": "SYCP1 TEX14 AKAP4 ODF1 ADAM2 TDRD1 HORMAD1 MAEL PIWIL2 DDX43",
    "Others": "CTCFL SPA17 SPAG9 PBK CEP55 TPTE CT45A1 CT83 PLAC1 CSAG1 LY6K NUTM1",
}
SWEEPS = {
    "max_off_target_tpm": [5.0, 2.0, 1.0, 0.5, 0.1],
    "max_off_target_ntpm": [5.0, 2.0, 1.0, 0.5, 0.1],
    "max_disposable": [10.0, 1.0],
    "rare_cell": ["no_dismissed"],
    "literature": ["no_conflict"],
    "mapping": ["no_risk"],
    "normal_function": ["no_vital", "no_protected", "no_manageable"],
    "readthrough": ["exclude"],
    "min_patients": [2, 5, 10, 50],
    "min_tumor_fold": [10.0, 30.0, 100.0],
    "min_reproductive": [5.0, 10.0],
    "membership": ["cancer_antigen", "reproductive_antigen", "cta"],
    "min_papers": [1, 2, 3, 4],
    "min_families": [2],
    "exact_members_only": [True],
}
CATALOG_RULES = {
    "CTA": "Passes the shared screen; cancer evidence and reproductive evidence both qualify",
    "Cancer Antigen only": "Passes the shared screen; only cancer evidence qualifies",
    "Reproductive Antigen only": "Passes the shared screen; only reproductive evidence qualifies",
    "Other": "Not excluded, but fails coverage, identity, expression bounds, sequence or "
    "both evidence routes",
    "Unsafe": "An accepted normal-tissue or clinical exclusion",
}
LAYER_ROWS = [
    (
        "1 `off_target`",
        "How low must expression be outside tumors and reproductive tissue?",
        "`max_off_target_tpm`, `max_off_target_ntpm`, `max_disposable`",
    ),
    (
        "2 `warnings`",
        "Should flags the catalog records, but does not exclude on, remove a protein?",
        "`rare_cell`, `literature`, `mapping`",
    ),
    (
        "3 `curation`",
        "Does cited literature place the protein in a normal tissue, or is it a read-through?",
        "`normal_function`, `readthrough`",
    ),
    (
        "4 `sets`",
        "How strong must each set's evidence be, and which set do you want?",
        "`min_patients`, `min_tumor_fold`, `min_reproductive`, `membership`",
    ),
    (
        "5 `papers`",
        "How much published support?",
        "`min_papers`, `families`, `min_families`, `exact_members_only`",
    ),
    (
        "6 `exclusions`",
        "Which gene families do you rule out by name?",
        "`exclude_prefixes`, `except_labels`",
    ),
]


def counts(d):
    kept = d[d.keep]
    by = kept.filtered_category.value_counts()
    return [len(kept), *(int(by.get(k, 0)) for k in SETS)]


def md(frame):
    def cell(v):
        if isinstance(v, float):
            return "" if pd.isna(v) else f"{v:,.1f}"
        return str(v)

    head = "| " + " | ".join(frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    body = ["| " + " | ".join(cell(v) for v in row) + " |" for row in frame.to_numpy()]
    return "\n".join([head, rule, *body])


def lookup(d, name):
    hit = d[d.protein_group_label.map(lambda x: name in (x, x.split("/")[0]))]
    return None if hit.empty else hit.iloc[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    table = tf.measures(target_selection.results())
    cand = table[table.is_target_candidate]
    cols = ["Kept", *SETS]

    applied = {name: tf.apply(f, table) for name, f in tf.FILTER_MODES.items()}
    verdict = {
        name: d.set_index("entity_id").filter_failed.replace("", "✓") for name, d in applied.items()
    }

    modes = pd.DataFrame(
        [(f"`{n}`", *counts(d)) for n, d in applied.items()], columns=["Filter mode", *cols]
    )
    modes.assign(rule=[tf.describe(f).replace("\n", " ") for f in tf.FILTER_MODES.values()]).to_csv(
        OUT / "filter-modes.csv", index=False
    )

    rows = []
    for knob, values in SWEEPS.items():
        for v in values:
            rows.append((f"`{knob}`", v, *counts(tf.apply(tf.TargetFilter(**{knob: v}), table))))
    sweep = pd.DataFrame(rows, columns=["Knob", "Value", *cols])
    sweep.to_csv(OUT / "knob-sweep.csv", index=False)

    keep = [
        "protein_group_label",
        "category",
        "off_target_tpm",
        "off_target_tpm_tissue",
        "off_target_ntpm",
        "disposable_max_tpm",
        "disposable_max_ntpm",
        "tumor_max_tpm",
        "independent_patients_above_normal",
        "reproductive_best",
        "papers_cta",
        "papers_oncofetal",
        "papers_tumor",
        "rare_cell_status",
        "mapping_risk",
    ]
    disposable = cand[cand[["disposable_max_tpm", "disposable_max_ntpm"]].max(axis=1).ge(1)]
    disposable = disposable.sort_values("disposable_max_tpm", ascending=False)
    disposable[keep].to_csv(OUT / "disposable-at-least-1.csv", index=False)
    low = cand[cand.off_target_tpm.lt(1) & cand.off_target_ntpm.lt(1)]
    low[keep].to_csv(OUT / "off-target-below-1.csv", index=False)

    papers = tf.paper_table()
    papers = papers[papers.entity_id.isin(cand.entity_id)]
    names = sources()[["source_id", "citation", "kind", "scope"]].drop_duplicates("source_id")
    paper_rows = (
        papers.groupby(["family", "source_id"])
        .size()
        .rename("candidates")
        .reset_index()
        .merge(names, on="source_id", how="left")
        .sort_values(["family", "candidates"], ascending=[True, False])
    )
    paper_rows.to_csv(OUT / "paper-families.csv", index=False)

    rows = []
    for group, members in FAMILIAR.items():
        for name in members.split():
            r = lookup(table, name)
            if r is None:
                rows.append({"Group": group, "Name": name, "Catalog": "not in catalog"})
                continue
            candidate = bool(r.is_target_candidate)
            reasons = r.Unsafe_reasons if isinstance(r.Unsafe_reasons, str) else ""
            row = {
                "Group": group,
                "Name": name,
                "Catalog": r.category,
                "Off-target TPM": r.off_target_tpm,
                "in": r.off_target_tpm_tissue if candidate else "",
                "Off-target nTPM": r.healthy_maximum_ntpm if not candidate else r.off_target_ntpm,
                "Tumors above normal": (
                    ""
                    if pd.isna(r.independent_patients_above_normal)
                    else int(r.independent_patients_above_normal)
                ),
                "CTA papers": int(r.papers_cta),
                "Flags": ", ".join(
                    flag
                    for flag, on in (
                        ("dismissed rare-cell", r.rare_cell_status == "dismissed"),
                        ("mapping risk", bool(r.mapping_risk)),
                        ("literature conflict", r.literature_verdict in tf.LITERATURE_CONFLICTS),
                    )
                    if on
                ),
            }
            for view in VIEWS:
                row[f"`{view}`"] = verdict[view][r.entity_id] if candidate else "—"
            row["Unsafe reasons"] = reasons
            rows.append(row)
    familiar = pd.DataFrame(rows).astype(object).fillna("")
    familiar.to_csv(OUT / "familiar-ctas.csv", index=False)

    mage = cand[cand.protein_group_label.str.startswith("MAGE")].sort_values("protein_group_label")
    unnamed = {
        view: tf.apply(replace(tf.FILTER_MODES[view], exclude_prefixes=(), except_labels=()), table)
        .set_index("entity_id")
        .filter_failed.replace("", "✓")
        for view in ("standard", "strict")
    }
    mage_rows = pd.DataFrame(
        {
            "MAGE": mage.protein_group_label,
            "Catalog": mage.category,
            "Off-target TPM": mage.off_target_tpm,
            "Off-target nTPM": mage.off_target_ntpm,
            "Tumors above normal": mage.independent_patients_above_normal.astype(int),
            "CTA papers": mage.papers_cta,
            "Mapping risk": mage.mapping_risk.map({True: "yes", False: ""}),
            "`standard`": mage.entity_id.map(unnamed["standard"]),
            "`strict` on its own data": mage.entity_id.map(unnamed["strict"]),
        }
    )
    mage_rows.to_csv(OUT / "mage-family.csv", index=False)

    from oncoref import target_curation

    loci = target_curation.readthrough()
    comps = target_curation.compartments()
    tier_rows = pd.DataFrame(
        [
            (
                tier.capitalize(),
                int(cand.curated_tier.eq(tier).sum()),
                ", ".join(
                    comps.set_index("compartment").label[
                        sorted(
                            {
                                k
                                for text in cand[cand.curated_tier.eq(tier)].curated_compartments
                                for k in text.split("; ")
                                if comps.set_index("compartment").tier[k] == tier
                            }
                        )
                    ]
                ),
                knob,
            )
            for tier, knob in zip(
                target_curation.TIERS,
                (
                    "`no_vital` and stricter",
                    "`no_protected` and stricter",
                    "`no_manageable`",
                    "never",
                ),
            )
        ],
        columns=["Curated tier", "Candidates", "Tissues", "Removed by `normal_function`"],
    )
    proposed = comps[comps.tier_source == "proposed"]
    proposed_rows = pd.DataFrame(
        {
            "Compartment": proposed.label,
            "Proposed tier": proposed.tier.str.capitalize(),
            "Confidence": proposed.confidence,
            "References": proposed.tier_references,
        }
    )
    catalog_counts = table.category.value_counts()
    catalog_rows = pd.DataFrame(
        [(k, f"{int(catalog_counts[k]):,}", v) for k, v in CATALOG_RULES.items()],
        columns=["Category", "Proteins", "Rule"],
    )
    layer_alone = []
    for (layer, question, knobs), view in zip(LAYER_ROWS, tf.LAYERS):
        strict = applied["strict"]
        removed = int((strict.is_target_candidate & ~strict[f"pass_{view}"]).sum())
        layer_alone.append((layer, question, knobs, removed))
    layers = pd.DataFrame(
        layer_alone, columns=["Layer", "Question", "Knobs", "Removed by `strict`, alone"]
    )

    funnels = []
    for name in tf.FILTER_MODES:
        if name == "all":
            continue
        steps = tf.funnel(name, table)
        chain = " → ".join(f"{r.layer} **{r.remaining}**" for r in steps.itertuples())
        rule = "\n".join(f"- {line}" for line in tf.describe(name).splitlines()[1:])
        funnels.append(f"### `{name}`\n\n{rule}\n\n{chain}")

    prame, ctag2 = lookup(table, "PRAME"), lookup(table, "CTAG2")
    both_low = int((cand.off_target_tpm.lt(1) & ~cand.off_target_ntpm.lt(1)).sum())
    nominated = int(cand.papers_cta.gt(0).sum())
    pass_all = familiar[(familiar[[f"`{v}`" for v in VIEWS]] == "✓").all(axis=1)].Name
    is_candidate = familiar.Catalog.isin(SETS)
    by_catalog = familiar.Catalog.value_counts()
    nominating = paper_rows[paper_rows.family == "cta"]
    top3 = nominating.head(3)
    in_top3 = papers[papers.source_id.isin(top3.source_id)].groupby("entity_id").size()
    table3 = tf.apply(tf.TargetFilter(min_papers=3), table)
    three = table3[table3.keep]
    only_big = int(three.entity_id.map(in_top3).eq(3).sum())
    dropped = set(mage.protein_group_label) - {"MAGEA4"}
    data_pass = mage_rows[mage_rows["`strict` on its own data"] == "✓"].MAGE
    data_pass = data_pass[data_pass != "MAGEA4"]

    report = f"""# Target filters: categories, knobs and where proteins land

**The catalog assigns five categories and they are not changed here. Filters only
remove proteins from the {len(cand)} Target candidates, through six layers that
follow the catalog's own shape.** Filter policy `{tf.FILTER_POLICY_ID}`, applied to
catalog `oncoref.target-selection.2026-11-10.1`. Every number below is computed by
`scripts/build_target_filter_audit.py`. [Rules reference](../../target-filters.md).

These are outputs of screening rules. They are not validated tumor antigens.

## 1. The categories

{md(catalog_rows)}

Target candidates are the first three rows: {len(cand)} proteins. A filter never moves
a protein into or out of Other or Unsafe.

## 2. The six filter layers

A protein is kept only if it passes all six. `filter_failed` names the first one it
fails. `max_*` knobs mean **below** the value, `min_*` knobs mean **at least** the
value, and a missing measurement never passes.

{md(layers)}

**Layer 4 works on the two sets separately.** `min_patients` and `min_tumor_fold`
raise the bar for being a Cancer Antigen. `min_reproductive` raises it for being a
Reproductive Antigen. A protein that loses one membership keeps the other, so a CTA
with weak tumor evidence becomes Reproductive Antigen only and stays in the list.
`membership` then picks the set: `membership="cta"` requires both.

## 3. Each knob alone

Target candidates kept with only that knob changed. Categories are the ones after
layer 4 has tightened the evidence.

{md(sweep)}

## 4. Filter modes

`all` → `loose` → `standard` → `strict`: each mode keeps a subset of the one before.
Every mode keeps the positive controls {", ".join(tf.POSITIVE_CONTROLS)}; a test fails
if a change to a mode drops one.

{md(modes)}

{chr(10).join(chr(10) + block for block in funnels).strip()}

## 5. Where the familiar CTAs land

{len(familiar)} well-known names. In the catalog, {int(is_candidate.sum())} are Target
candidates, {int(by_catalog.get("Unsafe", 0))} are Unsafe and
{int(by_catalog.get("Other", 0))} Other. A mode column shows ✓ or the first layer
that removes the protein. **{len(pass_all)} pass every mode:** {", ".join(pass_all)}.

{md(familiar[is_candidate].drop(columns=["Group", "Unsafe reasons"]))}

Not candidates in the catalog:

{md(familiar[~is_candidate][["Name", "Catalog", "Off-target TPM", "Off-target nTPM", "Unsafe reasons"]])}

## 6. The MAGE family

`strict` drops every label starting with `MAGE` except MAGEA4:
{len(dropped)} of {len(mage)} MAGE candidates. **This is a judgment by family name, not
a measurement.** The last column shows what `strict` would decide from each protein's
own data with the name rule switched off: {len(data_pass)} of the {len(dropped)} it drops
would pass ({", ".join(data_pass)}).

{md(mage_rows)}

The data-driven alternative is `mapping="no_risk"`, which drops proteins whose reads
could belong to a paralog or a host gene. It removes the MAGEs marked above and keeps
MAGEA4, but it also removes {int(cand.mapping_risk.sum())} of the {len(cand)}
candidates, including most SSX, GAGE and CT45 proteins.

## 7. Disposable tissue and the "below 1 everywhere" list

- **Disposable.** {len(disposable)} of {len(cand)} candidates reach 1 TPM or 1 nTPM in
  appendix, gallbladder or tonsil. The highest is
  {disposable.iloc[0].protein_group_label} at {disposable.iloc[0].disposable_max_tpm:.1f} TPM.
  [List](disposable-at-least-1.csv).
- **Below 1 everywhere off target.** `max_off_target_tpm=1.0` with
  `max_off_target_ntpm=1.0` keeps {len(low)} candidates.
  [List](off-target-below-1.csv). It is a filter, not a catalog category.

## 8. The curation layer

Layer 3 acts on cited literature. A finding acts only when a second reader confirmed
it as human, adult, normal-tissue and gene-specific. Its tissue is tiered by the
hazard registry, or by a cited proposal where the registry has no call.
[Findings and method](normal-function-sweep.md) ·
[every finding with its sources](normal-function-sweep.csv).

{md(tier_rows)}

Compartments the registry did not tier, with the tier proposed here:

{md(proposed_rows)}

`readthrough="exclude"` removes {len(loci)} loci that share exons with two or more other
protein-coding genes, including ATF7-NPFF, ENSG00000273155 and ENSG00000267618:
{", ".join(loci.protein_group_label)}.

## 9. What to keep in mind

- **The paper layer does most of the cutting.** {nominated} of {len(cand)} candidates
  are nominated by any CTA-list paper. `min_papers=3` keeps {len(three)}, and
  {only_big} of those are named by all of the three largest lists
  ({", ".join(top3.citation)}).
- **Paper counts are not independent replication.** The nominating sources are
  published catalogs and screens, and their lists overlap. Several lists naming a gene
  shows that several screens picked it, not that separate experiments confirmed it.
  [Sources and recorded scopes](paper-families.csv).
- **No mode reads the nTPM knob or the warning knobs.** PRAME is
  {prame.off_target_tpm:.1f} TPM at tissue level ({prame.off_target_tpm_tissue}) but
  {prame.off_target_ntpm:.1f} nTPM on one single-cell cluster, and carries a dismissed
  rare-cell signal and a conflicting literature verdict. CTAG2 is
  {ctag2.off_target_tpm:.1f} TPM ({ctag2.off_target_tpm_tissue}) and
  {ctag2.off_target_ntpm:.1f} nTPM. A mode using `max_off_target_ntpm`,
  `rare_cell="no_dismissed"` or `literature="no_conflict"` would drop PRAME, so those
  stay knobs you set yourself.
- **{both_low} candidates are below 1 TPM at tissue level but not below 1 nTPM.**
  `max_off_target_tpm=1.0` with `max_off_target_ntpm=1.0` is the "below 1 everywhere"
  filter; it removes PRAME by construction.
"""
    (OUT / "index.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
