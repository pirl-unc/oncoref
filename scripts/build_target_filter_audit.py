"""Rebuild docs/audits/target-filters-20261112: the filter report and its tables.

Every number in ``index.md`` is computed here from the packaged catalog.
Usage: PYTHONPATH=. python scripts/build_target_filter_audit.py
"""

import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import target_filters as tf  # noqa: E402
from oncoref import target_selection  # noqa: E402
from oncoref.cta_literature import sources  # noqa: E402

OUT = ROOT / "docs/audits/target-filters-20261112"
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
    "normal_function": ["no_vital", "no_protected", "no_manageable"],
    "min_tumor_fold": [20.0, 30.0, 50.0, 100.0],
    "min_reproductive_fold": [20.0, 30.0, 50.0, 100.0],
    "min_tumor_patients": [1, 2, 3, 5],
    "membership": ["cancer_antigen", "reproductive_antigen", "cta"],
    "min_papers": [1, 2, 3, 4],
}
#: Literature verdicts reported as a flag in the familiar-CTA table.
LITERATURE_CONFLICTS = ("conflicting_or_nonspecific", "specialized_tissue_expression")
CATALOG_RULES = {
    "CTA": "Passes the shared screen; cancer evidence and reproductive evidence both qualify",
    "Cancer Antigen only": "Passes the shared screen; only cancer evidence qualifies",
    "Reproductive Antigen only": "Passes the shared screen; only reproductive evidence qualifies",
    "Other": "Not excluded, but fails coverage, identity, expression bounds, sequence, the "
    "home-tissue rule or both evidence routes",
    "Unsafe": "An accepted normal-tissue or clinical exclusion",
}
LAYER_ROWS = [
    (
        "1 `off_target`",
        "How low must expression be outside tumors and reproductive tissue?",
        "`max_off_target_tpm`",
    ),
    (
        "2 `curation`",
        "Is the protein's recorded normal home in a tissue of this tier or above?",
        "`normal_function`",
    ),
    (
        "3 `host`",
        "Do protected compartments have sufficient coverage and bounded exposure?",
        "`host` (HostContext), or explicit profile/tissue choices",
    ),
    (
        "4 `sets`",
        "How much biological enrichment/replication, and which set?",
        "`min_tumor_fold`, `min_reproductive_fold`, `min_tumor_patients`, `membership`",
    ),
    (
        "5 `papers`",
        "How many papers nominate it or report it in tumors?",
        "`min_papers`",
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
    # Reported, not filtered on: HPA nTPM (single-cell clusters included) and mapping flags.
    table["off_target_ntpm"] = table[["healthy_maximum_ntpm", "disposable_max_ntpm"]].max(
        axis=1, skipna=False
    )
    table["mapping_risk"] = table.paralog_ambiguous.eq(True) | table.host_overlapped.eq(True)
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

    profiles = tf.profile_summary(table)
    profiles.to_csv(OUT / "host-profiles.csv", index=False)

    keep = [
        "protein_group_label",
        "category",
        "off_target_tpm",
        "off_target_tpm_tissue",
        "off_target_ntpm",
        "disposable_max_tpm",
        "disposable_max_ntpm",
        "tumor_max_tpm",
        "tumor_fold",
        "reproductive_fold",
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "anonymous_samples_at_tenfold",
        "cancer_types_at_tenfold",
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
                "Tumor fold": "" if pd.isna(r.tumor_fold) else round(float(r.tumor_fold)),
                "CTA papers": int(r.papers_cta),
                "Flags": ", ".join(
                    flag
                    for flag, on in (
                        ("dismissed rare-cell", r.rare_cell_status == "dismissed"),
                        ("mapping risk", bool(r.mapping_risk)),
                        ("literature conflict", r.literature_verdict in LITERATURE_CONFLICTS),
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
            "Tumor fold": mage.tumor_fold.round(),
            "CTA papers": mage.papers_cta,
            "Mapping risk": mage.mapping_risk.map({True: "yes", False: ""}),
            "`standard`": mage.entity_id.map(unnamed["standard"]),
            "`strict` on its own data": mage.entity_id.map(unnamed["strict"]),
        }
    )
    mage_rows.to_csv(OUT / "mage-family.csv", index=False)

    from oncoref import target_curation, tissue_tiers

    loci = target_curation.readthrough()
    loci = loci.drop_duplicates("entity_id").merge(
        table[["entity_id", "category", "classification"]], on="entity_id"
    )
    moved = loci[loci.classification.eq("Other:Sequence")]
    document = target_curation.literature()
    entries = target_curation.entries()
    excluded = entries[entries.effect.eq("exclude") & entries.review.eq("candidate_sweep")]
    excluded_rows = pd.DataFrame(
        {
            "Protein": excluded.protein,
            "Tissue": excluded.decision_tissue,
            "Tier": excluded.tissue_tier.str.capitalize(),
            "Basis": excluded.decision_reason,
        }
    )
    removal = ("`no_vital` and stricter", "`no_protected` and stricter", "`no_manageable`", "never")
    tier_rows = pd.DataFrame(
        [
            (
                tier.capitalize(),
                int(cand.curated_tier.eq(tier).sum()),
                ", ".join(
                    sorted(
                        {
                            text.split(": ")[0]
                            for text in cand[cand.curated_tier.eq(tier)].curated_compartments
                            for text in text.split("; ")
                        }
                    )
                ),
                knob,
            )
            for tier, knob in zip(target_curation.TIERS, removal)
        ],
        columns=["Recorded tier", "Candidates", "Tissues", "Removed by `normal_function`"],
    )
    added = target_curation.compartments()
    references = tissue_tiers.references().set_index("reference_id").url
    proposed_rows = pd.DataFrame(
        {
            "Tissue": added.label,
            "Tier": added.tier.str.capitalize(),
            "Confidence": added.confidence,
            "References": added.references.map(
                lambda keys: ", ".join(f"[{k}]({references[k]})" for k in keys)
            ),
        }
    )
    families = target_curation.families()
    family_rows = pd.DataFrame(
        {
            "Family": families.label,
            "Candidates": families.id.map(cand.literature_family.value_counts())
            .fillna(0)
            .astype(int),
            "Tissue": families.tissue,
            "Tier": families.tier.str.capitalize(),
            "Used by filters": families.filter_acts.map({True: "yes", False: "no"}),
            "Sources": families.sources,
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
    antigens = cand[cand.is_cancer_antigen]
    fold_bands = pd.cut(
        antigens.tumor_fold, [10, 20, 50, 100, 1000, float("inf")], right=False
    ).value_counts(sort=False)
    fold_text = ", ".join(
        f"{int(n)} at {int(b.left)}"
        + (f" to {int(b.right)}x" if b.right != float("inf") else "x or more")
        for b, n in fold_bands.items()
    )
    top_samples = antigens.tumor_top_sample.value_counts()
    other_tissue = table[table.classification.eq("Other:Unscreened tissue")]

    report = f"""# Target filters: categories, knobs and where proteins land

**The catalog assigns five categories and they are not changed here. Filters only
remove proteins from the {len(cand)} Target candidates, through six layers that
follow the catalog's own shape.** Filter policy `{tf.FILTER_POLICY_ID}`, applied to
catalog `{target_selection.POLICY_ID}`. Every number below is computed by
`scripts/build_target_filter_audit.py`. [Rules reference](../../target-filters.md).

These are outputs of screening rules. They are not validated tumor antigens.

## 1. The categories

{md(catalog_rows)}

Target candidates are the first three rows: {len(cand)} proteins. A filter never moves
a protein into or out of Other or Unsafe.

**One cancer rule.** A protein is a Cancer Antigen when one tumor sample, from any
QC-passing comparable TPM cohort including germ cell tumors, is at least 10 times its highest
non-reproductive normal value (normal tissue, immune cells and Disposable tissue; values
below 1 TPM count as 1). The {len(antigens)} Cancer Antigens: {fold_text}.
{int(antigens.samples_at_tenfold.eq(1).sum())} have exactly one qualifying sample;
{int(antigens.independent_patients_at_tenfold.ge(2).sum())} have qualifying samples from
at least two known patients. Anonymous samples are reported separately.
RA uses the same 10× denominator in reproductive tissue or a supported embryo-lineage mean.
{len(other_tissue)} proteins whose normal home is a tissue no panel resolves are
Other:Unscreened tissue unless reproductive evidence qualifies.

## 2. The six filter layers

A protein is kept only if it passes all six. `filter_failed` names the first one it
fails. `max_*` knobs mean **below** the value, `min_*` knobs mean **at least** the
value, and a missing measurement never passes.

{md(layers)}

**Layer 3 tightens the two sets independently.** Both fold knobs raise the catalog's 10x;
a protein that loses Cancer Antigen membership keeps its Reproductive Antigen
membership, so a CTA with a weaker tumor becomes Reproductive Antigen only and stays in
the list. `membership` then picks the set: `membership="cta"` requires both.

## 3. Each knob alone

Target candidates kept with only that knob changed. Categories are the ones after
layer 3 has tightened Cancer Antigens.

{md(sweep)}

## 4. Filter modes and host toxicity profiles

`all` → `loose` → `standard` → `strict`: each mode keeps a subset of the one before.
The reference proteins {", ".join(tf.POSITIVE_CONTROLS)} pass through the same thresholds as every other protein.

{md(modes)}

{chr(10).join(chr(10) + block for block in funnels).strip()}

### Host anatomy and tissue-tolerance profiles

Host protection is independent of biological membership. Presets assume anatomy and
protect no reproductive tissue by default. Explicit `protect`, `accept_loss` and `absent`
choices have different meanings.
Missing coverage and unresolved protein evidence fail the screen. Native TPM and nTPM
are screened separately; protected TPM also enters every off-target limit and cancer fold.
[Precise contract](../../host-selection.md). [CSV](host-profiles.csv).

{md(profiles[["profile", "protected_compartments", "CTA", "CA_only", "RA_only", "Union"]].rename(columns={"profile": "Profile", "protected_compartments": "Protected compartments", "CA_only": "CA only", "RA_only": "RA only", "Union": "CA | RA"}))}

## 5. Where the familiar CTAs land

{len(familiar)} well-known names. In the catalog, {int(is_candidate.sum())} are Target
candidates, {int(by_catalog.get("Unsafe", 0))} are Unsafe and
{int(by_catalog.get("Other", 0))} Other. A mode column shows ✓ or the first layer
that removes the protein. **{len(pass_all)} pass every mode:** {", ".join(pass_all)}.
Flags and the nTPM column are reported for reading; no knob acts on them.

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

## 7. Disposable tissue and the "below 1 everywhere" list

- **Disposable.** {len(disposable)} of {len(cand)} candidates reach 1 TPM or 1 nTPM in
  appendix, gallbladder or tonsil. The highest is
  {disposable.iloc[0].protein_group_label} at {disposable.iloc[0].disposable_max_tpm:.1f} TPM.
  [List](disposable-at-least-1.csv).
- **Below 1 everywhere off target.** {len(low)} candidates are below 1 TPM at tissue
  level and below 1 HPA nTPM, single-cell clusters included.
  [List](off-target-below-1.csv). It is a list, not a knob or a catalog category.

## 8. Curation: what the catalog decides and what the filters add

**The catalog and the filters read the same records.** The literature review
(`normal-expression-literature.json`) holds one entry for each of {len(entries)}
proteins. The tissue-tier registry tiers the tissue. Nothing below lives in a separate
filter file. [Sweep findings and method](normal-function-sweep.md).

The catalog acts on the review's rubric. A gene-specific report in normal adult tissue
tiered Vital or Protected makes the protein Unsafe, and needs a quote retrieved from
the source. The 2026-11-12 candidate sweep adds {len(excluded)} such proteins:

{md(excluded_rows)}

The catalog's sequence check moves read-through loci to Other. {len(loci)} loci in the
catalog have every coding transcript joining two separate genes. {len(moved)} are
Other:Sequence; the rest are Unsafe or Other for another reason. Other:Sequence:
{", ".join(sorted(moved.protein_group_label))}.

Layer 2 acts on the recorded normal home of the candidates that remain: a confirmed
finding that falls short of the rubric, a family rule, or the nose and skin reference.

{md(tier_rows)}

Family rules. A member whose family's home is a tissue no panel resolves, and that has
no reproductive evidence, is Other:Unscreened tissue in the catalog; the rows count the members
that remain candidates:

{md(family_rows)}

Tissues added to the tissue-tier registry for these findings
(`normal-tissue-tier-updates.yaml`). Olfactory epithelium is Manageable since
2026-11-15, on the oncology precedent ([tissue tiers](../../tissue-tiers.md#how-a-tier-is-chosen)):

{md(proposed_rows)}

Limits: {"; ".join(x[0].lower() + x[1:] for x in document["limitations"][-4:])}.

## 9. What to keep in mind

- **One tumor sample can admit a protein.** The rule asks for a single sample at 10x,
  as requested; it does not ask for replication. The most frequent top sample is
  `{top_samples.index[0]}`, the best tumor for {int(top_samples.iloc[0])} Cancer
  Antigens. `samples_at_tenfold`, `independent_patients_at_tenfold`,
  `anonymous_samples_at_tenfold` and `cancer_types_at_tenfold` record qualifying support.
  `min_tumor_patients` optionally requires known-patient replication at the selected
  `min_tumor_fold` (10× by default).
- **The paper layer does most of the cutting.** {nominated} of {len(cand)} candidates
  are nominated by any CTA-list paper. `min_papers=3` keeps {len(three)}, and
  {only_big} of those are named by all of the three largest lists
  ({", ".join(top3.citation)}).
- **Paper counts are not independent replication.** The nominating sources are
  published catalogs and screens, and their lists overlap. Several lists naming a gene
  shows that several screens picked it, not that separate experiments confirmed it.
  [Sources and recorded scopes](paper-families.csv).
"""
    (OUT / "index.md").write_text(report)
    write_guide(table)
    print(report)


def write_guide(table):
    """Render the public filter guide from the same catalog and filter definitions."""
    mode_rows = []
    for mode in tf.FILTER_MODES:
        kept = tf.apply(tf.FILTER_MODES[mode], table)
        mode_rows.append(
            f"| `{mode}` | {int(kept.keep.sum())} | `{tf.FILTER_MODES[mode].changed()}` |"
        )
    profiles = tf.profile_summary(table)
    profile_rows = [
        f"| `{r.profile}` | {r.protected_compartments} | **{r.CTA}** | {r.CA_only} | {r.RA_only} | {r.Union} |"
        for r in profiles.itertuples()
    ]
    from string import Template

    template = (ROOT / "scripts/cta_evidence_assets/target-filter-guide.md").read_text()
    (ROOT / "docs/target-filters.md").write_text(
        Template(template).substitute(
            catalog_policy=target_selection.POLICY_ID,
            filter_policy=tf.FILTER_POLICY_ID,
            profile_rows="\n".join(profile_rows),
            mode_rows="\n".join(mode_rows),
        )
    )


if __name__ == "__main__":
    main()
