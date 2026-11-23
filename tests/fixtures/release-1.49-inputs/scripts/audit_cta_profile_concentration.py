#!/usr/bin/env python3
"""Audit common-panel coverage, unassigned reasons and active profile diffuseness.

Uses pinned HPA v23 / GTEx v11 through their uniform cache readers. Reports the
selection replay; the 10% identity boundary is not a safety threshold.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from build_cta_alternative_evidence import anatomy  # noqa: E402
from cta_profile_data import profiles, summarize  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402
from render_cta_source_walkthrough import markdown  # noqa: E402

from oncoref import expression_diffuseness as metric  # noqa: E402
from oncoref import reference_data  # noqa: E402

CURRENT = ROOT / "docs/audits/cta-selection-20261028"
OUT = ROOT / "docs/audits/cta-unassigned-20261029"
BASELINE = OUT / "baseline-1.22.0.csv.gz"
ID = "protein_equivalence_group_id"
HPA_BRANCHES = metric.HPA_REPRODUCTIVE_BRANCHES


def reason_table(frame):
    """Count disjoint combinations; never combine missingness with measured failures."""
    labels = {
        "Reproductive / embryonic association": "Reproductive expression",
        "Required RNA coverage": "Required coverage",
        "Tumor sample ≥2 TPM/pTPM": "Cancer expression",
    }
    reasons = frame.exclusion_reasons.fillna("").map(
        lambda text: " + ".join(labels.get(s, s) for s in text.split(" | "))
    )
    return reasons.value_counts().rename_axis("Failed requirements").reset_index(name="Proteins")


def spotlight(facts, hp, gp, hp_panel, gp_panel, hp_reference, gp_reference):
    """Export every term of ZNF730's score and a common-anatomy sensitivity check."""
    entity = facts.loc[facts.protein_group_label.eq("ZNF730"), ID].item()
    tables, reproductive, organs = {}, {}, {}
    for name, matrix, panel, reference, unit in [
        ("HPA", hp, hp_panel, hp_reference, "nTPM"),
        ("GTEx", gp, gp_panel, gp_reference, "TPM"),
    ]:
        values = matrix.loc[entity]
        r = values[reference].max()
        p = panel.copy()
        p["RNA"] = p.compartment.map(values)
        p["Reference RNA"] = r
        p["Unit"] = unit
        p["Similarity (%)"] = 100 * p.RNA / (p.RNA + r)
        p["Weight"] = p.weight
        p["Contribution (pp)"] = p["Similarity (%)"] * p.weight
        tables[name] = p
        organs[name] = p.groupby("anatomy")["Similarity (%)"].mean()
        reproductive[name] = (
            values[reference]
            .rename(unit)
            .sort_values(ascending=False)
            .rename_axis("Tissue")
            .reset_index()
        )
        write_csv(p, OUT / f"znf730-{name.lower()}-calculation.csv")
    common = organs["HPA"].index.intersection(organs["GTEx"].index)
    summary = pd.DataFrame(
        [
            {
                "Source": name,
                "Normal contexts": len(table),
                "Anatomy units": table.anatomy.nunique(),
                "Reference": table["Reference RNA"].iloc[0],
                "Unit": table.Unit.iloc[0],
                "Full-panel score (%)": table["Contribution (pp)"].sum(),
                "Same 20 anatomy units (%)": organs[name][common].mean(),
            }
            for name, table in tables.items()
        ]
    )
    text = (
        """# ZNF730: every input to the diffuseness score

Protein identity: **ZNF730 / ENSG00000183850**. These are RNA summaries, not
protein measurements or a probability of toxicity. The HPA adult/placental identity
cutoff is 10%; ZNF730 remains Unsafe on independent normal-expression evidence.

For every row: **similarity = N / (N + R)**, then **contribution in percentage
points = 100 × weight × similarity**. Sum the last column to obtain the score.
R is the source's highest Reproductive value, **testis in both tables**.
Calculations use unrounded values; displayed values are rounded.

"""
        + markdown(summary)
        + """

**The panels differ.** HPA includes tonsil, lymph node, bone
marrow, retina and other tissues absent from this GTEx bulk panel. GTEx includes
blood, arteries and peripheral nerve absent from HPA's consensus panel. Restricting
both summaries to the same 20 anatomy units gives **6.36% HPA versus 1.71% GTEx**.
This is a sensitivity analysis, not a newly calibrated admission rule; corresponding
organs can still have different subdivisions and assay summaries.

Thus the 9.53% versus 1.65% difference reflects **both panel composition and source
measurements**. It is not a difference measured in an identical tissue panel.
The Reproductive references happen to be very similar; the normal values are not.

## HPA v23 consensus

[HPA v23 downloadable data](https://v23.proteinatlas.org/about/download) ·
[ZNF730 tissue page](https://v23.proteinatlas.org/ENSG00000183850-ZNF730/tissue) ·
[Full-precision calculation CSV](znf730-hpa-calculation.csv).

29 anatomy units each receive weight 1/29. Eight brain compartments share that
weight (1/232 each); colon/rectum and small intestine/duodenum share their respective
organ weights (1/58 each). Other rows have weight 1/29.

**Thymus is Permissive**: its 2.5 nTPM remains in the source profile, but it is
excluded from both the comparator and the Reproductive reference. The denominator
is 29 anatomy units, not 30 with a zero substituted for thymus.

### Reproductive reference candidates

"""
        + markdown(reproductive["HPA"])
        + "\n\n### Normal-tissue terms\n\n"
        + markdown(
            tables["HPA"][["compartment", "RNA", "Weight", "Similarity (%)", "Contribution (pp)"]]
        )
        + """

**Sum: 9.53156210%.** Excluding thymus reduces the prior 10.33996949% score.
Tonsil, lymph node and bone marrow still count. They are absent from this direct
GTEx bulk panel. Disposable tissues remain in the comparison.

## Direct GTEx v11 donor medians

[GTEx v11 data and release information](https://gtexportal.org/home/downloads/adult-gtex) ·
[Full-precision calculation CSV](znf730-gtex-calculation.csv).

23 anatomy units each receive weight 1/23, divided equally among that organ's
native subdivisions. The 12 brain regions share 1/23; spinal cord is a separate
unit. These are donor medians, not the highest individual-donor measurements.

### Reproductive reference candidates

"""
        + markdown(reproductive["GTEx"])
        + "\n\n### Normal-tissue terms\n\n"
        + markdown(
            tables["GTEx"][["compartment", "RNA", "Weight", "Similarity (%)", "Contribution (pp)"]]
        )
        + """

**Sum: 1.64599344%.** A zero here is a measured donor median; it does not mean every
donor or cell is negative. HPA also incorporates GTEx evidence, so the two sources
are not independent biological replication.

[Method, fixed panels and missingness rules](index.md).
"""
    )
    text = text.replace("9.53156210%", f"{tables['HPA']['Contribution (pp)'].sum():.8f}%")
    (OUT / "znf730.md").write_text(text)


def impact_report(facts, hs, hp):
    """Evaluate the chosen rule and nearby thresholds against its fixed baseline."""
    before = pd.read_csv(OUT / "before-diffuseness.csv.gz", low_memory=False)
    score = hs.set_index("entity_id")
    comparison = before[
        [
            ID,
            "protein_group_label",
            "category",
            "reproductive_association",
            "developmental_association",
            "reproductive_ntpm",
        ]
    ].copy()
    comparison["score"] = comparison[ID].map(score.diffuseness)
    comparison["complete"] = comparison[ID].map(score.reference_complete).fillna(False).astype(
        bool
    ) & np.isclose(comparison[ID].map(score.normal_coverage), 1)
    comparison["after"] = comparison[ID].map(facts.set_index(ID).category)
    changes = comparison[
        comparison.category.fillna("Unassigned").ne(comparison.after.fillna("Unassigned"))
    ].copy()
    write_csv(changes, OUT / "diffuseness-category-changes.csv")
    rows = []
    current_cta = comparison[comparison.category.eq("CTA")]
    for cutoff in [0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25]:
        adult = current_cta[
            current_cta.reproductive_association
            & current_cta.complete
            & current_cta.score.gt(cutoff)
        ]
        lost = adult[~adult.developmental_association]
        rows.append(
            {
                "Cutoff (%)": cutoff * 100,
                "Adult routes failing": len(adult),
                "CTA to Neither": len(lost),
                "Proteins moved": ", ".join(sorted(lost.protein_group_label)) or "None",
            }
        )
    sensitivity = pd.DataFrame(rows)
    write_csv(sensitivity, OUT / "diffuseness-threshold-sensitivity.csv")
    counts = (
        pd.DataFrame(
            {
                "Before": before.category.fillna("Unassigned").value_counts(),
                "After": facts.category.fillna("Unassigned").value_counts(),
            }
        )
        .reindex(["CTA", "Neither", "Unsafe", "Unassigned"])
        .rename_axis("Category")
        .reset_index()
    )
    named = facts[
        facts.protein_group_label.isin(
            [
                "MTCP1",
                "PRAME",
                "MAGEA4",
                "NY-ESO-1",
                "XAGE1A/B",
                "CTAG2",
                "CT83",
                "CT55",
                "FAM72D",
                "ZNF878",
                "ZNF730",
            ]
        )
    ][
        [
            ID,
            "protein_group_label",
            "category",
            "reproductive_ntpm",
            "developmental_mean_tpm",
            "adult_profile_status",
        ]
    ].copy()
    named["Diffuseness (%)"] = named[ID].map(score.diffuseness) * 100
    named = named.drop(columns=ID).rename(
        columns={
            "protein_group_label": "Protein",
            "category": "Category",
            "reproductive_ntpm": "HPA reference nTPM",
            "developmental_mean_tpm": "Embryonic TPM",
            "adult_profile_status": "Adult contrast",
        }
    )
    named["Adult contrast"] = named["Adult contrast"].replace(
        {
            "no_qualifying_adult_reference": "Weak adult reference",
            "unassessed": "Missing comparison",
            "concentrated": "Concentrated",
            "diffuse": "Diffuse",
        }
    )
    write_csv(named, OUT / "diffuseness-benchmarks.csv")
    adult = current_cta[current_cta.reproductive_association]
    missing = adult[~adult.complete]
    embryo_only = current_cta[
        ~current_cta.reproductive_association & current_cta.developmental_association
    ]
    mtcp = (
        hp.loc[facts.loc[facts.protein_group_label.eq("MTCP1"), ID].item()]
        .rename("nTPM")
        .rename_axis("Tissue")
        .reset_index()
    )
    mtcp = mtcp.merge(
        metric.hpa_bulk_roles()[["compartment", "role"]],
        left_on="Tissue",
        right_on="compartment",
        validate="one_to_one",
    ).drop(columns="compartment")
    write_csv(mtcp, OUT / "mtcp1-all-hpa-tissues.csv")
    text = f"""# Diffuseness: the chosen panel, rule and CTA impact

**Use HPA v23 consensus RNA and a 10% cutoff for adult/placental Reproductive concentration.**
The current local replay moves **MTCP1 from CTA to Neither**. Safety criteria and
Unsafe membership are unchanged. These are research identity criteria, not clinical
safety estimates. The boundary is applied to every protein, without named-target exceptions.

{markdown(counts)}

[Every changed assignment](diffuseness-category-changes.csv) ·
[All source-specific scores](protein-diffuseness.csv.gz) ·
[Complete tissue roles and weights](hpa-v23-context-roles.csv) ·
[Full curation report](../cta-source-walkthrough-20261028/index.md)

## The fixed panel

Use [HPA v23 consensus RNA](https://v23.proteinatlas.org/about/download) because it
covers placenta and the adult reproductive organs used by the positive-expression
rule. This avoids mixing a GTEx donor median with an HPA maximum or a single-cell
embryonic mean. HPA consensus is a curated cross-dataset summary, not a single donor
cohort; this score cannot resolve donor heterogeneity or antigen presentation.

There are **38 comparator contexts representing 29 anatomical units**, 11
Reproductive reference contexts and one excluded **Permissive** context: thymus.
Brain subdivisions share one organ weight, as do colon/rectum and small intestine/duodenum.
The full comparator comprises adipose, adrenal, appendix, bone marrow, brain,
colon, esophagus, gallbladder, heart, kidney, liver, lung, lymph node, pancreas,
parathyroid, pituitary, retina, salivary gland, skeletal muscle, skin, small intestine,
smooth muscle, spinal cord, spleen, stomach, thyroid, tongue, tonsil and urinary bladder.

Disposable tissues count for identity. Thymus counts in neither the comparator
nor the Reproductive reference; its RNA remains visible. This Permissive role is
separate from the safety tier and from the legacy scope alias named `permissive`.
Thymic expression of peripheral antigens has experimental support in
[Projection of an Immunological Self Shadow Within the Thymus by the Aire Protein](https://doi.org/10.1126/science.1075958).
It does not imply every thymic cell expresses every protein. Resolved thymic anatomy
uses the same profile role in human, dog and mouse; quantitative scores here are human-only.
Explicit embryonic/fetal stage retains Reproductive precedence.

## Continuous score and cutoff

For each comparator context, use **N / (N + R)**, where N is that context's nTPM
and R is the highest of the 11 Reproductive nTPM values in this same HPA matrix.
Average within each anatomical unit, then average the 29 unit scores.
Multiply by 100 to display a percentage. It is an average RNA similarity,
**not the percentage of tissues expressing the protein**.

There is no per-tissue on/off expression threshold. If all normals are tenfold
below R the score is 9.09%; if all equal R it is 50%. One arbitrarily high organ
contributes at most 1/29 = 3.45 percentage points. A broad profile accumulates
contributions across organs. Scores **above 10%** fail the adult/placental identity
route; equality passes. This custom boundary expresses our preference, not a
published biological or safety threshold. The [benchmark of tissue-specificity
metrics](https://doi.org/10.1093/bib/bbw008) motivates attention to tissue sampling;
our bounded Reproductive-anchored score is not its Tau statistic.

## Which route can this reject?

Apply the identity gate only with a **qualifying adult/placental reference of at
least 5 nTPM**, complete HPA comparator/reference measurements and exact same-source
protein-member coverage. Missing members are not zero and missing tissues are not
removed from the denominator. Scores for weaker anchors remain descriptive.

A supported embryonic route remains independent. An adult profile cannot establish
embryonic diffuseness without compatible developmental/adult measurements. Among the
523 proteins entering this change, **{len(embryo_only)} have embryonic-only qualifying
evidence**; their developmental concentration remains explicitly unassessed.
**{len(missing)} adult-positive proteins lack a complete HPA contrast**:
{", ".join(sorted(missing.protein_group_label))}. Their required GTEx core remains
complete, and absence of this additional contrast is disclosed rather than imputed
as either passing or failing. This limitation should be resolved with compatible
panels, not by borrowing a favorable ratio from a different source.

A failed identity route with no remaining qualifying route becomes **Neither**
when required coverage is complete and normal safety has no unresolved hold.
Independent safety exclusions still produce Unsafe. Missing safety evidence does
not become Neither. IHC confidence changes the safety assessment, not this RNA score.

## Why MTCP1 changes

MTCP1 has testis 5.1 nTPM, skin 4.1, marrow 2.8, pancreas and tonsil 2.5,
liver/retina/lymph node 2.4, and widespread values around 1-2 nTPM.
Its **24.53%** HPA score indicates diffuse expression despite passing the existing
normal-tissue safety gates. Direct GTEx independently yields a descriptive **16.82%**
score; the sources are not independent cohorts because HPA incorporates GTEx.
Its supported embryonic summary is 1.32 TPM, below the 5 TPM association floor.
[Every MTCP1 HPA measurement](mtcp1-all-hpa-tissues.csv).

## Established targets and cases needing more data

{markdown(named)}

PRAME, MAGEA4, NY-ESO-1, XAGE1A/B and CTAG2 remain CTA. FAM72D and ZNF878 retain
embryonic evidence; their weak adult anchors cannot settle developmental specificity.
ZNF730's HPA score falls from 10.34% to 9.53% after excluding thymus, but its independent
normal-expression exclusions keep it Unsafe. Diffuseness does not rescue unsafe proteins.

## Sensitivity around the chosen cutoff

{markdown(sensitivity)}

At 5%, CTAG2 and GDF9 lose the adult concentration route but retain supported
embryonic evidence; GOLGA8K and MTCP1 would move to Neither. Between 7.5% and 20%,
only MTCP1 moves. At 25%, it remains CTA. We choose 10% for an interpretable broad
contrast, not to force a requested number of removals. [Full sensitivity table](diffuseness-threshold-sensitivity.csv).

## GTEx and reproduction

[GTEx v11](https://gtexportal.org/home/downloads/adult-gtex) remains the mandatory
52-context coverage core and a separately reported donor-median comparison. It is
not an interchangeable score: it lacks placenta and has different reproductive
values for some familiar antigens. No score takes the minimum or maximum across sources.

Rebuild with `scripts/regenerate_cta_walkthrough.py`, followed by
`scripts/audit_cta_profile_concentration.py`. Native tables use the uniform cache
readers. The report manifest pins code, source checksums, before/after assignments
and every output. `oncoref.expression_diffuseness.assess()` exposes numeric terms;
`oncoref.cta_selection_policy.concentration_evidence()` applies the route policy.
"""
    (OUT / "diffuseness-impact.md").write_text(text)


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(CURRENT / "protein-group-classifications.csv.gz", low_memory=False)
    baseline = pd.read_csv(BASELINE, low_memory=False)
    hp, gp, mapping = profiles(facts)
    hp_panel = metric.hpa_bulk_panel()
    hp_reference = sorted(t for tissues in HPA_BRANCHES.values() for t in tissues)
    if set(hp_reference) != {t for tissues in HPA_BRANCHES.values() for t in tissues}:
        raise ValueError("HPA Reproductive branch mapping is incomplete")
    types = anatomy(pd.DataFrame({"tissue": sorted(set(mapping.values()))})).set_index("tissue")
    gp_reference = [c for c, t in mapping.items() if types.loc[t, "category"] == "reproductive"]
    brain = {
        "amygdala",
        "basal ganglia",
        "cerebellum",
        "cerebral cortex",
        "hippocampal formation",
        "hypothalamus",
        "midbrain",
    }
    gp_panel = pd.DataFrame(
        [
            {"compartment": c, "anatomy": "brain" if t in brain else t, "tissue": t}
            for c, t in mapping.items()
            if c not in gp_reference
        ]
    )
    gp_panel["weight"] = (
        1 / gp_panel.groupby("anatomy").compartment.transform("size") / gp_panel.anatomy.nunique()
    )
    hs = summarize(hp, hp_reference, hp_panel, "hpa_rna_consensus", "v23", "nTPM", "consensus")
    gs = summarize(gp, gp_reference, gp_panel, "gtex_tissue_rna", "v11", "TPM", "donor median")
    spotlight(facts, hp, gp, hp_panel, gp_panel, hp_reference, gp_reference)
    impact_report(facts, hs, hp)
    identity = facts.rename(columns={ID: "entity_id", "protein_group_label": "Protein"})
    scores = pd.concat([hs, gs], ignore_index=True).merge(
        identity[
            [
                "entity_id",
                "Protein",
                "category",
                "developmental_association",
                "adult_profile_status",
                "developmental_profile_status",
            ]
        ],
        on="entity_id",
        validate="many_to_one",
    )
    scores["selection_role"] = np.where(
        scores.source.eq("hpa_rna_consensus"),
        "adult_placental_identity_gate",
        "descriptive_comparison",
    )
    # A low adult anchor cannot adjudicate an embryonic/placental-only route.
    scores["developmental_contrast"] = "unassessed: no compatible developmental/adult panel"
    write_csv(scores, OUT / "protein-diffuseness.csv.gz")
    write_csv(hp.reset_index(), OUT / "hpa-v23-protein-profiles.csv.gz")
    write_csv(gp.reset_index(), OUT / "gtex-v11-protein-profiles.csv.gz")
    write_csv(hp_panel, OUT / "hpa-v23-normal-panel.csv")
    write_csv(metric.hpa_bulk_roles(), OUT / "hpa-v23-context-roles.csv")
    write_csv(gp_panel, OUT / "gtex-v11-normal-panel.csv")
    core = pd.DataFrame({"native_context": list(mapping), "tissue": list(mapping.values())})
    core["role"] = np.where(
        core.native_context.isin(gp_reference), "Reproductive reference", "Normal comparator"
    )
    write_csv(core, OUT / "gtex-v11-core-panel.csv")
    branches = pd.DataFrame(index=hp.index)
    for branch, tissues in HPA_BRANCHES.items():
        branches[branch + " nTPM"] = hp[tissues].max(axis=1)
    branches["Embryonic TPM"] = identity.set_index("entity_id").developmental_mean_tpm
    branches["Fetal"] = "unassessed"
    branches["Other extraembryonic"] = "unassessed"
    write_csv(branches.reset_index(), OUT / "reproductive-branches.csv.gz")
    current = facts.merge(
        baseline[[ID, "category", "assessment_status", "exclusion_reasons"]],
        on=ID,
        suffixes=("", "_1_22"),
        validate="one_to_one",
    )
    flags = current[current.category.isna() | current.category_1_22.isna()]
    write_csv(flags, OUT / "unassigned-protein-details.csv.gz")
    for status in ["incomplete_evidence", "outside_category_definitions"]:
        write_csv(
            reason_table(baseline[baseline.assessment_status.eq(status)]),
            OUT / (status + "-1.22.csv"),
        )
        write_csv(
            reason_table(facts[facts.assessment_status.eq(status)]), OUT / (status + "-current.csv")
        )
    complete = gs.set_index("entity_id").core_complete
    lost = current[current.category_1_22.eq("CTA") & ~current[ID].map(complete)].copy()
    write_csv(
        lost[
            [
                ID,
                "protein_group_label",
                "member_gene_ids",
                "category",
                "assessment_status",
                "exclusion_reasons",
            ]
        ],
        OUT / "gtex-coverage-losses.csv",
    )
    examples = [
        "ZNF730",
        "MTCP1",
        "PRAME",
        "MAGEA4",
        "NY-ESO-1",
        "XAGE1A/B",
        "CTAG2",
        "CT83",
        "FOXR2",
        "CT55",
        "TBXT",
    ]
    named = scores[scores.Protein.isin(examples)].copy()
    named["Score (%)"] = 100 * named.diffuseness
    named = named[
        ["Protein", "category", "source", "reference_value", "unit", "Score (%)", "review"]
    ]
    write_csv(named, OUT / "named-protein-scores.csv")
    old = baseline[baseline.category.isna()]
    now = facts[facts.category.isna()]
    counts = facts.category.fillna("Unassigned").value_counts().to_dict()
    report = f"""# Common coverage, Reproductive identity and unassigned proteins

**Current local replay: {counts["CTA"]} CTA proteins, {counts["Unsafe"]:,} Unsafe,
{counts["Neither"]} Neither and {counts["Unassigned"]:,} unassigned.** The maximum-versus-maximum
restriction rule has been removed. Complete GTEx core coverage is required. The
HPA diffuseness gate uses **10%** for adult/placental identity.
[Chosen method, changed proteins and cutoff sensitivity](diffuseness-impact.md).

[Current curation report](../cta-source-walkthrough-20261028/index.md) ·
[Paper titles and database versions](../cta-source-walkthrough-20261028/references.md) ·
[Every unassigned protein and its measurements](unassigned-protein-details.csv.gz) ·
[All protein scores](protein-diffuseness.csv.gz) · [Every ZNF730 calculation](znf730.md)

## One common GTEx core

The fixed [GTEx v11](https://gtexportal.org/home/downloads/adult-gtex) panel has **52 native
bulk contexts**: nine Reproductive and 43 other contexts. The latter represent 23
anatomical units after subdivisions share an organ's weight. Require a measured
value for **every member locus in every context**. A measured zero counts as coverage;
a missing or ambiguous identity does not. Read-through and retired identifiers need
identity repair, not substitution with a nearby gene.

**{int(complete.sum()):,} / {len(complete):,} proteins have complete GTEx coverage.**
In this pinned matrix, all other proteins lack a complete member profile throughout
the panel; requiring 50%, 80% or 100% therefore gives the same answer. Coverage
failure excludes a protein from CTA and yields unassigned/incomplete evidence unless
independent safety evidence already establishes Unsafe. It does not establish Neither.
GTEx does not cover placenta, fetal stages, every rare cell or every Disposable tissue.

[All 52 contexts](gtex-v11-core-panel.csv) · [All GTEx profiles](gtex-v11-protein-profiles.csv.gz) ·
[23 previously retained proteins missing GTEx](gtex-coverage-losses.csv).
Removing the old maximum veto alone would add TGIF2-RAB5IF and ATF7-NPFF; both also
lack GTEx core coverage, so remain unassigned. The required HPA tissue/immune panels
remain in force. Optional DICE/cornea coverage does not shrink or expand the score's
denominator; measured counterevidence from those sources still contributes to safety.

## What the numerator and denominator mean

For each non-Reproductive context, calculate:

```text
tissue similarity = N / (N + R)
diffuseness = sum(fixed tissue weight × tissue similarity)
```

**N** is that tissue's RNA. **R** is the highest measured Reproductive RNA within
the **same source, version, unit and summary**. The numerator is normal RNA; the
denominator is that normal RNA plus the Reproductive reference. The outer average
has fixed weights summing to one. It is **not “10% of tissues”**, a count of positive
tissues, a percentage of total body RNA, or a probability of harm.

| Uniform normal profile | Score |
|---|---:|
| Every normal tissue is zero, with a positive reference | 0% |
| Every normal tissue is 1/10 of the reference | 9.09% |
| Every normal tissue equals the reference | 50% |
| Every normal tissue greatly exceeds the reference | Approaches 100% |

The **adult/placental identity boundary is >10%**. In the uniform example this means normals above
1/9 of the reference. One extreme GTEx organ contributes at most 1/23 = 4.35%; one
HPA organ at most 1/29 = 3.45%. Thus a single outlier does not veto identity, while
many comparable tissues raise the score. These weights correct panel sampling;
they are not organ sizes, cell counts or safety priorities.

Missing normals retain their fixed weight: report a lower/upper bound, never a
zero-filled or smaller-denominator score. A reference below **2 native units** gets
a descriptive value but no concentrated/diffuse verdict. This guard is a proposed
research parameter, not an assay detection limit. Zero reference gives no score.

## Consequences of diffuseness

**A supported diffuse adult/placental profile with no other qualifying route and
clear normal safety belongs in Neither, not Unsafe.** The HPA gate requires its
existing ≥5 nTPM positive reference and a complete same-source comparison.
Independent safety exclusions retain Unsafe precedence. Required coverage failures
remain unassigned. Missing optional concentration measurements remain unassessed,
without becoming a pass/fail value. [Full policy and impact](diffuseness-impact.md).

An adult comparison cannot reject supported embryonic evidence. Compatible
embryonic/adult measurements are needed to resolve its concentration. GTEx is a
separate descriptive comparison, not a second veto or a fallback denominator.
IHC affects safety but does not mathematically reduce RNA diffuseness.

## Measurements in familiar proteins

{markdown(named)}

GTEx uses donor medians; HPA uses its consensus summary. These differ in sampling,
annotation and processing and are not independent replications (HPA incorporates
GTEx evidence). Do not average these scores or transfer their reference values.
NY-ESO-1 and XAGE1A/B illustrate why mandatory **coverage** must not silently become
a GTEx-only Reproductive expression floor. Both are well measured in GTEx, but their
GTEx Reproductive references are low. Their HPA evidence remains separately visible.

[HPA v23](https://v23.proteinatlas.org/about/download) provides a companion panel of
38 normal contexts / 29 anatomy units, including appendix, gallbladder and tonsil.
Thymus is separately **Permissive**, excluded from both sides of the identity
comparison. [All 50 context roles](hpa-v23-context-roles.csv). Disposable tissues
count in this identity assessment even though their expression
is permitted by the separate safety policy. [Panel weights](hpa-v23-normal-panel.csv) ·
[All HPA profiles](hpa-v23-protein-profiles.csv.gz).

## Reproductive is the umbrella

| Branch | Scope |
|---|---|
| Testicular | Testis and explicitly identified testicular cells |
| Other male | Other male reproductive organs and gametes |
| Ovarian | Ovary, follicles and oocytes |
| Female tract | Oviducts, uterus, cervix, vagina and other mapped tract tissues |
| Mammary | Mapped breast/mammary compartments |
| Placental | Mapped placenta and placental compartments |
| Other extraembryonic | Amnion, chorion, yolk sac, allantois and umbilical cord |
| Embryonic | Explicitly embryonic samples |
| Fetal | Explicitly fetal samples, including fetal organs |

These branches union to Reproductive. A protein can have evidence in multiple
branches. “Embryonic-only evidence” describes which measured route qualifies,
not proven absence from every other branch. Adult maternal organs do not change
category during pregnancy. The current developmental inputs cover embryonic stages;
they do not supply a comprehensive fetal panel. [Branch measurements](reproductive-branches.csv.gz).

## Open the original 617 / 649 buckets

The original release had **{len(old):,} unassigned proteins**. Its statuses combined
observed failures with evidence gaps. They did not mean 617 had only missing data
or 649 were unsafe. Here are disjoint combinations of actual failed requirements:

### Incomplete evidence: 617 in release 1.22.0

{markdown(reason_table(old[old.assessment_status.eq("incomplete_evidence")]))}

Missing developmental data or a multi-locus tumor profile can prevent the stricter
low-expression Neither determination without blocking an otherwise qualifying
adult-Reproductive CTA. The status tracked those gaps together. A failure of the
Reproductive expression floor means no qualifying RNA route in these inputs, not
proof that the protein is absent or that all reproductive IHC is negative.

### Outside definitions: 649 in release 1.22.0

{markdown(reason_table(old[old.assessment_status.eq("outside_category_definitions")]))}

Most fail only the association floor. For example FOXR2 has 1.2 HPA nTPM despite
measured cancer expression; this failure is not normal-tissue danger. MAGEA3, TBATA,
SPRR4 and LY6G6D have qualifying association but separate normal-expression holds.
Those holds must be inspected independently of missingness or diffuseness.

Current unassigned status counts after removing the maximum veto and requiring
GTEx: **{now.assessment_status.value_counts().to_dict()}**. Download the
[current incomplete combinations](incomplete_evidence-current.csv) and
[current outside-definition combinations](outside_category_definitions-current.csv).
Every record in the detailed download retains its release-1.22 status and current
reason, so a coverage change cannot masquerade as newly measured toxicity.

## Method basis and reproducibility

[A benchmark of gene expression tissue-specificity metrics](https://doi.org/10.1093/bib/bbw008)
compares continuous tissue-specificity measures and discusses sensitivity to tissue
sampling. Our Reproductive-anchored bounded average is a **custom metric**, not
the paper's Tau statistic and not a validated 10% clinical or biological boundary.

Rebuild with `python scripts/audit_cta_profile_concentration.py` after the current
selection replay. `CANCERDATA_DATA_DIR` selects the shared source-cache root.
`oncoref.expression_diffuseness.assess()` exposes the score, weights, missingness
bounds and reference-quality status for downstream use. CSV source values remain
numeric and retain explicit units; labels in the report link to explanatory pages.
"""
    (OUT / "index.md").write_text(report)
    (OUT / "manifest.json").write_text(
        json.dumps(
            {
                "status": "active HPA adult/placental 10% identity gate; GTEx descriptive",
                "metric": metric.POLICY_ID,
                "baseline_commit": "5d8dc1ca377c32c69548a6195540528b0c2e9aff",
                "inputs": input_hashes(
                    [
                        Path(__file__),
                        ROOT / "oncoref/expression_diffuseness.py",
                        ROOT / "scripts/cta_profile_data.py",
                        OUT / "before-diffuseness.csv.gz",
                        ROOT / "scripts/cta_source_links.py",
                        ROOT / "scripts/render_cta_source_walkthrough.py",
                        CURRENT / "protein-group-classifications.csv.gz",
                        BASELINE,
                        CURRENT / "hpa_rna_consensus-accepted-identities.csv.gz",
                        CURRENT / "alternative-evidence/manifest.json",
                    ]
                ),
                "sources": [
                    reference_data.provenance(s, v, verify_content=True)
                    for s, v in [("hpa_rna_consensus", "v23"), ("gtex_tissue_rna", "v11")]
                ],
                "outputs": {
                    p.name: sha256(p) for p in sorted(OUT.iterdir()) if p.name != "manifest.json"
                },
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    print(json.dumps({"core_complete": int(complete.sum()), "current_counts": counts}))


if __name__ == "__main__":
    build()
