#!/usr/bin/env python3
"""Build the all-cohort, optional-DICE/cornea, three-way CTA research funnel."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import plot_cta_geometric as style  # noqa: E402
from cta_geometric_sources_report import plot_routes  # noqa: E402
from cta_report_common import input_hashes  # noqa: E402
from cta_report_common import sha256 as digest  # noqa: E402
from cta_report_common import write_csv as write  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402

from oncoref import cta_tissues  # noqa: E402
from oncoref.cta_funnel import classify_layers  # noqa: E402

ASSAYS = [
    "Normal bulk",
    "HPA tissue clusters",
    "HPA immune",
    "Monaco immune",
    "DICE immune",
    "Adult cornea",
]
OPTIONAL = {"DICE immune", "Adult cornea"}


def markdown_table(frame: pd.DataFrame):
    """Render a small report table from the provided columns and rows."""
    return "\n".join(
        [
            "| " + " | ".join(map(str, frame.columns)) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |",
            *[
                "| " + " | ".join(map(str, row)) + " |"
                for row in frame.itertuples(index=False, name=None)
            ],
        ]
    )


def report(
    out: Path,
    result: pd.DataFrame,
    funnel: pd.DataFrame,
    independent: pd.DataFrame,
    scenarios: pd.DataFrame,
) -> None:
    counts = result.category.value_counts()
    previous_failures = result[result.positive_context & ~result.tumor_maximum_max_ptpm.ge(10)]
    previous_missing = previous_failures[previous_failures.tumor_maximum_max_ptpm.isna()]
    rescued = int(previous_failures.tumor_union_ge2.sum())
    remaining = int((~previous_failures.tumor_union_ge2).sum())
    cohort = json.loads((out / "tumor-union/summary.json").read_text())
    named = result[
        result.symbol.isin(
            [
                "PRAME",
                "CTAG2",
                "MAGEA4",
                "CTAG1A",
                "CTAG1B",
                "XAGE1A",
                "XAGE1B",
                "SSX2",
                "CT83",
                "TBXT",
                "SMC1B",
                "PAGE1",
                "FTMT",
            ]
        )
    ][["symbol", "category", "exclusion_reasons"]].copy()
    named["exclusion_reasons"] = named.exclusion_reasons.str.replace(
        " | ", "; ", regex=False
    ).replace("", "All required gates pass")
    content = f"""# Expanded tumor evidence and a three-way CTA funnel

The new research classification is **{counts["CTA"]:,} CTA / {counts["Unsafe"]:,} Unsafe /
{counts["Neither"]:,} Neither**, across all {len(result):,} reference genes. There are
**{int(result.quantitative_candidate.sum())} quantitative candidates** before accepted protein/safety
exclusions and unresolved antibody holds. These are newly computed categories,
not the old labels carried inside the previous 319 candidates. Released package
expression labels and forbidden sequence references are not overwritten.

## What changed

{markdown_table(scenarios.rename(columns={"scenario": "Scenario", "genes": "Genes"}))}

![Effects of each change](figures/policy-effects.png)

[Editable SVG](figures/policy-effects.svg), [scenario counts](policy-effects.csv).

Tumor admission now requires **at least one sample ≥2 in its source-native
TPM or HPA pTPM**. The shared normal-tissue/IHC grid is unchanged. The different
RNA scales remain labeled; they are not pooled, averaged or treated as calibrated
equivalents. DICE and cornea are optional for coverage: absent measurements do
not hold a gene, but supported measured expression still contributes exclusions.

## Cancer cohorts actually used

All **{cohort["used_cohorts"]} of the {cohort["registered_cohorts"]} published OncoRef matrices** marked
compatible with an absolute TPM threshold were processed, including adult and
pediatric cancers. The screen uses raw linear TPM from individual samples passing
the existing QC and absolute-TPM suitability checks, not clean-TPM, cohort medians
or representative samples. These are {cohort["sample_occurrences"]:,} sample occurrences,
or **{cohort["namespace_deduplicated_samples"]:,} identities after deduplicating shared source views**,
across {cohort["source_namespaces"]} source namespaces. This is not a unique-patient count;
the same patient can occur in independent studies or multiple biopsies.

The existing HPA 25.1 cancer reference adds its 21 TCGA and 10 validation cohorts
(8,384 source samples). HPA overlaps OncoRef, so those counts are not added to
claim a unique total. The eleven omitted matrices use microarray/CPM proxies,
single-cell pseudobulk scales or a known invalid legacy scale. Registry entries
without published sample matrices are not silently treated as acquired cohorts.

[Exact cohort inventory and exclusions](tumor-union/cohort-availability.csv),
[per-source sample identities](tumor-union/sample-identity.csv.gz),
[per-gene observed maxima and winning samples](tumor-union/tumor-maxima.csv.gz),
[source-matrix and QC hashes](tumor-union/manifest.json).

Of the previous **112 tumor failures, {rescued} now pass**: 42 through lowering the
HPA floor alone, then another {rescued - 42} through the wider cohort union. All 48 formerly
measured-below-10 genes now qualify; {int(previous_missing.tumor_union_ge2.sum())} of the 64 formerly unmeasured in HPA
now qualify. {remaining} remain below 2 in available, unambiguously matched source matrices. This still does
not establish that they are never expressed by any tumor.
[All previous 112 decisions](previous-112-tumor-failures.csv),
[remaining failures](remaining-tumor-failures.csv).

## The categories and every layer's contribution

- **Unsafe:** accepted nonpermitted normal protein, reliable nonconflicted HPA
  normal staining, an existing applicable clinical safety exclusion, or supported
  normal RNA at or above **twice** its applicable candidate ceiling.
- **Neither:** no Unsafe evidence, but at least one admission failure, mandatory
  coverage gap, uncertain/conflicted protein finding, unresolved antibody hold,
  or normal RNA above the candidate ceiling and below the upper boundary.
- **CTA:** quantitative reproductive/embryonic association and tumor admission,
  complete mandatory panels, normal RNA at or below each applicable candidate
  ceiling, and no protein/safety exclusion or unresolved antibody hold.

Unsafe is a **research policy exclusion**, not proof of protein presentation or
clinical toxicity. The upper RNA boundary previously appeared as a reference-review
annotation; this three-way proposal explicitly uses it for the Unsafe category.
It still does not automatically create or publish a forbidden-protein FASTA.

Every layer evaluates **all genes**, including genes held earlier. A later strong
RNA/protein finding can therefore change Neither to Unsafe. Final labels are
unchanged when layer order is reversed. Intermediate CTA counts mean *eligible
so far*, not final clearance. Uncertain HPA staining blocks CTA but does not alone
establish Unsafe. Antibody-concern flags prevent negative-IHC relaxation and remain
an explicit final review hold unless independent stronger evidence establishes Unsafe.

{markdown_table(funnel.rename(columns={"layer": "Layer", "new_candidate_removals": "Removed from eligibility here", "new_Unsafe": "New Unsafe", "CTA_eligible_so_far": "CTA eligible so far"}))}

![Three-way classification at every layer](figures/three-way-funnel.png)

[Editable SVG](figures/three-way-funnel.svg), [funnel counts](three-way-funnel.csv).

![Independent contributions from each source](figures/independent-evidence.png)

[Editable SVG](figures/independent-evidence.svg),
[overlapping evidence and leave-one-layer-out effects](independent-layer-effects.csv).
Independent source counts overlap and must not be added together.

## Normal thresholds and coverage

| Candidate ceiling | Missing/unusable IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 2 | 4 | 8 | 16 |
| Protected | 4 | 8 | 16 | 32 |
| Manageable | 8 | 16 | 32 | 64 |

Equality passes the candidate ceiling. Unsafe RNA evidence begins at **2×**
each entry. Bulk, supported HPA clusters and sorted immune assays use the linked
native TPM/nTPM values. Corneal CP10K uses the nominal table value divided by 100;
that rescaling is not an empirical conversion to HPA-equivalent abundance.
Immune populations default to Protected. Explicit reproductive/developmental
immune contexts are neutral. HPA clusters require ≥20 reads; DICE uses the third
highest reported value; cornea requires ≥5 positive cells and ≥20 counts in at
least two sample groups. Matching negative IHC retains the existing profile
control, exact context and antibody-conflict requirements.

Mandatory coverage remains 50 bulk tissues, 426 observed HPA clusters, 19 HPA
immune summaries and 30 Monaco summaries, plus the required numeric fields.
In the expanded funnel, **84** otherwise eligible genes lack mandatory coverage;
all lack complete bulk coverage, and some also lack the other required panels.
The old 122 coverage failures included **91 cornea-only gaps**, which are no
longer mandatory holds. This does not rescue all 91 through the remaining filters.
The final CTA set contains **10 genes with optional-panel gaps**, explicitly flagged.

[Mandatory gaps across the reference universe](mandatory-coverage-gaps.csv.gz),
[CTA genes with optional gaps](cta-optional-coverage-gaps.csv),
[normal measurements closest to or above each ceiling](normal/normal-peak-evidence.csv.gz).

## Bulk gaps: names and historical IDs were checked

None of the **31 previously flagged bulk gaps** can be filled through a registered
alias present in the pinned source. Each has 41 of 50 tissue entries, missing
amygdala, basal ganglia, cerebellum, hippocampal formation, hypothalamus, midbrain,
pituitary gland, spinal cord and vagina. They are absent source rows, not zero RNA.

| Example | Source identity result |
|---|---|
| POLR2A — ENSG00000181222 | Same stable ID and symbol; only 41 source tissues |
| EEF1AKMT4 — ENSG00000284753 | Same stable ID and symbol; only 41 source tissues |
| NOTCH2NLB — ENSG00000286019 | Same stable ID and symbol; only 41 source tissues |
| CIMIP3 — ENSG00000287363 | Source name GUCA1ANB, but already the same stable ID; renaming cannot restore missing rows |
| ACTL10 — ENSG00000288649 | Registered previous ID ENSG00000182584 is absent from this source |
| POLR2J3 — ENSG00000285437 | Another same-name ID, ENSG00000168255, coexists with 50 tissues; distinct values and no registered identity mapping mean it cannot be borrowed |

POLR2J3's two source IDs have maxima of 73.2 and 3.8 nTPM respectively, illustrating
why matching a symbol is not sufficient to merge measurements. A separate fully
absent canonical ID, **PAXX ENSG00000310560**, was recovered using the registered
historical ID ENSG00000148362: all 50 source rows are retained with their original
ID. This does not rescue any of the earlier 31 partial-coverage cases. This audit
does not claim that other GTEx/HPA releases lack the missing observations.

[The previous 31 cases](previous-31-bulk-gaps.csv),
[partial bulk-coverage identity audit](normal/bulk-missing-identity-audit.csv),
[registered-alias observations actually recovered](normal/bulk-alias-filled-observations.csv).

## Named genes and complete downloads

{markdown_table(named.rename(columns={"symbol": "Gene", "category": "New category", "exclusion_reasons": "Reason"}))}

Of 360 quantitative candidates, FTMT moves to Unsafe through accepted normal
protein and 54 move to Neither through unresolved antibody review. TBXT/brachyury
is in that latter group: it qualifies quantitatively, but its protein-evidence
uncertainty is not erased. PRAME and CTAG2 remain excluded under the unchanged
normal-expression policy. There are no named-gene admission exceptions.

[Final CTA list](cta-genes.csv), [all gene classifications and every reason](gene-classifications.csv.gz),
[Unsafe genes](unsafe-genes.csv.gz), [Neither genes](neither-genes.csv.gz),
[360 quantitative candidates](quantitative-candidates.csv),
[55 quantitative candidates requiring exclusion or review](quantitative-candidate-review-exclusions.csv),
[named-gene evidence](named-gene-decisions.csv).

![Association-source overlap before filtering and in final CTA](figures/quantitative-source-overlap.png)

[Editable SVG](figures/quantitative-source-overlap.svg), [source counts](source-route-overlaps.csv).

## Reproduction and checks

```bash
python scripts/cta_tumor_union.py docs/audits/cta-expanded-funnel-20261027/tumor-union
python scripts/cta_expanded_normal.py docs/audits/cta-expanded-funnel-20261027/normal
MPLCONFIGDIR=/tmp/oncoref-matplotlib-geometric python scripts/cta_expanded_funnel.py docs/audits/cta-expanded-funnel-20261027
```

The all-gene normal assessment reproduces the previous evaluated normal decisions
and the 319-gene baseline; the vectorized bulk match reproduces 115,440 prior
contexts. Final categories are tested for order independence and all stage totals
sum to 20,097. Input/source/output hashes are recorded in the [manifest](manifest.json)
and the two component manifests. Every figure has PNG and editable SVG versions.
"""
    (out / "index.md").write_text(content)


def coverage(facts: pd.DataFrame, *, optional: bool = True) -> pd.Series:
    """Return mandatory-panel completeness, retaining optional-panel policy."""
    required = (
        facts.normal_bulk_tissues.eq(50)
        & facts.cell_records.eq(426)
        & facts.hpa_immune_cell_records.eq(19)
        & facts.hpa_immune_monaco_records.eq(30)
        & facts.reproductive_ntpm.notna()
        & facts.adult_bulk_max_ntpm.notna()
        & facts.hpa_immune_cell_max_nTPM.notna()
        & facts.hpa_immune_monaco_max_TPM.notna()
        & (facts.cell_supported_records.eq(0) | facts.cell_max_ntpm.notna())
    )
    if not optional:
        required &= (
            facts.dice_sources.eq(15)
            & facts.cornea_records.eq(84)
            & facts.dice_third_max_tpm.notna()
        )
    return required


def plot_funnel(
    funnel: pd.DataFrame, independent: pd.DataFrame, scenarios: pd.DataFrame, out: Path
) -> None:
    style.setup()
    colors = [style.TEAL, style.RED, "#BEA477"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 10), gridspec_kw={"width_ratios": [1.8, 1]})
    fig.subplots_adjust(left=0.19, right=0.97, top=0.82, bottom=0.16, wspace=0.24)
    style.title(
        fig,
        "Every layer can contribute to the final category",
        "All 20,097 genes are evaluated at every layer—even after an earlier hold. Intermediate CTA = eligible so far.",
    )
    y = np.arange(len(funnel))
    left = np.zeros(len(funnel))
    for col, label, color in zip(
        ["CTA_eligible_so_far", "Unsafe", "Neither"],
        ["CTA / eligible so far", "Unsafe", "Neither"],
        colors,
        strict=True,
    ):
        axes[0].barh(y, funnel[col], left=left, color=color, label=label, height=0.65)
        left += funnel[col]
    axes[0].set_yticks(y, funnel.layer, fontsize=10)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Genes · disjoint categories at each stage")
    axes[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=3, fontsize=9)
    axes[0].set_xlim(0, 20097)
    axes[1].axis("off")
    axes[1].set_ylim(len(funnel) - 0.5, -0.5)
    axes[1].set_xlim(0, 1)
    for x, title in zip(
        [0.12, 0.43, 0.73, 0.98], ["CTA*", "Unsafe", "Neither", "New Unsafe"], strict=True
    ):
        axes[1].text(x, -0.85, title, ha="right", fontsize=10, weight="bold")
    for i, r in enumerate(funnel.itertuples()):
        for x, n in zip(
            [0.12, 0.43, 0.73, 0.98],
            [r.CTA_eligible_so_far, r.Unsafe, r.Neither, r.new_Unsafe],
            strict=True,
        ):
            axes[1].text(
                x,
                i,
                f"{n:,}",
                ha="right",
                va="center",
                fontsize=11,
                weight="bold" if i == len(funnel) - 1 else "normal",
            )
    style.footer(
        fig,
        "*Final CTA requires every layer. Unsafe = research-policy exclusion, not proven toxicity. Neither includes insufficient or conflicting evidence.\nTumor admission: ≥2 native TPM/pTPM in any qualifying sample. DICE and cornea: optional coverage; observed positive evidence still applies.",
    )
    style.save(fig, out, "three-way-funnel")
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(left=0.28, right=0.94, top=0.79, bottom=0.18)
    style.title(
        fig,
        "Evidence remains visible after an earlier exclusion",
        "Independent counts overlap. A gene can have evidence from bulk, clusters, immune cells and protein simultaneously.",
    )
    ax.barh(
        y - 0.17,
        independent.candidate_exclusions_anywhere,
        height=0.3,
        color="#BEA477",
        label="Any candidate hold / exclusion",
    )
    ax.barh(
        y + 0.17,
        independent.Unsafe_evidence_anywhere,
        height=0.3,
        color=style.RED,
        label="Unsafe evidence",
    )
    ax.set_yticks(y, independent.layer, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, independent.candidate_exclusions_anywhere.max() * 1.12)
    for i, row in enumerate(independent.itertuples()):
        for offset, count in [
            (-0.17, row.candidate_exclusions_anywhere),
            (0.17, row.Unsafe_evidence_anywhere),
        ]:
            ax.text(count + 90, i + offset, f"{count:,}", va="center", fontsize=8)
    ax.set_xlabel("Genes · full 20,097-gene denominator")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=2, fontsize=9)
    style.footer(
        fig,
        "Coverage and nomination failures alone never establish Unsafe. Reliable protein and the upper linked RNA thresholds can establish a research exclusion.\nThe companion table also reports how many categories change if each individual layer is removed.",
    )
    style.save(fig, out, "independent-evidence")
    fig, ax = plt.subplots(figsize=(12, 6.5))
    fig.subplots_adjust(left=0.46, right=0.92, top=0.76, bottom=0.2)
    style.title(
        fig,
        "What changed the candidate count",
        "Separate the lower tumor floor, broader cohort union, optional coverage and final evidence reconciliation.",
    )
    ax.barh(np.arange(len(scenarios)), scenarios.genes, color=[style.MUTED] * 4 + [style.TEAL])
    ax.set_yticks(np.arange(len(scenarios)), scenarios.scenario, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, scenarios.genes.max() * 1.2)
    for i, n in enumerate(scenarios.genes):
        ax.text(n + 4, i, str(n), va="center", weight="bold", fontsize=13)
    ax.set_xlabel("Genes")
    style.footer(
        fig,
        "The first four rows are quantitative candidates. The last row additionally applies accepted normal-protein, clinical-safety and antibody-review rules.\nOptional coverage does not discard observed DICE or corneal RNA evidence. No gene-specific positive-control admission rule is used.",
    )
    style.save(fig, out, "policy-effects")


def build(out: Path, baseline: Path) -> None:
    """Build the versioned audit tables, report and source provenance."""
    f = pd.read_csv(baseline / "gene-comparison.csv.gz", low_memory=False)
    bulk = pd.read_csv(
        out / "normal/bulk-contexts.csv.gz",
        usecols=["Ensembl_Gene_ID", "tissue", "rna_ntpm"],
    )
    bulk_fields = {
        "normal_bulk_tissues": bulk.groupby("Ensembl_Gene_ID").rna_ntpm.count(),
        "reproductive_ntpm": bulk[bulk.tissue.isin(cta_tissues.CORE_REPRODUCTIVE_TISSUES)]
        .groupby("Ensembl_Gene_ID")
        .rna_ntpm.max(),
        "adult_bulk_max_ntpm": bulk[~bulk.tissue.isin(cta_tissues.EXTENDED_REPRODUCTIVE_TISSUES)]
        .groupby("Ensembl_Gene_ID")
        .rna_ntpm.max(),
    }
    for column, values in bulk_fields.items():
        f["previous_" + column] = f[column]
        f[column] = f.gene_id.map(values)
    f["reproductive_route"] = f.reproductive_ntpm.ge(5) & f.reproductive_ntpm.gt(
        f.adult_bulk_max_ntpm
    )
    f["positive_context"] = f.reproductive_route | f.developmental_route
    t = pd.read_csv(out / "tumor-union/tumor-maxima.csv.gz")
    f = f.merge(
        t.rename(columns={c: "oncoref_" + c for c in t if c != "gene_id"}),
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    f = f.copy()
    normal = pd.read_csv(out / "normal/normal-assessments.csv.gz")
    empty = pd.Series(False, index=f.index)
    normal_hold, normal_strong = {}, {}
    for assay in ASSAYS:
        a = normal[normal.assay.eq(assay)].set_index("gene_id")
        normal_hold[assay] = (
            f.gene_id.map(a.candidate_exclusion).fillna(assay not in OPTIONAL).astype(bool)
        )
        normal_strong[assay] = f.gene_id.map(a.unsafe_evidence).fillna(False).astype(bool)
        # Existing 319-gene audit is an independent regression gate for the full
        # normal evidence re-evaluation, restricted to its evaluated universe.
        pool = f.ihc_contexts_evaluated
        if not normal_hold[assay][pool].eq(~f.loc[pool, "pass_" + assay]).all():
            raise ValueError(f"Normal assay does not replay earlier selection: {assay}")
    f["hpa_tumor_ge2"] = f.tumor_maximum_max_ptpm.ge(2)
    f["oncoref_tumor_ge2"] = f.oncoref_maximum_tpm.ge(2)
    f["tumor_union_ge2"] = f.hpa_tumor_ge2 | f.oncoref_tumor_ge2
    f["required_coverage"] = coverage(f)
    f["optional_DICE_complete"] = f.dice_sources.eq(15) & f.dice_third_max_tpm.notna()
    f["optional_cornea_complete"] = f.cornea_records.eq(84)
    layers = [
        ("Reproductive / embryonic association", ~f.positive_context, empty),
        ("Tumor sample ≥2 TPM/pTPM", ~f.tumor_union_ge2, empty),
        ("HPA normal protein", f.ihc_detected, f.ihc_reliable_detected & ~f["Antibody concern"]),
        ("Required RNA coverage", ~f.required_coverage, empty),
        *[(assay, normal_hold[assay], normal_strong[assay]) for assay in ASSAYS],
        (
            "Accepted normal-protein evidence",
            f["Accepted normal protein"],
            f["Accepted normal protein"],
        ),
        ("Clinical safety exclusion", f["Safety veto"], f["Safety veto"]),
        ("Unresolved antibody review", f["Antibody concern"], empty),
    ]
    result, funnel, independent = classify_layers(f, layers)
    result["policy"] = "oncoref.expanded-three-way.2026-10-27.1"
    normal_ok = ~pd.DataFrame(normal_hold).any(axis=1)
    common = f.positive_context & ~f.ihc_detected & normal_ok
    scenarios = pd.DataFrame(
        [
            {
                "scenario": "HPA ≥10; all panels required",
                "genes": int(
                    (common & f.tumor_maximum_max_ptpm.ge(10) & coverage(f, optional=False)).sum()
                ),
            },
            {
                "scenario": "HPA ≥2; all panels required",
                "genes": int((common & f.hpa_tumor_ge2 & coverage(f, optional=False)).sum()),
            },
            {
                "scenario": "HPA + OncoRef ≥2; all panels required",
                "genes": int((common & f.tumor_union_ge2 & coverage(f, optional=False)).sum()),
            },
            {
                "scenario": "HPA + OncoRef ≥2; DICE/cornea optional",
                "genes": int((common & f.tumor_union_ge2 & coverage(f)).sum()),
            },
            {
                "scenario": "Final CTA after evidence reconciliation",
                "genes": int(result.category.eq("CTA").sum()),
            },
        ]
    )
    assert scenarios.genes.iloc[0] == 319
    assert (funnel.CTA_eligible_so_far + funnel.Unsafe + funnel.Neither).eq(len(f)).all()
    reverse, _, _ = classify_layers(f, list(reversed(layers)))
    assert reverse.category.equals(result.category)
    quantitative = common & f.tumor_union_ge2 & coverage(f)
    result["quantitative_candidate"] = quantitative
    result["previous_geometric_candidate"] = f.two_axis_doubling_maximum10_selected
    write(result, out / "gene-classifications.csv.gz")
    for category in ["CTA", "Unsafe", "Neither"]:
        write(result[result.category.eq(category)], out / f"{category.lower()}-genes.csv.gz")
    write(
        result.loc[
            result.category.eq("CTA"),
            [
                "gene_id",
                "symbol",
                "category",
                "reproductive_route",
                "developmental_route",
                "reproductive_ntpm",
                "adult_bulk_max_ntpm",
                "developmental_mean_tpm",
                "tumor_maximum_max_ptpm",
                "tumor_maximum_cancer_code",
                "tumor_maximum_cohort",
                "oncoref_maximum_tpm",
                "oncoref_maximum_sample",
                "oncoref_cancer_code",
                "oncoref_source_cohort",
                "optional_DICE_complete",
                "optional_cornea_complete",
            ],
        ],
        out / "cta-genes.csv",
    )
    write(result[quantitative], out / "quantitative-candidates.csv")
    write(
        result[quantitative & ~result.category.eq("CTA")],
        out / "quantitative-candidate-review-exclusions.csv",
    )
    write(funnel, out / "three-way-funnel.csv")
    write(independent, out / "independent-layer-effects.csv")
    write(scenarios, out / "policy-effects.csv")
    write(
        result[f.positive_context & ~f.tumor_maximum_max_ptpm.ge(10)],
        out / "previous-112-tumor-failures.csv",
    )
    write(result[f.positive_context & ~f.tumor_union_ge2], out / "remaining-tumor-failures.csv")
    write(result[~f.required_coverage], out / "mandatory-coverage-gaps.csv.gz")
    write(
        result[
            result.category.eq("CTA") & ~(f.optional_DICE_complete & f.optional_cornea_complete)
        ],
        out / "cta-optional-coverage-gaps.csv",
    )
    write(
        result[
            result.symbol.isin(
                [
                    "PRAME",
                    "MAGEA4",
                    "CTAG1A",
                    "CTAG1B",
                    "XAGE1A",
                    "XAGE1B",
                    "CTAG2",
                    "SSX2",
                    "CT83",
                    "TBXT",
                    "SMC1B",
                    "PAGE1",
                    "FTMT",
                ]
            )
        ],
        out / "named-gene-decisions.csv",
    )
    figdir = out / "figures"
    figdir.mkdir(exist_ok=True)
    plot_funnel(funnel, independent, scenarios, figdir)
    routes = []
    for scope, mask in [
        ("association", f.positive_context),
        ("selected", result.category.eq("CTA")),
    ]:
        for name, route in [
            ("HPA only", f.reproductive_route & ~f.developmental_route),
            ("Tyser only", f.developmental_route & ~f.reproductive_route),
            ("Both", f.reproductive_route & f.developmental_route),
        ]:
            routes.append({"scope": scope, "route": name, "genes": int((mask & route).sum())})
    routes = pd.DataFrame(routes)
    write(routes, out / "source-route-overlaps.csv")
    plot_routes(routes, figdir)
    gaps = pd.read_csv(out / "normal/bulk-missing-identity-audit.csv")
    earlier_gaps = f[
        f.two_axis_doubling_first_failure.eq("Complete required RNA panels")
        & f.previous_normal_bulk_tissues.lt(50)
    ]
    write(
        earlier_gaps[["gene_id"]].merge(gaps, on="gene_id", validate="one_to_one"),
        out / "previous-31-bulk-gaps.csv",
    )
    report(out, result, funnel, independent, scenarios)
    summary = {
        "policy": result.policy.iloc[0],
        "categories": result.category.value_counts().to_dict(),
        "quantitative_candidates": int(quantitative.sum()),
        "tumor_failures_after_association": int((f.positive_context & ~f.tumor_union_ge2).sum()),
        "optional_gaps_in_CTA": int(
            (
                result.category.eq("CTA") & ~(f.optional_DICE_complete & f.optional_cornea_complete)
            ).sum()
        ),
        "historical_classes_changed": False,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    inputs = [
        Path(__file__),
        ROOT / "oncoref/cta_funnel.py",
        ROOT / "scripts/plot_cta_geometric.py",
        ROOT / "scripts/cta_geometric_sources_report.py",
        ROOT / "scripts/cta_cell_policy_report.py",
        ROOT / "scripts/cta_quantitative_threshold_report.py",
        ROOT / "oncoref/expression.py",
        ROOT / "oncoref/gene_ids.py",
        ROOT / "oncoref/cell_rna_policy.py",
        ROOT / "oncoref/cta_normal_screen.py",
        ROOT / "oncoref/normal_cell_data.py",
        ROOT / "oncoref/normal_cell_screen.py",
        ROOT / "oncoref/normal_tissues.py",
        ROOT / "oncoref/data/ensembl-id-aliases.csv",
        ROOT / "oncoref/data/canonical-gene-space.csv.gz",
        ROOT / "oncoref/data/ncbi-symbol-synonyms.csv.gz",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/normal-tissue-mappings.csv",
        baseline / "gene-comparison.csv.gz",
        out / "normal/manifest.json",
        out / "tumor-union/manifest.json",
    ]
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(inputs),
                "outputs": {
                    str(p.relative_to(out)): digest(p)
                    for p in [
                        *out.glob("*.csv*"),
                        out / "summary.json",
                        out / "index.md",
                        *figdir.glob("*"),
                    ]
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(scenarios.to_string(index=False))
    print(funnel.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--baseline", type=Path, default=ROOT / "docs/audits/cta-geometric-thresholds-20261027"
    )
    args = parser.parse_args()
    build(args.output, args.baseline)
