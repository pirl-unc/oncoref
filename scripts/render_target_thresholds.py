#!/usr/bin/env python3
"""Render the threshold comparison and an explicit recommended research policy."""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from audit_target_thresholds import OUT  # noqa: E402
from build_target_sample_ranks import sha  # noqa: E402

from oncoref import target_selection  # noqa: E402

PREFERRED = "fold_L1_R1_S10"
LABELS = {
    "fold_L1_R1_S5": "Current: two patients above normal, or one at 5×",
    PREFERRED: "Recommended: two patients above normal, or one at 10×",
    "fold_L1_R1_S20": "Two patients above normal, or one at 20×",
    "fold_L1_R2_S5": "Two patients at 2×, or one at 5×",
    "rank_P80_P95_R1": "Within-tumor p80 repeated / p95 single, plus current enrichment",
    "cohort_Q95_L1_R1_S0": "Across-patient p95 alone, in any eligible cancer cohort",
    "cohort_Q95_L1_R1_S5": "Across-patient p95, or one tumor at 5×",
}


def build():
    grid = pd.read_csv(OUT / "threshold-grid.csv")
    evidence = pd.read_csv(OUT / "protein-threshold-evidence.csv.gz").set_index("entity_id")
    decisions = pd.read_csv(OUT / "cancer-rule-decisions.csv.gz").set_index("entity_id")
    summary = json.loads((OUT / "summary.json").read_text())
    uncapped = grid[grid.disposable_cap_tpm.isna()].set_index("id")
    selected = uncapped.loc[list(LABELS)].copy()
    selected["rule"] = pd.Series(LABELS)
    selected.to_csv(OUT / "selected-rules.csv")
    # All maximum and 32-TPM cap checks are also applied to the proposal.
    eligible = (
        evidence.target_maximum_pass
        & evidence.disposable_max_tpm.le(32)
        & evidence.disposable_max_ntpm.le(32)
    )
    cancer = decisions[PREFERRED] & eligible
    reproductive = evidence.is_reproductive_antigen & eligible
    proposed = np.select(
        [cancer & reproductive, cancer, reproductive],
        ["CTA", "Cancer Antigen only", "Reproductive Antigen only"],
        default="Other",
    )
    possible_repeated = evidence.independent_patients_above_normal.add(
        evidence.unresolved_patient_samples_above_normal
    ).ge(2)
    cancer_status = np.select(
        [decisions[PREFERRED], possible_repeated],
        ["qualifies", "unassessed"],
        default="does_not_qualify",
    )
    replay = target_selection.classify(evidence.assign(unsafe=False, cancer_evidence=cancer_status))
    # Proposed cap/maximum rules do not remove any current candidate here.
    assert (replay.category.to_numpy() == proposed).all()
    evidence["proposed_cancer_evidence"] = cancer_status
    evidence["proposed_classification"] = replay.classification
    changes = evidence[evidence.category.ne(proposed)].copy()
    changes["proposed_category"] = pd.Series(proposed, index=evidence.index).loc[changes.index]
    changes.to_csv(OUT / "recommended-category-changes.csv")
    retained = evidence.assign(proposed_category=proposed)
    retained.to_csv(
        OUT / "recommended-eligible-inventory.csv.gz", compression={"method": "gzip", "mtime": 0}
    )

    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    heat = uncapped[
        (uncapped.kind.eq("fold"))
        & uncapped.floor.eq(1)
        & uncapped.repeated_fold.gt(0)
        & uncapped.single_fold.gt(0)
    ]
    for ax, value, title in zip(
        axes,
        ["Cancer Antigens", "cancer_lost_without_largest_sample"],
        [
            "Cancer Antigens passing",
            "Lose cancer qualification when their\nlargest tumor sample is removed",
        ],
    ):
        table = heat.pivot(index="repeated_fold", columns="single_fold", values=value).sort_index()
        ax.imshow(table.to_numpy(), cmap="Blues")
        ax.set_xticks(range(len(table.columns)), [f"{int(v)}×" for v in table.columns])
        ax.set_yticks(range(len(table.index)), [f"{int(v)}×" for v in table.index])
        ax.set_xlabel("Single-tumor enrichment")
        ax.set_ylabel("Enrichment in two patients")
        ax.set_title(title, fontsize=12, pad=12)
        for i, j in np.ndindex(table.shape):
            ax.text(
                j,
                i,
                str(table.iloc[i, j]),
                ha="center",
                va="center",
                color="white"
                if (table.iloc[i, j] - table.min().min())
                > 0.55 * (table.max().max() - table.min().min())
                else "#172b3a",
                fontweight="bold",
            )
        ax.add_patch(plt.Rectangle((0.5, -0.5), 1, 1, fill=False, edgecolor="#d56818", lw=3))
    fig.suptitle(
        "A stronger repeated-expression gate can increase reliance on single tumors", fontsize=14
    )
    fig.savefig(OUT / "enrichment-grid.png", dpi=180)
    fig.savefig(OUT / "enrichment-grid.svg")
    plt.close(fig)
    counts = pd.read_csv(OUT / "sample-rank-cutoffs.csv.gz")
    quantile_stats = counts[["p80_tpm", "p90_tpm", "p95_tpm"]].quantile([0, 0.1, 0.5, 0.9, 1])
    quantile_stats.to_csv(
        OUT / "within-sample-cutoff-summary.csv", index_label="quantile_across_samples"
    )
    fig, ax = plt.subplots(figsize=(9.4, 4.8), layout="constrained")
    labels = [
        "Current",
        "Recommended\nsingle 10×",
        "Within-tumor\np80 / p95",
        "Across-patient\np95 only",
        "Across-patient p95\n+ single 5×",
    ]
    ids = [
        "fold_L1_R1_S5",
        PREFERRED,
        "rank_P80_P95_R1",
        "cohort_Q95_L1_R1_S0",
        "cohort_Q95_L1_R1_S5",
    ]
    bars = ax.bar(
        labels,
        uncapped.loc[ids, "Cancer Antigens"],
        color=["#627d8a", "#126e82", "#aa7742", "#aa7742", "#aa7742"],
    )
    ax.bar_label(bars, padding=4)
    ax.set_ylim(0, 1000)
    ax.set_ylabel("Cancer Antigens, including CTAs")
    ax.set_title(
        "Percentile gates select abundance or prevalence, not the same biological property"
    )
    fig.savefig(OUT / "percentile-comparison.png", dpi=180)
    fig.savefig(OUT / "percentile-comparison.svg")
    plt.close(fig)

    table = "\n".join(
        f"| {LABELS[i]} | {r['Cancer Antigens']} | {r.CTA} | {r['Cancer Antigen only']} | {r['Reproductive Antigen only']} | {r.current_targets_to_other} |"
        for i, r in selected.iterrows()
    )
    paragraphs = f"""# A simpler expression rule: threshold comparison

**Recommended for the next policy draft: retain detection at 1 native TPM/nTPM; require two independent tumors above the normal maximum, or one tumor at least tenfold above it; add separate Disposable caps of 32 TPM and 32 nTPM in their respective panels. Keep percentiles as abundance and prevalence annotations.**

With the existing normal, translation, identity and sequence screen unchanged, this gives **849 Cancer Antigens: 617 CTAs and 232 Cancer Antigens only**. Another **43** are Reproductive Antigens only. The union is **892 Target candidates**, with **742 Other and 18,259 Unsafe**. This recommendation is a counterfactual; the rebuilt catalog still uses the current fivefold single-tumor rule.

## One decision flow

1. **Apply shared eligibility.** Retain the established normal-tissue exclusions, required coverage, protein identity, translation and sequence rules for both antigen sets. Unsafe takes precedence. Missing required evidence goes to Other:Review.
2. **Establish association.** Use 1 native TPM/nTPM with the existing read/cell support rules, or qualifying protein evidence. Reproductive includes adult reproductive, placental, extraembryonic, embryonic and fetal contexts.
3. **Establish cancer enrichment.** Let N be the maximum supported nonreproductive, non-Disposable normal value and E = T / max(N, 1). Require two independently identified patients with T ≥ 1 and T > N, **or one qualified tumor with E ≥ 10**. Preserve the exact qualifying sample and route. A single tumor remains sufficient for a rare strong signal.
4. **Bound Disposable expression.** Require max(cancer, reproductive/developmental) to exceed every measured nonreproductive normal/Disposable value, and Disposable expression ≤32 TPM and ≤32 nTPM, assessed separately on each source scale. The 32-TPM cap continues the 4/8/16/32 ordering of nominal tier ceilings as a screening convention; it is not a measured tolerance limit. These are two nominal source-scale caps, not a conversion between TPM and nTPM; existing IHC adjustments remain unchanged.
5. **Assign membership.** Cancer plus Reproductive = CTA; one qualifying branch gives its respective Antigen-only category; neither gives Other. Curated paper nominations remain a separate linked evidence annotation.

The recommendation changes only the proposed single-tumor threshold and adds explicit maximum/cap requirements. The maximum and separate 32-unit caps remove no current candidates in the measured panels. **15 current CTAs become Reproductive Antigens only; 18 Cancer Antigens only become Other:Review.** All 18 have multiple above-normal samples, but insufficiently resolved patient IDs. Their review reason is **patient identity unresolved**; resolving independent patients could restore membership through the repeated route. [All 33 changes](recommended-category-changes.csv) include proposed evidence status, values, tumor witnesses, transcripts and retained paper IDs. Losing cancer qualification does not establish normal expression or Unsafe status.

## What was compared

The sweep tested **{summary["rules"]} cancer-expression rules × six Disposable-cap settings = {summary["scenarios"]} scenarios** over **{summary["physical_samples"]:,} physical tumor samples** and **957 proteins passing the existing shared screen**. It varies detection floors (1, 2, 5, 10), repeated enrichment (1×, 2×, 3×, 5×), single-tumor enrichment (5×, 10×, 20×), within-tumor p80/p90/p95 ranks, across-patient p80/p90/p95 values, and Disposable caps (none, 4, 8, 16, 32, 128, applied independently to TPM and nTPM observations). Normal-exclusion thresholds are held fixed; this sweep does not optimize them against the proteins they already admitted.

| Cancer rule | Cancer Antigens total | CTA | Cancer Antigen only | Reproductive Antigen only | Current targets moving to Other |
|---|---:|---:|---:|---:|---:|
{table}

All rows also require the target maximum to exceed the measured nonreproductive maximum. These selected counts are identical with a 32-TPM Disposable cap. The [full filterable grid](report.html), [CSV grid](threshold-grid.csv) and [protein-by-rule decisions](cancer-rule-decisions.csv.gz) contain every scenario.

## Why this recommendation

**Keep the repeated route permissive about abundance.** Raising its enrichment from >1× to ≥2× while leaving the single-tumor rule at 5× reduces Cancer Antigens from 882 to 836. But it increases the number dependent on HL/Merkel from **40 to 86**, and the number lost after removing each protein's largest tumor sample from **58 to 101**. It discards repeated modest signals and leaves more proteins dependent on a high isolated value.

**Tighten the single-tumor route instead.** Moving from 5× to 10× gives 849 Cancer Antigens, retains **425 of the 430 current published cancer/CTA nominations**, and reduces the corresponding dependencies to **31 on HL/Merkel** and **49 on the largest sample**. A 20× rule gives 821 Cancer Antigens, retains 423 nominations, and reduces those dependencies to 13 and 27. Tenfold is a simple intermediate convention, not an empirically proven optimum. Published catalog nominations describe prior claims and are not validation labels used to fit the rule.

![Enrichment and dependence on one tumor](enrichment-grid.png)

Each largest-sample analysis removes one maximum separately for each protein, then reevaluates both routes. Repeated samples from the same known patient count once; another measurement from that patient can remain. The HL/Merkel analysis removes both cohorts. Dependence is not proof that a measurement is false. The [source audit](../scrutiny/index.md) still requires quantification and mapping review.

## The two different meanings of p95

**Within one tumor:** p95 means the protein ranks among the highest 5% of measured catalog proteins in that sample. The ranking denominator includes Unsafe and Other proteins, not just the candidate set. At the median sample, p80, p90 and p95 correspond to about **30, 62 and 117 TPM**, respectively. The p95 cutoff ranges from **2.89 to 410.55 TPM** across sources and samples. Rank is relative to the rest of the transcriptome; it is not a fixed abundance.

The tested p80/p95 gate requires p80 in the repeated route and p95 in the single-tumor route **in addition to** the existing enrichment rules. It leaves 300 Cancer Antigens and only 168 of the 430 current published nominations. For example, ESX1 peaks at 114.60 TPM but lies just below p95 in that sample; it loses cancer qualification under this rule. SSU72L2 remains, while TEDDM1 and TMEM269 lose cancer qualification. These are abundance filters, not evidence that excluded proteins are absent from tumors.

**Across patients within a cancer type:** p95 is the upper-tail expression value for that protein in that cohort. This tests prevalence. We use the observed lower order statistic, never interpolate a positive value between a zero and one extreme tumor. At least 20 measured, independently identified patients are required for each protein/cohort; repeated measurements from a patient contribute one maximum.

Only **{summary["cohorts_with_20_known_patients"]} of {summary["cohorts"]} cohort labels** have at least 20 documented patients. **{summary["unknown_patient_samples"]:,} samples lack resolved patient IDs** and cannot establish replication or patient percentiles. They remain eligible for the single-sample route. Unassessed cohorts are not negative cohorts. Percentiles are computed within cancer types, not across a mixture dominated by large common-cancer cohorts.

Across-patient p95 alone gives 400 Cancer Antigens. It removes SSU72L2's cancer membership despite a 990.91-TPM tumor maximum: its high RMS source lacks resolved patient IDs, and the highest p95 in eligible cohorts is 0.50 TPM. Adding the current rare-sample route yields 809 Cancer Antigens, **802 of which already qualify through that single-sample route**. The percentile requirement then adds complexity while doing little to control single-sample dependence.

![Percentile alternatives](percentile-comparison.png)

Use both percentile axes in the tables: **within-tumor abundance rank** and **within-cancer prevalence/percentiles**. Keep them out of the basic membership rule. The HPA likewise separates detection from enrichment in its [methods](https://v23.proteinatlas.org/about/assays%2Bannotation); its tissue-specificity definitions do not validate a therapeutic antigen threshold.

## Coverage and category boundaries

The Disposable TPM panel contains **appendix, gallbladder and tonsil**. It lacks isolated hair-follicle and melanocyte TPM. The HPA v23 consensus bulk nTPM panel independently peaks at 23.4 nTPM for LCE3C and 4.9 nTPM for NLRP4, both in tonsil. Adding these separate nTPM caps does not change any scenario count. [Every Disposable bulk nTPM observation](disposable-bulk-ntpm.csv) is retained. The separate [HPA melanocyte and anatomy audit](../expression-maxima/index.md) retains nTPM, reads and missing loci. A stricter number does not replace these missing compartments. Olfactory, bone and dental compartments still need explicit curation and expression coverage.

The normal comparison uses supported means and donor medians, not maxima across individual healthy donors. Same-unit source comparability remains under review. Within-tumor ranks also have variable denominators: {counts.measured_proteins.min():,.0f}–{counts.measured_proteins.max():,.0f} measured protein records per sample, from an inventory of 19,893. These limitations are reasons to complete source and coverage review before release, not to tune thresholds until questionable families disappear.

Cap failure removes target eligibility through Other:Expression. Required coverage/identity failures remain Other:Review. Cancer-only qualification lost under a tested rule means no qualifying cancer evidence under that rule; it does not establish absence of expression. Reproductive membership is retained when its own evidence still qualifies.

[All scenarios](threshold-grid.csv) · [All 957 eligible protein records](protein-threshold-evidence.csv.gz) · [Recommended inventory](recommended-eligible-inventory.csv.gz) · [Recommended changes](recommended-category-changes.csv) · [Patient coverage](patient-coverage.csv) · [Sample-specific rank thresholds](sample-rank-cutoffs.csv.gz) · [Summary](summary.json) · [Provenance](manifest.json).

Reproduce in order: `scripts/build_target_sample_ranks.py`, `scripts/audit_target_thresholds.py`, `scripts/render_target_thresholds.py`. The source matrices and candidate matrix hashes must match the pinned rebuild.
"""
    (OUT / "index.md").write_text(paragraphs)
    columns = [
        "id",
        "kind",
        "floor",
        "repeated_fold",
        "single_fold",
        "repeated_rank",
        "single_rank",
        "quantile",
        "disposable_cap_tpm",
        "disposable_cap_ntpm",
        "Cancer Antigens",
        "CTA",
        "Cancer Antigen only",
        "Reproductive Antigen only",
        "Target candidates",
        "current_targets_to_other",
        "published_cancer_nominations",
        "cancer_lost_without_hl_merkel",
        "cancer_lost_without_largest_sample",
    ]
    payload = json.dumps(
        grid[columns].astype(object).where(grid[columns].notna(), None).to_dict("records"),
        allow_nan=False,
    ).replace("<", "\\u003c")
    page = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>OncoRef threshold comparisons</title><style>body{font:15px system-ui;color:#173042;margin:2rem}h1{font-size:1.8rem}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:.6rem;border-bottom:1px solid #d6dfe4;text-align:right}th{position:sticky;top:0;background:#eaf0f3}td:first-child,th:first-child{text-align:left}input,select{padding:.5rem;margin:.5rem}a{color:#08677b}.scroll{overflow:auto;max-height:70vh}.note{max-width:80ch;line-height:1.5}</style><h1>Threshold comparisons</h1><p class="note">Counterfactual rules, with existing normal exclusions unchanged. The current catalog remains 632 CTA / 250 Cancer Antigen only / 28 Reproductive Antigen only. <a href="../thresholds/">Read the recommendation and definitions</a> · <a href="threshold-grid.csv">Download every scenario</a></p><label>Search <input id="search" placeholder="e.g. fold_L1_R1_S10"></label><label>Rule family <select id="kind"><option value="">All</option><option value="fold">TPM and fold enrichment</option><option value="rank">Within-tumor rank</option><option value="cohort">Across-patient percentile</option></select></label><label>Disposable cap <select id="cap"><option value="all">All</option><option value="none">None</option><option>4</option><option>8</option><option>16</option><option>32</option><option>128</option></select></label><p id="count"></p><div class="scroll"><table><thead><tr id="head"></tr></thead><tbody id="body"></tbody></table></div><script>const rows=__DATA__;const cols=Object.keys(rows[0]);const head=document.getElementById('head'),body=document.getElementById('body'),kind=document.getElementById('kind'),cap=document.getElementById('cap'),search=document.getElementById('search'),count=document.getElementById('count');const names={id:'Rule ID',kind:'Rule family',floor:'TPM floor',repeated_fold:'Two-patient fold',single_fold:'Single-tumor fold',repeated_rank:'Two-patient rank',single_rank:'Single-tumor rank',quantile:'Across-patient quantile',disposable_cap_tpm:'Disposable cap (TPM)',disposable_cap_ntpm:'Disposable cap (nTPM)',current_targets_to_other:'Move to Other',published_cancer_nominations:'Published nominations retained',cancer_lost_without_hl_merkel:'Lost without HL/Merkel',cancer_lost_without_largest_sample:'Lost without largest sample'};head.replaceChildren(...cols.map(c=>{const e=document.createElement('th');e.textContent=names[c]||c;return e}));function render(){const selected=rows.filter(r=>(!kind.value||r.kind===kind.value)&&(cap.value==='all'||(cap.value==='none'?r.disposable_cap_tpm===null:r.disposable_cap_tpm===Number(cap.value)))&&r.id.toLowerCase().includes(search.value.toLowerCase()));count.textContent=selected.length+' scenarios';body.replaceChildren(...selected.map(r=>{const tr=document.createElement('tr');for(const c of cols){const td=document.createElement('td');td.textContent=r[c]===null?'—':r[c];tr.append(td)}return tr}));}for(const id of ['kind','cap','search'])document.getElementById(id).addEventListener('input',render);render();</script></html>""".replace(
        "__DATA__", payload
    )
    (OUT / "report.html").write_text(page)
    manifest = {
        "implementation_sha256": {
            p: sha(ROOT / p)
            for p in [
                "scripts/build_target_sample_ranks.py",
                "scripts/audit_target_thresholds.py",
                "scripts/render_target_thresholds.py",
                "oncoref/gene_matching.py",
            ]
        },
        "inputs_sha256": summary["input_sha256"],
        "outputs_sha256": {
            p.name: sha(p)
            for p in sorted(OUT.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    build()
