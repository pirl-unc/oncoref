#!/usr/bin/env python3
"""Reproduce the continuous CTA pilot using checksum-pinned published evidence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_clinical_curation_report import read_verified  # noqa: E402
from cta_curation_report import readable_reason  # noqa: E402
from cta_recalibration_report import markdown_table  # noqa: E402
from cta_scope_report import sha256, write_csv  # noqa: E402

from oncoref import figure_style  # noqa: E402
from oncoref.cta_continuous_screen import (  # noqa: E402
    ContinuousPolicy,
    assess_joint_expression,
    continuous_tissue_screen,
)
from oncoref.normal_tissues import (  # noqa: E402
    SCREENED_NORMAL_TISSUE_TIERS,
    normal_tissue_screening_tier,
    resolve_normal_tissue,
)

BASE = ROOT / "docs/audits/cta-curation-20261021-r2"
CATEGORIES = ROOT / "docs/audits/normal-expression-categories-20261021-r2"
NUMERICAL_REASONS = {"rna_above_tissue_ceiling", "reproductive_selectivity_below_threshold"}


def eligible_before_rna(candidates):
    other_ok = candidates.failure_reasons.map(
        lambda x: not (set(filter(None, x.split(";"))) - NUMERICAL_REASONS)
    )
    return (
        other_ok
        & candidates.quality_flags.eq("")
        & candidates.reproductive_positive
        & candidates.family_eligible
        & candidates.nomination_review_eligible
        & ~candidates.external_safety_hold
    )


def sensitivity_counts(candidates, tissue_audit):
    base = set(candidates.loc[eligible_before_rna(candidates), "Ensembl_Gene_ID"])
    old = set(candidates.loc[candidates.retained, "Ensembl_Gene_ID"])
    rows = []
    for delta in [0, 0.05, 0.1, 0.25, 0.5, 1]:
        result = continuous_tissue_screen(
            tissue_audit, policy=ContinuousPolicy(allowance_ntpm=delta)
        )
        tested = result.loc[result.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)]
        fail = set(tested.loc[~tested.point_result.eq("Pass"), "Ensembl_Gene_ID"])
        passed = base - fail
        rows.append(
            {
                "Allowance (nTPM)": delta,
                "Pass point screen": len(passed),
                "Gained": len(passed - old),
                "Lost": len(old - passed),
            }
        )
    return pd.DataFrame(rows)


def panels():
    manifest = json.loads((CATEGORIES / "run-manifest.json").read_text())
    normal, reproductive = [], []
    for tissue, _ in manifest["routine_contexts"]["bulk_rna"]:
        mapped = resolve_normal_tissue(
            tissue, species="human", source_id="hpa_rna_consensus", source_version="v23"
        )
        if "reproductive" in mapped.scopes:
            reproductive.append(tissue)
        elif (
            normal_tissue_screening_tier(mapped, policy="ihc_adaptive_cta")
            in SCREENED_NORMAL_TISSUE_TIERS
        ):
            normal.append(tissue)
    return normal, reproductive


def candidate_results(candidates, tissue_audit, policy):
    normal, reproductive = panels()
    annotated = continuous_tissue_screen(tissue_audit, policy=policy)
    groups = dict(iter(tissue_audit.groupby("Ensembl_Gene_ID", sort=False)))
    current = dict(iter(annotated.groupby("Ensembl_Gene_ID", sort=False)))
    eligible = dict(zip(candidates.Ensembl_Gene_ID, eligible_before_rna(candidates)))
    results = []
    for gene in candidates.itertuples():
        rows = current.get(gene.Ensembl_Gene_ID, annotated.iloc[:0])
        tested = rows.loc[rows.screening_tier.isin(SCREENED_NORMAL_TISSUE_TIERS)]
        point_pass = bool(
            eligible[gene.Ensembl_Gene_ID] and len(tested) and tested.point_result.eq("Pass").all()
        )
        hard, review = [], list(filter(None, gene.quality_flags.split(";")))
        if gene.external_safety_hold:
            hard.append("external_safety_hold")
        if "positive_ihc_outside_permitted" in gene.failure_reasons.split(";"):
            hard.append("positive_ihc")
        if not gene.family_eligible:
            hard.append("family_exclusion")
        if not gene.nomination_review_eligible:
            review.append("nomination_review_unresolved")
        if not gene.reproductive_positive:
            review.append("no_measured_reproductive_expression")
        if rows.empty:
            assessment = {
                "result": "Excluded" if hard else "Review",
                "reason": ";".join(sorted(set(hard or [*review, "missing_rna_evidence"]))),
                "probability_mean_rules": None,
                "probability_with_donor_tail": None,
                "posterior_draws": 0,
            }
        else:
            assessment = assess_joint_expression(
                groups[gene.Ensembl_Gene_ID],
                normal_tissues=normal,
                reproductive_tissues=reproductive,
                policy=policy,
                hard_exclusions=hard,
                review_reasons=review,
                other_gates_passed=bool(eligible[gene.Ensembl_Gene_ID]),
            )
        findings = "; ".join(
            f"{r.tissue}: {r.rna_ntpm:g} nTPM; {r.point_reasons.replace('_', ' ')}"
            for r in tested.loc[~tested.point_result.eq("Pass")].itertuples()
        )
        unchanged = set(filter(None, gene.failure_reasons.split(";"))) - NUMERICAL_REASONS
        if not gene.family_eligible:
            unchanged.add("family_exclusion")
        if not gene.nomination_review_eligible:
            unchanged.add("existing_nomination_or_manual_review")
        if not gene.reproductive_positive:
            unchanged.add("insufficient_positive_reproductive_evidence")
        reason = "; ".join(
            filter(
                None,
                [
                    readable_reason(";".join(sorted(unchanged))),
                    readable_reason(gene.quality_flags),
                    findings,
                ],
            )
        )
        results.append(
            {
                "Gene": gene.Symbol,
                "Gene ID": gene.Ensembl_Gene_ID,
                "Previous point screen": "Pass" if gene.retained else "Does not pass",
                "Continuous point screen": "Pass" if point_pass else "Does not pass",
                "Point findings": reason or "All point-estimate screening gates pass",
                "Uncertainty result": assessment["result"],
                "Uncertainty reasons": assessment["reason"].replace("_", " ").replace(";", "; "),
                "Joint mean-rule probability": assessment["probability_mean_rules"],
                "Joint probability including donor tail": assessment["probability_with_donor_tail"],
                "Posterior draws": assessment["posterior_draws"],
                "Highest reproductive RNA (nTPM)": gene.reproductive_rna_max_ntpm,
                "Clinical evidence": getattr(gene, "Clinical_evidence", ""),
                "Change": "Gained"
                if point_pass and not gene.retained
                else "Lost"
                if gene.retained and not point_pass
                else "Unchanged",
            }
        )
    result = pd.DataFrame(results)
    result["Clinical evidence"] = result["Gene ID"].map(
        candidates.set_index("Ensembl_Gene_ID")["Clinical evidence"]
    )
    return result.sort_values(["Gene", "Gene ID"]), annotated


def figures(out, sensitivity, allowance):
    import matplotlib.pyplot as plt
    import numpy as np

    figure_style.apply()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, xmax in zip(axes, [2, 10]):
        r = np.linspace(0, xmax, 500)
        old = np.minimum(r / 3, 2)
        new = np.minimum.reduce([r, r / 3 + allowance, np.full_like(r, 2)])
        ax.fill_between(r, 0, new, color="#2c7b76", alpha=0.12)
        ax.plot(r, old, color="#777777", ls="--", label="Ratio only")
        ax.plot(r, new, color="#2c7b76", label=f"Continuous allowance: {allowance:g} nTPM")
        ax.set(xlabel="Highest reproductive RNA (nTPM)", ylabel="Maximum normal RNA (nTPM)")
        ax.set_title("Low-expression detail" if xmax == 2 else "Absolute ceiling remains binding")
    axes[0].legend(fontsize=8)
    fig.suptitle("Example: Vital tissue, no usable IHC — 75% selectivity and 2 nTPM ceiling")
    fig.tight_layout()
    for ext in ["png", "svg"]:
        figure_style.save(fig, out / f"continuous-boundary.{ext}", keep_titles=True)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(
        sensitivity["Allowance (nTPM)"],
        sensitivity["Pass point screen"],
        marker="o",
        color="#2c7b76",
    )
    for row in sensitivity.itertuples(index=False, name=None):
        ax.annotate(
            str(row[1]), (row[0], row[1]), xytext=(0, 8), textcoords="offset points", ha="center"
        )
    ax.set(
        xlabel="Fixed allowance (nTPM)",
        ylabel="Candidates passing the point screen",
        title="Allowance sensitivity — not uncertainty-qualified inclusion",
    )
    ax.margins(y=0.2)
    fig.tight_layout()
    for ext in ["png", "svg"]:
        figure_style.save(fig, out / f"allowance-sensitivity.{ext}", keep_titles=True)
    plt.close(fig)


def explorer(out, results):
    # Match the existing report's plain searchable table and expandable evidence.
    records = results.fillna("").to_dict(orient="records")
    body = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Continuous CTA screen — gene results</title>
<style>body{font:17px/1.5 system-ui,sans-serif;color:#172b38;max-width:1200px;margin:32px auto;padding:0 20px}a{color:#176665}input,select,button{font:inherit;padding:8px;margin:5px}table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:12px;border-bottom:1px solid #dce3e6;vertical-align:top}th{background:#eff5f5}dl{max-width:650px}dt{font-weight:600;margin-top:12px}dd{margin:0;overflow-wrap:anywhere}.scroll{overflow-x:auto}summary{cursor:pointer}small{display:block}</style>
<main><p><a href="./">Report and methods</a> · <a href="gene-results.csv">Download all results</a></p><h1>Continuous CTA screen</h1><p>Point-screen passes and uncertainty decisions are separate. No donor model has been fitted for this report; Review does not mean known unsafe.</p>
<label>Gene or ID <input type="search" id="search"></label><label>Point screen <select id="point"><option value="">All</option><option>Pass</option><option>Does not pass</option></select></label><label>Uncertainty <select id="decision"><option value="">All</option><option>Included</option><option>Excluded</option><option>Review</option></select></label><label><input type="checkbox" id="gained"> New point-screen passes</label><p id="status" aria-live="polite"></p>
<div class="scroll"><table><thead><tr><th>Gene / evidence</th><th>Point screen</th><th>Uncertainty</th><th>Change</th></tr></thead><tbody id="rows"></tbody></table></div><p><button id="prev">Previous</button><span id="page"></span><button id="next">Next</button></p></main>
<script type="application/json" id="data">__DATA__</script><script>
'use strict';const data=JSON.parse(document.getElementById('data').textContent);const controls=Object.fromEntries(['search','point','decision','gained'].map(id=>[id,document.getElementById(id)]));let page=0;const size=50;
function node(tag,text){const e=document.createElement(tag);e.textContent=String(text);return e}
function update(){const q=controls.search.value.trim().toLowerCase();const rows=data.filter(r=>(r.Gene+' '+r['Gene ID']).toLowerCase().includes(q)&&(!controls.point.value||r['Continuous point screen']===controls.point.value)&&(!controls.decision.value||r['Uncertainty result']===controls.decision.value)&&(!controls.gained.checked||r.Change==='Gained'));const pages=Math.max(1,Math.ceil(rows.length/size));page=Math.min(page,pages-1);const body=document.getElementById('rows');body.replaceChildren();for(const r of rows.slice(page*size,(page+1)*size)){const tr=document.createElement('tr');const gene=node('td',r.Gene);gene.append(node('small',r['Gene ID']));const details=document.createElement('details');details.append(node('summary','Evidence and reasons'));const dl=document.createElement('dl');for(const [k,v] of Object.entries(r)){if(v===''||['Gene','Gene ID'].includes(k))continue;dl.append(node('dt',k),node('dd',v))}details.append(dl);gene.append(details);tr.append(gene,node('td',r['Continuous point screen']),node('td',r['Uncertainty result']),node('td',r.Change));body.append(tr)}if(!rows.length){const tr=document.createElement('tr');const td=node('td','No matching genes. Clear a filter or try a shorter name.');td.colSpan=4;tr.append(td);body.append(tr)}document.getElementById('status').textContent=rows.length.toLocaleString()+' matching genes';document.getElementById('page').textContent='Page '+(page+1)+' of '+pages;document.getElementById('prev').disabled=page===0;document.getElementById('next').disabled=page===pages-1}
for(const c of Object.values(controls))c.addEventListener(c.type==='search'?'input':'change',()=>{page=0;update()});document.getElementById('prev').addEventListener('click',()=>{page--;update()});document.getElementById('next').addEventListener('click',()=>{page++;update()});update();
</script></html>"""
    (out / "explore.html").write_text(
        body.replace("__DATA__", json.dumps(records, allow_nan=False).replace("<", "\\u003c"))
    )


def build(out, allowance=0.1):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    policy = ContinuousPolicy(allowance_ntpm=allowance)
    candidates = read_verified(BASE, "candidate-decisions.csv.gz")
    audit = read_verified(BASE, "tissue-evidence.csv.gz")
    clinical = read_verified(BASE, "clinical-target-audit.csv")
    result, evidence = candidate_results(candidates, audit, policy)
    sensitivity = sensitivity_counts(candidates, audit)
    write_csv(result, out / "gene-results.csv")
    write_csv(evidence, out / "tissue-evidence.csv.gz")
    write_csv(result.loc[result.Change.ne("Unchanged")], out / "changes.csv")
    write_csv(sensitivity, out / "sensitivity.csv")
    clinical_result = clinical[
        ["Symbol", "Ensembl_Gene_ID", "Clinical evidence", "source_url", "limitations"]
    ].merge(
        result.drop(columns="Clinical evidence"),
        left_on="Ensembl_Gene_ID",
        right_on="Gene ID",
        validate="one_to_one",
    )
    write_csv(clinical_result, out / "clinical-results.csv")
    counts = (
        result["Uncertainty result"]
        .value_counts()
        .reindex(["Included", "Excluded", "Review"], fill_value=0)
    )
    summary = pd.DataFrame({"Uncertainty result": counts.index, "Genes": counts.values})
    write_csv(summary, out / "uncertainty-summary.csv")
    normal, reproductive = panels()
    base_manifest = json.loads((BASE / "run-manifest.json").read_text())
    base_policy = json.loads((BASE / "screening-policy.json").read_text())
    if sha256(BASE / "screening-policy.json") != base_manifest["outputs"]["screening-policy.json"]:
        raise ValueError("Base policy checksum mismatch")
    (out / "screening-policy.json").write_text(
        json.dumps(
            {
                **policy.manifest(),
                "base_policy": base_policy,
                "normal_tissues": normal,
                "reproductive_tissues": reproductive,
                "donor_model_fitted": False,
                "uncertainty_input_coverage": 0,
            },
            indent=2,
        )
        + "\n"
    )
    figures(out, sensitivity, allowance)
    explorer(out, result)
    passed = int(result["Continuous point screen"].eq("Pass").sum())
    changes = result.loc[
        result.Change.eq("Gained"),
        ["Gene", "Highest reproductive RNA (nTPM)", "Uncertainty result"],
    ]
    spotlight = result.loc[
        result.Gene.isin(
            [
                "CT83",
                "CTAG1A",
                "CTAG1B",
                "CTAG2",
                "SAGE1",
                "CTCFL",
                "GAGE1",
                "SSX4",
                "SSX5",
                "XAGE1A",
                "XAGE1B",
                "PRAME",
                "MAGEA3",
                "MAGEA12",
                "SPAG9",
            ]
        ),
        ["Gene", "Continuous point screen", "Uncertainty result"],
    ]
    text = f"""# Continuous CTA expression screen

Updated **22 October 2026**. **{passed:,} candidates pass the point-estimate screen** ({int(result.Change.eq("Gained").sum())} gained, {int(result.Change.eq("Lost").sum())} lost versus the published 746). **No candidate has an uncertainty-qualified inclusion yet:** the pinned input bundle contains tissue summaries, not a fitted donor model.

## Start here

- **[Find a gene](explore.html)** — search, filter, and open its reasons.
- **[What changed](#what-changed)** · **[Named CTAs](#named-ctas)** · **[Filtering rules](#filtering-rules)**.
- **[Uncertainty results](#uncertainty-results)** explains why a point-screen pass can still be Review.
- **[Downloads](#downloads-and-reproduction)** includes every gene, every measured tissue and checksums.
- [Publication sources and full intake](../cta-curation-20261021-r2/index.md#publication-sources) · [Existing disjoint gene/protein sets](../normal-sequence-sets-20261021/index.md).

## What changed

We replaced the percentage cliff with a continuous expression margin, and required reproductive RNA to be strictly higher than each screened normal tissue. The new allowance is **{allowance:g} nTPM**, fixed across genes. It is a provisional project parameter, **not an assay detection limit, measurement error or a confidence interval**. No gene-specific rescue is used.

The published absolute ceilings, local IHC requirements, positive-IHC veto, reproductive admission, coverage checks, family reviews and external holds are preserved. The previously discussed tighter 730-gene absolute grid is not silently substituted. Historical gene/protein exports remain tied to their original policy; this pilot does not create new k-mer exclusion sets from Review candidates.

### Newly passing the point screen

{markdown_table(changes)}

These candidates need evidence review; passing revised parameters does not establish tumor antigenicity or clinical safety. [Full changes](changes.csv).

## Filtering rules

For each Vital, Protected or Manageable tissue:

```text
R = highest reproductive tissue RNA
N = this normal tissue's RNA
C = its existing category/IHC absolute ceiling
s = its existing category/IHC selectivity minimum

N <= C
N <= ((1 - s) / s) * R + allowance
R > N
```

All three must hold. The first two boundaries are inclusive; equality of reproductive and normal expression fails the third. Disposable tissue is ignored for these comparisons and never contributes to R. Reproductive admission remains RNA >0 or positive reproductive IHC; the RNA comparison needs a positive reproductive RNA reference. Missing measurements are not zero. Any positive normal IHC, including Low, still excludes.

| Tissue | Absolute ceilings, nTPM | Selectivity minima |
| --- | --- | --- |
| Vital | 2 / 4 / 6 / 12 | 75 / 70 / 65 / 60% |
| Protected | 6 / 8 / 12 / 24 | 70 / 65 / 60 / 55% |
| Manageable | 10 / 20 / 30 / 60 | 65 / 60 / 55 / 50% |

Order: missing/uncertain IHC → Approved negative → Supported negative → Enhanced negative. Relaxation still requires complete routine negative IHC in the same tissue, using its weakest confidence. These are research parameters, not validated safety limits.

![Continuous expression boundary](continuous-boundary.png)

The example shows a 75% minimum, equivalent to N ≤ R/3 before the allowance. The shaded region also enforces the absolute ceiling and R > N; the equality boundary R = N is excluded. At low expression the allowance matters more, without a jump at 1 nTPM. At high expression the original ratio or absolute ceiling dominates.

### Allowance sensitivity

{markdown_table(sensitivity)}

![Allowance sensitivity](allowance-sensitivity.png)

Zero allowance passes 747 candidates: the old 746 plus **OR2A4**, whose exact 0.3 / (0.3 + 0.1) = 75% boundary previously failed because of floating-point rounding. Inclusive boundaries now accommodate one representable numerical step; this is separate from the biological allowance. Thus 1 of this run's {int(result.Change.eq("Gained").sum())} new passes is a numerical correction and {int(result.Change.eq("Gained").sum()) - 1} come from the {allowance:g} allowance. The starting allowance was selected as a small, inspectable pilot tolerance, not fitted to clinical targets. This table is sensitivity analysis, not validation. Independent donor data and controls are needed to calibrate it.

## Named CTAs

{markdown_table(spotlight)}

All 16 previously retained loci in the clinical-evidence panel still pass the point screen. A trial registration is evidence of a nominated treatment target, not proof of safety. [All 30 clinical-panel loci, references and limitations](clinical-results.csv).

**CTAG2:** pinned reproductive maximum 10.1 and heart 5.1 nTPM still pass the point rules. Existing donor evidence reports two positive HPA heart samples and two zeros; that heterogeneity is a reason to fit a donor model, not assign a narrow interval around 5.1. Its mean-rule pass is not uncertainty-qualified inclusion. [Detailed donor/clinical review](../ctag2-heart-evidence.md).

**BAGE and GAGE2B** remain source nominations without an assessed canonical RNA row. This change supplies neither missing mappings nor donor measurements. **SPAG9** remains assessed with its original counterevidence. [Identity audit](../cta-curation-20261021-r2/index.md#bage-and-gage2b-what-is-missing).

## Uncertainty results

{markdown_table(summary)}

This second result is deliberately separate from the point screen. Excluded currently reflects explicit protein, family or safety exclusions. Review covers unavailable donor uncertainty and other evidence gaps. RNA point failures alone do not become statistically confident failures without an uncertainty model.

The implementation evaluates **joint posterior draws** for all {len(normal)} screened normal tissues and {len(reproductive)} reproductive reference tissues in the independent HPA panel:

1. Recompute the reproductive maximum in each draw, keeping tissue correlations and reference-selection uncertainty.
2. Test all three expression rules across all normal tissues in that draw. The passing fraction is the joint mean-rule probability; marginal probabilities are never multiplied.
3. Separately check posterior draws of each tissue's **95th donor percentile** against its absolute ceiling. This prevents a precise mean from concealing donor variability. Tail concerns yield Review rather than a clinical-toxicity claim.
4. **Included:** at least 95% joint support including the donor-tail check, with the other gates satisfied. **Excluded:** at least 95% support for failing the mean rules, or a hard protein/curation/safety exclusion. **Review:** everything else, including missing or degenerate uncertainty inputs.

The 95% decision probability and donor percentile are explicit provisional modeling choices. The result is conditional on the model and measured panel, not the probability that treatment is safe. At least 1,000 independent posterior draws and passed model checks are required. IHC labels are not converted into invented antibody sensitivities or numeric probabilities.

### What remains unavailable

**No donor model was fitted in this run.** All posterior probability fields are blank, not zero. The imported HPA v23 consensus uses aggregated expression and the maximum across its underlying sources; these are not biological replicates from which a sampling interval can be recovered. [HPA v23 download definitions](https://v23.proteinatlas.org/about/download).

HPA also describes sample-level transcript TPM and estimated counts, which could support a future fit. Their normalization, gene aggregation, donor identities and coverage must be reconciled with the pinned gene-level nTPM scale first. Counts/TPM from another release or GTEx cannot simply be treated as HPA v23 nTPM. [HPA sample data](https://v24.proteinatlas.org/humanproteome/tissue/data) · [Expression units](../../cta-expression-units.md).

A count-aware model can estimate and shrink unstable low-count effects; the published apeglm method provides relevant statistical precedent. It does not validate this project's tolerances or turn tissue summaries into confidence intervals. [Zhu et al., 2019](https://pubmed.ncbi.nlm.nih.gov/30395178/). A future fit also needs cohort heterogeneity, donor variation, zero-count uncertainty, mapping ambiguity and held-out calibration. Repeating consensus values or bootstrapping only observed zeros is explicitly rejected as sufficient uncertainty evidence.

## Downloads and reproduction

- [Gene results](gene-results.csv) · [Changes](changes.csv) · [Clinical evidence/results](clinical-results.csv).
- [All tissue margins and IHC context](tissue-evidence.csv.gz) · [Allowance sensitivity](sensitivity.csv).
- [Uncertainty counts](uncertainty-summary.csv) · [Policy and independent tissue panels](screening-policy.json).
- [Run manifest and hashes](run-manifest.json) · [Boundary SVG](continuous-boundary.svg) · [Sensitivity SVG](allowance-sensitivity.svg).

```bash
python scripts/cta_continuous_report.py --output /tmp/cta-continuous-audit
```

The input archive, output tables and implementation are hashed. The public Python API is documented in [continuous CTA screening](../../cta-continuous-screen.md). Frozen reproductive/normal sequence backgrounds remain unchanged until a reviewed policy explicitly regenerates them.
"""
    (out / "index.md").write_text(text)
    inputs = [
        BASE / x
        for x in [
            "run-manifest.json",
            "candidate-decisions.csv.gz",
            "tissue-evidence.csv.gz",
            "clinical-target-audit.csv",
            "screening-policy.json",
        ]
    ]
    inputs += [CATEGORIES / "run-manifest.json", ROOT / "oncoref/data/normal-tissue-taxonomy.yaml"]
    implementation = [Path(__file__), ROOT / "oncoref/cta_continuous_screen.py"]
    manifest = {
        "schema": "oncoref.cta-continuous-report.v1",
        **policy.manifest(),
        "sources": base_manifest["sources"],
        "donor_model_fitted": False,
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "implementation_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in implementation},
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir())},
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"point_pass": passed, "uncertainty": counts.to_dict(), "output": str(out)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allowance", type=float, default=0.1)
    args = parser.parse_args()
    build(args.output, allowance=args.allowance)


if __name__ == "__main__":
    main()
