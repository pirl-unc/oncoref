#!/usr/bin/env python3
"""Build the ordered CTA screen and excluded Variable review set offline."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_clinical_curation_report import read_verified  # noqa: E402
from cta_recalibration_report import markdown_table  # noqa: E402
from cta_scope_report import sha256, write_csv  # noqa: E402

from oncoref.cta_screen import POLICY_ID, reviewed_donor_rna, screen_candidates  # noqa: E402

BASE = ROOT / "docs/audits/cta-curation-20261021-r2"
PREVIOUS = ROOT / "docs/audits/cta-continuous-20261022"


def explorer(out, result):
    data = json.dumps(result.fillna("").to_dict(orient="records"), allow_nan=False).replace(
        "<", "\\u003c"
    )
    template = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CTA screening results</title>
<style>body{font:17px/1.5 system-ui,sans-serif;max-width:1200px;margin:32px auto;padding:0 20px;color:#172b38}a{color:#176665}.note{padding:15px;background:#fff2d9;border-left:4px solid #9b6d18}input,select,button{font:inherit;padding:8px;margin:5px}table{border-collapse:collapse;width:100%}td,th{text-align:left;padding:12px;border-bottom:1px solid #dce3e6;vertical-align:top}th{background:#eef5f4}small{display:block}summary{cursor:pointer}dt{font-weight:600;margin-top:8px}dd{margin:0;overflow-wrap:anywhere}.scroll{overflow-x:auto}.Variable{color:#885800;font-weight:700}.Excluded{color:#9c3434}.Retained{color:#176665;font-weight:700}</style>
<main><p><a href="./">Rules and report</a> · <a href="gene-results.csv">All results</a> · <a href="variable.csv">Variable candidates</a></p><h1>CTA screening results</h1>
<p class="note"><b>Variable means excluded now.</b> Retained means a research screen pass. Individual donor RNA is currently reviewed for CTAG2 only; all other donor profiles remain unassessed.</p>
<label>Gene or ID <input id="search" type="search"></label><label>Category <select id="category"><option value="">All</option><option>Retained</option><option>Variable</option><option>Excluded</option><option>Unassessed</option></select></label><p id="status" aria-live="polite"></p>
<div class="scroll"><table><thead><tr><th>Gene and evidence</th><th>Category</th><th>Donor review</th></tr></thead><tbody id="rows"></tbody></table></div><p><button id="prev">Previous</button><span id="page"></span><button id="next">Next</button></p></main>
<script type="application/json" id="data">__DATA__</script><script>
'use strict';const data=JSON.parse(document.getElementById('data').textContent),search=document.getElementById('search'),category=document.getElementById('category');let page=0;const size=50;
function node(tag,text){const e=document.createElement(tag);e.textContent=String(text);return e}
function update(){const q=search.value.trim().toLowerCase(),filtered=data.filter(r=>(r.Gene+' '+r['Gene ID']).toLowerCase().includes(q)&&(!category.value||r.Category===category.value)),pages=Math.max(1,Math.ceil(filtered.length/size));page=Math.min(page,pages-1);const body=document.getElementById('rows');body.replaceChildren();for(const r of filtered.slice(page*size,(page+1)*size)){const tr=document.createElement('tr'),gene=node('td',r.Gene||r['Gene ID']);gene.append(node('small',r['Gene ID']));const details=document.createElement('details');details.append(node('summary','Reasons and evidence'));const dl=document.createElement('dl');for(const [k,v] of Object.entries(r)){if(v===''||['Gene','Gene ID'].includes(k))continue;dl.append(node('dt',k),node('dd',v))}details.append(dl);gene.append(details);const cat=node('td',r.Category);cat.className=r.Category;tr.append(gene,cat,node('td',r['Donor review']));body.append(tr)}if(!filtered.length){const tr=document.createElement('tr'),td=node('td','No matching genes.');td.colSpan=3;tr.append(td);body.append(tr)}document.getElementById('status').textContent=filtered.length.toLocaleString()+' matching genes';document.getElementById('page').textContent='Page '+(page+1)+' of '+pages;document.getElementById('prev').disabled=page===0;document.getElementById('next').disabled=page===pages-1}
search.addEventListener('input',()=>{page=0;update()});category.addEventListener('change',()=>{page=0;update()});document.getElementById('prev').addEventListener('click',()=>{page--;update()});document.getElementById('next').addEventListener('click',()=>{page++;update()});update();</script></html>"""
    (out / "explore.html").write_text(template.replace("__DATA__", data))


def build(out, dataset=None):
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use an empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    candidates = read_verified(BASE, "candidate-decisions.csv.gz")
    audit = read_verified(BASE, "tissue-evidence.csv.gz")
    previous = read_verified(PREVIOUS, "gene-results.csv")
    clinical = read_verified(BASE, "clinical-target-audit.csv")
    donor_rna = reviewed_donor_rna()
    result, tissues, donor_results = screen_candidates(candidates, audit, donors=donor_rna)
    old_ids = set(previous.loc[previous["Continuous point screen"].eq("Pass"), "Gene ID"])
    result["Previous screen"] = result["Gene ID"].map(
        lambda gid: "Pass" if gid in old_ids else "Does not pass"
    )
    result["Change"] = [
        "Removed"
        if gid in old_ids and not selected
        else "Added"
        if gid not in old_ids and selected
        else "Unchanged"
        for gid, selected in zip(result["Gene ID"], result.Selected)
    ]
    write_csv(result, out / "gene-results.csv")
    write_csv(tissues, out / "tissue-evidence.csv.gz")
    write_csv(donor_results, out / "donor-results.csv")
    write_csv(donor_rna, out / "donor-rna.csv")
    shutil.copyfile(ROOT / "oncoref/data/cta-donor-ihc.csv", out / "donor-ihc.csv")
    write_csv(result.loc[result.Category.eq("Variable")], out / "variable.csv")
    write_csv(result.loc[result.Selected], out / "retained.csv")
    changed = result.loc[result.Change.ne("Unchanged")]
    write_csv(changed, out / "changes.csv")
    clinical_results = clinical[["Symbol", "Ensembl_Gene_ID", "source_url", "limitations"]].merge(
        result, left_on="Ensembl_Gene_ID", right_on="Gene ID", validate="one_to_one"
    )
    write_csv(clinical_results, out / "clinical-results.csv")
    if dataset:
        write_csv(result, dataset)
    counts = result.Category.value_counts().reindex(
        ["Retained", "Variable", "Excluded", "Unassessed"], fill_value=0
    )
    write_csv(pd.DataFrame({"Category": counts.index, "Genes": counts.values}), out / "counts.csv")
    policy = {
        "policy_id": POLICY_ID,
        "basis": "HPA v23 human research screen",
        "negative_ihc_relaxation": "requires reliable reproductive positive in curated gene profile",
        "positive_control_scope": "curated gene profile; antibody/run/specimen matching not established",
        "selectivity_increment": 0.05,
        "allowance_ntpm": 0.1,
        "absolute_ceilings": {
            "Vital": [2, 4, 6, 12],
            "Protected": [6, 8, 12, 24],
            "Manageable": [10, 20, 30, 60],
        },
        "selectivity_minima": {
            "Vital": [0.80, 0.75, 0.70, 0.65],
            "Protected": [0.75, 0.70, 0.65, 0.60],
            "Manageable": [0.70, 0.65, 0.60, 0.55],
        },
        "confidence_order": [
            "Missing/uncertain or no reproductive control",
            "Approved",
            "Supported",
            "Enhanced",
        ],
        "donor_rule": "apply the same absolute, selectivity and R>N rules to each reviewed vital donor",
        "variable": "reported zero and failing donor in same cohort/tissue; complete local negative IHC; no competing hard exclusion/quality flag",
        "variable_selected": False,
        "unreviewed_donors": "Not assessed",
        "donor_model_fitted": False,
        "clinical_safety_probability": None,
    }
    (out / "screening-policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    explorer(out, result)
    spotlight = result.loc[
        result.Gene.isin(
            [
                "CTAG2",
                "ZNF730",
                "CTAG1A",
                "CTAG1B",
                "XAGE1A",
                "XAGE1B",
                "PRAME",
                "MAGEA4",
                "CT83",
                "GAGE1",
                "SAGE1",
                "CTCFL",
                "SSX4",
                "SSX5",
                "MAGEA12",
            ]
        ),
        ["Gene", "Category", "Donor review"],
    ]
    text = f"""# CTA screening: retained, variable and excluded

Updated **23 October 2026**. **{counts["Retained"]} retained**, **{counts["Variable"]} Variable (excluded now)** and **{counts["Excluded"]} other exclusions** among {len(result):,} mapped candidates. {len(changed)} former point-screen passes are removed; no new candidates are admitted.

**Retained is a research-policy pass, not clinical clearance. Only CTAG2 currently has imported individual-donor RNA evidence. All other donor profiles remain Not assessed.**

- **[Find a gene](explore.html)** — filter by category and open its evidence.
- [Simple rules](#the-rules-in-order) · [CTAG2 and Variable](#variable-excluded-now-reconsiderable-later) · [Named targets](#named-targets).
- [All gene results](gene-results.csv) · [Retained set](retained.csv) · [Variable set](variable.csv) · [Changes](changes.csv).
- [Publication intake and unmapped nominations](../cta-curation-20261021-r2/index.md) · [Previous 758-gene pilot](../cta-continuous-20261022/index.md).

## The rules, in order

| Rule | Action |
| --- | --- |
| No reproductive expression, unresolved identity/coverage, family exclusion, or a separate safety hold | Do not select. Preserve the reason. |
| Positive IHC, including Low, in a Vital, Protected or Manageable compartment | Exclude. A positive reproductive control cannot override it. |
| Missing/uncertain normal IHC, incomplete local staining, or no reliable positive reproductive IHC | Use baseline RNA limits. Missing staining is never a negative result. |
| Complete reliable negative local IHC **and** reliable positive reproductive staining | Allow the existing confidence-dependent RNA relaxation. |
| Any screened tissue fails its RNA ceiling, selectivity, or reproductive-higher rule | Exclude. |
| A reviewed individual donor fails those same RNA rules in a Vital tissue | Exclude even when the tissue average passes. Classify the qualifying zero/high, IHC-negative pattern as **Variable**. |

Reproductive positivity is checked at the mapped cell compartment. Disposable expression is ignored. The positive control is currently established **within the curated HPA gene profile**: the export does not prove identical antibody, batch or specimen across every tissue. It does not calibrate antibody sensitivity or a proportional RNA-to-protein relationship.

### The numerical limits

| Tissue category | Absolute nTPM ceilings | Selectivity minima |
| --- | --- | --- |
| Vital | 2 / 4 / 6 / 12 | 80 / 75 / 70 / 65% |
| Protected | 6 / 8 / 12 / 24 | 75 / 70 / 65 / 60% |
| Manageable | 10 / 20 / 30 / 60 | 70 / 65 / 60 / 55% |

Order: baseline → Approved-negative → Supported-negative → Enhanced-negative. Every relaxation requires the reproductive positive control above. Selectivity is tightened by **five percentage points throughout**, as tested against the named panel; absolute ceilings are unchanged. These are research settings, not validated safety limits.

For each tissue, let `R` be the highest reproductive nTPM, `N` its normal nTPM, `C` its absolute ceiling and `s` its selectivity minimum. Require **N ≤ C**, **N ≤ R × (1−s)/s + 0.1**, and **N < R**. The small continuous allowance stays fixed at 0.1 nTPM. Reviewed donor values replace `N` in **exactly these same comparisons**. There is no separate 10 nTPM donor cutoff.

## Variable: excluded now, reconsiderable later

Membership requires:

1. At least one **reported-zero RNA donor** and a different **RNA donor failing the applicable limits**, within the same gene, Vital tissue, cohort, source version and units.
2. Complete reliable negative IHC in the affected tissue's routinely measured compartments; no competing normal-protein exclusion or separate safety hold.
3. Measured reproductive expression and usable candidate identity/coverage.

**Variable never means selected.** It records a reason to revisit the exclusion. No measurement becomes zero because it is missing, and zeros/high values from different cohorts are not combined to manufacture membership. Negative IHC in other volunteers does not clear an RNA-positive donor. A patient with an RNA-zero measurement is not automatically eligible.

### CTAG2 / LAGE-1

- Pinned v23 heart RNA donors: **13.5, 0, 0 and 7.1 nTPM**. The consensus heart estimate is 5.1 nTPM; reproductive maximum is testis 10.1 nTPM.
- Heart has Enhanced-negative cardiomyocyte IHC and testis has positive reproductive staining. The applicable absolute ceiling is **12 nTPM** and the selectivity-derived ceiling is **5.538 nTPM**. Thus **both nonzero donors fail**, while the mean passes.
- Three separately sampled heart IHC donors (women aged 52/62, man aged 37) have Not detected annotations. Their metadata are preserved in [donor IHC](donor-ihc.csv). These live-source records establish the sample mismatch; they do not replace the pinned v23 numerical data. Current portal RNA values differ slightly (13.1/0/0/6.9) and are not silently mixed into this run.
- Result: **Variable, excluded from the retained set**. This is an unresolved donor-pattern exclusion, not a conclusion that CTAG2 has caused cardiac injury.

Before reconsideration: identify the expressing cell types and isoforms, obtain paired RNA/protein evidence in RNA-positive and RNA-zero specimens, investigate donor/physiological/technical determinants, and independently validate any proposed subgroup. Additional relevant antigen-presentation evidence may inform that review. A new explicit review is required; clinical nomination never overrides the hold.

HPA's **curated** annotations outside testis are negative, but its primary antibody data include staining that HPA disregarded as presumed off-target. Keep the raw staining and curated interpretation separate. The donor rule here uses the pinned curated layer and does not pretend that every raw image is blank. [Primary staining](https://www.proteinatlas.org/ENSG00000126890-CTAG2/tissue/primary%2Bdata) · [Antibody validation](https://www.proteinatlas.org/ENSG00000126890-CTAG2/summary/antibody) · [Heart evidence review](../ctag2-heart-evidence.md).

## Named targets

{markdown_table(spotlight)}

CTAG2 is the only loss from the previously retained 22-gene clinical/requested benchmark. **The other 21 remain retained**, including XAGE1A/B, CTAG1A/B, PRAME and MAGEA4. ZNF730 fails in cerebral cortex, bone marrow, lymph node and thymus. The 17 other removals come from the tighter selectivity rules. [Every changed gene and reason](changes.csv) · [Clinical panel, including previously excluded targets](clinical-results.csv).

The earlier frozen core/extended CTA APIs already exclude CTAG2. Historical reports and sequence exports retain their dated policies; the new explicit `oncoref.cta_screen` selection excludes Variable. This update does not certify reproductive-specific proteins or regenerate k-mer backgrounds.

## Downloads and reproduction

- [Gene results](gene-results.csv) · [Tissue evidence](tissue-evidence.csv.gz) · [Counts](counts.csv).
- [Donor decisions](donor-results.csv) · [Pinned donor RNA](donor-rna.csv) · [IHC donor metadata](donor-ihc.csv).
- [Policy](screening-policy.json) · [Run manifest and checksums](run-manifest.json).

```bash
python scripts/cta_screen_report.py --output /tmp/cta-screen-audit
```

```python
from oncoref import cta_screen
cta_screen.results()                 # full decisions, provenance and reasons
cta_screen.retained_gene_ids()       # all 740 loci; excludes Variable
cta_screen.retained_gene_names()     # known symbols; unnamed loci remain in IDs
cta_screen.variable_gene_names()     # {{"CTAG2"}}
```

No donor model is fitted, no false-negative probabilities are invented, and an unreviewed gene is never labelled donor-negative. Human nTPM ceilings are not applied to dog/mouse measurements. [API and limitations](../../cta-screen.md).
"""
    (out / "index.md").write_text(text)
    inputs = [
        BASE / "run-manifest.json",
        BASE / "candidate-decisions.csv.gz",
        BASE / "tissue-evidence.csv.gz",
        BASE / "clinical-target-audit.csv",
        PREVIOUS / "run-manifest.json",
        PREVIOUS / "gene-results.csv",
        ROOT / "docs/audits/normal-expression-categories-20261021-r2/run-manifest.json",
    ]
    inputs += [
        ROOT / "oncoref/data" / name
        for name in [
            "cta-donor-rna.csv",
            "cta-donor-ihc.csv",
            "cta-reviewed-evidence.csv",
            "normal-tissue-taxonomy.yaml",
            "normal-tissue-mappings.csv",
        ]
    ]
    implementation = [
        Path(__file__),
        ROOT / "oncoref/cta_screen.py",
        ROOT / "oncoref/cta_continuous_screen.py",
        ROOT / "oncoref/normal_tissues.py",
    ]
    manifest = {
        "policy_id": POLICY_ID,
        "donor_model_fitted": False,
        "donor_reviewed_genes": sorted(set(donor_rna.Symbol)),
        "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
        "implementation_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in implementation},
        "outputs": {p.name: sha256(p) for p in sorted(out.iterdir())},
        "published_dataset_sha256": sha256(dataset) if dataset else None,
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"counts": counts.to_dict(), "removed": len(changed), "output": str(out)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--dataset", type=Path, help="Also publish the result as the package dataset"
    )
    args = parser.parse_args()
    build(args.output, args.dataset)
