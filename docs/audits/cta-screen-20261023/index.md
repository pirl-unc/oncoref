# CTA screening: retained, variable and excluded

Updated **23 October 2026**. **740 retained**, **1 Variable (excluded now)** and **1799 other exclusions** among 2,540 mapped candidates. 18 former point-screen passes are removed; no new candidates are admitted.

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

| Gene | Category | Donor review |
| --- | --- | --- |
| CT83 | Retained | Not assessed |
| CTAG1A | Retained | Not assessed |
| CTAG1B | Retained | Not assessed |
| CTAG2 | Variable | Reviewed |
| CTCFL | Retained | Not assessed |
| GAGE1 | Retained | Not assessed |
| MAGEA12 | Excluded | Not assessed |
| MAGEA4 | Retained | Not assessed |
| PRAME | Retained | Not assessed |
| SAGE1 | Retained | Not assessed |
| SSX4 | Retained | Not assessed |
| SSX5 | Retained | Not assessed |
| XAGE1A | Retained | Not assessed |
| XAGE1B | Retained | Not assessed |
| ZNF730 | Excluded | Not assessed |

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
cta_screen.variable_gene_names()     # {"CTAG2"}
```

No donor model is fitted, no false-negative probabilities are invented, and an unreviewed gene is never labelled donor-negative. Human nTPM ceilings are not applied to dog/mouse measurements. [API and limitations](../../cta-screen.md).
