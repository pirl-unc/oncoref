# OncoRef 1.24.0: TCR-T clinical evidence and serious events

The human selection remains **522 CTA proteins (558 member gene loci), 18,031
Unsafe, 20 Neither and 1,290 unassigned**. There are **no category changes** from
the pinned 1.23 selection. This release adds a reviewed clinical-evidence layer
and shows how that evidence changes individual RNA ceilings.

[Full curation report and funnel](audits/cta-source-walkthrough-20261028/index.md) ·
[TCR-T serious-event review](audits/cta-tcr-trials-20261029/index.md) ·
[Nine-page clinical PDF](audits/cta-tcr-trials-20261029/cta-tcr-safety-1.24.0.pdf) ·
[Logic slides](audits/cta-source-walkthrough-20261028/cta-selection-1.24.0.pdf) ·
[Rules and API](cta-selection.md)

## What the clinical review covers

Eight TCR-T products, nine treated cohorts and eleven sentinel events are linked
to their source papers, exact epitopes, HLA restrictions and protein origins.
Fifty-five peptide/isoform links were checked against the pinned Ensembl 115
protein sequences. Twelve original paper/protocol artifacts use the existing
checksum-pinned reference-data API and cache layout.

This is a targeted review, not an exhaustive list of trials or all SAEs. SAE
seriousness and CTCAE severity are recorded separately. The report distinguishes
treatment-related injury from evidence that a receptor recognized healthy tissue:

- The MAGE-A3 a3a cardiac deaths implicated unintended recognition of titin.
- The separate NCI MAGE-A3/A9 receptor neurological cases implicated the
  cross-recognized MAGE-A12 peptide.
- c259's late cardiac death and Guillain–Barré cases remain explicit caveats
  without an established peptide-specific injury mechanism.
- TBI-1301's grade 3 lung injury retains unresolved specificity and blocks
  credit for that product.

## How the RNA screen changes

A reviewed treated cohort without an identified grade 4–5 specificity-related
signal can support **one existing RNA ladder step**, capped at Enhanced-negative.
Trials do not stack. Missing or unresolved attribution does not count as a
favorable result. Protein and clinical-injury exclusions remain independent.

| Tissue tier | Baseline RNA ceilings | With qualifying clinical credit |
|---|---|---|
| Vital | 4 / 8 / 16 / 32 | 8 / 16 / 32 / 32 |
| Protected | 8 / 16 / 32 / 64 | 16 / 32 / 64 / 64 |
| Manageable | 16 / 32 / 64 / 128 | 32 / 64 / 128 / 128 |

Each ladder follows missing/unusable, Approved-negative, Supported-negative,
Enhanced-negative IHC. Values are nominal TPM/nTPM research parameters, not
measured clinical toxicity limits; corneal CP10K keeps its separate /100 scale.
The credit does not invent or upgrade IHC evidence, waive coverage, or establish
cancer or Reproductive expression.

Five products currently qualify: c259/lete-cel, IMA203, afami-cel, IMA202 and
ADP-A2M10. Their evidence covers six protein profiles and seven member gene loci.
Every member of a collapsed protein must be supported; there is no family,
alias-prefix or cross-species transfer. Exact peptide/product limits remain visible.

**CTAG2 was already CTA.** HPA v23 heart RNA is 5.1 nTPM with Enhanced-negative
IHC and a 32 nTPM ceiling. The new credit cannot raise an already capped ceiling.
c259 recognizes SLLMWITQC shared by NY-ESO-1 and LAGE-1A; this does not establish
safety for every CTAG2 isoform, epitope or receptor.

## API and reproducibility

```python
from oncoref import tcr_trials

tcr_trials.products()  # product-level credit, exact targets and limitations
tcr_trials.cohorts()   # treated denominators and reviewed attribution
tcr_trials.events()    # sentinel events; unknown grades remain missing
tcr_trials.protein_credit(["ENSG00000126890"])  # CTAG2
```

`oncoref tcr-trials --table all` exports the complete ledger offline. Existing
selection explanations now show baseline and effective RNA ceilings and the
supporting product IDs. Pin package **1.24.0**, policy
**`oncoref.cta-selection.2026-10-29.4`**, and the returned source snapshot.

The [clinical report](audits/cta-tcr-trials-20261029/index.md) includes the
baseline snapshot, exact membership comparison, context-level threshold effects,
source pins and reproduction scripts. These are research classifications with
epitope-specific clinical evidence, not clinical clearance of whole proteins.
