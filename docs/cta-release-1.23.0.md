# OncoRef 1.23.0: reproductive concentration and a common tissue panel

The human selection contains **522 CTA proteins (558 member gene loci), 18,031
Unsafe, 20 Neither and 1,290 unassigned** across 19,863 proteins. PRAME, MAGE-A4,
NY-ESO-1, XAGE1A/B and CTAG2 remain CTA under the same general rules. These are
research classifications, not demonstrations of clinical safety.

[Full report and funnel](audits/cta-source-walkthrough-20261028/index.md) ·
[30-slide PDF](audits/cta-source-walkthrough-20261028/cta-selection-1.23.0.pdf) ·
[Diffuseness report PDF](audits/cta-unassigned-20261029/cta-diffuseness-1.23.0.pdf) ·
[ZNF730 calculation tables](audits/cta-unassigned-20261029/znf730.md) ·
[Rules and API](cta-selection.md)

## What changed

| Outcome | 1.22 human selection | 1.23.0 | Explanation |
|---|---:|---:|---|
| CTA | 546 | 522 | 23 lack the required GTEx core; MTCP1 has a diffuse reproductive profile |
| Unsafe | 18,031 | 18,031 | Exact membership unchanged |
| Neither | 20 | 20 | MTCP1 enters; NPHP3-ACAD11 moves to unassigned for missing GTEx coverage |
| Unassigned | 1,266 | 1,290 | 24 proteins lack the required GTEx core |

The former rule requiring reproductive or developmental expression to exceed the
highest normal value is removed. The new questions are independent:

1. **Normal-tissue safety:** the existing tissue/IHC RNA ceilings and normal-protein
   rules remain. Only this axis can establish Unsafe. Manageable excess blocks CTA
   without automatically making a protein Unsafe.
2. **Cancer expression:** at least one measured cancer sample reaches 2 TPM or HPA
   pTPM. No cancer-versus-normal ratio is imposed.
3. **Reproductive expression:** adult/placental HPA RNA reaches 5 nTPM, or a supported
   embryonic lineage mean reaches 5 TPM with the existing cell-support requirements.
   Reproductive includes placental, extraembryonic, embryonic and fetal branches.
4. **Reproductive concentration:** a qualifying adult/placental profile should be
   much higher than most normal tissues, rather than diffusely similar to them.

Complete **GTEx v11 coverage across 52 native contexts for every member locus** is
required alongside the existing normal bulk, cell and immune panels. Missing
measurements remain missing. DICE and corneal coverage is optional; observed
signals still count. [Coverage losses and unassigned reasons](audits/cta-unassigned-20261029/index.md).

## How diffuseness works

Use **HPA v23 consensus RNA**, with 38 normal contexts representing **29 equally
weighted anatomical units**. For each context compute `N / (N + R)`, where N is
normal RNA and R is the largest Reproductive value in that same matrix. Average
within an organ, then across organs. Disposable tissues count. **Thymus is
Permissive** and is excluded from both the comparator and the reference; its
normal-safety tier is unchanged.

This is an average expression similarity, not the fraction of tissues expressing
RNA. If every normal value is tenfold below R, the score is 9.09%; if every normal
value equals R, it is 50%. A single extreme organ contributes at most 3.45
percentage points. [Every tissue, weight and calculation](audits/cta-unassigned-20261029/diffuseness-impact.md).

A score **above 10%** fails the adult/placental route only when R is at least
5 nTPM and all required HPA comparison/reference measurements are present.
With complete required normal coverage and no unresolved normal-safety hold, a
protein with no remaining qualifying route becomes **Neither**, never Unsafe.
The numeric boundary is a project identity criterion, not a clinical threshold.

This gate alone moves **MTCP1** from CTA to Neither: testis 5.1 nTPM, widespread
normal values around 1–3 nTPM and a 24.53% score. Its embryonic summary is below
the qualifying floor. The same single protein moves at cutoffs from 7.5% through
20%. GTEx is retained as a separate descriptive profile, not a second veto or a
source chosen opportunistically for the better score.

## Explicit limits

Supported embryonic evidence retains an independent route because the available
adult and embryonic panels do not establish a comparable developmental contrast.
Of the 523 CTAs entering the diffuseness step, 154 qualify through embryonic
expression alone. Seven adult-positive proteins lack a complete HPA contrast:
CIMIP3, H2AL1Q, PTTG1IP2, SMIM42, TMEM247, TSPY9 and ZC3H11C. Their required GTEx
coverage is complete; their adult concentration remains **unassessed**, not passed.
Both limitations appear in the public evidence flags.

HPA consensus is a cross-dataset summary, not a matched donor experiment. RNA
summaries cannot establish antigen presentation, absence from rare cells, or
pediatric/pregnancy safety. Protein identities collapse identical selected
reference sequences, not every alternative isoform.

The merged release preserves dog managed protein/peptide datasets and the mouse
P1A native-evidence seed from 1.22.1–1.22.3. Dog's 14 nominations and mouse's single
P1A nomination remain unassigned/unassessed; there is no quantitative native
selection policy and no transfer of human safety evidence through orthology.

## API and downstream migration

```python
from oncoref import cta_selection

proteins = cta_selection.results(species="human", category="CTA")
reasons = cta_selection.audit(species="human").explain("MTCP1")
```

Pin package **1.23.0**, policy **`oncoref.cta-selection.2026-10-29.3`**, and the
returned source snapshot. The result schema is unchanged. Preserve member gene
identities, evidence flags and null categories. Unassigned is not Neither.

Older `cta_genes`, exact-protein backgrounds and
`normal_expression_categories.classify_normal_expression()` remain separately
versioned products. The last API retains **Other → Disposable/Silent** under its
older RNA/IHC policy; current Neither has no dedicated Disposable-only branch.
[The migration gap is tracked in #647](https://github.com/pirl-unc/oncoref/issues/647).
This release does not silently rewrite downstream k-mer backgrounds.

The analysis, PDF and slide scripts are committed, with pinned measurements,
calculation tables and content receipts. The release was checked with a complete
and cached replay producing identical classifications, independent score
recalculation, retained-target checks and package/API tests.
