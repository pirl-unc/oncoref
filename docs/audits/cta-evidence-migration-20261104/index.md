# Historical evidence and clinical attribution

All **196 historical findings** have an explicit disposition and linked
current record: 178 revised, 14 retained, 3 unresolved, 1 rejected.
The [complete ledger](historical-dispositions.csv) retains each original finding,
prior status, access depth, limitations and antibody adjudication. This inventory
does not complete the wider literature reviews in issues #636 and #638.

## Material unresolved cases

MAGEC2, LGALS13 and PLCZ1 enter **Other with a sequence-review requirement**. Their
reviews are accepted as complete at the stated access depth; the interpretations
remain unresolved. No new positive normal-protein exclusion is asserted.

| Case | Evidence and unresolved question |
|---|---|
| MAGEC2 | [Li 2003](https://pubmed.ncbi.nlm.nih.gov/12920247/) reports Purkinje staining with TC-1. The abstract was rechecked; full reagent/control validation was unavailable. A negative bulk profile cannot resolve specificity. |
| LGALS13 | [1999 PP13 work](https://pubmed.ncbi.nlm.nih.gov/10527825/) reports extra-placental immunoreactivity, while [2009 work](https://pubmed.ncbi.nlm.nih.gov/19497882/) supports placental restriction. Adult/fetal anatomy and galectin identity remain unresolved; only abstracts were rechecked. |
| PLCZ1 | [Zhu 2003](https://www.cuhk.edu.hk/ecbrc/cliparts/Publication%20list/Lancet2003.pdf), Methods/Figure 1: qualitative pancreatic NYD-SP27 RNA remains unresolved quantitatively and by isoform. Immunostaining measured CFTR, so the rejected PLCZ1 protein interpretation stays rejected. |

The other 178 historical antibody conflicts retain their observations and prior
adjudications. Their blanket historical hold is **revised by the declared current
HPA level/reliability, prior reagent-eligibility and compartment rules**, not by a
claim that every antibody conflict was experimentally resolved. Current assay
decisions remain reproducible and visible in the candidate audit; consumers must
not apply a second, private historical IHC classifier.

## Product injury and competing peptide explanations

The [original neurological report](https://pubmed.ncbi.nlm.nih.gov/23377668/)
recorded three neurological cases, including two deaths, among nine recipients
of the HLA-A*02 MAGE-A3/A9/A12-reactive TCR. The
[2021 re-examination](https://pmc.ncbi.nlm.nih.gov/articles/PMC7946352/) challenges
the MAGEA12 attribution and supplies EPS8L2 cross-reactivity evidence. Its RNA FISH,
peptide assays, HLA immunopeptidomics and cell-line killing experiments support an
alternative without proving the initiating patient-brain peptide. Rare MAGEA12
cells and other targets remain possible. Full-text Methods, Results, Figures 1–4,
Tables 1–3 and Discussion were reviewed; supplementary experiments were not reanalysed.

Fatal product-specific injury still blocks favorable clinical credit. Unique
gene-level causality is unresolved. **MAGEA12 remains Unsafe through independent
rare Vital-cell RNA evidence**, whose limitations remain separate. The HLA-A*01
a3a/titin cardiac cases involve a different receptor and peptide.

## Reproducible consequences

Relative to the decision-contract PR: **216
classification changes, 197 sequence-role
changes and 704 reason changes**.
[All changes](changes-from-contract-pr.csv) include unchanged categories with
changed reasons. The [series-wide delta](../cta-definition-20261104/changes-from-1.29.0.csv)
uses the frozen 1.29.0 baseline.

```python
from oncoref import gene_knowledge, cta_selection_audit, tcr_trials
gene_knowledge.historical_dispositions("PLCZ1")
cta_selection_audit.audit().findings("PLCZ1")  # includes rejected/context findings
tcr_trials.events()  # causal_findings retains both peptide hypotheses
```

[Current findings](current-findings.csv.gz), [clinical ledger](clinical-counterevidence.json),
and [receipt](manifest.json) preserve the source and policy boundaries.
