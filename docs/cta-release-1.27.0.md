# OncoRef 1.27.0: simpler decisions and adult Disposable evidence

**701 CTA / 18,062 Unsafe / 453 Other / 647 Unassessed** representative proteins.
All 701 previous CTAs and all 18,062 previous Unsafe assignments are preserved.

[Simple rules, funnel and four target reviews](audits/cta-simplification-20261031/index.md) ·
[Full source report](audits/cta-source-walkthrough-20261028/index.md) ·
[Rules and API](cta-selection.md) ·
[PDF logic slides](audits/cta-source-walkthrough-20261028/cta-selection-1.27.0.pdf)

## One decision, separate evidence fields

The main decision is CTA, Unsafe, Other or Unassessed. Other now records only
its unmet identity requirement: **448 Nonreproductive / 5 Diffuse**; no current
protein fails only cancer presence. The former Low, Weak, Normal, Tumor-supported
and Disposable labels mixed abundance, confidence and anatomy.

`cta_selection.profiles()` exposes those evidence dimensions independently:
separate native TPM/nTPM/pTPM/CP10K summaries, replicated tumor RNA, normal protein,
reproductive presence, developmental measurements and adult Disposable expression.
The top-level result and identity fields remain available through `results()`.
Update old subtype filters explicitly; the new subtype names describe reasons.

## Disposable is an adult tissue disposition

An adult compartment's permission never depends on embryonic/fetal/placental
expression. The adult Disposable-only expression attribute now recognizes
TCHHL1, LCE3B and LCE3C among Other proteins; developmental measurements remain
separate, including missing data. The existing adult anatomy mapping and safety
ceilings are preserved. Whole skin remains Manageable; mapped follicles and
cutaneous melanocytes have conditional Disposable allowances.

## One additional normal-evidence hold

OR10A4 moves from Other to Unassessed. Original human olfactory-mucosa data show
RNA across 26 adult donors, and a separate study detects normal nasal-fluid RNA.
These native qPCR/copy measurements are not TPM or protein. They prompt review of
a poorly covered normal compartment, rather than a claimed clinical toxicity.
The two source spreadsheets are pinned through the shared `catalog` / `reference_data`
API at version `2014-05-06`. [Measurements, citations and limitations](audits/cta-simplification-20261031/index.md#or10a4-has-additional-normal-human-nasal-evidence).

DBX1, OR2B2 and KRTAP20-2 retain strong tumor evidence but lack qualifying
reproductive identity. Their atlas normal peaks and early-developmental values
are in the report. High tumor expression alone does not make a protein a CTA.

## Downstream adoption

```python
from oncoref import cta_selection, catalog

results = cta_selection.results()
profiles = cta_selection.profiles()
roles = cta_selection.sequence_dispositions()
normal_source = catalog.ensure("olfactory_verbeurgt_cq", version="2014-05-06")
```

Retain policy/reference identity, raw evidence and sequence QC. OR10A4 now requires
review for shared-sequence allowance. Historical FASTA/k-mer bundles remain
separately versioned; see the unresolved sequence-bundle and translation-QC
contracts in [#641](https://github.com/pirl-unc/oncoref/issues/641) and
[#657](https://github.com/pirl-unc/oncoref/issues/657). Human policy does not transfer
to dog/mouse; their profile measurements stay missing.
