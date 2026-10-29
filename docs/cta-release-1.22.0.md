# OncoRef 1.22.0: safety and CTA identity

The human research selection contains **546 CTA protein-equivalence groups
(582 loci), 18,031 Unsafe, 20 Neither and 1,266 unassigned**, across 19,863 groups.
MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B pass the uniform rules.

[Full report and funnel](audits/cta-source-walkthrough-20261028/index.md) ·
[28-slide PDF](audits/cta-source-walkthrough-20261028/cta-selection-1.22.0.pdf) ·
[One-page logic PDF](images/cta-evidence-axes.pdf) ·
[Rules](cta-evidence-axes.md) · [API](cta-selection.md)

## Four independent requirements

1. **Safety:** absolute normal RNA limits by tissue tier and usable local IHC,
   qualifying normal protein and applicable clinical safety findings. Only this
   assessment can establish Unsafe. Manageable excess blocks CTA under the
   current policy without automatically establishing Unsafe.
2. **Cancer expression:** one measured cancer sample ≥2 TPM or HPA pTPM. There
   is no tumor-versus-normal ratio and RNA does not prove protein abundance.
3. **Reproductive/developmental association:** reproductive ≥5 nTPM or a supported
   developmental lineage mean ≥5 TPM, with the existing cell-support requirements.
4. **Adult restriction:** a qualifying reproductive signal exceeds supported
   adult normal nTPM, or a qualifying developmental signal exceeds supported
   adult normal TPM. Vital, Protected and Manageable observations form the adverse
   comparator. Reproductive and Disposable do not. CP10K remains safety-only.

The last three can block CTA admission but cannot supply Unsafe evidence. Neither
requires a measured low-expression or diffuse-Manageable profile. Other incomplete
or nonqualifying cases remain unassigned. These are research heuristics; source
normalization, sampling and antibody sensitivity limit interpretation.

The 29 October restriction change moves **TGIF2-RAB5IF and ATF7-NPFF** from the
preceding 548-group draft to unassigned. Their developmental means are slightly
below their supported normal TPM values. **No safety exclusion changed.**
[Individual changes](audits/cta-selection-20261028/adult-restriction-changes.csv) ·
[Comparison receipt](audits/cta-selection-20261028/adult-restriction-comparison.json).

## Source coverage and provenance

This release adds the source-version API and expression-coverage repairs described
in the [full report](audits/cta-source-walkthrough-20261028/index.md). It combines
HPA normal panels with direct GTEx v11 and HPA brain coverage, adds Petropoulos
preimplantation measurements alongside Tyser gastrulation, and resolves source
identities while quarantining conflicting aliases. The saved selection pins its
normal HPA measurements to v23; current default RNA loaders can resolve v25.1.
Source versions, annotations, processing and units remain distinct fields.

The public audit records four assessment rows per human group, including passes:

```python
from oncoref import cta_selection

rows = cta_selection.results(species="human", category="CTA")
evidence = cta_selection.audit(species="human").explain("PRAME")
assessments = evidence[evidence.layer.str.startswith("Assessment: ")]
```

Dog's 14 native nominations remain unassigned/unassessed; mouse selection remains
unavailable. Human thresholds or clinical evidence never transfer through orthology.
Neither a missing observation nor an unavailable species becomes zero expression.

## Downstream migration

Pin the package version and selection policy
`oncoref.cta-selection.2026-10-29`. Use `cta_selection.results()` for this selection
and keep its group/member identity and evidence flags. Unassigned has a null
category, not Neither or Unsafe. MAGE-A4, PRAME and CTAG2 should not inherit older
selection dispositions when adopting this policy.

Historical `cta_genes`/624-gene compatibility APIs and the v1.20 exact-protein
near-self background remain separately versioned. This release does not silently
rewrite k-mer backgrounds or product-specific blacklists. A gene-group admission
is not clearance of every isoform, peptide, receptor or unmeasured tissue.
