# Numerical audit of Neither and Unassigned

Evidence: OncoRef 1.24.0, `oncoref.cta-selection.2026-10-29.4`. Reproduce with `.venv/bin/python scripts/audit_cta_neither_boundaries.py`. This report does not change the released assignments.

## Structural fallback counterfactual

| Category | Current | Fallback only |
| --- | --- | --- |
| Unsafe | 18031 | 18031 |
| Unassigned | 1290 | 527 |
| CTA | 522 | 522 |
| Neither | 20 | 783 |

763 Unassigned proteins pass required normal coverage, have no unresolved evidence and have no normal-screen hold. 757 fail only reproductive association; six additionally lack a qualifying cancer sample. This is a mechanical counterfactual, not a claim that these proteins are clinically safe or all belong in Neither.

## Important audit counts

| Finding | Proteins |
| --- | --- |
| Current Neither with a normal-screen hold | 7 |
| Current Neither with non-permitted normal >30 TPM or nTPM | 6 |
| Normal-pass Unassigned with positive reproductive IHC in every member | 52 |
| Subset with Medium/High, Supported/Enhanced reproductive IHC | 38 |
| Unsafe with RNA-only exclusions <=25% over every Unsafe ceiling | 138 |

The IHC rows are review evidence, not a new admission rule. Boundary proteins may still have Manageable excess or other admission failures; removing one RNA exclusion would not necessarily make them candidates.

## All current Neither proteins

| protein_group_label | neither_profile | reproductive_ntpm | healthy_maximum_ntpm | healthy_maximum_tpm | oncoref_maximum_tpm | normal_screen_hold |
| --- | --- | --- | --- | --- | --- | --- |
| CRYGN | diffuse_manageable | 0.7 | 15 | 18.04 | 334.6 | True |
| KRT34 | diffuse_manageable | 0.8 | 64.5 | 0.9323 | 472.4 | True |
| KRT85 | diffuse_manageable | 1.1 | 167.2 | 0.2186 | 413.2 | True |
| KRT71 | diffuse_manageable | 0.1 | 96.9 | 0.9216 | 569 | True |
| OXGR1 | diffuse_manageable | 0.6 | 5.6 | 3.5 | 162 | False |
| KRT9 | diffuse_manageable | 0.2 | 11.8 | 2.519 | 1866 | False |
| LEP | diffuse_manageable | 4.4 | 296.7 | 194.3 | 907.8 | True |
| LCE2A | diffuse_manageable | 0 | 114.9 | 119.2 | 770.5 | True |
| KRT25 | diffuse_manageable | 0.8 | 108.8 | 0.263 | 1352 | True |
| MTCP1 | diffuse_reproductive | 5.1 | 4.1 | 7.319 | 492.3 | False |
| FAM90A5 | low_expression | 0 | 0.1 | 0 | 0.4293 | False |
| USP17L5 | low_expression | 0 | 0.1 | 0 | 0.165 | False |
| USP17L12 | low_expression | 0 | 0.5 | 0 | 1.649 | False |
| FAM90A3 | low_expression | 0 | 0.1 | 0 | 0.5907 | False |
| USP17L18 | low_expression | 0 | 0.2 | 0.2 | 1.575 | False |
| TAF11L4 | low_expression | 0 | 0 | 0 | 1.777 | False |
| TMDD1 | diffuse_manageable | 0.3 | 13.5 | 1.272 | 24.85 | False |
| FAM90A22 | low_expression | 0 | 0 | 0 | 0.1795 | False |
| FAM90A14 | low_expression | 0 | 0.2 | 0 | 0.4801 | False |
| FAM90A7 | low_expression | 0 | 0 | 0 | 0.5501 | False |

Normal maxima are supported non-permitted normal compartments, not maxima over every donor or every tissue in the body. Source units remain separate. The tumor column is the OncoRef tumor maximum; HPA pTPM is retained separately in the CSV.

## Tumor-positive, low-normal examples

| protein_group_label | healthy_maximum_ntpm | healthy_maximum_tpm | reproductive_ntpm | developmental_all_lineage_max_tpm | oncoref_maximum_tpm | oncoref_cancer_code | hpa_tumor_group_lower_bound_ptpm | optional_DICE_complete | optional_cornea_complete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| USP17L3 | 0.1 | 0 | 0 | 0.0006052 | 3.82 | HNSC | 5.984 | True | True |
| FAM90A9 | 0.1 | 0 | 0 | 0.221 | 1.75 | STAD | 15.1 | False | False |

Neither example has a measured HPA tissue-IHC row in the pinned group-IHC table. No positive IHC is therefore not equivalent to negative IHC. Tumor maxima do not establish prevalence, malignant-cell origin, or protein production.

## Disposable bulk examples

| protein_group_label | tonsil | appendix | gallbladder | healthy_maximum_ntpm | healthy_maximum_tpm |
| --- | --- | --- | --- | --- | --- |
| IL36B | 3.7 | 0 | 0 | 7.9 | 7.163 |
| IL1F10 | 2.7 | 0 | 0 | 12.1 | 13.94 |
| LCE3B | 3 | 0 | 0 | 0 | 0 |
| PLA2G2E | 2.4 | 0 | 0 | 0.7 | 0.2677 |
| LCE3C | 23.4 | 0 | 0 | 1.7 | 0.1582 |

These are examples with permitted bulk RNA >=2 nTPM, not an exhaustive Disposable-only set: hair-follicle and melanocyte compartments require their specific cell data. Skin as a whole is Manageable.

## Interpretation and limitations

- `low_expression` unnecessarily requires low tumor expression and measured developmental profiles; those requirements do not establish adult-normal safety.
- `diffuse_manageable` counts >=2 Manageable tissue labels with RNA >=2. It can allow above-ceiling expression and counts thymus, unlike the separate organic diffuseness score.
- Current reproductive admission is RNA-only (>=5); positive reproductive IHC is not an independent admission route.
- Required-panel completeness is not exhaustive biological coverage. Rare cells, activation states and paralog-specific read assignment remain material limitations.
- Low observed expression is not proof of absence. Optional-source missingness remains explicit.
- Review queues are overlapping diagnostics, not disjoint replacement categories.
