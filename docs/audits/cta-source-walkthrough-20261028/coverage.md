# Which sources have missing evidence?

These counts describe the **475 unassigned proteins with incomplete evidence**, not every gene in the reference. Sources overlap: their counts must not be added.

| Evidence gap | Proteins |
| --- | --- |
| normal panel missing | 453 |
| tumor measurement missing | 34 |
| tumor group profile incomplete | 3 |
| embryonic measurement missing | 86 |
| unresolved evidence | 28 |

The earlier **71 tumor-measurement gaps** combined 41 groups without a mapped measurement and 30 measured groups with an incomplete low-expression profile. After the identity repair, **34** groups lack a mapped measurement, while **3** have measured tumor values but incomplete group coverage for the low-expression profile. Multi-locus group lower bounds cannot establish low expression. This restriction does not erase observed values and does not block tumor admission when a sample exceeds the floor.

The repair also fixes a source-name sentinel bug: the text `NaN` was interpreted as a historical SCN11A synonym instead of a missing name, contradicting valid ENSG IDs. Seven affected cohort views were reprocessed from their pinned sample matrices and original passing-sample QC: ACINIC, ADCC, CRANIO, DIPG, MENINGIOMA, SARC_MMNST and VSCC. Their corrected measurements replace the older cohort summaries; source samples and units are unchanged. [Repaired cohort evidence and hashes](../cta-selection-20261028/tumor-identity-repair/manifest.json).

## Required normal panels

| Source | Proteins | Examples |
| --- | --- | --- |
| HPA v23 consensus bulk | 374 | RTEL1-TNFRSF6B; ENSG00000100101; ENSG00000125695; ZNF670-ZNF695; ENSG00000141979 |
| Combined normal bulk: HPA consensus + HPA brain + direct GTEx | 366 | RTEL1-TNFRSF6B; ENSG00000100101; ENSG00000125695; ZNF670-ZNF695; ENSG00000141979 |
| HPA v23 tissue-cell clusters | 335 | RTEL1-TNFRSF6B; ENSG00000100101; ENSG00000125695; ZNF670-ZNF695; ENSG00000141979 |
| HPA v23 immune cells | 304 | RTEL1-TNFRSF6B; ENSG00000100101; ENSG00000125695; ZNF670-ZNF695; ENSG00000141979 |
| Monaco immune cells via HPA v23 | 304 | RTEL1-TNFRSF6B; ENSG00000100101; ENSG00000125695; ZNF670-ZNF695; ENSG00000141979 |
| Required GTEx v11 bulk panel (52 contexts) | 413 | BORCS8-MEF2B; ENSG00000100101; CORO7-PAM16; ENSG00000125695; ZNF670-ZNF695 |

Bulk completeness now requires every nonreproductive tissue concept in the HPA consensus panel, using compatible measurements from HPA consensus, the separate HPA brain atlas, or direct GTEx v11. The fixed GTEx v11 core additionally requires all 52 native bulk contexts, including its nine adult Reproductive contexts, for every member locus. The separate requirements remain 426 HPA tissue/cell contexts, 19 HPA immune populations and 30 Monaco immune populations for every group member. The table distinguishes single-export gaps from gaps remaining after the source union. A zero-valued row is measured coverage. A missing tissue row, unresolved identity or absent group member is not. DICE and adult cornea remain optional and do not contribute to this mandatory-panel gap.

[Every missing group/source/context](normal-panel-coverage-gaps.csv.gz) · [Every member, current/historical IDs, synonyms and observed/missing contexts](normal-panel-member-coverage.csv.gz) · [Source-row identity matches and unresolved alternatives](coverage-identity-matches.csv.gz).

## Cancer and embryonic sources

| Source | Cohort views | Groups without any mapped measurement |
| --- | --- | --- |
| HPA 25.1 | 31 | 332 |
| OncoRef | 132 | 34 |

[Per-protein, per-cohort coverage](tumor-cohort-coverage.csv.gz) identifies each OncoRef cohort view and HPA cancer/cohort, the measured-member count and missing member IDs. Cohort views can share biological samples. [Tumor identity alternatives](tumor-missing-identity-candidates.csv.gz) retains source rows matching the missing groups through any registered current/historical ID or known name candidate, including conflicts that cannot be substituted.

Embryonic coverage now combines [Tyser et al. (2021), E-MTAB-9388](https://doi.org/10.1038/s41586-021-04158-y) (gastrulation) and [Petropoulos et al. (2016), E-MTAB-3929](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results) (preimplantation days 3–7). Both supply filtered gene-level TPM, with stages and source-specific lineages kept separate. Petropoulos contributes 1,519 measured cells from 88 embryos; the metadata-only ten additional cells are not inserted into the matrix. It recovers **683 of the previous 1,006 Tyser-missing groups**. The current **86** combined gaps are after reclassification and therefore have a changed denominator. [Coverage and missing member IDs for each atlas and their union](embryonic-coverage-gaps.csv.gz), [Petropoulos source](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-3929/Results), [download hashes and scope](../cta-selection-20261028/alternative-evidence/manifest.json). No later fetal atlas contributes yet. The author raw-count table remains a separate diagnostic, not a TPM substitute.

The new identity pass recovers 19 gene rows omitted by the earlier exact-ID filter. It keeps 5 conflicting alternate-row cases unresolved. [Every attempted recovery and source ID](developmental-identity-recoveries.csv). Canonical rows take precedence; multiple alternative rows are accepted only when their entire cell-expression vectors agree. Different values are neither added nor arbitrarily selected.

## Unresolved evidence

The **28** current unresolved groups have contributing **HPA v23 normal-tissue IHC rated Uncertain**, with Medium or High staining. These are measured but unresolved observations. [Exact gene, tissue, cell type, staining and reliability](unresolved-protein-contexts.csv.gz).

## Identity matching

The shared matcher checks canonical ENSG IDs, registered historical and alternate IDs, and unique official/synonym names. It preserves the source identifier and matching method. A conflicting explicit ENSG ID is not overwritten with a symbol guess, and ambiguous aliases do not get assigned to the first matching paralog. An empty history list means no replacement is registered in this snapshot, not that a gene has never had another identifier.

For example, **POLR2A / ENSG00000181222** and **CIMIP3 / ENSG00000287363** each have 41 of 50 HPA bulk tissue rows. CIMIP3 already matches the source name **GUCA1ANB** through its stable ID. Both lack amygdala, basal ganglia, cerebellum, hippocampal formation, hypothalamus, midbrain, pituitary gland, spinal cord and vagina rows. Eight missing tissues are Vital (the seven brain/spinal regions plus pituitary); vagina is reproductive. The separate HPA brain file supplies all seven missing brain/spinal regions for both genes. Direct GTEx v11 also supplies POLR2A pituitary (**109.319 TPM**) and the other missing regions, completing its normal bulk coverage. GTEx v11 also contains CIMIP3 as **ENSG00000287363.4 / GUCA1ANB**, with a measured pituitary median of **0 TPM**, completing its normal bulk coverage. The zero is an observation; no GUCA1A measurement is borrowed. [Combined anatomical coverage](combined-rna-coverage.csv.gz). Exact whole-organ GTEx measurements receive the same reviewed bulk IHC grade as their HPA anatomical equivalent; finer brain subregions do not inherit parent-organ negative staining. HPA nTPM and direct GTEx median TPM remain separate measurements with the same nominal cutoffs; they are not pooled or converted. Conversely, a tumor row named NOX5 but carrying the SPESP1-NOX5 readthrough ID remains an explicit ID/name conflict; it is not assigned to either locus by name alone.

This audit distinguishes absent accepted rows, partial tissue coverage and unresolved identity alternatives. An omitted source row can be supplemented by another compatible source; no missing source row is filled with zero. The registry also contains an explicit unresolved protein-group member, ENSG00000275869 in the TP53TG3E/F group; its absence from the canonical catalog remains visible rather than being silently deleted from the group.

[Source hashes and audit code receipt](coverage-audit-receipt.json).
