# Which sources have missing evidence?

These counts describe the **1,325 unassigned groups with incomplete evidence**, not every gene in the reference. Sources overlap: their counts must not be added.

| Evidence gap | Groups |
| --- | --- |
| normal panel missing | 574 |
| tumor measurement missing | 37 |
| tumor group profile incomplete | 30 |
| embryonic measurement missing | 1006 |
| unresolved evidence | 33 |

The earlier **71 tumor-measurement gaps** combined 41 groups without a mapped measurement and 30 measured groups with an incomplete low-expression profile. After the identity repair, **37** groups lack a mapped measurement, while **30** have measured tumor values but incomplete group coverage for the low-expression profile. Multi-locus group lower bounds cannot establish low expression. This restriction does not erase observed values and does not block tumor admission when a sample exceeds the floor.

The repair also fixes a source-name sentinel bug: the text `NaN` was interpreted as a historical SCN11A synonym instead of a missing name, contradicting valid ENSG IDs. Seven affected cohort views were reprocessed from their pinned sample matrices and original passing-sample QC: ACINIC, ADCC, CRANIO, DIPG, MENINGIOMA, SARC_MMNST and VSCC. Their corrected measurements replace the older cohort summaries; source samples and units are unchanged. [Repaired cohort evidence and hashes](../cta-selection-20261028/tumor-identity-repair/manifest.json).

## Required normal panels

| Source | Groups | Examples |
| --- | --- | --- |
| HPA v23 consensus bulk | 574 | RTEL1-TNFRSF6B; GSTT2; HMGB1P1; OR11H1; ZNF670-ZNF695 |
| HPA v23 tissue-cell clusters | 352 | RTEL1-TNFRSF6B; ZNF670-ZNF695; LINC03040; ZNF840P; MEIS3P2 |
| HPA v23 immune cells | 312 | RTEL1-TNFRSF6B; ZNF670-ZNF695; LINC03040; ZNF840P; MEIS3P2 |
| Monaco immune cells via HPA v23 | 312 | RTEL1-TNFRSF6B; ZNF670-ZNF695; LINC03040; ZNF840P; MEIS3P2 |

Completeness requires all 50 HPA consensus bulk tissues, 426 HPA tissue/cell contexts, 19 HPA immune populations and 30 Monaco immune populations for every protein-equivalence-group member. A zero-valued row is measured coverage. A missing tissue row, unresolved identity or absent group member is not. DICE and adult cornea remain optional and do not contribute to this mandatory-panel gap.

[Every missing group/source/context](normal-panel-coverage-gaps.csv.gz) · [Every member, current/historical IDs, synonyms and observed/missing contexts](normal-panel-member-coverage.csv.gz) · [Source-row identity matches and unresolved alternatives](coverage-identity-matches.csv.gz).

## Cancer and embryonic sources

| Source | Cohort views | Groups without any mapped measurement |
| --- | --- | --- |
| HPA 25.1 | 31 | 406 |
| OncoRef | 132 | 37 |

[Per-group, per-cohort coverage](tumor-cohort-coverage.csv.gz) identifies each OncoRef cohort view and HPA cancer/cohort, the measured-member count and missing member IDs. Cohort views can share biological samples. [Tumor identity alternatives](tumor-missing-identity-candidates.csv.gz) retains source rows matching the missing groups through any registered current/historical ID or known name candidate, including conflicts that cannot be substituted.

The embryonic input is specifically **Tyser et al. (2021), E-MTAB-9388**, EBI Single Cell Expression Atlas's filtered gene-TPM matrix using Ensembl 106. The **1,006** gaps mean at least one group member has no usable mapped row in this input. No later fetal atlas contributes to this coverage flag. [Missing embryonic member IDs](embryonic-coverage-gaps.csv.gz). The broader author raw-count table is a separate diagnostic and is not a TPM substitute.

The new identity pass recovers 19 gene rows omitted by the earlier exact-ID filter. It keeps 5 conflicting alternate-row cases unresolved. [Every attempted recovery and source ID](developmental-identity-recoveries.csv). Canonical rows take precedence; multiple alternative rows are accepted only when their entire cell-expression vectors agree. Different values are neither added nor arbitrarily selected.

## Unresolved evidence

All **33** current unresolved groups have contributing **HPA v23 normal-tissue IHC rated Uncertain**, with Medium or High staining. These are measured but unresolved observations. [Exact gene, tissue, cell type, staining and reliability](unresolved-protein-contexts.csv.gz).

## Identity matching

The shared matcher checks canonical ENSG IDs, registered historical and alternate IDs, and unique official/synonym names. It preserves the source identifier and matching method. A conflicting explicit ENSG ID is not overwritten with a symbol guess, and ambiguous aliases do not get assigned to the first matching paralog. An empty history list means no replacement is registered in this snapshot, not that a gene has never had another identifier.

For example, **POLR2A / ENSG00000181222** and **CIMIP3 / ENSG00000287363** each have 41 of 50 HPA bulk tissue rows. CIMIP3 already matches the source name **GUCA1ANB** through its stable ID. Both lack amygdala, basal ganglia, cerebellum, hippocampal formation, hypothalamus, midbrain, pituitary gland, spinal cord and vagina rows. Renaming the gene cannot restore those omitted measurements. Conversely, a tumor row named NOX5 but carrying the SPESP1-NOX5 readthrough ID remains an explicit ID/name conflict; it is not assigned to either locus by name alone.

This audit distinguishes absent accepted rows, partial tissue coverage and unresolved identity alternatives. It cannot recover a measurement that the downloaded source omitted. The registry also contains an explicit unresolved protein-group member, ENSG00000275869 in the TP53TG3E/F group; its absence from the canonical catalog remains visible rather than being silently deleted from the group.

[Source hashes and audit code receipt](coverage-audit-receipt.json).
