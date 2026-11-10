# Thresholds, patient identity and named CTA placements

**Corrected draft: 638 CTA, 260 Cancer Antigen only, 22 Reproductive Antigen only, 714 Other and 18,259 Unsafe.** These use the existing fivefold single-tumor rule. [Every protein](../rebuild/report.html) · [Named CTA families and studied targets](named-ctas.md) · [All threshold values as CSV](thresholds.csv).

## Samples and patients

OncoRef's cancer matrices have one column per sample/library. The previous comparison resolved only TCGA, TARGET and MMRF patient barcodes, omitting patient metadata already available in OncoRef. `expression.sample_metadata` now joins actual matrix columns to the central `samples.sample_identities` policy. The rebuild and threshold analysis use that same mapping. Source-specific patient IDs, physical sample keys, status, basis and references stay separate.

The correction resolves **2,471 additional samples**. **13,164 of 14,440 sample columns** have documented identities, representing **13,034 patient IDs**; **1,276 samples** remain unresolved. Identities use each source's documented patient namespace; unavailable cross-study linkages remain unknown. A copied sample label in a GEO case field is not accepted as a patient mapping.

| Source namespace | Additional resolved samples |
|---|---:|
| BEATAML_OHSU_2022 | 600 |
| CGCI_BLGSP | 182 |
| CLLMAP_2022 | 708 |
| OPENPBTA_V23_CBTN_CRANIO | 29 |
| OPENPBTA_V23_DIPG_H3K27 | 32 |
| TREEHOUSE_POLYA_25_01 | 905 |
| TREEHOUSE_RIBOD_25_01 | 15 |

The [CLL-map source](https://www.cllmap.org/downloads.html) specifies six duplicate exclusions plus one suspected MCL exclusion, already applied to OncoRef's 708-sample cohort. Treehouse uses curated donor IDs where present and the existing builder's donor-family rule for older identifiers. Multiple S01/S02 libraries from one family count once. [Complete source coverage](patient-coverage.csv).

**Multiple samples are not mandatory.** One qualified tumor ≥5 × max(N,1) still qualifies. Patient identity is needed only to establish independent support for the repeated route and patient-based percentiles. Unknown IDs do not erase measured expression.

| Protein | Before correction | Corrected category | Patients above normal |
|---|---|---|---:|
| DEFB106A/B | Reproductive Antigen only | CTA | 2 |
| TKTL2 | Reproductive Antigen only | CTA | 34 |
| ARL13A | Reproductive Antigen only | CTA | 14 |
| OR4X1 | Other | Cancer Antigen only | 2 |
| OR56B2 | Other | Cancer Antigen only | 10 |
| OR4C15 | Other | Cancer Antigen only | 2 |
| OR4A15 | Other | Cancer Antigen only | 3 |
| KRTAP19-2 | Other | Cancer Antigen only | 4 |
| MINDY4B | Reproductive Antigen only | CTA | 2 |
| LACTBL1 | Other | Cancer Antigen only | 3 |
| FAM90A26 | Other | Cancer Antigen only | 2 |
| OR4F3 | Other | Cancer Antigen only | 3 |
| OR4C13 | Other | Cancer Antigen only | 4 |
| DPEP2NB | Reproductive Antigen only | CTA | 3 |
| PRAMEF13 | Reproductive Antigen only | CTA | 2 |
| FAM246B | Other | Cancer Antigen only | 2 |

The expression values and normal exclusions did not change. [Current-rule changes](patient-mapping-category-changes.csv) and [tenfold-alternative changes](patient-mapping-proposed-changes.csv) retain complete evidence. The baseline is [commit 2e24b126](https://github.com/pirl-unc/oncoref/commit/2e24b126); the [threshold comparison](../thresholds/index.md) has been regenerated with the corrected identities.

## Threshold inventory

These are the rules used in the current draft or inherited normal screen. Proposed alternatives are labeled separately. Historical 1.36.0 used a ≥5 reproductive RNA threshold, tumor presence >0 and the >10% concentration veto; those are not the current draft's association rules. Older compatibility taxonomy parameters do not replace the active `cta_selection_policy` ceilings.

| Rule | Threshold | Status | Meaning |
|---|---|---|---|
| Cancer: repeated | ≥2 documented patients; each tumor ≥1 TPM and strictly above N | Current draft | N is the largest supported nonreproductive, non-Disposable normal context statistic; repeated libraries count once. |
| Cancer: single tumor | ≥5 × max(N, 1) TPM in ≥1 sample | Current draft | One strong sample suffices, including one without a resolved patient ID. |
| Reproductive bulk RNA | ≥1 native nTPM | Current draft | Includes reproductive, placental and developmental scope; normal exclusions still apply. |
| Reproductive-cell RNA | ≥1 nTPM and ≥20 reads | Current draft | Indirect sites also require complete outside-cluster coverage, every outside value <1, and no curated normal-induction conflict. |
| Developmental RNA | Mean ≥1 TPM; ≥20 cells, ≥5 positive and ≥10% positive; positive means ≥1 TPM | Current draft | Qualifying developmental lineage/stage; bulk fetal/reproductive scope is retained separately. |
| Reproductive protein | Reliable mapped positive IHC or reviewed human reproductive protein for every member | Current draft | An alternative to RNA association; antibody/identity conflicts remain unresolved. |
| Normal RNA ceilings | Vital 4/8/16/32; Protected 8/16/32/64; Manageable 16/32/64/128 | Retained normal screen | Columns are missing/uncertain, Approved, Supported, Enhanced matched negative IHC. Units remain native TPM or nTPM. |
| Reviewed clinical RNA credit | One IHC rank increment, capped at rank 3 | Retained normal screen | Only qualifying product-specific evidence; no stacking or override of injury/protein exclusions. |
| Normal RNA consequence | Above Vital/Protected ceiling → Unsafe; above Manageable ceiling → hold; Manageable >128 → Unsafe | Retained normal screen | Equality to the ceiling passes. Disposable and reproductive scope have no current numerical ceiling in this gate. |
| HPA tissue-cell support | ≥20 reads | Retained normal screen | Cell-specific tier; only exact/declared compartment-matched IHC can supply credit. |
| HPA/Monaco immune | Cell-specific 4/8/16 bases on native nTPM/TPM | Retained normal screen | No matched negative-IHC credit; a reproductive location does not exempt shared immune cells. |
| DICE immune | Third-highest reported TPM; ≥3 observations | Retained normal screen | Apply cell-specific ceiling; unpaired distributions keep their bounds. |
| Adult cornea | ≥5 positive cells and ≥20 counts in ≥2 sample groups; RNA ceilings divided by 100 | Retained normal screen | Native CP10K screening scale, not a TPM conversion. |
| Rare-cell detection | ≥1% detected, ≥3,000 cells and ≥10 expressing cells | Retained normal screen | Existing brain/spinal cord/heart/lung panel: Vital → Unsafe; Protected → review hold. |
| Normal protein | Reliable non-Low positive in Vital/Protected → Unsafe; Manageable → hold | Retained normal screen | Low is flagged; uncertain is unresolved. Reviewed/external findings retain reagent, cell and identity scope. |
| Required RNA coverage | Every nonreproductive combined-bulk tissue concept for every member | Retained normal screen | The comparison uses 36 concepts. Missing optional cell/immune/cornea panels alone does not block; observed exclusions still apply. |
| Translation entry | Reference coding CDS, or supported ncORF: MS peptide or ≥2 independent Ribo-seq studies | Current draft | All proteins still undergo the normal screen. |
| Representative translation | Complete primary-locus coding sequence; no LoF/NMD, incomplete tags, internal stops or non-AA20 residues | Current draft | MANE Select first, Ensembl canonical second; conflicting evidence requires review. X is not allowed. |
| Nested coding sequence | <50 amino acids inside a same-strand coding gene without independent ncORF evidence | Current draft | Other:Sequence. |
| Unsafe sequence overlap | ≥95% of 9-mers occur in Unsafe proteins | Current draft | Other:Sequence; this does not create expression-based Unsafe evidence. |
| Reproductive concentration | Recorded normal/(normal + reproductive maximum); previous >10% veto removed | Descriptive only | Developmental and adult reproductive association now share the same admission flow. |
| Published nomination | Exact-member claim in a curated CA/CTA/oncofetal target source | Annotation | Retain source, table and claim; neither automatic admission nor exemption from normal findings. |
| Proposed single-tumor rule | ≥10 × max(N, 1) TPM | Counterfactual; not adopted | Repeated route remains ≥2 patients at ≥1 TPM and above N. |
| Proposed target maximum | max(cancer, reproductive/developmental) ≥1 and strictly above every measured nonreproductive value | Counterfactual; not adopted | Source-compatible values only; missing compartments are not zeros. |
| Proposed Disposable caps | ≤32 TPM and ≤32 nTPM, assessed separately | Counterfactual; not adopted | No unit conversion or IHC-derived numerical allowance for these proposed caps. |
| Percentile alternatives | Within-tumor p80/p90/p95; across-patient p80/p90/p95 with ≥20 measured patients | Sensitivity analysis | Within-tumor ranks use the full measured catalog; across-patient values use the lower observed order statistic. |

The base normal ceilings are **4 for Vital, 8 for Protected and 16 for Manageable**, before qualifying IHC or clinical adjustments. They are screening parameters, not measured injury thresholds. [Normal-screen implementation](https://github.com/pirl-unc/oncoref/blob/feat/antigen-set-policy/oncoref/cta_selection_policy.py), [antigen membership implementation](https://github.com/pirl-unc/oncoref/blob/feat/antigen-set-policy/oncoref/target_selection.py) and the [tissue registry with cited rationales](../../../tissue-hazards.md) identify the distinct policy layers.

## Hair follicles

**Hair follicles are already classified Disposable; whole skin is Protected.** The missing component in the numerical cap comparison is isolated hair-follicle RNA, not a missing ontology bucket. HPA v23's tissue-cell RNA file has no hair-follicle cluster. It does have cutaneous melanocytes on its separate nTPM scale; the [cell audit](../expression-maxima/index.md) retains read support and missing loci.

HPA v23 also contains hair protein staining: cortex/medulla, cuticle, internal/external root sheath, and hair-follicle cells in skin sections. Thus “no hair data” was too broad. A separate mapping defect was found: the skin cell-label regex `hair follic` does not match `hair follicles` and falls back to Protected skin. [Issue #687](https://github.com/pirl-unc/oncoref/issues/687) tracks correction together with a replay of the inherited normal-protein gates and affected sequence containment; the patient-only correction here does not revise those gates.

Human hair-follicle RNA datasets are available, including [GSE129611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129611) and [the human scalp single-cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9126928/). They have not been added to the pinned cap panel. Missing hair RNA cannot count as zero Disposable expression.
