# Linked tissue, cell-cluster and immune RNA thresholds

This research comparison uses **one tissue/IHC ladder** for bulk RNA, HPA cell
clusters and sorted immune RNA. Source-specific numerical cutoffs are not
independently tuned. The original reproductive/developmental association,
normal-protein exclusion, RNA coverage and measurement-support requirements
remain. This does not change the released expression classes or forbidden
sequence reference.

## Results

| Tumor admission | Shared normal-RNA ladder, including corneal count scaling | Sensitivity retaining the old 1 CP10K corneal cutoff |
|---|---:|---:|
| ≥1 pTPM in ≥3 samples and ≥5% of one cohort | 244 | 248 |
| **≥10 pTPM in any one measured sample** | **337** | 341 |
| ≥30 pTPM in any one measured sample | 197 | 201 |

The ≥10 result gains 99 and loses six compared with the recurrence result;
it gains 140 compared with the ≥30 result. Relative to the old 201-gene screen,
the linked recurrence result gains 46 and loses three. These comparisons keep
their tumor rules explicit rather than attributing every difference to normal
RNA thresholds. All 20,097 genes were reprocessed; 2,960 qualified for detailed
local assessment, producing 793,280 cell-context records.

**SMC1B and PAGE1 pass all three linked tumor scenarios.** SMC1B's largest DICE
third-highest value is 3.4997 TPM, below the shared Protected ceiling of 6 TPM.
Its stimulated-CD4 population has 58/88 reported values ≥2 TPM, but none ≥5;
the measured finding is retained even though the arbitrary 2 TPM veto is gone.

The eight requested identities MAGEA4, CTAG1A, CTAG1B, XAGE1A, XAGE1B, SSX2,
CT83 and TBXT pass. PRAME and CTAG2 still fail independent Vital tissue RNA
and bulk gates. There are no named-gene overrides.

The corneal sensitivity differs by DKKL1 (0.041542 CP10K), INSL4 (0.112708),
DCAF4L1 (0.021790) and RBAK-RBAKDN (0.078576). These are supported recurrent
statistics, not individual-cell maxima.

The 337 carry **141 historical CTA, 9 Unsafe and 187 Neither** annotations.
Those are existing expression-policy labels, not new classifications. The nine
historical Unsafe flags are retained for TNP1, OBP2A, SMCP, XAGE3, HMGB4, FTMT,
TMEM95, LELP1 and BOD1L2. Their presence must not be read as removal of a previous
protein/antibody concern or as clinical clearance.

[Selected 337-gene table with retained caveats](selected-genes.csv).

## Where the conflicting numbers came from

The 21 October registry defined inclusive candidate ceilings of Vital
**2/4/6/12**, Protected **6/8/12/24**, and Manageable **10/20/30/60**, ordered by
Missing-or-Uncertain, Approved, Supported and Enhanced negative IHC.
Those are user-directed project parameters, not clinical ground truth.

The 26 October illustrative cell screen independently introduced lower/upper
pairs of **2/5**, **5/10**, and **10/20**, and classified all sorted immune assays
as Vital. That is why DICE acquired a **2 TPM** candidate exclusion. Neither
DICE data nor a validated single-cell purity model established this reduction.
The old 201-gene comparison inherited that separate grid.

The protected **6 and 12** in the earlier registry are its Missing and Supported
IHC ceilings, respectively. They were not originally two validated candidate
and toxicity boundaries. In this revision, stronger RNA reference review starts
at **2× the candidate ceiling**, uniformly across tiers and assays. Thus a
Protected context without usable IHC now has a candidate ceiling of **6** and
a stronger-review boundary of **12**. That factor of two is an explicit
provisional policy choice, not evidence of protein production or toxicity.

## One shared ladder

| Tissue or cell tier | Missing/unusable IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 2 | 4 | 6 | 12 |
| Protected | 6 | 8 | 12 | 24 |
| Manageable | 10 | 20 | 30 | 60 |

Values at the candidate ceiling pass. Values above it exclude from candidate
selection; values at or above twice the ceiling additionally enter stronger
RNA reference review. RNA review does not automatically add a sequence to the
forbidden reference. Reliable positive protein evidence remains separate.

Bulk and cell code both obtain candidate values from
`normal_tissue_screening_policy("ihc_adaptive_cta")["negative_ihc_rna_limits"]`.
`cell_rna_policy.linked_thresholds()` generates every cell-assay row from that
definition. The old independent bulk cap of 10, cluster cap of 30, immune cap
of 5, cell IHC multipliers and corneal cap of 1 are not additional gates in the
linked result. Their presence as historical columns does not affect selection.

HPA bulk and cell data retain **nTPM**, DICE and Monaco retain **TPM**. Equal
numerical settings are a provisional policy linkage, not a claim of calibrated
assay equivalence. A source-specific deviation needs ground-truth calibration
and a versioned justification; none are configured in this policy.

Negative IHC only changes a ceiling when there is a reliable reproductive-positive
profile and no recorded antibody concern. Bulk relaxation requires all expected
IHC contexts in the same anatomical tissue to be negative. A cell population
requires the exact tissue and cell label. Missing, partial, unmatched or
Uncertain IHC receives the baseline. Detected nonpermitted IHC, including
Uncertain staining, still excludes selection.

## Immune cells and reproductive context

Identified immune cells default to **Protected**, including DICE, HPA/Monaco
sorted populations and mapped immune populations in HPA tissue clusters.
Their nonimmune neighbors retain the organ-specific tier. No immune subset has
a reviewed Manageable exception in this run; an anatomical location or presumed
replenishment alone does not create one. This is a selection-policy overlay;
the anatomical registry and measured values are preserved.

Explicit reproductive/placental/embryonic/fetal immune contexts stay neutral
for exclusion and do not supply positive association. All 426 source HPA
clusters still contribute to coverage: 399 screened and 27 newly neutral.
The policy does not exempt normal childhood expression.

DICE still uses the **third-highest reported value** per population, so the
candidate ceiling is exceeded only with at least three reported values above
6 TPM. HPA/Monaco use their source summaries. These are not invented independent
donors or individual-donor maximum gates. No equivalent IHC match is assumed
for sorted or activated immune states.

## Corneal counts and the meaning of pTPM

**pTPM means protein-coding transcripts per million**, an RNA normalization.
HPA removes noncoding transcripts and rescales the remaining abundance; nTPM
adds its between-sample normalization. [HPA methods](https://v23.proteinatlas.org/about/assays%2Bannotation).

Our adult Collin corneal statistic is the mean of per-cell counts divided by
whole-cell library size and multiplied by 10,000, including zeros. Therefore
**CPM = 100 × CP10K** exactly for that same statistic. With this 10x 3′ UMI
assay, a bulk-fragment gene-length correction is not the way to convert it.
[Primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343164/).

The result is still not HPA-equivalent nTPM: protein-coding denominators, TMM,
gene inclusion and pooled-cluster versus mean-per-cell aggregation differ.
The linked policy uses the shared nominal per-million ceiling divided by 100
in native CP10K, without relabeling observations. For Vital corneal populations
without IHC, **2 per million corresponds to 0.02 CP10K**, with stronger review
at 0.04. The original **1 CP10K = 100 CPM** rule is exported as a separate
sensitivity so this substantial numerical assumption is visible.

Corneal support remains ≥5 positive cells and ≥20 counts in each of ≥2 distinct
sample groups. HPA tissue clusters still need ≥20 gene reads. Coverage is
reported separately; low support and missing data are not interchangeable with
a measured zero.

## Tumor admission and reproduction

The main alternative asks for **at least one measured tumor sample ≥10 pTPM**
in HPA 25.1, without a minimum prevalence or three-sample requirement. The same
normal-expression rules are compared against the old recurrence rule and the
one-sample ≥30 pTPM rule. These use 21 TCGA and ten validation cohorts, not all
OncoRef cancer datasets. They are tumor nomination thresholds, not an inversion
of a validated normal-tissue safety limit.

All 20,097 reference genes are reprocessed for context neutrality. The complete
union passing association, protein, coverage and either tumor rule receives
local bulk/cell assessment. Old RNA exclusions do not restrict that union.

```bash
python scripts/cta_quantitative_threshold_report.py \
  --neutral-reproductive-immune --linked-tiers --tumor-floor 10 \
  --output /path/to/empty-directory
```

[Every gene and selection scenario](gene-comparison.csv.gz),
[cell observations and decisions](cell-contexts.csv.gz),
[bulk observations and decisions](bulk-contexts.csv.gz),
[generated shared cell thresholds](cell-thresholds.csv),
[machine-readable counts](summary.json),
[input/source/output hashes](run-manifest.json).

These remain adult-scope research parameters. No ground-truth toxicity,
peptide-presentation calibration or pediatric normal-panel validation is
claimed by a passing candidate or passing software tests.
