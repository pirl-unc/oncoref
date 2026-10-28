# Structured external staining evidence

6,236 observations across 193 matched genes are represented as
**detected, not detected, mixed, reported cross-reactivity or not assessed**.
Intensity, localization, explicitly reported percentages, H-scores, specimen counts,
age/stage and specimen state have separate fields. Missing numbers remain missing.
Source statements and table/figure locations accompany each observation.

The 201-product MSVA collection supplies additional panels for **CD70, CLDN6, CTLA4, DNTT, INSM1, MYOD1, TBXT, TIGIT**
among proteins with missing or uncertain HPA evidence. Paragraph observations have
categorical results even when anatomical mapping still requires review.
10 additional findings were reviewed against publication text or product images.

## How an external result can supersede an Uncertain HPA observation

Use the same rule for every gene: reviewed human adult normal-tissue IHC, independent
identified reagent, a positive control, specificity validation in IHC, the exact
tissue and cell compartment, and at least two normal donors for a negative result.
Two donors is a provisional review-policy minimum, not a validated absence cutoff.
Documented positives require at least one normal donor and the same context and
specificity checks. Validated positive/negative disagreement remains a conflict.
Western-blot specificity does not automatically validate tissue IHC. Supplier
panels do not inherit HPA's Approved, Supported or Enhanced grades.

**0 observations currently meet every requirement.**
Consequently this evidence collection produces no automatic negative-IHC overrides.
Each failure is machine-readable in `ineligible_reasons`; collection is not clearance.

## Findings that matter to the current review

- TBXT: 54-1 reports negative normal lung, while HMV321 reports negative organ panels.
  Missing donor/cell-level documentation prevents a context-specific override.
  Adult disc staining remains recorded separately, including surgical specimen state
  and the companion flow assay (5 specimens; mean 76.2% positive cells).
- CT83: CL4762 staining was reported in 76/84 tumor-adjacent gastric specimens;
  median H-score 50, range 0–150 on a 0–300 scale. This is a review flag for a named
  control, not a healthy-donor panel and not grounds for pretending the signal is zero.
- ADAM30: cortical neuronal staining with a custom Interchim reagent is a positive
  finding; the separate RNA donor counts are not borrowed for IHC.
- CDX4 and GOT1L1 have additional positive tissue images, with specificity and reagent
  independence unresolved. RNASE11 and FOXR2 findings are reproductive/placental.

The detailed rows link each finding to its original source.
[Reviewed findings](reviewed-publication-and-image-evidence.csv),
[all assessed evidence](assessed-evidence.csv.gz),
[catalog and gene identity mapping](external-panel-inventory.csv),
[original categorical extraction](external-panel-observations.csv.gz),
[source receipts](source-receipts.json), [summary](evidence-summary.json).

Reproduce with `python scripts/collect_external_ihc.py` followed by
`python scripts/curate_external_ihc.py`. These are research evidence annotations;
released classifications and sequence references are unchanged.
