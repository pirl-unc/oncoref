# Tissue specific CTA screening and reproductive immune contexts

Immune cells sampled in an explicitly reproductive, placental, fetal or embryonic
context are **neutral for normal-expression exclusion** in this revision. The
same immune lineage in blood, heart, lung or another nonreproductive context
remains screened. The observations keep their measured values and source labels;
neutral observations do not supply positive CTA admission or a claim of protein
uptake. This is a selection convention, not a claim that local immune injury
would be harmless.

All **20,097** reference genes were reprocessed before candidate selection.
Of the previous 426 screened HPA clusters, **27 become neutral and 399 remain
screened**. All 426 still count toward source coverage. Three placental Hofbauer
clusters were already reproductive under the earlier anatomy mapping.

| Newly neutral sampling site | Clusters |
|---|---:|
| Breast | 8 |
| Endometrium | 2 |
| Fallopian tube | 9 |
| Ovary | 1 |
| Placenta | 3 |
| Prostate | 3 |
| Testis | 1 |

The fetal/embryonic rule requires an explicit stage label and supported immune
identity. It does not add a fetal atlas to the current normal panel or extend
the exemption to normal childhood tissues. Broad organ names do not invent
immune identity. Source immune-group metadata or the reviewed immune anatomy
mapping is required.

## Candidate counts

| Normal-expression scenario after the context revision | Original tumor recurrence rule | At least one tumor sample ≥30 pTPM | Requested identities |
|---|---:|---:|---:|
| Original uniform RNA limits | 457 | 364 | 10/10 |
| Tissue-specific limits with local IHC certainty | 201 | 167 | 8/10 |
| Same tissue-specific screen, but protected cell clusters allowed below 30 nTPM | 203 | 169 | 8/10 |

The original recurrence rule is ≥1 pTPM in ≥3 samples and ≥5% of one HPA cancer
cohort. The broader uniform set gains ANKRD7, TNP1, OBP2A, PRM2, SMCP, XAGE3,
PRM1, SLC25A21, C17orf99 and LELP1 compared with the archived 447. Membership
still requires reproductive or developmental association, tumor expression,
specified RNA coverage, and the remaining normal-expression gates.

The tissue-specific screen excludes PRAME and CTAG2/LAGE-1. The other eight
requested identities pass it. The 201 retain historical expression labels of
82 CTA, 6 Unsafe and 113 Neither; the 167 retain 71 CTA, 7 Unsafe and 89 Neither.
Those are annotations from the older expression policy, not recalculated
categories under the new context rule. In particular, a historical RNA veto
may have come from a newly neutral context.

The archived 447-gene API and forbidden-sequence reference are not overwritten.
This revision is available in [gene-comparison.csv.gz](gene-comparison.csv.gz),
with all gates, old membership and context-specific results retained.

## Heart and airway versus colon thresholds

Heart, bronchus and lung are in the vital/critical tier. Colon is in the
protected tier. The main tissue-cluster rule is:

| Tissue-cluster tier | Missing or unmatched IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital or critical | <2 | <2.5 | <3 | <4 |
| Protected | <5 | <6.25 | <7.5 | <10 |
| Manageable | <10 | <12.5 | <15 | <20 |

Values are native HPA nTPM, with ≥20 gene reads required. IHC adjustment requires
the same tissue and cell label, a measured reliable reproductive-positive
profile, and no recorded antibody concern. The bulk, immune, corneal and protein
rules remain those in the [complete threshold table](../cta-threshold-details-20261027/index.md).
All numerical choices remain exploratory research parameters.

The separate protected-30 sensitivity allows every protected-tier cluster below
30 nTPM, without gene-name exceptions. It keeps heart, airway, bulk, immune,
protein and coverage rules unchanged. It adds only two genes to either main
tissue-specific result. Reference-review cutoffs are not changed by this
candidate-admission sensitivity.

**PRAME is not a colon-only signal.** Its colon undifferentiated cluster is
20.3 nTPM with 96 reads, but heart fibroblasts independently reach **12.0 nTPM
with 72 reads**, and vascular fibroblasts reach 2.8 nTPM with 22 reads. Negative
cardiomyocyte staining does not cover fibroblasts. Relaxing the protected/colon
threshold therefore cannot make PRAME pass the stricter heart criterion.
Colon is a different relative risk tier, not a harmless compartment; its
undifferentiated cluster label does not by itself establish which progenitor
or differentiated population is expressing the transcript.

## CTAG2 and CTAG1A or CTAG1B

The pinned HPA v23 measurements distinguish these genes:

| Gene | Heart bulk consensus nTPM | Highest heart cluster nTPM | Highest bronchial cluster nTPM | Highest lung cluster nTPM |
|---|---:|---:|---:|---:|
| CTAG2 / LAGE-1 | 5.1 | 0.0 | 24.9 | 0.0 |
| CTAG1A | 0.0 | 0.0 | 0.0 | 0.0 |
| CTAG1B | 0.0 | 0.0 | 0.0 | 0.0 |

CTAG1A and CTAG1B also have zero gene reads in the imported heart, bronchial,
lung and colon clusters. This is source-specific non-detection, not proof of
absence in all donors or physiological states. CTAG2's heart bulk and single-cell
results are different assays/samples and should not be collapsed into a claim
that the bulk signal is absent. Reviewed heart donors reach 13.5 nTPM.
[HPA v23 source files](https://v23.proteinatlas.org/about/download),
[raw named-gene measurements](named-cell-measurements.csv),
[heart evidence review](../ctag2-heart-evidence.md).

CTAG2 has supported bronchial signals in **basal respiratory cells (24.9 nTPM,
58 reads), club cells (17.0, 197 reads), and ciliated cells (8.0, 39 reads)**.
Airway basal cells contribute to epithelial maintenance and repair. If a
therapy recognized an antigen presented by these normal cells, injury could
impair airway integrity or regeneration. That is a plausible targeting risk,
not demonstrated CTAG2-mediated toxicity from the RNA result alone.
[Primary airway basal-cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2714281/).

HPA reports Enhanced-confidence negative CTAG2 IHC in the broader respiratory
epithelial category. This weighs against readily detectable protein by that
assay, but does not establish absence of the relevant peptide-HLA complex in
the RNA-positive basal or club populations. Source read assignment, transcript
identity and relevant normal-cell protein/presentation remain unresolved.

Low CTAG1A/B expression does not clear every NY-ESO-1 targeting receptor:
**lete-cel recognizes an epitope shared by NY-ESO-1 and LAGE-1A**. Consequently,
CTAG2's normal-tissue evidence matters to receptors with that recognition scope.
This is a receptor/target-identity consideration, not a claim that every NY-ESO-1
receptor recognizes every LAGE-1 isoform.
[Primary lete-cel report](https://pmc.ncbi.nlm.nih.gov/articles/PMC11788651/).

A separate nine-patient NY-ESO-1 TCR study reported one grade-3 lung injury with
TCR-T-cell infiltration. The authors proposed cytokine-associated capillary
leak or recognition of pulmonary tumor micrometastases as explanations; they
did not demonstrate CTAG2 expression as the cause. It therefore cannot be used
as proof of basal-cell targeting, or as evidence that pulmonary toxicity has
never occurred in this therapeutic class.
[Ishihara et al. 2022](https://jitc.bmj.com/content/10/6/e003811).

## CT83 and MAGE-A4 under the neutral-context rule

- **CT83:** testicular macrophage RNA at 8.1 nTPM/13 reads is neutral. Its largest
  remaining raw tissue-cluster value is 0.9 nTPM/6 reads in bronchial ciliated
  cells; its supported maximum is 0.6 nTPM/29 reads in salivary mucus glandular
  cells. Its largest nonreproductive bulk value remains 2.6 nTPM in salivary gland.
- **MAGEA4:** placental T-cell RNA at 0.8 nTPM/1 read is neutral, as is the
  testicular macrophage observation. Remaining screened clusters report 0.0
  nTPM. The largest nonreproductive bulk value remains 0.6 nTPM in lymph node.

Phagocytosis or ambient RNA could contribute to an immune-cell RNA observation,
but the source measurements do not identify that mechanism. Neutrality is
applied from explicit location and lineage without asserting why the RNA was
detected or erasing the measurement.

## Reproduction

Run:

```bash
python scripts/cta_quantitative_threshold_report.py \
  --neutral-reproductive-immune --output /path/to/empty-directory
```

The full-genome context reprocessing precedes the IHC comparison; it is not
restricted to the previous 447 genes. The expensive local IHC comparison then
covers the 490-gene union that passes the remaining gates under either tumor
rule. [Source and output hashes](run-manifest.json),
[all context-reprocessed facts](context-reprocessed-facts.csv.gz),
[cell-level decisions](cell-contexts.csv.gz),
[named-gene decisions](named-gene-comparison.csv).
