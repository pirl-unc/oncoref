# Current target decisions: PRAME, CTAG2 and TBXT

Reviewed **2026-10-27**, against pinned **HPA v23** measurements and the
[current complete-reference funnel](../cta-harmonized-20261027/index.md).
That funnel selects **343 protein-equivalence groups / 358 gene loci**, with
**17,824 Unsafe** and **1,696 Neither** among 19,863 groups. Previous CTA labels
are comparison annotations, not admission criteria or the candidate universe.

Here **Unsafe means an expression-policy exclusion, not demonstrated clinical
toxicity**. A gene-level expression screen and a peptide/HLA/receptor-specific
clinical judgment answer different questions. The latter must retain its scope.

## Measurements responsible for PRAME and CTAG2 flags

The first boundary is the inclusive candidate ceiling; reaching twice that
boundary supplies the stronger expression exclusion. RNA units below are native
HPA nTPM. Supported cluster observations require at least 20 reads.

| Gene | Normal context | RNA | Reads | Candidate / stronger boundary | Current contribution |
|---|---|---:|---:|---|---|
| PRAME | Heart muscle fibroblasts | 12.0 | 72 | 2 / 4 | Strong RNA exclusion |
| PRAME | Colon undifferentiated cells | 20.3 | 96 | 4 / 8 | Strong RNA exclusion |
| PRAME | Adrenal bulk | 3.6 | — | 2 / 4 | Candidate hold only |
| PRAME | Vascular fibroblasts | 2.8 | 22 | 2 / 4 | Candidate hold only |
| PRAME | Heart bulk | 0.3 | — | 16 / 32 | Pass |
| CTAG2 / LAGE-1 | Bronchial basal respiratory cells | 24.9 | 58 | 2 / 4 | Strong RNA exclusion |
| CTAG2 / LAGE-1 | Bronchial club cells | 17.0 | 197 | 2 / 4 | Strong RNA exclusion |
| CTAG2 / LAGE-1 | Bronchial ciliated cells | 8.0 | 39 | 2 / 4 | Strong RNA exclusion |
| CTAG2 / LAGE-1 | Kidney bulk | 3.4 | — | 2 / 4 | Candidate hold only |
| CTAG2 / LAGE-1 | Adipose macrophages | 4.6 | 66 | 4 / 8 | Candidate hold only |
| CTAG2 / LAGE-1 | Adipose NK cells | 5.4 | 22 | 4 / 8 | Candidate hold only |
| CTAG2 / LAGE-1 | Heart bulk | 5.1 | — | 16 / 32 | Pass |

PRAME has Enhanced-negative heart IHC for **cardiomyocytes**, not fibroblasts.
Negative colon glandular-cell staining does not establish a negative in the
undifferentiated cluster. CTAG2 has Enhanced-negative broad respiratory-epithelial
IHC, but the matcher does not certify coverage of each basal/club/ciliated subset.
The high heart-bulk ceiling already reflects matched negative cardiomyocyte IHC.
**CTAG2's heart measurement does not cause its current exclusion.**

PRAME also has an imported Accepted normal-protein finding for **adrenal cortex**.
The underlying primary abstract reports low protein expression. Its stored review
explicitly lacks full methods and donor-level assessment; it is not evidence of
treatment injury. This older acceptance path is less structured than the new
external-IHC validation criteria and needs a uniform quality review across all
accepted findings before being treated as an adjudicated protein threshold.
[Lezcano et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33890816/).

## How IMA203 changes the interpretation

The published phase 1 trial included **40 infused patients**. Its specific
PRAME–HLA-A*02:01 target was unquantifiable in **41 normal tissue types**.
Normal-cell testing found a very weak renal epithelial response in **1/2 donors**;
the paper reported no specific renal toxicity at that cutoff. These are substantial
context-specific data, not a guarantee for every PRAME-directed therapy.
[Wermke et al., 2025](https://www.nature.com/articles/s41591-025-03650-6).

The June 2026 sponsor update describes **33 melanoma patients** at the recommended
dose with manageable tolerability. It is a subset/follow-up report, not 33 independent
additional patients; CRS, neurotoxicity and conditioning toxicities still occurred.
[ASCO 2026 sponsor report](https://www.sec.gov/Archives/edgar/data/1809196/000095010326008113/dp247360_ex9903.htm).

The appropriate interpretation is **PRAME: clinically supported targeting in the
IMA203 context, with normal-expression flags**. The expression-policy label alone
understates that evidence. Clinical evidence is stored categorically and numerically
in [clinical-context-evidence.csv](clinical-context-evidence.csv), separately from
the expression decision; unreviewed products/genes are not coded as negative.

## TBXT and the corrected antibody rule

TBXT now passes the uniform strict-RNA fallback. An uncertain antibody disables
negative-IHC threshold increases but no longer adds a second gene-wide veto.
Removing that veto admits 52 groups, including TBXT; reliable/context-qualified
positive protein findings and the RNA exclusions still apply to every gene.

The 54-1 study specifically challenges the HPA macrophage staining and reports
negative normal lung. This supports a reagent-specific explanation, not universal
adult absence. Independent adult intervertebral-disc findings remain a review flag:
the 2017 series used RNA (60 specimens), IHC (56) and flow cytometry (5).
[Hamilton et al., 2015](https://www.oncotarget.com/article/3086/text/),
[Richardson et al., 2017](https://www.nature.com/articles/s41598-017-01567-w).
The [full antibody comparison](../brachyury-antibody-20261027/index.md) retains
reagent identity, compartment disagreement and missing coverage. No external
antibody has automatically been assigned an HPA Enhanced-negative grade.

## Which rule changes are justified

1. **Applied:** eliminate the redundant gene-wide antibody veto; retain the
   lowest RNA ceilings when IHC is unusable. Keep review flags in selected rows.
2. **Required before protein adjudication:** use the same explicit reagent,
   specimen, donor, cellular localization and specificity criteria for older
   accepted findings and newly collected IHC. A negative in one cell type cannot
   erase a finding in another cell type of the same organ.
3. **Clinical interpretation:** record peptide, HLA, product, exposure, follow-up,
   tissue assay coverage and attributable organ toxicities. Favorable evidence
   can qualify that target context without declaring the entire protein safe.
   This audit adds evidence, not a new per-gene admission override.

For scale only: PRAME's largest RNA-to-ceiling ratio is **6**. Multiplying the
entire linked grid by six would clear its RNA gates, but would leave the adrenal
protein veto. The next power-of-two grid starts at vital/missing **16 nTPM**
(8× the current 2), and still leaves CTAG2's 24.9-nTPM bronchial observation above
the candidate ceiling. Neither choice is clinically validated, and neither has
been substituted into the current funnel merely to admit PRAME. A numeric
threshold fitted to an outcome cannot itself establish that outcome's safety.

## Reproducible evidence

Run `python scripts/audit_cta_target_evidence.py` after the current funnel builder.
[Raw HPA clusters](hpa-cluster-observations.csv),
[raw curated HPA IHC](hpa-ihc-observations.csv),
[assessed cluster contexts](assessed-cluster-contexts.csv),
[assessed bulk contexts](assessed-bulk-contexts.csv),
[all assay peaks](all-assay-peaks.csv),
[current decisions](current-decisions.csv), and [receipts](manifest.json)
retain the underlying values and the rule that consumed them. The source review
does not change released classifications or the forbidden sequence reference.
