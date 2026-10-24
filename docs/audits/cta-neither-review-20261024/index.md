# CTA and Neither review: high normal RNA

**399 CTA genes remain. 253 previous Neither genes and 37 previous CTAs now enter the forbidden background. No remaining Neither gene has a measured non-permitted/unsited value above 30 TPM or nTPM in the audited sources.**

[Search genes](explore.html) · [Review findings and fixes](review.md) · [Candidate follow-up](follow-up.md) · [All changes](changed-genes.csv) · [Previous report](../cta-normal-cells-20261023/index.md)

## What changed?

The v1.17 rules put normal-cell RNA conflicts into Neither regardless of magnitude. The tissue-cluster rule also compared normal RNA with the reproductive maximum, allowing large normal signals to escape review when reproductive expression was still greater. The new absolute rule is independent of that ratio and of IHC confidence.

**Any observed normal RNA >30 TPM or >30 nTPM outside explicitly mapped Reproductive/Disposable contexts yields Unsafe / High RNA.** Exactly 30 does not trigger this additional rule. It applies to both former CTA and Neither genes, including Variable. Existing stricter protein, RNA and safety exclusions remain. No candidate is promoted and no clinical-target exception is made.

Unsafe here means a research exclusion and membership in the forbidden background, not demonstrated clinical toxicity. Unresolved anatomy does not establish permission: 1 newly excluded genes have only unsited high findings in this audit. Low-read clusters, ambiguous loci and possible contamination remain matters for adjudication; they do not silently waive the requested hard ceiling.

| Scope | CTA | Unsafe | Neither |
|---|---:|---:|---:|
| 2,540 mapped CTA nominations | 399 | 1,580 | 561 |
| 20,097 canonical human genes | 399 | 17,699 | 1,999 |

Across remaining Neither genes the largest audited non-permitted/unsited observations are **29.9398 TPM** and **30 nTPM**. Missing observations remain missing. High reproductive/disposable values are permitted and therefore can still occur in Neither.

## Important targets

**MAGEA3 and ROPN1 leave the retained benchmark.** MAGEA3 has 53.1 nTPM and 251 reads in a colon undifferentiated-cell cluster. ROPN1 has 32.9 nTPM and 496 reads in salivary serous cells. Both cross the general rule; these are RNA-based research exclusions, not established target-mediated toxicity.

The other **19 of the previous 21** benchmark targets remain, including XAGE1A/B, CTAG1A/B, PRAME, MAGEA4, CT83, SAGE1, CTCFL, GAGE1 and SSX4/5. CT55 remains retained. ACRBP, TEX101, GAGE12F, PAGE2B and TRIP13 now leave Neither because the expanded audit finds high normal RNA. HORMAD1, TDRD1, BRDT, DAZL, FOXR2 and CTAG2 remain Neither with explicit reconsideration questions.

| Watched gene | Current result | Maximum non-permitted/unsited nTPM | Maximum non-permitted/unsited TPM |
|---|---|---:|---:|
| [ACRBP](follow-up.md#acrbp) | Unsafe / High RNA | 5316.5 | 56.543 |
| [HORMAD1](follow-up.md#hormad1) | Neither / Cellular | 14.8 | 19.733 |
| [TDRD1](follow-up.md#tdrd1) | Neither / Cellular | 10.8 | 6.434 |
| [BRDT](follow-up.md#brdt) | Neither / Cellular | 17.7 | 6 |
| [TEX101](follow-up.md#tex101) | Unsafe / High RNA | 119.722 | 74.9 |
| [DAZL](follow-up.md#dazl) | Neither / Cellular | 2.7 | 8.913 |
| [FOXR2](follow-up.md#foxr2) | Neither / Antibody | 0 | 0.513 |
| [CTAG2](follow-up.md#ctag2) | Neither / Variable | 24.9 | 5.075 |
| [GAGE12F](follow-up.md#gage12f) | Unsafe / High RNA | 83.07 | 71.665 |
| [PAGE2B](follow-up.md#page2b) | Unsafe / High RNA | 35 | 8.043 |
| [TRIP13](follow-up.md#trip13) | Unsafe / High RNA | 88.8 | 8.639 |
| [CTAGE1](follow-up.md#ctage1) | Neither / Antibody | 1.9 | 0.388 |
| [ZNF730](follow-up.md#znf730) | Neither / Cellular | 11.3 | 7 |
| [MAGEA3](follow-up.md#magea3) | Unsafe / High RNA | 53.1 | 0.166 |
| [ROPN1](follow-up.md#ropn1) | Unsafe / High RNA | 32.9 | 0.595 |

The [follow-up record](follow-up.md) contains measurements, uncertainties and evidence needed for reconsideration for each gene. It is a review queue, not an allowlist. [All important genes](important-ctas.csv) · [Important complete-protein decisions](important-proteins.csv).

## Coverage and units

The audit checks **5,448,388 RNA observations across 2,705 genes**: all previous CTA and Neither genes plus the established important-target audit. Other Unsafe assignments remain inherited; they are not newly cleared. Sources are HPA bulk consensus, cell-type and tissue-cluster RNA, sorted immune summaries, individual HPA immune specimens, original DICE distributions, and the four reviewed CTAG2 heart RNA specimens. Both original TPM and nTPM columns from the HPA immune-sample table are checked separately.

- No averaging across donors, cells, tissues, studies or units can hide a value >30. A single reported high value suffices under this conservative policy.
- The older HPA cell-type table preserves its unsited contexts; it is not assigned an organ by guessing. Immune cells sampled in reproductive/disposable organs retain their existing hematopoietic mapping.
- HPA summaries and their original data overlap. They are not counted as independent biological replication.
- DICE slots are not stable donor identifiers across genes. RNA and IHC specimens are not assumed paired.
- Corneal CP10K is **not TPM** and does not enter the numerical ceiling. The existing corneal rule and uncertainty remain; developmental/culture data are not pooled into adult-normal evidence.
- This is measured-source coverage, not a claim that every normal tissue, physiological state or donor has been sampled. Human thresholds are not transferred to mouse or dog.

The following source counts overlap and must not be added as independent genes or cohorts:

| Source with newly excluding observations | Native unit | Changed genes | Observations >30 |
|---|---|---:|---:|
| dice_b_cell_naive | TPM | 16 | 695 |
| dice_cd4_naive | TPM | 12 | 561 |
| dice_cd4_stim | TPM | 18 | 719 |
| dice_cd8_naive | TPM | 15 | 668 |
| dice_cd8_stim | TPM | 11 | 571 |
| dice_m2 | TPM | 9 | 147 |
| dice_monocytes | TPM | 7 | 202 |
| dice_nk | TPM | 24 | 1,115 |
| dice_tfh | TPM | 17 | 780 |
| dice_th1 | TPM | 17 | 690 |
| dice_th17 | TPM | 18 | 765 |
| dice_th2 | TPM | 16 | 767 |
| dice_thstar | TPM | 14 | 740 |
| dice_treg_mem | TPM | 17 | 762 |
| dice_treg_naive | TPM | 17 | 659 |
| hpa_immune_cell | nTPM | 17 | 39 |
| hpa_immune_cell_samples | TPM | 26 | 184 |
| hpa_immune_cell_samples | nTPM | 29 | 212 |
| hpa_immune_monaco | TPM | 20 | 72 |
| hpa_immune_schmiedel | TPM | 21 | 93 |
| hpa_single_cell | nTPM | 113 | 289 |
| hpa_single_cell_tissue | nTPM | 236 | 2,325 |

## Files and API

[All genes](gene-results.csv.gz) · [CTA](cta-genes.csv.gz) · [Neither](neither-genes.csv.gz) · [Unsafe](unsafe-genes.csv.gz) · [All high observations](high-normal-rna.csv.gz) · [Context maxima including permitted tissues](rna-context-maxima.csv.gz) · [Source provenance](source-provenance.json) · [Run manifest](run-manifest.json)

```python
from oncoref import target_sets, normal_rna_review
target_sets.results()                 # current decisions and native-unit maxima
normal_rna_review.followup_reviews()   # watched genes, questions and evidence
```

Complete protein sets are **723 CTA, 126,368 Unsafe and 2,894 Neither**. Shared exact sequences use Unsafe > Neither > CTA. The high-RNA exclusions also update incomplete translations in the additional forbidden background. [Protein assets and checksums](protein-assets.json) · [Exclusion recipe](exclusion-basis.json). Expression is gene-level and does not establish isoform-specific expression.

Reproduce from the same shared cache:

```bash
python scripts/review_neither_ctas.py --assets /path/to/new-assets --previous-sequences /path/to/pinned-ensembl115-assets
python scripts/render_neither_cta_review.py --assets /path/to/new-assets --publish-dataset
```

Policy: `oncoref.target-sets.2026-10-24.1`. Baseline: `oncoref.target-sets.2026-10-23.2`. The prior published report remains unchanged.
