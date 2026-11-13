# Published CTAs, exclusive antigen sets and normal-expression peaks

**MAGE-A4, PRAME, NY-ESO-1, XAGE-1 and CT83 remain CTA.** These are the current draft's classifications, after the patient-identity correction, using the fivefold single-tumor rule.

[Search and filter every protein](report.html) · [113 familiar CTA-family proteins and studied targets](familiar-ctas.md) · [All exclusive-set members](exclusive-sets.md) · [The highest 20 in every category and unit](normal-top20.md).

## All curated published target nominations

The OncoRef ledger contains **2,233 protein records with a curated published CA/CTA/oncofetal target nomination from 10 sources**. This includes broad catalog nominations and experimental target papers. It does not label every entry as an experimentally validated CTA. Reproductive-only background lists are excluded from this nomination flag, while their citations remain in the main catalog.

| Current category | Published nominations |
|---|---:|
| CTA | 408 |
| Cancer Antigen only | 24 |
| Reproductive Antigen only | 4 |
| Other | 52 |
| Unsafe | 1745 |

[Every published protein](published-targets.csv.gz) · [Every exact-member paper claim and table receipt](publication-claims.csv.gz). The [familiar-name table](familiar-ctas.md) shows why MAGEC2 is Other:Review, CT45A3 is Other:Sequence, and MAGEA3/MAGEA12, CTCFL/BORIS and HORMAD1 are Unsafe. Published nomination does not override normal-expression evidence or an identity/translation hold.

The 113 familiar CTA-family proteins and studied targets divide into **65 CTA, 43 Unsafe, three Other, one Cancer Antigen only and one Reproductive Antigen only**.

These sources supply the nomination flag. Counts are current protein records with an exact-member claim, not the original papers' gene totals; sources overlap. The claim download retains each source's role and table locations.

| Paper | Matching protein records |
|---|---:|
| [Bradley et al. (2020)](https://doi.org/10.1038/s41467-020-19141-w) | 9 |
| [Bruggeman et al. (2018)](https://doi.org/10.1038/s41388-018-0357-2) | 692 |
| [Carter et al. (2023)](https://doi.org/10.1136/jitc-2023-007935) | 103 |
| [Chang et al. (2019)](https://doi.org/10.1002/cam4.2223) | 985 |
| [Garg et al. (2007)](https://doi.org/10.1158/1078-0432.CCR-06-2340) | 1 |
| [Jamin et al. (2021)](https://doi.org/10.1002/1878-0261.12900) | 469 |
| [Loriot et al. (2025)](https://doi.org/10.1371/journal.pgen.1011734) | 202 |
| [Seager et al. (2024)](https://doi.org/10.1186/s12967-024-04918-0) | 15 |
| [Wang et al. (2016)](https://doi.org/10.1038/ncomms10499) | 988 |
| [da Silva et al. (2017)](https://doi.org/10.18632/oncotarget.21715) | 724 |

## Cancer-only and reproductive-only

**260 Cancer Antigen only protein records cover 261 genes.** They include 147 olfactory receptors, 32 keratin-associated proteins and 81 other records. PPIAL4E/F is one protein record with two member genes. [Complete list and values](exclusive-sets.md#cancer-antigen-only-260-proteins-261-genes).

**22 proteins are Reproductive Antigen only.** Fifteen have unresolved patient independence for potentially qualifying repeated tumor measurements; seven do not meet either cancer RNA rule. Their category therefore does not mean “never expressed in cancer.” [Complete list](exclusive-sets.md#reproductive-antigen-only-22-proteins-22-genes).

## Largest normal values in the three antigen sets

These peaks **include Disposable tissues** and exclude reproductive/developmental contexts. The primary table uses native TPM. The complete ranking has separate nTPM and CP10K tables, plus Other and Unsafe. RNA values are attributed to protein records through their member genes; they are not protein-abundance measurements.

| Category | Protein | Highest normal TPM | Context | Source |
|---|---|---:|---|---|
| CTA | DEFB132 | 24.7 | adipose tissue | HPA tissue mean |
| Cancer Antigen only | LCE3C | 21.4 | tonsil | HPA tissue mean |
| Reproductive Antigen only | VN1R2 | 7.794 | blood / follicular helper T cells | DICE third-highest |

The normal screen permits different RNA ceilings by cell/tissue tier and matched IHC evidence, so a high normal value need not itself create Unsafe. “Highest” here means the highest supported source statistic: tissue means, GTEx donor medians, cell summaries or DICE third-highest TPM. It is not the largest value in any individual healthy donor. DICE unpaired-member upper bounds retain that label in the download. HPA cell RNA requires ≥20 reads; low-read raw signals do not enter this supported ranking. No conversion or pooled ranking is made between TPM, nTPM and CP10K.

PRAME has **20.3 nTPM in undifferentiated colon cells**, and CTAG2/LAGE-1 has **24.9 nTPM in bronchial basal cells**. Both retain CTA status under a **32-nTPM ceiling with Enhanced matched negative IHC**. These are concrete consequences of the current IHC-adjusted rule. The rankings display the source-scale ceiling and IHC grade where assigned.

The full-inventory expansion quarantines conflicting duplicate normal locus/context vectors and collapses identical duplicates; [every affected source row](duplicate-normal-observations.csv.gz) is retained. [Incomplete protein-member observations](incomplete-normal-contexts.csv.gz), including unbounded DICE group estimates, do not enter a quantitative maximum. Infinite upper bounds represent missing coverage, not measured expression. The resulting non-Disposable native-TPM maxima reproduce the existing comparator for all shared-eligible proteins. Classification outputs are unchanged. The hair-follicle mapping and coverage work remains tracked in [#687](https://github.com/pirl-unc/oncoref/issues/687).

Reproduce with `scripts/audit_target_category_profiles.py`, then `scripts/render_target_category_profiles.py`. [Provenance](manifest.json) · [All normal peaks](normal-peaks.csv.gz) · [Separate normal/Disposable source peaks](normal-and-disposable-source-peaks.csv.gz).
