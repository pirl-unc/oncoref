# CTA expression sources and threshold comparisons

The frozen quantitative proposal contains **447 genes**. Replacing its tumor
prevalence rule with **at least one measured tumor sample at ≥30 pTPM** yields
**352 genes**, including all ten requested identities. This replacement removes
126 and adds 31; it is not a nested subset. The IHC-dependent exclusion scenario
implemented here retains **187** under the original tumor rule, or **150** under
the ≥30 pTPM alternative. Both stricter scenarios lose PRAME and CTAG2/LAGE-1.

These are research selection rules, not validated toxicity thresholds. The
447-gene artifact and released expression/sequence-reference policies remain
unchanged. The comparison applies existing tissue/IHC rules consistently;
it does not fit new per-gene exceptions to retain the requested controls.

| Normal RNA policy | Tumor admission | Genes | Requested identities | Existing CTA / Unsafe / Neither annotations |
|---|---|---:|---:|---|
| Original uniform limits | ≥1 pTPM in ≥3 samples and ≥5% of one cohort | 447 | 10/10 | 127 / 110 / 210 |
| Original uniform limits | ≥30 pTPM in at least one sample | 352 | 10/10 | 101 / 88 / 163 |
| Additional local IHC-dependent limits | Original recurrence rule | 187 | 8/10 | 80 / 1 / 106 |
| Additional local IHC-dependent limits | ≥30 pTPM in at least one sample | 150 | 8/10 | 67 / 0 / 83 |

Existing risk labels are annotations from the separate expression policy; zero
genes labeled Unsafe does not establish clinical safety. No reference-review
finding in this audit adds a sequence to the forbidden-protein reference.

## RNA units and cancer cohorts

**pTPM means protein-coding transcripts per million.** HPA removes noncoding
transcripts from the normalization denominator and rescales the retained RNA
abundances to one million per sample. This is RNA abundance, not a protein
measurement. **nTPM** additionally uses between-sample TMM normalization. The
corneal **CP10K** statistic is count-normalized RNA per 10,000 counts and is not
TPM. [HPA normalization methods](https://www.proteinatlas.org/humanproteome/cancer/method/transcriptomics).

The tumor gate currently queries **HPA 25.1**, imported into OncoRef as
`hpa_cancer_rna`, version `v25.1-1`. Its 8,384 nominal samples comprise 6,918
TCGA samples in 21 cohorts and 1,466 samples in ten validation cohorts. The
recurrence rule is evaluated within one cancer type and cohort; TCGA and
validation are separate. It does not query all OncoRef expression sources.
[HPA cancer source](https://www.proteinatlas.org/humanproteome/cancer/data),
[exact cohort labels and denominators](tumor-cohorts.csv).

- TCGA: BLCA, BRCA, CESC, COAD, GBM, HNSC, KICH, KIRC, KIRP, LIHC, LUAD, LUSC,
  OV, PAAD, PRAD, READ, SKCM, STAD, TGCT, THCA and UCEC.
- Validation: BRCA, COAD, GBM, KIRC, LIHC, LUAD, LUSC, OV, PAAD and READ.

Other registered OncoRef sources include TARGET, Treehouse and disease-specific
cohorts. They are not inputs to either tumor comparison here. An expanded tumor
admission rule could accept valid patient-tumor evidence from those sources,
while retaining the sample, cancer, assay, native unit and reference identity.
Cell lines, computed aggregate cohorts and noncomparable pseudobulk units must
not silently become independent patient observations.

The original 1-pTPM/3-sample/5% rule asks whether tumor expression is reproducible
within a cohort. Its prevalence floor is a design choice, not a safety
requirement, and can reject rare but strongly expressing tumors. The ≥30-pTPM
alternative directly supports that rare-expression use case. A single maximum
remains vulnerable to an outlier or a sample/identity error; recurrence can be
reported without making it a hard admission condition.

It is not a literal inversion of a normal-expression safety threshold:
30 tumor pTPM, 30 cluster nTPM and 30 bulk TPM are different measurements.
Neither the existing normal cutoff nor the proposed tumor cutoff has been
validated as a universal protein or peptide-presentation threshold.

## Developmental source and qualifying lineages

The implemented source is **Tyser et al. 2021, E-MTAB-9388**, through the
**EBI Single Cell Expression Atlas filtered gene-TPM matrix**: 1,195 cells
from **one human Carnegie-stage-7 embryo**, approximately **16–19 days after
fertilization**. Quantification uses Ensembl 106 and Kallisto 0.46.2; filtered
TPM is used without renormalizing the retained genes.
[Primary study](https://www.nature.com/articles/s41586-021-04158-y),
[EBI matrix](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-9388/downloads),
[pinned provenance](../cta-quantitative-20261026/developmental-source.json).

“Lineage” means the study's author-assigned cell-type label, not a prediction
that an adult tumor is embryonic. All eleven groups are processed identically:

| Author label | Cells |
|---|---:|
| Advanced mesoderm | 164 |
| Axial mesoderm | 23 |
| Ectodermal cell | 29 |
| Emergent mesoderm | 185 |
| Endodermal cell | 135 |
| Epiblast cell | 133 |
| Erythrocyte | 32 |
| Hemogenic endothelial progenitor | 111 |
| Nascent mesoderm | 98 |
| Primitive streak | 202 |
| Yolk sac mesoderm | 83 |

The mean includes every cell in the group, including zeros. To qualify, a group
must have ≥20 cells, ≥5 cells with ≥1 TPM, ≥10% such positive cells, and a mean
≥5 TPM. There is no minimum of two embryos: cells from this single embryo do
not provide donor replication. TBXT qualifies in axial mesoderm at **245.09 TPM**,
with **22/23 positive cells**.

Cao/DESCARTES fetal organs, later organogenesis, pediatric normal tissues,
organoids and animal datasets do not contribute to this developmental gate.
Registered fetal corneal data is also not used for admission here. The source
is an embryonic-expression atlas, not a curated oncofetal-antigen list.

## Every threshold in the original 447-gene proposal

The source downloads below are pinned to **HPA v23 for normal tissues** and
**HPA 25.1 for cancer RNA**. A current live HPA page can contain different values.

| Measurement | Actual source | Passing rule and support |
|---|---|---|
| Reproductive association | HPA v23 `rna_tissue_consensus.tsv`, HPA + GTEx consensus | Maximum in testis, ovary or placenta ≥5 nTPM **and greater than** the largest nonreproductive bulk value. Alternative to the developmental route. |
| Developmental association | E-MTAB-9388 EBI filtered TPM and author labels | Mean ≥5 TPM with the four support requirements above. Alternative to reproductive association. |
| Tumor RNA | HPA 25.1 `rna_cancer_sample.tsv.gz`, summarized by OncoRef | ≥1 pTPM in ≥3 measured samples and ≥5% of one cancer/cohort. |
| Normal bulk | HPA v23 `rna_tissue_consensus.tsv` | Every one of 39 nonreproductive tissue estimates ≤10 nTPM. These are tissue consensus estimates, not maximum individual donors. |
| Normal tissue clusters | HPA v23 `rna_single_cell_type_tissue.tsv` | Maximum supported nonpermitted cluster <30 nTPM; support is ≥20 gene reads per cluster. |
| HPA immune RNA | HPA v23 `rna_immune_cell.tsv` | Every one of 19 summaries <5 nTPM: 18 immune-cell types plus total PBMC. |
| Monaco immune RNA | Monaco et al. 2019, distributed as HPA v23 `rna_immune_cell_monaco.tsv` | Every one of 30 summaries <5 **TPM**, the actual native column in the pinned file. |
| DICE immune RNA | Original DICE DB1 distributions, release 2022-02-23, 15 cell types/states | Third-highest reported value in each distribution <5 TPM; exclusion requires ≥3 values at or above 5. The HPA DICE summary is not independent replication. |
| Corneal RNA | Collin et al. 2021, GSE155683 **adult** counts, HCA cell annotations | No cell type with mean ≥1 CP10K in ≥2 sample groups, each supported by ≥5 positive cells and ≥20 counts. Groups are not asserted to be independent donors. |
| Normal protein | HPA v23 `normal_tissue.tsv` curated tissue/cell IHC | Any nonpermitted Low/Medium/High detection excludes from selection, including Uncertain-grade staining. Missing IHC is flagged but does not exclude in the original proposal. |
| RNA coverage | Those exact pinned panels | 50 bulk tissues, 426 nonpermitted clusters, 19 HPA immune summaries, 30 Monaco summaries, 15 DICE distributions and 84 corneal cell/sample groups; required numerical inputs must be present. |

[HPA v23 downloads](https://v23.proteinatlas.org/about/download),
[original DICE](https://dice-database.org/),
[Collin study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8343164/),
[every source URL and hash](../cta-quantitative-20261026/cell-source-receipts.json).

Exactly 10 nTPM passes the bulk gate; exactly 30 nTPM fails the cluster gate.
No supported cluster signal is distinct from a measured zero. Missing source
rows cannot be interpreted as absence. The 11 bulk reproductive exemptions are
breast, cervix, endometrium, epididymis, fallopian tube, ovary, placenta, prostate,
seminal vesicle, testis and vagina. Immune cells in reproductive organs remain
screened immune compartments. Conditional disposable-cell mappings are also
retained in the cell/IHC screen.

There is no publication minimum, named-control override, fixed target count,
continuous-breadth cutoff, donor-percentile gate, or negative-IHC waiver in the
447-gene proposal. Existing donor and breadth annotations remain available.

## Adult versus pediatric interpretation

“Adult” describes the intended normal-expression screening scope. The current
pipeline labels the imported normal-cell observations as adult; it does not
establish comprehensive donor-age coverage or a pediatric age-specific normal
reference. HPA and GTEx bulk tissue estimates, HPA's aggregated normal-cell
datasets and the adult corneal subset cannot establish low expression throughout
childhood. The developmental atlas is a positive nomination source, not a
negative normal reference for children.

For pediatric evaluation, normal childhood expression must remain exclusion
evidence rather than being excused as developmental expression. Age-specific
normal data are needed, particularly for developmental remnants, growing and
maturing tissues and immune development. The current numbers have not been
validated for children. TBXT's adult disc/nucleus-pulposus coverage gap also
remains outside these routine panels; see the
[developmental source assessment](../../developmental-targets.md).

## IHC-dependent RNA comparison

HPA's `Reliability` field describes the **curated gene protein-expression
profile**, informed by antibody validation and other evidence. It is not a
numerical antibody-specific false-negative probability. Enhanced, Supported,
Approved and Uncertain are therefore evidence tiers, not calibrated RNA-to-
protein conversions. [HPA profile reliability](https://v23.proteinatlas.org/about/assays%2Bannotation).

The new comparison keeps the original caps and adds local IHC-dependent gates.
Negative IHC can increase the local allowance only when there is also a reliable
measured reproductive-positive profile and no recorded antibody concern. This
is a quantitative evidence condition applied to every gene, not the list of
ten requested controls. For bulk RNA, all expected IHC contexts in the same
anatomical tissue must be negative; the weakest reliability governs. The
inventory is the full set of tissue/cell annotations in the pinned normal IHC
export. Broad and detailed annotations are not silently treated as equivalent.

| Bulk anatomical tier | Missing, partial or Uncertain IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | ≤2 | ≤4 | ≤6 | ≤10 |
| Protected | ≤6 | ≤8 | ≤10 | ≤10 |
| Manageable | ≤10 | ≤10 | ≤10 | ≤10 |

All values are **nTPM**. This is the existing absolute ladder
(2/4/6/12, 6/8/12/24, 10/20/30/60) capped at the proposal's 10 nTPM. The previous
reproductive-selectivity ratio ladder is not applied in this comparison.

For tissue-cluster RNA, negative staining must match the **same tissue and cell
label**, normalized only for case and whitespace. An organ-level result cannot
clear an unmeasured cell population. The existing illustrative cell model gives:

| HPA tissue-cluster tier | Missing, unmatched or Uncertain IHC | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | <2 | <2.5 | <3 | <4 |
| Protected | <5 | <6.25 | <7.5 | <10 |
| Manageable | <10 | <12.5 | <15 | <20 |

All values are **nTPM**, with ≥20 reads still required. Equality excludes.
For sorted/activated immune cells, no reviewed equivalent HPA IHC matching is
available: the added gate is <2 in the source's native TPM or nTPM. The DICE
three-value requirement remains. Adult cornea retains <1 CP10K with its original
cell/count/recurrence guards. Those assays receive no organ-level IHC relaxation.

The exported stronger **reference-review** cutoffs are separate from candidate
admission: vital/protected/manageable tissue clusters use 5/10/20 nTPM, immune
sources use 5 native units, and cornea uses 2 CP10K. Exact negative IHC multiplies
both lower and upper cell cutoffs by 1, 1.25, 1.5 or 2 for
Missing-or-Uncertain/Approved/Supported/Enhanced. These upper cutoffs only
annotate evidence for review; they do not construct a forbidden sequence set.

Uncertain positive IHC remains a selection hold. Reliable positive IHC is
stronger protein evidence. Missing IHC is unknown and receives the stricter
RNA allowance. Only 1,499 of the 128,104 assessed cell contexts match an IHC
record exactly, and 1,451 qualify for the negative-IHC adjustment.

PRAME fails the added cell limits at colon undifferentiated cells (20.3 nTPM),
heart fibroblasts (12.0) and vascular fibroblasts (2.8). CTAG2 fails at bronchial
basal cells (24.9), club cells (17.0), ciliated cells (8.0), and several immune
populations. HPA does report Enhanced negative CTAG2 IHC in **bronchial
respiratory epithelial cells**, but it does not separately score those basal
and club subtypes. Additional conservative bulk failures are PRAME adrenal
3.6 nTPM and CTAG2 kidney 3.4 nTPM, each against 2 nTPM because detailed local
IHC coverage is incomplete under the full-inventory rule.

## Normal expression for the ten requested identities

These HPA values come from the pinned v23 files. Immune cells in placenta or
testis count as normal immune cells, not as exempt germline expression.

| Gene | Highest nonreproductive bulk consensus, nTPM | Highest normal tissue-cluster RNA, nTPM and reads | Highest supported cluster, nTPM |
|---|---|---|---|
| PRAME | 3.6, adrenal gland | 20.3, colon undifferentiated cells, 96 reads | 20.3 |
| MAGEA4 | 0.6, lymph node | 0.8, placental T cells, 1 read | None with ≥20 reads |
| CTAG1A | 0.1, skin | 0.1, adipose macrophages, 1 read | None with ≥20 reads |
| CTAG1B | 0.2, cerebellum | 0.0 in all screened clusters | None with ≥20 reads |
| XAGE1A | 0.1, adipose/lung/lymph node/skin | 2.2, skin macrophages, 2 reads | None with ≥20 reads |
| XAGE1B | 1.5, liver | 1.2, skin Langerhans cells, 1 read | None with ≥20 reads |
| CTAG2 / LAGE-1 | 5.1, heart muscle | 24.9, bronchial basal respiratory cells, 58 reads | 24.9 |
| SSX2 | 0.2, thyroid gland | Rounded 0.0; some neuronal clusters have 1 read | None with ≥20 reads |
| CT83 | 2.6, salivary gland | 8.1, testicular macrophages, 13 reads | 0.6, salivary mucus glandular cells, 29 reads |
| TBXT / brachyury | 2.5, parathyroid gland | 2.2, skeletal-muscle fibroblasts, 18 reads | 1.2, adipose fibroblasts, 35 reads |

The 24.9 figure is **RNA in HPA bronchus cluster c-14**, which contains 311
cells; 58 is its gene read count, not 58 positive cells. HPA reprocessed the
Lukassen et al. bronchial dataset (`fig11981034`). PRAME's 20.3 is colon cluster
c-8, from the HPA reprocessing of Parikh et al. (`GSE116222`).

The broader imported normal evidence includes higher individual immune values
than the tissue/cluster summaries: **XAGE1B 3.8312 nTPM in plasmacytoid dendritic
cells**, and **TBXT 4.6117 nTPM in naive B cells**, both HPA immune sample data.
They are retained in the [normal-context maxima](named-normal-maxima.csv).
The IHC-dependent comparison still gates HPA and Monaco immune **summaries**,
not individual HPA immune donors; it therefore does not exclude XAGE1B or TBXT
on these two individual values. Replacing a summary with the largest individual
normal measurement would be another policy change, separate from the tumor
maximum comparison.
CTAG2's highest original DICE value is **5.07497 TPM in TH1**; the maximum
within-distribution third-highest value used by the recurrence gate is
**1.646174 TPM**. An isolated high DICE value is therefore not hidden, but it
does not cross the recurrent-expression gate.

CTAG2's **5.1 nTPM heart consensus is not the highest individual donor**.
The reviewed HPA v23 heart values are **13.5, 7.1, 0 and 0 nTPM**. Its 10.1-nTPM
reproductive measurement, 5.1-nTPM bulk maximum and 24.9-nTPM supported cluster
maximum describe three different summaries. The original 447-gene rule uses
the latter two cutoffs; it does not apply an additional donor-maximum exclusion.

## Audit files and reproduction

- [All gene comparisons](gene-comparison.csv.gz)
- [Ten-gene tumor and IHC comparison](named-gene-comparison.csv)
- [All ten-gene normal evidence contexts](named-normal-contexts.csv.gz)
- [Bulk IHC matching, missing contexts and effective limits](bulk-contexts.csv.gz)
- [Cell IHC matching and limits](cell-contexts.csv.gz)
- [Cell exclusion findings](cell-exclusions.csv.gz)
- [Source and output checksums](run-manifest.json)

Run `python scripts/cta_quantitative_threshold_report.py --output /path/to/empty-directory`.
All 20,097 genes retain their baseline gates. The additional IHC evaluation
covers the 478-gene union eligible under either tumor rule; genes already
excluded by unchanged gates cannot be rescued by this exclusion-only layer.
The report records which genes were evaluated, rather than filling unevaluated
normal evidence with passing values.
