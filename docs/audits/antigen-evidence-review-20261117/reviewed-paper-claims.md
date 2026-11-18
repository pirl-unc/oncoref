# Centrally reviewed papers and gene claims

Generated from `oncoref/data/cta-literature-followup.json` (schema 2). Edit that central ledger, not this report. Original nomination tables are immutable source inputs. Reviewed claims override inherited scope for the same paper/gene in the current catalog, API and publication filters.

Gene lists show reviewed claims, not every gene mentioned by a paper. Reproductive expression, translation and functional experiments are separate. A publication use is neither automatic tenfold qualification nor clinical clearance.

## [Moyer et al., 2009](https://doi.org/10.1371/journal.pone.0007682)

Expression of Genes Encoding Multi-Transmembrane Proteins in Specific Primate Taste Cell Populations

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| CALHM3 | CALHM3 | Exact human gene, normal circumvallate taste-bud RNA; Figure 10D | Human circumvallate taste buds | normal_tissue_review |

**CALHM3** (`ENSG00000183128`): Microdissection and semi-quantitative PCR. Source: Human taste-bud results and Figure 10D; Methods: human tissue collection and PCR. Limitations: Human data are RNA from a microdissected compartment. Cell-subtype in situ experiments elsewhere in this paper used macaques. This is not human protein detection or a therapeutic depletion experiment.

## [Böhm et al., 2012](https://doi.org/10.3892/or.2012.1849)

Comparison of tear protein levels in breast cancer patients and healthy controls using a de novo proteomic approach

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| TRIM75 | TRIM75 | Exact human protein label in a pooled tear-MS identification table; peptide and clinical-group attribution unresolved | Human tear fluid, pooled cancer and control samples | normal_tissue_review, translation_context |

**TRIM75** (`ENSG00000250374`): MALDI mass spectrometry. Source: Table II; pooled sampling and MALDI identification methods; Table III does not quantify TRIM75. Limitations: One pooled tear sample per clinical group. Table II merges identifications and does not assign TRIM75 to a pool or supply its peptide sequences. Tear-fluid detection does not localize the producing tissue or demonstrate a tumor origin.

## [Zhu et al., 2021](https://doi.org/10.1093/molehr/gaab012)

Histone demethylase KDM4A overexpression improved the efficiency of corrected human tripronuclear zygote development

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| UBTFL1 | UBTFL1 | Reported human assay rejected for exact-copy coverage: primers map to mouse Sod1, not human UBTFL1 | Reported human zygote assay; primers match mouse Sod1 | identity_review |

**UBTFL1** (`ENSG00000255009`): Primer/reference sequence comparison. Source: Table I Ubtfl1 primers; Figure 5; NCBI XM_021185540.2 and exact forward/reverse-complement cDNA comparisons. Limitations: This finding rejects attribution of this specific assay. It does not reject the other experiments in the paper or establish normal UBTFL1 expression. Mouse Hmgpi and human UBTFL2/3 results cannot be substituted for this copy.

## [Shi et al., 2020/2021](https://doi.org/10.1136/jclinpath-2020-206979)

Identification of potential crucial genes associated with the pathogenesis and prognosis of liver hepatocellular carcinoma

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| TMEM269 | TMEM269 | Exact human gene in a significant differential-expression result table; no focused gene experiment or CTA nomination | Human liver cancer versus comparison tissue | tumor_support |

**TMEM269** (`ENSG00000274386`): TCGA differential-expression reanalysis. Source: Supplementary Table S3, PDF page 101; differential-expression table header on page 1. Limitations: The supplementary RNA comparison is a reanalysis of TCGA, not an independent tumor cohort. Low average RNA expression and a differential-expression statistic do not establish an accessible membrane protein.

## [Qi et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7868764/)

CDCA8, targeted by MYBL2, promotes malignant progression and olaparib insensitivity in ovarian cancer

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| SSU72L1 | SSU72P7 | Exact previous human symbol SSU72P7 in a significant RNA perturbation result table; no selective protein validation | Human ovarian cancer cells following CDCA8 knockdown | tumor_support |

**SSU72L1** (`ENSG00000284438`): RNA differential-expression supplementary table. Source: Figure 4 and Table S2; combined publisher PDF page 45 (supplement page 18). Limitations: The table reports a gene-level RNA result, not unique-read or peptide validation. SSU72 phosphatase experiments are not transferable to SSU72L1 or to SSU72L2-5.

## [Zou et al., 2022](https://pubmed.ncbi.nlm.nih.gov/36074823/)

Translatome and transcriptome co-profiling reveals a role of TPRXs in human zygotic genome activation

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| TPRX1 | TPRX1 | Named human gene in early-embryo translation study; combined TPRX1/2/L functional perturbation, not a single-copy functional validation | Human oocytes and early embryos; combined TPRX perturbation | reproductive_support, translation_support, functional_context |
| TPRX2 | TPRX2 | Named human gene in early-embryo translation study; combined TPRX1/2/L functional perturbation, not a single-copy functional validation | Human oocytes and early embryos; combined TPRX perturbation | reproductive_support, translation_support, functional_context |

**TPRX1** (`ENSG00000178928`): Ribo-lite, mRNA sequencing and combined knockdown. Source: PubMed abstract, PMID 36074823; human oocyte/early-embryo Ribo-lite and mRNA co-profiling. Limitations: Developmental and translation evidence, not tumor expression, antigen presentation or adult tissue safety. The abstract does not resolve every paralog-specific ribosome footprint or assign the combined knockdown phenotype to one copy. No TPM enrichment value is inferred.

**TPRX2** (`ENSG00000259009`): Ribo-lite, mRNA sequencing and combined knockdown. Source: PubMed abstract, PMID 36074823; human oocyte/early-embryo Ribo-lite and mRNA co-profiling. Limitations: Developmental and translation evidence, not tumor expression, antigen presentation or adult tissue safety. The abstract does not resolve every paralog-specific ribosome footprint or assign the combined knockdown phenotype to one copy. No TPM enrichment value is inferred.

## [Characterization and target genes of nine human PRD-like homeobox domain genes expressed exclusively in early embryos](https://pmc.ncbi.nlm.nih.gov/articles/PMC4944136/)

Characterization and target genes of nine human PRD-like homeobox domain genes expressed exclusively in early embryos

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| DPRX | DPRX | Embryo-expressed homeobox regulator | Human preimplantation embryos; hESC overexpression | reproductive_support, functional_context |
| TPRX2 | TPRX2 | Embryo-expressed homeobox regulator | Human preimplantation embryos; hESC overexpression | reproductive_support, functional_context |

**DPRX** (`ENSG00000204595`): Embryo cDNA cloning and transcriptomic perturbation. Source: Cloning methods; Figure 3 regulatory responses. Limitations: Overexpression supports regulatory activity; no tumor validation or endogenous protein quantification.

**TPRX2** (`ENSG00000259009`): Embryo cDNA cloning and transcriptomic perturbation. Source: Cloning methods; Figure 3 regulatory responses. Limitations: Overexpression supports regulatory activity; no tumor validation or endogenous protein quantification.

## [The human PRD-like homeobox gene LEUTX has a central role in embryo genome activation](https://pmc.ncbi.nlm.nih.gov/articles/PMC5087614/)

The human PRD-like homeobox gene LEUTX has a central role in embryo genome activation

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| DPRX | DPRX | Early-embryo-expressed transcriptional repressor | Human embryos, hESC overexpression and HEK-293 reporter assays | reproductive_support, functional_context |

**DPRX** (`ENSG00000204595`): RNA profiling and luciferase reporter assays. Source: Figures 4 and 5G. Limitations: Cancer-cell-line discussion concerns LEUTX; it is not DPRX tumor evidence. Embryo function is inferred from expression and regulatory assays.

## [Comprehensive characterization of the embryonic factor LEUTX](https://pmc.ncbi.nlm.nih.gov/articles/PMC9978639/)

Comprehensive characterization of the embryonic factor LEUTX

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| TPRX2 | TPRX2 / TPRX2P | LEUTX-responsive transcript in an embryonic regulatory program | Human cells with experimentally induced LEUTX | reproductive_support |

**TPRX2** (`ENSG00000259009`): RT-qPCR after LEUTX induction. Source: Figure 4B; TPRX2/TPRX2P results. Limitations: This paper cites the separate combined-TPRX knockdown study; that cited experiment is not a second independent functional validation.

## [Dunwell TL et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28446706/)

A sister of NANOG regulates genes expressed in pre-implantation human development.

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| NANOGNB | NANOGNB | Preimplantation-expressed regulator of an early-embryo transcription program | Human 8-cell/morula RNA; ectopic expression in adult fibroblasts | reproductive_support, functional_context |

**NANOGNB** (`ENSG00000205857`): RNA-seq and transgene perturbation. Source: Results 2.2/2.3; Figures 2/3. Limitations: Fibroblast transfection is an assay context, not endogenous adult expression. No human embryo loss-of-function experiment or tumor validation.

## [Gong et al. (2021)](https://doi.org/10.1038/s41467-021-22695-y)

The RNA landscape of the human placenta in health and disease

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| FOXL3 | WI2-2373I1.1 | Testis-enriched RNA in a cross-tissue comparison | Human testis; GTEx comparison within a placenta-focused paper | reproductive_support |
| SSU72L6 | SSU72P8 | Testis-enriched RNA in a cross-tissue comparison | Human testis; GTEx comparison within a placenta-focused paper | reproductive_support |

**FOXL3** (`ENSG00000248767`): Tissue-enriched transcript table. Source: Supplementary Data 5, Data 5 - tissue-enriched PC, Excel row 1274, Tissue in query = Testis. Limitations: RNA association only; no demonstrated reproductive function, tumor validation or locus-specific protein detection.

**SSU72L6** (`ENSG00000230268`): Tissue-enriched transcript table. Source: Supplementary Data 5, Data 5 - tissue-enriched PC, Excel row 1242, Tissue in query = Testis. Limitations: RNA association only; no demonstrated reproductive function, tumor validation or locus-specific protein detection.

## [Choi MS et al. (2004)](https://pubmed.ncbi.nlm.nih.gov/15034037/)

IFN-gamma stimulates the expression of a novel secretoglobin that regulates chemotactic cell migration and invasion.

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| SCGB1D4 | IIS | Interferon-responsive secretoglobin with reported chemotactic activity | Human lymphoid tissues/cells, including activated states | normal_tissue_review, functional_context |

**SCGB1D4** (`ENSG00000197745`): Expression and cell-migration experiments described in the abstract. Source: PMID 15034037 abstract. Limitations: Full probe specificity was not adjudicated. Reported tissue breadth conflicts with resting HPA; retain induced-expression concern without automatic Unsafe classification.

## [Lu X et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21385388/)

The cytokine-driven regulation of secretoglobins in normal human upper airway and their expression, particularly that of uteroglobin-related protein 1, in chronic rhinosinusitis.

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| SCGB1D4 | SCGB1D4 / IIS | Secretoglobin RNA detected in human sinonasal mucosa and increased by inflammatory cytokines | Human control sinonasal tissue, rhinosinusitis and cytokine-treated nasal explants | normal_tissue_review |

**SCGB1D4** (`ENSG00000197745`): qRT-PCR with sequenced PCR products. Source: Results, Figures 1-3; expression-regulation experiments. Limitations: Cellular location and protein production are unresolved; a normal-expression concern, not RA support or a tumor nomination.

## [Kurji et al., 2010](https://doi.org/10.1167/iovs.09-3622)

Microarray Analysis Identifies Changes in Inflammatory Gene Expression in Response to Amyloid-β Stimulation of Cultured Human Retinal Pigment Epithelial Cells

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| PSG8 | PSG8 | The PSG8-labelled probe has a reported 1.64-fold decrease after 24 hours of amyloid-β treatment (P=0.0123). This is an additional ocular RNA observation requiring assay-specific attribution. | Cultured human retinal pigment epithelial cells, amyloid-β treatment and control | normal_tissue_review |

**PSG8** (`ENSG00000124467`): Agilent oligonucleotide microarray. Source: Table 2, PSG8 / NM_182707; microarray and qRT-PCR methods. Limitations: Probe specificity against other PSG genes has not been established in this review. PSG8 was not in the reported 17 upregulated genes selected for qRT-PCR validation. Cultured RPE is not evidence of expression in intact healthy retina or of PSG8 protein. This does not resolve the retinal single-cell mapping hold.

## [Zhang et al., 2020](https://doi.org/10.1186/s12885-020-07155-y)

Suppressing Dazl modulates tumorigenicity and stemness in human glioblastoma cells

| Gene | Original name | Reported as | Context | Allowed uses |
|---|---|---|---|---|
| DAZL | DAZL | The study reports higher DAZL RNA in glioblastoma and glioblastoma cell lines than the normal comparators. It supplies tumor-expression support, not a test of normal retinal astrocytes. | Human glioblastoma specimens and cell lines, with adjacent brain and commercial normal human astrocyte comparators | tumor_support, normal_tissue_review |

**DAZL** (`ENSG00000092345`): qRT-PCR; the study also reports immunoblotting and perturbation experiments. Source: Figure 1; cell culture and qRT-PCR methods. Limitations: Relative expression in commercial normal human astrocytes does not establish absence in all normal astrocytes. The normal comparator is not the ocular astrocyte population with the Census signal. This review does not validate antibody specificity or use the paper to dismiss that signal.

