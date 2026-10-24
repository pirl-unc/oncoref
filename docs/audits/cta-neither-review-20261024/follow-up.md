# Candidate follow-up record

Reviewed 24 October 2026 by Codex during the follow-up to PR #618. This is an implementation/evidence review, not an independent clinical adjudication. These records never override the general classification rules.

[Report](index.md) · [Full follow-up CSV](follow-up.csv) · [Measured high RNA](high-normal-rna.csv.gz) · [All context maxima](rna-context-maxima.csv.gz) · [Pinned source URLs and checksums](source-provenance.json)

## ACRBP

**Neither → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **5316.5 nTPM** and **56.543 TPM** (separate units). Reproductive bulk maximum: 275.6 nTPM.

**Open question:** Does the broad immune and platelet RNA represent endogenous ACRBP expression?

**Evidence needed:** Resolve normal-cell RNA identity and reliable local protein evidence before changing the general exclusion.

**Limits:** The retracted PMID 35463688 is excluded; the decision uses pinned RNA. Platelet RNA does not establish peptide presentation.

**High normal RNA findings:** dice_m2: blood/M2, max 39.1126 TPM, 12 observations >30; hpa_immune_cell: blood/non-classical monocyte, max 31.2 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/non-classical monocyte, max 56.5427 TPM, 1 observations >30; hpa_immune_cell_samples: blood/non-classical monocyte, max 77.8462 nTPM, 1 observations >30; hpa_immune_monaco: blood/intermediate monocyte, max 30.6 TPM, 1 observations >30; hpa_immune_monaco: blood/non-classical monocyte, max 37.5 TPM, 1 observations >30; hpa_single_cell_tissue: pbmc/platelets, max 5316.5 nTPM, 1 observations >30

**Immune RNA findings:** DICE B_CELL_NAIVE (unstimulated): 77/91 values ≥5 TPM, max 13.4406; DICE MONOCYTES (unstimulated): 91/91 values ≥5 TPM, max 20.75; DICE M2 (unstimulated): 90/90 values ≥5 TPM, max 39.1126; DICE NK (unstimulated): 35/90 values ≥5 TPM, max 8.60777; DICE TREG_MEM (unstimulated): 64/89 values ≥5 TPM, max 8.61494; DICE CD4_NAIVE (unstimulated): 72/88 values ≥5 TPM, max 11.1159; DICE CD4_STIM (activated): 4/88 values ≥5 TPM, max 7.53764; DICE TREG_NAIVE (unstimulated): 83/89 values ≥5 TPM, max 12.1914; DICE TFH (unstimulated): 24/89 values ≥5 TPM, max 8.47375; DICE TH1 (unstimulated): 22/81 values ≥5 TPM, max 8.18315; DICE TH2 (unstimulated): 18/89 values ≥5 TPM, max 7.37337; DICE CD8_NAIVE (unstimulated): 86/89 values ≥5 TPM, max 14.0219; DICE CD8_STIM (activated): 20/88 values ≥5 TPM, max 9.34729; hpa_immune_cell classical monocyte: 12.5 nTPM; hpa_immune_cell intermediate monocyte: 15.8 nTPM; hpa_immune_cell myeloid DC: 17.7 nTPM; hpa_immune_cell non-classical monocyte: 31.2 nTPM; hpa_immune_cell total PBMC: 10.6 nTPM; hpa_immune_monaco basophil: 5 TPM; hpa_immune_monaco classical monocyte: 17.6 TPM; hpa_immune_monaco Exhausted memory B-cell: 5.9 TPM; hpa_immune_monaco intermediate monocyte: 30.6 TPM; hpa_immune_monaco myeloid DC: 20.3 TPM; hpa_immune_monaco naive B-cell: 5.6 TPM; hpa_immune_monaco neutrophil: 5.2 TPM; hpa_immune_monaco non-classical monocyte: 37.5 TPM; hpa_immune_monaco Non-switched memory B-cell: 5 TPM; hpa_immune_monaco plasmacytoid DC: 10.2 TPM; hpa_immune_monaco Progenitor cell: 28.5 TPM; hpa_immune_monaco Switched memory B-cell: 5.1 TPM; hpa_immune_monaco total PBMC: 10.8 TPM

**Cell RNA findings:** pbmc/platelets cluster c-9: 5316.5 nTPM, 559 reads; reproductive-cell max 1326.2

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000111644-ACRBP) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## HORMAD1

**Neither → Neither / Cellular.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **14.8 nTPM** and **19.733 TPM** (separate units). Reproductive bulk maximum: 114.3 nTPM.

**Open question:** Why is RNA repeatedly detected in monocytes and naive B cells despite reproductive-positive curated IHC?

**Evidence needed:** Review gene-specific read mapping and obtain reliable protein evidence in the implicated immune populations.

**Limits:** DICE observations are not paired RNA/IHC donors; negative bulk-organ IHC does not assess every sorted immune population.

**Immune RNA findings:** DICE B_CELL_NAIVE (unstimulated): 33/91 values ≥5 TPM, max 11.3605; DICE MONOCYTES (unstimulated): 90/91 values ≥5 TPM, max 19.7327; DICE M2 (unstimulated): 48/90 values ≥5 TPM, max 14.6589; hpa_immune_monaco neutrophil: 5.7 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000143452-HORMAD1) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## TDRD1

**Neither → Neither / Cellular.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **10.8 nTPM** and **6.434 TPM** (separate units). Reproductive bulk maximum: 15.2 nTPM.

**Open question:** Are the few activated-CD4 values and Monaco plasmacytoid-DC summary reproducible?

**Evidence needed:** Resolve the borderline immune RNA with source-specific read review and reliable local protein evidence.

**Limits:** Three of 88 activated-CD4 DICE values reach 5 TPM; this is a review threshold rather than demonstrated toxicity.

**Immune RNA findings:** DICE CD4_STIM (activated): 3/88 values ≥5 TPM, max 6.43364; hpa_immune_monaco plasmacytoid DC: 6.4 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000095627-TDRD1) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## BRDT

**Neither → Neither / Cellular.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **17.7 nTPM** and **6 TPM** (separate units). Reproductive bulk maximum: 127.6 nTPM.

**Open question:** Does the Monaco-only borderline signal reflect genuine immune expression?

**Evidence needed:** Reconcile Monaco CD8/basophil RNA with the much lower DICE values and obtain relevant protein evidence.

**Limits:** Monaco summary values 5 and 6 TPM are weaker evidence than recurrent high original-sample RNA; DICE does not sample basophils.

**Immune RNA findings:** hpa_immune_monaco basophil: 5 TPM; hpa_immune_monaco Effector memory CD8 T-cell: 6 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000137948-BRDT) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## TEX101

**Neither → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **119.722 nTPM** and **74.9 TPM** (separate units). Reproductive bulk maximum: 114.4 nTPM.

**Open question:** What explains high basophil RNA in HPA and Monaco?

**Evidence needed:** Resolve basophil RNA identity and local protein evidence before reconsideration.

**Limits:** TPM and nTPM are separate measurements; DICE has no basophil population.

**High normal RNA findings:** hpa_immune_cell: blood/basophil, max 55.7 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/basophil, max 66.4757 TPM, 2 observations >30; hpa_immune_cell_samples: blood/basophil, max 119.722 nTPM, 2 observations >30; hpa_immune_monaco: blood/basophil, max 74.9 TPM, 1 observations >30

**Immune RNA findings:** hpa_immune_cell basophil: 55.7 nTPM; hpa_immune_monaco basophil: 74.9 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000131126-TEX101) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## DAZL

**Neither → Neither / Cellular.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **2.7 nTPM** and **8.913 TPM** (separate units). Reproductive bulk maximum: 77.3 nTPM.

**Open question:** Is the naive-B-cell RNA translated in normal B cells?

**Evidence needed:** Review gene-specific mapping and reliable protein evidence in naive B cells.

**Limits:** Reproductive-positive Enhanced IHC does not establish absence in the particular sorted immune cells measured by RNA.

**Immune RNA findings:** DICE B_CELL_NAIVE (unstimulated): 16/91 values ≥5 TPM, max 8.91336; hpa_immune_monaco naive B-cell: 7.5 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000092345-DAZL) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## FOXR2

**Neither → Neither / Antibody.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **0 nTPM** and **0.513 TPM** (separate units). Reproductive bulk maximum: 1.2 nTPM.

**Open question:** Which primary antibody staining is specific to FOXR2?

**Evidence needed:** Resolve discordance between HPA034487 and HPA057358 using independently validated protein evidence.

**Limits:** Both profiles are Uncertain; reproductive/tumor publication support does not validate normal-organ staining or waive counterevidence.

**Antibody findings:** Cerebral cortex/Endothelial cells: Low, HPA034487; conflicting or uncurated primary profile; Colon/Peripheral nerve/ganglion: Low, HPA034487; conflicting or uncurated primary profile; Lymph node/Non-germinal center cells: Low, HPA034487; conflicting or uncurated primary profile; Pancreas/Exocrine glandular cells: Medium, HPA034487; conflicting or uncurated primary profile; Skeletal muscle/Myocytes: Medium, HPA034487; conflicting or uncurated primary profile; Skin 1/Keratinocytes: Low, HPA034487; conflicting or uncurated primary profile; Skin 2/Epidermal cells: Low, HPA034487; conflicting or uncurated primary profile; Thyroid gland/Glandular cells: Low, HPA034487; conflicting or uncurated primary profile; Adipose tissue/Adipocytes: Medium, HPA057358; conflicting or uncurated primary profile; Adrenal gland/Glandular cells: Medium, HPA057358; conflicting or uncurated primary profile; Bone marrow/Hematopoietic cells: Low, HPA057358; conflicting or uncurated primary profile; Bronchus/Respiratory epithelial cells: Low, HPA057358; conflicting or uncurated primary profile; Caudate/Glial cells: Medium, HPA057358; conflicting or uncurated primary profile; Cerebellum/Purkinje cells: Medium, HPA057358; conflicting or uncurated primary profile; Cerebral cortex/Endothelial cells: Low, HPA057358; conflicting or uncurated primary profile; Cerebral cortex/Neuronal cells: Medium, HPA057358; conflicting or uncurated primary profile; Cerebral cortex/Neuropil: Low, HPA057358; conflicting or uncurated primary profile; Colon/Endothelial cells: Medium, HPA057358; conflicting or uncurated primary profile; Heart muscle/Cardiomyocytes: Medium, HPA057358; conflicting or uncurated primary profile; Lung/Alveolar cells: Low, HPA057358; conflicting or uncurated primary profile; Lung/Macrophages: High, HPA057358; conflicting or uncurated primary profile; Nasopharynx/Respiratory epithelial cells: Low, HPA057358; conflicting or uncurated primary profile; Skin 1/Fibroblasts: Medium, HPA057358; conflicting or uncurated primary profile; Skin 2/Epidermal cells: Low, HPA057358; conflicting or uncurated primary profile; Stomach 2/Glandular cells: Low, HPA057358; conflicting or uncurated primary profile; Thyroid gland/Glandular cells: Medium, HPA057358; conflicting or uncurated primary profile

**Publication sources:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9437574/

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000189299-FOXR2) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## CTAG2

**Neither → Neither / Variable.** Donor-variable vital RNA with negative local IHC; excluded from targets and forbidden background pending review

Maximum outside explicitly permitted contexts: **24.9 nTPM** and **5.075 TPM** (separate units). Reproductive bulk maximum: 10.1 nTPM.

**Open question:** When and why do some normal heart specimens contain CTAG2 RNA?

**Evidence needed:** Resolve donor variability and obtain paired RNA/protein information from the implicated cardiac context.

**Limits:** RNA donors and negative-IHC donors differ; neither absence of established target-mediated toxicity nor separate negative specimens clears the positive RNA donors.

**Donor findings:** 149: 13.5 nTPM (fails); 150: 0 nTPM (passes); 151: 0 nTPM (passes); 152: 7.1 nTPM (fails)

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000126890-CTAG2) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## GAGE12F

**Neither → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **83.07 nTPM** and **71.665 TPM** (separate units). Reproductive bulk maximum: 103.3 nTPM.

**Open question:** Is the immune signal correctly assigned within the GAGE family?

**Evidence needed:** Review uniquely assigning reads and exact gene/protein identity before reconsidering local expression.

**Limits:** Highly related paralogs can complicate attribution; no paralog substitution or RNA-to-protein inference is permitted.

**High normal RNA findings:** hpa_immune_cell_samples: blood/NK-cell, max 51.8037 TPM, 1 observations >30; hpa_immune_cell_samples: blood/T-reg, max 71.6653 TPM, 1 observations >30; hpa_immune_cell_samples: blood/classical monocyte, max 38.8727 TPM, 1 observations >30; hpa_immune_cell_samples: blood/intermediate monocyte, max 32.0759 TPM, 1 observations >30; hpa_immune_cell_samples: blood/memory CD4 T-cell, max 46.1675 TPM, 1 observations >30; hpa_immune_cell_samples: blood/naive CD4 T-cell, max 35.4421 TPM, 1 observations >30; hpa_immune_cell_samples: blood/naive CD8 T-cell, max 33.0234 TPM, 1 observations >30; hpa_immune_cell_samples: blood/plasmacytoid DC, max 41.922 TPM, 1 observations >30; hpa_immune_cell_samples: blood/total PBMC, max 40.9921 TPM, 1 observations >30; hpa_immune_cell_samples: blood/NK-cell, max 53.5002 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/T-reg, max 83.0695 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/classical monocyte, max 46.6655 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/intermediate monocyte, max 47.1668 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/memory B-cell, max 33.8287 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/memory CD4 T-cell, max 49.0849 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/naive B-cell, max 32.1517 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/naive CD4 T-cell, max 37.4272 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/naive CD8 T-cell, max 32.3725 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/non-classical monocyte, max 31.2073 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/plasmacytoid DC, max 50.2688 nTPM, 1 observations >30; hpa_immune_cell_samples: blood/total PBMC, max 46.4602 nTPM, 1 observations >30

**Immune RNA findings:** hpa_immune_cell classical monocyte: 7.7 nTPM; hpa_immune_cell intermediate monocyte: 7.8 nTPM; hpa_immune_cell memory B-cell: 5.6 nTPM; hpa_immune_cell memory CD4 T-cell: 8.1 nTPM; hpa_immune_cell memory CD8 T-cell: 5.5 nTPM; hpa_immune_cell naive B-cell: 5.3 nTPM; hpa_immune_cell naive CD4 T-cell: 6.2 nTPM; hpa_immune_cell naive CD8 T-cell: 5.3 nTPM; hpa_immune_cell NK-cell: 8.9 nTPM; hpa_immune_cell non-classical monocyte: 6.2 nTPM; hpa_immune_cell plasmacytoid DC: 8.3 nTPM; hpa_immune_cell T-reg: 13.8 nTPM; hpa_immune_cell total PBMC: 7.7 nTPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000236362-GAGE12F) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## PAGE2B

**Neither → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **35 nTPM** and **8.043 TPM** (separate units). Reproductive bulk maximum: 55.0 nTPM.

**Open question:** Are borderline DICE B-cell RNA and other cell-context signals genuine?

**Evidence needed:** Review exact loci and source-specific normal-cell signals; reliable protein evidence must address implicated cells.

**Limits:** A borderline immune summary must not hide a higher signal elsewhere; complete RNA audit accompanies this record.

**High normal RNA findings:** hpa_single_cell_tissue: liver/erythroid cells, max 35 nTPM, 1 observations >30

**Immune RNA findings:** DICE B_CELL_NAIVE (unstimulated): 3/91 values ≥5 TPM, max 8.04323

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000238269-PAGE2B) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## TRIP13

**Neither → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **88.8 nTPM** and **8.639 TPM** (separate units). Reproductive bulk maximum: 29.0 nTPM.

**Open question:** Does its broader normal-cell expression permit a reproductive-restricted classification?

**Evidence needed:** Resolve the documented non-permitted cell expression before any reproductive-restricted designation.

**Limits:** Tumor relevance and reproductive enrichment do not establish normal-cell restriction.

**High normal RNA findings:** hpa_single_cell: /Ciliated cells, max 46.8 nTPM, 1 observations >30; hpa_single_cell_tissue: adipose tissue/t-cells, max 33.6 nTPM, 1 observations >30; hpa_single_cell_tissue: bone marrow/plasma cells, max 37.8 nTPM, 1 observations >30; hpa_single_cell_tissue: bronchus/ciliated cells, max 51.6 nTPM, 1 observations >30; hpa_single_cell_tissue: liver/erythroid cells, max 35 nTPM, 1 observations >30; hpa_single_cell_tissue: lung/ciliated cells, max 88.8 nTPM, 1 observations >30; hpa_single_cell_tissue: pancreas/fibroblasts, max 30.1 nTPM, 1 observations >30; hpa_single_cell_tissue: placenta/mixed immune cells, max 60.2 nTPM, 1 observations >30; hpa_single_cell_tissue: spleen/plasma cells, max 45.7 nTPM, 1 observations >30; hpa_single_cell_tissue: thymus/plasma cells, max 48.4 nTPM, 2 observations >30

**Immune RNA findings:** DICE CD4_STIM (activated): 8/88 values ≥5 TPM, max 7.21408; DICE CD8_STIM (activated): 41/88 values ≥5 TPM, max 8.63936; hpa_immune_monaco Plasmablast: 7.4 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000071539-TRIP13) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## CTAGE1

**Neither → Neither / Antibody.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **1.9 nTPM** and **0.388 TPM** (separate units). Reproductive bulk maximum: 10.0 nTPM.

**Open question:** Can the uncertain primary-antibody normal staining be resolved?

**Evidence needed:** Use independently validated protein evidence to adjudicate primary staining conflicts.

**Limits:** An uncertain antibody is not reliable proof of either expression or absence.

**Antibody findings:** Adrenal gland/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Bone marrow/Hematopoietic cells: Medium, HPA000387; conflicting or uncurated primary profile; Caudate/Glial cells: Low, HPA000387; conflicting or uncurated primary profile; Caudate/Neuronal cells: High, HPA000387; conflicting or uncurated primary profile; Cerebellum/Cells in granular layer: Medium, HPA000387; conflicting or uncurated primary profile; Cerebellum/Cells in molecular layer: Medium, HPA000387; conflicting or uncurated primary profile; Cerebellum/Purkinje cells: Medium, HPA000387; conflicting or uncurated primary profile; Cerebral cortex/Endothelial cells: Low, HPA000387; conflicting or uncurated primary profile; Cerebral cortex/Glial cells: Low, HPA000387; conflicting or uncurated primary profile; Cerebral cortex/Neuronal cells: High, HPA000387; conflicting or uncurated primary profile; Colon/Endothelial cells: Low, HPA000387; conflicting or uncurated primary profile; Colon/Glandular cells: Medium, HPA000387; conflicting or uncurated primary profile; Duodenum/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Esophagus/Squamous epithelial cells: Low, HPA000387; conflicting or uncurated primary profile; Hippocampus/Glial cells: Low, HPA000387; conflicting or uncurated primary profile; Hippocampus/Neuronal cells: High, HPA000387; conflicting or uncurated primary profile; Kidney/Cells in tubules: Medium, HPA000387; conflicting or uncurated primary profile; Liver/Cholangiocytes: Low, HPA000387; conflicting or uncurated primary profile; Liver/Hepatocytes: Medium, HPA000387; conflicting or uncurated primary profile; Lung/Alveolar cells: Medium, HPA000387; conflicting or uncurated primary profile; Lung/Macrophages: High, HPA000387; conflicting or uncurated primary profile; Lymph node/Germinal center cells: Medium, HPA000387; conflicting or uncurated primary profile; Nasopharynx/Respiratory epithelial cells: High, HPA000387; conflicting or uncurated primary profile; Oral mucosa/Squamous epithelial cells: Medium, HPA000387; conflicting or uncurated primary profile; Pancreas/Exocrine glandular cells: High, HPA000387; conflicting or uncurated primary profile; Parathyroid gland/Glandular cells: Low, HPA000387; conflicting or uncurated primary profile; Rectum/Glandular cells: Low, HPA000387; conflicting or uncurated primary profile; Salivary gland/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Skin 1/Keratinocytes: Low, HPA000387; conflicting or uncurated primary profile; Skin 2/Epidermal cells: Medium, HPA000387; conflicting or uncurated primary profile; Small intestine/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Smooth muscle/Smooth muscle cells: Medium, HPA000387; conflicting or uncurated primary profile; Spleen/Cells in red pulp: Low, HPA000387; conflicting or uncurated primary profile; Spleen/Cells in white pulp: Low, HPA000387; conflicting or uncurated primary profile; Stomach 1/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Stomach 2/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Thyroid gland/Glandular cells: High, HPA000387; conflicting or uncurated primary profile; Urinary bladder/Urothelial cells: Medium, HPA000387; conflicting or uncurated primary profile; Adrenal gland/Glandular cells: Low, HPA000922; conflicting or uncurated primary profile; Bronchus/Respiratory epithelial cells: Medium, HPA000922; conflicting or uncurated primary profile; Caudate/Neuronal cells: Low, HPA000922; conflicting or uncurated primary profile; Cerebellum/Purkinje cells: Low, HPA000922; conflicting or uncurated primary profile; Cerebral cortex/Neuronal cells: Low, HPA000922; conflicting or uncurated primary profile; Colon/Glandular cells: Low, HPA000922; conflicting or uncurated primary profile; Duodenum/Glandular cells: Medium, HPA000922; conflicting or uncurated primary profile; Hippocampus/Neuronal cells: High, HPA000922; conflicting or uncurated primary profile; Kidney/Cells in tubules: Low, HPA000922; conflicting or uncurated primary profile; Liver/Hepatocytes: Low, HPA000922; conflicting or uncurated primary profile; Lung/Macrophages: Medium, HPA000922; conflicting or uncurated primary profile; Lymph node/Germinal center cells: Low, HPA000922; conflicting or uncurated primary profile; Nasopharynx/Respiratory epithelial cells: Medium, HPA000922; conflicting or uncurated primary profile; Pancreas/Exocrine glandular cells: Medium, HPA000922; conflicting or uncurated primary profile; Pancreas/Pancreatic endocrine cells: Medium, HPA000922; conflicting or uncurated primary profile; Rectum/Glandular cells: Low, HPA000922; conflicting or uncurated primary profile; Salivary gland/Glandular cells: Medium, HPA000922; conflicting or uncurated primary profile; Skin 2/Epidermal cells: Low, HPA000922; conflicting or uncurated primary profile; Small intestine/Glandular cells: Low, HPA000922; conflicting or uncurated primary profile; Stomach 1/Glandular cells: Medium, HPA000922; conflicting or uncurated primary profile; Stomach 2/Glandular cells: Medium, HPA000922; conflicting or uncurated primary profile; Thyroid gland/Glandular cells: Low, HPA000922; conflicting or uncurated primary profile; Urinary bladder/Urothelial cells: Low, HPA000922; conflicting or uncurated primary profile; Adrenal gland/Glandular cells: Low, CAB016204; conflicting or uncurated primary profile; Bone marrow/Hematopoietic cells: High, CAB016204; conflicting or uncurated primary profile; Bronchus/Respiratory epithelial cells: High, CAB016204; conflicting or uncurated primary profile; Cerebellum/Cells in molecular layer: Low, CAB016204; conflicting or uncurated primary profile; Colon/Glandular cells: Low, CAB016204; conflicting or uncurated primary profile; Esophagus/Squamous epithelial cells: High, CAB016204; conflicting or uncurated primary profile; Hippocampus/Glial cells: Low, CAB016204; conflicting or uncurated primary profile; Kidney/Cells in tubules: High, CAB016204; conflicting or uncurated primary profile; Liver/Cholangiocytes: Medium, CAB016204; conflicting or uncurated primary profile; Liver/Hepatocytes: High, CAB016204; conflicting or uncurated primary profile; Lung/Alveolar cells: Medium, CAB016204; conflicting or uncurated primary profile; Lung/Macrophages: Low, CAB016204; conflicting or uncurated primary profile; Lymph node/Non-germinal center cells: Medium, CAB016204; conflicting or uncurated primary profile; Nasopharynx/Respiratory epithelial cells: High, CAB016204; conflicting or uncurated primary profile; Oral mucosa/Squamous epithelial cells: High, CAB016204; conflicting or uncurated primary profile; Pancreas/Exocrine glandular cells: Medium, CAB016204; conflicting or uncurated primary profile; Rectum/Glandular cells: Medium, CAB016204; conflicting or uncurated primary profile; Salivary gland/Glandular cells: High, CAB016204; conflicting or uncurated primary profile; Skin 1/Langerhans: High, CAB016204; conflicting or uncurated primary profile; Spleen/Cells in red pulp: High, CAB016204; conflicting or uncurated primary profile; Stomach 1/Glandular cells: Low, CAB016204; conflicting or uncurated primary profile; Stomach 2/Glandular cells: Medium, CAB016204; conflicting or uncurated primary profile; Urinary bladder/Urothelial cells: High, CAB016204; conflicting or uncurated primary profile

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000212710-CTAGE1) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## ZNF730

**Neither → Neither / Cellular.** Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins

Maximum outside explicitly permitted contexts: **11.3 nTPM** and **7 TPM** (separate units). Reproductive bulk maximum: 4.9 nTPM.

**Open question:** Does the full normal-cell distribution support reproductive restriction?

**Evidence needed:** Review diffuse RNA and missing reliable protein evidence rather than selecting the largest reproductive value.

**Limits:** A high reproductive maximum does not erase normal-cell expression; no clinical safety conclusion is recorded.

**Immune RNA findings:** hpa_immune_monaco Progenitor cell: 7 TPM

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000183850-ZNF730) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## MAGEA3

**CTA → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **53.1 nTPM** and **0.166 TPM** (separate units). Reproductive bulk maximum: 7.8 nTPM.

**Open question:** Is the colonic undifferentiated-cell RNA correctly attributed and endogenous?

**Evidence needed:** Resolve the pinned high normal-cell signal and exact locus/protein identity before reconsideration.

**Limits:** Clinical interest does not create a gene exception; RNA does not by itself establish protein or peptide display.

**High normal RNA findings:** hpa_single_cell_tissue: colon/undifferentiated cells, max 53.1 nTPM, 1 observations >30

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000221867-MAGEA3) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)

## ROPN1

**CTA → Unsafe / High RNA.** Observed RNA above 30 TPM or nTPM outside explicitly permitted reproductive/disposable contexts; research exclusion, not proven clinical toxicity

Maximum outside explicitly permitted contexts: **32.9 nTPM** and **0.595 TPM** (separate units). Reproductive bulk maximum: 154.4 nTPM.

**Open question:** Is the salivary serous-cell RNA above 30 nTPM a genuine local signal?

**Evidence needed:** Review the underlying cluster and gene-specific reads with reliable local protein evidence before reconsideration.

**Limits:** The value is close to the hard ceiling; a general rule is applied rather than a benchmark exception.

**High normal RNA findings:** hpa_single_cell_tissue: salivary gland/serous glandular cells, max 32.9 nTPM, 1 observations >30

[HPA v23 gene record](https://v23.proteinatlas.org/ENSG00000065371-ROPN1) · [Previous decision and measured antibodies](../cta-normal-cells-20261023/explore.html)
