# Completed CTA adjudications — 25 October 2026

**All 178 antibody cases remaining in the local 391-gene draft have a completed decision: 177 Neither and 1 Unsafe.** No member of this discrepancy queue remains CTA. The released list has **213 CTA genes**.

[Current rules and funnel](../cta-policy-20261025/index.md) · [All case records](adjudications.csv) · [Raw positive observations](primary-observations.csv.gz) · [Reviewed findings](../cta-policy-20261025/reviewed-findings.csv)

## Decision boundary

- Accepted endogenous protein in non-permitted normal cells excludes unrestricted gene-level use: **Unsafe**.
- A material raw normal stain whose molecular specificity is unresolved gives **Neither**. Low intensity is not, by itself, proof of nonspecific staining. A generic off-target annotation does not identify which specimens are false positives.
- Curated negative IHC and passing RNA remain useful evidence, but do not silently erase a conflicting raw stain. Enhanced profile reliability does not validate every individual primary stain.
- A clinical product recognizes a particular peptide/HLA/receptor combination. Product-specific efficacy or tolerability does not clear every protein or sequence encoded by its target gene.

These are conservative research classifications. Neither is not proven toxicity; Unsafe is not a claim that an existing therapy caused injury. Future changes require source-specific evidence revisions, not an allowlist of familiar genes.

## What was actually reviewed

The 178-case review covers the pinned mapped positive primary-antibody rows, their named normal cells, intensity, primary reagents, HPA profile notes and staining summaries, with the RNA context preserved. This is a review of annotations and reagent evidence, **not an independent histopathology reread of every image** or new experimental validation. All 178 decisions are complete; molecular specificity remains unresolved in the raw-antibody findings. PRAME has additional accepted independent normal-protein evidence.

The eight initial cases were rechecked against their primary studies at the access depth recorded in the findings table. Full original methods were available for FTMT, YBX2, ACTRT1 and PLCZ1; other records explicitly retain abstract/excerpt and reagent-access limitations. A search result is not a full-paper adjudication. The broader 2,540-gene literature review and a continuous breadth metric remain open work.

## Eight initial cases and the important clinical targets

| Gene | Decision | Adjudication and limits |
|---|---|---|
| FTMT | Unsafe | Normal human nigral neurons with antibody specificity controls; small postmortem series. [Study](https://pubmed.ncbi.nlm.nih.gov/28386150/) |
| ABCC11 | Unsafe | Reported normal neural protein with multimodal support. Full reagent methods remain an access limitation. [Study](https://pubmed.ncbi.nlm.nih.gov/16359813/) |
| YBX2 | Unsafe | Vascular smooth-muscle/adrenal protein; not evidence of cardiomyocyte expression. Polyclonal and specimen limitations remain. [Study](https://pubmed.ncbi.nlm.nih.gov/16479255/) |
| CCDC158 | Unsafe | Reported normal renal tubular protein/function. Renal versus testicular isoforms remain unresolved; no isoform-specific clearance is claimed. [Study](https://pubmed.ncbi.nlm.nih.gov/39319391/) |
| ACTRT1 | Unsafe | Endogenous protein in differentiating normal human keratinocytes. Cultured-cell state is explicit; an RPE cell line is not normal retinal tissue evidence. [Study](https://pubmed.ncbi.nlm.nih.gov/33972689/) |
| MAGEC2 | Neither | Purkinje staining with TC-1 requires gene/family specificity resolution. [Study](https://pubmed.ncbi.nlm.nih.gov/12920247/) |
| LGALS13 | Neither | Historical extra-placental protein reactivity versus placental RNA restriction; galectin identity and anatomy are unresolved. [Early study](https://pubmed.ncbi.nlm.nih.gov/10527825/), [placental study](https://pubmed.ncbi.nlm.nih.gov/19497882/) |
| PLCZ1 | Neither | **Assay correction:** NYD-SP27 was measured by RT-PCR; the paper's pancreatic immunostaining measured CFTR. Reject the old PLCZ1 protein interpretation. Qualitative pancreatic RNA and isoform identity still require quantitative reconciliation. [Original paper](https://www.cuhk.edu.hk/ecbrc/cliparts/Publication%20list/Lancet2003.pdf) |
| PRAME | Unsafe | Independent monoclonal work reports low normal adrenal-cortex protein. A later negative panel using other reagents does not negate it. The initial positive study was reviewed at abstract level; clinical injury is not established. [Positive study](https://pubmed.ncbi.nlm.nih.gov/33890816/), [negative panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC9588667/) |
| MAGEA4 | Neither | HPA normal neuronal/epidermal stains remain unresolved. An independent OTI1F9 assay has reassuring normal-panel data but does not directly resolve the older HPA reagents. Tumor diagnostic assay performance is not whole-gene safety. [Validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11170170/) |

**Preserved therapeutic evidence:** [PRAME/IMA203](https://pubmed.ncbi.nlm.nih.gov/40205198/) and [MAGE-A4/afami-cel](https://pubmed.ncbi.nlm.nih.gov/38554725/) remain product-specific clinical evidence. Neither publication is discarded, and neither supplies an unrestricted gene-level rescue. CTAG1A/B and XAGE1A/B remain CTA. CTAG2 remains Neither/Variable, with the donor evidence preserved.

## Every antibody case

Each heading is searchable by gene. The CSV carries exact IDs, source-native RNA maxima, reagent summaries and all contexts. “Evidence resolved = false” is a completed hold decision, not an unfinished review.

### ACR — Neither

`ENSG00000100312` · 5 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA048687.

**Contexts:** Bronchus / Respiratory epithelial cells; Colon / Glandular cells; Kidney / Cells in tubules; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000100312-ACR/tissue/primary+data). Review records: ihc-20261025-ENSG00000100312.

### ACSBG2 — Neither

`ENSG00000130377` · 20 positive raw rows; 11 Medium/High; 6 Vital Medium/High. Reagents: HPA043421; HPA051013.

**Contexts:** Bronchus / Respiratory epithelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Nasopharynx / Respiratory epithelial cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 11 Medium/High normal annotations (6 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000130377-ACSBG2/tissue/primary+data). Review records: ihc-20261025-ENSG00000130377.

### ACTL7B — Neither

`ENSG00000148156` · 4 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA021803; HPA024664.

**Contexts:** Bronchus / Respiratory epithelial cells; Cerebral cortex / Neuronal cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000148156-ACTL7B/tissue/primary+data). Review records: ihc-20261025-ENSG00000148156.

### ADAM2 — Neither

`ENSG00000104755` · 28 positive raw rows; 16 Medium/High; 11 Vital Medium/High. Reagents: CAB022390; HPA024621; HPA026581.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Neuronal cells; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Lung / Macrophages; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** 16 Medium/High normal annotations (11 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000104755-ADAM2/tissue/primary+data). Review records: ihc-20261025-ENSG00000104755.

### AKAP4 — Neither

`ENSG00000147081` · 14 positive raw rows; 3 Medium/High; 2 Vital Medium/High. Reagents: HPA005949; HPA020046.

**Contexts:** Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Cells in molecular layer; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Hippocampus / Glial cells; Hippocampus / Neuronal cells; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** 3 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000147081-AKAP4/tissue/primary+data). Review records: ihc-20261025-ENSG00000147081.

### ASZ1 — Neither

`ENSG00000154438` · 8 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA020558.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Kidney / Cells in tubules; Rectum / Glandular cells; Salivary gland / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000154438-ASZ1/tissue/primary+data). Review records: ihc-20261025-ENSG00000154438.

### BOLL — Neither

`ENSG00000152430` · 12 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA018678; HPA048813.

**Contexts:** Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Purkinje cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Kidney / Cells in tubules; Lung / Macrophages; Lymph node / Non-germinal center cells; Spleen / Cells in red pulp; Thyroid gland / Glandular cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000152430-BOLL/tissue/primary+data). Review records: ihc-20261025-ENSG00000152430.

### C10orf53 — Neither

`ENSG00000178645` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA037951.

**Contexts:** Kidney / Cells in tubules; Salivary gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000178645-C10orf53/tissue/primary+data). Review records: ihc-20261025-ENSG00000178645.

### C10orf62 — Neither

`ENSG00000203942` · 6 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA072947.

**Contexts:** Caudate / Neuronal cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Hippocampus / Neuronal cells; Nasopharynx / Respiratory epithelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000203942-C10orf62/tissue/primary+data). Review records: ihc-20261025-ENSG00000203942.

### C11orf42 — Neither

`ENSG00000180878` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA063404.

**Contexts:** Adrenal gland / Glandular cells; Smooth muscle / Smooth muscle cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000180878-C11orf42/tissue/primary+data). Review records: ihc-20261025-ENSG00000180878.

### C20orf144 — Neither

`ENSG00000149609` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA019029.

**Contexts:** Salivary gland / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000149609-C20orf144/tissue/primary+data). Review records: ihc-20261025-ENSG00000149609.

### CABS1 — Neither

`ENSG00000145309` · 3 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA055352.

**Contexts:** Duodenum / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000145309-CABS1/tissue/primary+data). Review records: ihc-20261025-ENSG00000145309.

### CALR3 — Neither

`ENSG00000269058` · 4 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: CAB025749.

**Contexts:** Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Lymph node / Non-germinal center cells; Spleen / Cells in red pulp.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000269058-CALR3/tissue/primary+data). Review records: ihc-20261025-ENSG00000269058.

### CAPZA3 — Neither

`ENSG00000177938` · 4 positive raw rows; 2 Medium/High; 2 Vital Medium/High. Reagents: HPA038688.

**Contexts:** Colon / Glandular cells; Duodenum / Glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000177938-CAPZA3/tissue/primary+data). Review records: ihc-20261025-ENSG00000177938.

### CATSPER4 — Neither

`ENSG00000188782` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA045257.

**Contexts:** Bronchus / Respiratory epithelial cells; Small intestine / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000188782-CATSPER4/tissue/primary+data). Review records: ihc-20261025-ENSG00000188782.

### CBY2 — Neither

`ENSG00000174015` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA039359.

**Contexts:** Skin 1 / Keratinocytes; Stomach 2 / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000174015-CBY2/tissue/primary+data). Review records: ihc-20261025-ENSG00000174015.

### CETN1 — Neither

`ENSG00000177143` · 41 positive raw rows; 23 Medium/High; 13 Vital Medium/High. Reagents: CAB037244; HPA028956; HPA069157.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebellum / Cells in granular layer; Cerebellum / Cells in molecular layer; Cerebellum / Purkinje cells; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Endothelial cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Glial cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Liver / Cholangiocytes; Lung / Alveolar cells; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Fibroblasts; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 23 Medium/High normal annotations (13 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000177143-CETN1/tissue/primary+data). Review records: ihc-20261025-ENSG00000177143.

### CFAP97D1 — Neither

`ENSG00000231256` · 17 positive raw rows; 6 Medium/High; 3 Vital Medium/High. Reagents: HPA053028.

**Contexts:** Bronchus / Respiratory epithelial cells; Cerebral cortex / Neuropil; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 6 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000231256-CFAP97D1/tissue/primary+data). Review records: ihc-20261025-ENSG00000231256.

### CGB3 — Neither

`ENSG00000104827` · 22 positive raw rows; 14 Medium/High; 8 Vital Medium/High. Reagents: CAB010884; HPA038925.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Alveolar cells; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skin 1 / Keratinocytes; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 14 Medium/High normal annotations (8 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000104827-CGB3/tissue/primary+data). Review records: ihc-20261025-ENSG00000104827.

### CHCT1 — Neither

`ENSG00000141371` · 4 positive raw rows; 2 Medium/High; 1 Vital Medium/High. Reagents: HPA044415.

**Contexts:** Pancreas / Exocrine glandular cells; Salivary gland / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000141371-CHCT1/tissue/primary+data). Review records: ihc-20261025-ENSG00000141371.

### CIMAP1A — Neither

`ENSG00000177947` · 15 positive raw rows; 3 Medium/High; 3 Vital Medium/High. Reagents: HPA038919; HPA039241.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Hippocampus / Glial cells; Liver / Hepatocytes; Oral mucosa / Squamous epithelial cells; Salivary gland / Glandular cells; Skin 2 / Epidermal cells; Smooth muscle / Smooth muscle cells; Stomach 2 / Glandular cells.

**Decision:** 3 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000177947-CIMAP1A/tissue/primary+data). Review records: ihc-20261025-ENSG00000177947.

### CIMIP2A — Neither

`ENSG00000188163` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA023307.

**Contexts:** Parathyroid gland / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000188163-CIMIP2A/tissue/primary+data). Review records: ihc-20261025-ENSG00000188163.

### CIMIP4 — Neither

`ENSG00000185264` · 15 positive raw rows; 4 Medium/High; 1 Vital Medium/High. Reagents: HPA005998.

**Contexts:** Adrenal gland / Glandular cells; Cerebellum / Purkinje cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Lung / Alveolar cells; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** 4 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000185264-CIMIP4/tissue/primary+data). Review records: ihc-20261025-ENSG00000185264.

### CT45A10 — Neither

`ENSG00000269586` · 38 positive raw rows; 12 Medium/High; 10 Vital Medium/High. Reagents: HPA044735; HPA044757; HPA046872; HPA046987.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 12 Medium/High normal annotations (10 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000269586-CT45A10/tissue/primary+data). Review records: ihc-20261025-ENSG00000269586.

### CT45A2 — Neither

`ENSG00000271449` · 38 positive raw rows; 12 Medium/High; 10 Vital Medium/High. Reagents: HPA044735; HPA044757; HPA046872; HPA046987.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 12 Medium/High normal annotations (10 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000271449-CT45A2/tissue/primary+data). Review records: ihc-20261025-ENSG00000271449.

### CT45A7 — Neither

`ENSG00000273696` · 38 positive raw rows; 12 Medium/High; 10 Vital Medium/High. Reagents: HPA044735; HPA044757; HPA046872; HPA046987.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 12 Medium/High normal annotations (10 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000273696-CT45A7/tissue/primary+data). Review records: ihc-20261025-ENSG00000273696.

### CT47A1 — Neither

`ENSG00000236371` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000236371-CT47A1/tissue/primary+data). Review records: ihc-20261025-ENSG00000236371.

### CT47A10 — Neither

`ENSG00000224089` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000224089-CT47A10/tissue/primary+data). Review records: ihc-20261025-ENSG00000224089.

### CT47A11 — Neither

`ENSG00000226929` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000226929-CT47A11/tissue/primary+data). Review records: ihc-20261025-ENSG00000226929.

### CT47A12 — Neither

`ENSG00000226685` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000226685-CT47A12/tissue/primary+data). Review records: ihc-20261025-ENSG00000226685.

### CT47A2 — Neither

`ENSG00000242362` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000242362-CT47A2/tissue/primary+data). Review records: ihc-20261025-ENSG00000242362.

### CT47A3 — Neither

`ENSG00000236126` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000236126-CT47A3/tissue/primary+data). Review records: ihc-20261025-ENSG00000236126.

### CT47A4 — Neither

`ENSG00000230594` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000230594-CT47A4/tissue/primary+data). Review records: ihc-20261025-ENSG00000230594.

### CT47A5 — Neither

`ENSG00000237957` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000237957-CT47A5/tissue/primary+data). Review records: ihc-20261025-ENSG00000237957.

### CT47A6 — Neither

`ENSG00000226023` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000226023-CT47A6/tissue/primary+data). Review records: ihc-20261025-ENSG00000226023.

### CT47A7 — Neither

`ENSG00000228517` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000228517-CT47A7/tissue/primary+data). Review records: ihc-20261025-ENSG00000228517.

### CT47A8 — Neither

`ENSG00000230347` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000230347-CT47A8/tissue/primary+data). Review records: ihc-20261025-ENSG00000230347.

### CT47A9 — Neither

`ENSG00000226600` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000226600-CT47A9/tissue/primary+data). Review records: ihc-20261025-ENSG00000226600.

### CT47B1 — Neither

`ENSG00000236446` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042483.

**Contexts:** Bone marrow / Hematopoietic cells; Lymph node / Non-germinal center cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000236446-CT47B1/tissue/primary+data). Review records: ihc-20261025-ENSG00000236446.

### CTCFL — Neither

`ENSG00000124092` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA001472.

**Contexts:** Rectum / Glandular cells; Small intestine / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000124092-CTCFL/tissue/primary+data). Review records: ihc-20261025-ENSG00000124092.

### CYLC1 — Neither

`ENSG00000183035` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA070355.

**Contexts:** Parathyroid gland / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000183035-CYLC1/tissue/primary+data). Review records: ihc-20261025-ENSG00000183035.

### DAZ3 — Neither

`ENSG00000187191` · 5 positive raw rows; 3 Medium/High; 3 Vital Medium/High. Reagents: HPA047113.

**Contexts:** Colon / Glandular cells; Duodenum / Glandular cells; Kidney / Cells in tubules; Rectum / Glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 3 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000187191-DAZ3/tissue/primary+data). Review records: ihc-20261025-ENSG00000187191.

### DDI1 — Neither

`ENSG00000170967` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA079511.

**Contexts:** Lymph node / Non-germinal center cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000170967-DDI1/tissue/primary+data). Review records: ihc-20261025-ENSG00000170967.

### DDX4 — Neither

`ENSG00000152670` · 34 positive raw rows; 15 Medium/High; 9 Vital Medium/High. Reagents: CAB026170; HPA037763; HPA037764.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebellum / Cells in molecular layer; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 15 Medium/High normal annotations (9 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000152670-DDX4/tissue/primary+data). Review records: ihc-20261025-ENSG00000152670.

### DKKL1 — Neither

`ENSG00000104901` · 7 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: CAB025410.

**Contexts:** Caudate / Neuronal cells; Kidney / Cells in glomeruli; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Smooth muscle / Smooth muscle cells; Urinary bladder / Urothelial cells.

**Decision:** 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000104901-DKKL1/tissue/primary+data). Review records: ihc-20261025-ENSG00000104901.

### DMRT1 — Neither

`ENSG00000137090` · 6 positive raw rows; 2 Medium/High; 1 Vital Medium/High. Reagents: HPA027850.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000137090-DMRT1/tissue/primary+data). Review records: ihc-20261025-ENSG00000137090.

### DMRTB1 — Neither

`ENSG00000143006` · 15 positive raw rows; 4 Medium/High; 2 Vital Medium/High. Reagents: HPA036490; HPA058553.

**Contexts:** Adipose tissue / Adipocytes; Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Heart muscle / Cardiomyocytes; Kidney / Cells in glomeruli; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Skin 1 / Fibroblasts; Skin 1 / Langerhans; Skin 2 / Epidermal cells.

**Decision:** 4 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000143006-DMRTB1/tissue/primary+data). Review records: ihc-20261025-ENSG00000143006.

### DNAJC5G — Neither

`ENSG00000163793` · 8 positive raw rows; 2 Medium/High; 1 Vital Medium/High. Reagents: HPA041445.

**Contexts:** Bone marrow / Hematopoietic cells; Esophagus / Squamous epithelial cells; Oral mucosa / Squamous epithelial cells; Salivary gland / Glandular cells; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Spleen / Cells in red pulp; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000163793-DNAJC5G/tissue/primary+data). Review records: ihc-20261025-ENSG00000163793.

### EDDM3B — Neither

`ENSG00000181552` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA001323.

**Contexts:** Bronchus / Respiratory epithelial cells; Nasopharynx / Respiratory epithelial cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000181552-EDDM3B/tissue/primary+data). Review records: ihc-20261025-ENSG00000181552.

### EQTN — Neither

`ENSG00000120160` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA015089.

**Contexts:** Adrenal gland / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000120160-EQTN/tissue/primary+data). Review records: ihc-20261025-ENSG00000120160.

### EZHIP — Neither

`ENSG00000187690` · 8 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA004003; HPA006128.

**Contexts:** Bronchus / Respiratory epithelial cells; Cerebral cortex / Neuronal cells; Kidney / Cells in tubules; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000187690-EZHIP/tissue/primary+data). Review records: ihc-20261025-ENSG00000187690.

### FAM170A — Neither

`ENSG00000164334` · 16 positive raw rows; 8 Medium/High; 4 Vital Medium/High. Reagents: HPA037900; HPA037902.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Neuronal cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Glandular cells; Esophagus / Squamous epithelial cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Salivary gland / Glandular cells; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 8 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000164334-FAM170A/tissue/primary+data). Review records: ihc-20261025-ENSG00000164334.

### FAM187B — Neither

`ENSG00000177558` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA079535.

**Contexts:** Kidney / Cells in tubules.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000177558-FAM187B/tissue/primary+data). Review records: ihc-20261025-ENSG00000177558.

### FAM24A — Neither

`ENSG00000203795` · 2 positive raw rows; 2 Medium/High; 1 Vital Medium/High. Reagents: HPA052490.

**Contexts:** Bronchus / Respiratory epithelial cells; Nasopharynx / Respiratory epithelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000203795-FAM24A/tissue/primary+data). Review records: ihc-20261025-ENSG00000203795.

### FBXO24 — Neither

`ENSG00000106336` · 4 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA050951.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Hippocampus / Neuronal cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000106336-FBXO24/tissue/primary+data). Review records: ihc-20261025-ENSG00000106336.

### FKBP6 — Neither

`ENSG00000077800` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA061109.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Cerebellum / Cells in granular layer; Kidney / Cells in tubules.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000077800-FKBP6/tissue/primary+data). Review records: ihc-20261025-ENSG00000077800.

### FMR1NB — Neither

`ENSG00000176988` · 18 positive raw rows; 2 Medium/High; 1 Vital Medium/High. Reagents: CAB026403; HPA011284.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Skin 1 / Keratinocytes; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 2 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000176988-FMR1NB/tissue/primary+data). Review records: ihc-20261025-ENSG00000176988.

### FNDC8 — Neither

`ENSG00000073598` · 8 positive raw rows; 5 Medium/High; 3 Vital Medium/High. Reagents: HPA026521.

**Contexts:** Bronchus / Respiratory epithelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Kidney / Cells in tubules; Lung / Macrophages; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 5 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000073598-FNDC8/tissue/primary+data). Review records: ihc-20261025-ENSG00000073598.

### FSCB — Neither

`ENSG00000189139` · 14 positive raw rows; 7 Medium/High; 5 Vital Medium/High. Reagents: HPA003546.

**Contexts:** Caudate / Neuronal cells; Cerebellum / Cells in molecular layer; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lung / Macrophages; Oral mucosa / Squamous epithelial cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (5 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000189139-FSCB/tissue/primary+data). Review records: ihc-20261025-ENSG00000189139.

### GARIN3 — Neither

`ENSG00000170613` · 7 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA039478; HPA040771.

**Contexts:** Colon / Endothelial cells; Colon / Glandular cells; Esophagus / Squamous epithelial cells; Rectum / Glandular cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000170613-GARIN3/tissue/primary+data). Review records: ihc-20261025-ENSG00000170613.

### GDF9 — Neither

`ENSG00000164404` · 3 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA069146.

**Contexts:** Adrenal gland / Glandular cells; Heart muscle / Cardiomyocytes; Pancreas / Exocrine glandular cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000164404-GDF9/tissue/primary+data). Review records: ihc-20261025-ENSG00000164404.

### GLIPR1L1 — Neither

`ENSG00000173401` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA042095.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Lung / Macrophages.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000173401-GLIPR1L1/tissue/primary+data). Review records: ihc-20261025-ENSG00000173401.

### GOLGA6A — Neither

`ENSG00000159289` · 15 positive raw rows; 5 Medium/High; 2 Vital Medium/High. Reagents: HPA044978; HPA050087.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Lung / Macrophages; Lymph node / Non-germinal center cells; Parathyroid gland / Glandular cells; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 5 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000159289-GOLGA6A/tissue/primary+data). Review records: ihc-20261025-ENSG00000159289.

### GOLGA6C — Neither

`ENSG00000167195` · 15 positive raw rows; 5 Medium/High; 2 Vital Medium/High. Reagents: HPA044978; HPA050087.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Lung / Macrophages; Lymph node / Non-germinal center cells; Parathyroid gland / Glandular cells; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 5 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000167195-GOLGA6C/tissue/primary+data). Review records: ihc-20261025-ENSG00000167195.

### GOLGA6D — Neither

`ENSG00000140478` · 15 positive raw rows; 5 Medium/High; 2 Vital Medium/High. Reagents: HPA044978; HPA050087.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Lung / Macrophages; Lymph node / Non-germinal center cells; Parathyroid gland / Glandular cells; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 5 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000140478-GOLGA6D/tissue/primary+data). Review records: ihc-20261025-ENSG00000140478.

### H2AP — Neither

`ENSG00000187516` · 10 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA003356.

**Contexts:** Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Kidney / Cells in tubules; Oral mucosa / Squamous epithelial cells; Rectum / Glandular cells; Small intestine / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000187516-H2AP/tissue/primary+data). Review records: ihc-20261025-ENSG00000187516.

### HDGFL1 — Neither

`ENSG00000112273` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA063707.

**Contexts:** Small intestine / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000112273-HDGFL1/tissue/primary+data). Review records: ihc-20261025-ENSG00000112273.

### HHIPL2 — Neither

`ENSG00000143512` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA059673.

**Contexts:** Cerebral cortex / Neuronal cells; Hippocampus / Neuronal cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000143512-HHIPL2/tissue/primary+data). Review records: ihc-20261025-ENSG00000143512.

### HORMAD2 — Neither

`ENSG00000176635` · 1 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA003074.

**Contexts:** Cerebral cortex / Neuropil.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000176635-HORMAD2/tissue/primary+data). Review records: ihc-20261025-ENSG00000176635.

### HSFX2 — Neither

`ENSG00000268738` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA063040.

**Contexts:** Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000268738-HSFX2/tissue/primary+data). Review records: ihc-20261025-ENSG00000268738.

### IHO1 — Neither

`ENSG00000173421` · 18 positive raw rows; 4 Medium/High; 3 Vital Medium/High. Reagents: HPA045690.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Spleen / Cells in red pulp; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 4 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000173421-IHO1/tissue/primary+data). Review records: ihc-20261025-ENSG00000173421.

### IRGC — Neither

`ENSG00000124449` · 19 positive raw rows; 8 Medium/High; 4 Vital Medium/High. Reagents: HPA046769; HPA060064.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebellum / Purkinje cells; Colon / Endothelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 8 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000124449-IRGC/tissue/primary+data). Review records: ihc-20261025-ENSG00000124449.

### KIAA1210 — Neither

`ENSG00000250423` · 20 positive raw rows; 10 Medium/High; 4 Vital Medium/High. Reagents: HPA048322.

**Contexts:** Caudate / Neuronal cells; Cerebellum / Cells in molecular layer; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Liver / Hepatocytes; Lung / Alveolar cells; Lung / Macrophages; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skin 1 / Fibroblasts; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** 10 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000250423-KIAA1210/tissue/primary+data). Review records: ihc-20261025-ENSG00000250423.

### KLF17 — Neither

`ENSG00000171872` · 4 positive raw rows; 3 Medium/High; 1 Vital Medium/High. Reagents: HPA024629.

**Contexts:** Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Spleen / Cells in red pulp; Spleen / Cells in white pulp.

**Decision:** 3 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000171872-KLF17/tissue/primary+data). Review records: ihc-20261025-ENSG00000171872.

### LDHAL6B — Neither

`ENSG00000171989` · 12 positive raw rows; 5 Medium/High; 4 Vital Medium/High. Reagents: HPA053962.

**Contexts:** Adipose tissue / Adipocytes; Bone marrow / Hematopoietic cells; Cerebellum / Purkinje cells; Cerebral cortex / Endothelial cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Nasopharynx / Respiratory epithelial cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 5 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000171989-LDHAL6B/tissue/primary+data). Review records: ihc-20261025-ENSG00000171989.

### LRRD1 — Neither

`ENSG00000240720` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA079440.

**Contexts:** Stomach 1 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000240720-LRRD1/tissue/primary+data). Review records: ihc-20261025-ENSG00000240720.

### LYZL4 — Neither

`ENSG00000157093` · 18 positive raw rows; 9 Medium/High; 4 Vital Medium/High. Reagents: HPA038806.

**Contexts:** Adipose tissue / Adipocytes; Caudate / Glial cells; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Alveolar cells; Lung / Macrophages; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Fibroblasts; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 9 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000157093-LYZL4/tissue/primary+data). Review records: ihc-20261025-ENSG00000157093.

### LYZL6 — Neither

`ENSG00000275722` · 14 positive raw rows; 3 Medium/High; 0 Vital Medium/High. Reagents: HPA025734; HPA053073.

**Contexts:** Adrenal gland / Glandular cells; Cerebellum / Purkinje cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Colon / Endothelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Lung / Macrophages; Rectum / Glandular cells; Skeletal muscle / Myocytes; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 3 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000275722-LYZL6/tissue/primary+data). Review records: ihc-20261025-ENSG00000275722.

### MAGEA1 — Neither

`ENSG00000198681` · 6 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: CAB002620.

**Contexts:** Adipose tissue / Adipocytes; Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Hippocampus / Neuronal cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000198681-MAGEA1/tissue/primary+data). Review records: ihc-20261025-ENSG00000198681.

### MAGEA10 — Neither

`ENSG00000124260` · 19 positive raw rows; 8 Medium/High; 4 Vital Medium/High. Reagents: HPA003333.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Glial cells; Cerebellum / Purkinje cells; Cerebral cortex / Glial cells; Duodenum / Glandular cells; Hippocampus / Glial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 8 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000124260-MAGEA10/tissue/primary+data). Review records: ihc-20261025-ENSG00000124260.

### MAGEA4 — Neither

`ENSG00000147381` · 26 positive raw rows; 3 Medium/High; 2 Vital Medium/High. Reagents: CAB040535; HPA021942.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Colon / Endothelial cells; Colon / Glandular cells; Esophagus / Squamous epithelial cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 3 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting. The independent OTI1F9 study supports a validated tumor assay and reports a normal panel negative outside reproductive tissues; it does not directly resolve the different HPA021942/CAB040535 neural and epidermal signals. Preserve afami-cel evidence at product scope.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000147381-MAGEA4/tissue/primary+data). Review records: ihc-20261025-ENSG00000147381.

### MAGEA8 — Neither

`ENSG00000156009` · 19 positive raw rows; 2 Medium/High; 1 Vital Medium/High. Reagents: HPA003998.

**Contexts:** Adipose tissue / Adipocytes; Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebellum / Cells in granular layer; Cerebellum / Cells in molecular layer; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Glandular cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Liver / Cholangiocytes; Nasopharynx / Respiratory epithelial cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000156009-MAGEA8/tissue/primary+data). Review records: ihc-20261025-ENSG00000156009.

### MAGEB1 — Neither

`ENSG00000214107` · 12 positive raw rows; 7 Medium/High; 2 Vital Medium/High. Reagents: HPA001193; HPA002820.

**Contexts:** Colon / Glandular cells; Duodenum / Glandular cells; Oral mucosa / Squamous epithelial cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** 7 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000214107-MAGEB1/tissue/primary+data). Review records: ihc-20261025-ENSG00000214107.

### MAGEB3 — Neither

`ENSG00000198798` · 12 positive raw rows; 5 Medium/High; 4 Vital Medium/High. Reagents: HPA072203.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Nasopharynx / Respiratory epithelial cells; Salivary gland / Glandular cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 5 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000198798-MAGEB3/tissue/primary+data). Review records: ihc-20261025-ENSG00000198798.

### MAGEC1 — Neither

`ENSG00000155495` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA004622.

**Contexts:** Kidney / Cells in tubules.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000155495-MAGEC1/tissue/primary+data). Review records: ihc-20261025-ENSG00000155495.

### MEIOC — Neither

`ENSG00000180336` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA027266.

**Contexts:** Bone marrow / Hematopoietic cells; Spleen / Cells in red pulp.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000180336-MEIOC/tissue/primary+data). Review records: ihc-20261025-ENSG00000180336.

### MGAT4D — Neither

`ENSG00000205301` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA036523.

**Contexts:** Oral mucosa / Squamous epithelial cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000205301-MGAT4D/tissue/primary+data). Review records: ihc-20261025-ENSG00000205301.

### MORC1 — Neither

`ENSG00000114487` · 12 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA036413.

**Contexts:** Adrenal gland / Glandular cells; Colon / Glandular cells; Lung / Alveolar cells; Oral mucosa / Squamous epithelial cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000114487-MORC1/tissue/primary+data). Review records: ihc-20261025-ENSG00000114487.

### MOV10L1 — Neither

`ENSG00000073146` · 16 positive raw rows; 5 Medium/High; 3 Vital Medium/High. Reagents: HPA019034.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Glial cells; Cerebellum / Cells in granular layer; Cerebral cortex / Neuronal cells; Heart muscle / Cardiomyocytes; Hippocampus / Glial cells; Liver / Hepatocytes; Lung / Macrophages; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Skin 2 / Epidermal cells; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 5 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000073146-MOV10L1/tissue/primary+data). Review records: ihc-20261025-ENSG00000073146.

### MROH2B — Neither

`ENSG00000171495` · 15 positive raw rows; 6 Medium/High; 5 Vital Medium/High. Reagents: HPA059457.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Cerebellum / Purkinje cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Alveolar cells; Pancreas / Exocrine glandular cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skin 1 / Keratinocytes; Small intestine / Glandular cells; Stomach 1 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 6 Medium/High normal annotations (5 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000171495-MROH2B/tissue/primary+data). Review records: ihc-20261025-ENSG00000171495.

### NAA11 — Neither

`ENSG00000156269` · 15 positive raw rows; 4 Medium/High; 2 Vital Medium/High. Reagents: HPA035922.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Cells in granular layer; Colon / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Skeletal muscle / Myocytes; Skin 1 / Fibroblasts; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 4 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000156269-NAA11/tissue/primary+data). Review records: ihc-20261025-ENSG00000156269.

### NANOS2 — Neither

`ENSG00000188425` · 27 positive raw rows; 10 Medium/High; 9 Vital Medium/High. Reagents: HPA014190.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Glial cells; Caudate / Neuronal cells; Cerebellum / Cells in granular layer; Cerebellum / Cells in molecular layer; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Glial cells; Hippocampus / Neuronal cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Lung / Macrophages; Lymph node / Non-germinal center cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 10 Medium/High normal annotations (9 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000188425-NANOS2/tissue/primary+data). Review records: ihc-20261025-ENSG00000188425.

### NUP210L — Neither

`ENSG00000143552` · 47 positive raw rows; 18 Medium/High; 11 Vital Medium/High. Reagents: HPA064245; HPA066707.

**Contexts:** Adipose tissue / Adipocytes; Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebellum / Purkinje cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lung / Alveolar cells; Lung / Macrophages; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 18 Medium/High normal annotations (11 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000143552-NUP210L/tissue/primary+data). Review records: ihc-20261025-ENSG00000143552.

### NUTM1 — Neither

`ENSG00000184507` · 4 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA040421.

**Contexts:** Adrenal gland / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Oral mucosa / Squamous epithelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000184507-NUTM1/tissue/primary+data). Review records: ihc-20261025-ENSG00000184507.

### ODF1 — Neither

`ENSG00000155087` · 13 positive raw rows; 8 Medium/High; 7 Vital Medium/High. Reagents: HPA064295.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Pancreas / Pancreatic endocrine cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 8 Medium/High normal annotations (7 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000155087-ODF1/tissue/primary+data). Review records: ihc-20261025-ENSG00000155087.

### ODF4 — Neither

`ENSG00000184650` · 7 positive raw rows; 2 Medium/High; 0 Vital Medium/High. Reagents: HPA072129.

**Contexts:** Cerebral cortex / Glial cells; Cerebral cortex / Neuropil; Lymph node / Germinal center cells; Lymph node / Non-germinal center cells; Skeletal muscle / Myocytes; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000184650-ODF4/tissue/primary+data). Review records: ihc-20261025-ENSG00000184650.

### OOEP — Neither

`ENSG00000203907` · 4 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA055165.

**Contexts:** Cerebral cortex / Neuropil; Kidney / Cells in glomeruli; Rectum / Glandular cells; Smooth muscle / Smooth muscle cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000203907-OOEP/tissue/primary+data). Review records: ihc-20261025-ENSG00000203907.

### PASD1 — Neither

`ENSG00000166049` · 25 positive raw rows; 3 Medium/High; 2 Vital Medium/High. Reagents: HPA011122; HPA011152.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Purkinje cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Lymph node / Germinal center cells; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Rectum / Glandular cells; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Spleen / Cells in red pulp; Spleen / Cells in white pulp; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 3 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000166049-PASD1/tissue/primary+data). Review records: ihc-20261025-ENSG00000166049.

### PDCL2 — Neither

`ENSG00000163440` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA048260.

**Contexts:** Bronchus / Respiratory epithelial cells; Nasopharynx / Respiratory epithelial cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000163440-PDCL2/tissue/primary+data). Review records: ihc-20261025-ENSG00000163440.

### PDHA2 — Neither

`ENSG00000163114` · 93 positive raw rows; 57 Medium/High; 32 Vital Medium/High. Reagents: HPA047487; HPA047864; HPA063053.

**Contexts:** Adipose tissue / Adipocytes; Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Cells in granular layer; Cerebellum / Cells in molecular layer; Cerebellum / Purkinje cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Glial cells; Hippocampus / Neuronal cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Cholangiocytes; Liver / Hepatocytes; Lung / Alveolar cells; Lung / Macrophages; Lymph node / Germinal center cells; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Fibroblasts; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp; Spleen / Cells in white pulp; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 57 Medium/High normal annotations (32 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000163114-PDHA2/tissue/primary+data). Review records: ihc-20261025-ENSG00000163114.

### PGK2 — Neither

`ENSG00000170950` · 31 positive raw rows; 10 Medium/High; 5 Vital Medium/High. Reagents: HPA045385; HPA073644; HPA073656.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebral cortex / Glial cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Endothelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skin 1 / Fibroblasts; Skin 1 / Keratinocytes; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 10 Medium/High normal annotations (5 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000170950-PGK2/tissue/primary+data). Review records: ihc-20261025-ENSG00000170950.

### PKDREJ — Neither

`ENSG00000130943` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA034587.

**Contexts:** Bronchus / Respiratory epithelial cells; Duodenum / Glandular cells; Lung / Macrophages; Rectum / Glandular cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000130943-PKDREJ/tissue/primary+data). Review records: ihc-20261025-ENSG00000130943.

### PLAC1 — Neither

`ENSG00000170965` · 6 positive raw rows; 2 Medium/High; 2 Vital Medium/High. Reagents: HPA063425.

**Contexts:** Adrenal gland / Glandular cells; Caudate / Neuronal cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules.

**Decision:** 2 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000170965-PLAC1/tissue/primary+data). Review records: ihc-20261025-ENSG00000170965.

### POM121L12 — Neither

`ENSG00000221900` · 11 positive raw rows; 3 Medium/High; 3 Vital Medium/High. Reagents: HPA027836.

**Contexts:** Adrenal gland / Glandular cells; Cerebral cortex / Neuropil; Duodenum / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 3 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000221900-POM121L12/tissue/primary+data). Review records: ihc-20261025-ENSG00000221900.

### PRAME — Unsafe

`ENSG00000185686` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA045153.

**Contexts:** Colon / Glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting. A separate gene-resolved publication reports adrenal-cortical protein; that accepted finding takes precedence and makes the final gene category Unsafe.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000185686-PRAME/tissue/primary+data). Review records: ihc-20261025-ENSG00000185686, review-20261025-PRAME-adrenal.

### PRDM7 — Neither

`ENSG00000126856` · 3 positive raw rows; 2 Medium/High; 0 Vital Medium/High. Reagents: HPA059944.

**Contexts:** Colon / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000126856-PRDM7/tissue/primary+data). Review records: ihc-20261025-ENSG00000126856.

### PRKACG — Neither

`ENSG00000165059` · 15 positive raw rows; 4 Medium/High; 3 Vital Medium/High. Reagents: CAB004530.

**Contexts:** Adipose tissue / Adipocytes; Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Colon / Peripheral nerve/ganglion; Kidney / Cells in tubules; Lung / Macrophages; Lymph node / Non-germinal center cells; Pancreas / Exocrine glandular cells; Parathyroid gland / Glandular cells; Spleen / Cells in red pulp; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 4 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000165059-PRKACG/tissue/primary+data). Review records: ihc-20261025-ENSG00000165059.

### PRSS37 — Neither

`ENSG00000165076` · 4 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA020541.

**Contexts:** Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Salivary gland / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000165076-PRSS37/tissue/primary+data). Review records: ihc-20261025-ENSG00000165076.

### PRSS54 — Neither

`ENSG00000103023` · 3 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA041549.

**Contexts:** Adrenal gland / Glandular cells; Duodenum / Glandular cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000103023-PRSS54/tissue/primary+data). Review records: ihc-20261025-ENSG00000103023.

### PSG1 — Neither

`ENSG00000231924` · 19 positive raw rows; 7 Medium/High; 1 Vital Medium/High. Reagents: HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000231924-PSG1/tissue/primary+data). Review records: ihc-20261025-ENSG00000231924.

### PSG11 — Neither

`ENSG00000243130` · 24 positive raw rows; 11 Medium/High; 3 Vital Medium/High. Reagents: HPA042240; HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 11 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000243130-PSG11/tissue/primary+data). Review records: ihc-20261025-ENSG00000243130.

### PSG3 — Neither

`ENSG00000221826` · 19 positive raw rows; 7 Medium/High; 1 Vital Medium/High. Reagents: HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000221826-PSG3/tissue/primary+data). Review records: ihc-20261025-ENSG00000221826.

### PSG6 — Neither

`ENSG00000170848` · 19 positive raw rows; 7 Medium/High; 1 Vital Medium/High. Reagents: HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000170848-PSG6/tissue/primary+data). Review records: ihc-20261025-ENSG00000170848.

### PSG7 — Neither

`ENSG00000221878` · 19 positive raw rows; 7 Medium/High; 1 Vital Medium/High. Reagents: HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000221878-PSG7/tissue/primary+data). Review records: ihc-20261025-ENSG00000221878.

### PSG8 — Neither

`ENSG00000124467` · 19 positive raw rows; 7 Medium/High; 1 Vital Medium/High. Reagents: HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000124467-PSG8/tissue/primary+data). Review records: ihc-20261025-ENSG00000124467.

### PSG9 — Neither

`ENSG00000183668` · 19 positive raw rows; 7 Medium/High; 1 Vital Medium/High. Reagents: HPA046327.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Glial cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Esophagus / Squamous epithelial cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Skin 1 / Fibroblasts; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000183668-PSG9/tissue/primary+data). Review records: ihc-20261025-ENSG00000183668.

### RAD21L1 — Neither

`ENSG00000244588` · 23 positive raw rows; 7 Medium/High; 3 Vital Medium/High. Reagents: HPA053282.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Liver / Hepatocytes; Oral mucosa / Squamous epithelial cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Fibroblasts; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 7 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000244588-RAD21L1/tissue/primary+data). Review records: ihc-20261025-ENSG00000244588.

### RBMXL2 — Neither

`ENSG00000170748` · 23 positive raw rows; 15 Medium/High; 11 Vital Medium/High. Reagents: HPA038912.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Colon / Endothelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Liver / Cholangiocytes; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 15 Medium/High normal annotations (11 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000170748-RBMXL2/tissue/primary+data). Review records: ihc-20261025-ENSG00000170748.

### RBMY1A1 — Neither

`ENSG00000234414` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA053147.

**Contexts:** Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000234414-RBMY1A1/tissue/primary+data). Review records: ihc-20261025-ENSG00000234414.

### RBMY1B — Neither

`ENSG00000242875` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA053147.

**Contexts:** Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000242875-RBMY1B/tissue/primary+data). Review records: ihc-20261025-ENSG00000242875.

### RBMY1D — Neither

`ENSG00000244395` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA053147.

**Contexts:** Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000244395-RBMY1D/tissue/primary+data). Review records: ihc-20261025-ENSG00000244395.

### RBMY1E — Neither

`ENSG00000242389` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA053147.

**Contexts:** Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000242389-RBMY1E/tissue/primary+data). Review records: ihc-20261025-ENSG00000242389.

### RBMY1F — Neither

`ENSG00000169800` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA053147.

**Contexts:** Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000169800-RBMY1F/tissue/primary+data). Review records: ihc-20261025-ENSG00000169800.

### RBMY1J — Neither

`ENSG00000226941` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA053147.

**Contexts:** Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000226941-RBMY1J/tissue/primary+data). Review records: ihc-20261025-ENSG00000226941.

### RHOXF2 — Neither

`ENSG00000131721` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: CAB020794.

**Contexts:** Skeletal muscle / Myocytes.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000131721-RHOXF2/tissue/primary+data). Review records: ihc-20261025-ENSG00000131721.

### RIMBP3 — Neither

`ENSG00000275793` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA001183.

**Contexts:** Cerebellum / Purkinje cells; Cerebral cortex / Glial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000275793-RIMBP3/tissue/primary+data). Review records: ihc-20261025-ENSG00000275793.

### RIMBP3C — Neither

`ENSG00000183246` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA001183.

**Contexts:** Cerebellum / Purkinje cells; Cerebral cortex / Glial cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000183246-RIMBP3C/tissue/primary+data). Review records: ihc-20261025-ENSG00000183246.

### RNASE12 — Neither

`ENSG00000258436` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA044175.

**Contexts:** Small intestine / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000258436-RNASE12/tissue/primary+data). Review records: ihc-20261025-ENSG00000258436.

### RSPH6A — Neither

`ENSG00000104941` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA045382.

**Contexts:** Adrenal gland / Glandular cells; Duodenum / Glandular cells; Small intestine / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000104941-RSPH6A/tissue/primary+data). Review records: ihc-20261025-ENSG00000104941.

### SAGE1 — Neither

`ENSG00000181433` · 30 positive raw rows; 9 Medium/High; 4 Vital Medium/High. Reagents: HPA003033; HPA003208.

**Contexts:** Adipose tissue / Adipocytes; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Purkinje cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skin 1 / Keratinocytes; Smooth muscle / Smooth muscle cells; Spleen / Cells in red pulp; Spleen / Cells in white pulp; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 9 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000181433-SAGE1/tissue/primary+data). Review records: ihc-20261025-ENSG00000181433.

### SEL1L2 — Neither

`ENSG00000101251` · 1 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA059121.

**Contexts:** Cerebellum / Purkinje cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000101251-SEL1L2/tissue/primary+data). Review records: ihc-20261025-ENSG00000101251.

### SLCO6A1 — Neither

`ENSG00000205359` · 4 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA054126.

**Contexts:** Kidney / Cells in tubules; Lymph node / Non-germinal center cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000205359-SLCO6A1/tissue/primary+data). Review records: ihc-20261025-ENSG00000205359.

### SMC1B — Neither

`ENSG00000077935` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA001500.

**Contexts:** Liver / Hepatocytes.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000077935-SMC1B/tissue/primary+data). Review records: ihc-20261025-ENSG00000077935.

### SMIM21 — Neither

`ENSG00000206026` · 7 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA041686.

**Contexts:** Bone marrow / Hematopoietic cells; Heart muscle / Cardiomyocytes; Parathyroid gland / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Spleen / Cells in red pulp.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000206026-SMIM21/tissue/primary+data). Review records: ihc-20261025-ENSG00000206026.

### SOX30 — Neither

`ENSG00000039600` · 40 positive raw rows; 8 Medium/High; 5 Vital Medium/High. Reagents: CAB009382; HPA006159.

**Contexts:** Adipose tissue / Adipocytes; Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Caudate / Neuronal cells; Cerebellum / Cells in molecular layer; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Glandular cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Liver / Hepatocytes; Lung / Macrophages; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Skin 1 / Keratinocytes; Skin 1 / Langerhans; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 8 Medium/High normal annotations (5 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000039600-SOX30/tissue/primary+data). Review records: ihc-20261025-ENSG00000039600.

### SPACA1 — Neither

`ENSG00000118434` · 12 positive raw rows; 5 Medium/High; 4 Vital Medium/High. Reagents: HPA026744; HPA043297.

**Contexts:** Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Colon / Glandular cells; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Lung / Macrophages; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Smooth muscle / Smooth muscle cells.

**Decision:** 5 Medium/High normal annotations (4 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000118434-SPACA1/tissue/primary+data). Review records: ihc-20261025-ENSG00000118434.

### SPACA4 — Neither

`ENSG00000177202` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA041927.

**Contexts:** Pancreas / Exocrine glandular cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000177202-SPACA4/tissue/primary+data). Review records: ihc-20261025-ENSG00000177202.

### SPACA5 — Neither

`ENSG00000171489` · 4 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA044747.

**Contexts:** Caudate / Neuronal cells; Heart muscle / Cardiomyocytes; Kidney / Cells in tubules; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000171489-SPACA5/tissue/primary+data). Review records: ihc-20261025-ENSG00000171489.

### SPAM1 — Neither

`ENSG00000106304` · 5 positive raw rows; 2 Medium/High; 2 Vital Medium/High. Reagents: HPA017984.

**Contexts:** Bronchus / Respiratory epithelial cells; Cerebral cortex / Neuropil; Duodenum / Glandular cells; Skin 2 / Epidermal cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000106304-SPAM1/tissue/primary+data). Review records: ihc-20261025-ENSG00000106304.

### SPATA12 — Neither

`ENSG00000186451` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA075578.

**Contexts:** Kidney / Cells in tubules.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000186451-SPATA12/tissue/primary+data). Review records: ihc-20261025-ENSG00000186451.

### SPATA16 — Neither

`ENSG00000144962` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA037884.

**Contexts:** Cerebellum / Cells in molecular layer; Cerebral cortex / Neuronal cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000144962-SPATA16/tissue/primary+data). Review records: ihc-20261025-ENSG00000144962.

### SPATA19 — Neither

`ENSG00000166118` · 2 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA040082.

**Contexts:** Salivary gland / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000166118-SPATA19/tissue/primary+data). Review records: ihc-20261025-ENSG00000166118.

### SPATA3 — Neither

`ENSG00000173699` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA018254.

**Contexts:** Lung / Macrophages; Pancreas / Pancreatic endocrine cells; Skin 2 / Epidermal cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000173699-SPATA3/tissue/primary+data). Review records: ihc-20261025-ENSG00000173699.

### SPATS1 — Neither

`ENSG00000249481` · 22 positive raw rows; 6 Medium/High; 5 Vital Medium/High. Reagents: HPA049722.

**Contexts:** Cerebellum / Cells in granular layer; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Cerebral cortex / Neuropil; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells; Rectum / Glandular cells; Salivary gland / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 6 Medium/High normal annotations (5 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000249481-SPATS1/tissue/primary+data). Review records: ihc-20261025-ENSG00000249481.

### SPDYE4 — Neither

`ENSG00000183318` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA051750.

**Contexts:** Kidney / Cells in tubules; Thyroid gland / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000183318-SPDYE4/tissue/primary+data). Review records: ihc-20261025-ENSG00000183318.

### SPMIP9 — Neither

`ENSG00000172073` · 10 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA043987.

**Contexts:** Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Duodenum / Glandular cells; Esophagus / Squamous epithelial cells; Lymph node / Germinal center cells; Oral mucosa / Squamous epithelial cells; Rectum / Glandular cells; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000172073-SPMIP9/tissue/primary+data). Review records: ihc-20261025-ENSG00000172073.

### STPG4 — Neither

`ENSG00000239605` · 18 positive raw rows; 6 Medium/High; 5 Vital Medium/High. Reagents: HPA052691.

**Contexts:** Bone marrow / Hematopoietic cells; Bronchus / Respiratory epithelial cells; Cerebellum / Purkinje cells; Cerebral cortex / Endothelial cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Lymph node / Non-germinal center cells; Nasopharynx / Respiratory epithelial cells; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Thyroid gland / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 6 Medium/High normal annotations (5 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000239605-STPG4/tissue/primary+data). Review records: ihc-20261025-ENSG00000239605.

### SUN5 — Neither

`ENSG00000167098` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA048529.

**Contexts:** Rectum / Glandular cells; Smooth muscle / Smooth muscle cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000167098-SUN5/tissue/primary+data). Review records: ihc-20261025-ENSG00000167098.

### SYCP1 — Neither

`ENSG00000198765` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA021083.

**Contexts:** Bronchus / Respiratory epithelial cells; Hippocampus / Neuronal cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000198765-SYCP1/tissue/primary+data). Review records: ihc-20261025-ENSG00000198765.

### TBC1D21 — Neither

`ENSG00000167139` · 8 positive raw rows; 2 Medium/High; 2 Vital Medium/High. Reagents: HPA043319.

**Contexts:** Adrenal gland / Glandular cells; Colon / Glandular cells; Heart muscle / Cardiomyocytes; Oral mucosa / Squamous epithelial cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000167139-TBC1D21/tissue/primary+data). Review records: ihc-20261025-ENSG00000167139.

### TCP11 — Neither

`ENSG00000124678` · 2 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA048311.

**Contexts:** Duodenum / Glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000124678-TCP11/tissue/primary+data). Review records: ihc-20261025-ENSG00000124678.

### TEX13A — Neither

`ENSG00000268629` · 1 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA072257.

**Contexts:** Skeletal muscle / Myocytes.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000268629-TEX13A/tissue/primary+data). Review records: ihc-20261025-ENSG00000268629.

### TEX14 — Neither

`ENSG00000121101` · 3 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA077020.

**Contexts:** Kidney / Cells in tubules; Salivary gland / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000121101-TEX14/tissue/primary+data). Review records: ihc-20261025-ENSG00000121101.

### TEX15 — Neither

`ENSG00000133863` · 8 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA036800.

**Contexts:** Caudate / Neuronal cells; Cerebral cortex / Neuronal cells; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Nasopharynx / Respiratory epithelial cells; Rectum / Glandular cells; Small intestine / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000133863-TEX15/tissue/primary+data). Review records: ihc-20261025-ENSG00000133863.

### TEX19 — Neither

`ENSG00000182459` · 4 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA054577.

**Contexts:** Adrenal gland / Glandular cells; Colon / Glandular cells; Duodenum / Glandular cells; Small intestine / Glandular cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000182459-TEX19/tissue/primary+data). Review records: ihc-20261025-ENSG00000182459.

### TEX44 — Neither

`ENSG00000177673` · 13 positive raw rows; 6 Medium/High; 2 Vital Medium/High. Reagents: HPA049917; HPA056433.

**Contexts:** Adrenal gland / Glandular cells; Cerebral cortex / Neuronal cells; Colon / Glandular cells; Duodenum / Glandular cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Parathyroid gland / Glandular cells; Rectum / Glandular cells; Skeletal muscle / Myocytes; Small intestine / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 6 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000177673-TEX44/tissue/primary+data). Review records: ihc-20261025-ENSG00000177673.

### TGIF2LX — Neither

`ENSG00000153779` · 2 positive raw rows; 2 Medium/High; 2 Vital Medium/High. Reagents: HPA034543.

**Contexts:** Pancreas / Exocrine glandular cells; Pancreas / Pancreatic endocrine cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 2 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000153779-TGIF2LX/tissue/primary+data). Review records: ihc-20261025-ENSG00000153779.

### TMCO5A — Neither

`ENSG00000166069` · 4 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA056530.

**Contexts:** Bronchus / Respiratory epithelial cells; Kidney / Cells in tubules; Stomach 1 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000166069-TMCO5A/tissue/primary+data). Review records: ihc-20261025-ENSG00000166069.

### TMEM210 — Neither

`ENSG00000185863` · 2 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA066907.

**Contexts:** Adrenal gland / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000185863-TMEM210/tissue/primary+data). Review records: ihc-20261025-ENSG00000185863.

### TSBP1 — Neither

`ENSG00000204296` · 8 positive raw rows; 6 Medium/High; 6 Vital Medium/High. Reagents: HPA057426.

**Contexts:** Adrenal gland / Glandular cells; Bone marrow / Hematopoietic cells; Duodenum / Glandular cells; Kidney / Cells in glomeruli; Kidney / Cells in tubules; Rectum / Glandular cells; Skin 2 / Epidermal cells; Small intestine / Glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 6 Medium/High normal annotations (6 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000204296-TSBP1/tissue/primary+data). Review records: ihc-20261025-ENSG00000204296.

### TSPY1 — Neither

`ENSG00000258992` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA049384.

**Contexts:** Adrenal gland / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000258992-TSPY1/tissue/primary+data). Review records: ihc-20261025-ENSG00000258992.

### TSPY10 — Neither

`ENSG00000236424` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA049384.

**Contexts:** Adrenal gland / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000236424-TSPY10/tissue/primary+data). Review records: ihc-20261025-ENSG00000236424.

### TSPY2 — Neither

`ENSG00000168757` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA049384.

**Contexts:** Adrenal gland / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000168757-TSPY2/tissue/primary+data). Review records: ihc-20261025-ENSG00000168757.

### TSPY3 — Neither

`ENSG00000228927` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA049384.

**Contexts:** Adrenal gland / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000228927-TSPY3/tissue/primary+data). Review records: ihc-20261025-ENSG00000228927.

### TSPY4 — Neither

`ENSG00000233803` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA049384.

**Contexts:** Adrenal gland / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000233803-TSPY4/tissue/primary+data). Review records: ihc-20261025-ENSG00000233803.

### TSPY8 — Neither

`ENSG00000229549` · 3 positive raw rows; 0 Medium/High; 0 Vital Medium/High. Reagents: HPA049384.

**Contexts:** Adrenal gland / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. Only Low raw normal annotations remain; faint intensity alone is not evidence that staining is nonspecific. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000229549-TSPY8/tissue/primary+data). Review records: ihc-20261025-ENSG00000229549.

### TSPYL6 — Neither

`ENSG00000178021` · 14 positive raw rows; 7 Medium/High; 2 Vital Medium/High. Reagents: HPA034699; HPA034700.

**Contexts:** Bronchus / Respiratory epithelial cells; Cerebellum / Purkinje cells; Colon / Glandular cells; Duodenum / Glandular cells; Lung / Macrophages; Lymph node / Germinal center cells; Lymph node / Non-germinal center cells; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Skin 1 / Keratinocytes; Skin 2 / Epidermal cells; Small intestine / Glandular cells; Spleen / Cells in red pulp; Thyroid gland / Glandular cells.

**Decision:** 7 Medium/High normal annotations (2 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000178021-TSPYL6/tissue/primary+data). Review records: ihc-20261025-ENSG00000178021.

### TSSK1B — Neither

`ENSG00000212122` · 1 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA027827.

**Contexts:** Bone marrow / Hematopoietic cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000212122-TSSK1B/tissue/primary+data). Review records: ihc-20261025-ENSG00000212122.

### TTC23L — Neither

`ENSG00000205838` · 20 positive raw rows; 3 Medium/High; 3 Vital Medium/High. Reagents: HPA037807.

**Contexts:** Adrenal gland / Glandular cells; Bronchus / Respiratory epithelial cells; Cerebellum / Purkinje cells; Cerebral cortex / Neuronal cells; Colon / Endothelial cells; Colon / Glandular cells; Colon / Peripheral nerve/ganglion; Duodenum / Glandular cells; Heart muscle / Cardiomyocytes; Hippocampus / Neuronal cells; Kidney / Cells in tubules; Liver / Hepatocytes; Nasopharynx / Respiratory epithelial cells; Oral mucosa / Squamous epithelial cells; Pancreas / Exocrine glandular cells; Rectum / Glandular cells; Small intestine / Glandular cells; Stomach 1 / Glandular cells; Stomach 2 / Glandular cells; Urinary bladder / Urothelial cells.

**Decision:** 3 Medium/High normal annotations (3 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000205838-TTC23L/tissue/primary+data). Review records: ihc-20261025-ENSG00000205838.

### TTLL13 — Neither

`ENSG00000213471` · 1 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA041305.

**Contexts:** Adrenal gland / Glandular cells.

**Decision:** 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000213471-TTLL13/tissue/primary+data). Review records: ihc-20261025-ENSG00000213471.

### UBQLN3 — Neither

`ENSG00000175520` · 3 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA020547; HPA062018.

**Contexts:** Colon / Endothelial cells; Heart muscle / Cardiomyocytes; Skin 2 / Epidermal cells.

**Decision:** 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000175520-UBQLN3/tissue/primary+data). Review records: ihc-20261025-ENSG00000175520.

### USP26 — Neither

`ENSG00000134588` · 2 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA027310.

**Contexts:** Lymph node / Non-germinal center cells; Pancreas / Exocrine glandular cells.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000134588-USP26/tissue/primary+data). Review records: ihc-20261025-ENSG00000134588.

### VCX — Neither

`ENSG00000182583` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA044996.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Pancreas / Exocrine glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000182583-VCX/tissue/primary+data). Review records: ihc-20261025-ENSG00000182583.

### VCX3A — Neither

`ENSG00000169059` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA044996.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Pancreas / Exocrine glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000169059-VCX3A/tissue/primary+data). Review records: ihc-20261025-ENSG00000169059.

### VCX3B — Neither

`ENSG00000205642` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA044996.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Pancreas / Exocrine glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000205642-VCX3B/tissue/primary+data). Review records: ihc-20261025-ENSG00000205642.

### VCY — Neither

`ENSG00000129864` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA044996.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Pancreas / Exocrine glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000129864-VCY/tissue/primary+data). Review records: ihc-20261025-ENSG00000129864.

### VCY1B — Neither

`ENSG00000129862` · 5 positive raw rows; 1 Medium/High; 1 Vital Medium/High. Reagents: HPA044996.

**Contexts:** Adrenal gland / Glandular cells; Kidney / Cells in tubules; Liver / Hepatocytes; Pancreas / Exocrine glandular cells; Small intestine / Glandular cells.

**Decision:** The source warns that the reagent targets multiple genes; a family signal cannot prove which locus produces the normal protein. 1 Medium/High normal annotations (1 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000129862-VCY1B/tissue/primary+data). Review records: ihc-20261025-ENSG00000129862.

### ZCCHC13 — Neither

`ENSG00000187969` · 3 positive raw rows; 1 Medium/High; 0 Vital Medium/High. Reagents: HPA003382.

**Contexts:** Lymph node / Germinal center cells; Lymph node / Non-germinal center cells; Spleen / Cells in white pulp.

**Decision:** The profile discusses presumed off-target binding, but does not provide a specimen-and-reagent-specific resolution of every retained normal signal. 1 Medium/High normal annotations (0 in Vital compartments) remain specificity conflicts. The curated negative profile and passing RNA screen do not identify the molecule responsible for these raw stains. This evidence does not establish clinical toxicity, but is insufficient to clear unrestricted gene-level targeting.

**To reconsider:** Resolve the named normal-cell signal using identity-specific endogenous protein or reagent validation; quantify any confirmed normal expression. A trial or another negative bulk atlas alone is not a rebuttal.

[Pinned primary data](https://v23.proteinatlas.org/ENSG00000187969-ZCCHC13/tissue/primary+data). Review records: ihc-20261025-ENSG00000187969.

## Reproduce

Run `python scripts/build_cta_adjudication_report.py` after rebuilding the current policy. This renders reviewed records; it does not perform new adjudication. [Checksums](run-manifest.json).
