# Current policy and IHC matching

The baseline RNA grid is Vital 4/8/16/32, Protected 8/16/32/64, Manageable 16/32/64/128 across missing, Approved-negative, Supported-negative and Enhanced-negative IHC. There is one boundary; equality passes. Above-ceiling Vital/Protected evidence is Unsafe. Manageable excess blocks CTA but does not automatically establish Neither.

Reviewed TCR-T evidence can add one RNA step, capped at Enhanced: Vital 8/16/32/32, Protected 16/32/64/64, Manageable 32/64/128/128. The IHC label stays unchanged. Clinical credit carries the exact peptide/HLA/product limits and does not override normal-protein or clinical-injury exclusions. [Trial and SAE review](../cta-tcr-trials-20261029/index.md).

## RNA/IHC correspondence

The matcher records **exact_cell**, **reviewed_compartment**, **same_tissue_unresolved_cell**, or **unmatched**. Its source/version-specific crosswalk is gene-independent. Exact adult tissue/cell matches and declared histological compartment matches can receive negative-IHC credit; unrelated same-tissue cells, stimulated immune cells, fetal samples and other species cannot.

The new compartment mapping relates intestinal epithelial populations to HPA glandular cells, and bronchial basal/club/ciliated cells to respiratory epithelial cells. This is an anatomical interpretation of the histological compartment, not direct proof that every rare/stem cell was observed. The source grade is retained, and the coverage limitation is explicit. Undifferentiated intestinal populations retain a progenitor caution. This assumption affects selection and requires empirical validation; a crosswalk does not improve antibody sensitivity.

[HPA colon histology](https://www.proteinatlas.org/learn/dictionary/normal/colon), [HPA bronchus histology](https://www.proteinatlas.org/learn/dictionary/normal/bronchus), [all matched measurements](../cta-selection-20261028/compartment-ihc-matches.csv.gz).

PRAME colon RNA is 20.3 nTPM with 96 reads. It receives the declared epithelial-compartment match to Enhanced-negative glandular IHC. Its Manageable ceiling stays capped at 128. PRAME's adrenal bulk 3.6 nTPM has a missing-IHC baseline ceiling of 4, raised to 8 by the IMA203 clinical credit. Glandular-cell staining does not supply complete coverage of the specified adrenal medullary and cortical subcompartments; the clinical credit does not change that IHC finding. Low reported adrenal protein remains a positive flag, not negative evidence.

CTAG2/LAGE-1 has Enhanced-negative heart cardiomyocyte IHC and bronchial respiratory-epithelial IHC. The revised matcher applies the latter to the basal respiratory compartment; the report retains the parent-compartment coverage assumption. Negative IHC is counterevidence, not proof of clinical safety.

## Anatomy and protein

Identified fibroblasts and large-bowel epithelial cells use Manageable in both RNA and HPA protein assessment. They are not called disposable. Other cell types retain their anatomical tiers; organ-level bulk is not assigned the tier of one constituent cell type. Accepted non-Low Vital/Protected protein can independently supply Unsafe evidence. Manageable protein blocks CTA only. Uncertain staining is unresolved; explicitly Low protein remains flagged under the existing policy.

Selective fibroblast depletion in mouse contexts is not evidence that all human tissue-resident fibroblasts are dispensable. [Kuwabara et al. (2022)](https://elifesciences.org/articles/69854). Epithelial targeting can cause severe colitis: [Parkhurst et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3048186/). The tier choices describe this screening policy's tolerance assumptions, not validated toxicity predictions.

## Curated developmental nominations

The current developmental inputs are Tyser gastrulation and Petropoulos preimplantation expression. Newly identified literature sources remain outside current counts: [Shaw et al. (2024)](https://doi.org/10.1038/s41467-024-47649-y) reports 157 cancer-specific targets with six oncofetal annotations, including exon-specific FN1/TNC evidence; [the TRIM71 liver study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11388079/) contains a liver oncofetal signature. [Cardoso-Moreira et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6658352/) supplies developmental organ trajectories with child/adult comparisons. Their identities, feature scope and developmental/postnatal measurements require ingestion before contributing nominations or exclusion evidence.
