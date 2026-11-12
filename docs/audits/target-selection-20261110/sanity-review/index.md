# Gene-by-gene sanity review of the candidate sets

Every one of the 899 Target candidates was checked against four sources that the classification does not use in this form: GTEx v11 tissue medians, HPA v23 single-cell clusters, HPA v23 immunohistochemistry and all 1,093,048 Tabula Sapiens UMI cells (23 donors, 30 tissue groups, immune and stromal cells included). Review notes are hand-written, one per protein, in [manual-notes.csv](manual-notes.csv). [Full table](candidates.csv). These are outputs of screening rules, not validated tumor antigens.

## Result

- **660 of 899** candidates show nothing in any of the four sources: under 1 TPM in every non-reproductive GTEx tissue, under 3 nTPM in every non-reproductive HPA cell cluster, no positive normal IHC and no Tabula Sapiens population with 1% of cells carrying two molecules.
- **108** proteins have a review note; 95 of them are candidates in the current catalog.
- **116** Cancer Antigens rest on a tumor maximum under 5 TPM or on at most two patients.

## Families whose normal tissue no panel samples

Absence from the sampled tissues is not absence from the body. These families pass because their tissue is missing, and most qualify as Cancer Antigens on sparse tumor counts.

| Family | Candidates | CTA | Cancer Antigen only | Median tumor max (TPM) | Median patients above normal | Normal tissue |
|---|---|---|---|---|---|---|
| Cornified-envelope proteins | 2 | 0 | 2 | 144.1 | 54.5 | Upper epidermis. |
| Enamel and dentin proteins | 3 | 1 | 2 | 1695.6 | 75.0 | Tooth; adult ameloblasts are gone, odontoblasts are not. |
| Hair keratins and keratin-associated proteins | 42 | 8 | 33 | 35.5 | 9.5 | Hair follicle; treated as Disposable, replay still open (#687). |
| Male reproductive-tract secretions | 44 | 42 | 1 | 26.2 | 16.0 | Epididymis, seminal vesicle and prostate. |
| Olfactory receptors | 164 | 25 | 137 | 16.6 | 5.0 | Olfactory epithelium is not sampled. |
| Taste, vomeronasal and trace-amine receptors | 13 | 1 | 11 | 21.0 | 14.0 | Taste buds and nasal epithelium are not sampled. |
| Unnamed gene records | 24 | 10 | 8 | 30.8 | 0.0 | No symbol; several lie inside another gene. |

## Reviewed proteins

### Normal non-reproductive cells express it in an independent atlas

39 reviewed; 32 are candidates in the current catalog.

| Protein | Current class | What I saw |
|---|---|---|
| APOL5 | Cancer Antigen only | Hepatocytes 11% in 2 donors |
| ARMH2 | CTA | Lung multiciliated cells 17% in 2 donors; normal reference already 9.1 TPM |
| C5orf52 | Unsafe | Colon crypt stem cells 30% in 3 donors; HPA Paneth 7.7 nTPM; GTEx colon 3 TPM |
| CCDC116 | CTA | GTEx >= 1 TPM in 14 tissues |
| CCDC172 | Unsafe | Tabula Sapiens photoreceptors 57% (32% multi-molecule) in 4 donors; HPA single-cell rods 3.1 nTPM |
| CCN6 | CTA | Cartilage gene (loss causes progressive pseudorheumatoid dysplasia); GTEx 3.5 TPM in 4 tissues |
| CPSF4L | CTA | Glandular basal cells 31% (13% multi-molecule) in 8 donors |
| CSTL1 | CTA | Airway multiciliated cells 10% in 2 donors |
| DCST1 | CTA | Photoreceptors 28% in 6 donors; HPA cones 8.5 nTPM; GTEx skin 3.2 TPM |
| EPPIN-WFDC6 | CTA | Read-through record; lung multiciliated cells 11% in 2 donors |
| EVX1 | CTA | Bladder urothelium 12% in 2 donors; GTEx bladder 2.3 TPM |
| FKBP6 | CTA | Airway multiciliated cells 9% in 3 donors |
| FSCN3 | CTA | Paneth cells 8% in 3 donors |
| HOXC12 | Cancer Antigen only | GTEx skin 5.4 TPM in 4 tissues; posterior HOX gene of adult skin and limb |
| IL36B | Cancer Antigen only | Keratinocyte cytokine: GTEx skin 7.2 TPM, tongue squamous cells 5% in 2 donors |
| IZUMO2 | CTA | Corneal epithelium 12% in 3 donors |
| KRTAP5-10 | CTA | Bladder urothelium 13%; GTEx bladder 1.7 TPM |
| KRTAP5-9 | CTA | GTEx >= 1 TPM in 18 tissues; shares sequence with KRTAP5 paralogs |
| MMP20 | Unsafe | Paneth cells 13% in 5 donors; enamel protease otherwise |
| MROH2B | CTA | Hepatocytes 14% in 5 donors; HPA hepatocytes 2.4 nTPM |
| NKX2-6 | CTA | Cardiomyocytes 61% in one Tabula Sapiens donor; GTEx atrial appendage 0.79 TPM |
| OR8D4 | CTA | Pancreatic acinar cells 23% in 5 donors |
| POTEF | CTA | Ventricular cardiomyocytes 21% in one donor; shares sequence with actin-domain POTE paralogs |
| RAD21L1 | CTA | Photoreceptors 14% in 4 donors; exon inside another gene |
| RIPPLY3 | CTA | Alveolar type 2 cells 6% in 4 donors |
| SCGB1D1 | Unsafe | Ocular keratocytes 67% (43% multi-molecule) in 3 donors; a tear and lacrimal secretoglobin |
| SMIM9 | CTA | Cardiac endothelium 13% in 6 donors |
| SPACA3 | Unsafe | Pancreatic acinar cells 31% in 5 donors; GTEx pancreas 4.7 TPM |
| SPINT4 | Unsafe | Photoreceptors 19% in 2 donors; low IHC in 6 normal tissues |
| SRMS | CTA | GTEx >= 1 TPM in 15 tissues (small intestine 4.6); urothelium 10% in 3 donors |
| STPG3 | CTA | GTEx >= 1 TPM in 16 tissues (thyroid 5.7) |
| TBX22 | CTA | GTEx thyroid 14.5 TPM; cardiac fibroblasts 11% in 2 donors |
| TDRD15 | Unsafe | Hepatocytes 32% (13% multi-molecule) in 7 donors |
| TEX46 | CTA | Airway multiciliated cells 9% in 2 donors |
| TMEM270 | CTA | Bladder urothelium 8% (4% multi-molecule) in 4 donors |
| TUBB8 | CTA | Lung multiciliated cells 10% in 4 donors; shares sequence with other tubulins |
| WFDC8 | CTA | Glandular basal cells 11% in 2 donors |
| ZACN | CTA | GTEx colon 13.2 TPM and >= 1 TPM in 4 tissues |
| ZC3H11C | CTA | GTEx >= 1 TPM in 35 tissues; shares sequence with ZC3H11A |

### Blood or immune cells express it in Tabula Sapiens (outside the donor panel by design)

14 reviewed; 12 are candidates in the current catalog.

| Protein | Current class | What I saw |
|---|---|---|
| C18orf63 | Unsafe | Neutrophils 15% in 5 donors |
| CAGE1 | CTA | Plasma cells 20% in 7 donors; exon inside another gene |
| CBY2 | CTA | Plasma cells 17% in 7 donors |
| CHCT1 | CTA | Monocytes 21% in 6 donor populations |
| CXorf58 | CTA | Monocytes 15% in 10 donor populations |
| NUTM1 | CTA | Bone-marrow neutrophils 40% (17% multi-molecule) in 4 donors; immune cells are outside the donor panel |
| PRR5-ARHGAP8 | CTA | Read-through record; plasma cells 31% in 10 donors |
| PSG8 | CTA | Neutrophils 29% in 7 donors; shares sequence with other PSG genes |
| PSMA8 | CTA | Plasma cells 14% in 9 donors; HPA plasma cells 2.8 nTPM |
| SLCO6A1 | Unsafe | Monocytes 17% in 8 donor populations |
| TXNDC8 | CTA | Monocytes 15% in 6 donor populations |
| XKR3 | CTA | Macrophages 5% in 4 donors |
| ZAR1L | CTA | Plasma cells 8% in 4 donors |
| ZNF492 | CTA | Erythroid progenitors 40% (23% multi-molecule) in 3 donors; HPA erythroid 7.2 nTPM |

### Counts are likely those of a parent or host gene, in tumors as much as in normal tissue

20 reviewed; 20 are candidates in the current catalog.

| Protein | Current class | What I saw |
|---|---|---|
| AK4P3 | CTA | Processed copy of AK4; fibroblasts 23% in 12 donors |
| ATF7-NPFF | CTA | Read-through record; absent from GTEx |
| ENSG00000248235 | Cancer Antigen only | Unnamed record inside another gene; absent from GTEx; no resolved patient |
| ENSG00000249773 | Reproductive Antigen only | Unnamed record inside another gene; fibroblasts in 3 donors |
| ENSG00000255432 | CTA | Unnamed record inside another gene; macrophages in 3 donors |
| ENSG00000255835 | CTA | Unnamed record inside another gene; fibroblasts in 6 donors; no resolved patient |
| ENSG00000256349 | Cancer Antigen only | Unnamed record inside another gene; absent from GTEx; no resolved patient |
| ENSG00000256825 | Cancer Antigen only | Unnamed record inside another gene; myofibroblasts in 8 donors; no resolved patient |
| ENSG00000258691 | CTA | Unnamed record inside another gene; myofibroblasts 24% in 9 donors; no resolved patient |
| ENSG00000266202 | Cancer Antigen only | Unnamed record inside another gene; monocytes in 8 donor populations; no resolved patient |
| ENSG00000267140 | CTA | Unnamed record inside another gene; absent from GTEx; no resolved patient |
| ENSG00000267618 | CTA | Unnamed record inside another gene; absent from GTEx; no resolved patient |
| ENSG00000270299 | Cancer Antigen only | Unnamed record inside another gene; absent from GTEx; no resolved patient |
| ENSG00000279983 | CTA | Unnamed record inside another gene; fibroblasts in 6 donors |
| ENSG00000284934 | CTA | Unnamed record inside another gene; no resolved patient |
| HMGB1P1 | CTA | Processed copy of HMGB1; cardiomyocytes 16% in 5 donors and 1,797 tumor patients are what HMGB1 reads would produce |
| PABPC3 | CTA | Retrocopy of PABPC1: low IHC in 22 normal tissues and macrophages in 3 donors fit cross-reaction and shared reads |
| RPS4Y2 | CTA | Y-linked ribosomal paralog; normal reference 14.3 TPM, B cells 7% in 3 donors |
| URGCP-MRPS24 | CTA | Read-through record; absent from GTEx; tumor counts match host genes |
| ZNF816-ZNF321P | CTA | Read-through record; 7,261 tumor patients and urothelium in 5 donors match host-gene expression |

### Known normal tissue is absent or diluted in every panel

20 reviewed; 17 are candidates in the current catalog.

| Protein | Current class | What I saw |
|---|---|---|
| AWAT1 | Unsafe | Sebaceous-gland wax synthase; low IHC in skin |
| BPIFB3 | Cancer Antigen only | Nasal BPI-fold protein; GTEx pituitary 1.3 TPM |
| BPIFB6 | Cancer Antigen only | Nasal and salivary BPI-fold protein; GTEx minor salivary gland 1.6 TPM |
| C7orf33 | Unsafe | HPA cone photoreceptors 10.2 nTPM |
| DSPP | Cancer Antigen only | Dentin and inner ear; mutations cause dentinogenesis imperfecta with hearing loss |
| GALP | CTA | Neuropeptide of a small hypothalamic arcuate population |
| GFRA4 | Cancer Antigen only | Thyroid C cells, a population too rare for bulk thyroid |
| GRXCR1 | Unsafe | Inner-ear hair-cell gene; loss causes human deafness DFNB25, the same pattern that made OTOG Unsafe |
| OBP2A | CTA | Nasal and lacrimal odorant-binding protein; HPA secretory cells 44 nTPM |
| OCM2 | CTA | Oncomodulin of cochlear outer hair cells |
| OCSTAMP | Cancer Antigen only | Osteoclast fusion protein; bone is not sampled |
| OTOP1 | Cancer Antigen only | Inner-ear and taste-cell proton channel |
| PANX3 | Cancer Antigen only | Cartilage and bone channel; neither tissue is sampled |
| SCGB1C1 | Cancer Antigen only | Nasal mucosa secretoglobin |
| SCGB1D4 | CTA | HPA ciliated cells 1,029 nTPM (fallopian tube); airway ciliated cells not excluded |
| SP7 | CTA | Osterix, required in adult osteoblasts; bone is not in any normal panel |
| SPATA31F3 | CTA | HPA cone photoreceptors 8.5 nTPM |
| TBXT | CTA | Notochord remnants of the adult intervertebral disc; GTEx pituitary 1.7 TPM |
| TMEM95 | CTA | HPA Schwann cells 12.2 nTPM |
| TSGA10IP | CTA | HPA melanocytes 24 nTPM |

### Normal IHC without any RNA; probably antibody cross-reaction

7 reviewed; 7 are candidates in the current catalog.

| Protein | Current class | What I saw |
|---|---|---|
| ACTRT2 | CTA | Low IHC in small intestine and stomach with no RNA; actin-related |
| CT45A2/8/9 | CTA | Low IHC in heart muscle with no RNA; CT45 paralogs share epitopes |
| DCAF8L1 | CTA | Low IHC in bronchus and nasopharynx with no RNA |
| MAGEB4 | CTA | Low IHC in 13 normal tissues with no RNA anywhere: antibody cross-reaction is the likelier reading |
| OOSP2 | CTA | Low IHC in 8 normal tissues with no RNA |
| TOMM20L | CTA | Low IHC in duodenum and small intestine; TOMM20 cross-reaction likely |
| TSPY9 | CTA | Low IHC in adrenal and stomach with no RNA |

### Explained by a permitted tissue or a known technical effect

8 reviewed; 7 are candidates in the current catalog.

| Protein | Current class | What I saw |
|---|---|---|
| BTN1A1 | CTA | Milk-fat globule protein; plasma cells 10% in 3 donors is unexplained |
| LALBA | CTA | Milk protein of lactating breast; breast tumors reach 79,700 TPM |
| PRAME | Unsafe | Known low adrenal and endometrial expression (GTEx adrenal 2.5 TPM); clinically targeted |
| PRM1 | CTA | Protamine; GTEx >= 1 TPM in 10 tissues is spill-over from very high testis signal in sequencing batches |
| PRM2 | CTA | Protamine; as PRM1 |
| RLN1 | CTA | Prostate relaxin |
| SEMG2 | CTA | Seminal-vesicle protein |
| TGM4 | CTA | Prostate transglutaminase; prostate is treated as a permitted tissue |

## Limits

- Tabula Sapiens is one study; its donors share libraries and a reference, so a signal there is replicated across donors but not across studies.
- A monocyte population of one donor (TSP21 lymph node, 237 cells) shows low counts for many testis genes; single-population signals there were not treated as findings.
- The notes use general knowledge of each gene. They are a reviewer's reading, not a literature search, and they change no category.
