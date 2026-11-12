# Rare normal-cell RNA by adult donor

**1,804 proteins were graded from 13,581,391 adult cells and nuclei of
1,779 donors in 93 collections, plus the directly analysed
Siletti cerebellar granule nuclei. 132 have established normal-cell RNA and are
Unsafe; 264 have a signal that is dismissed with its reason.** Category
consequences are in the [catalog rebuild](../../target-selection-20261110/rebuild/index.md).

A count matrix cannot show where a read aligned. One mis-assigned read per cell, from a
paralog, a pseudogene or an overlapping gene, looks exactly like sparse expression. The
grade therefore ignores single-molecule cells and asks five questions. Every outcome is
terminal; nothing waits for review.

| Question | Requirement |
|---|---|
| Is there a signal? | >= 1% of a donor's cells and >= 10 cells carry >= 2 molecules. Protected populations need 2 times these frequencies. |
| Is it replicated? | The signal in >= 2 independent studies of the same organ. Donors of one study share its libraries and reference. Without replication, a signal counts only when >= 5% of a donor's cells carry multiple molecules. |
| Could the reads belong elsewhere? | Shared transcript sequence with another locus, or an exon inside another same-strand gene. Counts for a flagged gene cannot be attributed to it. |
| Does another assay agree? | Bulk RNA >= 1 native TPM in the same organ in both the HPA and GTEx sample sets, a different assay with its own read assignment. |
| Is it published? | Where a reviewed publication reports the gene in the same organ by a gene-specific assay, 1% detection in one donor is established, single molecules included: two assay types then agree. Published protein in a Vital or Protected tissue excludes through the knowledge layer without any counts. |

| Outcome | Rule |
|---|---|
| Established | A replicated signal at a locus without mapping risk, or any signal with bulk confirmation. Unsafe. |
| Dismissed | A signal that is neither: counts at a flagged locus, or from one study, that no other assay confirms. Recorded with its reason; no category effect. |
| Below the screen | No replicated or frequent multi-molecule signal. No finding. |

Below the screen is not a demonstration of absence or safety. The
[literature review](../../target-selection-20261110/rebuild/literature.md) covers every
candidate whose signal is dismissed or who left the pooled flag.
Donors are never pooled, frequencies are never added across paralogs, and unmeasured
genes stay unassessed.

## Grades

| Grade | Tier | Catalog action | Proteins |
|---|---|---|---:|
| below_screen | — | none | 1,372 |
| dismissed | protected | none | 12 |
| dismissed | vital | none | 252 |
| established | protected | unsafe | 5 |
| established | vital | unsafe | 127 |
| unassessed | — | none | 36 |

## Mapping risk

Of 1,862 extracted genes, 612 share at least
10% of their transcript 31-mers, or of their 3′ kilobase,
with another locus, and 638 have an exon inside another same-strand
gene; 1,113 carry either flag. Ensembl 115, sense strand,
primary chromosomes. [Flags, shared fractions and partner loci for every gene](gene-mapping-risk.csv.gz).
A flag marks a possible
confounder, not a demonstrated mis-assignment, and an unflagged gene can still collect
intronic reads.

## Sensitivity to the detection cutoff

Reported comparisons only, with at least ten positive cells in the donor population.
The first two columns count any molecule; they show why single-molecule detection is
not graded.

| Tier | Cutoff | Proteins above it in any donor | In two or more studies | With ≥2 molecules per cell, any donor | With ≥2 molecules per cell, two or more studies |
|---|---:|---:|---:|---:|---:|
| vital | 1% | 1,102 | 886 | 647 | 396 |
| vital | 2% | 993 | 771 | 589 | 323 |
| vital | 5% | 839 | 564 | 470 | 231 |
| protected | 1% | 672 | 356 | 249 | 95 |
| protected | 2% | 577 | 275 | 200 | 70 |
| protected | 5% | 378 | 172 | 128 | 53 |

## Established

| Protein | Organ: strongest population | Studies / donors with multi-molecule cells | Strongest donor: cells with ≥2 molecules | Bulk TPM reached in both HPA and GTEx | Basis |
|---|---|---:|---:|---:|---|
| SPAG9 | small intestine: enterocyte of epithelium proper of ileum | 8 / 27 | 2,255/2,263 (99.6%) | 23 | replicated signal at a locus without mapping risk |
| TSGA10 | brain: dopaminergic neuron | 23 / 432 | 1,757/1,912 (91.9%) | 1.9 | replicated signal at a locus without mapping risk |
| REDIC1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 19 / 209 | 1,722/1,938 (88.9%) | 0 | replicated signal at a locus without mapping risk |
| SPAG1 | brain: choroid plexus epithelial cell | 22 / 265 | 4,740/6,146 (77.1%) | 7.9 | replicated signal at a locus without mapping risk |
| SPA17 | lung: multiciliated columnar cell of tracheobronchial tree | 2 / 5 | 1,762/2,359 (74.7%) | 5.2 | replicated signal at a locus without mapping risk |
| ROPN1L | lung: ciliated cell | 2 / 5 | 946/1,287 (73.5%) | 4.9 | signal confirmed by bulk RNA (exon inside another same-strand gene) |
| UBE2U | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 9 | 745/1,021 (73.0%) | 0.1 | replicated signal at a locus without mapping risk |
| SYCP2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 20 / 393 | 944/1,362 (69.3%) | 4.3 | replicated signal at a locus without mapping risk |
| CSNK2A2IP | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 11 / 152 | 10,449/15,320 (68.2%) | 0.088 | replicated signal at a locus without mapping risk |
| GRXCR2 | heart: regular atrial cardiac myocyte | 5 / 30 | 3,539/5,328 (66.4%) | 0.2 | replicated signal at a locus without mapping risk |
| MTCP1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 13 / 309 | 823/1,362 (60.4%) | 1.9 | signal confirmed by bulk RNA (exon inside another same-strand gene) |
| LIPI | brain: glutamatergic neuron | 12 / 143 | 3,975/6,853 (58.0%) | 0.43 | replicated signal at a locus without mapping risk |
| OOEP | lung: pulmonary alveolar type 2 cell | 5 / 19 | 587/1,181 (49.7%) | 0.2 | replicated signal at a locus without mapping risk |
| C7orf33 | eye: retinal cone cell | 2 / 41 | 1,079/2,436 (44.3%) | — | replicated signal at a locus without mapping risk |
| PBK | small intestine: transit amplifying cell of small intestine | 1 / 1 | 2,105/4,928 (42.7%) | 3.3 | signal confirmed by bulk RNA |
| ADAM7 | eye: OFF midget ganglion cell | 2 / 46 | 606/1,440 (42.1%) | — | replicated signal at a locus without mapping risk |
| MORC1 | brain: hippocampal pyramidal neuron | 2 / 6 | 585/1,451 (40.3%) | 0.1 | replicated signal at a locus without mapping risk |
| MMP26 | eye: non-pigmented ciliary epithelial cell | 3 / 51 | 443/1,110 (39.9%) | — | replicated signal at a locus without mapping risk |
| MAGEC3 | brain: dopaminergic neuron | 19 / 366 | 731/1,912 (38.2%) | 1.5 | replicated signal at a locus without mapping risk |
| GDPD4 | liver: hepatocyte | 2 / 3 | 1,989/5,256 (37.8%) | 0.1 | replicated signal at a locus without mapping risk |
| CEP55 | small intestine: transit amplifying cell of small intestine | 1 / 1 | 1,775/4,928 (36.0%) | 2.9 | signal confirmed by bulk RNA |
| PIWIL4 | brain: L6b glutamatergic cortical neuron | 17 / 364 | 503/1,471 (34.2%) | 5.6 | signal confirmed by bulk RNA (exon inside another same-strand gene) |
| HORMAD1 | eye: basal cell | 3 / 7 | 966/2,848 (33.9%) | — | replicated signal at a locus without mapping risk |
| LDHC | eye: endothelial cell of venule | 9 / 36 | 367/1,148 (32.0%) | — | replicated signal at a locus without mapping risk |
| SYCP3 | brain: astrocyte | 17 / 179 | 617/1,951 (31.6%) | 1.6 | replicated signal at a locus without mapping risk |
| LINC03040 | eye: retinal pigment epithelial cell | 2 / 12 | 573/1,818 (31.5%) | — | replicated signal at a locus without mapping risk |
| RGS21 | eye: Mueller cell | 5 / 50 | 345/1,104 (31.2%) | — | replicated signal at a locus without mapping risk |
| PIWIL1 | eye: retinal cone cell | 13 / 120 | 430/1,398 (30.8%) | — | replicated signal at a locus without mapping risk |
| METTL21EP | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 5 / 74 | 631/2,201 (28.7%) | — | replicated signal at a locus without mapping risk |
| SH2D4B | brain: glutamatergic neuron | 9 / 105 | 740/2,606 (28.4%) | 0.1 | replicated signal at a locus without mapping risk |
| TDRD1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 9 / 91 | 527/1,938 (27.2%) | 0.17 | replicated signal at a locus without mapping risk |
| RBM46 | eye: retinal pigment epithelial cell | 3 / 22 | 485/1,818 (26.7%) | — | replicated signal at a locus without mapping risk |
| ZYG11A | kidney: epithelial cell of proximal tubule | 5 / 26 | 656/2,606 (25.2%) | 0.7 | replicated signal at a locus without mapping risk |
| EFCAB8 | eye: keratocyte | 6 / 19 | 1,569/6,549 (24.0%) | — | replicated signal at a locus without mapping risk |
| TEX101 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 43 | 272/1,209 (22.5%) | 0.31 | replicated signal at a locus without mapping risk |
| C5orf47 | eye: retinal pigment epithelial cell | 3 / 23 | 364/1,818 (20.0%) | — | replicated signal at a locus without mapping risk |
| PKD1L3 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 10 / 145 | 196/1,047 (18.7%) | 0.16 | replicated signal at a locus without mapping risk |
| ADAM29 | brain: dopaminergic neuron | 11 / 152 | 354/1,912 (18.5%) | 0.1 | replicated signal at a locus without mapping risk |
| FAM217A | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 6 / 181 | 247/1,362 (18.1%) | 0.062 | replicated signal at a locus without mapping risk |
| HORMAD2 | liver: hepatocyte | 2 / 4 | 871/5,256 (16.6%) | 0.61 | replicated signal at a locus without mapping risk |
| WFDC6 | lung: multiciliated columnar cell of tracheobronchial tree | 1 / 2 | 386/2,359 (16.4%) | 1.3 | signal confirmed by bulk RNA (exon inside another same-strand gene) |
| ARL5C | eye: diffuse bipolar 2 cell | 4 / 60 | 461/2,849 (16.2%) | — | replicated signal at a locus without mapping risk |
| SPINT4 | eye: retinal rod cell | 4 / 61 | 959/6,080 (15.8%) | — | replicated signal at a locus without mapping risk |
| TMPRSS12 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 18 | 192/1,258 (15.3%) | 0.026 | replicated signal at a locus without mapping risk |
| LY6K | eye: conjunctival epithelial cell | 3 / 9 | 184/1,246 (14.8%) | — | replicated signal at a locus without mapping risk |
| MROH2A | liver: hepatocyte | 2 / 4 | 761/5,256 (14.5%) | 0.3 | replicated signal at a locus without mapping risk |
| PLCZ1 | brain: L6b glutamatergic cortical neuron | 4 / 18 | 356/2,487 (14.3%) | 0 | replicated signal at a locus without mapping risk |
| CCDC172 | eye: retinal rod cell | 4 / 34 | 3,777/27,262 (13.9%) | — | replicated signal at a locus without mapping risk |
| XAGE2 | lung: pulmonary alveolar type 1 cell | 8 / 42 | 1,254/10,205 (12.3%) | 2.2 | replicated signal at a locus without mapping risk |
| SMIM23 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 24 | 145/1,190 (12.2%) | 0 | replicated signal at a locus without mapping risk |
| OIP5 | small intestine: transit amplifying cell of small intestine | 1 / 1 | 590/4,928 (12.0%) | 1.8 | signal confirmed by bulk RNA |
| DDX43 | eye: endothelial cell of venule | 7 / 19 | 133/1,148 (11.6%) | — | replicated signal at a locus without mapping risk |
| MMP20 | eye: OFF midget ganglion cell | 4 / 30 | 622/5,400 (11.5%) | — | replicated signal at a locus without mapping risk |
| SPACA3 | pancreas: pancreatic acinar cell | 2 / 2 | 540/4,757 (11.4%) | 3.1 | replicated signal at a locus without mapping risk |
| SAXO1 | brain: choroid plexus epithelial cell | 16 / 194 | 697/6,146 (11.3%) | 0.05 | replicated signal at a locus without mapping risk |
| PIWIL2 | brain: dopaminergic neuron | 15 / 146 | 649/5,945 (10.9%) | 1.7 | replicated signal at a locus without mapping risk |
| TDRD15 | liver: hepatocyte | 2 / 4 | 539/5,256 (10.3%) | 0.1 | replicated signal at a locus without mapping risk |
| OR51D1 | brain: pvalb GABAergic cortical interneuron | 2 / 25 | 217/2,146 (10.1%) | 0 | replicated signal at a locus without mapping risk |
| MDFIC2 | brain: oligodendrocyte | 5 / 118 | 157/1,559 (10.1%) | 0.6 | replicated signal at a locus without mapping risk |
| DGAT2L6 | eye: retinal cone cell | 3 / 17 | 190/1,958 (9.7%) | — | replicated signal at a locus without mapping risk |
| C5orf52 | colon: intestinal crypt stem cell of colon | 1 / 1 | 122/1,305 (9.3%) | 1.5 | signal confirmed by bulk RNA |
| PLAC1 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 54 | 1,134/12,360 (9.2%) | 0.0051 | replicated signal at a locus without mapping risk |
| OOSP3 | brain: glutamatergic neuron | 4 / 57 | 238/2,606 (9.1%) | 0.21 | replicated signal at a locus without mapping risk |
| MYO1H | eye: S cone cell | 6 / 35 | 100/1,109 (9.0%) | — | replicated signal at a locus without mapping risk |
| ERICH6 | brain: glutamatergic neuron | 10 / 103 | 583/6,853 (8.5%) | 0.13 | replicated signal at a locus without mapping risk |
| SCGB1D1 | eye: acinar cell | 1 / 3 | 100/1,178 (8.5%) | — | detection confirmed by a published gene-specific assay |
| TOPAZ1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 8 / 44 | 105/1,258 (8.3%) | 0.046 | replicated signal at a locus without mapping risk |
| SMC1B | brain: glutamatergic neuron | 9 / 81 | 548/6,853 (8.0%) | 0.12 | replicated signal at a locus without mapping risk |
| DUXB | brain: glutamatergic neuron | 3 / 23 | 510/6,853 (7.4%) | 0.056 | replicated signal at a locus without mapping risk |
| SLCO6A1 | eye: Mueller cell | 5 / 39 | 99/1,579 (6.3%) | — | replicated signal at a locus without mapping risk |
| COX7B2 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 6 / 16 | 941/15,320 (6.1%) | 0.043 | replicated signal at a locus without mapping risk |
| LACTBL1 | eye: flat midget bipolar cell | 2 / 36 | 87/1,437 (6.1%) | — | replicated signal at a locus without mapping risk |
| SELENOV | eye: trabecular meshwork cell | 3 / 5 | 255/4,340 (5.9%) | — | replicated signal at a locus without mapping risk |
| CCDC83 | eye: endothelial cell of venule | 2 / 10 | 62/1,148 (5.4%) | — | replicated signal at a locus without mapping risk |
| GRXCR1 | eye: OFF midget ganglion cell | 2 / 14 | 578/10,959 (5.3%) | — | replicated signal at a locus without mapping risk |
| WNT8A | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 5 / 70 | 783/15,132 (5.2%) | 0.75 | replicated signal at a locus without mapping risk |
| SPATA45 | brain: lamp5 GABAergic cortical interneuron | 2 / 5 | 80/1,548 (5.2%) | 0.33 | replicated signal at a locus without mapping risk |
| BOLL | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 18 | 70/1,362 (5.1%) | 0.017 | replicated signal at a locus without mapping risk |
| SPMIP11 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 4 | 65/1,362 (4.8%) | 0 | replicated signal at a locus without mapping risk |
| SPMIP3 | eye: endothelial cell of venule | 2 / 5 | 49/1,032 (4.7%) | — | replicated signal at a locus without mapping risk |
| C9orf57 | brain: oligodendrocyte | 3 / 21 | 74/1,608 (4.6%) | 0 | replicated signal at a locus without mapping risk |
| AADACL4 | eye: retinal pigment epithelial cell | 2 / 8 | 126/2,751 (4.6%) | — | replicated signal at a locus without mapping risk |
| ASZ1 | brain: oligodendrocyte | 2 / 53 | 50/1,141 (4.4%) | 0 | replicated signal at a locus without mapping risk |
| LRRIQ4 | brain: astrocyte | 2 / 2 | 65/1,508 (4.3%) | 0.052 | replicated signal at a locus without mapping risk |
| GARIN1B | brain: dopaminergic neuron | 3 / 8 | 82/1,912 (4.3%) | 0.73 | replicated signal at a locus without mapping risk |
| C16orf96 | eye: Mueller cell | 8 / 25 | 57/1,361 (4.2%) | — | replicated signal at a locus without mapping risk |
| LRRC72 | brain: GABAergic neuron | 3 / 8 | 46/1,174 (3.9%) | 0.04 | replicated signal at a locus without mapping risk |
| ASB17 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 10 | 582/15,132 (3.8%) | 0 | replicated signal at a locus without mapping risk |
| OTOL1 | brain: corticothalamic-projecting glutamatergic cortical neuron | 2 / 4 | 244/6,470 (3.8%) | 0.097 | replicated signal at a locus without mapping risk |
| SLC25A31 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 5 / 14 | 531/15,320 (3.5%) | 0.12 | replicated signal at a locus without mapping risk |
| MARCOL | brain: astrocyte | 2 / 39 | 64/1,951 (3.3%) | 0.089 | replicated signal at a locus without mapping risk |
| MSLNL | eye: anterior lens cell | 2 / 8 | 40/1,239 (3.2%) | — | replicated signal at a locus without mapping risk |
| ERICH6B | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 14 | 41/1,285 (3.2%) | 0.1 | replicated signal at a locus without mapping risk |
| OTOG | eye: Mueller cell | 4 / 20 | 71/2,328 (3.0%) | — | replicated signal at a locus without mapping risk |
| KLHL10 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 24 | 41/1,362 (3.0%) | 0.72 | replicated signal at a locus without mapping risk |
| MROH5 | eye: pigmented epithelial cell | 2 / 4 | 42/1,544 (2.7%) | — | replicated signal at a locus without mapping risk |
| PXT1 | brain: astrocyte | 2 / 3 | 27/1,029 (2.6%) | 0.4 | replicated signal at a locus without mapping risk |
| CBY3 | brain: dopaminergic neuron | 5 / 17 | 50/1,912 (2.6%) | 1.4 | replicated signal at a locus without mapping risk |
| SPMAP2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 6 | 27/1,034 (2.6%) | 0.41 | replicated signal at a locus without mapping risk |
| C1orf185 | brain: astrocyte | 5 / 9 | 39/1,508 (2.6%) | 0 | replicated signal at a locus without mapping risk |
| OR2T6 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 12 | 396/15,320 (2.6%) | 0.14 | replicated signal at a locus without mapping risk |
| MANSC4 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 34 | 35/1,362 (2.6%) | 2 | signal confirmed by bulk RNA (exon inside another same-strand gene) |
| TEX36 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 2 / 3 | 28/1,093 (2.6%) | 0 | replicated signal at a locus without mapping risk |
| A3GALT2 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 6 | 27/1,093 (2.5%) | 0.43 | replicated signal at a locus without mapping risk |
| HHIPL2 | eye: astrocyte | 2 / 3 | 27/1,094 (2.5%) | — | replicated signal at a locus without mapping risk |
| HTR3C | lung: capillary endothelial cell | 2 / 3 | 43/1,801 (2.4%) | 1.9 | signal confirmed by bulk RNA (sequence shared with another locus) |
| FNDC7 | brain: astrocyte | 6 / 8 | 38/1,612 (2.4%) | 0.29 | replicated signal at a locus without mapping risk |
| NPBWR1 | brain: GABAergic neuron | 2 / 7 | 39/1,718 (2.3%) | 0.78 | replicated signal at a locus without mapping risk |
| CIMIP7 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 18 | 22/1,038 (2.1%) | 0 | replicated signal at a locus without mapping risk |
| TMCO5A | brain: oligodendrocyte | 2 / 22 | 48/2,281 (2.1%) | 0 | replicated signal at a locus without mapping risk |
| TP53TG3/B/C/D/E/F | eye: retinal ganglion cell | 2 / 2 | 101/4,961 (2.0%) | — | replicated signal at a locus without mapping risk |
| NLRP14 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 2 / 9 | 349/17,556 (2.0%) | 0.096 | replicated signal at a locus without mapping risk |
| OR2AJ1 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 12 | 300/15,246 (2.0%) | 0.073 | replicated signal at a locus without mapping risk |
| ZAN | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 9 | 20/1,093 (1.8%) | 0 | replicated signal at a locus without mapping risk |
| FRMD7 | brain: pvalb GABAergic cortical interneuron | 3 / 5 | 45/2,483 (1.8%) | 0.069 | replicated signal at a locus without mapping risk |
| SOX30 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 4 | 18/1,034 (1.7%) | 0.088 | replicated signal at a locus without mapping risk |
| UCP1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 4 | 21/1,209 (1.7%) | 0.14 | replicated signal at a locus without mapping risk |
| DCAF4L1 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 5 | 260/15,320 (1.7%) | 0.82 | replicated signal at a locus without mapping risk |
| ZNF648 | brain: lamp5 GABAergic cortical interneuron | 2 / 5 | 26/1,608 (1.6%) | 0.061 | replicated signal at a locus without mapping risk |
| MS4A5 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 7 | 33/2,072 (1.6%) | 0 | replicated signal at a locus without mapping risk |
| WDR87 | brain: dopaminergic neuron | 3 / 7 | 30/1,912 (1.6%) | 0.11 | replicated signal at a locus without mapping risk |
| ISM2 | brain: indirect pathway medium spiny neuron | 3 / 4 | 27/1,738 (1.6%) | 2.3 | replicated signal at a locus without mapping risk |
| ROPN1B | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 8 | 30/1,938 (1.5%) | 1.7 | signal confirmed by bulk RNA (sequence shared with another locus) |
| LRRC74A | heart: cardiac endothelial cell | 2 / 3 | 44/2,891 (1.5%) | 0 | replicated signal at a locus without mapping risk |
| C18orf63 | lung: pulmonary alveolar type 2 cell | 2 / 4 | 49/3,297 (1.5%) | 0.055 | replicated signal at a locus without mapping risk |
| PRAME | eye: melanocyte | 2 / 2 | 15/1,026 (1.5%) | — | replicated signal at a locus without mapping risk |
| MGAT4D | brain: astrocyte | 3 / 3 | 22/1,612 (1.4%) | 0.02 | replicated signal at a locus without mapping risk |
| C1orf94 | brain: oligodendrocyte | 2 / 3 | 16/1,198 (1.3%) | 0.58 | replicated signal at a locus without mapping risk |
| ACRV1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 3 | 18/1,362 (1.3%) | 1.6 | replicated signal at a locus without mapping risk |
| EXD1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 2 | 16/1,285 (1.2%) | 0.16 | replicated signal at a locus without mapping risk |
| AWAT1 | eye: retinal cone cell | 2 / 2 | 17/1,447 (1.2%) | — | replicated signal at a locus without mapping risk |
| MAGEA12 | brain: pericyte | 0 / 0 | 5/1,263 (0.4%) | 0.3 | detection confirmed by a published gene-specific assay (sequence shared with another locus | exon inside another same-strand gene) |

## Dismissed

| Protein | Organ: strongest population | Studies / donors with multi-molecule cells | Strongest donor: cells with ≥2 molecules | Bulk TPM reached in both HPA and GTEx | Basis |
|---|---|---:|---:|---:|---|
| LINC03051 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 16 / 375 | 1,041/1,045 (99.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000214265 | brain: indirect pathway medium spiny neuron | 5 / 27 | 1,970/1,979 (99.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000272162 | eye: retinal rod cell | 3 / 19 | 15,411/16,080 (95.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000260729 | brain: dopaminergic neuron | 8 / 308 | 1,782/1,912 (93.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000266302 | brain: glutamatergic neuron | 9 / 188 | 1,088/1,168 (93.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000269237 | brain: dopaminergic neuron | 12 / 322 | 1,759/1,912 (92.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000255872 | eye: retinal cone cell | 5 / 30 | 1,295/1,447 (89.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000268797 | brain: glutamatergic neuron | 11 / 318 | 2,252/2,606 (86.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000284461 | brain: glutamatergic neuron | 6 / 231 | 1,261/1,621 (77.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| C1QTNF3-AMACR | brain: oligodendrocyte | 12 / 334 | 892/1,198 (74.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000274322 | brain: glutamatergic neuron | 8 / 295 | 1,798/2,606 (69.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000276490 | small intestine: enterocyte of epithelium proper of ileum | 1 / 4 | 1,496/2,263 (66.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000284772 | brain: glutamatergic neuron | 7 / 105 | 1,055/1,621 (65.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000284341 | eye: Mueller cell | 1 / 4 | 739/1,140 (64.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000262633 | brain: glutamatergic neuron | 5 / 27 | 923/1,440 (64.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000285932 | small intestine: enterocyte of epithelium proper of ileum | 1 / 4 | 1,403/2,263 (62.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000266076 | brain: glutamatergic neuron | 10 / 125 | 994/1,621 (61.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273398 | eye: retinal rod cell | 4 / 14 | 9,754/16,080 (60.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000269755 | brain: dopaminergic neuron | 9 / 128 | 1,137/1,912 (59.5%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| SLCO1B3-SLCO1B7 | liver: hepatocyte | 1 / 2 | 3,206/5,393 (59.4%) | 0.5 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000266086 | brain: dopaminergic neuron | 10 / 142 | 1,132/1,912 (59.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000281593 | brain: choroid plexus epithelial cell | 10 / 286 | 3,625/6,146 (59.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000260851 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 7 / 205 | 1,214/2,072 (58.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000205236 | eye: retinal ganglion cell | 1 / 6 | 754/1,308 (57.6%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| IGF2BP3 | heart: regular atrial cardiac myocyte | 1 / 1 | 3,236/5,692 (56.9%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| TSNAX-DISC1 | eye: retinal rod cell | 5 / 18 | 6,017/10,611 (56.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000283977 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 8 / 305 | 630/1,176 (53.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000264324 | eye: retinal cone cell | 5 / 29 | 757/1,447 (52.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000289700 | eye: corneal epithelial cell | 1 / 2 | 771/1,474 (52.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000266202 | brain: glutamatergic neuron | 9 / 289 | 3,424/6,853 (50.0%) | 0.2 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285635 | brain: glutamatergic neuron | 5 / 278 | 3,367/6,853 (49.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| PSG8 | eye: rod bipolar cell | 8 / 103 | 2,791/6,071 (46.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000262302 | lung: pulmonary alveolar type 2 cell | 4 / 7 | 1,889/4,263 (44.3%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000270099 | brain: astrocyte | 7 / 268 | 929/2,129 (43.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| OR56A1 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 1 / 2 | 472/1,106 (42.7%) | 0 | one study only; no bulk RNA |
| HHLA1 | brain: dopaminergic neuron | 4 / 37 | 810/1,912 (42.4%) | 0.018 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000262304 | liver: hepatocyte | 1 / 3 | 2,207/5,256 (42.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000278817 | eye: retinal ganglion cell | 6 / 23 | 2,871/7,102 (40.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000276418 | eye: retinal rod cell | 1 / 4 | 9,403/23,420 (40.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000287725 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 86 | 453/1,176 (38.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| NKX2-6 | heart: regular atrial cardiac myocyte | 1 / 1 | 2,093/5,692 (36.8%) | 0 | one study only; no bulk RNA |
| DCST1 | eye: retinal cone cell | 5 / 73 | 783/2,183 (35.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000255292 | eye: corneal epithelial cell | 2 / 7 | 522/1,474 (35.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000287694 | brain: GABAergic neuron | 5 / 272 | 603/1,718 (35.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273088 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 78 | 468/1,362 (34.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ZNF670-ZNF695 | brain: dopaminergic neuron | 8 / 81 | 656/1,912 (34.3%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| POLR2J2-UPK3BL1 | brain: glutamatergic neuron | 5 / 27 | 551/1,621 (34.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| PLD5P1 | brain: GABAergic neuron | 8 / 186 | 707/2,086 (33.9%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000255330 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 9 / 122 | 394/1,176 (33.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267168 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 136 | 362/1,093 (33.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285329 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 13 / 369 | 440/1,362 (32.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000239920 | eye: Mueller cell | 10 / 228 | 812/2,556 (31.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000270149 | brain: glutamatergic neuron | 8 / 80 | 510/1,621 (31.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ARMH2 | lung: multiciliated columnar cell of tracheobronchial tree | 2 / 4 | 329/1,068 (30.8%) | 0.4 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000283228 | eye: amacrine cell | 1 / 4 | 555/1,805 (30.7%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000288000 | brain: glutamatergic neuron | 4 / 166 | 795/2,606 (30.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273734 | brain: glutamatergic neuron | 5 / 58 | 2,031/6,853 (29.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000279983 | brain: glutamatergic neuron | 6 / 127 | 572/1,955 (29.3%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| MROH7-TTC4 | eye: Mueller cell | 4 / 21 | 332/1,140 (29.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000282218 | lung: pulmonary alveolar type 2 cell | 3 / 8 | 1,239/4,263 (29.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000249967 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 135 | 395/1,362 (29.0%) | 0.1 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000269242 | heart: regular ventricular cardiac myocyte | 1 / 17 | 900/3,136 (28.7%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| TRIM59-IFT80 | brain: oligodendrocyte | 9 / 280 | 427/1,529 (27.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273291 | eye: Mueller cell | 4 / 18 | 397/1,427 (27.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000286221 | brain: glutamatergic neuron | 4 / 140 | 721/2,606 (27.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258677 | eye: amacrine cell | 2 / 7 | 498/1,805 (27.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267699 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 8 / 282 | 365/1,362 (26.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000269693 | brain: glutamatergic neuron | 8 / 108 | 434/1,621 (26.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000285991 | tongue: stratified squamous epithelial cell | 1 / 3 | 952/3,616 (26.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267618 | eye: Mueller cell | 2 / 6 | 517/2,017 (25.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| SYS1-DBNDD2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 8 / 285 | 298/1,176 (25.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285953 | brain: GABAergic neuron | 3 / 177 | 294/1,174 (25.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273025 | eye: retinal rod cell | 2 / 7 | 4,003/16,080 (24.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285444 | brain: indirect pathway medium spiny neuron | 4 / 39 | 430/1,738 (24.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000256566 | brain: dopaminergic neuron | 9 / 123 | 466/1,912 (24.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| FKBP6 | lung: multiciliated columnar cell of tracheobronchial tree | 2 / 2 | 258/1,068 (24.2%) | 0.1 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000254732 | eye: Mueller cell | 1 / 4 | 339/1,427 (23.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285117 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 132 | 323/1,362 (23.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285404 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 5 / 86 | 346/1,468 (23.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285713 | brain: oligodendrocyte | 1 / 4 | 33,058/140,968 (23.5%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| SATL1 | brain: L6b glutamatergic cortical neuron | 13 / 150 | 344/1,471 (23.4%) | 0.42 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000251357 | eye: retinal rod cell | 2 / 7 | 5,418/23,420 (23.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| CIMIP3 | eye: retinal cone cell | 3 / 26 | 231/1,028 (22.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000100101 | brain: dopaminergic neuron | 4 / 14 | 426/1,912 (22.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ROPN1 | eye: melanocyte | 2 / 4 | 226/1,026 (22.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000273167 | eye: Mueller cell | 1 / 4 | 430/2,017 (21.3%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| GSX2 | eye: retinal pigment epithelial cell | 2 / 4 | 386/1,818 (21.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000284779 | pancreas: type B pancreatic cell | 1 / 15 | 294/1,395 (21.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267385 | liver: hepatocyte | 2 / 3 | 1,102/5,393 (20.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| NPHP3-ACAD11 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 125 | 223/1,093 (20.4%) | 0.3 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000272305 | brain: glutamatergic neuron | 8 / 106 | 551/2,702 (20.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000279686 | lung: vein endothelial cell | 2 / 5 | 293/1,459 (20.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273259 | pancreas: pancreatic acinar cell | 2 / 2 | 237/1,187 (20.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| BLOC1S5-TXNDC5 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 10 / 234 | 233/1,176 (19.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000280893 | eye: amacrine cell | 1 / 4 | 350/1,805 (19.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000268870 | brain: dopaminergic neuron | 7 / 113 | 370/1,912 (19.4%) | 0.1 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000283515 | eye: amacrine cell | 1 / 4 | 293/1,529 (19.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| CPHXL | brain: glutamatergic neuron | 3 / 59 | 1,284/6,853 (18.7%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000285000 | brain: oligodendrocyte | 4 / 265 | 357/1,937 (18.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285064 | eye: retinal rod cell | 1 / 2 | 205/1,137 (18.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258052 | eye: Mueller cell | 1 / 4 | 357/2,017 (17.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000284686 | brain: Bergmann glial cell | 5 / 37 | 242/1,393 (17.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000283563 | eye: amacrine cell | 3 / 9 | 301/1,805 (16.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| OR5V1 | brain: ependymal cell | 1 / 3 | 171/1,029 (16.6%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000257062 | liver: hepatocyte | 2 / 4 | 893/5,393 (16.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000255339 | eye: amacrine cell | 1 / 4 | 294/1,805 (16.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| DCAF8L2 | eye: invaginating midget bipolar cell | 1 / 2 | 246/1,513 (16.3%) | — | one study only; no bulk RNA |
| ENSG00000159239 | brain: oligodendrocyte | 5 / 65 | 533/3,287 (16.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ATF7-NPFF | eye: Mueller cell | 1 / 4 | 327/2,017 (16.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000271949 | eye: amacrine cell | 3 / 8 | 290/1,805 (16.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| EFCAB3 | lung: pulmonary alveolar type 2 cell | 2 / 10 | 277/1,738 (15.9%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000285791 | brain: astrocyte | 5 / 142 | 203/1,297 (15.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| OR4C3 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 2 | 247/1,582 (15.6%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000254536 | small intestine: enterocyte of epithelium proper of ileum | 1 / 4 | 351/2,263 (15.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000167807 | brain: glutamatergic neuron | 2 / 6 | 299/1,955 (15.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| RGPD3 | eye: diffuse bipolar 2 cell | 3 / 17 | 166/1,102 (15.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| TSGA10IP | urinary bladder: umbrella cell of urothelium | 1 / 1 | 219/1,489 (14.7%) | 0 | one study only; no bulk RNA |
| ENSG00000285942 | brain: Bergmann glial cell | 3 / 42 | 203/1,393 (14.6%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000281938 | eye: amacrine cell | 1 / 4 | 257/1,805 (14.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285777 | brain: ependymal cell | 1 / 3 | 500/3,573 (14.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000259753 | heart: regular atrial cardiac myocyte | 4 / 12 | 796/5,692 (14.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285558 | eye: non-pigmented ciliary epithelial cell | 2 / 8 | 155/1,110 (14.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285542 | brain: glutamatergic neuron | 4 / 106 | 955/6,853 (13.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| PSG11 | eye: rod bipolar cell | 3 / 44 | 816/6,071 (13.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000285471 | eye: astrocyte | 1 / 5 | 147/1,094 (13.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| PRDM7 | eye: melanocyte | 4 / 6 | 179/1,389 (12.9%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000255835 | brain: glutamatergic neuron | 8 / 87 | 200/1,621 (12.3%) | 0.3 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| KBTBD11-OT1 | eye: Mueller cell | 2 / 6 | 248/2,017 (12.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258465 | brain: glutamatergic neuron | 8 / 130 | 251/2,043 (12.3%) | 0.7 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| POTEC | eye: corneal endothelial cell | 2 / 2 | 150/1,243 (12.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| RGPD4 | eye: diffuse bipolar 2 cell | 3 / 16 | 337/2,849 (11.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000264058 | eye: retinal bipolar neuron | 1 / 4 | 592/5,162 (11.5%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ARHGAP19-SLIT1 | brain: dopaminergic neuron | 1 / 2 | 216/1,912 (11.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| SPECC1L-ADORA2A | eye: Mueller cell | 1 / 4 | 225/2,017 (11.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000277611 | brain: oligodendrocyte | 7 / 118 | 357/3,210 (11.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258539 | eye: retinal rod cell | 1 / 4 | 1,785/16,080 (11.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| KRTAP10-4 | eye: retinal rod cell | 1 / 10 | 678/6,161 (11.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000260272 | eye: retinal bipolar neuron | 1 / 4 | 285/2,655 (10.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| WFDC8 | eye: flat midget bipolar cell | 1 / 6 | 141/1,314 (10.7%) | — | one study only; no bulk RNA |
| WFDC13 | eye: retinal rod cell | 1 / 5 | 648/6,080 (10.7%) | — | one study only; no bulk RNA |
| ENSG00000279956 | eye: Mueller cell | 3 / 14 | 148/1,427 (10.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000269026 | brain: glutamatergic neuron | 7 / 119 | 710/6,853 (10.4%) | 0.1 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000260234 | eye: Mueller cell | 3 / 10 | 147/1,427 (10.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285854 | brain: GABAergic neuron | 3 / 65 | 120/1,174 (10.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000258864 | brain: dopaminergic neuron | 5 / 50 | 190/1,912 (9.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000269476 | brain: dopaminergic neuron | 4 / 20 | 188/1,912 (9.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| IZUMO2 | brain: dopaminergic neuron | 3 / 8 | 182/1,912 (9.5%) | 0.072 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000258790 | eye: Mueller cell | 1 / 4 | 191/2,017 (9.5%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000264187 | brain: dopaminergic neuron | 4 / 26 | 181/1,912 (9.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| URGCP-MRPS24 | eye: amacrine cell | 1 / 4 | 169/1,805 (9.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000268790 | eye: Mueller cell | 1 / 4 | 188/2,017 (9.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000281348 | eye: amacrine cell | 1 / 4 | 168/1,805 (9.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258466 | eye: Mueller cell | 1 / 4 | 105/1,140 (9.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285901 | brain: glutamatergic neuron | 3 / 70 | 239/2,606 (9.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000266997 | eye: retinal rod cell | 2 / 5 | 1,457/16,080 (9.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| TUBB8 | brain: glutamatergic neuron | 6 / 34 | 177/1,955 (9.1%) | 0.3 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000286231 | brain: glutamatergic neuron | 3 / 87 | 235/2,606 (9.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000263809 | eye: retinal rod cell | 1 / 4 | 1,119/12,412 (9.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| PPP3R2 | brain: dopaminergic neuron | 3 / 49 | 169/1,912 (8.8%) | 0.025 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000286001 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 72 | 137/1,582 (8.7%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000258461 | lung: pulmonary alveolar type 2 cell | 1 / 3 | 125/1,473 (8.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| UBE2F-SCLY | brain: GABAergic neuron | 9 / 93 | 99/1,174 (8.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000276087 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 66 | 114/1,362 (8.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000255663 | brain: dopaminergic neuron | 4 / 9 | 160/1,912 (8.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| POTEF | heart: ventricular cardiac muscle cell | 1 / 1 | 182/2,180 (8.3%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000256825 | eye: amacrine cell | 3 / 6 | 148/1,805 (8.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| SPESP1-NOX5 | brain: hippocampal pyramidal neuron | 1 / 6 | 273/3,388 (8.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| TECTB | eye: retinal pigment epithelial cell | 1 / 5 | 572/7,623 (7.5%) | — | one study only; no bulk RNA |
| ENSG00000251012 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 2 / 29 | 80/1,093 (7.3%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000288684 | eye: conjunctival epithelial cell | 1 / 3 | 127/1,740 (7.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| HNRNPUL2-BSCL2 | eye: amacrine cell | 1 / 4 | 130/1,805 (7.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000286070 | lung: pulmonary alveolar type 2 cell | 1 / 3 | 72/1,005 (7.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000284431 | eye: Mueller cell | 1 / 4 | 102/1,427 (7.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000141979 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 63 | 96/1,362 (7.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000272921 | brain: glutamatergic neuron | 3 / 23 | 180/2,606 (6.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ATP6V1G2-DDX39B | eye: retinal bipolar neuron | 1 / 4 | 182/2,655 (6.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ZNF280A | eye: conjunctival epithelial cell | 2 / 6 | 84/1,227 (6.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| FBXO47 | eye: keratocyte | 1 / 1 | 443/6,549 (6.8%) | — | one study only; no bulk RNA |
| TENT5D | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 24 | 87/1,318 (6.6%) | 0 | one study only; no bulk RNA |
| ENSG00000281883 | brain: oligodendrocyte | 7 / 180 | 78/1,198 (6.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267426 | brain: glutamatergic neuron | 2 / 5 | 104/1,621 (6.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| MYMK | eye: anterior lens cell | 1 / 1 | 79/1,239 (6.4%) | — | one study only; no bulk RNA |
| ENSG00000125695 | eye: amacrine cell | 1 / 4 | 115/1,805 (6.4%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258653 | eye: amacrine cell | 1 / 4 | 97/1,529 (6.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| POTEJ | heart: regular atrial cardiac myocyte | 2 / 6 | 88/1,400 (6.3%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ZC3H11B | brain: glutamatergic neuron | 2 / 5 | 101/1,621 (6.2%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000257390 | eye: amacrine cell | 1 / 4 | 112/1,805 (6.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| BORCS7-ASMT | eye: amacrine cell | 3 / 9 | 112/1,805 (6.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000213204 | brain: dopaminergic neuron | 1 / 3 | 117/1,912 (6.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000271793 | eye: amacrine cell | 3 / 8 | 110/1,805 (6.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000285589 | lung: pulmonary alveolar type 2 cell | 1 / 2 | 89/1,473 (6.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| RAD21L1 | eye: ON midget ganglion cell | 4 / 35 | 123/2,036 (6.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000288520 | eye: non-pigmented ciliary epithelial cell | 3 / 6 | 66/1,110 (5.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| MINDY4B | eye: non-pigmented ciliary epithelial cell | 2 / 38 | 191/3,229 (5.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| DAZL | eye: astrocyte | 1 / 4 | 163/2,782 (5.9%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000272741 | eye: amacrine cell | 1 / 4 | 88/1,529 (5.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000262560 | lung: lung secretory cell | 1 / 1 | 72/1,252 (5.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| TRIM77 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 12 | 68/1,190 (5.7%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000256349 | eye: retinal rod cell | 1 / 4 | 703/12,412 (5.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000288643 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 44 | 77/1,362 (5.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| KRTAP5-11 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 6 / 90 | 984/17,434 (5.6%) | 0.6 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| TMEM269 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 37 | 859/15,320 (5.6%) | 0.86 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| BTNL2 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 1 / 16 | 61/1,093 (5.6%) | 0.27 | one study only; no bulk RNA |
| ENSG00000259066 | brain: dopaminergic neuron | 5 / 45 | 106/1,912 (5.5%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000198211 | eye: amacrine cell | 1 / 4 | 82/1,529 (5.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000254553 | brain: dopaminergic neuron | 6 / 78 | 100/1,912 (5.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| POM121L2 | brain: dopaminergic neuron | 1 / 5 | 100/1,912 (5.2%) | 0.048 | one study only; no bulk RNA |
| ENSG00000267477 | eye: amacrine cell | 1 / 3 | 94/1,805 (5.2%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000269179 | eye: conjunctival epithelial cell | 2 / 3 | 64/1,249 (5.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| RGSL1 | brain: lamp5 GABAergic cortical interneuron | 4 / 26 | 102/2,004 (5.1%) | 0.012 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000258830 | eye: retinal bipolar neuron | 1 / 4 | 272/5,348 (5.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267881 | colon: intestinal crypt stem cell of colon | 1 / 1 | 66/1,305 (5.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000286239 | eye: melanocyte | 2 / 6 | 149/2,953 (5.0%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000286022 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 65 | 68/1,362 (5.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267179 | brain: glutamatergic neuron | 3 / 19 | 79/1,621 (4.9%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000249773 | brain: glutamatergic neuron | 4 / 35 | 79/1,621 (4.9%) | 0.5 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| DUX4 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 14 | 63/1,326 (4.8%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ZNF99 | kidney: kidney loop of Henle thick ascending limb epithelial cell | 4 / 8 | 54/1,165 (4.6%) | 0.1 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| OOSP1 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 16 | 687/15,320 (4.5%) | 0.11 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| VCX3B | eye: retinal ganglion cell | 2 / 3 | 315/7,102 (4.4%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000260371 | brain: glutamatergic neuron | 4 / 23 | 70/1,621 (4.3%) | 0.1 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| CTAG2 | lung: pulmonary alveolar type 2 cell | 2 / 3 | 166/3,865 (4.3%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000272442 | brain: glutamatergic neuron | 4 / 40 | 110/2,606 (4.2%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| FAM90A26 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 13 | 46/1,093 (4.2%) | 0.026 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| TMEM256-PLSCR3 | brain: glutamatergic neuron | 3 / 8 | 66/1,621 (4.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285880 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 29 | 47/1,285 (3.7%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ENSG00000267552 | brain: dopaminergic neuron | 3 / 4 | 68/1,912 (3.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| GAGE10 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 17 | 534/15,132 (3.5%) | 0.6 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| PSG3 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 21 | 38/1,093 (3.5%) | 0.013 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| PRR5-ARHGAP8 | eye: pigmented epithelial cell | 2 / 6 | 233/6,957 (3.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285816 | brain: oligodendrocyte | 4 / 77 | 59/1,766 (3.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| CAGE1 | brain: dopaminergic neuron | 2 / 3 | 63/1,912 (3.3%) | 0.013 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000243501 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 34 | 42/1,285 (3.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000285505 | brain: GABAergic neuron | 5 / 24 | 68/2,086 (3.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000261732 | eye: Mueller cell | 3 / 7 | 63/2,017 (3.1%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| OOSP4B | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 2 / 4 | 471/15,132 (3.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| OR2T2 | kidney: epithelial cell of proximal tubule | 3 / 4 | 94/3,189 (2.9%) | 0.033 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ZNF286A-TBC1D26 | brain: glutamatergic neuron | 3 / 22 | 32/1,168 (2.7%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| TAS2R1 | brain: direct pathway medium spiny neuron | 3 / 7 | 59/2,186 (2.7%) | 0.37 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000256206 | brain: glutamatergic neuron | 3 / 27 | 183/6,853 (2.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| MYH13 | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 5 / 14 | 40/1,525 (2.6%) | 0.061 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| OR4N2 | brain: GABAergic neuron | 5 / 15 | 28/1,090 (2.6%) | 0.012 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| EEF1AKMT4-ECE2 | brain: GABAergic neuron | 2 / 15 | 36/1,409 (2.6%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| RIPPLY3 | lung: pulmonary alveolar type 2 cell | 2 / 4 | 107/4,263 (2.5%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| OR2M2 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 15 | 422/18,270 (2.3%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000256407 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 7 | 27/1,285 (2.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000267173 | brain: astrocyte | 3 / 3 | 31/1,508 (2.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| TBCEL-TECTA | brain: GABAergic neuron | 3 / 12 | 24/1,174 (2.0%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| CYP26C1 | eye: Mueller cell | 2 / 2 | 26/1,434 (1.8%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| OR2F1 | eye: melanocyte | 3 / 3 | 19/1,067 (1.8%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| GSTT2 | brain: dopaminergic neuron | 2 / 3 | 26/1,497 (1.7%) | 0.11 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000268133 | brain: dopaminergic neuron | 3 / 3 | 33/1,912 (1.7%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000273269 | brain: astrocyte | 4 / 6 | 26/1,508 (1.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000284695 | eye: basal cell | 2 / 2 | 41/2,439 (1.7%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ZNF492 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 4 | 340/20,998 (1.6%) | 0.13 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000287839 | eye: endothelial cell of venule | 2 / 3 | 18/1,148 (1.6%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| C12orf50 | brain: L2/3-6 intratelencephalic projecting glutamatergic neuron | 2 / 2 | 17/1,093 (1.6%) | 0 | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000285446 | eye: conjunctival epithelial cell | 2 / 2 | 27/1,740 (1.6%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus | exon inside another same-strand gene) |
| ZNF705G | eye: Mueller cell | 2 / 4 | 18/1,221 (1.5%) | — | counts cannot be attributed to this locus; no bulk RNA (sequence shared with another locus) |
| ENSG00000258064 | brain: dopaminergic neuron | 2 / 4 | 28/1,912 (1.5%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000256966 | brain: glutamatergic neuron | 3 / 5 | 35/2,606 (1.3%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000248993 | brain: dopaminergic neuron | 2 / 3 | 24/2,117 (1.1%) | 0 | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000283580 | brain: glutamatergic neuron | 2 / 4 | 13/1,168 (1.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |
| ENSG00000272647 | brain: astrocyte | 2 / 2 | 16/1,508 (1.1%) | — | counts cannot be attributed to this locus; no bulk RNA (exon inside another same-strand gene) |

## Scope and limits

The panel is every named Vital or Protected cell population of 19
non-reproductive tissues (adipose tissue, adrenal gland, blood, brain, colon, esophagus, eye, heart, kidney, liver, lung, lymph node, pancreas, small intestine, spinal cord, stomach, tongue, ureter, urinary bladder), without immune, stromal, generic,
unknown or developmental labels. Abundant populations are included. Reproductive organs are
left out because germ-cell RNA contaminates their somatic cells. Bulk confirmation is
available only for organs sampled by both HPA and GTEx. Adult age, primary tissue and donor
identity come from metadata. Only UMI protocols
are graded. Census release 2025-11-08; filter `is_primary_data == True and disease == 'normal' and tissue_type == 'tissue'`.

- A gene-count matrix cannot check alignment, strand or exonic support; the mapping-risk flags are reference properties, not read-level checks.
- No empty-droplet or ambient-RNA correction is applied to the pinned counts.
- Donor overlap across collections has not been reviewed.
- Proteins that were Unsafe for another reason were not extracted; their other exclusions stand.

[Every finding](findings.csv) · [Donor populations with at least ten positive cells](donor-populations.csv.gz) ·
[Panel](panel.csv) · [Cutoff comparison](threshold-sensitivity.csv) ·
[Proteins selected for extraction](selected-proteins.csv.gz)

```sh
python scripts/build_rare_cell_evidence.py --counts tmp/normal-atlas/wholebody-2025-11-08 --prepare
python scripts/extract_rare_cell_panel_counts.py --counts tmp/normal-atlas/wholebody-2025-11-08
python scripts/build_gene_mapping_risk.py --proteins tmp/normal-atlas/counts-2025-11-08/selected-proteins.csv.gz
python scripts/build_rare_cell_evidence.py --counts tmp/normal-atlas/wholebody-2025-11-08
```
