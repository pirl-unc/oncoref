# Rare normal-cell RNA by adult donor

**1,804 proteins were graded from 3,526,769 adult cells and nuclei of
955 donors in 49 collections, plus the directly analysed
Siletti cerebellar granule nuclei. 13 have established Vital-cell RNA and
282 need review.** Category consequences are in the
[catalog rebuild](../../target-selection-20261110/rebuild/index.md).

A count matrix cannot show where a read aligned. One mis-assigned read per cell, from a
paralog, a pseudogene or an overlapping gene, looks exactly like sparse expression. The
grade therefore ignores single-molecule cells and asks four separate questions.

| Question | Requirement |
|---|---|
| Is there a signal? | >= 1% of a donor's cells and >= 10 cells carry >= 2 molecules. |
| Is it replicated? | The signal in >= 2 independent studies of the same organ. Donors of one study share its libraries and reference. Without replication, a signal counts only when >= 5% of a donor's cells carry multiple molecules. It can reach review, never Unsafe. |
| Could the reads belong elsewhere? | Shared transcript sequence with another locus, or an exon inside another same-strand gene. A flagged gene cannot be established from counts. |
| Does another assay agree? | Bulk RNA >= 1 native TPM in the same organ in both the HPA and GTEx sample sets. Accepted protein evidence is assessed by its own rules. |

| Grade | Meaning | Vital population | Protected population |
|---|---|---|---|
| Established | Replicated signal, no mapping risk, and confirmation. | Unsafe | Other:Review |
| Unresolved | A replicated signal lacking a clean mapping or confirmation, or a frequent signal from one study. | Other:Review | Recorded |
| Below the screen | No replicated or frequent multi-molecule signal | No finding | No finding |

Below the screen is not a demonstration of absence or safety. Read-level review settles an unresolved signal: unique alignment and exonic or splice-spanning reads establish it; their absence dismisses it with the evidence.
Donors are never pooled, frequencies are never added across paralogs, and unmeasured
genes stay unassessed.

## Grades

| Grade | Tier | Catalog action | Proteins |
|---|---|---|---:|
| below_screen | — | none | 1,468 |
| established | vital | unsafe | 13 |
| unassessed | — | none | 36 |
| unresolved | protected | none | 5 |
| unresolved | vital | review | 282 |

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
| vital | 1% | 982 | 811 | 468 | 340 |
| vital | 2% | 870 | 701 | 409 | 278 |
| vital | 5% | 695 | 508 | 309 | 183 |
| protected | 1% | 438 | 239 | 123 | 65 |
| protected | 2% | 414 | 210 | 114 | 55 |
| protected | 5% | 266 | 141 | 77 | 37 |

## Established in a Vital population

| Protein | Organ: strongest population | Studies / donors with multi-molecule cells | Strongest donor: cells with ≥2 molecules | Bulk TPM reached in both HPA and GTEx | Bulk source |
|---|---|---:|---:|---:|---|
| SPAG9 | brain: dopaminergic neuron | 22 / 232 | 1,855/1,912 (97.0%) | 61 | gtex_tissue_rna:v11: Brain_Cerebellar_Hemisphere |
| TSGA10 | brain: dopaminergic neuron | 22 / 231 | 1,757/1,912 (91.9%) | 1.9 | gtex_tissue_rna:v11: Brain_Cerebellar_Hemisphere |
| SPAG1 | brain: choroid plexus epithelial cell | 22 / 227 | 4,740/6,146 (77.1%) | 7.9 | gtex_tissue_rna:v11: Brain_Cerebellar_Hemisphere |
| SYCP2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 20 / 219 | 944/1,362 (69.3%) | 4.3 | gtex_tissue_rna:v11: Brain_Cerebellum |
| MAGEC3 | brain: dopaminergic neuron | 13 / 126 | 731/1,912 (38.2%) | 1.5 | rna_brain_gtex:v23: hippocampal formation |
| SYCP3 | brain: astrocyte | 17 / 176 | 617/1,951 (31.6%) | 1.6 | gtex_tissue_rna:v11: Brain_Cerebellum |
| SPA17 | brain: dopaminergic neuron | 11 / 143 | 355/1,912 (18.6%) | 9.7 | rna_brain_hpa:v23: choroid plexus |
| XAGE2 | lung: pulmonary alveolar type 1 cell | 3 / 8 | 1,254/10,205 (12.3%) | 2.2 | rna_tissue_gtex:v23: lung |
| PIWIL2 | brain: dopaminergic neuron | 14 / 124 | 649/5,945 (10.9%) | 1.7 | gtex_tissue_rna:v11: Brain_Cerebellar_Hemisphere |
| DDX43 | lung: pulmonary alveolar type 1 cell | 2 / 4 | 106/2,752 (3.9%) | 1.5 | rna_tissue_hpa:v23: lung |
| CEP55 | lung: basal cell | 2 / 2 | 49/1,430 (3.4%) | 1.2 | rna_tissue_gtex:v23: lung |
| CBY3 | brain: dopaminergic neuron | 5 / 15 | 50/1,912 (2.6%) | 1.4 | rna_brain_gtex:v23: midbrain |
| ISM2 | brain: indirect pathway medium spiny neuron | 2 / 3 | 27/1,738 (1.6%) | 2.3 | rna_brain_hpa:v23: nucleus accumbens |

## Needing review

| Protein | Organ: strongest population | Studies / donors with multi-molecule cells | Strongest donor: cells with ≥2 molecules | Bulk TPM reached in both HPA and GTEx | Unresolved because |
|---|---|---:|---:|---:|---|
| LINC03051 | brain: astrocyte | 16 / 187 | 1,558/1,565 (99.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000214265 | brain: indirect pathway medium spiny neuron | 4 / 22 | 1,970/1,979 (99.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000260729 | brain: dopaminergic neuron | 8 / 132 | 1,782/1,912 (93.2%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000266302 | brain: glutamatergic neuron | 9 / 120 | 1,088/1,168 (93.2%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000269237 | brain: dopaminergic neuron | 11 / 144 | 1,759/1,912 (92.0%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| REDIC1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 19 / 189 | 1,722/1,938 (88.9%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000268797 | brain: glutamatergic neuron | 11 / 142 | 2,252/2,606 (86.4%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000284461 | brain: glutamatergic neuron | 6 / 108 | 1,261/1,621 (77.8%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| C1QTNF3-AMACR | brain: dopaminergic neuron | 12 / 155 | 1,353/1,912 (70.8%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000255872 | brain: dopaminergic neuron | 10 / 134 | 1,336/1,912 (69.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| GRXCR2 | heart: cardiac muscle myoblast | 4 / 20 | 4,504/6,478 (69.5%) | 0.2 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000274322 | brain: glutamatergic neuron | 8 / 124 | 1,798/2,606 (69.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| UBE2U | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 6 | 712/1,047 (68.0%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000284772 | brain: glutamatergic neuron | 7 / 102 | 1,055/1,621 (65.1%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000262633 | brain: glutamatergic neuron | 5 / 24 | 923/1,440 (64.1%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000266076 | brain: glutamatergic neuron | 10 / 123 | 994/1,621 (61.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| MTCP1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 13 / 145 | 823/1,362 (60.4%) | 1.9 | exon inside another same-strand gene |
| ENSG00000269755 | brain: dopaminergic neuron | 9 / 124 | 1,137/1,912 (59.5%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000266086 | brain: dopaminergic neuron | 10 / 131 | 1,132/1,912 (59.2%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000281593 | brain: choroid plexus epithelial cell | 10 / 134 | 3,625/6,146 (59.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000260851 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 7 / 118 | 1,214/2,072 (58.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| LIPI | brain: glutamatergic neuron | 12 / 139 | 3,975/6,853 (58.0%) | 0.43 | bulk RNA not confirmed in both HPA and GTEx |
| IGF2BP3 | heart: regular atrial cardiac myocyte | 1 / 1 | 3,236/5,692 (56.9%) | 0 | frequent multi-molecule signal in one study only | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000283977 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 8 / 133 | 630/1,176 (53.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000266202 | brain: glutamatergic neuron | 9 / 135 | 3,424/6,853 (50.0%) | 0.2 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285635 | brain: glutamatergic neuron | 4 / 82 | 3,367/6,853 (49.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000272162 | brain: glutamatergic neuron | 9 / 130 | 711/1,621 (43.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000270099 | brain: astrocyte | 7 / 124 | 929/2,129 (43.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| HHLA1 | brain: dopaminergic neuron | 4 / 30 | 810/1,912 (42.4%) | 0.018 | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OR56A1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 2 | 912/2,165 (42.1%) | 0 | frequent multi-molecule signal in one study only | bulk RNA not confirmed in both HPA and GTEx |
| ADAM7 | brain: dopaminergic neuron | 7 / 107 | 791/1,912 (41.4%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| MORC1 | brain: hippocampal pyramidal neuron | 2 / 6 | 585/1,451 (40.3%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| CSNK2A2IP | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 11 / 146 | 486/1,209 (40.2%) | 0.088 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000264324 | brain: dopaminergic neuron | 9 / 129 | 752/1,912 (39.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000287725 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 84 | 453/1,176 (38.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| NKX2-6 | heart: regular atrial cardiac myocyte | 1 / 1 | 2,093/5,692 (36.8%) | 0 | frequent multi-molecule signal in one study only | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000287694 | brain: GABAergic neuron | 5 / 87 | 603/1,718 (35.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000273088 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 76 | 468/1,362 (34.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ZNF670-ZNF695 | brain: dopaminergic neuron | 8 / 79 | 656/1,912 (34.3%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| PIWIL4 | brain: L6b glutamatergic cortical neuron | 16 / 189 | 503/1,471 (34.2%) | 5.6 | exon inside another same-strand gene |
| POLR2J2-UPK3BL1 | brain: glutamatergic neuron | 5 / 23 | 551/1,621 (34.0%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| PLD5P1 | brain: GABAergic neuron | 8 / 132 | 707/2,086 (33.9%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000255330 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 9 / 102 | 394/1,176 (33.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285329 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 12 / 182 | 440/1,362 (32.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| PSG8 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 10 / 135 | 619/1,938 (31.9%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000270149 | brain: glutamatergic neuron | 8 / 77 | 510/1,621 (31.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000255292 | brain: glutamatergic neuron | 10 / 124 | 446/1,440 (31.0%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000239920 | brain: astrocyte of the cerebral cortex | 12 / 157 | 406/1,312 (30.9%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| PIWIL1 | brain: near-projecting glutamatergic cortical neuron | 5 / 12 | 1,388/4,513 (30.8%) | 0.029 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000288000 | brain: glutamatergic neuron | 4 / 100 | 795/2,606 (30.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000273734 | brain: glutamatergic neuron | 5 / 54 | 2,031/6,853 (29.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000279983 | brain: glutamatergic neuron | 6 / 115 | 572/1,955 (29.3%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285932 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 105 | 398/1,362 (29.2%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000249967 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 96 | 395/1,362 (29.0%) | 0.1 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| METTL21EP | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 5 / 72 | 631/2,201 (28.7%) | — | bulk RNA not confirmed in both HPA and GTEx |
| SH2D4B | brain: glutamatergic neuron | 9 / 93 | 740/2,606 (28.4%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000286221 | brain: glutamatergic neuron | 4 / 101 | 721/2,606 (27.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| TDRD1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 7 / 38 | 527/1,938 (27.2%) | 0.17 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267699 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 8 / 127 | 365/1,362 (26.8%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000269693 | brain: glutamatergic neuron | 8 / 106 | 434/1,621 (26.8%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SYS1-DBNDD2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 7 / 100 | 298/1,176 (25.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285953 | brain: GABAergic neuron | 3 / 93 | 294/1,174 (25.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285444 | brain: indirect pathway medium spiny neuron | 4 / 37 | 430/1,738 (24.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000256566 | brain: dopaminergic neuron | 9 / 112 | 466/1,912 (24.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285117 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 66 | 323/1,362 (23.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285404 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 5 / 84 | 346/1,468 (23.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SATL1 | brain: L6b glutamatergic cortical neuron | 13 / 148 | 344/1,471 (23.4%) | 0.42 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| TEX101 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 22 | 272/1,209 (22.5%) | 0.31 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000205236 | brain: glutamatergic neuron | 5 / 19 | 687/3,063 (22.4%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000100101 | brain: dopaminergic neuron | 4 / 14 | 426/1,912 (22.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267168 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 76 | 427/1,938 (22.0%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000269242 | heart: pericyte | 1 / 5 | 221/1,024 (21.6%) | — | frequent multi-molecule signal in one study only | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000272305 | brain: glutamatergic neuron | 8 / 95 | 551/2,702 (20.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OOEP | lung: basal cell | 4 / 7 | 290/1,430 (20.3%) | 0.2 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000279686 | lung: vein endothelial cell | 2 / 5 | 293/1,459 (20.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| BLOC1S5-TXNDC5 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 10 / 132 | 233/1,176 (19.8%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000268870 | brain: dopaminergic neuron | 7 / 108 | 370/1,912 (19.4%) | 0.1 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000282218 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 73 | 370/1,938 (19.1%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000273025 | brain: glutamatergic neuron | 4 / 79 | 309/1,621 (19.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| MROH7-TTC4 | brain: astrocyte | 10 / 104 | 285/1,508 (18.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| MMP26 | brain: glutamatergic neuron | 4 / 53 | 509/2,702 (18.8%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| CPHXL | brain: glutamatergic neuron | 3 / 48 | 1,284/6,853 (18.7%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| PKD1L3 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 10 / 132 | 196/1,047 (18.7%) | 0.16 | bulk RNA not confirmed in both HPA and GTEx |
| ADAM29 | brain: dopaminergic neuron | 11 / 142 | 354/1,912 (18.5%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| FAM217A | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 6 / 107 | 247/1,362 (18.1%) | 0.062 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000284686 | brain: Bergmann glial cell | 5 / 32 | 242/1,393 (17.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000262302 | lung: basal cell | 2 / 4 | 235/1,430 (16.4%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| NPHP3-ACAD11 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 85 | 216/1,362 (15.9%) | 0.3 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285791 | brain: astrocyte | 5 / 101 | 203/1,297 (15.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OR4C3 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 1 | 247/1,582 (15.6%) | 0 | frequent multi-molecule signal in one study only | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000167807 | brain: glutamatergic neuron | 2 / 6 | 299/1,955 (15.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| TMPRSS12 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 12 | 192/1,258 (15.3%) | 0.026 | bulk RNA not confirmed in both HPA and GTEx |
| LINC03040 | lung: pulmonary alveolar type 1 cell | 2 / 8 | 1,532/10,205 (15.0%) | — | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285942 | brain: Bergmann glial cell | 3 / 41 | 203/1,393 (14.6%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| DCAF8L2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 46 | 280/1,938 (14.4%) | 0 | frequent multi-molecule signal in one study only | bulk RNA not confirmed in both HPA and GTEx |
| PLCZ1 | brain: L6b glutamatergic cortical neuron | 4 / 18 | 356/2,487 (14.3%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000259753 | heart: regular atrial cardiac myocyte | 4 / 12 | 796/5,692 (14.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285542 | brain: glutamatergic neuron | 4 / 74 | 955/6,853 (13.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000273291 | brain: glutamatergic neuron | 8 / 112 | 370/2,702 (13.7%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285000 | brain: glutamatergic neuron | 4 / 106 | 355/2,606 (13.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000255835 | brain: glutamatergic neuron | 8 / 86 | 200/1,621 (12.3%) | 0.3 | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000258465 | brain: glutamatergic neuron | 8 / 98 | 251/2,043 (12.3%) | 0.7 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SMIM23 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 12 | 145/1,190 (12.2%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| OIP5 | heart: cardiac muscle myoblast | 4 / 10 | 457/3,880 (11.8%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| ROPN1L | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 11 / 132 | 253/2,201 (11.5%) | 3.6 | exon inside another same-strand gene |
| SAXO1 | brain: choroid plexus epithelial cell | 15 / 144 | 697/6,146 (11.3%) | 0.05 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000273398 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 6 / 72 | 154/1,362 (11.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ARHGAP19-SLIT1 | brain: dopaminergic neuron | 1 / 2 | 216/1,912 (11.3%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285777 | brain: choroid plexus epithelial cell | 3 / 20 | 690/6,146 (11.2%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000283563 | brain: glutamatergic neuron | 5 / 98 | 125/1,122 (11.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000159239 | brain: dopaminergic neuron | 5 / 53 | 209/1,912 (10.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285713 | brain: Bergmann glial cell | 1 / 3 | 147/1,393 (10.6%) | — | frequent multi-molecule signal in one study only | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000269026 | brain: glutamatergic neuron | 7 / 86 | 710/6,853 (10.4%) | 0.1 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285854 | brain: GABAergic neuron | 3 / 62 | 120/1,174 (10.2%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OR51D1 | brain: pvalb GABAergic cortical interneuron | 2 / 25 | 217/2,146 (10.1%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| TRIM59-IFT80 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 6 / 72 | 195/1,938 (10.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000258864 | brain: dopaminergic neuron | 5 / 50 | 190/1,912 (9.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000269476 | brain: dopaminergic neuron | 4 / 19 | 188/1,912 (9.8%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000284341 | brain: astrocyte | 1 / 5 | 148/1,508 (9.8%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| IZUMO2 | brain: dopaminergic neuron | 3 / 8 | 182/1,912 (9.5%) | 0.072 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000264187 | brain: dopaminergic neuron | 4 / 26 | 181/1,912 (9.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285901 | brain: glutamatergic neuron | 3 / 67 | 239/2,606 (9.2%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OOSP3 | brain: glutamatergic neuron | 4 / 37 | 238/2,606 (9.1%) | 0.21 | bulk RNA not confirmed in both HPA and GTEx |
| TUBB8 | brain: glutamatergic neuron | 6 / 25 | 177/1,955 (9.1%) | 0.3 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| TSNAX-DISC1 | brain: dopaminergic neuron | 5 / 42 | 173/1,912 (9.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000286231 | brain: glutamatergic neuron | 3 / 85 | 235/2,606 (9.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| PPP3R2 | brain: dopaminergic neuron | 3 / 49 | 169/1,912 (8.8%) | 0.025 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000286001 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 71 | 137/1,582 (8.7%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ERICH6 | brain: glutamatergic neuron | 10 / 88 | 583/6,853 (8.5%) | 0.13 | bulk RNA not confirmed in both HPA and GTEx |
| EFCAB3 | brain: L6b glutamatergic cortical neuron | 10 / 48 | 125/1,471 (8.5%) | 0.017 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| SLCO1B3-SLCO1B7 | brain: choroid plexus epithelial cell | 1 / 2 | 521/6,146 (8.5%) | 0 | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| UBE2F-SCLY | brain: GABAergic neuron | 9 / 91 | 99/1,174 (8.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000276087 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 64 | 114/1,362 (8.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000255663 | brain: dopaminergic neuron | 3 / 6 | 160/1,912 (8.4%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| POTEF | heart: ventricular cardiac muscle cell | 1 / 1 | 182/2,180 (8.3%) | 0 | frequent multi-molecule signal in one study only | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| TOPAZ1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 7 / 30 | 105/1,258 (8.3%) | 0.046 | bulk RNA not confirmed in both HPA and GTEx |
| SPESP1-NOX5 | brain: hippocampal pyramidal neuron | 1 / 6 | 273/3,388 (8.1%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SMC1B | brain: glutamatergic neuron | 9 / 68 | 548/6,853 (8.0%) | 0.12 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000254536 | lung: basal cell | 1 / 3 | 113/1,430 (7.9%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| DUXB | brain: glutamatergic neuron | 3 / 19 | 510/6,853 (7.4%) | 0.056 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000289700 | lung: basal cell | 1 / 1 | 106/1,430 (7.4%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000255339 | brain: glutamatergic neuron | 6 / 63 | 117/1,621 (7.2%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000251012 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 27 | 137/1,938 (7.1%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000141979 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 62 | 96/1,362 (7.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| MDFIC2 | heart: mesothelial cell | 5 / 13 | 93/1,325 (7.0%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000272921 | brain: glutamatergic neuron | 3 / 13 | 180/2,606 (6.9%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| LY6K | lung: pulmonary alveolar type 1 cell | 3 / 3 | 185/2,752 (6.7%) | 0.6 | bulk RNA not confirmed in both HPA and GTEx |
| TENT5D | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 20 | 87/1,318 (6.6%) | 0 | frequent multi-molecule signal in one study only | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000279956 | brain: Bergmann glial cell | 3 / 18 | 264/4,013 (6.6%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267426 | brain: glutamatergic neuron | 2 / 5 | 104/1,621 (6.4%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| POTEJ | heart: regular atrial cardiac myocyte | 2 / 6 | 88/1,400 (6.3%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ZC3H11B | brain: glutamatergic neuron | 2 / 4 | 101/1,621 (6.2%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000276490 | lung: basal cell | 1 / 2 | 88/1,430 (6.2%) | — | frequent multi-molecule signal in one study only | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000213204 | brain: dopaminergic neuron | 1 / 3 | 117/1,912 (6.1%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000258677 | brain: dopaminergic neuron | 4 / 44 | 116/1,912 (6.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267618 | brain: dopaminergic neuron | 1 / 5 | 114/1,912 (6.0%) | 0.1 | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000271949 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 6 / 57 | 68/1,176 (5.8%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| HORMAD1 | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 7 / 35 | 127/2,201 (5.8%) | 0.054 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000262560 | lung: lung secretory cell | 1 / 1 | 72/1,252 (5.8%) | — | frequent multi-molecule signal in one study only | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| TRIM77 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 9 | 68/1,190 (5.7%) | 0 | frequent multi-molecule signal in one study only | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| PLAC1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 38 | 69/1,209 (5.7%) | 0.0051 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000278817 | lung: pulmonary alveolar epithelial cell | 2 / 4 | 109/1,917 (5.7%) | — | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000288643 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 43 | 77/1,362 (5.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000286070 | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 4 / 17 | 85/1,525 (5.6%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000259066 | brain: dopaminergic neuron | 5 / 44 | 106/1,912 (5.5%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000277611 | brain: GABAergic neuron | 6 / 34 | 64/1,174 (5.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| MYO1H | brain: dopaminergic neuron | 6 / 48 | 111/2,117 (5.2%) | 0.15 | bulk RNA not confirmed in both HPA and GTEx |
| POM121L2 | brain: dopaminergic neuron | 1 / 5 | 100/1,912 (5.2%) | 0.048 | frequent multi-molecule signal in one study only | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000254553 | brain: dopaminergic neuron | 5 / 61 | 100/1,912 (5.2%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| BTNL2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 1 / 15 | 71/1,362 (5.2%) | 0.27 | frequent multi-molecule signal in one study only | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000271793 | brain: glutamatergic neuron | 6 / 73 | 58/1,122 (5.2%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SPATA45 | brain: lamp5 GABAergic cortical interneuron | 2 / 5 | 80/1,548 (5.2%) | 0.33 | bulk RNA not confirmed in both HPA and GTEx |
| BOLL | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 16 | 70/1,362 (5.1%) | 0.017 | bulk RNA not confirmed in both HPA and GTEx |
| HORMAD2 | brain: astrocyte | 12 / 73 | 58/1,134 (5.1%) | 0.33 | bulk RNA not confirmed in both HPA and GTEx |
| RGSL1 | brain: lamp5 GABAergic cortical interneuron | 4 / 24 | 102/2,004 (5.1%) | 0.012 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| DGAT2L6 | brain: astrocyte | 3 / 33 | 53/1,042 (5.1%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| BORCS7-ASMT | brain: dopaminergic neuron | 6 / 59 | 96/1,912 (5.0%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000286022 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 62 | 68/1,362 (5.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000249773 | brain: glutamatergic neuron | 4 / 29 | 79/1,621 (4.9%) | 0.5 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267179 | brain: glutamatergic neuron | 3 / 19 | 79/1,621 (4.9%) | 0 | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SPMIP11 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 4 | 65/1,362 (4.8%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| DUX4 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 14 | 63/1,326 (4.8%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| COX7B2 | brain: lamp5 GABAergic cortical interneuron | 6 / 15 | 95/2,004 (4.7%) | 0.043 | bulk RNA not confirmed in both HPA and GTEx |
| KRTAP5-11 | brain: glutamatergic neuron | 6 / 51 | 126/2,702 (4.7%) | 0.6 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| RAD21L1 | brain: glutamatergic neuron | 4 / 22 | 118/2,606 (4.5%) | 0.063 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000260371 | brain: glutamatergic neuron | 4 / 22 | 70/1,621 (4.3%) | 0.1 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| LRRIQ4 | brain: astrocyte | 2 / 2 | 65/1,508 (4.3%) | 0.052 | bulk RNA not confirmed in both HPA and GTEx |
| ZNF99 | brain: dopaminergic neuron | 9 / 17 | 82/1,912 (4.3%) | 0.078 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| GARIN1B | brain: dopaminergic neuron | 3 / 7 | 82/1,912 (4.3%) | 0.73 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000272442 | brain: glutamatergic neuron | 4 / 39 | 110/2,606 (4.2%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| C5orf47 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 7 / 44 | 48/1,176 (4.1%) | 0.32 | bulk RNA not confirmed in both HPA and GTEx |
| TMEM256-PLSCR3 | brain: glutamatergic neuron | 3 / 7 | 66/1,621 (4.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| WNT8A | brain: pvalb GABAergic cortical interneuron | 5 / 54 | 103/2,532 (4.1%) | 0.75 | bulk RNA not confirmed in both HPA and GTEx |
| LDHC | brain: dopaminergic neuron | 8 / 11 | 76/1,912 (4.0%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| FAM90A26 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 11 | 43/1,086 (4.0%) | 0.026 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| LRRC72 | brain: GABAergic neuron | 3 / 8 | 46/1,174 (3.9%) | 0.04 | bulk RNA not confirmed in both HPA and GTEx |
| MROH2A | brain: near-projecting glutamatergic cortical neuron | 2 / 6 | 106/2,768 (3.8%) | 0.035 | bulk RNA not confirmed in both HPA and GTEx |
| KBTBD11-OT1 | brain: glutamatergic neuron | 3 / 10 | 62/1,621 (3.8%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OTOL1 | brain: corticothalamic-projecting glutamatergic cortical neuron | 2 / 4 | 244/6,470 (3.8%) | 0.097 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285880 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 29 | 47/1,285 (3.7%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285991 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 27 | 70/1,938 (3.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267552 | brain: dopaminergic neuron | 3 / 4 | 68/1,912 (3.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OOSP1 | brain: lamp5 GABAergic cortical interneuron | 3 / 9 | 71/2,004 (3.5%) | 0.11 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285558 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 23 | 41/1,176 (3.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| MARCOL | brain: astrocyte | 2 / 16 | 64/1,951 (3.3%) | 0.089 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000243501 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 34 | 42/1,285 (3.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285505 | brain: GABAergic neuron | 5 / 24 | 68/2,086 (3.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ASB17 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 6 | 63/1,938 (3.3%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| PRR5-ARHGAP8 | brain: astrocyte | 4 / 16 | 69/2,129 (3.2%) | 0.1 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ERICH6B | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 13 | 41/1,285 (3.2%) | 0.1 | bulk RNA not confirmed in both HPA and GTEx |
| RGPD3 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 6 / 32 | 40/1,298 (3.1%) | 0.31 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| KLHL10 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 18 | 41/1,362 (3.0%) | 0.72 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000256825 | brain: glutamatergic neuron | 3 / 35 | 58/1,955 (3.0%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| CCDC83 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 38 | 57/1,938 (2.9%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| GDPD4 | brain: L6b glutamatergic cortical neuron | 4 / 14 | 127/4,415 (2.9%) | 0.046 | bulk RNA not confirmed in both HPA and GTEx |
| PSG3 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 15 | 54/1,938 (2.8%) | 0.013 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| PBK | lung: respiratory tract suprabasal cell | 2 / 2 | 38/1,380 (2.8%) | 0.9 | bulk RNA not confirmed in both HPA and GTEx |
| ZNF286A-TBC1D26 | brain: glutamatergic neuron | 3 / 22 | 32/1,168 (2.7%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| C16orf96 | brain: glutamatergic neuron | 6 / 38 | 185/6,853 (2.7%) | 0.64 | bulk RNA not confirmed in both HPA and GTEx |
| TAS2R1 | brain: direct pathway medium spiny neuron | 3 / 7 | 59/2,186 (2.7%) | 0.37 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000256206 | brain: glutamatergic neuron | 3 / 27 | 183/6,853 (2.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| PSG11 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 12 | 51/1,938 (2.6%) | 0.018 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| RGS21 | brain: corticothalamic-projecting glutamatergic cortical neuron | 3 / 6 | 63/2,399 (2.6%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| PXT1 | brain: astrocyte | 2 / 3 | 27/1,029 (2.6%) | 0.4 | bulk RNA not confirmed in both HPA and GTEx |
| MYH13 | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 5 / 14 | 40/1,525 (2.6%) | 0.061 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| SPMAP2 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 6 | 27/1,034 (2.6%) | 0.41 | bulk RNA not confirmed in both HPA and GTEx |
| C1orf185 | brain: astrocyte | 5 / 9 | 39/1,508 (2.6%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| MANSC4 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 31 | 35/1,362 (2.6%) | 2 | exon inside another same-strand gene |
| OR4N2 | brain: GABAergic neuron | 4 / 14 | 28/1,090 (2.6%) | 0.012 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| EEF1AKMT4-ECE2 | brain: GABAergic neuron | 2 / 15 | 36/1,409 (2.6%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| EFCAB8 | brain: glutamatergic neuron | 5 / 14 | 65/2,606 (2.5%) | 0.054 | bulk RNA not confirmed in both HPA and GTEx |
| C9orf57 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 17 | 47/1,938 (2.4%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| TMEM269 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 25 | 47/1,938 (2.4%) | 0.86 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000251357 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 4 / 17 | 31/1,285 (2.4%) | 0.1 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| HTR3C | lung: capillary endothelial cell | 2 / 3 | 43/1,801 (2.4%) | 1.9 | sequence shared with another locus |
| FNDC7 | brain: astrocyte | 6 / 8 | 38/1,612 (2.4%) | 0.29 | bulk RNA not confirmed in both HPA and GTEx |
| ZYG11A | brain: choroid plexus epithelial cell | 3 / 5 | 143/6,146 (2.3%) | 0.48 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267385 | lung: pulmonary alveolar type 1 cell | 2 / 5 | 237/10,205 (2.3%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| NPBWR1 | brain: GABAergic neuron | 2 / 7 | 39/1,718 (2.3%) | 0.78 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000286239 | brain: glutamatergic neuron | 2 / 6 | 57/2,606 (2.2%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| GAGE10 | brain: indirect pathway medium spiny neuron | 4 / 13 | 38/1,738 (2.2%) | 0.6 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| CIMIP7 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 16 | 22/1,038 (2.1%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000256407 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 7 | 27/1,285 (2.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000267173 | brain: astrocyte | 3 / 3 | 31/1,508 (2.1%) | 0 | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| TBCEL-TECTA | brain: GABAergic neuron | 3 / 12 | 24/1,174 (2.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000258539 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 18 | 24/1,178 (2.0%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| A3GALT2 | brain: astrocyte | 3 / 4 | 31/1,565 (2.0%) | 0.43 | bulk RNA not confirmed in both HPA and GTEx |
| NLRP14 | brain: L5 extratelencephalic projecting glutamatergic cortical neuron | 2 / 8 | 30/1,525 (2.0%) | 0.096 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000264058 | brain: glutamatergic neuron | 3 / 9 | 51/2,606 (2.0%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OR2T6 | brain: pvalb GABAergic cortical interneuron | 3 / 8 | 61/3,159 (1.9%) | 0.14 | bulk RNA not confirmed in both HPA and GTEx |
| FRMD7 | brain: pvalb GABAergic cortical interneuron | 3 / 5 | 45/2,483 (1.8%) | 0.069 | bulk RNA not confirmed in both HPA and GTEx |
| RIPPLY3 | lung: pulmonary alveolar type 1 cell | 2 / 2 | 48/2,752 (1.7%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| SOX30 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 4 | 18/1,034 (1.7%) | 0.088 | bulk RNA not confirmed in both HPA and GTEx |
| UCP1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 3 | 21/1,209 (1.7%) | 0.14 | bulk RNA not confirmed in both HPA and GTEx |
| GSTT2 | brain: dopaminergic neuron | 2 / 3 | 26/1,497 (1.7%) | 0.11 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000268133 | brain: dopaminergic neuron | 3 / 3 | 33/1,912 (1.7%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000273269 | brain: astrocyte | 4 / 6 | 26/1,508 (1.7%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| DCAF4L1 | brain: glutamatergic neuron | 3 / 4 | 44/2,606 (1.7%) | 0.82 | bulk RNA not confirmed in both HPA and GTEx |
| SPMIP3 | brain: glutamatergic neuron | 4 / 14 | 43/2,606 (1.7%) | 0.49 | bulk RNA not confirmed in both HPA and GTEx |
| ZNF648 | brain: lamp5 GABAergic cortical interneuron | 2 / 5 | 26/1,608 (1.6%) | 0.061 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000285589 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 2 | 22/1,362 (1.6%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| TMCO5A | brain: glutamatergic neuron | 2 / 4 | 110/6,853 (1.6%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| MS4A5 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 6 | 33/2,072 (1.6%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| WDR87 | brain: dopaminergic neuron | 3 / 5 | 30/1,912 (1.6%) | 0.11 | bulk RNA not confirmed in both HPA and GTEx |
| OR2M2 | brain: glutamatergic neuron | 3 / 6 | 42/2,702 (1.6%) | 0 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |
| ROPN1B | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 7 | 30/1,938 (1.5%) | 1.7 | sequence shared with another locus |
| ENSG00000285446 | brain: glutamatergic neuron | 2 / 2 | 40/2,606 (1.5%) | — | sequence shared with another locus | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000258052 | brain: GABAergic neuron | 2 / 5 | 18/1,174 (1.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| LRRC74A | heart: cardiac endothelial cell | 2 / 3 | 44/2,891 (1.5%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| ZAN | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 3 / 8 | 20/1,362 (1.5%) | 0 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000258064 | brain: dopaminergic neuron | 2 / 3 | 28/1,912 (1.5%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| SLC25A31 | brain: lamp5 GABAergic cortical interneuron | 5 / 11 | 29/2,004 (1.4%) | 0.12 | bulk RNA not confirmed in both HPA and GTEx |
| ARL5C | brain: hippocampal pyramidal neuron | 2 / 3 | 33/2,325 (1.4%) | 0.15 | bulk RNA not confirmed in both HPA and GTEx |
| MGAT4D | brain: astrocyte | 3 / 3 | 22/1,612 (1.4%) | 0.02 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000256966 | brain: glutamatergic neuron | 3 / 5 | 35/2,606 (1.3%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| OR2AJ1 | brain: glutamatergic neuron | 2 / 2 | 36/2,702 (1.3%) | 0.073 | bulk RNA not confirmed in both HPA and GTEx |
| EXD1 | brain: L2/3 intratelencephalic projecting glutamatergic neuron | 2 / 2 | 16/1,285 (1.2%) | 0.16 | bulk RNA not confirmed in both HPA and GTEx |
| RBM46 | brain: choroid plexus epithelial cell | 2 / 2 | 74/6,146 (1.2%) | 0.035 | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000248993 | brain: dopaminergic neuron | 2 / 3 | 24/2,117 (1.1%) | 0 | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000283580 | brain: glutamatergic neuron | 2 / 4 | 13/1,168 (1.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ENSG00000272647 | brain: astrocyte | 2 / 2 | 16/1,508 (1.1%) | — | exon inside another same-strand gene | bulk RNA not confirmed in both HPA and GTEx |
| ZNF492 | brain: astrocyte | 2 / 2 | 35/3,456 (1.0%) | 0.13 | sequence shared with another locus | bulk RNA not confirmed in both HPA and GTEx |

## Scope and limits

The panel keeps the scope of the earlier pooled screen: brain, spinal cord, heart and
selected lung populations, without immune, stromal, generic, unknown or developmental
labels, and without populations above 10% of a tissue's sampled
cells. Adult age, primary tissue and donor identity come from metadata. Only UMI protocols
are graded. Census release 2025-11-08; filter `is_primary_data == True and disease == 'normal' and tissue_type == 'tissue'`.

- A gene-count matrix cannot check alignment, strand or exonic support; the mapping-risk flags are reference properties, not read-level checks.
- No empty-droplet or ambient-RNA correction is applied to the pinned counts.
- Donor overlap across collections has not been reviewed.
- Proteins that were Unsafe for another reason were not extracted; their other exclusions stand.

[Every finding](findings.csv) · [Donor populations with at least ten positive cells](donor-populations.csv.gz) ·
[Panel](panel.csv) · [Cutoff comparison](threshold-sensitivity.csv) ·
[Proteins selected for extraction](selected-proteins.csv.gz)

```sh
python scripts/extract_rare_cell_panel_counts.py --counts tmp/normal-atlas/counts-2025-11-08
python scripts/build_gene_mapping_risk.py --proteins tmp/normal-atlas/counts-2025-11-08/selected-proteins.csv.gz
python scripts/build_rare_cell_evidence.py --counts tmp/normal-atlas/counts-2025-11-08
```
