# Rebuilt antigen catalog

Policy `oncoref.target-selection.2026-11-18.2`. **203 CTA, 64 Cancer Antigen only, 211 Reproductive Antigen only, 1162 Other and 18,253 Unsafe.** These are outputs of screening rules, not validated tumor antigens.

[Filterable protein table](report.html) · [Changes from 1.39.0](changes-from-1.39.csv) · [Transition counts](changes-from-1.39-summary.csv) · [Changes from release 1.36.0](category-changes.csv) · [Exact comparison scope](comparison.json)

## End-to-end architecture and CTA tiers

![End-to-end curation architecture](summary/source-to-tiers-architecture.png)

![CTA tiers and host profiles](summary/cta-tiers-and-profiles.png)

## Curation funnel

![Curation funnel](funnel.png)

| Filter | Result | Proteins removed | Remaining |
|---|---|---:|---:|
| Bulk, cluster, immune, protein or reviewed normal exclusion | Unsafe | 18,128 | 1,765 |
| Established rare-cell RNA or published normal expression only | Unsafe | 125 | 1,640 |
| Coverage, identity, translation or unresolved evidence | Other:Review | 731 | 909 |
| Nested, Unsafe-contained or read-through sequence | Other:Sequence | 48 | 861 |
| No qualifying expression, target maximum or Disposable cap | Other:Expression | 320 | 541 |
| Home tissue no panel resolves, without reproductive evidence | Other:Unscreened tissue | 63 | 478 |

The 478 remaining Candidate antigen proteins are 203 CTA, 64 Cancer Antigen only and 211 Reproductive Antigen only. [Funnel table](funnel.csv) · [Curation summary with figures](summary.md).

![Rebuilt counts](counts.png)

## The filters

| Filter | Rule | Failure |
|---|---|---|
| Normal expression | Existing bulk, cluster, immune, cornea, protein and reviewed-knowledge rules, unchanged | Unsafe for accepted exclusions; otherwise Other:Review |
| Rare normal-cell RNA | >= 1% of a donor's cells and >= 10 cells carry >= 2 molecules. The signal in >= 2 independent studies of the same organ, which agree on a cell type or include >= 3 donors. Donors of one study share its libraries and reference; one donor in each of two studies, in different cell types, agrees on nothing but the organ. | Established: Unsafe. Unresolved: Other:Review. Dismissed: no effect |
| Published normal expression | Gene-specific report in normal adult Vital or Protected tissue, by a fixed rubric | Unsafe |
| Coverage, identity and translation | Unchanged | Other:Review |
| Sequence | Fewer than 95% of 9-mers inside Unsafe proteins; no unsupported nested protein or read-through locus | Other:Sequence |
| Expression bounds | Target maximum above every measured nonreproductive and Disposable value; Disposable at most 32 TPM and 32 nTPM | Other:Expression |
| Cancer evidence | One tumor sample ≥10 × max(N, Disposable, 1), N the highest non-reproductive normal value | No Cancer Antigen membership |
| Reproductive evidence | One reproductive sample or supported embryo-lineage mean ≥10 × max(N, Disposable, 1), on the same scale and reference as CA | No Reproductive Antigen membership |

## Rare normal-cell RNA

A count matrix cannot show where a read aligned, and one mis-assigned read per cell from a paralog, pseudogene or overlapping gene looks like sparse expression. Single-molecule detection alone does not establish expression. Ambient RNA and read assignment are possible confounders, not demonstrated explanations. Strong unresolved signals are held for review; mapping flags do not prove artifact.

- **Established, Unsafe:** A replicated signal at a locus without mapping risk, or any signal with bulk confirmation. Unsafe.
- **Unresolved, Other:Review:** A replicated signal with unresolved locus attribution, or a frequent single-study signal without independent confirmation. Other:Review pending adjudication.
- **Dismissed, no effect:** A weak or inconsistent signal below replication and frequency requirements. Recorded with its reason; no category effect.
- **Mapping risk:** Shared transcript sequence with another locus, or an exon inside another same-strand gene. Attribution needs read-level review; a flag does not prove artifact.
- **Bulk confirmation:** Bulk RNA >= 1 native TPM in the same organ in both the HPA and GTEx sample sets, a different assay with its own read assignment.
- **Literature:** Where a reviewed publication reports the gene in the same organ by a gene-specific assay, 1% detection in one donor is established, single molecules included: two assay types then agree. Published protein in a Vital or Protected tissue excludes through the knowledge layer without any counts.

80 draft non-Unsafe proteins become Unsafe, among them ADAM7, RBM46, C5orf47, SPINT4, MROH2A, CCDC172, SMIM23, MMP20, SPACA3, OR51D1, DGAT2L6, C5orf52, OOSP3, ERICH6. 4 candidates carry a dismissed signal, among them LRRIQ4, C12orf50, PRAME, EXD1.

**17 proteins that were Unsafe only through the pooled 1% detection flag are candidates again**: ADAD1, ADAM18, C12orf50, CHCT1, CSAG1, DPPA2, EXD1, GFY, NLRP8, PABPC3, PPP1R2B, SPAG11A, SPAG11B, SPATA3, TPTE, TRIM60, VRTN. Each was taken through the literature review below and none has a gene-specific report of expression in vital normal tissue. That is a search result, not proof of absence. [Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md).

## Literature review

912 proteins were reviewed: every protein that would be a Target candidate while its rare-cell signal is dismissed, or after release from the pooled 1% detection flag, and every Target candidate of the 2026-11-10 catalog. 19 are Unsafe on published evidence: ADAM30 (brain, R2); AMBN (small intestine, R2); ASCL4 (skin, R2); BECN2 (brain, R2); BTNL2 (small intestine, R2); CCN6 (cartilage, R2); DEFB106A/B (lung, R2); GALR3 (adrenal gland, R2); GDF9 (bone marrow, R2); ITIH6 (kidney, R2); MRGPRX1 (dorsal root ganglion, R2); MYH13 (skeletal muscle, R2); OBP2A (lung, R2); OTOG (inner ear, R3); SCGB1C1 (nasal mucosa, R2); SLC25A2 (liver, R2); SRY (brain, R2); SULT6B1 (kidney, R2); TCHHL1 (skin, R2). [Every protein, rubric, sources and quotes](literature.md).

## Changes from the 2026-11-09 draft

| Draft | Now | Proteins |
|---|---|---:|
| CTA | Cancer Antigen only | 49 |
| CTA | Other:Expression | 95 |
| CTA | Other:Review | 37 |
| CTA | Other:Sequence | 7 |
| CTA | Other:Unscreened tissue | 11 |
| CTA | Reproductive Antigen only | 194 |
| CTA | Unsafe | 54 |
| Cancer Antigen only | CTA | 1 |
| Cancer Antigen only | Other:Expression | 174 |
| Cancer Antigen only | Other:Review | 7 |
| Cancer Antigen only | Other:Sequence | 2 |
| Cancer Antigen only | Other:Unscreened tissue | 51 |
| Cancer Antigen only | Unsafe | 11 |
| Other:Review | CTA | 2 |
| Other:Review | Other:Expression | 20 |
| Other:Review | Reproductive Antigen only | 1 |
| Other:Review | Unsafe | 15 |
| Other:Sequence | Other:Expression | 1 |
| Other:Sequence | Other:Review | 1 |
| Reproductive Antigen only | Other:Expression | 11 |
| Reproductive Antigen only | Other:Review | 1 |
| Reproductive Antigen only | Other:Sequence | 1 |
| Unsafe | CTA | 9 |
| Unsafe | Cancer Antigen only | 1 |
| Unsafe | Other:Expression | 4 |
| Unsafe | Other:Review | 64 |
| Unsafe | Other:Unscreened tissue | 1 |
| Unsafe | Reproductive Antigen only | 7 |

The tenfold single-tumor rule accounts for the moves to Reproductive Antigen only and for Cancer Antigens only whose reproductive compartment does not reach the shared threshold. Known patient identity is a replication descriptor, not a single-sample admission gate.

Cancer Antigens include 267 proteins; Reproductive Antigens include 414. Their intersection is CTA and their union contains 478 Candidate antigen proteins.

Other contains 731 review, 48 sequence, 320 expression, 63 unscreened tissue. A curated published target nomination is retained for 143 CTAs and 33 Cancer Antigens only. Broader paper coverage also includes biology, normal-expression findings and other evidence.

## Changes from release 1.36.0

| Released category | Rebuilt category | Proteins |
|---|---|---:|
| CTA | CTA | 188 |
| CTA | Cancer Antigen only | 37 |
| CTA | Other | 131 |
| CTA | Reproductive Antigen only | 202 |
| CTA | Unsafe | 49 |
| Other | CTA | 6 |
| Other | Cancer Antigen only | 26 |
| Other | Other | 962 |
| Other | Reproductive Antigen only | 2 |
| Other | Unsafe | 31 |
| Unsafe | CTA | 9 |
| Unsafe | Cancer Antigen only | 1 |
| Unsafe | Other | 69 |
| Unsafe | Reproductive Antigen only | 7 |
| Unsafe | Unsafe | 18173 |

New CTAs: MAGEC2, ANKRD60, VRTN, TTLL8, DPPA2, SPAG11B, C12orf50, ADAM18, TRIM60, SPAG11A, SPATA31A6, CSAG1, TXNDC8, PPP1R2B, MAGEB6B.

Former CTAs now requiring review: BPY2/B/C, CDY2A/B, DEFB130A/B, DAZL, EVX1, TTLL2, RGSL1, TBX22, CSTL1, CST11, CTAG2, CDX4, HHLA1, USP26, GOLGA6D, POM121L2, WFDC8, IZUMO2, LRRC52, ZSWIM2, DCST1, CAGE1, SSX5, GOLGA6C, WFDC13, OR10G3, OR6F1, ZNF280A, TAS2R1, C8orf74, OR1N1, CDY1, CDY1B, XKR3, CCDC63, TENT5D, ARL13A, TSGA10IP, CCDC197, OR4N2, OR4M1, MAGEB10, SLC35G5, GOLGA8J, OR1R1P, CALHM3, RIPPLY3, POTEC, FAM9A, OR56A3, OR51M1, PRR23B, ZNF479, MAGEA11, DEFB132, ZACN, KRTAP21-2, PPP3R2, NANOS2, KRTAP24-1, DCAF8L2, SPATA31D4, RGPD4, ZNF705A, STPG3, OR8G1, RTP2, CT83, LIPJ, OR11A1, OR2J3, FBXO47, VCX3B, SMIM21, PRR23A, ZNF99, MINDY4B, TRIM77, C12orf71, CDRT15L2, GOLGA6B, GAGE10, ZNF705G, GAB4, PSG3, PSG7, GOLGA6L25, OR5H1, AK4P3, FAM237A, OR7E24, ENSG00000241690, PSG11, RBMY1D, EPPIN-WFDC6, TAF11L12, TRIM75, KIAA1210, TAF11L14, EXOC1L, OOSP4B, ZNF878, ENSG00000258417, GOLGA8T, DPEP2NB, KLF14, ATF7-NPFF, CTAGE15, ENSG00000273155, TMEM269, ENSG00000274944, XKR5, GSTT4, POTEB3, PRAMEF18, CCDC201, CTXND2, HSFX4, LITAFD, C19orf85, CIROP, HSFX3, TAF11L11, CPHXL, PMIS2, TAF11L13, OR8B3, OR8B2, ENSG00000284934, CIMIP3, ENSG00000289809.

## Reproductive Antigens only

No tumor sample reaches 10 × the protein's non-reproductive normal maximum.

| Protein | Tumor maximum TPM | Normal comparison TPM | Cancer assessment |
|---|---:|---:|---|
| DEFB104A/B | 5.94 | 5.61 | Below cancer rule |
| DEFB107A/B | 2.4 | 0.383 | Below cancer rule |
| GSC2 | 22.4 | 2.4 | Below cancer rule |
| FNDC8 | 11.4 | 26.9 | Below cancer rule |
| DEFB127 | 8.53 | 1.1 | Below cancer rule |
| CST9L | 38.6 | 15.3 | Below cancer rule |
| PRSS54 | 6.41 | 6.55 | Below cancer rule |
| ZNRF4 | 13.3 | 34.3 | Below cancer rule |
| SPAM1 | 3.06 | 3 | Below cancer rule |
| FSCN3 | 136 | 24.4 | Below cancer rule |
| HDGFL1 | 121 | 22.9 | Below cancer rule |
| SPACA1 | 23.6 | 11.7 | Below cancer rule |
| LYZL1 | 29 | 7.62 | Below cancer rule |
| PRM2 | 1.36e+03 | 3.05e+03 | Below cancer rule |
| CCDC70 | 21.2 | 26.8 | Below cancer rule |
| ATXN3L | 5.96 | 0.868 | Below cancer rule |
| IRGC | 16.5 | 48.1 | Below cancer rule |
| TCP11 | 56.3 | 57.8 | Below cancer rule |
| CST8 | 4.25 | 7.12 | Below cancer rule |
| DEFB129 | 30.7 | 14.1 | Below cancer rule |
| C20orf173 | 10.2 | 5.86 | Below cancer rule |
| RFPL3 | 29.2 | 3.2 | Below cancer rule |
| BMP15 | 3.33 | 1.78 | Below cancer rule |
| BPIFA3 | 12.3 | 13.5 | Below cancer rule |
| SCP2D1 | 7.07 | 18 | Below cancer rule |
| C4orf17 | 16.4 | 2.87 | Below cancer rule |
| CHCT1 | 41.8 | 28.2 | Below cancer rule |
| SLC25A52 | 5.58 | 2.97 | Below cancer rule |
| ZIM3 | 3.58 | 1.07 | Below cancer rule |
| CABS1 | 28.4 | 76.2 | Below cancer rule |
| AKAP4 | 16 | 84.9 | Below cancer rule |
| CPXCR1 | 3.33 | 1.19 | Below cancer rule |
| ACTL7B | 73.5 | 56.5 | Below cancer rule |
| C20orf144 | 109 | 14.4 | Below cancer rule |
| TKTL2 | 20.2 | 6.76 | Below cancer rule |
| LYZL2 | 25 | 14 | Below cancer rule |
| PABPC3 | 69.7 | 7.03 | Below cancer rule |
| GOLGA8F | 8.73 | 0.726 | Below cancer rule |
| ODF1 | 41.2 | 99.3 | Below cancer rule |
| CYLC2 | 11.7 | 13 | Below cancer rule |
| TRIM42 | 16.4 | 5.54 | Below cancer rule |
| CCDC116 | 48.6 | 8.37 | Below cancer rule |
| CCDC42 | 114 | 16 | Below cancer rule |
| CCDC27 | 13.2 | 8.41 | Below cancer rule |
| PDHA2 | 9.43 | 13.3 | Below cancer rule |
| SMCP | 77.3 | 186 | Below cancer rule |
| TEX55 | 18.5 | 21 | Below cancer rule |
| ADAD1 | 18.3 | 18.2 | Below cancer rule |
| FAM170A | 3.45 | 2.84 | Below cancer rule |
| PRKACG | 19.9 | 5.47 | Below cancer rule |
| PRSS37 | 44.1 | 7.14 | Below cancer rule |
| SSMEM1 | 7.4 | 5.94 | Below cancer rule |
| SPATA19 | 13.1 | 36.8 | Below cancer rule |
| C16orf78 | 5.59 | 16.5 | Below cancer rule |
| CCDC182 | 14.7 | 6.19 | Below cancer rule |
| SUN5 | 4.94 | 5.31 | Below cancer rule |
| TBC1D21 | 26.8 | 9.25 | Below cancer rule |
| GOT1L1 | 7.06 | 2.13 | Below cancer rule |
| ACTRT2 | 21.3 | 48.5 | Below cancer rule |
| GARIN3 | 3.07 | 11.3 | Below cancer rule |
| RBMXL2 | 47.9 | 7.08 | Below cancer rule |
| TPD52L3 | 42.2 | 14.7 | Below cancer rule |
| PGK2 | 18.1 | 50.8 | Below cancer rule |
| DDI1 | 30.8 | 6.01 | Below cancer rule |
| SPMIP9 | 71.3 | 31.9 | Below cancer rule |
| FAM170B | 10.6 | 3.73 | Below cancer rule |
| IQCF1 | 12.9 | 21.9 | Below cancer rule |
| RNASE11 | 20 | 10.6 | Below cancer rule |
| NLRP13 | 6.04 | 3.13 | Below cancer rule |
| SPATA3 | 24.1 | 11.4 | Below cancer rule |
| CBY2 | 24.5 | 28.7 | Below cancer rule |
| GOLGA6L2 | 16.1 | 9.09 | Below cancer rule |
| HIGD2B | 10.1 | 1.45 | Below cancer rule |
| UBQLN3 | 29.7 | 45.8 | Below cancer rule |
| PRM1 | 1.71e+03 | 2.86e+03 | Below cancer rule |
| RBMXL3 | 12.9 | 1.45 | Below cancer rule |
| TMEM270 | 29.6 | 21.1 | Below cancer rule |
| HMGB4 | 111 | 72.6 | Below cancer rule |
| CETN1 | 22.9 | 29.9 | Below cancer rule |
| FAM187B | 48.1 | 20.3 | Below cancer rule |
| TEX44 | 29.8 | 46.2 | Below cancer rule |
| DEFB114 | 7.64 | 0.484 | Below cancer rule |
| CAPZA3 | 129 | 62.9 | Below cancer rule |
| CIMAP1A | 15 | 1.95 | Below cancer rule |
| TNP2 | 9.04 | 4.38 | Below cancer rule |
| CCDC185 | 20.5 | 13.3 | Below cancer rule |
| H1-8 | 7.68 | 1.33 | Below cancer rule |
| EXD1 | 4.69 | 6.28 | Below cancer rule |
| SIGLECL1 | 16.1 | 3.07 | Below cancer rule |
| TMEM31 | 85.1 | 45 | Below cancer rule |
| DNAJB8 | 16.1 | 17.8 | Below cancer rule |
| RNF151 | 34.5 | 47.2 | Below cancer rule |
| NLRP8 | 3.76 | 0.4 | Below cancer rule |
| GARIN5B | 41.1 | 10.3 | Below cancer rule |
| WFDC9 | 4.39 | 4.07 | Below cancer rule |
| WFDC10A | 12.2 | 24 | Below cancer rule |
| DEFB123 | 50.5 | 26.4 | Below cancer rule |
| C3orf22 | 11.8 | 6.69 | Below cancer rule |
| SPEM1 | 9.16 | 35.6 | Below cancer rule |
| ACTL9 | 36.4 | 35.9 | Below cancer rule |
| ADIG | 33.8 | 11.4 | Below cancer rule |
| ZNF716 | 10.3 | 1.4 | Below cancer rule |
| CCDC87 | 23.3 | 4.54 | Below cancer rule |
| TMEM95 | 18.3 | 5.65 | Below cancer rule |
| CYLC1 | 2.57 | 13 | Below cancer rule |
| C10orf120 | 12.1 | 10.3 | Below cancer rule |
| GOLGA8G | 6.24 | 1.75 | Below cancer rule |
| TRIML1 | 30.1 | 3.03 | Below cancer rule |
| DEFB108B | 10.4 | 1.44 | Below cancer rule |
| IQCF2 | 20.6 | 22.4 | Below cancer rule |
| SPEM2 | 33.5 | 39.9 | Below cancer rule |
| ODF4 | 12.4 | 2.58 | Below cancer rule |
| CIMIP4 | 14.1 | 17.8 | Below cancer rule |
| SPPL2C | 17.7 | 19.2 | Below cancer rule |
| DEFB128 | 4.76 | 2.58 | Below cancer rule |
| TEX38 | 47.5 | 47 | Below cancer rule |
| PRR30 | 13.9 | 36.9 | Below cancer rule |
| SPATA31D3 | 12.6 | 1.73 | Below cancer rule |
| ACTL7A | 29.3 | 62.7 | Below cancer rule |
| TMEM262 | 47.4 | 7.4 | Below cancer rule |
| H2AP | 112 | 16.4 | Below cancer rule |
| PRAMEF10 | 142 | 17.4 | Below cancer rule |
| C2orf78 | 13.2 | 5.06 | Below cancer rule |
| ZCCHC13 | 32.6 | 6.39 | Below cancer rule |
| MAGEB5 | 2.09 | 0.486 | Below cancer rule |
| CATSPER4 | 8.09 | 2.45 | Below cancer rule |
| TMCO2 | 29 | 30.8 | Below cancer rule |
| DUSP21 | 41.6 | 4.9 | Below cancer rule |
| FAM47B | 2.25 | 2.71 | Below cancer rule |
| FSCB | 33.6 | 10.3 | Below cancer rule |
| ZAR1L | 12.1 | 1.43 | Below cancer rule |
| SPANXN3 | 20.3 | 17 | Below cancer rule |
| SPANXN4 | 33 | 4.87 | Below cancer rule |
| TBC1D28 | 4.07 | 2.84 | Below cancer rule |
| PFN3 | 19.8 | 7.58 | Below cancer rule |
| TOMM20L | 46.6 | 6.59 | Below cancer rule |
| SPMIP10 | 14.1 | 7.27 | Below cancer rule |
| CCER1 | 106 | 10.7 | Below cancer rule |
| TMEM239 | 19.3 | 6.51 | Below cancer rule |
| CCT8L2 | 30.8 | 4.34 | Below cancer rule |
| MAGEB3 | 7.55 | 0.732 | Below cancer rule |
| TEDDM1 | 15.2 | 3 | Below cancer rule |
| LELP1 | 119 | 212 | Below cancer rule |
| FAM24A | 4.13 | 2.83 | Below cancer rule |
| CXorf66 | 14.5 | 5.14 | Below cancer rule |
| C10orf62 | 27.4 | 87.4 | Below cancer rule |
| DEFB110 | 5.4 | 1.1 | Below cancer rule |
| GLT6D1 | 7.41 | 3.42 | Below cancer rule |
| TMEM225 | 18.3 | 2.82 | Below cancer rule |
| SPANXN5 | 7.93 | 18.8 | Below cancer rule |
| PRAMEF17 | 5.69 | 1.33 | Below cancer rule |
| PRAMEF19 | 9.62 | 1.16 | Below cancer rule |
| DEFB121 | 2.15 | 1.06 | Below cancer rule |
| OR2H1 | 3.45 | 2.77 | Below cancer rule |
| SPATA31A1 | 6.83 | 1.22 | Below cancer rule |
| PABPN1L | 58.5 | 25.9 | Below cancer rule |
| SPATA31F1 | 8.67 | 6.84 | Below cancer rule |
| RNASE13 | 26.6 | 8.13 | Below cancer rule |
| TSSK1B | 6.49 | 2.05 | Below cancer rule |
| TSGA13 | 5.95 | 1.77 | Below cancer rule |
| WEE2 | 10.7 | 2.23 | Below cancer rule |
| EFCAB9 | 8.56 | 3.73 | Below cancer rule |
| NOTO | 4.26 | 1.7 | Below cancer rule |
| DEFB113 | 2.27 | 0.539 | Below cancer rule |
| IQCF5 | 13 | 1.6 | Below cancer rule |
| SPATA31D1 | 25.5 | 4.14 | Below cancer rule |
| TCP11X2 | 15.4 | 1.91 | Below cancer rule |
| C16orf90 | 56.7 | 14.5 | Below cancer rule |
| POM121L12 | 7.05 | 4.15 | Below cancer rule |
| ZNF735 | 5.29 | 0.353 | Below cancer rule |
| GPX5 | 7.92 | 0.496 | Below cancer rule |
| KIF4B | 6.33 | 0.948 | Below cancer rule |
| TEX46 | 16 | 13.9 | Below cancer rule |
| BOD1L2 | 35.7 | 96.3 | Below cancer rule |
| TEX50 | 12.9 | 4.23 | Below cancer rule |
| C16orf82 | 26.2 | 94.7 | Below cancer rule |
| FREY1 | 34.9 | 9.69 | Below cancer rule |
| ENSG00000237378 | 6.52 | 1.26 | Below cancer rule |
| TEX51 | 3.36 | 7.54 | Below cancer rule |
| GOLGA6L24 | 2.27 | 0.36 | Below cancer rule |
| TSPY9 | 5.93 | 1.2 | Below cancer rule |
| MPC1L | 33.2 | 6.82 | Below cancer rule |
| PRR20G | 7.77 | 15.2 | Below cancer rule |
| PRAMEF11 | 7.24 | 0.91 | Below cancer rule |
| OR14A2 | 6.35 | 1.54 | Below cancer rule |
| H2AL1Q | 4.02 | 5.06 | Below cancer rule |
| P3R3URF | 16.5 | 29.6 | Below cancer rule |
| PTTG1IP2 | 3.23 | 2.23 | Below cancer rule |
| PRSS58 | 38.2 | 4.25 | Below cancer rule |
| RNASE12 | 13.9 | 2.17 | Below cancer rule |
| C20orf141 | 56.9 | 34.1 | Below cancer rule |
| DUXA | 38.2 | 9.24 | Below cancer rule |
| TLE7 | 1.34 | 0.971 | Below cancer rule |
| IGBP1C | 7.1 | 2.58 | Below cancer rule |
| TEX13A | 3.77 | 3.29 | Below cancer rule |
| SPANXN2 | 1.81 | 1.16 | Below cancer rule |
| CALR3 | 14.6 | 15.7 | Below cancer rule |
| LYPD4 | 73.7 | 22.1 | Below cancer rule |
| TPTE | 121 | 13 | Below cancer rule |
| SPMAP1 | 21.7 | 3.6 | Below cancer rule |
| LYZL6 | 11.9 | 8.06 | Below cancer rule |
| GOLGA6L6 | 5.37 | 0.747 | Below cancer rule |
| CT47C1 | 2.17 | 0.941 | Below cancer rule |
| GOLGA6L22 | 5.33 | 0.794 | Below cancer rule |
| PRAMEF13 | 4.53 | 0.289 | Below cancer rule |
| SPEM3 | 5.22 | 15.1 | Below cancer rule |
| C1orf202 | 12.1 | 1.64 | Below cancer rule |
| ENSG00000285566 | 8.23 | 0.545 | Below cancer rule |
| MSANTD5 | 3.43 | 3.83 | Below cancer rule |
| SPADH | 11.7 | 12 | Below cancer rule |
| ETDA/B | 26.2 | 2.93 | Below cancer rule |

## Evidence and reproduction

The comparison uses source-native TPM, individual GTEx v11 and HPA v23 tissue samples and supported tissue/cell summaries. It is not a maximum over every healthy donor. It does not correct cross-study batch effects. NUTM's stored three-sample matrix fails the nonnegative-TPM check and is listed in the [coverage ledger](tumor-source-coverage.csv).

[Full protein decisions](protein-decisions.csv.gz) · [All member genes](gene-decisions.csv.gz) · [HPA sample QC](normal-library-qc.csv) · [GTEx sample QC and donors](gtex-library-qc.csv.gz) · [Normal comparison values](normal-comparison.csv.gz) · [Tumor witnesses](tumor-witnesses.csv.gz) · [Reproductive observations](reproductive-comparison.csv.gz) · [Reproductive witnesses](reproductive-witnesses.csv.gz) · [Paper claims](paper-attestations.json.gz) · [Input hashes](inputs.json).

The [pinned HPA download receipt](downloaded-hpa-sources.json) supplies exact URLs, paths and SHA-256 hashes for the newly retrieved native-TPM tables. With the pinned source files in the paths listed in `inputs.json`, run:

```sh
python scripts/build_rare_cell_evidence.py --regrade-pinned
python scripts/rebuild_target_catalog.py
python scripts/render_target_catalog.py
python scripts/build_target_selection_report.py
```

The builder verifies the pinned tumor matrices, re-resolves source identities and assesses all proteins. Temporary tumor matrices stay outside the published report. The new catalog is available through `target_selection.results()`; the `cta_selection` API and existing sequence downloads retain the 1.36.0 snapshot.
