# Rebuilt antigen catalog

Policy `oncoref.target-selection.2026-11-10.1`. **616 CTA, 244 Cancer Antigen only, 30 Reproductive Antigen only, 873 Other and 18,130 Unsafe.** These are outputs of screening rules, not validated tumor antigens.

[Filterable protein table](report.html) · [Every change from the 2026-11-09 draft](changes-from-draft.csv) · [Changes from release 1.36.0](category-changes.csv) · [Exact comparison scope](comparison.json)

## Curation funnel

![Curation funnel](funnel.png)

| Filter | Result | Proteins removed | Remaining |
|---|---|---:|---:|
| Accepted normal-tissue exclusion | Unsafe | 18,130 | 1,763 |
| Unresolved rare-cell RNA | Other:Review | 270 | 1,493 |
| Other normal evidence, coverage, identity or translation review | Other:Review | 548 | 945 |
| Nested or Unsafe-contained sequence | Other:Sequence | 40 | 905 |
| No qualifying expression, target maximum or Disposable cap | Other:Expression | 15 | 890 |

The 890 remaining Target candidates are 616 CTA, 244 Cancer Antigen only and 30 Reproductive Antigen only. [Funnel table](funnel.csv).

![Rebuilt counts](counts.png)

## The filters

| Filter | Rule | Failure |
|---|---|---|
| Normal expression | Existing bulk, cluster, immune, cornea, protein and reviewed-knowledge rules, unchanged | Unsafe for accepted exclusions; otherwise Other:Review |
| Rare normal-cell RNA | >= 1% of a donor's cells and >= 10 cells carry >= 2 molecules. The signal in >= 2 independent studies of the same organ. Donors of one study share its libraries and reference. | Established in a Vital population: Unsafe. Unresolved: Other:Review |
| Coverage, identity and translation | Unchanged | Other:Review |
| Sequence | Fewer than 95% of 9-mers inside Unsafe proteins, recomputed for the new Unsafe set | Other:Sequence |
| Expression bounds | Target maximum above every measured nonreproductive and Disposable value; Disposable at most 32 TPM and 32 nTPM | Other:Expression |
| Cancer evidence | Two independent patients each ≥1 TPM and above the normal maximum N, or one tumor ≥10 × max(N, 1) | No Cancer Antigen membership |
| Reproductive evidence | ≥1 native TPM/nTPM with read and cell support, or reliable protein | No Reproductive Antigen membership |

## Rare normal-cell RNA

A count matrix cannot show where a read aligned, and one mis-assigned read per cell from a paralog, pseudogene or overlapping gene looks like sparse expression. Cells with one molecule never change a category. A mis-assigned read from a paralog, pseudogene or overlapping gene produces exactly this pattern. Without replication, a signal counts only when >= 5% of a donor's cells carry multiple molecules. It can reach review, never Unsafe.

- **Established:** Replicated signal, no mapping risk, and confirmation. Bulk RNA >= 1 native TPM in the same organ in both the HPA and GTEx sample sets. Accepted protein evidence is assessed by its own rules.
- **Unresolved:** A replicated signal lacking a clean mapping or confirmation, or a frequent signal from one study.
- **Mapping risk:** Shared transcript sequence with another locus, or an exon inside another same-strand gene. A flagged gene cannot be established from counts.
- **Leaving review:** Read-level review settles an unresolved signal: unique alignment and exonic or splice-spanning reads establish it; their absence dismisses it with the evidence.

60 draft candidates move to Other:Review on this rule, including HHLA1, ADAM7, ENSG00000279983, CPHXL, OR4C3, DCAF8L2, ENSG00000255835, SMIM23, OR51D1, IZUMO2, OOSP3, PPP3R2. 1 becomes Unsafe: CBY3.

**34 proteins that were Unsafe only through the pooled 1% detection flag are candidates again.** Their donor-level signal is single-molecule, or multi-molecule in one study below 5% of cells. They are ADAD1, ADAM18, AMBN, C12orf50, C7orf33, CDCP2, CHCT1, CSAG1, CXorf58, CYP26C1, DPPA2, GDF9, GFY, ITIH6, MAGEA12, NLRP8, OR2F1, OR5V1, OR8D4, OTOG, PABPC3, PPP1R2B, PRDM7, RSPH6A, SLC25A2, SMIM9, SPAG11A, SPAG11B, SPATA3, TDRD15, TPTE, TRIM60, URGCP-MRPS24, VRTN. This is the rule working as written, not a demonstration that these proteins are absent from normal tissue. [Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md).

## Changes from the 2026-11-09 draft

| Draft | Now | Proteins |
|---|---|---:|
| CTA | Other:Review | 51 |
| CTA | Reproductive Antigen only | 7 |
| CTA | Unsafe | 1 |
| Cancer Antigen only | Other:Review | 18 |
| Other:Review | CTA | 5 |
| Other:Sequence | CTA | 1 |
| Unsafe | CTA | 31 |
| Unsafe | Cancer Antigen only | 2 |
| Unsafe | Other:Review | 95 |
| Unsafe | Other:Sequence | 1 |
| Unsafe | Reproductive Antigen only | 1 |

The tenfold single-tumor rule accounts for the moves to Reproductive Antigen only and for Cancer Antigens only that now lack a qualifying route. The target-maximum and Disposable-cap checks remove no candidate.

Cancer Antigens include 860 proteins; Reproductive Antigens include 646. Their intersection is CTA and their union contains 890 Target candidates.

Other contains 818 review, 40 sequence, 15 expression. A curated published target nomination is retained for 401 CTAs and 24 Cancer Antigens only. Broader paper coverage also includes biology, normal-expression findings and other evidence.

## Changes from release 1.36.0

| Released category | Rebuilt category | Proteins |
|---|---|---:|
| CTA | CTA | 530 |
| CTA | Other | 50 |
| CTA | Reproductive Antigen only | 26 |
| CTA | Unsafe | 1 |
| Other | CTA | 55 |
| Other | Cancer Antigen only | 242 |
| Other | Other | 727 |
| Other | Reproductive Antigen only | 3 |
| Unsafe | CTA | 31 |
| Unsafe | Cancer Antigen only | 2 |
| Unsafe | Other | 96 |
| Unsafe | Reproductive Antigen only | 1 |
| Unsafe | Unsafe | 18129 |

New CTAs: FKBP6, AMELY, H2BW2, ITIH6, PRSS54, RSPH6A, KRT37, SLC25A2, ANKRD60, SRMS, PRDM7, OR1E2, GALR3, VRTN, OCM2, TTLL8, CHCT1, CIB3, PABPC3, DPPA2, ADAD1, GDF9, SPAG11B, CXorf58, C12orf50, POTED, OR2B2, ADAM18, PROKR1, HOXD12, C7orf33, TRIM60, SPAG11A, AMBN, DEFB125, NLRP8, OR8D4, OR52D1, MAGEB17, SPATA31A6, OTOG, AADACL3, HUS1B, OR14A16, OR51C1, KRTAP4-1, CSAG1, SMIM9, GLYATL3, TXNDC8, KRTAP5-10, MAS1L, ARSH, KRTAP19-8, KRTAP20-4, OR2F1, MAGEA12, ZC3H11C, TRIM64C, ZNF705D, DEFB115, TDRD15, TRIM51G, ZNF816-ZNF321P, SUPT20HL1, PPP1R2B, RAB6D, PRR23C, WFDC6, OR5V1, KRTAP5-7, KRTAP5-9, ENSG00000255432, ENSG00000258691, ARMH2, GFY, ENSG00000267140, URGCP-MRPS24, TPTE, SPATA31A3, SPATA31A7, SPATA31A5, RPS4Y2, FAM240B, ENSG00000285566, POTEB/2.

Former CTAs now requiring review: WNT8A, ADAM7, SPMAP2, RGSL1, HHLA1, CCDC83, RBM46, BOLL, ASB17, NLRP14, POM121L2, KLHL10, IZUMO2, ERICH6, ERICH6B, TMCO5A, TAS2R1, COX7B2, TOPAZ1, TENT5D, OR4N2, DCAF4L1, CALHM3, RIPPLY3, C5orf47, SPATA45, SMIM23, LRRIQ4, PPP3R2, DCAF8L2, C1orf185, C9orf57, MGAT4D, MANSC4, LRRC72, ZNF99, TRIM77, GAGE10, PSG3, CIMIP7, PSG11, TAF11L12, TRIM75, TAF11L14, TMEM269, PRAMEF18, DUXB, TAF11L11, CPHXL, TAF11L13.

## Reproductive Antigens only

21 have potentially qualifying tumor measurements but unresolved patient independence; none reaches the single-sample threshold. The remaining 9 do not meet either cancer RNA rule.

| Protein | Tumor maximum TPM | Normal comparison TPM | Cancer assessment |
|---|---:|---:|---|
| CDY2A/B | 8.64 | 0 | Patient identity unresolved |
| SLC25A52 | 5.58 | 2.97 | Patient identity unresolved |
| ZSWIM2 | 2.74 | 1.52 | Patient identity unresolved |
| GARIN3 | 6.5 | 2.07 | Patient identity unresolved |
| CDY1 | 1.6 | 0 | Below cancer rule |
| CDY1B | 0.3 | 0.1 | Below cancer rule |
| SPATA3 | 24.1 | 3 | Below cancer rule |
| TMEM262 | 47.4 | 7.4 | Patient identity unresolved |
| MAGEB5 | 2.09 | 0.0392 | Patient identity unresolved |
| KRTAP24-1 | 38 | 6 | Patient identity unresolved |
| FAM47B | 2.25 | 0 | Patient identity unresolved |
| VN1R2 | 12.3 | 7.79 | Patient identity unresolved |
| OR8G1 | 2.69 | 2.5 | Below cancer rule |
| DEFB113 | 2.27 | 0 | Patient identity unresolved |
| ZNF705B | 3.78 | 0 | Patient identity unresolved |
| OR5H1 | 1.4 | 2.3 | Below cancer rule |
| FAM237A | 7.8 | 2.26 | Patient identity unresolved |
| ENSG00000237378 | 8.16 | 0.7 | Patient identity unresolved |
| ENSG00000241690 | 5.7 | 0.5 | Below cancer rule |
| H2AL1Q | 4.02 | 2.6 | Patient identity unresolved |
| P3R3URF | 16.5 | 3.4 | Patient identity unresolved |
| PTTG1IP2 | 3.23 | 0.3 | Patient identity unresolved |
| KDM4F | 5.49 | 0 | Patient identity unresolved |
| ENSG00000258417 | 6.76 | 1.2 | Patient identity unresolved |
| ENSG00000274944 | 13.5 | 2.8 | Patient identity unresolved |
| CCDC201 | 5.38 | 1.1 | Patient identity unresolved |
| HSFX3 | 4.88 | 1.8 | Patient identity unresolved |
| MSANTD5 | 3.43 | 2.2 | Below cancer rule |
| CIMIP3 | 1.2 | 0.2 | Below cancer rule |
| ENSG00000289809 | 0.0679 | 4.5 | Below cancer rule |

## Evidence and reproduction

The comparison uses source-native TPM and maxima of supported tissue/cell summaries. It does not correct cross-study batch effects. NUTM's stored three-sample matrix fails the nonnegative-TPM check and is listed in the [coverage ledger](tumor-source-coverage.csv).

[Full protein decisions](protein-decisions.csv.gz) · [All member genes](gene-decisions.csv.gz) · [Normal comparison values](normal-comparison.csv.gz) · [Tumor witnesses](tumor-witnesses.csv.gz) · [Paper claims](paper-attestations.json.gz) · [Input hashes](inputs.json).

The [pinned HPA download receipt](downloaded-hpa-sources.json) supplies exact URLs, paths and SHA-256 hashes for the newly retrieved native-TPM tables. With the pinned source files in the paths listed in `inputs.json`, run:

```sh
python scripts/build_rare_cell_evidence.py --counts tmp/normal-atlas/counts-2025-11-08
python scripts/rebuild_target_catalog.py
python scripts/render_target_catalog.py
python scripts/build_target_selection_report.py
```

The builder verifies the pinned tumor matrices, re-resolves source identities and assesses all proteins. Temporary tumor matrices stay outside the published report. The new catalog is available through `target_selection.results()`; the `cta_selection` API and existing sequence downloads retain the 1.36.0 snapshot.
