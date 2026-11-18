# Recovering evidence for former CTAs

Audit of released OncoRef 1.42.0. The catalog remains **290 CTAs**. Every scenario below is a sensitivity calculation, not a change to membership. They are evaluated independently and their gains must not be added.

## The two large drops, evaluated against current evidence

| Historical drop | Left CTA then | Already restored | Still absent | Tumor ≥1 TPM | Tumor ≥10 TPM | Tumor ≥10× reference |
|---|---:|---:|---:|---:|---:|---:|
| 1.38 | 173 | 3 | 170 | 150 | 56 | 2 |
| 1.40 | 191 | 44 | 147 | 147 | 147 | 124 |

Tumor counts in the final three columns concern only proteins still absent. Unassessed proteins are not called tumor-negative. Counts refer to selected-reference protein groups. One passing tumor sample already suffices; independent patient counts and paper counts are optional filters. The 1.38 and 1.40 historical populations can overlap.

## Reproductive scope is already broad

All measured breast, prostate, testis, ovary, placenta, epididymis, seminal-vesicle, fallopian-tube, endometrial and cervical rows are reproductive, not normal. They contribute zero rows to the default normal denominator. Both mammary states are included. Neither anatomy preset protects reproductive tissues by default. [Exact scope counts](reproductive-tissue-scope.csv).

## Evidence recovery versus changing the rule

| Counterfactual, holding other stated gates fixed | Additional CTA | Total | Interpretation |
|---|---:|---:|---|
| Use existing native-TPM reproductive tissue summaries | +10 | 300 | Retrieve and verify individual donors; summary values are currently excluded from the positive witness route. |
| Resolve all normal-review holds that otherwise pass | +15 | 305 | An upper bound: resolution can confirm normal expression and retain/exacerbate exclusion. |
| Remove each protein’s highest individual normal library | +8 | 298 | Dependence on one sample; not evidence that the library is erroneous. |
| Omit all individual HPA normal observations | +17 | 307 | Sensitivity only; no existing normal-library QC flag justifies blanket removal. |
| Add the previously audited reliable High-IHC route | +9 | 299 | Requires antibody/cell review; missing tissues remain a limitation, including SP7 bone cells. |
| Allow ≥3 cells in early embryos, still ≥50% positive and two embryos | +2 | 292 | A support-policy change, not the already-fixed 20-cell bug. |
| Accept seven TGCT libraries rejected for the testis signature | +13 | 303 | Requires malignant-cell attribution; every added protein peaks in TCGA-2G-AAFL-01. |
| Combine normal-review resolution and reproductive summaries | +25 | 315 | All 25 require the corresponding evidence resolution. |
| Also omit individual HPA normals | +47 | 337 | Optimistic sensitivity with simultaneous unaccepted assumptions; not a forecast. |

Simply reducing the tumor fold while preserving RA 10× would retain **315 at 8×**, **358 at 5×**, **394 at 2×**, **401 at 1×**. Reaching 400 this way changes cancer specificity; it is not correction of a tissue label or a missing sample. No such change is activated.

Requiring only a QC-passing tumor measurement ≥1 TPM while preserving RA 10× would retain **403 proteins**. This is a tumor-detected RA set, not the released cancer-enriched CTA set. The distinction explains how a 400-plus research list can coexist with 290 tenfold CTAs.

## Prioritized work

Tracked follow-ups: [reproductive donors #710](https://github.com/pirl-unc/oncoref/issues/710), [normal-cell attribution #711](https://github.com/pirl-unc/oncoref/issues/711), [normal-sample peaks #712](https://github.com/pirl-unc/oncoref/issues/712), and [TGCT attribution #713](https://github.com/pirl-unc/oncoref/issues/713).

1. **Native reproductive donor data.** Ten measured leads already exceed 10× in a native-TPM summary. Verify the current gene/sequence identity, individual sample, donor, tissue and QC; for multi-locus protein groups require within-sample co-occurrence instead of adding separate gene medians. Then evaluate the unchanged 10× rule.
2. **Normal-cell attribution.** Start with CAGE1, RAD21L1, PPP3R2, RGSL1 and DAZL. Distinguish uniquely spliced/exonic signal from host pre-mRNA, paralogs and study-specific counting. Use public spliced/unspliced counts or targeted reads where necessary. Negative bulk/IHC and a mapping-risk flag alone cannot establish artifact. Record gene-specific adjudications centrally. CTAG2 also fails reproductive 10× (3.15×); clearing its review alone is insufficient.
3. **Normal-library attribution.** Review the eight single-library-dependent cases before any broad statistic change. Check provider identity, cross-tissue signatures, technical support and replication in the same cell/tissue. Preserve the raw observation and an explicit supersession if rejected.
4. **Tumor evidence.** Recheck the 25 currently RA-only proteins at 8–10× using additional QC-passing cohorts and donor-resolved malignant-cell data. Thirteen additional leads depend on excluded TGCT libraries; require tumor attribution before restoring those observations. Do not relax the cutoff or reintroduce a contaminated library merely to recover a label.
5. **Broader gaps.** [#638](https://github.com/pirl-unc/oncoref/issues/638) provides individual literature adjudication; [#687](https://github.com/pirl-unc/oncoref/issues/687) requires a direct hair-follicle gate replay and downstream sequence recomputation; [#707](https://github.com/pirl-unc/oncoref/issues/707) covers the concrete IHC/SP7 gap. [#647](https://github.com/pirl-unc/oncoref/issues/647) concerns Disposable/Silent policy integration. [#665](https://github.com/pirl-unc/oncoref/issues/665) could recover valid tumor witnesses by rebuilding the currently invalid UNC NUTM source; its gain cannot be estimated before quantification. None currently supplies an established path to 400 CTAs. Paper-count fixes change optional prioritization, not this broad catalog.

## Gene lists

### Use existing native-TPM reproductive tissue summaries

TTLL8, TMEM95, SPANXD, CTAGE1, CT47B1, CCDC166, SLC35G6, CTAGE6, GOLGA6L1, C1orf202

### Resolve all normal-review holds that otherwise pass

DAZL, RGSL1, PSG8, CAGE1, ZNF280A, TAS2R1, POTEC, PPP3R2, DCAF8L2, FBXO47, VCX3B, PSG3, RAD21L1, OOSP4B, TUBB8

### Remove each protein’s highest individual normal library

C4orf17, TTLL8, GOLGA6L2, CCDC87, RNASE13, GAB4, PSG7, GSTT4

### Omit all individual HPA normal observations

CDX4, C4orf17, TTLL8, GOLGA6C, GOLGA6L2, CCDC197, ADIG, CCDC87, PFN3, STPG3, RNASE13, GAB4, PSG7, AK4P3, CSAG3, SPATA31A7, GSTT4

### Add the previously audited reliable High-IHC route

PKDREJ, LCN9, SSX5, SP7, MAGEB18, OTUD6A, DEFB136, RBMY1E, RBMY1B

### Allow ≥3 cells in early embryos, still ≥50% positive and two embryos

MOS, PRR20A/B/C/D/E

### Accept seven TGCT libraries rejected for the testis signature

ZNRF4, PRM2, CCDC70, BPIFA3, CYLC2, SPATA19, PGK2, CETN1, TEX38, PRR30, TSSK2, IQCF5, BOD1L2

### Combine normal-review resolution and reproductive summaries

DAZL, RGSL1, PSG8, TTLL8, CAGE1, ZNF280A, TAS2R1, TMEM95, POTEC, PPP3R2, DCAF8L2, SPANXD, FBXO47, VCX3B, CTAGE1, PSG3, CT47B1, RAD21L1, CCDC166, OOSP4B, SLC35G6, TUBB8, CTAGE6, GOLGA6L1, C1orf202

### Also omit individual HPA normals

FKBP6, DAZL, RGSL1, PSG8, CDX4, HHLA1, C4orf17, TTLL8, LRRC52, CAGE1, GOLGA6C, ZNF280A, TAS2R1, GOLGA6L2, TSGA10IP, CCDC197, OR4N2, ADIG, CCDC87, TMEM95, POTEC, PPP3R2, DCAF8L2, SPANXD, PFN3, STPG3, FBXO47, VCX3B, RNASE13, CTAGE1, GAB4, PSG3, PSG7, AK4P3, CT47B1, PSG11, RAD21L1, CCDC166, OOSP4B, SLC35G6, TUBB8, CSAG3, CTAGE6, GOLGA6L1, SPATA31A7, GSTT4, C1orf202

## Reproduce and inspect

`python scripts/audit_cta_recovery.py` reads pinned catalog observations and the recorded TGCT matrix; `--tgct-matrix PATH` permits another local copy with the exact expected checksum. It never writes catalog files. The standard embryo rule must replay with zero additions, the baseline must reproduce the exact CTA set, and no scenario may admit an Unsafe protein.

[Candidate measurements and every scenario](candidate-evidence.csv.gz) · [TGCT witnesses](tgct-qc-counterfactual.csv.gz) · [Early embryo sensitivity](early-embryo-counterfactual.csv.gz) · [Summary](summary.json) · [Input/output receipt](manifest.json)
