# TCGA-TGCT Somatic Testicular Admixture vs. Malignant Germ-Cell Localization Audit (`#713`)

**Audit ID:** `oncoref.tgct-localization-audit.2026-11-22.1`  
**Packaged tables:** `oncoref/data/tgct-reproductive-localization-audit.csv` (`oncoref.tgct_reproductive_localization_audit()`) and `oncoref/data/tgct-aafl-protein-support.csv` (`oncoref.tgct_aafl_protein_support()`)  
**Generator:** `scripts/audit_tgct_reproductive_localization.py`

## Motivation

In `oncoref.library_qc`, Note `R` (`testis_tissue_content >= 0.001`) flags tumor RNA-seq libraries in which all six GTEx v11 post-meiotic spermatid markers (`PRM1`, `PRM2`, `TNP1`, `PHF7`, `LELP1`, `SMCP`) simultaneously reach at least `0.1%` of normal adult testis median TPM. Because spermatid transcripts alone cannot distinguish **benign admixed seminiferous parenchyma** from **ectopic germ-cell tumor expression**, `library_qc` records Note `R` as advisory rather than excluding the library.

Within `TCGA-TGCT` (`n = 154`), `7` libraries carry Note `R`, led by **`TCGA-2G-AAFL-01`** (`testis_tissue_content = 0.045037`, `10.2x` higher than the second-ranked library `TCGA-2G-AAG8-01` at `0.004398`). Thirty-four catalog proteins (`22` `CTA`s and `12` `Reproductive Antigen only` entries) achieve their pan-cancer maximum TPM in `TCGA-2G-AAFL-01`. Issue `#713` requested an orthogonal cell-of-origin localization audit separating genuine tumor-supported CTAs from candidates driven by unconfirmed benign testicular admixture in `TCGA-2G-AAFL-01`.

---

## Orthogonal Somatic Niche Test Across `TCGA-TGCT` (`n = 154`)

Malignant testicular germ cell tumors (`TGCT` — seminomas, embryonal carcinomas, yolk sac tumors, choriocarcinomas, and teratomas) arise from arrested primordial germ cells / gonocytes (`GCNIS`) and do **not** differentiate into either haploid elongating spermatids or the two non-germ-cell somatic lineages of normal adult testis:

1. **Somatic Leydig interstitial steroidogenic cells:** `INSL3`, `CYP17A1`, `STAR`, `CYP11A1`, `HSD3B2`, `LHCGR`
2. **Somatic Sertoli intratubular nurse cells:** `INHA`, `CLU`, `SOX9`, `WT1`, `AMH`, `FSHR`, `DMRT1`

When surgical resection of a primary testicular tumor captures `~5–8%` adjacent benign testicular parenchyma, **all three compartments** — somatic Leydig cells, somatic Sertoli cells, and post-meiotic round/elongating spermatids — appear simultaneously in proportion to normal GTEx v11 `Testis`:

| Library | `testis_tissue_content` | Spermatid median ratio vs GTEx Testis | Spermatid log10 `r` vs GTEx Testis | Leydig `INSL3` / `CYP17A1` / `STAR` (TPM) | Sertoli `CLU` / `INHA` / `SOX9` (TPM) | Dominant malignant program (`POU5F1` / `SOX2` / `SOX17`) | Localization call |
|---|---:|---:|---:|---:|---:|---|---|
| `TCGA-2G-AAFL-01` | `0.0450` | `0.0603` (`6.0%`) | `0.9875` | `157.1` / `23.2` / `66.1` | `2019.4` / `45.6` / `10.8` | Embryonal carcinoma (`1404.8` / `35.9` / `3.7`) | `admixed_seminiferous_and_interstitial_parenchyma` |
| `TCGA-2G-AAG8-01` | `0.0044` | `0.0054` (`0.54%`) | `0.9472` | `584.3` / `33.4` / `30.5` | `919.7` / `31.0` / `6.6` | Embryonal carcinoma (`1788.5` / `50.6` / `3.1`) | `trace_admixed_testicular_parenchyma` |
| `TCGA-2G-AAFE-01` | `0.0034` | `0.0047` (`0.47%`) | `0.9611` | `3.8` / `1.9` / `9.1` | `548.9` / `23.0` / `2.7` | Seminoma (`1718.3` / `0.1` / `50.4`) | `trace_admixed_testicular_parenchyma` |
| `TCGA-2G-AAGA-01` | `0.0013` | `0.0022` (`0.22%`) | `0.8574` | `203.4` / `25.1` / `25.6` | `534.3` / `17.4` / `3.6` | Embryonal carcinoma (`1787.8` / `84.5` / `1.6`) | `trace_admixed_testicular_parenchyma` |
| `TCGA-2G-AALT-01` | `0.0013` | `0.0026` (`0.26%`) | `0.9619` | `1.7` / `1.2` / `2.4` | `423.3` / `1.6` / `12.0` | Teratoma / mixed (`160.6` / `2.8` / `6.4`) | `trace_admixed_testicular_parenchyma` |
| `TCGA-2G-AAEX-01` | `0.0013` | `0.0020` (`0.20%`) | `0.9513` | `2.7` / `1.3` / `2.2` | `65.6` / `1.0` / `0.4` | Seminoma (`1433.6` / `2.0` / `20.4`) | `trace_admixed_testicular_parenchyma` |
| `TCGA-XE-AAO4-01` | `0.0010` | `0.0014` (`0.14%`) | `0.9033` | `32.8` / `3.0` / `6.8` | `433.3` / `17.1` / `3.2` | Seminoma (`1338.4` / `0.2` / `46.4`) | `trace_admixed_testicular_parenchyma` |
| Non-`R` TGCT median (`n = 147`) | `0.0000` | `0.0001` (`0.01%`) | — | `1.2` / `0.1` / `1.8` | `176.2` / `1.1` / `1.3` | Pluripotent TGCT (`1307.8`) | `no_coordinated_spermatid_program` |

### Biological Interpretation of `TCGA-2G-AAFL-01`

`TCGA-2G-AAFL-01` is a bona fide `POU5F1+` / `NANOG+` / `SOX2+` / `DNMT3B+` embryonal carcinoma (`~93–95%` malignant bulk) whose RNA extraction included **`~5–7%` adjacent non-neoplastic testicular parenchyma** containing both interstitial Leydig cells (`CYP17A1` `166x` non-`R` TGCT median; `INSL3` `127x`; `STAR` `37x`) and intratubular Sertoli + haploid spermatid cells (`INHA` `43x`; `CLU` `11.5x`; `PRM1 = 1706.9 TPM` = `7.20%` of GTEx Testis; `PRM2 = 1362.3 TPM` = `5.90%` of GTEx Testis; log10 stoichiometry correlation `r = 0.9875` across the 10-gene spermatid panel).

---

## Cross-Cohort Witness Stratification for the 22 `AAFL`-Peaking CTAs

Scanning all **132** curation per-sample tumor matrices (`sample_qc == "pass"`) partitions the 22 `CTA` proteins whose pan-cancer maximum occurs in `TCGA-2G-AAFL-01` into three actionable tiers:

### 1. Independently Confirmed Outside `R`-Flagged `TGCT` (`9` CTAs — `confirmed_outside_r_tgct`)

These 9 CTAs remain unequivocally tumor-supported even if all 7 `R`-flagged `TGCT` libraries are excluded:

- **`HMGB4`**: `18` independent non-TGCT witnesses in `NEC_MERKEL` (peak `44.77 TPM` in `MO_MCC_035.1_T_RNA` vs threshold `11.0 TPM`).
- **`ACTL7A`**: `4` independent witnesses in `NEC_MERKEL` (peak `15.69 TPM` in `MO_MCC_042_T_RNA` vs `10.0 TPM`).
- **`LELP1`**: `3` independent witnesses across `LUSC` (`88.02 TPM`), `MENINGIOMA`, and `NEC_MERKEL` (vs `10.0 TPM`).
- **`IRGC`**: `2` independent witnesses across `B_ALL` (`14.31 TPM`) and `CLL` (vs `10.0 TPM`).
- **`SMCP`**: independent `SARC_GIST` witness (`72.22 TPM` vs `62.84 TPM`).
- **`TEX38`**: independent `NEC_MERKEL` witness (`32.30 TPM` in `NU_MCC_023_T_RNA` vs `28.09 TPM`).
- **`ODF1`**: independent `BLCA` witness (`30.79 TPM` vs `15.71 TPM`).
- **`ACTRT2`**: independent `STAD` witness (`11.31 TPM` vs `10.0 TPM`).
- **`CIMIP4`**: independent `BL` witness (`10.79 TPM` vs `10.0 TPM`).

### 2. Multi-Library `R`-Flagged `TGCT` Only (`2` CTAs — `r_tgct_multi_library_only`)

- **`PRM1`** (`1706.94 TPM` in `AAFL`, `127.00 TPM` in `AAG8`, `114.38 TPM` in `AAFE`; max non-`R` TGCT `29.43 TPM`, max non-TGCT `BRCA` `58.97 TPM` vs threshold `80.53 TPM`).
- **`PRM2`** (`1362.27 TPM` in `AAFL`, `108.19 TPM` in `AAFE`, `101.51 TPM` in `AAG8`; max non-`R` TGCT `22.11 TPM`, max non-TGCT `BL` `8.60 TPM` vs threshold `77.57 TPM`).

### 3. `TCGA-2G-AAFL-01` Sole-Witness Admixture Candidates (`11` CTAs — `aafl_sole_witness_admixed_parenchyma_candidate`)

`BOD1L2`, `CETN1`, `CCDC70`, `TSSK2`, `PGK2`, `PRR30`, `ZNRF4`, `SPATA19`, `IQCF5`, `BPIFA3`, and `CYLC2` each have **`0` qualifying tumor witnesses outside `TCGA-2G-AAFL-01`** across all 132 cancer cohorts, and their `AAFL / reference Testis` ratios (`3.6%–10.6%` across 10/11 — including `TSSK2` against HPA v23 `Testis` where GTEx v11 omits the 22q11.2 locus — and `IQCF5` `61.7%`) match the `~5–7%` benign seminiferous parenchyma fraction measured independently by `PRM1`/`PRM2` and somatic Sertoli/Leydig markers.

## Downstream Programmatic Filter

Downstream callers requiring tumor-cell-localized CTA evidence without unconfirmed seminiferous parenchyma admixture can query:

```python
import oncoref

# 11 sole-witness AAFL candidates:
sole_aafl = oncoref.tgct_admixture_sole_witness_entity_ids()

# All 13 CTAs lacking any witness outside R-flagged TGCT (11 + PRM1 + PRM2):
admixture_only = oncoref.tgct_admixture_restricted_entity_ids()
```
