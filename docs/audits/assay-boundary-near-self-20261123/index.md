# Assay-Boundary Adjudications & 1-Mismatch Near-Self Audit (#768, #769, #779)

Policy ID: `oncoref.assay-boundary-and-near-self.2026-11-23.1`

## 1. Assay-Boundary Adjudications (`#768`, `#779`)

1. **`#713` Non-admixed tumor witness gating (`TCGA-TGCT` Note `R` libraries):**
   `11` post-meiotic spermatid markers whose `>= 10x` tumor evidence rested solely on `library_qc` Note `R` parenchymal-admixture `TCGA-TGCT` libraries (`CETN1`, `PGK2`, `PRM2`, `BOD1L2`, `CCDC70`, `TSSK2`, `PRR30`, `ZNRF4`, `IQCF5`, `BPIFA3`, `CYLC2`) reclassify from `CTA -> Reproductive Antigen only`, whereas **`PRM1` stays `CTA`** via its independent female `BRCA` witness (`58.97 TPM >= 57.11 TPM` under `#764` bleed QC) and **`SPATA19` stays `CTA`** via its `#737` `cSCC` somatic tumor witness (`GSM3567709 = 13.28 TPM >= 10.0 TPM`, matching `CANCER_ONLY` MS evidence).
2. **RNA-zero single-antibody HPA IHC reconciliation (`cta-ihc-unreliable.csv`):**
   **`MAGEB6`** (`0.0118 TPM` somatic normal max, `0.0 nTPM`, `100%` required RNA coverage complete, `22.14 TPM` Testis, `10.44 TPM` `MM`) advances `Unsafe -> CTA`, while `NLRP9` (`2.32 TPM` somatic bulk) and `XAGE2` remain `Unsafe`.
3. **`Adult cornea` (`GSE155683`) multi-UMI cell support floor (`>= 10` cells with `>= 2` UMIs):**
   **`DKKL1`** (`8` multi-UMI cells across `21,343` adult cornea cells; `3.39 TPM` normal max, `1,346.2 TPM` Testis, `61.62 TPM` `BRCA`) advances `Unsafe -> CTA`, while `INSL4` (`64`), `GPHB5` (`24`), `QRFP` (`40`), and `C20orf202` (`42`) remain `Unsafe`.
4. **Protamine-diluted germline transcription-factor qualification:**
   **`FOXR2`** (`135` Testis donors `>= 3 TPM`, `5.2 nTPM` spermatogonia), **`SPATA31A3`** (`23` Testis donors `>= 3 TPM`, `1.2 nTPM` spermatids), and **`PPP4R3C`** (`230` Testis donors `>= 3 TPM`, `22.7 nTPM` germ cells) advance `Cancer Antigen only -> CTA`, while single-donor `DEFB125` (`1` donor) remains `Cancer Antigen only`.
5. **Droplet ambient milk soup & `<= 100 kb` cis-chromatin spillover (`#779`):**
   **`CSN2`** (`ENSG00000135222` on `4q13.3`; `breast c-3` macrophages `683.2 -> 0.0 nTPM` after `single_cell_tissues(qc="all", cleanup="all")`; `salivary gland` serous clusters `30.2 kb` from `STATH`/`HTN3`/`HTN1` at `321,577–896,409 nTPM` with spillover ratio `<= 0.000050 <= 1e-4`, HPA surgical salivary bulk `0.90 TPM < 1.0 TPM`, salivary IHC `Not detected`, `0` CELLxGENE rare-cell flags, `2.99 TPM` clean somatic normal max, `28,460.4 TPM` lactating breast, `99,761.8 TPM` `BRCA`, `987.6 TPM` `BRCA_TNBC`) advances `Unsafe -> CTA` (`CTA-F` once coordinated `XY` galactorrhea biopsies are separated from resting male breast), while `CSN1S1` (`16.9 nTPM` brain neurons), `CSN3` (`261.4 nTPM` + `395.9 TPM` bulk), `ANKRD30A` (`4.9 nTPM` brain neurons), and `MAGEA3` remain strictly `Unsafe`.

## 2. 1-Mismatch Hamming Near-Self & Tryptic-MS Paralog Audit (`#769`)

- **`MAGEA6`:** `306` total 9-mers -> `86` exact-unique -> **`16` safe `>= 2-mm`** (`81.40%` 1-mm near-self against `Unsafe` `MAGEA3`/`MAGEA12`).
- **`MAGEA2/B`:** `306` total -> `165` exact-unique -> **`58` safe `>= 2-mm`** (`64.85%` 1-mm near-self).
- **`PAGE2`:** `103` total -> `21` exact-unique -> **`12` safe `>= 2-mm`** (`42.86%` 1-mm near-self).
- **`CTAGE1` (`ENSG00000212710`) vs. `MIA2` (`ENSG00000150527`):** `55 / 138` (`39.86%`) canonical tryptic peptides (`I == L`) are shared directly with ubiquitous `Unsafe` `MIA2 (49); CTAGE4 (8); CTAGE9 (8); ACR (1); ANKRD11 (1); CEP135 (1)`, explaining healthy multi-organ bottom-up MS attributions without requiring `CTAGE1`-specific translation.

### Summary of Controversial Target Dispositions & Near-Self Footprints

| Protein | Promoted Category | Sex Subset | Exact Unique 9-mers | Safe `>= 2-mm` 9-mers | `1-mm` Near-Self % | Shared `Unsafe` Tryptic Peptides | Key Rule / Note |
|---|---|---|---:|---:|---:|---:|---|
| `PLAC1` | `CTA` | `CTA-M` | `204` | `204` | `0.00%` | `0` | `#765` nuclear 5'-UTR intron & `PABIR2` readthrough cleared |
| `DKKL1` | `CTA` | `CTA-MF` | `234` | `234` | `0.00%` | `0` | `#768` sparse adult cornea (`8 / 21,343` multi-UMI cells) cleared |
| `MAGEB6` | `CTA` | `CTA-M` | `367` | `358` | `2.45%` | `0` | `#768` single-antibody IHC with `0.012 TPM` somatic RNA cleared |
| `CSN2` | `CTA` | `CTA-F` | `218` | `209` | `4.13%` | `0` | `#779` breast ambient milk soup + `4q13.3` SCPP cis-spillover + `XY` galactorrhea separated |
| `CTAG2` | `CTA` | `CTA-MF` | `172` | `172` | `0.00%` | `0` | `#766` `LAGE-1a` (`ENST00000369585`) `SLLMWITQC` clinical clearance |
| `DUX4` | `CTA` | `CTA-E` | `413` | `412` | `0.24%` | `0` | `#765 + #767` symmetric `D4Z4` repeat-array silence |
| `FOXR2` | `CTA` | `CTA-M` | `303` | `303` | `0.00%` | `0` | `#768` protamine-diluted spermatogonial TF (`135` Testis donors) |
| `SPATA31A3` | `CTA` | `CTA-M` | `1,303` | `1,254` | `3.76%` | `14` | `#768` germline TF; shares `14` tryptic peptides with `SPATA31C1/C2` |
| `PRM1` | `CTA` | `CTA-M` | `43` | `43` | `0.00%` | `0` | `#764`/`#776` Thyroid & opposite-sex bleed cleaned (`5.71 TPM`); female `BRCA` `58.97 TPM` |
| `PRM2` | `Reproductive Antigen only` | `RA-M-only` | `94` | `94` | `0.00%` | `0` | `#768`/`#776` `TCGA-TGCT` Note `R` parenchymal admixture demoted; female-safe `RA-M-only` |
| `SPATA19` | `CTA` | `CTA-M` | `159` | `159` | `0.00%` | `0` | `#768 + #737 + #776` somatic `cSCC` witness (`13.28 TPM`) + female breast `0.41 TPM` |
| `MAGEA6` | `CTA` | `CTA-MF` | `86` | `16` | `81.40%` | `9` | `#769` high 1-mm cross-reactivity risk vs. `Unsafe` `MAGEA3`/`MAGEA12` |
| `CTAGE1` | `CTA` | `CTA-M` | `539` | `520` | `3.53%` | `55` | `#769` `55/138` tryptic peptides shared with `Unsafe` `MIA2` |

## Programmatic Inspection

```python
from oncoref import assay_boundary_adjudications, peptide_subtraction, target_selection

adj = assay_boundary_adjudications.assay_boundary_adjudications()
near_self = assay_boundary_adjudications.candidate_near_self_1mm()
ms_audit = assay_boundary_adjudications.tryptic_ms_paralog_audit()
policy_1mm = assay_boundary_adjudications.hamming_1mm_policy()
```
