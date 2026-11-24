# Assay-Boundary Adjudications & 1-Mismatch Near-Self Audit (#768, #769)

Policy ID: `oncoref.assay-boundary-and-near-self.2026-11-23.1`

## 1. Assay-Boundary Adjudications (`#768`)

1. **`#713` Non-admixed tumor witness gating (`TCGA-TGCT` Note `R` libraries):**
   `11` post-meiotic spermatid markers whose `>= 10x` tumor evidence rested solely on `library_qc` Note `R` parenchymal-admixture `TCGA-TGCT` libraries (`CETN1`, `PGK2`, `PRM2`, `BOD1L2`, `CCDC70`, `TSSK2`, `PRR30`, `ZNRF4`, `IQCF5`, `BPIFA3`, `CYLC2`) reclassify from `CTA -> Reproductive Antigen only`, whereas **`PRM1` stays `CTA`** via its independent female `BRCA` witness (`58.97 TPM >= 57.11 TPM` under `#764` bleed QC) and **`SPATA19` stays `CTA`** via its `#737` `cSCC` somatic tumor witness (`GSM3567709 = 13.28 TPM >= 10.0 TPM`, matching `CANCER_ONLY` MS evidence).
2. **RNA-zero single-antibody HPA IHC reconciliation (`cta-ihc-unreliable.csv`):**
   **`MAGEB6`** (`0.0118 TPM` somatic normal max, `0.0 nTPM`, `100%` required RNA coverage complete, `22.14 TPM` Testis, `10.44 TPM` `MM`) advances `Unsafe -> CTA`, while `NLRP9` (`2.32 TPM` somatic bulk) and `XAGE2` remain `Unsafe`.
3. **`Adult cornea` (`GSE155683`) multi-UMI cell support floor (`>= 10` cells with `>= 2` UMIs):**
   **`DKKL1`** (`8` multi-UMI cells across `21,343` adult cornea cells; `3.39 TPM` normal max, `1,346.2 TPM` Testis, `61.62 TPM` `BRCA`) advances `Unsafe -> CTA`, while `INSL4` (`64`), `GPHB5` (`24`), `QRFP` (`40`), and `C20orf202` (`42`) remain `Unsafe`.
4. **Protamine-diluted germline transcription-factor qualification:**
   **`FOXR2`** (`135` Testis donors `>= 3 TPM`, `5.2 nTPM` spermatogonia), **`SPATA31A3`** (`23` Testis donors `>= 3 TPM`, `1.2 nTPM` spermatids), and **`PPP4R3C`** (`230` Testis donors `>= 3 TPM`, `22.7 nTPM` germ cells) advance `Cancer Antigen only -> CTA`, while single-donor `DEFB125` (`1` donor) remains `Cancer Antigen only`.

## 2. 1-Mismatch Hamming Near-Self & Tryptic-MS Paralog Audit (`#769`)

- **`MAGEA6`:** `306` total 9-mers -> `86` exact-unique -> **`16` safe `>= 2-mm`** (`81.40%` 1-mm near-self against `Unsafe` `MAGEA3`/`MAGEA12`).
- **`MAGEA2/B`:** `306` total -> `165` exact-unique -> **`58` safe `>= 2-mm`** (`64.85%` 1-mm near-self).
- **`PAGE2`:** `103` total -> `21` exact-unique -> **`12` safe `>= 2-mm`** (`42.86%` 1-mm near-self).
- **`CTAGE1` (`ENSG00000212710`) vs. `MIA2` (`ENSG00000150527`):** `55 / 138` (`39.86%`) canonical tryptic peptides (`I == L`) are shared directly with ubiquitous `Unsafe` `MIA2 (49); CTAGE4 (8); CTAGE9 (8); ACR (1); ANKRD11 (1); CEP135 (1)`, explaining healthy multi-organ bottom-up MS attributions without requiring `CTAGE1`-specific translation.
