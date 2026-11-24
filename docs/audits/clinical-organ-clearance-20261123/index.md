# Clinical pMHC-Directed Organ-Safety Evidence Tier (`#766`)

**Policy ID:** `oncoref.clinical-organ-clearance.2026-11-23.1`  
**Ensembl Release:** `115` (`GRCh38`)

## Guardrails & `tsarina` Isoform Contract

- **Isoform- & epitope-scoped clearance (`CTAG2` `LAGE-1a` vs. `LAGE-1b`):**
  `ENST00000369585` (`LAGE-1a`, `180 aa`) contains `SLLMWITQC` (`157–165`), whereas
  `ENST00000247306` (`LAGE-1b`, `210 aa`) lacks `SLLMWITQC` and is explicitly listed in
  `uncleared_coding_transcripts`.
- **Vascular perfusion & BBB protection:** Systemic IV `TCR-T` clears `lung` and `urinary bladder`
  (`CTAG2` also exceeds its `6.92 TPM` GTEx Urinary Bladder baseline by `87.3x` at `604.33 TPM` in `OV`),
  while never clearing `brain` clinically (resolved separately via `#711` `Xq28` attribution) or `heart`
  (`#544` confounded day-161 cardiac event preserved).

## Reviewed Entries

| Clearance ID | Target | Product | Status | Patients | Cleared Organs | Cleared Transcripts | Uncleared Transcripts |
| --- | --- | --- | --- | ---: | --- | --- | --- |
| `COC-20261123-CTAG2-LETE-CEL` | `CTAG2` | `lete_cel_c259` | `cleared` | `20` | `lung,urinary bladder` | `1` | `1` |
| `COC-20261123-NYESO1-LETE-CEL` | `NY-ESO-1` | `lete_cel_c259` | `cleared` | `20` | `lung,urinary bladder` | `2` | `2` |
| `COC-20261123-PRAME-IMA203` | `PRAME` | `ima203` | `cleared` | `40` | `kidney,lung` | `5` | `1` |
| `COC-20261123-MAGEA4-AFAMI-CEL` | `MAGEA4` | `afami_cel_c1032` | `cleared` | `90` | `lung,liver,kidney` | `26` | `1` |
| `COC-20261123-NEG-MAGEA3-TITIN` | `MAGEA3` | `magea3_a3a` | `blocked_clinical_toxicity` | `2` | `-` | `0` | `0` |
| `COC-20261123-NEG-MAGEA12-NEURO` | `MAGEA12` | `magea3_a2_neuro` | `blocked_clinical_toxicity` | `9` | `-` | `0` | `0` |
| `COC-20261123-NEG-TBI1301-LUNG` | `CTAG2` | `tbi1301` | `blocked_unresolved_event` | `9` | `-` | `0` | `0` |

## `CTAG2` (`LAGE-1`) Resolution & Downstream `tsarina` Schema

`CTAG2` reaches `604.33 TPM` in ovarian carcinoma (`OV`, `375` known independent patients `≥ 10×`) and `102.0 TPM` in testis over a `6.92 TPM` GTEx Urinary Bladder baseline (`87.3×` tumor fold, `14.7×` reproductive fold), but was previously trapped in `Other:Review` by two sparse `snRNA-seq` flags (`4.3%` lung alveolar type 2 cells and cortical neurons) combined with `90.2%` 31-mer sharing across the closed `Xq28` `CTAG1A`/`CTAG1B`/`CTAG2` family:

1. **Lung AT2 & bladder clearance (`#766`):** `lete-cel` (`c259`, `HLA-A*02:01 / SLLMWITQC`) infuses billions of high-affinity engineered T cells intravenously—whose first capillary bed is the pulmonary microvasculature—across synovial sarcoma, MRCLS, NSCLC, and myeloma cohorts without on-target pulmonary or urothelial toxicity, clearing `lung` and `urinary bladder` specifically for `ENST00000369585` (`LAGE-1a`, `180 aa`).
2. **Brain neuronal `Xq28` attribution (`#711` / `#767`):** `CTAG2`'s external `Unsafe` sequence overlap is `< 1%` (`0.0%` outside `CTAG1A/B`), bulk brain polyA RNA is `0.15 TPM`, and systemic IV TCR-T is forbidden from claiming clinical clearance across the blood–brain barrier (`validate()` raises `ValueError` on any BBB-protected organ).
3. **Machine-readable columns on `target_selection.results()` and `sequence_bundle`:**
   - `clinical_clearance_scope`: `"epitope"` (never `"gene"`)
   - `clinical_cleared_epitopes`: `["SLLMWITQC"]`
   - `clinical_hla_alleles`: `["HLA-A*02:01", "HLA-A*02:05", "HLA-A*02:06"]`
   - `clinical_cleared_organs`: `["lung", "urinary bladder"]`
   - `clinical_cleared_transcripts`: `["ENST00000369585"]`
   - `clinical_unverified_isoforms`: `["ENST00000247306"]` (`LAGE-1b`, `210 aa`, intron-2 retention lacking `SLLMWITQC`)

## Programmatic Inspection

```python
import json
from oncoref import clinical_organ_clearance, target_selection

ledger = clinical_organ_clearance.clearances_table()
ctag2 = target_selection.results("CTAG2").iloc[0]
print(ctag2["classification"], json.loads(ctag2["clinical_unverified_isoforms"]))
```
