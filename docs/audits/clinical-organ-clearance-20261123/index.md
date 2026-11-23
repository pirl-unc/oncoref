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
