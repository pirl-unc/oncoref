# Nuclear snRNA-seq Compartment & Transcriptional Spillover QC (`#765`)

**Policy ID:** `oncoref.nuclear-spillover-qc.2026-11-23.1`  
**Ensembl Release:** `115` (`GRCh38`)  
**CELLxGENE Census Release:** `2025-11-08` (`93` primary adult collections)

## Overview

pMHC-directed T-cell cytotoxicity requires spliced cytoplasmic `CDS` mRNA, translation,
and surface HLA presentation. Intron-inclusive single-nucleus RNA-seq (`snRNA-seq`,
`assay_compartment == "nucleus"`) captures unspliced pre-mRNA across long non-coding
5'-UTR introns (`PLAC1`, `TENT5D`), termination readthrough downstream of active
same-strand somatic neighbors (`PABIR2 -> PLAC1`, `EPO -> ZAN`), and divergent
promoter-upstream transcripts (`NOP10 <-> NUTM1`, `135 bp` head-to-head CpG island).

## Guardrails Enforced (`tsarina` Cross-Repo Contract)

1. **Bulk-dilution & literature/IHC guardrails (`MAGEA12` regression anchor):**
   `MAGEA12` (`finding_bulk_tpm = 0.30 TPM` in brain, `literature_effect = confirm_detection`)
   remains strictly `Unsafe` (`retained_literature_or_ihc_supported`). Clearance requires
   all four conditions simultaneously: (a) `nucleus_only_signal == True` (zero cytoplasmic
   `scRNA-seq` multi-UMI populations), (b) `finding_bulk_tpm < 0.05 TPM`, (c) negative
   somatic IHC and no `exclude`/`confirm_detection` literature finding, and (d) verified
   Ensembl 115 spillover geometry (`< 1 kb` bidirectional promoter `>= 10 TPM`, `< 10 kb`
   same-strand upstream partner `>= 5 TPM`, or `>= 50 kb` non-coding 5'-UTR span with
   `CDS`/span ratio `<= 5%`).
2. **Intronless / single-CDS-exon coding support (`DUX4`, `PLAC1`):**
   Records `is_intronless_cds` and explicitly accepts `cds_exonic_reads` alongside
   `spliced_exon_junction` for coding-sequence attribution.
3. **Downstream `nuclear_only` disposition:**
   Cleared loci emit `effective_rare_cell_status = "nuclear_only"` and
   `nuclear_rna_caution = "snrnaseq_intronic_spillover_cleared"`, advancing loci that
   satisfy all `10x` gates (`PLAC1`: `152.6 TPM` `WILMS`, `187.0 TPM` placenta vs `0.82 TPM`
   normal; `ZAN`: `10.5 TPM` `TGCT` / `10.2 TPM` `B_ALL`, `16.7 TPM` testis vs `0.51 TPM`
   normal) directly into `CTA` rather than trapping them in `Unsafe` or `Other:Review`.

## Cleared Nuclear / Promoter Spillover Loci

| Protein | Baseline | Finding Organ | Bulk TPM | Spillover Mechanism | Key Neighbor / Span | Disposition |
| --- | --- | --- | ---: | --- | --- | --- |
| `METTL21EP` | `Unsafe` | `brain` | `0.0000` | `same_strand_upstream_readthrough` | `ERCC5:ENSG00000134899` | `nuclear_spillover_cleared` |
| `PLAC1` | `Unsafe` | `brain` | `0.0051` | `same_strand_upstream_readthrough;long_noncoding_5prime_intron_span` | `PABIR2:ENSG00000156504` | `nuclear_spillover_cleared` |
| `TENT5D` | `Other:Review` | `brain` | `0.0000` | `long_noncoding_5prime_intron_span` | `5'-UTR 107,036 bp (1.07% CDS)` | `nuclear_spillover_cleared` |
| `ZAN` | `Unsafe` | `brain` | `0.0000` | `same_strand_upstream_readthrough` | `EPO:ENSG00000130427` | `nuclear_spillover_cleared` |
| `NUTM1` | `CTA` | `bone marrow (Tabula Sapiens 10x neutrophil audit)` | `0.0000` | `bidirectional_promoter_prompt` | `NOP10:ENSG00000182117` | `promoter_spillover_immune_explained` |

## Locus-Specific Biological & Orthogonal Proteomic Evidence

- **`PLAC1` (`Unsafe → CTA`, `CTA-M`):** Canonical trophoblast surface protein whose alternative placental promoter spans a `197,673 bp` non-coding 5'-UTR intron (`CDS` is a single intronless exon comprising `0.32%` of locus span) located `5,211 bp` downstream of ubiquitous same-strand `PABIR2` (`54.1 TPM` in brain). Cortical glutamatergic `snRNA-seq` collections capture nuclear readthrough/intronic pre-mRNA (`9.2%` of nuclei), whereas cytoplasmic polyA bulk brain RNA is `0.0051 TPM`, `0/4` somatic HPA IHC tissues stain positive, and immunopeptidomics/MS datasets detect `PLAC1` exclusively in tumors (`CANCER_ONLY`). With peripheral baseline `0.82 TPM`, `PLAC1` clears `10×` in both Wilms tumor (`152.6 TPM`) and placenta (`187.0 TPM`).
- **`NUTM1` (retained `CTA` with `promoter_spillover_immune_explained`):** Shares a `135 bp` bidirectional CpG-island promoter with core H/ACA snoRNP subunit `NOP10` (`4,501 TPM` in bone-marrow neutrophils). Sparse `10x` neutrophil droplet UMIs originate from divergent promoter initiation near `NOP10`'s TSS, whereas sorted peripheral neutrophil bulk polyA RNA (`HPA`/`Monaco`) is `0.00 TPM` and MS proteomics across 7 studies is strictly `CANCER_ONLY` (`56.2 TPM` in `NUTM1`-rearranged `LAML_ELNadv`).
- **`DUX4` (joint `#765 + #767 → CTA`, `CTA-E`):** Clears its single-study nuclear `snRNA-seq` cortical neuron hold (`0.0000 TPM` bulk brain, intronless `cds_exonic_reads` contract) under `#765`, and simultaneously clears `D4Z4` macrosatellite `paralog_ambiguous` under `#767`.

## Programmatic Inspection

```python
from oncoref import rare_cell_compartment_qc, target_selection

compartments = rare_cell_compartment_qc.census_collection_compartments()
adjudications = rare_cell_compartment_qc.nuclear_spillover_adjudications()

target_selection.results().query(
    "protein_group_label in ['PLAC1', 'ZAN', 'TENT5D', 'NUTM1', 'DUX4', 'MAGEA12']"
)[
    [
        "protein_group_label",
        "classification",
        "rare_cell_status",
        "rare_cell_action",
        "nuclear_spillover_765_disposition",
        "nuclear_rna_caution",
    ]
]
```
