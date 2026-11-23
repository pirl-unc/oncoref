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
