# Closed Homologous Repeat-Family Sequence Equivalence (#767)

Policy ID: `oncoref.repeat-family-equivalence.2026-11-23.1`

## 1. Sequence-Partner Partitioning (`P_self` vs. `P_ext`)

Across all `209` `Other:Review` candidates carrying `paralog_ambiguous = True`, sequence partners are partitioned into:

- **`P_self(g)` (closed self-family homologs / repeat copies):** primary-assembly `protein_coding` and `unprocessed_pseudogene` loci sharing `>= 85%` `CDS` 31-mers or `>= 90%` coding amino-acid sequence with `g` and exhibiting zero `Unsafe` somatic expression (`unsafe_sequence_overlap_pct < 5.0%`, somatic family donor `P95 < 5.0 TPM`).
- **`P_ext(g)` (external / heterologous paralogs):** `Unsafe` coding paralogs (e.g., `PSG8 <-> PSG1` `74.07%`, `POTEC <-> POTED/E` `46.47%`, `TUBB8B <-> TUBB*` `50.0%`, `GAGE2E <-> GAGE12F/G` `100%`), `processed_pseudogene` retrocopies, or somatic `lncRNA`/chimeric hosts.

Only `2` candidates satisfy `P_self != empty` and `P_ext == empty` with somatic family ceiling `< 5.0 TPM`.

## 2. `DUX4` & `D4Z4` (`4q35` / `10q26`) Symmetric Multimapper Evaluation

Summing natively across `DUX4` (`ENSG00000260596`) and all `20` primary-assembly `D4Z4` repeat units (`>= 85%` `CDS` 31-mer recall: `DUX4L1..DUX4L8` on `4q35` and `DUX4L10..15, DUX4L20..25` on `10q26`):

- **GTEx v11 Skeletal Muscle (`FSHD` locus):** `0.0000 TPM` donor `P95` (confirming normal adult skeletal muscle silence).
- **GTEx v11 Skin (`Sun_Exposed` / `Not_Sun_Exposed`):** `0.0000 TPM` donor `P95`.
- **Highest somatic organ (`rna_tissue_hpa:somatic_max`):** `0.3000 TPM` donor `P95` (HPA genome-aligned `rna_tissue_hpa` somatic max = `0.40 TPM`; note `hpa_tissue_samples:v23` transcript-only `Kallisto` EM inflates coding `DUX4` because `unprocessed_pseudogene` `DUX4L*` decoys are absent from its transcript index).
- **Joint `#765 + #767` dependency:** Neither `#765` alone (`paralog_ambiguous` still blocks) nor `#767` alone (brain `snRNA-seq` nuclear flag still blocks) admits `DUX4`; enabling both `#765` and `#767` resolves both holds.

## 3. Downstream Coding Closed-Family Contract (`tsarina`)

`closed-homolog-families.csv` (`48` rows across `12` closed families) populates `closed_homolog_family_id`, `co_silenced_coding_sibling_ids`, and `closed_self_paralog_ids` on `target_selection.results()`, `antigen_sets`, and `sequence_bundle.target_selection_dispositions()` so `tsarina`'s `cta_exclusive_peptides()` retains epitopes shared with somatic-safe germline siblings (`CT47C1`, `TSPY9`, `SPANXN2/5`, `CSAG2/3`) while strictly subtracting `Unsafe` and `Other:Review` somatic paralogs.
