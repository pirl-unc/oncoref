# Release 1.46.0

Policy `oncoref.target-selection.2026-11-20.2` · Sequence-role policy `oncoref.target-selection-sequence-roles.2026-11-21.1` · Species policy `oncoref.species-screening-policy.2026-11-21.1`.

Across **19,893** human protein-equivalence groups (**20,099** Ensembl 115 coding loci), the catalog classifies **327 CTA, 69 Cancer Antigen only, 128 Reproductive Antigen only (524 candidates), 1,120 Other, and 18,249 Unsafe**, backed by the refreshed sequence bundle `oncoref-sequences-1.46.0.zip` (**1,087** allowed reference proteins across all **524/524** candidate groups, **0** `Unsafe` allowed).

## Highlights

1. **`target_selection`-backed sequence bundle (`#725`, `#726`):**
   - Rebuilds the downloadable protein sequence bundle (`oncoref-sequences-1.46.0.zip`, `oncoref:sequence-bundle:v1:sha256:938da7bdfa3903f582e086797288199d42d6f449b561fad7fc4780eab4badb9e`) directly from `target_selection` (`oncoref.target-selection.2026-11-20.2`) under `oncoref.target-selection-sequence-roles.2026-11-21.1`.
   - Eliminates all previously allowed `Unsafe` and `Other` records, covers **327/327** `CTA`, **69/69** `Cancer Antigen only`, and **128/128** `Reproductive Antigen only` groups (**570** `direct_target_records`, **650** `cancer_antigen_records`, **731** `reproductive_antigen_records`, and **1,087** `shared_sequence` records), and archives provenance under [sequence-bundle-20261121](audits/sequence-bundle-20261121/index.md).

2. **Gene-identifier `'nan'` corruption fix, complete `data_manifest` registration, and multi-species API consistency (`#727`):**
   - Fixes `gene_ids` (`_canonical_gene_index` / `_canonical_symbol_index`) so missing symbols in `canonical-gene-space.csv.gz` return `None` rather than literal `'nan'` and no longer shadow the `SCN11A` `'NaN'` synonym (`canonical_gene_id("NaN") == "ENSG00000168356"`).
   - Cleans all 51 literal `'nan'` `gene_symbol` cells in `rna-protein-hpa-priors.csv.gz` with read-time schema validation in `oncoref.rna_protein`.
   - Registers every bundled dataset in `oncoref.data_manifest` (`0` unregistered allowlist entries) and aligns `antigen_sets` (`species=` support across `proteins`, `gene_ids`, `sequences`), `oncoref antigens subtract --species`, and top-level `oncoref` submodule exports.

3. **Canonical `target_selection` CTA wiring across coverage, plots, and peptide subtraction (`#728`):**
   - Switches default `cta_gene_ids()` / `cta_gene_id_to_name()` in `coverage.py`, `plots.py`, and `peptides.py` from the frozen legacy `cta` snapshot to canonical `antigen_sets.gene_ids("CTA")`, while subtracting `peptides.py` 9-mers against all `Unsafe` and `Other` proteins (`CA|RA` candidate universe).
   - Extends `coverage.resolve_gene_set()` to all `antigen_sets.CANDIDATE_SETS` tokens (`cta`, `ca`, `ra`, `ca|ra`, `mra`, `fra`, `era`, `sra`, `cta-m`, `cta-f`, `cta-e`, `cta-mf`) plus `legacy-cta`, and preserves row 0 in headerless one-gene-per-line CSV/TSV/TXT panels.

4. **`pirlygenes` plotting/coverage parity and clarified CTA curation & active-policy gate exports (`#607`, `#729`):**
   - Adds `oncoref.figure_export` (`queue_scatter_labels`, `save_figure`) writing `>=300` dpi PNGs paired with vector PDFs (`pdf.fonttype=42`) and deterministic collision-free scatter label placement with data-coordinate leader lines across `figure_style.save`, `plots._save`, and `coverage`.
   - Brings `oncoref.coverage` and `oncoref plot` to full `pirlygenes` parity (`Threshold`, owner `linear_tpm_comparable` auto/TPM/percentile contracts, honest zero-weight `burden_weights(codes) -> (weights, unmapped)`, cDNA-collapsed `covering_set_matrix`, `greedy_covering_set`, `render_covering_set`, and `oncoref plot covering-set`).
   - Clarifies `cta_curation_plots` stage labels and RNA-vs-IHC scope headers, and exports `cta-stage-definitions.json`, `cta-final-stage-attrition.csv`, `active-policy-gate-matrix.csv` (explicit per-gate pass/fail/unassessed/not_applied states over all 19,893 human protein groups), and `active-policy-transitions.csv`.
