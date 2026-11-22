# Release 1.47.0

Policy `oncoref.target-selection.2026-11-20.2` · Sequence-role policy `oncoref.target-selection-sequence-roles.2026-11-21.1` · Near-self subtraction contract `oncoref.near-self-subtraction.2026-11-21.1` · Species policy `oncoref.species-screening-policy.2026-11-21.1`.

Across **19,893** human protein-equivalence groups (**20,099** Ensembl 115 coding loci), the catalog classifies **327 CTA, 69 Cancer Antigen only, 128 Reproductive Antigen only (524 candidates), 1,120 Other, and 18,249 Unsafe**.

## Highlights

1. **Multi-species Ensembl ID helpers, canonical schema normalization, and audited partial technical-compartment handling (`#732`, `#733`, `#737`, `#739`):**
   - Adds `is_ensembl_gene_id(..., species=...)` and `ENSEMBL_GENE_PREFIXES` (`ENSG`, `ENSCAFG`, `ENSMUSG`) in `oncoref.gene_ids`, decoupling `coverage`, `drivers`, `expression`, and `protein_annotations` from hardcoded human `ENSG*` assumptions and preventing human symbol fallback from contaminating canine/murine gene panels in `coverage`.
   - Adds `gene_id_column`, `symbol_column`, `cancer_code_column`, and `canonicalize_gene_columns(style='expression' | 'lower', include_cancer_code=...)`, plus `driver_cancer_codes()` and `driver_gene_evidence(..., cancer_code=..., include_cancer_code=...)`.
   - Audits `SCLC` (`sclc_ucologne_2015`) in `source_matrices.partial_technical_compartment_cohorts()` and adds `on_empty_compartment='warn' | 'record'` with `attrs['unfilled_compartments']` in `normalization.clean_tpm`.

2. **Structured legacy CTA deprecation and `Disposable` / `Silent` reconciliation (`#730`, `#647`, `#740`):**
   - Moves top-level `oncoref.cta_*` / `oncoref.CTA_*` exports to lazy `__getattr__` accessors that emit a structured once-per-process `DeprecationWarning` pointing to `oncoref.antigen_sets.gene_ids('CTA')` / `oncoref.target_selection.results()`, while preserving explicit `oncoref.cta.*` imports warning-free for historical provenance reproduction.
   - Adds `oncoref.legacy_api_migrations()` / `legacy_api_migration(surface)` mapping all historical CTA/target API generations onto current surfaces.
   - Adds `normal_expression_categories.reconcile_with_target_selection()` and `adult_disposable_only_proteins()` cross-tabulating legacy `Normal / Reproductive / Other(Disposable, Silent) / Unassessed` concepts onto `target_selection.results()` (`DISPOSABLE_CAP = 32`) without relabeling `Other`, weakening `Unsafe` precedence, or dropping `Other` from default subtraction backgrounds.

3. **Offline `sequence_bundle` ↔ `PeptideIndex` / `antigen_sets` bridge and near-self subtraction contract (`#731`, `#640`, `#741`):**
   - Adds `SequenceBundle.to_sequences_table()` / `sequence_bundle.sequences_table()`, `sequence_source='auto' | 'sequence_bundle' | 'pyensembl'` across `antigen_sets` and `oncoref antigens subtract`, and `PeptideIndex.from_sequence_bundle(...)` so exact and bounded edit-distance subtraction works offline straight from `oncoref-sequences` without `pyensembl` (`0` missing genes across all 19,893 human `target_selection` protein groups).
   - Formalizes the downstream near-self / forbidden-reference contract via `peptide_subtraction.near_self_contract()` and `antigen_sets.forbidden_background(...)` with explicit provenance and uncertainty metadata (`omitted_categories`, `unmeasured_home_entity_ids`, `review_hold_entity_ids`, `host_excluded_entity_ids`, `cta_exemption_applied=False`, and fail-closed species guardrails).

4. **Bounded multi-family cache inventory and hard-link-aware retention pruning (`#591`, `#742`):**
   - Adds `oncoref.cache_retention` (`cache_inventory`, `prune_caches`) and extends `oncoref data prune [all|bundle|source|reference|hpa|sequence|species] [--keep-version V] [--include-legacy-pirlygenes] [--json] [--yes]`.
   - Deduplicates `(st_dev, st_ino)` hard links so `logical_bytes`, `unique_bytes`, and `reclaimable_bytes` report true physical disk savings across overlay bundles while transitively retaining `base_bundle` dependencies, per-cohort `source_matrices` pins + newest local versions, `DEFAULT_HPA_VERSION` + `curation.required_reference_inputs()`, symlinks, custom files, and active `.lock` files.
