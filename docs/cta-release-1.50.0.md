# Release 1.50.0

Policy `oncoref.target-selection.2026-11-20.2` · Sequence-role policy `oncoref.target-selection-sequence-roles.2026-11-21.1` · Species policy `oncoref.species-screening-policy.2026-11-21.1` · Murine comparison panel `oncoref.mouse-grcm39-tpm.2026-11-22.1` · Integrated human preview `oncoref.target-selection.integrated-preview.2026-11-23.1`.

Across **19,893** human protein-equivalence groups (**20,099** Ensembl 115 coding loci), `evaluate_species_catalog` reproduces **100% (`19,893 / 19,893`)** of the active human catalog classifications (**327 CTA, 69 Cancer Antigen only, 128 Reproductive Antigen only [524 candidates], 1,120 Other, and 18,249 Unsafe**) and **100% (`455 / 455`)** sex-restriction parity (`MRA`/`CTA-M`, `FRA`/`CTA-F`, `ERA`/`CTA-E`, `SRA`/`CTA-MF`) through the same declarative `SpeciesScreeningSpec` evaluator used for mouse (**21,433** protein groups) and dog, while keeping the built PyPI wheel under **20 MB** (`~15.9 MB`).

## Highlights

1. **Decoupled historical audit manifest SHA-256 checks from live runtime sources (`#753`, `#757`):**
   - Snapshots the **137** frozen source and script inputs referenced by historical `docs/audits/*/run-manifest.json` ledgers under `tests/fixtures/release-1.49-inputs/`, updating `oncoref.curation` (`_matches_live_or_snapshot`) and audit tests so historical provenance digests remain verified without freezing active runtime modules against refactoring.

2. **Unified runtime species selection across human, mouse, and dog (`#754`, `#758`):**
   - Eliminates `taxon_id == 9606` runtime branching across `species_target_selection`, `host_context`, `target_filters`, `antigen_sets`, and `antigen_cli`, backed by `oncoref/data/species/human/target-selection.json` so every species resolves catalog tables, patient/host context filters, sex-restricted antigen subsets, and CLI projections through a single declarative runtime path.

3. **Unified `datacache` + `OPENVAX_DATA_CACHE` cache roots, offline preflight APIs, and `pandas 3.1+` safety (`#747`, `#761`, `#762`, `#759`):**
   - Consolidates cache-root resolution across `_downloads`, `data_bundle`, `reference_data`, `source_matrices`, `sequence_bundle`, `species_assets`, `load_dataset`, and `cache_retention` onto `datacache.get_cache_root("oncoref", ...)` and `OPENVAX_DATA_CACHE` (ignoring empty or lock-only directories during legacy cache fallback), migrating bundle installation and pruning onto `datacache.install_archive` and `datacache.bundle_store`.
   - Adds `allow_download: bool = True` across `reference_data.ensure`, `reference_data.bundled_or_managed_path`, `species.target_selection`, `species.target_selection_coverage`, and `species.host_exposure`, plus `reference_data.is_cached(..., include_bundled=...)` and `reference_data.preflight(...)` (`#762`), and eliminates `Pandas4Warning` on empty `sex_class_by_eid` mappings under `pandas >= 3.1` (`#761`).

4. **Unified human catalog evaluation onto `evaluate_species_catalog` + integrated v1.49 evidence preview (`#755`, `#760`):**
   - Adds `oncoref.human_screening_inputs` and `scripts/evaluate_human_species_catalog.py`, achieving exact **`19,893 / 19,893`** baseline category/subtype parity and **`455 / 455`** sex-restriction parity when evaluating the full human Ensembl 115 proteome through `evaluate_species_catalog(human_spec, ...)`.
   - Publishes `human-737-tumor-observations.{csv.gz,json}` (**5,494** `THYMCA`, linear `GSE125285` `BCC`/`cSCC` observation rows) and `human-integrated-policy-preview.{csv,json}` (**71** protein transitions across `#711`, `#712`, `#735`, and `#737` evidence layers).

5. **Quarantined frozen legacy `cta*.py` modules under `oncoref._legacy` and modularized `oncoref.expression` (`#756`, `#763`):**
   - Moves **22** frozen pre-`target_selection` CTA modules into `oncoref._legacy/` behind compact top-level `sys.modules` compatibility shims (`oncoref.cta*`) while keeping active shared evidence, literature, sequence-overlap, peptide, and coverage modules at top level.
   - Extracts `expression_universe_deltas`, `expression_qc`, `expression_availability`, and `expression_pan_cancer` from `oncoref.expression` (`7,051` → `~4,170` lines) with full public and test-monkeypatch re-export compatibility.
