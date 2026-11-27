# OncoRef 1.55.0 — Multi-Species Sequence & Subtraction Speedups, `HostContext.for_host`, HPA v25.1 Single-Cell `nCPM`, and Species Catalog Polish

OncoRef `1.55.0` delivers multi-species `pyensembl` sequence resolution and `>15×` background-build speedups, promotes `HostContext.for_host(...)` and cross-species sex-restriction helpers with explicit pregnancy scope controls, adds in-process HPA cancer prevalence caching with canonical `symbol` attachment plus `oncoref.dog` and canine/murine `representative_transcripts`, and restores HPA v25.1 single-cell `nCPM` compatibility.

## Highlights

1. **Multi-species sequence resolution, strict offline `auto` fallback, and `>15×` build speedups:**
   - `antigen_sets._sequences()` now resolves each species' pinned Ensembl release dynamically from `packaged_species_metadata(species)["reference"]` (`Ensembl 115` for human/mouse; `Ensembl 98` / `CanFam3.1` for dog), and `PeptideIndex.from_sequence_bundle(..., species=...)` supports `"human"`, `"dog"`, `"mouse"`, and custom non-human `SequenceBundle` archives without `pyensembl`.
   - `sequence_source="auto"` probes `allow_download=False` before falling back to local `pyensembl` (never triggering implicit asset downloads), backed by `sequence_bundle.is_cached()`.
   - Replaced ~18.4k per-gene `genome.gene_by_id(gid).transcripts` SQLite round-trips with a single bulk `CDS` query joined against `protein_sequences.fasta_dictionary` (**76.4s → 4.7s** on `230,575` human `Unsafe` isoforms) and replaced pandas `groupby("sequence")` + `.to_dict("records")` in `PeptideIndex` with single-pass index grouping (**>35s → 1.1s**).

2. **`HostContext.for_host`, cross-species `compartments()`, `SEX_RESTRICTION_CODES`, and pregnancy scope controls:**
   - Promoted `host_context.from_legacy(...)` to `HostContext.for_host(...)` and `host_context.for_host(...)` (supporting `max_tpm`, `max_ntpm`, `negative_ihc_multiplier`, and `allow_incomplete_scope`, with `from_legacy` preserved as a warning-free alias).
   - Exported public `host_context.compartments(species="human")` and `SEX_RESTRICTION_CODES = ("CTA-M", "CTA-F", "CTA-MF", "CTA-P", "RA-M-only", "RA-F-only", "RA-MF-only", "RA-P-only")` (also re-exported on `oncoref.antigens`).
   - Kept `pregnant=True` fail-closed (`0` eligible) by default across human, dog, and mouse because `embryo` carries `scope_complete=False` (integrated early-embryo atlases lack 1st–3rd trimester fetal organogenesis coverage), documented maternal-fetal interface screening via `protect=("placenta", "lactating_breast")` (`scope_complete=True`, `306 / 322` human CTAs eligible), and added explicit `allow_incomplete_scope=True` (`--allow-incomplete-scope`) across `HostContext`, `TargetFilter`, `antigens.*`, `antigen_sets.*`, and `oncoref antigens` CLI to screen measured early-embryo + placenta + lactating RNA while still failing closed on unmeasured or `>= 8` TPM/nTPM exposure.

3. **HPA cancer prevalence caching + canonical `symbol`, `oncoref.dog` symmetry, and dog/mouse `representative_transcripts`:**
   - `hpa_cancer_rna_prevalence()`, `hpa_cancer_ihc_prevalence()`, and `hpa_cancer_rna_ihc_comparison()` now cache parsed tables in-process keyed by `(path, mtime_ns, size)` (returning `.copy()` and registering with `load_dataset._clear_cache()`) and attach canonical OncoRef `symbol` immediately after `gene_id` (`include_symbol: bool = True`).
   - Exposed `oncoref.dog` (`import oncoref.dog` / `from oncoref import dog`) as a symmetric module alias for `oncoref.canine` alongside `oncoref.human` and `oncoref.mouse`.
   - Deterministically populated `representative_transcripts` across all `20,082` dog (`ENSCAFT*`) and `21,433` mouse (`ENSMUST*`) `target_selection()` rows offline from compact SHA-256-verified wheel sidecars and `group_reference_proteins()`.

4. **HPA v25.1 single-cell `nCPM` compatibility & provenance unit:**
   - Updated `hpa_cell_type_expression(version=None)`, `gene_cell_type_ntpm(gene_id, *, version=None)`, and `normal_cell_data.observations("hpa_single_cell")` to resolve both HPA v25.1's `nCPM` column and historical v23's `nTPM` column (recording `attrs["value_column"]` and returning defensive `.copy()` frames), and updated `hpa_single_cell` `v25.1` `"unit"` to `"nCPM"` in `reference-releases.json`.
