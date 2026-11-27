# OncoRef 1.54.0 — Unified `oncoref.antigens` API, Live Multi-Species Explorer, and Zero-Sidecar Catalog Architecture

OncoRef `1.54.0` simplifies programmatic, command-line, and interactive exploration of the antigen catalogs across **human (`9606`)**, **domestic dog (`9615`)**, and **laboratory mouse (`10090`)**, eliminates runtime catalog sidecar patching and hardcoded entity exceptions, and removes internal issue-number and acronym clutter across the documentation and research whitepapers.

## Highlights

1. **Single, self-explanatory `oncoref.antigens` Python API & `oncoref antigens` CLI:**
   - `oncoref.antigens` exposes seven plain-keyword functions (`table`, `gene_ids`, `explain`, `summary`, `funnel`, `sequences`, `subtract`) across `"human"`, `"dog"`, and `"mouse"` without requiring callers to construct internal `HostContext`, `TargetFilter`, `AntigenSelection`, `SpeciesScreeningSpec`, or `MatchPolicy` objects.
   - `sex="female"` (`--sex female`) protects all non-pregnant, non-lactating adult female reproductive and resting mammary organs of the species and enforces `>=10×` tumor fold above protected-organ exposure (`258` retained of `270` `CTA-M + CTA-P` human CTAs), and `sex="male"` (`--sex male`) protects all adult male reproductive and resting mammary organs (`58` retained of `60` `CTA-F + CTA-P` human CTAs), resolving the legacy `--profile female-patient` / `male-patient` trap (which left `protect=()` empty by default).
   - First-class `ectopic_only=True` (`--ectopic-only`), `replicated_only=True` (`--replicated`), and `oncoref antigens summary` make ectopic activation, multi-patient replication (`tumor_support`), and cross-species counts one-line queries.

2. **Live multi-species [Interactive Antigen Explorer](antigen-explorer.html):**
   - Replaces the frozen `2026-11-10` human-only HTML link on the site front door with a deterministic, zero-dependency multi-species interactive explorer (`docs/antigen-explorer.html`, built by `scripts/render_antigen_explorer.py` and verified by `oncoref.curation`) covering Human (`322` CTAs + `CA only` + `RA only` + clinical TCR benchmarks), Dog (`118` CTAs + candidates + 14 literature nominations), and Mouse (`8` CTAs + candidates + TIL/clinical anchors) with one-click physiological presets (`Female-safe`, `Male-safe`, `Adult-absent CTA-P`, `Female-restricted CTA-F`), ectopic/replication/cell-line toggles, expandable evidence cards, and CSV export.

3. **Zero-sidecar catalog materialization & general biological rules:**
   - All clinical organ-clearance (`clinical_clearance_*`), closed homologous repeat-family (`closed_homolog_family_id`, `co_silenced_coding_sibling_ids`), 1-mismatch near-self 9-mer (`exact_unique_9mers`, `safe_2mm_9mers`, `unsafe_near_self_1mm_pct`, `unsafe_shared_tryptic_peptides`), and `tumor_support` (`"replicated"` vs. `"single_sample"`) columns are pre-baked directly into `oncoref/data/target-selection.csv.gz` at build time, eliminating runtime `annotate_catalog()` sidecar mutations on `target_selection.results()` and `packaged_species_results()`.
   - Eliminated hardcoded entity-ID checks and post-hoc `target_selection.results()` reconciliation in `human_screening_inputs.py` and `scripts/promote_integrated_human_catalog.py`, replacing them with general mammary protein-coding threshold and multi-cohort homologous-universe rules across somatic and reproductive tissues.

4. **De-slopped navigation, guides, and research whitepapers:**
   - Cleaned primary site navigation (`mkdocs.yml`) and stripped internal GitHub ticket tags, `"Note R"` jargon, `"sidecar"` framing, and repetitive legacy alias lists across [Home](index.md), [API guide](api.md), [Classification rules](target-selection.md), [Antigen whitepaper](antigen-whitepaper.md), and [Multi-species whitepaper](species-antigen-whitepaper.md).
