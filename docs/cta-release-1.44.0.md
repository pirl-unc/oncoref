# Release 1.44.0

Policy `oncoref.target-selection.2026-11-20.1` · Tissue-tier update `2026-11-20.1` · Comparison `oncoref.native-gene-tpm-normal-reference.2026-11-20.1`.

Across **19,893** human protein groups (**20,099** Ensembl 115 coding loci), the catalog classifies **327 CTA, 70 Cancer Antigen only, 128 Reproductive Antigen only, 1,121 Other, and 18,247 Unsafe**.

## Curation and data-source upgrades

1. **Consortium and integrated-atlas donor provenance deduplication (#715):**
   `target_selection.resolve_rare_cell_provenance()` resolves CELLxGENE Census collections belonging to shared consortia (Baylor/Chen eye, Sanes eye, Allen reference brain, Michigan kidney, Teichmann heart) and donors embedded inside integrated atlases (HLCA, Baylor/Chen eye, SEA-AD) back to their original study group and biological donor ID before grading multi-study replication.
   - **Catalog impact (3 proteins):** `RBM46` (`Unsafe` → `Other:Review`), `SELENOV` (`Unsafe` → `Other:Review`), and `AADACL4` (`Unsafe` → `Other:Expression`). Strong single-consortium signals stay held in `Other:Review` pending independent confirmation; zero candidates are auto-rescued into `CTA`, `CA`, or `RA`.
2. **Hair-follicle cell-pattern revision and normal-protein IHC replay (#687):**
   Repaired the truncated `skin_follicular` regex in `oncoref/data/normal-tissue-tier-updates.yaml` (`\b(hair follicles?(?: cells?)?|hair matrix|hair shaft|hair root|cells? in (?:cortex/medulla|cuticle|(?:external|internal) root sheath))\b`) and extended `scripts/rebuild_revised_tissue_gates.py` to replay cell-pattern revisions against `hpa-protein-evidence.csv.gz` (`620` positive IHC rows cleanly reclassified to `disposable`; zero unsafe-gene reversals).
3. **Upstream anatomical normalization:**
   Normalized `rna_tissue_gtex:v23` `"hippocampus"` to `"hippocampal formation"` in `scripts/rebuild_target_catalog.py`.
4. **Uniform, portable curation data fetch and normal-tissue summary option:**
   Added `python -m oncoref.curation fetch [--include-tumors]` (`curation.required_reference_inputs()`, `curation.required_tumor_matrices()`, `curation.fetch_inputs()`) and updated `scripts/build_host_exposure.py` and `scripts/rebuild_target_catalog.py` so every external reference file (82 files) and per-sample tumor matrix (132 cohorts) resolves and SHA-256-verifies through `oncoref.catalog` / `oncoref.reference_data` with no machine-local paths or out-of-band downloads. Added `--normal-stat mean_plus_2std` (and `target_selection.summarize_normal_tissue_tpm(..., method="mean_plus_2std")`) as an optional alternative to the default within-tissue `p95` normal baseline.

## Whitepaper and documentation redesign

- Rewrote the [antigen catalog whitepaper](antigen-whitepaper.md) with natural cancer-immunology terminology, crisp opening exposition on the five categories (`Unsafe`, `Other`, `CA`, `RA`, `CTA`) and why OncoRef combines bulk RNA, sorted immune subsets, single-cell/rare-cell UMI, normal IHC, and curated anatomy, plus successive-depth organization.
- Added two standalone conceptual diagrams (`classification-overview.png` / streamlined `source-to-tiers-architecture.png` and `normal-evidence-layers.png`), two sex-restriction figures (`sex-restricted-ra-cta.png`, `sex-restricted-landscape.png`, plus `sex-restricted-antigens.csv` partitioning `RA` and `CTA` into `MRA`/`CTA-M`, `FRA`/`CTA-F`, `ERA`/`CTA-E`, and `SRA`/`CTA-MF`), and split every cramped side-by-side whitepaper figure into full-width single-panel figures (`cancer-enrichment.png`, `reproductive-enrichment.png`, `literature-sources-counts.png`, `literature-sources-retention.png`, `host-profiles-matrix.png`, `host-profiles-outcomes.png`, `tpm-coverage-somatic.png`, `tpm-coverage-neural-repro.png`).
- Replaced boxy card/callout styling in `docs/stylesheets/extra.css` with quiet scientific typography and subtle palette-matched inline category tags.
