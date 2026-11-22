# Releases and Archive

OncoRef pins every dataset, classification rule, and output artifact by SHA-256 digest so analyses remain reproducible across package versions.

## Recent releases

- [Release 1.47.0](cta-release-1.47.0.md) — Multi-species Ensembl ID & schema helpers and audited partial technical-compartment matrices (#732, #733, #737), structured legacy `oncoref.cta_*` deprecations and `Disposable`/`Silent` reconciliation (#647, #730), offline `sequence_bundle` ↔ `PeptideIndex`/`antigen_sets` bridge and near-self subtraction contract (#640, #731), and hard-link-aware bounded cache inventory & pruning (#591).
- [Release 1.46.0](cta-release-1.46.0.md) — `target_selection`-backed `sequence_bundle` (`oncoref-sequences-1.46.0.zip`, #725), `gene_ids` `'nan'` / `SCN11A` `'NaN'` synonym fix and `rna-protein-hpa-priors` cleanup (#727), complete `data_manifest` registration and multi-species `antigen_sets` API consistency (#727), canonical `target_selection` CTA wiring across `coverage`, `plots`, and `peptides` (#728), and full `pirlygenes` plotting/coverage parity plus clarified CTA curation & active-policy gate exports (#607, #729).
- [Release 1.45.0](cta-release-1.45.0.md) — Declarative multi-species target selection, first-class sex-restricted subsets (`MRA`/`CTA-M`, `FRA`/`CTA-F`, `ERA`/`CTA-E`, `SRA`/`CTA-MF`) and per-species whitepapers (#574), complete Ensembl 115 mouse host reference (#651), `SP7` adult bone coverage gate (#707) and `#711` rare-cell Ensembl 115 attribution reviews, GLOBOCAN 2024 cancer-burden refresh and multi-bucket/disjointness APIs (#542, #543, #578), and historical CTA recovery audit (#714).
- [Release 1.44.0](cta-release-1.44.0.md) — Consortium and integrated-atlas donor deduplication in rare-cell RNA, hair-follicle IHC cell-pattern rule, portable reference/tumor fetching (`python -m oncoref.curation fetch`), optional `mean + 2 SD` normal-tissue statistic, and sex-restricted CTA subsets (`MRA`/`CTA-M`, `FRA`/`CTA-F`, `ERA`/`CTA-E`).
- [Release 1.43.1](cta-release-1.43.1.md) — GTEx v11 primary tissue matrix provenance and sample-level reproductive comparisons.
- [Release 1.43.0](cta-release-1.43.0.md) — Adult thymus compartment scoping and organ-loss hierarchy updates.
- [Release 1.42.0](cta-release-1.42.0.md) — Symmetric 10× cancer and reproductive enrichment and explicit patient-profile separation.
- [Release 1.41.0](cta-release-1.41.0.md) · [Release 1.40.0](cta-release-1.40.0.md)

## Technical audits and frozen compatibility snapshots

- [Catalog & curation funnel ledger](audits/target-selection-20261110/rebuild/index.md) · [Filter sweep report](audits/target-filters-20261112/index.md) · [Web summary](audits/target-filters-20261112/report.html)
- [Literature adjudication of normal-tissue expression](audits/target-selection-20261110/rebuild/literature.md) · [Normal-function sweep](audits/target-filters-20261112/normal-function-sweep.md)
- [Reference transcript & protein selection](audits/cta-transcript-selection-20261108/index.md) · [Sequence bundle](sequence-bundle.md)
- [TMB / ICI curation audit](curation-audit.md) · [Tissue injury ontology report](audits/tissue-hazards-1.29.0/index.md)

Frozen legacy modules (`oncoref.cta_selection`, `oncoref.cta`, `oncoref.target_sets`) and earlier dated audit directories under `docs/audits/` remain preserved in the package and site build for backward compatibility without affecting current APIs.
