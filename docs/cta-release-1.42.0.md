# OncoRef 1.42.0: CTA compartment tiers, lactating mammary separation, and host toxicity profiles

The default `CTA` catalog (**290 CTA, 106 Cancer Antigen only, 113 Reproductive
Antigen only, 1,131 Other and 18,253 Unsafe proteins** across **19,893** human
protein groups) intentionally spans all 12 reproductive, developmental, and
mammary compartments so no sex-specific, lactational, or surgically absent
target is discarded upstream. OncoRef 1.42.0 formalizes the biological hierarchy
inside those 12 compartments and lets callers slice candidate sets directly by
patient anatomy and tissue-loss tolerance.

## Highlights

- **Lactating vs. resting mammary separation (`lactating_breast` vs. `resting_breast`):**
  Individual HPA v23 breast biopsies expressing a coordinated milk-protein program
  (`>= 3` of `LALBA, CSN1S1, CSN2, CSN3` at `>= 100` native TPM; `breast.410` in HPA v23)
  are separated into `lactating_breast`, distinct from `resting_breast` (`breast.373`,
  `breast.390`, `breast.405`). Because men possess rudimentary mammary ducts,
  `resting_breast` is protected in both sexes, while milk-program targets (`LALBA`,
  `BTN1A1`) remain available in non-lactating women and men and are immediately
  excluded when `--lactating` (`lactating=True`) is enabled.
- **Three-way biological `cta_tier` decomposition:** Every CTA record carries an
  explicit `cta_tier` (`257 core + 9 core_with_extended_exposure + 24 extended_only = 290`)
  alongside per-compartment maxima (`reproductive_compartment_max_tpm`) and qualifying
  origins (`reproductive_qualifying_compartments`, `reproductive_qualifying_tpm`).
- **Unified host / patient toxicity profiles:** `antigen_sets.proteins`,
  `antigen_sets.gene_ids`, `target_filters.TargetFilter`, `target_filters.profile_summary`,
  and `oncoref antigens list` / `explain` / `funnel` support named profiles (`all` **290**,
  `female-patient` **285**, `female-patient --lactating` **283**, `male-patient` **260**,
  `core` **257**) plus `--pregnant`, `--allow-tissue` (e.g., `prostate`, `seminal vesicle`
  post-prostatectomy -> **264** CTAs), and `--protect-tissue` (e.g., `testis`, `ovary`).
  Every profile is governed by one rule over its protected compartment set **P**:
  each compartment in **P** must stay `<= 32` TPM, tumor expression must reach
  `>= min_tumor_fold * max(H0, max_P TPM)`, and `RA` must reach
  `>= min_reproductive_fold * H0` in at least one allowed compartment outside **P**.
- **End-to-end source-to-tiers architecture & compartment visualizations:**
  Generated figures (`source-to-tiers-architecture` and `cta-tiers-and-profiles`)
  trace all upstream sources (Ensembl 115 + MANE Select, 16 nomination catalogs,
  7 normal somatic panels, sequence/coverage/bounds checks, 130 QC-passing tumor
  cohorts, and 54 HPA biopsies + Petropoulos/Tyser embryos) through the gradual
  funnel into `CA`, `RA`, `CTA`, and host profiles.
- **Documentation & build DAG cleanup:** Reorganized the front-page architecture
  and navigation hierarchy (`docs/index.md`, `README.md`, `mkdocs.yml`) and added
  `filter_web_report` (`scripts/build_target_filter_web_report.py`) to the
  checksum-verified `oncoref.curation` build graph (`210` verified outputs).

Expression bundle **5.23.29** and raw source matrices **5.22.14** are unchanged.
Screening profiles describe anatomical and physiological compartment rules over the
declared panels; they do not establish clinical safety or absence of off-target
receptor recognition.

[Current rules](target-selection.md) ·
[Host profiles & candidate filters](target-filters.md) ·
[Commands and build workflow](curation.md) ·
[Architecture & summary figures](audits/target-selection-20261110/rebuild/summary.md) ·
[Filterable catalog](audits/target-selection-20261110/rebuild/report.html)
