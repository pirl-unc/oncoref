# Release 1.45.0

Policy `oncoref.target-selection.2026-11-20.2` · Species policy `oncoref.species-screening-policy.2026-11-21.1` · Tissue-tier update `2026-11-20.1`.

Across **19,893** human protein-equivalence groups (**20,099** Ensembl 115 coding loci), the catalog classifies **327 CTA, 69 Cancer Antigen only, 128 Reproductive Antigen only (524 candidates), 1,120 Other, and 18,249 Unsafe**.

## Highlights

1. **Multi-species target-selection framework, sex-restricted subsets, and per-species whitepapers (#574, #724):**
   - `oncoref.species_policies` and `oncoref.species_target_selection` separate the universal mathematical definitions of `CA`, `RA`, `CTA = CA ∩ RA`, and sex-restricted subsets (`MRA`/`CTA-M`, `FRA`/`CTA-F`, `ERA`/`CTA-E`, `SRA`/`CTA-MF`) from declarative per-species anatomy, dataset, and unit-threshold parameterizations across **human (`9606`)**, **domestic dog (`9615`)**, and **laboratory mouse (`10090`)**.
   - Promotes `MRA`, `FRA`, `ERA`, `SRA`, `CTA-M`, `CTA-F`, `CTA-E`, and `CTA-MF` (`classify_sex_restriction`) into first-class runtime (`antigen_sets.select`, `human.sex_restricted_antigens`) and CLI (`oncoref antigens list --set CTA-M`) sets parameterized by comparative mammalian anatomy (`normal-tissue-taxonomy.yaml`), failing closed when opposite-sex organs are unmeasured.
   - Exposes `screening_policy()`, `screening_datasets()`, `sex_restricted_antigens()`, `target_selection()`, `target_selection_coverage()`, `host_compartments()`, and `evaluate_target_catalog()` uniformly on `oncoref.human`, `oncoref.canine`, and `oncoref.mouse`, with bidirectional `HostContext` species validation and uniform `policy` / `target_policy_id` columns.
   - Publishes per-species sections in the [main antigen whitepaper](antigen-whitepaper.md) and the companion [multi-species antigen whitepaper](species-antigen-whitepaper.md).

2. **Complete pinned Ensembl 115 mouse host reference and `P1A` qualification (#651, #722):**
   - Registers SHA-256-pinned mouse Ensembl 115 GRCm39 (`GCA_000001635.9`, **66,668** annotated proteins across **57,314** distinct sequences) under `mouse.references.registered_references()` (`ensembl115`) and maps `P1A` (`Trap1a`) to `ENSMUSG00000051257` (`ENSMUST00000054534.4` / `ENSMUSP00000059448.4`, 224 aa) while retaining original DBA/2 `P815` (`M36387.1` / `AAA40475.1`) and NCBI Gene `22037` (`NM_011635.2` / `NP_035765.1`) provenance.

3. **Unscreened adult bone home gate (`SP7`, #707) and rare normal-cell attribution ledger (#711, #723):**
   - Adds `"bone"` to `UNMEASURED_HOME_TISSUES` with primary adult osteoblast (`10.18632/aging.203452`) and human `SP7` bone-turnover (`10.1038/s41467-022-28318-4`) provenance, moving `SP7` from `Cancer Antigen only` (`70 → 69`) to `Other:Unscreened tissue`.
   - Adds `oncoref/data/rare-cell-attribution-reviews.json` adjudicating all 15 prioritized `Other:Review` candidates against Ensembl 115 exon overlap, 3'-polyA distance, and 31-mer paralog sharing (`RGSL1` and `OOSP4B` → `Unsafe`; `DAZL` host confounder rejected while single-study review rule holds; 12 genuine 3'-host / paralog / single-study holds documented).

4. **GLOBOCAN 2024 cancer-burden refresh, multi-bucket scope, and disjointness APIs (#542, #543, #578, #721):**
   - Refreshes world incidence and mortality reference data to GLOBOCAN 2024 (Sung et al. 2026, `PMID:42417444` / `DOI:10.3322/caac.70090`, `20,636,778` cases and `9,762,507` deaths) with extracted counts, ICD-10 site provenance, and deterministic rebuild script `scripts/build_cancer_burden_reference.py`.
   - Adds `burden_scope()`, `burden_categories()`, `cohorts_for_burden_category(..., include_mixed=...)`, `burden_category_overlaps()`, and `burden_categories_disjoint()` for mixed-site cohorts (`CHOL`, `BTC`, salivary `head_and_neck`) and overlapping world site rows (`anus` vs. `colorectal`, `other_and_unknown_primary` vs. minor sites), and removes `os.path.commonprefix` symbol-string calls (#578).

5. **Historical CTA recovery sensitivity audit (#714):**
   - Archives the checksum-pinned sensitivity audit of historical `1.38`/`1.40` transitions and bounded recovery counterfactuals under [cta-recovery-20261117](audits/cta-recovery-20261117/index.md).
