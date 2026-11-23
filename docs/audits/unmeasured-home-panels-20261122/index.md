# Quantitative Adult Human Normal Panels for Unmeasured Home Tissues (`#735`)

Pinned to `oncoref.target-selection.2026-11-20.2` (**19,893** human protein-equivalence groups) and `oncoref:tissue-tiers:v1`. Generated reproducibly by [`build_unmeasured_home_panels.py`](https://github.com/pirl-unc/oncoref/blob/main/scripts/build_unmeasured_home_panels.py).

## 1. Scope & Unit-Separation Contract

Standard pan-tissue bulk atlases (`GTEx v11`, `HPA v23`) omit five specialized adult human compartments (`UNMEASURED_HOME_TISSUES`: `olfactory epithelium`, `hair follicle`, `taste bud`, `dental pulp`, `bone`) plus `nasal mucosa`. Across the human proteome, **568** protein-equivalence groups carry one of these compartments as their curated `home_tissue`, of which **95** (`89` `Other:Unscreened tissue` + `6` tenfold-tumor-enriched `Reproductive Antigen only`) are held back from `Cancer Antigen` / `CTA` admission solely by unmeasured home-tissue coverage.

To replace qualitative holdbacks with empirical measurements while preserving strict scale separation:

- **Linear bulk RNA-seq (`TPM`)** is evaluated within-unit against each compartment's `tissue_tiers` ceiling (`4.0 TPM` Vital `bone`, `8.0 TPM` Protected `taste bud` / `nasal mucosa`, `16.0 TPM` Manageable `olfactory epithelium` / `dental pulp`, `32.0 TPM` Disposable `hair follicle`) and against `tumor_max_tpm` for `10×` tumor-over-home enrichment (`tumor_max_tpm >= 10 × max(normal_reference_tpm, home_bulk_p95_tpm, 1.0)`).
- **Donor-resolved single-cell RNA-seq (`UMI`)** is filtered to QC-passing droplets (`>= 500 UMI`) and evaluated strictly on multi-molecule (`>= 2 UMI`) donor detection percentages — never converted or pooled into bulk `TPM`.
- **Targeted receptor qPCR copy numbers (`Verbeurgt et al. 2014`)** are retained as an auxiliary 26-donor olfactory receptor presence layer and never converted to `TPM`.

---

## 2. Checksum-Pinned Adult Human Normal Panels (`11` Panels)

| Home Tissue | Injury Tier | Panel ID | Accession | Modality & Unit | Donors | QC Cells | Ensembl 115 Genes |
|---|---|---|---|---|---:|---:|---:|
| `olfactory epithelium` | `manageable` | `gse80249_olfactory_bulk_tpm` | `GSE80249` | Bulk RNA-seq (`TPM`, GRCh38.p13) | `4` | — | `19,297` |
| `olfactory epithelium` | `manageable` | `gse139522_olfactory_scrnaseq_umi` | `GSE139522` | 10x scRNA-seq (`UMI`, Patient 2 & 3) | `2` | `23,535` | `19,658` |
| `olfactory epithelium` | `manageable` | `verbeurgt_2014_olfactory_qpcr_copies` | `pone.0096333.s004` | TaqMan qPCR (`TaqMan_copies`) | `26` | — | `335` |
| `nasal mucosa` | `protected` | `cellxgene_adult_nose_umi` | `CELLxGENE-2025-11-08` | Adult nasal scRNA-seq (`UMI`) | `26` | `110,020` | `366` |
| `hair follicle` | `disposable` | `gse129611_hair_follicle_scrnaseq_umi` | `GSE129611` | Anagen scalp Drop-seq + 10x (`UMI`) | `5` | `6,091` | `19,338` |
| `hair follicle` | `disposable` | `gse193269_hair_follicle_scrnaseq_umi` | `GSE193269` | Normal scalp 10x `HF1`–`HF3` (`UMI`) | `3` | `19,259` | `19,338` |
| `taste bud` | `protected` | `gse221229_fungiform_papillae_bulk_tpm` | `GSE221229` | Fungiform papillae bulk (`TPM`) | `7` | — | `19,297` |
| `dental pulp` | `manageable` | `gse305125_dental_pulp_bulk_tpm` | `GSE305125` | Central & odontoblast pulp bulk (`TPM`) | `6` | — | `19,844` |
| `dental pulp` | `manageable` | `gse164157_dental_pulp_scrnaseq_umi` | `GSE164157` | Healthy molar pulp 10x `Pulp1`–`5` (`UMI`) | `5` | `25,680` | `19,662` |
| `bone` | `vital` | `gse344036_mandible_bone_bulk_tpm` | `GSE344036` | Healthy adult mandible bone (`TPM`) | `6` | — | `17,474` |
| `bone` | `vital` | `gse147390_primary_osteoblast_scrnaseq_umi` | `GSE147390` | Primary `ALPL+`/`CD45-` osteoblasts (`UMI`) | `1` | `9,801` | `19,658` |

Canonical lineage markers confirm anatomical specificity and row alignment across panels:

- **Hair shaft / cuticle keratins (`KRT85`, `KRT31`):** detected at `>= 2 UMI` across `2` studies and `7/8` scalp hair-follicle donors (`21.74%` and `4.73%` max donor multi-molecule cells), vs `0.00%` in olfactory, dental pulp, and bone single cells.
- **Osteoblast & bone matrix markers (`ALPL`, `SP7`):** `600.74 TPM` and `192.99 TPM` in healthy human bone bulk (`GSE344036`), and `31.36%` and `1.60%` multi-molecule cells in primary human osteoblasts (`GSE147390`).
- **Olfactory neuroepithelium markers (`GFY`, `OMP`):** `394.82 TPM` and `139.40 TPM` in adult olfactory mucosa bulk (`GSE80249`), with `GFY` at `2.55%` multi-molecule cells in `GSE139522` olfactory neuroepithelium and `<= 0.10 TPM` in bone and taste papillae.

---

## 3. Quantitative Adjudications Across the 95 `CA`/`CTA`-Gated Proteins

1. **Active non-Disposable expression exceeding tier safety ceiling (`5` proteins — `home_expressed_above_tier_ceiling`):**
   - `SP7` (`bone`, Vital ceiling `4.0 TPM`): **`192.99 TPM`** in healthy adult human bone (`GSE344036`) and `1.60%` multi-molecule osteoblasts (`GSE147390`), leaving only `3.06×` osteosarcoma (`SARC_OS`) fold over normal bone.
   - `PANX3` (`dental pulp` home, cross-panel Vital `bone` violation): while low in extracted soft dental pulp (`0.20 TPM`), `PANX3` is expressed at **`54.68 TPM`** in healthy adult human bone (`GSE344036`, `13.7×` above the `4.0 TPM` Vital bone ceiling).
   - `BPIFB3` (`2,695.35 TPM`, `3.62%` scUMI), `GFY` (`394.82 TPM`, `2.55%` scUMI), and `GLYATL3` (`41.51 TPM`) in adult `olfactory epithelium` (`> 16.0 TPM` Manageable ceiling).
2. **Active home expression below ceiling that fails `10×` tumor-over-home enrichment (`8` proteins — `home_expressed_insufficient_tumor_fold`):**
   - `PABPC3` (`Reproductive Antigen only`, `13.25 TPM` in adult olfactory mucosa vs `69.66 TPM` tumor peak = `5.26× < 10×`), plus seven olfactory receptors expressed across adult olfactory biopsies (`OR11A1` `14.15 TPM` [`2.06×`], `OR7E24` `7.67 TPM` [`3.27×`], `OR2J3` `7.55 TPM` [`2.70×`], `OR6F1` `5.44 TPM` [`4.15×`], `OR56A3` `3.12 TPM` [`9.14×`], `OR4E2` `2.94 TPM` [`6.10×`], `OR11H4` `1.42 TPM` [`9.89×`]), all independently confirmed across `12–26 / 26` whole-olfactory-mucosa donors in Verbeurgt et al. 2014.
3. **Disposable hair-follicle compartment localization confirmed (`2` proteins — `home_disposable_follicle_expressed_tenfold`):**
   - `KRTAP21-2` (`0.43%` multi-molecule follicle cells, `33.77 TPM` tumor peak) and `KRTAP9-7` (`247.12 TPM` tumor peak) are directly observed in normal human scalp hair-follicle UMI libraries while clearing non-follicular somatic skin.
4. **Detected in home panel below ceiling and passing `10×` tumor-over-home enrichment (`50` proteins — `home_cleared_panel_passes_tenfold`):**
   - `49` olfactory receptors (including `CTA`-candidate `Reproductive Antigen only` receptors `OR10J1` [`22.9×`], `OR1N1`, and `OR4M2` [`30.2×`]) and `1` sweet-taste receptor (`TAS1R2`, `0.58 TPM` in fungiform papillae vs `12.23 TPM` tumor peak = `11.2×`).
5. **Mapped but undetected in finite sampled biopsies (`30` proteins — `home_undetected_in_sampled_panel`):**
   - `12` ultra-high-sulfur hair-shaft cortex `KRTAP*` proteins (depleted by enzymatic single-cell dissociation of cornified hair shafts), `10` sparse monogenic `OR*`/`TAAR9` receptors, `6` bitter-taste `TAS2R*` receptors (diluted in whole fungiform papillae), and `DSPP` (mineralized dentin-front odontoblast transcript, `0.06 TPM` in extracted soft pulp). Kept explicitly distinct from confirmed clearance.

## Programmatic Access

- `oncoref.unmeasured_home_tissue_panels(home_tissue=..., modality=...)`
- `oncoref.unmeasured_home_protein_adjudications(home_tissue=..., disposition=..., only_ca_cta_gated=...)`
- `oncoref.unmeasured_home_gene_expression(gene=...)` (externalized on `reference-unmeasured-homes-20261122` via `oncoref.reference_data`)
- `oncoref.unmeasured_home_tissue_summary()`
