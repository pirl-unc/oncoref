# Multiplex Index-Hopping & Cross-Sample Bleed QC (`#764`)

Policy `oncoref.multiplex-bleed-qc.2026-11-23.1` · Reference panel `gtex_individual_tissue_tpm:v11` · Orthogonal checks `hpa_single_cell_tissue:v23` + `rna_tissue_hpa:v23`.

## Method

Ultra-abundant tissue-specific transcripts (`>= 1,000 TPM`, such as spermatid `PRM1`/`PRM2`/`TNP1` in testis or `KLK3`/`MSMB`/`KLK2`/`ACP3` in prostate) bleed at `10^-5` to `5 x 10^-3` into unrelated biopsies sequenced on shared multiplexed Illumina flowcells ([Nieuwenhuis et al., *Genome Res* 2020](https://doi.org/10.1101/gr.255380.119)).

A biopsy $b$ in non-source organ $T$ is flagged for source-$S$ bleed if and only if:
1. **Multi-sentinel co-detection:** `>= 3` sentinel markers from $\mathcal{M}_S$ are simultaneously detected at `>= 0.25 TPM` in biopsy $b$.
2. **Stoichiometric proportionality:** Median $\hat{\alpha}_{b,g} = \text{TPM}_{b,g} / \text{median\_TPM}_{S,g} \in [10^{-5}, 5 \times 10^{-3}]$ with $\text{IQR}(\log_{10} \hat{\alpha}_{b,g}) \le 0.5$.
3. **Orthogonal somatic veto:** Independent `HPA v23` single-cell clusters (`< 0.5 nTPM`), non-GTEx `rna_tissue_hpa` (`< 0.5 nTPM`), and sorted immune cells (`< 0.5 TPM`) confirm no true expression in organ $T$.

## Top GTEx Sites with Coordinated Testis Sentinel Bleed

| GTEx Site | Organ | Samples | Flagged | Flagged % | Median Bleed Fraction |
|---|---|---:|---:|---:|---:|
| `Adrenal_Gland` | adrenal gland | 295 | 154 | 52.20% | 0.00011 |
| `Liver` | liver | 262 | 136 | 51.91% | 0.00011 |
| `Esophagus_Muscularis` | esophagus | 561 | 260 | 46.35% | 0.00012 |
| `Colon_Transverse` | colon | 479 | 210 | 43.84% | 0.00012 |
| `Stomach` | stomach | 407 | 174 | 42.75% | 0.00013 |
| `Colon_Sigmoid` | colon | 419 | 176 | 42.01% | 0.00013 |
| `Esophagus_Gastroesophageal_Junction` | esophagus | 403 | 150 | 37.22% | 0.00010 |
| `Artery_Aorta` | blood vessel | 472 | 173 | 36.65% | 0.00016 |
| `Small_Intestine_Terminal_Ileum` | small intestine | 207 | 75 | 36.23% | 0.00017 |
| `Pancreas` | pancreas | 362 | 131 | 36.19% | 0.00008 |
| `Esophagus_Mucosa` | esophagus | 614 | 206 | 33.55% | 0.00015 |
| `Thyroid` | thyroid gland | 684 | 222 | 32.46% | 0.00026 |

## Proteins with Bleed-Inflated GTEx Peripheral Baselines (`n = 2`)

| Protein | Source | Source Peak TPM | Raw Baseline TPM | Top Raw Site | Clean Baseline TPM | Top Non-TGCT Tumor TPM |
|---|---|---:|---:|---|---:|---|
| `PRM1` | testis | 58,750.8 | 8.05 | `Thyroid` | **5.71** | 58.97 (`BRCA`) |
| `PRM2` | testis | 50,820.5 | 7.76 | `Thyroid` | **6.65** | 8.60 (`BL`) |

## Denominator Integrity & Anatomic Boundary Carryover Guardrails

1. **Drop-from-denominator (never zero-fill):** Masked bleed-contaminated biopsies are removed from the donor-replicate pool before computing `gtex_donor_p95` (`compute_bleed_cleaned_donor_summary`). Zero-filling flagged libraries would artificially depress the 95th percentile below true biological background.
2. **Flowcell index-hopping vs. anatomic boundary carryover (`SEMG2`, `TGM4`):** Seminal-vesicle and prostate transcripts (`SEMG2` at `26,738 TPM`, `TGM4` at `22,107 TPM`) appear in GTEx `Bladder` (`9.52 TPM` raw) from physical urogenital sampling/dissection proximity rather than uniform flowcell-wide hopping (`bleed_qc_applied = False`), and `TGM4` carries independent healthy tissue MS support (`TGM2` shared peptides audited in `#769`), so neither `SEMG2` nor `TGM4` receives bleed masking.
3. **Interaction with `#768` non-admixed tumor witness rule:** Cleaning `PRM1`'s GTEx `Thyroid` bleed ceiling from `8.05 TPM` (`10×` bar `80.5 TPM`) to `5.71 TPM` (`10×` bar `57.1 TPM`, anchored by `Small_Intestine_Terminal_Ileum`) allows its genuine female breast cancer witness (`TCGA-AO-A0JL-01 = 58.97 TPM` in `BRCA`) to clear `10×` independently of admixed `TCGA-TGCT` biopsies, retaining `PRM1` as `CTA` while `PRM2` (`8.60 TPM` max non-TGCT in `BL` vs. `66.5 TPM` bar) reclassifies cleanly to `Reproductive Antigen only`.
4. **Opposite-sex `host-exposure` bleed masking & GTEx `resting_breast` donor-sex stratification (`#776`):** Because `336 / 514` (`65.4%`) GTEx `Breast_Mammary_Tissue` biopsies come from `XY` male cadavers (`RPS4Y1 + DDX3Y + UTY + KDM5D ≥ 5 TPM`, including male autopsy testis contamination up to `961.6 TPM` `PRM1` and `176.7 TPM` `TGM4`), `host-exposure` stratifies `resting_breast` into `female_max_tpm` (`n = 178` `XX` GTEx + female HPA biopsies) and `male_max_tpm` (`n = 336` `XY` GTEx biopsies) and applies `HOST_EXPOSURE_SENTINEL_PANELS` (including the 4-gene post-meiotic chromatin core `PRM1`/`PRM2`/`TNP1`/`TNP2`) across opposite-sex GTEx shards whenever HPA `scRNA-seq`, HPA normal `IHC`, and coding HPA bulk biopsies are negative—restoring `79` CTAs (`PRM1`, `SPATA19`, `CT83`, `TGM4`, `SEMG2`, `ODF1`, `ADAM2`, `AKAP4`, etc.) from false `CTA-MF` to **`CTA-M`** (`225 CTA-M · 39 CTA-E · 3 CTA-F · 64 CTA-MF`).

## Programmatic Inspection

```python
from oncoref import normal_library_qc, target_selection

# Site-level GTEx stoichiometric sentinel bleed incidence
site_summary = normal_library_qc.multiplex_bleed_tissue_summary()

# Per-protein baseline adjustments and active catalog provenance
deltas = normal_library_qc.multiplex_bleed_gene_deltas()
prm = target_selection.results().query("protein_group_label in ['PRM1', 'PRM2', 'SEMG2', 'TGM4']")[
    [
        "protein_group_label",
        "classification",
        "normal_reference_tpm",
        "tumor_max_tpm",
        "tumor_top_cancer_code",
        "multiplex_bleed_764_disposition",
    ]
]
```
