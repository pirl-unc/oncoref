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
