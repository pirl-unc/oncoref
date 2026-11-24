# Holistic GTEx Tissue Deconvolution & Donor Quantile Estimation (`#778`)

**Policy ID:** `oncoref.gtex-holistic-deconvolution.2026-11-24.3`  
**Evaluated Scope:** `52` GTEx v11 shards · `18,637` RNA-seq libraries · `946` `GTEX-XXXXX` donors · `904` Broad `SM-XXXX` 96-well prep-plate batches · `25` disjoint organ archetypes (`185` anchor genes, `342` 52-shard NNLS basis genes)

## Four Contamination Regimes Separated from Histological Lineage & Surgical Boundaries

1. **Patterned-flowcell index hopping (`flowcell_index_hopping`, `2×10⁻⁵ – 5×10⁻³`):** Pervasive multiplexed pool bleed across co-sequenced tissues; conservatively subtracted via `alpha_lower_q25 * source_median_tpm` when bleed accounts for `≥25%` of observed TPM.
2. **Cross-donor 96-well plate / extraction-batch pipetting splash (`plate_pipetting_splash`, `≥0.5%`):** Requires cross-donor `SM-XXXX` co-plating (without same-donor co-plating), verified cross-donor `SMNABTCH`/`SMGEBTCH` batch sharing, or chromosomal sex impossibility (`XX` recipient carrying `testis`/`prostate`, or `XY` recipient carrying `ovary`/`uterus`). Heavy splashes (`≥0.5%` from expressing sources) are dropped from the gene's donor denominator rather than zero-filled.
3. **Same-donor prosector carryover vs. surgical/anatomical boundaries (`same_donor_carryover` vs. `anatomical_boundary`):** Separates non-adjacent autopsy bench carryover within `GTEX-XXXXX` from true contiguous surgical planes (`prostate → bladder`, `epidermis → subcutaneous adipose`, `gastrocnemius → tibial nerve`, `pectoralis → breast`, `cardia → stomach`), which are **never** subtracted from somatic safety `N`.
4. **Gross biopsy tube swaps & cross-sex chimeras (`tube_swap`, `≥50%` foreign chimerism, `XX`/`XY` sample-donor sex discordance at `≥20%`, or `≥20%` foreign + joint 52-shard NNLS `native_retention < 0.25`):** Quarantined globally across all genes.

## Quarantined Foreign-Organ Swaps

| Sample ID | Donor Sex | Labeled Shard | Detected Organ | Plate Batch | Foreign $\hat{\alpha}$ | NNLS Native Retention | Detected Anchors |
| --- | --- | --- | --- | --- | ---: | ---: | --- |
| `GTEX-1GZHY-2426-SM-9JGGK` | `XX` | `ovary` | `lung` | `SM-9JGG` | `1.335` | `0.392` | `SCGB3A2;SFTPA1;SFTPA2;SFTPB;SFTPC` |
| `GTEX-11GS4-2826-SM-5HL56` | `XY` | `nerve_tibial` | `pancreas` | `SM-5HL5` | `0.609` | `0.504` | `CELA3A;CLPS;CPA1;CTRB2;PLA2G1B;PNLIP;PRSS1;PRSS2;REG1A;REG1B` |

## Top Cross-Donor Wet-Lab Pipetting Splashes (`plate_pipetting_splash`)

| Sample ID | Donor Sex | Target Shard | Splashed Organ | Plate Batch | SM-Plate Peer | Batch-Matched | Splash % | Detected Anchors |
| --- | --- | --- | --- | --- | --- | --- | ---: | --- |
| `GTEX-1H1DG-0011-R11a-SM-C1YQK` | `XY` | `brain_cerebellar_hemisphere` | `adrenal_gland` | `SM-C1YQ` | `True` | `False` | `43.36%` | `CYP11B1;CYP11B2;CYP17A1` |
| `GTEX-17GQL-2326-SM-731EJ` | `XY` | `artery_aorta` | `prostate` | `SM-731E` | `True` | `False` | `39.09%` | `ENSG00000305268;KLK2;KLK3;MSMB` |
| `GTEX-ZDYS-0426-SM-5IJEQ` | `XY` | `lung` | `stomach` | `SM-5IJE` | `True` | `False` | `32.85%` | `ATP4A;ATP4B;CBLIF;GKN1;LIPF;PGA3;PGA5;TFF1;TFF2` |
| `GTEX-ZPIC-1026-SM-59HKZ` | `XX` | `artery_aorta` | `pancreas` | `SM-59HK` | `True` | `False` | `25.75%` | `CELA3A;CLPS;CPA1;CTRB2;PLA2G1B;PNLIP;PRSS1;PRSS2;REG1A;REG1B` |
| `GTEX-ZDXO-2426-SM-5S2NJ` | `XY` | `adipose_subcutaneous` | `pancreas` | `SM-5S2N` | `False` | `True` | `24.88%` | `CELA3A;CLPS;CPA1;CTRB2;PLA2G1B;PNLIP;PRSS1;PRSS2;REG1A;REG1B` |
| `GTEX-145LU-0126-SM-5TDCB` | `XY` | `skin_sun_exposed_lower_leg` | `pituitary` | `SM-5TDC` | `True` | `False` | `23.78%` | `CGA;FSHB;GH1;GHRHR;LHB;MIR7-3HG;POMC;POU1F1;PRL;TSHB` |
| `GTEX-14JG6-0005-SM-9JGFV` | `XX` | `whole_blood` | `lung` | `SM-9JGF` | `True` | `False` | `12.49%` | `SFTPA1;SFTPA2;SFTPB;SFTPC` |
| `GTEX-1E2YA-3126-SM-9JGI6` | `XY` | `pituitary` | `esophagus_mucosa` | `SM-9JGI` | `True` | `False` | `11.45%` | `CNFN;CRNN;KRT13;KRT4;KRT6A;SPRR1A;SPRR1B;SPRR2A;SPRR3;TGM3` |
| `GTEX-14A5H-0526-SM-5TDCN` | `XY` | `adipose_subcutaneous` | `pituitary` | `SM-5TDC` | `True` | `False` | `11.39%` | `CGA;FSHB;GH1;GHRHR;LHB;MIR7-3HG;POMC;POU1F1;PRL;TSHB` |
| `GTEX-1H3VY-1826-SM-9JGI4` | `XY` | `adrenal_gland` | `minor_salivary_gland` | `SM-9JGI` | `True` | `False` | `11.03%` | `BPIFB2;C6orf58;FDCSP;LPO;MUC5B;MUC7;STATH;ZG16B` |
| `GTEX-1I1GR-1426-SM-CY8H8` | `XY` | `esophagus_muscularis` | `testis` | `SM-CY8H` | `False` | `True` | `10.72%` | `BOD1L2;C16orf82;LELP1;ODF1;PRM1;PRM2;SMCP;TNP1;TSACC;TUBA3C` |
| `GTEX-1J8EW-1826-SM-C1YQW` | `XX` | `adrenal_gland` | `brain_cerebellar_hemisphere` | `SM-C1YQ` | `True` | `True` | `10.54%` | `BARHL2;GABRA6;GRM4;KCNJ9;NEUROD1;NEUROD2;ZIC2;ZIC4` |
| `GTEX-XV7Q-2626-SM-4BRVA` | `XX` | `adipose_subcutaneous` | `vagina` | `SM-4BRV` | `True` | `False` | `10.50%` | `SERPINB3;SERPINB4;TMPRSS11D` |
| `GTEX-1RB15-0126-SM-G9ERT` | `XY` | `spleen` | `brain_nucleus_accumbens_basal_ganglia` | `SM-G9ER` | `True` | `False` | `9.78%` | `CAMKV;NTSR2;PDYN;PRKCG;SYNPR` |
| `GTEX-1GF9W-2726-SM-9JGGW` | `XY` | `brain_cerebellum` | `lung` | `SM-9JGG` | `True` | `False` | `9.77%` | `SCGB3A2;SFTPA1;SFTPB;SFTPC` |

## Programmatic Inspection

```python
from oncoref import gtex_deconvolution

receipt = gtex_deconvolution.gtex_deconvolution_receipt()
swaps = gtex_deconvolution.gtex_quarantined_swaps()
summary = gtex_deconvolution.gtex_deconvolution_tissue_summary()
print(receipt['anchor_source_archetypes'], receipt['quarantined_swap_samples'], len(summary))
```
