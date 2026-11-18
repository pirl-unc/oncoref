# Reproductive evidence and normal-cell holds

The requested cohorts are frozen at release 1.42.0: 78 CA-only proteins and 31 Other:Review proteins with tumor ≥10×. Current measurements are regenerated from pinned sources; this report does not change admission.

Bulk TPM, HPA cell-cluster nTPM and donor UMI frequencies are separate assays. Maxima include all measured non-reproductive normal tissues on the normal side and all reproductive tissues on the reproductive side. The normal maximum includes Disposable tissues; their separate D value is also retained in the CSV. A cell-cluster maximum is descriptive and cannot be divided by bulk TPM to qualify RA. Whole-organ bulk does not measure every cell population. Embryo TPM is a separate within-embryo/lineage mean; its observed maximum can fail cell/donor support.

Sources: HPA v23 tissue-resolved single-cell clusters and bulk consensus, native-TPM HPA individual tissues and bulk/brain summaries, GTEx v11 and HPA-hosted GTEx, Petropoulos/Tyser embryo observations, and the pinned Census donor UMI panel. The native normal reference also includes sorted immune panels. No single-cell value is silently substituted for missing bulk coverage.

## Codes

- **H**: Same-strand host overlap: resolve spliced/exonic locus attribution.
- **P**: Paralog ambiguity: validate gene-unique sequence or a specific protein assay.
- **S**: Selected finding qualifies in one catalog collection; verify independent study support.
- **F**: Accepted reproductive enrichment remains below 10×.
- **E**: Eye bulk comparison covers retina only; no matched HPA bulk or reliable eye staining.
- **C**: Integrated retina atlas contributes; audit original-study/donor overlap before counting replication.
- **I**: Reliable positive IHC in the implicated organ; inspect cell and antibody specificity.
- **—**: Not measured or not attributable; never a measured zero.

## The 78 CA-only proteins

[Full numeric table, exact contexts and sources](ca-only-78.csv). The accepted RA fold uses the existing native-TPM witness/support rules, not the largest value across these assays.

**A source-selection gap affects 10 of these proteins:** C1orf202, CCDC166, CT47B1, CTAGE1, CTAGE6, GOLGA6L1, SLC35G6, SPANXD, TMEM95, TTLL8. Their reproductive bulk TPM summaries clear the same normal/Disposable/1-TPM reference, but current RA admission uses individual HPA samples and supported embryo means, not GTEx tissue summaries. The CSV exposes this separate `bulk_R_fold`. This is not a cross-unit comparison. Before using these summaries for admission, verify the exact-locus quantification and underlying reproductive samples; a tissue summary does not supply an identified donor witness. [Follow-up](https://github.com/pirl-unc/oncoref/issues/716).

| Protein | Normal bulk TPM | Reproductive bulk TPM | Normal cells nTPM | Reproductive cells nTPM | Observed embryo TPM | Accepted R fold |
|---|---:|---:|---:|---:|---:|---:|
| AADACL3 | 3.77 · skin.226 | 7.25 · placenta.203 | 0.1 · brain / excitatory neurons / c-29 | 2.4 · placenta / extravillous trophoblasts / c-8 | 7.86 · 5 day / trophectoderm | 1.92× |
| ANKRD60 | 1.9 · nucleus ambiguus | 5.08 · Testis | 0.6 · thymus / t-cells / c-18 | 11 · testis / late spermatids / c-2 | 3.31 · 7 day / primitive endoderm | 2.05× |
| C1orf202 | 0.8 · dorsomedial nucleus | 12.4 · Testis | 6.1 · kidney / distal tubular cells / c-10 | 229 · testis / early spermatids / c-3 | 11.1 · 7 day / primitive endoderm | 2.46× |
| CCDC166 | 0.2 · nucleus ambiguus | 11.4 · Testis | 0.8 · eye / rod photoreceptor cells / c-4 | 39.6 · testis / early spermatids / c-8 | 0.408 · Carnegie stage 7 / yolk sac mesoderm | 1.15× |
| CIB3 | 5.19 · bone marrow.97 | 3.66 · testis.282 | 8.9 · colon / distal enterocytes / c-9 | 14.4 · testis / late spermatids / c-11 | 2.94 · 7 day / primitive endoderm | 0.556× |
| CLDN34 | 1.4 · temporo-occipital transitional zone | 3.84 · Testis | 4.7 · liver / hepatocytes / c-18 | 18.6 · testis / early spermatids / c-8 | 1.71 · 3 day / unspecified lineage | 2.23× |
| CPSF4L | 2.2 · cerebellar cortex | 3.04 · ovary.3662 | 3 · eye / bipolar cells / c-3 | 7.4 · ovary / oocytes / c-15 | 21.2 · 3 day / unspecified lineage | 7.26× |
| CSAG2 | 5.69 · lymph node.173 | 4.08 · testis.283 | 1 · lymph node / mixed cell types / c-6 | 0.9 · testis / spermatocytes / c-1 | 18.8 · 3 day / unspecified lineage | 1.96× |
| CSAG3 | 4.28 · lymph node.177 | 10.4 · testis.244 | 0.2 · colon / undifferentiated cells / c-8 | 1.8 · testis / spermatocytes / c-1 | 84 · 4 day / unspecified lineage | 5.25× |
| CT47B1 | 0 · Adipose_Subcutaneous | 14.6 · Testis | 0 · adipose tissue / adipocytes / c-12 | 200 · testis / spermatocytes / c-1 | 9.83 · 3 day / unspecified lineage | 6.49× |
| CTAGE1 | 0.0728 · small intestine.234 | 25 · Testis | 1.3 · brain / inhibitory neurons / c-26 | 227 · testis / early spermatids / c-8 | 0.123 · 7 day / epiblast | 7.7× |
| CTAGE15 | 1.6 · thyroid gland | 10.5 · testis | 0.8 · vascular / t-cells / c-10 | 1.5 · prostate / endothelial cells / c-11 | 0.389 · 6 day / epiblast | 5.33× |
| CTAGE6 | 0.8 · somatomor cortex, precentral gyrus (white matter) | 18.7 · testis | 2.2 · pancreas / exocrine glandular cells / c-10 | 9.4 · testis / early spermatids / c-8 | 0.0914 · 5 day / trophectoderm | 7.7× |
| CTXND2 | 0.9 · paraventricular nucleus | 1.77 · fallopian tube.306 | 0.9 · brain / microglial cells / c-43 | 7.8 · placenta / extravillous trophoblasts / c-25 | 1.68 · 7 day / epiblast | 1.77× |
| DCAF8L1 | 0.0513 · Bladder | 3.11 · testis.349 | 0.4 · skeletal muscle / t-cells / c-5 | 12.5 · testis / spermatocytes / c-1 | 5.81 · 3 day / unspecified lineage | 3.11× |
| DEFB125 | 0.1 · retina | 3.53 · seminal vesicle.386 | 0 · adipose tissue / adipocytes / c-12 | 1.8 · testis / late spermatids / c-11 | — | 3.53× |
| DEFB136 | 1.1 · supramarginal gyrus | 1.46 · epididymis.407 | 0.2 · colon / undifferentiated cells / c-8 | 2.4 · testis / late spermatids / c-9 | 0.195 · 6 day / epiblast | 1.33× |
| ENSG00000284934 | 6.65 · adrenal gland.89 | 4.32 · testis.279 | 3.7 · small intestine / proximal enterocytes / c-7 | 0.3 · endometrium / endometrial stromal cells / c-19 | 19 · 3 day / unspecified lineage | 1.25× |
| EVX2 | 3.33 · colon.67 | 7.6 · Testis | 11.7 · rectum / enteroendocrine cells / c-9 | 22.2 · testis / late spermatids / c-13 | 0.721 · 4 day / unspecified lineage | 1.52× |
| FAM240B | 3.4 · angular gyrus (white matter) | 1.53 · testis.278 | 2.7 · salivary gland / fibroblasts / c-19 | 4.1 · testis / spermatogonia / c-14 | 11.1 · 5 day / unspecified lineage | 2.05× |
| FAM47A | 1.3 · transversal temporal gyrus, anterior | 7.14 · testis.282 | 0 · adipose tissue / adipocytes / c-12 | 17.3 · testis / early spermatids / c-8 | 2.22 · 5 day / epiblast | 5.49× |
| FAM47C | 0.28 · Pituitary | 6.75 · Testis | 0.1 · brain / inhibitory neurons / c-23 | 15.3 · testis / early spermatids / c-8 | 0.354 · 6 day / epiblast | 6.5× |
| FOXR2 | 0.0287 · adipose tissue.86 | 3.53 · testis.244 | 0 · adipose tissue / adipocytes / c-12 | 5.2 · testis / spermatogonia / c-4 | 7.92 · 7 day / primitive endoderm | 3.53× |
| FTHL17 | 0.122 · skin.225 | 7.62 · testis.279 | 0.1 · salivary gland / salivary duct cells / c-12 | 27.5 · testis / spermatocytes / c-1 | 9.48 · 3 day / unspecified lineage | 7.62× |
| GAB4 | 1.96 · spleen.378 | 13.6 · testis.282 | 1 · vascular / schwann cells / c-18 | 50.3 · testis / late spermatids / c-11 | 1.55 · 3 day / unspecified lineage | 6.95× |
| GALP | 1.46 · liver.359 | 0.24 · testis.279 | 0.8 · liver / hepatocytes / c-0 | 9 · endometrium / glandular and luminal cells / c-12 | 1.52 · 3 day / unspecified lineage | 0.164× |
| GFY | 1.8 · vermis | 3 · Testis | 3.1 · pbmc / b-cells / c-8 | 16.2 · ovary / oocytes / c-15 | 10 · 3 day / unspecified lineage | 2.64× |
| GOLGA6L1 | 0.3 · cerebral cortex | 36.5 · testis | 1.1 · eye / muller glia cells / c-9 | 15.1 · testis / early spermatids / c-8 | 0.756 · 3 day / unspecified lineage | 0.0991× |
| GOLGA8K | 5.01 · thyroid gland.250 | 15.8 · testis.278 | 1 · spleen / t-cells / c-5 | 21.7 · testis / early spermatids / c-8 | 2.16 · 3 day / unspecified lineage | 3.15× |
| H2BW1 | 0.4 · paracentral lobule, posterior | 1.28 · Testis | 0.1 · brain / excitatory neurons / c-20 | 8.1 · testis / spermatocytes / c-1 | 1.12 · 4 day / unspecified lineage | 0.893× |
| HOXD12 | 6.32 · colon.68 | 2.35 · Vagina | 18.4 · rectum / enteroendocrine cells / c-9 | 0 · breast / adipocytes / c-23 | 1.36 · 5 day / unspecified lineage | 0.339× |
| HUS1B | 2.2 · nucleus ambiguus | 4.69 · Testis | 5.1 · skin / endothelial cells / c-17 | 30 · testis / spermatocytes / c-1 | 8.02 · 7 day / primitive endoderm | 0.169× |
| IQCA1L | 0.367 · skeletal muscle.340 | 16.7 · testis.282 | 0.7 · lymph node / b-cells / c-4 | 91.7 · testis / early spermatids / c-8 | 0.762 · 3 day / unspecified lineage | 7.26× |
| KLF14 | 2.03 · adipose tissue.86 | 2.01 · Testis | 3.9 · salivary gland / serous glandular cells / c-16 | 8 · placenta / smooth muscle cells / c-23 | 2.06 · 5 day / unspecified lineage | 0.95× |
| LCN9 | 2.46 · cerebral cortex.105 | 2.96 · epididymis.407 | 0.5 · brain / oligodendrocyte precursor cells / c-42 | 0.1 · endometrium / endometrial stromal cells / c-7 | 0.0308 · 7 day / trophectoderm | 1.2× |
| LRRC52 | 5.48 · parathyroid gland.1600 | 40.3 · Testis | 1 · skin / endothelial cells / c-17 | 295 · testis / early spermatids / c-3 | 0.682 · Carnegie stage 7 / emergent mesoderm | 6.61× |
| LRRIQ4 | 1.12 · colon.69 | 7.76 · testis.282 | 2.4 · esophagus / t-cells / c-8 | 75.5 · testis / late spermatids / c-5 | 2.8 · 4 day / unspecified lineage | 6.92× |
| MAGEB17 | 1.85 · Esophagus_Gastroesophageal_Junction | 4.55 · testis.282 | 2.8 · brain / inhibitory neurons / c-28 | 1 · testis / spermatocytes / c-1 | 8.51 · 7 day / primitive endoderm | 2.46× |
| MAGEB18 | 0 · Adipose_Subcutaneous | 2.71 · testis.279 | 0.4 · colon / undifferentiated cells / c-8 | 6.2 · testis / spermatocytes / c-16 | 0.486 · 5 day / epiblast | 2.71× |
| MOS | 0.5 · nucleus of the diagonal band | 2.06 · testis.279 | 0.1 · brain / excitatory neurons / c-20 | 10.8 · testis / spermatocytes / c-16 | 27.3 · 3 day / unspecified lineage | 7.42× |
| OCM2 | 0.39 · liver.367 | 1.51 · testis.354 | 0.1 · brain / excitatory neurons / c-0 | 1.5 · ovary / oocytes / c-15 | 4.88 · 3 day / unspecified lineage | 1.51× |
| OTUD6A | 2.8 · cerebellar cortex | 15.2 · testis.279 | 2.7 · eye / rod photoreceptor cells / c-0 | 110 · testis / spermatocytes / c-1 | 12.3 · 4 day / unspecified lineage | 5.42× |
| PGAM4 | 1.1 · retrosplenial cortex | 0.4 · prostate.412 | 4.2 · small intestine / undifferentiated cells / c-9 | 1.1 · fallopian tube / macrophages / c-14 | 22 · 6 day / epiblast | 5.16× |
| PKDREJ | 1.5 · occipital cortex (white matter) | 9.69 · testis.349 | 2.2 · bone marrow / t-cells / c-6 | 19.2 · testis / spermatocytes / c-16 | 0.0506 · 6 day / trophectoderm | 3.08× |
| PPP4R3C | 0 · Adipose_Subcutaneous | 5.45 · testis.349 | 0 · adipose tissue / adipocytes / c-12 | 22.7 · testis / spermatocytes / c-1 | 6.69 · 4 day / unspecified lineage | 5.45× |
| PRDM9 | 2.4 · angular gyrus (white matter) | 10.6 · testis.349 | 0.2 · brain / inhibitory neurons / c-25 | 61.8 · testis / spermatocytes / c-1 | 4.03 · 5 day / primitive endoderm | 4.4× |
| PROKR1 | 2.78 · adipose tissue.137 | 4.21 · epididymis.407 | 1.2 · stomach / t-cells / c-4 | 4.4 · breast / adipocytes / c-7 | 0.953 · 5 day / epiblast | 1.52× |
| PRR20A/B/C/D/E | 0.0172 · thymus.1599 | 2.44 · testis.282 | 0 · adipose tissue / adipocytes / c-12 | 0 · breast / adipocytes / c-23 | 37.9 · 3 day / unspecified lineage | 7.28× |
| PRR23C | 1.3 · flocculonodular lobe | 2.65 · testis.282 | 0.8 · adipose tissue / dendritic cells / c-20 | 5.9 · testis / late spermatids / c-5 | 8.59 · 4 day / unspecified lineage | 4.84× |
| PRSS38 | 0.3 · nuclei of the trapezoid body | 9.44 · Testis | 2.8 · adipose tissue / t-cells / c-7 | 51.1 · testis / spermatocytes / c-1 | 0.647 · 3 day / unspecified lineage | 7.84× |
| PSMA8 | 3.19 · lymph node.173 | 55.4 · testis.349 | 8.5 · liver / plasma cells / c-3 | 244 · testis / spermatocytes / c-17 | 6.89 · 5 day / trophectoderm | 9.96× |
| RBMY1B | 0.453 · kidney.377 | 3.59 · testis.280 | 0 · adipose tissue / adipocytes / c-12 | 0 · breast / adipocytes / c-23 | 11.8 · 3 day / unspecified lineage | 8.79× |
| RBMY1E | 0.228 · kidney.58 | 8.92 · testis.244 | 0 · adipose tissue / adipocytes / c-12 | 2.2 · testis / spermatocytes / c-1 | 12.1 · 3 day / unspecified lineage | 8.92× |
| RPS4Y2 | 0.411 · spleen.238 | 9.96 · prostate.388 | 14.5 · liver / b-cells / c-16 | 4.3 · prostate / prostatic glandular cells / c-14 | 28.7 · 5 day / primitive endoderm | 0.696× |
| SLC35G3 | 0.376 · lymph node.174 | 6.35 · Testis | 4.3 · heart muscle / cardiomyocytes / c-1 | 56 · testis / early spermatids / c-8 | 21.6 · 5 day / primitive endoderm | 3.02× |
| SLC35G5 | 0.812 · Spleen | 5.26 · Testis | 8.7 · bone marrow / b-cells / c-9 | 17.5 · testis / early spermatids / c-8 | 4.34 · 5 day / unspecified lineage | 0.774× |
| SLC35G6 | 0.7 · paracentral lobule, posterior | 13.6 · Testis | 2 · lymph node / t-cells / c-5 | 42.4 · testis / late spermatids / c-13 | 0.968 · 6 day / primitive endoderm | 5.29× |
| SP7 | 1.6 · temporal white matter | 4.07 · testis.279 | 2 · lymph node / t-cells / c-5 | 1.9 · testis / late spermatids / c-13 | 18.8 · 7 day / primitive endoderm | 2.63× |
| SPANXD | 2.6 · preoptic area | 32.1 · Testis | 0.1 · salivary gland / serous glandular cells / c-10 | 2.43e+03 · testis / late spermatids / c-2 | 18.4 · 6 day / epiblast | 9.13× |
| SPANXN1 | 0 · Adipose_Subcutaneous | 6.01 · testis.244 | 2.2 · pancreas / exocrine glandular cells / c-10 | 37.6 · testis / late spermatids / c-13 | — | 6.01× |
| SPATA31A3 | 0.0295 · kidney.58 | 8.86 · testis.283 | 0 · adipose tissue / adipocytes / c-12 | 1.2 · testis / late spermatids / c-9 | 0.278 · 5 day / unspecified lineage | 8.86× |
| SPATA31A6 | 0.3 · hippocampus | 2.95 · testis.278 | 0.1 · brain / excitatory neurons / c-31 | 1 · testis / late spermatids / c-13 | 0.229 · 7 day / trophectoderm | 2.95× |
| SRMS | 4.57 · stomach.241 | 8.97 · prostate.207 | 2.9 · heart muscle / smooth muscle cells / c-3 | 0.2 · breast / t-cells / c-0 | 0.504 · 6 day / primitive endoderm | 1.96× |
| SSX5 | 3.5 · ventromedial nucleus | 1.84 · testis.244 | 0.6 · lung / endothelial cells / c-5 | 5.8 · testis / spermatogonia / c-4 | 17.4 · 4 day / unspecified lineage | 2.98× |
| SUPT20HL1 | 0.897 · Colon_Sigmoid | 2.01 · Fallopian_Tube | 3.7 · esophagus / fibroblasts / c-0 | 3.9 · testis / spermatocytes / c-1 | 1.27 · 3 day / unspecified lineage | 1.07× |
| TCP11X1 | 0.783 · lung.164 | 2.45 · testis.283 | 0.9 · brain / inhibitory neurons / c-23 | 0.1 · testis / spermatogonia / c-14 | 8.97 · Carnegie stage 7 / axial mesoderm | 8.97× |
| TEX13C | 0.3 · putamen | 4.21 · testis.282 | 0.1 · adipose tissue / fibroblasts / c-8 | 13.6 · testis / spermatocytes / c-1 | 4.26 · 7 day / primitive endoderm | 4.21× |
| TEX13D | 0.218 · smooth muscle.330 | 8.4 · testis.349 | 0.6 · pbmc / t-cells / c-2 | 26.3 · testis / spermatocytes / c-1 | 25.1 · 4 day / unspecified lineage | 9.83× |
| TGIF2LY | 0.105 · small intestine.233 | 9.33 · testis.244 | 0 · adipose tissue / adipocytes / c-12 | 36.8 · testis / early spermatids / c-8 | 2.17 · 4 day / unspecified lineage | 9.33× |
| TMEM95 | 1.2 · corpus callosum, splenium (white matter) | 25.7 · Testis | 35.9 · esophagus / schwann cells / c-14 | 86.4 · testis / early spermatids / c-8 | 0.708 · 3 day / unspecified lineage | 9.11× |
| TRIM51G | 0.1 · perirhinal gyrus | 0.121 · testis.278 | 0 · adipose tissue / adipocytes / c-12 | 0.1 · testis / spermatogonia / c-14 | 7.97 · 3 day / unspecified lineage | 2.72× |
| TRIM64C | 0.0342 · rectum.333 | 0.128 · testis.278 | 0 · adipose tissue / adipocytes / c-12 | 0 · breast / adipocytes / c-23 | 10.4 · 3 day / unspecified lineage | 5.37× |
| TTLL8 | 1.5 · lymph node.175 | 16 · Testis | 0.8 · bone marrow / b-cells / c-3 | 9 · testis / late spermatids / c-11 | 0.125 · 6 day / trophectoderm | 8.68× |
| XKR3 | 0.331 · adipose tissue.331 | 25.4 · testis.371 | 1.3 · salivary gland / serous glandular cells / c-16 | 297 · testis / early spermatids / c-8 | — | 8.48× |
| XKR5 | 1.6 · dorsomedial nucleus | 1.55 · prostate.365 | 1.2 · rectum / mixed immune cells / c-8 | 7.2 · endometrium / endothelial cells / c-23 | 1.16 · 6 day / primitive endoderm | 0.971× |
| ZACN | 24.9 · skin.227 | 22 · prostate.206 | 4.8 · thymus / t-cells / c-18 | 1.8 · fallopian tube / ciliated cells / c-12 | 18.1 · 7 day / primitive endoderm | 0.883× |
| ZNF560 | 2.8 · retina | 15.5 · Testis | 3.9 · bronchus / basal respiratory cells / c-14 | 140 · testis / early spermatids / c-8 | 26.4 · 4 day / unspecified lineage | 7.23× |
| ZNF705D | 0.3 · choroid plexus | 3.88 · epididymis.407 | 0.6 · kidney / proximal tubular cells / c-5 | 6.7 · testis / late spermatids / c-13 | 9.79 · 3 day / unspecified lineage | 9.49× |

## Review of the 31 holds

Eye and retinal findings remain Vital-tissue concerns. All 31 remain excluded from antigen sets. None is rescued by low bulk RNA or a mapping-risk flag. Established gene-specific normal-cell expression would move a hold to Unsafe; the current review does not establish that attribution.

### What the deeper review found

- **Ocular evidence:** HPA's retina values under its GTEx label are actually [EyeGEx RNA processed by HPA](https://v23.proteinatlas.org/about/assays%2Bannotation). The older rare-cell confirmation summary omitted eye; the table below includes retinal bulk, including RAD21L1 (1.2 TPM), RGPD3 (1.3) and RGPD4 (1.4). These are additional concerns, not proof that the same cell or exact locus produced both signals. There is no separate HPA retinal bulk measurement here. Retinal RNA cannot assess corneal, conjunctival or optic-nerve populations, and duplicated HPA exports of EyeGEx are one source.
- **Replication:** collection IDs are not guaranteed independent studies. The [integrated retinal atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC13050529/) combines public studies and new samples. Its final publication postdates the pinned Census snapshot, so original-study and donor mapping must be checked against that snapshot. Counts below are collections and source-specific donor labels, not verified unique patients. Possible overlap does not dismiss the signals.
- **PSG8:** the strongest ocular observation reaches 46% multi-molecule cells; the paralog attribution question remains consequential. [Kurji et al.](https://doi.org/10.1167/iovs.09-3622) also report a PSG8-labelled microarray signal in cultured RPE (Table 2). The probe has not been checked for PSG-family specificity and the reported validation panel did not validate PSG8. This is additional RNA evidence, not exact-protein confirmation.
- **Host-overlap hypotheses:** CAGE1, PPP3R2, RAD21L1 and the other H-coded genes need locus-specific attribution. Same-strand overlap does not demonstrate host pre-mRNA contamination. The previous bulk confirmation score was the weaker source maximum, not the highest normal expression: for example, brain HPA reaches 0.5 TPM for CAGE1 and 2.7 for PPP3R2, despite much lower GTEx values. Negative bulk cannot rule out a rare essential population.
- **Single-collection findings:** DCAF8L2 (16% retinal bipolar cells), FBXO47 (6.8% keratocytes) and NKX2-6 (37% atrial cardiomyocytes) lack a mapping-risk flag but still need independent confirmation. Their frequency warrants holding them; it does not justify calling them artifacts.
- **DAZL and CTAG2:** [Zhang et al.](https://doi.org/10.1186/s12885-020-07155-y) provide DAZL glioblastoma-expression support, but their normal astrocyte comparator does not resolve the ocular astrocyte signal. CTAG2 has low but nonzero lung bulk measurements and remains below reproductive 10× even if its review hold were resolved. Both newly reviewed papers are recorded centrally with exact gene IDs and permitted uses.

No selected organ has reliable positive HPA staining in this extraction. Usually the assay is missing; only measured negative rows count as negative staining. Neither missing staining nor a negative bulk measurement clears a cell-specific concern. Use Unsafe|Other for conservative peptide subtraction so unresolved proteins remain in the exclusion background.

| Protein | Organ / strongest population | Collections / donor labels | Peak ≥2-UMI cells | Same-organ HPA / GTEx–EyeGEx TPM | Same-organ cell nTPM | Codes |
|---|---|---:|---:|---:|---:|---|
| ARMH2 | lung / multiciliated columnar cell of tracheobronchial tree | 2 / 4 | 30.8% | 6.06 / 0.4 | 4.4 | H F |
| CAGE1 | brain / dopaminergic neuron | 2 / 3 | 3.29% | 0.5 / 0.0125 | 2.6 | H |
| CTAG2 | lung / pulmonary alveolar type 2 cell | 2 / 3 | 4.29% | 0.537 / 0.3 | 0 | P F |
| DAZL | eye / astrocyte | 1 / 4 | 5.86% | — / 0 | 0 | H S E |
| DCAF8L2 | eye / invaginating midget bipolar cell | 1 / 2 | 16.3% | — / 0 | 0 | S E C |
| EFCAB3 | lung / pulmonary alveolar type 2 cell | 2 / 10 | 15.9% | 1.95 / 0 | 0.1 | P F |
| FBXO47 | eye / keratocyte | 1 / 1 | 6.76% | — / 0 | 0 | S E |
| GAGE10 | brain / L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 17 | 3.53% | 1.9 / 0.6 | 1.2 | P F |
| NKX2-6 | heart / regular atrial cardiac myocyte | 1 / 1 | 36.8% | 0.197 / 0.786 | 0.2 | S F |
| OOSP4B | brain / L2/3-6 intratelencephalic projecting glutamatergic neuron | 2 / 4 | 3.11% | 0.0223 / 0 | 0 | H |
| OR4N2 | brain / GABAergic neuron | 5 / 15 | 2.57% | 3.06 / 0.0118 | 2.3 | P F |
| POTEC | eye / corneal endothelial cell | 2 / 2 | 12.1% | — / 0 | 0 | P E |
| POTEF | heart / ventricular cardiac muscle cell | 1 / 1 | 8.35% | 0.0415 / 0.0195 | 0 | P S F |
| PPP3R2 | brain / dopaminergic neuron | 3 / 49 | 8.84% | 2.7 / 0.0253 | 0 | H |
| PRDM7 | eye / melanocyte | 4 / 6 | 12.9% | — / 0 | 0 | P F E |
| PSG3 | brain / L2/3-6 intratelencephalic projecting glutamatergic neuron | 4 / 21 | 3.48% | 0.1 / 0.0133 | 1.1 | P |
| PSG8 | eye / rod bipolar cell | 8 / 103 | 46% | — / 0 | 0 | P E C |
| RAD21L1 | eye / ON midget ganglion cell | 4 / 35 | 6.04% | — / 1.2 | 0.3 | H E C |
| RGPD3 | eye / diffuse bipolar 2 cell | 3 / 17 | 15.1% | — / 1.3 | 1 | P F E C |
| RGPD4 | eye / diffuse bipolar 2 cell | 3 / 16 | 11.8% | — / 1.4 | 3 | P F E C |
| RGSL1 | brain / lamp5 GABAergic cortical interneuron | 4 / 26 | 5.09% | 0.313 / 0.0124 | 0.9 | H |
| RIPPLY3 | lung / pulmonary alveolar type 2 cell | 2 / 4 | 2.51% | 0.72 / 0.624 | 2.5 | P F |
| SATL1 | brain / L6b glutamatergic cortical neuron | 13 / 150 | 23.4% | 1.3 / 0.42 | 1.5 | P F |
| TAS2R1 | brain / direct pathway medium spiny neuron | 3 / 7 | 2.7% | 3.2 / 0.373 | 0.1 | H |
| TUBB8 | brain / glutamatergic neuron | 6 / 34 | 9.05% | 0.4 / 0.3 | 1.8 | P |
| VCX3B | eye / retinal ganglion cell | 2 / 3 | 4.44% | — / 0.5 | 0.8 | P E C |
| ZC3H11B | brain / glutamatergic neuron | 2 / 5 | 6.23% | 0.107 / 0.2 | 0 | P F |
| ZNF280A | eye / conjunctival epithelial cell | 2 / 6 | 6.85% | — / 0 | 0 | P E |
| ZNF492 | brain / L2/3-6 intratelencephalic projecting glutamatergic neuron | 3 / 4 | 1.62% | 2.5 / 0.13 | 2.8 | P F |
| ZNF705G | eye / Mueller cell | 2 / 4 | 1.47% | — / 0.2 | 2.6 | P F E C |
| ZNF99 | kidney / kidney loop of Henle thick ascending limb epithelial cell | 4 / 8 | 4.64% | 0.279 / 0.109 | 0 | P F |

[31-gene review details](review-31.csv) · [Every organ finding](all-organ-findings.csv) · [Donor observations](donor-evidence.csv.gz) · [Primary source registry](census-sources.csv) · [HPA assay measurements](hpa-assay-measurements.csv.gz) · [Protein staining](protein-staining.csv.gz) · [All 124 comparisons](expression-comparison.csv) · [Receipt](manifest.json)
