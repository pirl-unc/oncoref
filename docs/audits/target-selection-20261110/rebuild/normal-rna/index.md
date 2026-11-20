# Normal RNA sources and antigen membership

Every native-TPM reference is evaluated within its anatomical context. The highest supported within-site statistic across sources sets N. Disposable tissue supplies a separate D. Both antigen routes use **H₀ = max(N, D, 1 TPM)**. The source panels overlap; measurements are never added across studies.

[Cell and isoform evidence](../../../../curation.md#cell-and-isoform-evidence) describes expression in olfactory neurons, developing teeth, bone, and peri-implantation embryos, with the assay and age for each observation.

![Tissue coverage and role by native TPM source](tpm-coverage.png)

A colored cell indicates measured coverage for at least one assessed group. Reproductive context denotes a tissue summary; qualifying RA evidence requires an individual reproductive sample or supported embryo-lineage mean. Gray denotes no qualifying coverage in that source and role. The inventory lists every source-native site and cell type used in these comparison records.

[Download source, tissue, native context, cell type and role](source-tissues.csv)

| Source | Measurement and use |
|---|---|
| GTEx P95 | GTEx v11. Linear 95th percentile within each anatomical site across donors. |
| GTEx mean+2SD | GTEx v11. Sample mean plus two sample standard deviations within each anatomical site across donors. |
| GTEx samples | GTEx v11. Individual reproductive samples can witness RA. |
| GTEx median | GTEx v11. Published site-level donor medians, with exact member-locus mapping. |
| HPA biopsies | HPA v23. Within-biopsy transcript TPM summed by locus. Within-tissue P95 sets the healthy biopsy reference and each reproductive biopsy can witness RA. Donor IDs are unavailable. |
| HPA tissue | HPA v23 tissue means in TPM. |
| GTEx via HPA | GTEx tissue summaries distributed in HPA v23. These overlap the GTEx resource. |
| HPA brain | HPA brain subregion means in TPM, mapped to anatomical regions. |
| GTEx brain via HPA | GTEx brain region summaries distributed by HPA v23. |
| HPA sorted immune | P95 of donor-level native TPM within each sorted blood population. |
| Monaco | Mean TPM in sorted blood immune populations. The immune screen also evaluates these measurements. |
| DICE | Supported blood immune-cell TPM, using within-population P95 for single loci and the conservative bound for unpaired protein-group members. The immune screen also uses DICE. |
| Petropoulos | Preimplantation embryo-lineage means in TPM; stage, cell and independent embryo support apply. |
| Tyser | Carnegie stage 7 embryo-lineage means in TPM; cell and embryo support apply. |
| HPA consensus | Normal tissue consensus nTPM defines anatomical coverage and contributes to the shared expression screen. |
| HPA cell clusters | Tissue cell-type cluster nTPM contributes to shared normal-expression and reproductive-context evidence. |
| HPA immune | Blood immune-cell nTPM evaluated against the tissue-tier limits. |
| HPA brain nTPM | Brain subregion nTPM contributes to the independent normal-expression screen. |
| Collin cornea | Adult cornea/conjunctiva cell populations evaluated in counts per 10,000. |
| Census donors | CELLxGENE Census donor-resolved single-cell and single-nucleus counts. The rare-cell gate uses within-donor frequency and study replication. |

## Tissues by source

Tissue names follow each panel's anatomical mapping. The downloadable inventory preserves finer site, subregion and cell-type labels.

| Source | Tissues represented |
|---|---|
| GTEx P95 | adipose tissue, adrenal gland, amygdala, basal ganglia, blood, blood vessel, cerebellum, cerebral cortex, colon, esophagus, heart muscle, hippocampal formation, hypothalamus, kidney, liver, lung, midbrain, pancreas, peripheral nerve, pituitary gland, salivary gland, skeletal muscle, skin, small intestine, spinal cord, spleen, stomach, thyroid gland, urinary bladder |
| GTEx mean+2SD |  |
| GTEx samples | breast, cervix, fallopian tube, ovary, prostate, testis, uterus, vagina |
| GTEx median | adipose tissue, adrenal gland, amygdala, basal ganglia, blood, blood vessel, breast, cerebellum, cerebral cortex, cervix, colon, esophagus, fallopian tube, heart muscle, hippocampal formation, hypothalamus, kidney, liver, lung, midbrain, ovary, pancreas, peripheral nerve, pituitary gland, prostate, salivary gland, skeletal muscle, skin, small intestine, spinal cord, spleen, stomach, testis, thyroid gland, urinary bladder, uterus, vagina |
| HPA biopsies | adipose tissue, adrenal gland, appendix, bone marrow, breast, cerebral cortex, cervix, choroid plexus, colon, duodenum, endometrium, epididymis, esophagus, fallopian tube, gallbladder, heart muscle, kidney, liver, lung, lymph node, ovary, pancreas, parathyroid gland, placenta, prostate, rectum, salivary gland, seminal vesicle, skeletal muscle, skin, small intestine, smooth muscle, spleen, stomach, testis, thymus, thyroid gland, tongue, tonsil, urinary bladder |
| HPA tissue | adipose tissue, adrenal gland, appendix, bone marrow, breast, cerebral cortex, cervix, choroid plexus, colon, duodenum, endometrium, epididymis, esophagus, fallopian tube, gallbladder, heart muscle, kidney, liver, lung, lymph node, ovary, pancreas, parathyroid gland, placenta, prostate, rectum, salivary gland, seminal vesicle, skeletal muscle, skin, small intestine, smooth muscle, spleen, stomach, testis, thymus, thyroid gland, tongue, tonsil, urinary bladder |
| GTEx via HPA | adipose tissue, adrenal gland, amygdala, basal ganglia, breast, caudate, cerebellum, cerebral cortex, cervix, colon, endometrium, esophagus, fallopian tube, heart muscle, hippocampal formation, hypothalamus, kidney, liver, lung, ovary, pancreas, pituitary gland, prostate, retina, salivary gland, skeletal muscle, skin, small intestine, spinal cord, spleen, stomach, substantia nigra, testis, thyroid gland, urinary bladder, vagina |
| HPA brain | amygdala, basal ganglia, cerebellum, cerebral cortex, choroid plexus, hippocampal formation, hypothalamus, medulla oblongata, midbrain, pons, spinal cord, thalamus, white matter |
| GTEx brain via HPA | amygdala, basal ganglia, cerebellum, cerebral cortex, hippocampal formation, hypothalamus, midbrain, pituitary gland, retina, spinal cord |
| HPA sorted immune | blood |
| Monaco | blood |
| DICE | blood |
| Petropoulos | embryo |
| Tyser | embryo |
| HPA consensus | adipose tissue, adrenal gland, amygdala, appendix, basal ganglia, bone marrow, breast, cerebellum, cerebral cortex, cervix, choroid plexus, colon, duodenum, endometrium, epididymis, esophagus, fallopian tube, gallbladder, heart muscle, hippocampal formation, hypothalamus, kidney, liver, lung, lymph node, midbrain, ovary, pancreas, parathyroid gland, pituitary gland, placenta, prostate, rectum, retina, salivary gland, seminal vesicle, skeletal muscle, skin, small intestine, smooth muscle, spinal cord, spleen, stomach, testis, thymus, thyroid gland, tongue, tonsil, urinary bladder, vagina |
| HPA cell clusters | adipose tissue, bone marrow, brain, breast, bronchus, colon, endometrium, esophagus, eye, fallopian tube, heart muscle, kidney, liver, lung, lymph node, ovary, pancreas, pbmc, placenta, prostate, rectum, salivary gland, skeletal muscle, skin, small intestine, spleen, stomach, testis, thymus, tongue, vascular |
| HPA immune | blood |
| HPA brain nTPM | amygdala, basal ganglia, cerebellum, cerebral cortex, choroid plexus, hippocampal formation, hypothalamus, medulla oblongata, midbrain, pons, spinal cord, thalamus, white matter |
| Collin cornea | cornea/conjunctiva |
| Census donors | adipose tissue, adrenal gland, blood, brain, colon, esophagus, eye, heart, kidney, liver, lung, lymph node, pancreas, small intestine, spinal cord, stomach, tongue, ureter, urinary bladder |

## Cross-tissue RNA in individual libraries

P denotes a pancreatic enzyme signature (PRSS1, CLPS, PNLIP, CELA3A). E denotes a squamous signature observed in esophageal mucosa (KRT4, KRT13). These markers come from [Nieuwenhuis et al.](https://www.nature.com/articles/s41467-020-15821-9).

OncoRef compares each library with samples from the same anatomical site or sorted immune population. The panel requires at least three complete libraries. Every marker must exceed its within-site interpolated P95 and exceed the site median by at least 0.1% of its native GTEx median. The source organ is skipped. A blank flag means the panel could not be assessed. The thresholds here are OncoRef choices.

These are advisory sample flags. Normal epithelial content can produce E, particularly in tissues with squamous epithelium. Shared marker expression alone cannot identify a transfer event. Tissue means, unpaired DICE ranks and single-cell UMI counts cannot support this within-library bulk TPM comparison.

| Source | Libraries | P flagged / assessed | E flagged / assessed | Either |
|---|---:|---:|---:|---:|
| GTEx v11 tissue | 18,637 | 209 / 18,275 | 181 / 17,059 | 381 |
| HPA v23 tissue | 186 | 2 / 176 | 1 / 173 | 3 |
| HPA v23 sorted immune | 109 | 0 / 109 | 0 / 109 | 0 |

[GTEx libraries](../gtex-library-qc.csv.gz) · [HPA tissue biopsies](../normal-library-qc.csv) · [HPA sorted immune libraries](immune-library-qc.csv) · [GTEx batch matches and percentile sensitivity](../../../gtex-cross-tissue-signals/index.md)

## Independent normal RNA screens

![Independent RNA screen coverage](shared-coverage.png)

The tissue-tier, immune, corneal and rare-cell screens can prevent shared eligibility before enrichment is evaluated. Some native TPM sources also inform these screens. Cluster nTPM, CP10K, donor UMIs and protein staining retain their own evidence rules. A bulk percentile cannot override an exclusion.

## From measurements to memberships

For GTEx v11 normal sites, sum the measured loci encoding each identical protein within a library. Give each donor one value per site (median of replicate libraries), then take the linearly interpolated 95th percentile of the measured donors. Interpolation produces a tissue statistic without a sample or donor ID. HPA biopsies contribute the linear P95 of measured biopsy TPM per healthy tissue. Missing required tissue coverage blocks qualification.

![Sources setting the reference for each assessed protein](denominator-sources.png)

| Route after shared eligibility | Definition | Protein groups |
|---|---|---:|
| CA | At least one QC-passing tumor reaches 10 × H₀ | 397 |
| RA | At least one reproductive sample or supported embryo-lineage observation reaches 10 × H₀ | 455 |
| CTA | Both CA and RA | 327 |
| CA ∪ RA | Either route | 525 |

![Cancer and reproductive antigen overlap](../summary/venn.png)

[Set overlap and eligibility funnel](../summary.md) · [Exact comparison contract](../comparison.json) · [Census study inventory](../../../normal-cell-atlas-20261110/index.md)
