# Other proteins: what the proposed subcategories actually contain

This is a diagnostic audit of the pending 1.25.0 data, **not a released taxonomy change**. The package still calls this top-level category Neither. The proposed user-facing name is **Other**. Counts refer to identical-reference-sequence proteins, not gene loci. Missing required normal evidence remains Unassigned; Other requires a normal-screen pass.

[All 725 current proteins and reasons](all-other-proteins.csv) · [All 400 current generic-bucket proteins](current-generic-bucket.csv) · [66 tumor-supported proteins](tumor-supported.csv) · [25 bulk corrections](bulk-correction-CTA.csv) · [107 further cell-evidence cases](cell-evidence-reconsideration.csv)

## The three names in question

| Proposed description | Count under current rules | Meaning |
|---|---:|---|
| Other:Reproductive | 1 | CDY1/B: qualifying reproductive IHC but no qualifying recorded tumor RNA. This is not the count of all proteins with any biological reproductive expression. |
| Other:Tumor-supported | 66 | Two distinct tumors ≥10 native TPM, low measured normals and no qualifying reproductive evidence; rare-normal-cell/state restriction is not yet established. |
| Existing generic bucket | 400 | This is a legacy fallback, not a biological class. It includes CDY1/B and 25 proteins affected by the stale Reproductive input. It must not simply become Other:Other. |

These three rows are an audit of existing flags, **not a disjoint partition**. The current 725 rows comprise 400 generic, 321 low-normal, three legacy Manageable-expression profiles and one diffuse profile. The 66 tumor-supported leads are contained in the 321 low-normal rows. The three legacy Manageable profiles are OXGR1, KRT9 and TMDD1; they must not be relabelled Disposable.

## Reproductive-only: CDY1/B

The identical-reference-sequence protein is supported by High/Supported testicular IHC at both loci. Its full HPA v23 reproductive bulk maximum is **3.1 nTPM**, compared with the stale input 2.7. OncoRef tumor RNA reaches **1.597 TPM** and the separate HPA tumor lower bound reaches **1.208 pTPM**, below the current one-tumor ≥2 rule. It therefore remains Other under the present cancer-admission rule.

[HPA v23 CDY1B pathology evidence](https://v23.proteinatlas.org/ENSG00000172352-CDY1B/pathology) also describes liver-cancer staining. That is a tumor-protein follow-up, not an automatic RNA-cutoff override: antibody family attribution and the staining localization need review. Reliable tumor-protein admission is currently implemented only in the Tumor-restricted prototype, another asymmetry to resolve before a broad admission change.

## Tumor-supported does not mean reproductive or tumor-exclusive

All 66 fail current Reproductive admission. Their full adult/placental bulk maximum is at most **0.9 nTPM** and available developmental maxima are below **2 TPM**. None has qualifying reproductive IHC or a ≥5 nTPM/≥20-read reproductive cell observation in this audit. Missing optional measurements remain missing, and early developmental coverage is not exhaustive fetal-organ coverage.

Forty labels are olfactory receptors. Other examples are **GFRAL, SCN10A, CNGA2, IFNA4, IFNW1, OTOR and OC90**. The relevant normal contexts include specialized sensory cells, brainstem nuclei and induced immune states. These are tumor-supported leads, not established tumor-exclusive targets. See the [full normal-coverage follow-up](https://github.com/pirl-unc/oncoref/issues/653).

## The generic bucket exposes two admission gaps

**1. Stale bulk tissue definition.** Reproductive admission inherited an older, narrower maximum although the current full HPA profile already includes epididymis and the other mapped branches. MAGEA11 was assessed at 4.3 nTPM even though epididymis is **20.5**. Correcting only this input, while retaining every normal hold, tumor rule and concentration gate, gives **598 CTA / 700 Other / 18,031 Unsafe / 534 Unassigned**. Every existing CTA remains CTA and Unsafe membership is identical.

| Protein | Old input | Full Reproductive nTPM | Tissue |
| --- | --- | --- | --- |
| DEFB130A/B | 0.0 | 5.0 | epididymis |
| DEFB127 | 0.0 | 1050.8 | epididymis |
| RLN1 | 0.3 | 22.7 | prostate |
| CCN6 | 0.5 | 9.9 | epididymis |
| DEFB126 | 0.1 | 1251.1 | epididymis |
| MMP26 | 0.1 | 8.8 | endometrium |
| WFDC13 | 0.1 | 87.3 | epididymis |
| DEFB114 | 0.0 | 51.4 | epididymis |
| WFDC9 | 0.0 | 964.7 | epididymis |
| WFDC10A | 0.1 | 280.5 | epididymis |
| DEFB112 | 0.0 | 24.0 | epididymis |
| DEFB108B | 0.0 | 12.5 | epididymis |
| MAGEA11 | 4.3 | 20.5 | epididymis |
| DEFB128 | 0.0 | 1163.1 | epididymis |
| RNASE9 | 0.0 | 96.2 | epididymis |
| SCGB1D4 | 1.0 | 52.2 | fallopian tube |
| TEDDM1 | 0.4 | 562.1 | epididymis |
| DEFB110 | 0.0 | 848.5 | epididymis |
| DEFB135 | 0.1 | 55.3 | epididymis |
| RNASE13 | 0.5 | 299.6 | epididymis |
| DEFB113 | 0.0 | 46.1 | epididymis |
| GPX5 | 1.4 | 9.0 | epididymis |
| PATE3 | 0.2 | 154.4 | epididymis |
| SERPINA2 | 0.0 | 104.5 | epididymis |
| SPADH | 0.0 | 115.4 | epididymis |

**2. Supported reproductive-cell RNA is not an admission route.** The current rule uses bulk RNA, supported developmental RNA or reliable reproductive IHC. The existing source cache nevertheless contains **113 current Other proteins** with mapped reproductive-cell RNA ≥5 nTPM and the existing ≥20-read support requirement. Six overlap the 25 bulk corrections, leaving **107 further cases**. This diagnostic evidence flag does not automatically promote them: it needs a consistent positive-association route, cell identity/attribution review and an explicit treatment of concentration when bulk is diluted. Endothelial or other generic cells found inside a reproductive organ must not automatically establish reproductive-specific identity.

| Protein | Reproductive bulk nTPM | Reproductive cell nTPM (≥5; ≥20 reads) | Tumor TPM | Non-permitted normal nTPM | Non-permitted normal TPM |
| --- | --- | --- | --- | --- | --- |
| FOXR2 | 1.2 | 5.2 | 239.421 | 0.0 | 0.0 |
| CAGE1 | 4.9 | 208.9 | 17.402 | 3.6 | 1.056 |
| PRDM9 | 3.9 | 61.8 | 70.59 | 2.6 | 0.5 |
| KRT82 | 0.1 | unmeasured | 54.537 | 12.3 | 0.121 |
| KRT28 | 0.7 | unmeasured | 1887.549 | 11.3 | 0.11 |
| IL1F10 | 0.5 | unmeasured | 37.48 | 12.1 | 13.939 |
| IFNB1 | 0.0 | unmeasured | 55.667 | 6.0 | 1.594 |
| OR2T7 | 0.1 | 22.8 | 1.184 | 3.2 | 0.761 |

After removing the 25 bulk-input errors and separating CDY1/B, **374 of the original 400 generic rows remain**. The three legacy Manageable profiles sit outside this count. The residual includes weak bulk reproductive signals, overlooked cell evidence and expression in other normal tissues. 373 of the 374 pass the current cancer gate; **OR2T7** fails both cancer and Reproductive admission. Full row-level explanations and native units are in the linked tables.

## Naming recommendation

Use **CTA / Unsafe / Other**, with Unassigned retained as an assessment status for incomplete or unresolved normal evidence. Within Other:

- **Low** means low measured expression on a stated panel, not never expressed. The existing 321 flag is low *normal* expression and says nothing about cancer abundance.
- **Disposable** requires affirmative Disposable expression and the agreed below-threshold conditions elsewhere. Skin/Manageable expression is not automatically Disposable.
- **Diffuse** means a measured dispersed profile, not a synonym for missing reproductive evidence.
- **Reproductive** means the Reproductive route qualifies but cancer evidence does not.
- **Tumor-supported** means strong cancer evidence without qualifying Reproductive identity or established tumor restriction.

Do not invent an Other:Other biological subtype. Leave Other without a suffix when the evidence does not establish one of these profiles, and show the specific failed or missing admission evidence. Display secondary evidence flags as well as a disjoint primary subtype; a low-normal profile and tumor support can coexist. Final subtype counts require the two admission gaps and the Disposable definition to be resolved first.

## Sources and reproduction

[Human Protein Atlas v23 downloadable datasets](https://v23.proteinatlas.org/about/download) defines separate consensus bulk and tissue-cluster RNA tables. This audit uses the same pinned cache, accepted gene identities, anatomy mappings and supported cell rows as the existing screen; it does not pool bulk and single-cell values into one numerical assay. The ≥5/≥20 diagnostic flag reuses project conventions, not an HPA clinical-safety cutoff.

Run `python scripts/audit_cta_other_subtypes.py`. The bulk counterfactual reruns association, concentration and classification while preserving every original normal-safety mask. Source hashes and output hashes are recorded in [manifest.json](manifest.json). No packaged selection table, report PDF or released artifact is rewritten by this script.
