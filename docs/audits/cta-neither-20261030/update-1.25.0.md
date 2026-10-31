# OncoRef 1.25: categorical CTA, Unsafe and Other decisions

**698 CTA, 600 Other, 18,031 Unsafe and 534 unassessed proteins.** All 522 previous CTAs remain included; Unsafe membership is unchanged. This release adds 176 CTA proteins. Counts collapse identical selected reference protein sequences; they are not gene-locus counts.

| Category | 1.24.0 | 1.25.0 |
| --- | --- | --- |
| CTA | 522 | 698 |
| Other | 20 | 600 |
| Unsafe | 18031 | 18031 |
| Unassigned | 1290 | 534 |

[Full report and funnel](../cta-source-walkthrough-20261028/index.md) · [PDF logic slides](../cta-source-walkthrough-20261028/cta-selection-1.25.0.pdf) · [Every protein, category and measurement](protein-categories-1.25.0.csv) · [Every transition](transitions-1.25.0.csv) · [Added CTAs](restored-CTA-1.25.0.csv)

## One category per assessed protein

The categorical label is `classification`: CTA, Unsafe, or one of the named Other categories below. `category` retains the broad family and `subtype` the Other name. Every Other protein receives one subtype; there is no Other:Other or unnamed remainder. Insufficient normal evidence is explicitly Unassessed, not a normal-pass category.

| Category | Proteins | Definition | Examples |
| --- | --- | --- | --- |
| Other:Diffuse | 8 | Complete HPA profile is dispersed (score >10% with a reference ≥2 nTPM), or the active concentration gate fails; no CTA route remains. | H2BW2, HUS1B, KRTAP5-9, MTCP1, OR5AU1 |
| Other:Reproductive | 2 | Qualifying Reproductive evidence, but no qualifying cancer expression; the normal screen passes. | CDY1/B, OR2T7 |
| Other:Inducible | 5 | Cited expression in activated nonreproductive normal cells; a localized indirect atlas signal cannot establish Reproductive identity. | IFNA14, IFNA17, IFNA2, IFNA5, IFNA7 |
| Other:Tumor-supported | 66 | Two distinct tumors ≥10 native TPM or reviewed reliable tumor protein, with low measured normal/developmental expression; not established tumor restriction. | CLDN25, CNGA2, DBX1, GFRAL, HES3 |
| Other:Disposable | 2 | Positive Disposable RNA/protein; measured RNA elsewhere below the descriptive 2-unit floor and no detected protein elsewhere, including low measured developmental RNA. | PLA2G2E, TCHHL1 |
| Other:Low | 254 | All measured normal RNA below the descriptive 2-unit floor and no detected normal protein. Cancer expression need not be low. | AADACL4, AMELY, CLRN2, DEFB115, DEFB116 |
| Other:Normal | 171 | Supported nonreproductive normal RNA reaches the descriptive 2-unit floor or reliable protein is detected; no preceding subtype or CTA definition qualifies. | A3GALT2, AMELX, ANKRD60, ANP32D, APOL5 |
| Other:Weak | 92 | Reproductive signals are below admission requirements or lack support; no preceding subtype qualifies and no supported appreciable nonreproductive expression is established. | AADACL3, ARSH, ASCL4, BPIFB3, BPY2/B/C |

These categories are assigned in the table's order after safety and CTA admission. Secondary evidence flags remain available. For example, a protein can have both low normal expression and strong tumor evidence, but receives Other:Tumor-supported. Descriptive 2-unit floors mean native TPM/nTPM, with the pre-existing separate nominal CP10K floor; the assays are not pooled or converted into equivalent biological measurements. Low does not mean never expressed. Normal and Inducible do not mean demonstrated clinical injury. All Other proteins pass every recorded normal gate; supported non-permitted maxima remain below 30 TPM/nTPM.

Raw maxima are not always supported expression. **OR2B6 and SCGB1C2** each have a 47.6 nTPM platelet-cluster value supported by only 5 reads (the rule requires 20). **OR6C2** has a 30.92 TPM DICE maximum, but the source-specific third-highest value is 1.14 TPM. These remain **Other:Weak**, with the raw observations retained, rather than being described as never expressed. This support distinction applies consistently across proteins.

## The repaired admission rules

1. **Normal safety:** unchanged tissue/IHC RNA ceilings, protein rules and clinical exclusions. Required normal coverage and evidence resolution remain mandatory. No reproductive or cancer signal overrides a normal exclusion.
2. **Cancer presence:** one sample ≥2 native TPM or HPA pTPM, or reviewed reliable tumor protein. The protein route is now shared by CTA admission and tumor-support assessment. Its current review ledger has no accepted rows; aggregate cancer-IHC counts alone do not qualify.
3. **Reproductive presence:** full mapped HPA v23 reproductive bulk maximum ≥5 nTPM, supported reproductive-cell RNA ≥5 nTPM with ≥20 reads, reliable positive reproductive IHC, or supported developmental RNA ≥5 TPM. Bulk, cells and developmental values retain separate scales and provenance. All matched protein-member loci must be present before summing a normal cell context.
4. **Reproductive concentration:** the existing complete HPA bulk contrast remains active, including Disposable and excluding Permissive thymus. Cell evidence cannot override a measured concentration failure. A missing/zero bulk reference remains explicitly unassessed rather than being replaced by cell RNA.

The full bulk tissue definition fixes the MAGEA11 example: **20.5 nTPM in epididymis**, not the stale 4.3 input. The cell route recovers direct signals such as **CAGE1 208.9 nTPM in early spermatids**, **PRDM9 61.8 in spermatocytes**, and **FOXR2 5.2 in spermatogonia**. [Supporting cell rows and read counts](restored-reproductive-cell-1.25.0.csv) · [Supporting IHC rows](restored-reproductive-IHC-1.25.0.csv).

### Generic cells and macrophages

A generic, immune or mixed-cell signal inside reproductive tissue may provide **indirect sample-site support** when it reaches ≥5 nTPM with ≥20 reads, the same-source outside-reproductive cluster panel is complete, and every outside cluster is <5 nTPM. This includes all measured outside cell types, not only the matching cell label. All normal-source and concentration gates still apply. It establishes the sampling site, not which cell produced the transcript. [Single-cell RNA sequencing uncovers heterogenous transcriptional signatures in macrophages during efferocytosis (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7459098/) documents passenger transcripts in phagocytes; this is an interpretation limit, not proof of the mechanism for any gene here.

Published normal induction prevents treating an indirect-only atlas signal as reproductive restriction. **IFNA2 and IFNA14** have high breast-macrophage atlas signals, but [Plasmacytoid dendritic cells respond to Epstein-Barr virus infection with a distinct type I interferon subtype profile (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/) reports expression in activated human nonreproductive immune cells. They therefore become **Other:Inducible**, not CTA. The curated source table also records IFNA5, IFNA7 and IFNA17. This qualitative counterevidence is not fabricated TPM and does not itself supply Unsafe evidence. Unreviewed family members do not inherit it by name or sequence similarity.

## Named targets and boundary cases

| protein_group_label | classification | reproductive_ntpm | reproductive_cell_ntpm | reproductive_cell_indirect_ntpm | oncoref_maximum_tpm |
| --- | --- | --- | --- | --- | --- |
| CDY1/B | Other:Reproductive | 3.1 | unmeasured | unmeasured | 1.597 |
| CTAG2 | CTA | 10.1 | 84.0 | unmeasured | 604.332 |
| MAGEA4 | CTA | 21.5 | 203.6 | unmeasured | 3201.684 |
| PRDM9 | CTA | 3.9 | 61.8 | unmeasured | 70.59 |
| CAGE1 | CTA | 4.9 | 208.9 | 9.1 | 17.402 |
| MAGEA11 | CTA | 20.5 | 58.4 | unmeasured | 109.401 |
| PRAME | CTA | 102.2 | 318.0 | unmeasured | 3692.464 |
| IFNA2 | Other:Inducible | 0.0 | unmeasured | 9.8 | 17.203 |
| FOXR2 | CTA | 1.2 | 5.2 | unmeasured | 239.421 |
| MTCP1 | Other:Diffuse | 5.1 | 0.8 | unmeasured | 492.267 |
| USP17L3 | Other:Low | 0.0 | unmeasured | unmeasured | 3.82 |
| OR2T7 | Other:Reproductive | 0.1 | 22.8 | unmeasured | 1.184 |
| IFNA14 | Other:Inducible | 0.0 | unmeasured | 9.5 | 12.515 |
| XKR5 | CTA | 0.7 | 1.5 | 5.0 | 22.912 |
| FAM90A9 | Other:Low | 0.0 | unmeasured | unmeasured | 1.75 |
| NY-ESO-1 | CTA | 15.5 | unmeasured | unmeasured | 944.32 |
| XAGE1A/B | CTA | 98.8 | 9.5 | unmeasured | 4068.639 |

MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B remain CTA. CDY1/B has reliable reproductive IHC but fails the cancer-RNA floor; OR2T7 now has reproductive-cell support but also lacks qualifying cancer RNA. HPA tumor staining for CDY1/B remains a protein-evidence follow-up because localization and family attribution need adjudication; the ledger does not silently accept it. [HPA v23 CDY1B pathology](https://v23.proteinatlas.org/ENSG00000172352-CDY1B/pathology).

## Tumor-supported and normal-coverage caveats

66 Other:Tumor-supported proteins have two qualifying ≥10 TPM tumors. They lack qualifying Reproductive identity and are not declared tumor-exclusive. Rare sensory cells and activated immune states remain relevant coverage gaps. [Every qualifying tumor observation](tumor-restricted-lead-observations-1.25.0.csv) · [Coverage follow-up #653](https://github.com/pirl-unc/oncoref/issues/653). Missing optional matrices are not measured zeros.

The seven previous Neither proteins with normal holds remain unassessed: CRYGN, KRT34, KRT85, KRT71, LEP, LCE2A, KRT25. Their older narrow profile cannot override the normal gates. Historical audits and the separate exact-protein/k-mer policy remain frozen.

## API and reproducibility

```python
from oncoref import cta_selection

cta_selection.results(category="CTA")
cta_selection.results(category="Other", subtype="Inducible")
cta_selection.results()[["label", "classification", "reasons"]]
cta_selection.candidates()  # CTA only
cta_selection.audit().explain("FOXR2")
```

The former `category="Neither"` query is a deprecated alias for Other; returned labels use Other. The proposed but unreleased Tumor-restricted top-level category is removed. Dog and mouse keep the same schema with explicit Unassessed status; human admission is not transferred by orthology. Pin policy **oncoref.cta-selection.2026-10-30.2** and the source snapshot. The existing legacy CTA and exact-sequence exclusion APIs remain separate products, so downstream code must not silently change k-mer backgrounds.

Run `scripts/regenerate_cta_walkthrough.py --verify-replay`, build the PDF/PPTX slides, finalize the receipts, and run this script. Inputs and output checksums are in [update-manifest.json](update-manifest.json). The [pre-fix diagnostic audit](../cta-other-20261030/index.md) is historical; its counts are not this release.
