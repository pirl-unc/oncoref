# Other proteins: literature and category audit

This audit reviews the **600 Other proteins / 611 member genes in 1.25.0**. Every member was searched. Search retrieval and annotation triage are not equivalent to reading and adjudicating every paper. The per-protein table states the depth of review; absence of a search hit never proves normal absence.

[Full decisions](protein-decisions.csv) · [Changed proteins](changed-proteins.csv) · [Reviewed findings and sources](reviewed-findings.csv) · [PDF report](other-protein-review.pdf) · [Current complete CTA report](../cta-source-walkthrough-20261028/index.md)

## Current complete universe

| Classification | Proteins |
| --- | --- |
| Unsafe | 18062 |
| CTA | 701 |
| Unassessed | 646 |
| Other:Low | 210 |
| Other:Normal | 124 |
| Other:Weak | 67 |
| Other:Tumor-supported | 46 |
| Other:Diffuse | 6 |
| Other:Disposable | 1 |

## What changed among the 600

| baseline | classification | proteins |
| --- | --- | --- |
| Other:Diffuse | Other:Diffuse | 6 |
| Other:Diffuse | Unassessed | 2 |
| Other:Disposable | Other:Disposable | 1 |
| Other:Disposable | Unassessed | 1 |
| Other:Inducible | Unsafe | 5 |
| Other:Low | Other:Low | 210 |
| Other:Low | Unassessed | 39 |
| Other:Low | Unsafe | 5 |
| Other:Normal | Other:Normal | 124 |
| Other:Normal | Unassessed | 38 |
| Other:Normal | Unsafe | 9 |
| Other:Reproductive | CTA | 2 |
| Other:Tumor-supported | Other:Tumor-supported | 46 |
| Other:Tumor-supported | Unassessed | 14 |
| Other:Tumor-supported | Unsafe | 6 |
| Other:Weak | CTA | 1 |
| Other:Weak | Other:Weak | 67 |
| Other:Weak | Unassessed | 18 |
| Other:Weak | Unsafe | 6 |

## Rules and consequences

1. **Normal safety comes first.** Existing tissue/IHC RNA ceilings remain unchanged. Reviewed adult human normal protein, or localized RNA plus orthogonal normal-cell function, can establish exclusion in Vital/Protected. Rare cells and normal activation states count. Literature findings retain assay, species, stage, exact identity and limitations.
2. **Unresolved normal evidence stays out of Other.** Sparse sensory-cell RNA, unresolved species/isoform attribution, and high raw normal observations without sufficient support receive a review hold. This does not label them clinically toxic.
3. **Cancer presence means >0 measured TPM or HPA pTPM in one compatible sample, or reliable tumor protein.** Abundance is reported separately. No cancer/normal ratio is required. Microarray proxies, pseudobulk and invalid legacy scales cannot be compared numerically with TPM.
4. **Reproductive presence accepts RNA or reliable protein.** The RNA-only floor remains 5 nTPM, with the existing read/cell support rules. Reviewed reproductive protein is an independent route and must cover the selected identical-sequence members. Concentration, safety and coverage still apply.
5. **CTA / Unsafe / Other / Unassessed have distinct sequence roles.** CTA is a candidate direct target. Unsafe supplies an exclusion background. Other supplies no exclusion under this policy but is not a direct-target category. Unassessed requires review and cannot silently be treated as allowed. These are research roles, not epitope/HLA or all-isoform clearance.

## Not all Other:Normal proteins are skin proteins

Among the 170 baseline Normal proteins with supported nonreproductive RNA at or above the descriptive floor, **41 have skin measurements; 33 reach the floor only in skin**. Overlapping counts include **81 Vital**, **59 Protected**, **50 Manageable**, and **4 Disposable**. These counts overlap and must not be added. Most observations remain below the existing safety ceiling. The full table retains cell names and native scales.

[Every qualifying tissue/cell measurement](normal-category-tissues.csv) · [Tissue summary](normal-tissue-summary.csv) · [Raw and supported observations](normal-contexts.csv.gz)

## Reconsideration and limitations

- The human olfactory supplement identifies 68 baseline Other receptor genes in 1–20 positive neurons. Sparse detection is a review hold, not proof of protein abundance. Unlisted receptor genes are not thereby negative.
- Species, stage and isoforms matter: rat kidney IHC cannot become human IHC; fetal cochlear cloning does not prove adult protein; DUX4 full-length and shorter healthy-muscle isoforms need sequence-specific treatment.
- Three raw normal observations exceed 30 TPM/nTPM despite failing recurrence/read support. They are retained as review holds with the raw values visible, not relabeled zero.
- Published normal production is not a demonstrated therapy-related adverse event. Functional receptor inhibition also differs from cytotoxic targeting of the expressing cell.
- The baseline universe is representative selected protein sequences, not every splice isoform. Normal atlases still miss some rare cells, activation states and donor variability.

## Cancer measurements for CDY1/B and OR2T7

[All 143 registered cohort views, with native-scale and missingness labels](requested-cancer-sources.csv). Compatible source maxima remain approximately 1.597 TPM for CDY1/B and 1.184 TPM for OR2T7. HPA independently reports 1.208 pTPM for CDY1/B. Positive low values now attest to cancer presence; they are not claims of high tumor abundance.

## Reusable API and pinned sources

```python
from oncoref import cta_selection, gene_knowledge, reference_data

findings = gene_knowledge.findings("GFRAL")
roles = cta_selection.sequence_dispositions()
source = reference_data.ensure("olfactory_durante_supplement", "2020-02-18")
receipt = reference_data.provenance("olfactory_durante_supplement", "2020-02-18")
```

The knowledge reference ships in the wheel. The olfactory workbook uses the same versioned, checksummed reference-data cache as HPA/DICE/cornea. Native counts/fractions stay in their own units. Human findings never transfer to dog or mouse by orthology. Legacy exact-protein backgrounds remain separately versioned; downstream code must explicitly adopt these roles.

## Search accounting

1,350 search queries; 611 genes; 1,346 fully paginated queries. Four broad IFNA1/NPS queries reached the retrieval cap; aliases also create unrelated hits. The bibliography is a lead list, not accepted evidence. Accepted findings are listed separately with their limitations.

## Every baseline Other protein

| Protein | Before | Now | Reason |
| --- | --- | --- | --- |
| A3GALT2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. A normal human dendritic-cell glycosphingolipid product is not direct measurement of the encoded enzyme; review pathway and gene attribution. |
| AADACL3 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| AADACL4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| AMELX | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| AMELY | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| ANKRD60 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ANP32D | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| APOL5 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ARL14EPL | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ARL5C | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ARSH | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| ASCL4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| AWAT1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| BHLHE23 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| BPIFB3 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| BPIFB6 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| BPY2/B/C | Other:Weak | CTA | Independent cancer presence, reproductive presence, concentration, normal safety and required coverage pass. VCY2/BPY2 protein is detected in normal human male germ cells. |
| C1orf232 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| CCDC195 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| CDCP2 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| CDY1/B | Other:Reproductive | CTA | Independent cancer presence, reproductive presence, concentration, normal safety and required coverage pass. Positive cancer RNA below 2 TPM now establishes presence; it does not establish therapeutic abundance. |
| CHRNA9 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expressed in cochlea, keratinocytes, pituitary gland, B-cells and T-cells |
| CIB3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| CLDN25 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported enriched RNA in a defined adult human olfactory population. |
| CLRN2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Plays a key role to hearing function. Required for normal organization and maintenance of the stereocilia bundle and for mechano-electrical transduction Deafness, autosomal recessive, 117: A form of non-syndromic deafness characterized by prelingual, moderate-to-profound sensorineural hearing loss. Sensorineural hearing loss results from damage to the neural receptors of the inner ear, the nerve pathways to the brain, or the area of the brain that receives sound information. |
| CNGA2 | Other:Tumor-supported | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Normal human arterial endothelium and smooth muscle show localized CNGA2 RNA and protein. |
| CST9 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expressed in heart, placenta, lung, liver, skeletal muscle and pancreas (PubMed:12535658). Not expressed in brain (PubMed:12535658). Expressed in epididymis, kidney, testis, spinal cord, and thymus with a strong expression in epididymis and kidney and a weak expression in the spinal cord and thymus (PubMed:20565543) |
| CYP26C1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| DBX1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| DCANP1 | Other:Weak | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Human DCNP1/DCANP1 protein localizes in dendritic cells; normal-tissue immunostaining supports production. |
| DEFB115 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| DEFB116 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| DEFB125 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| DGAT2L6 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| DSG4 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. In scalp follicles, present in the inner root sheath (IRS) and all layers of the matrix and precortex (at protein level) (PubMed:12705872, PubMed:26173648). Expressed in the brain, muscle, kidney, pancreas, spleen, thymus and weakly expressed in the lung, liver, heart and placenta (PubMed:12705872). Highly expressed in skin, testis and prostate; less in salivary gland (PubMed:12648213). Expressed in the esophageal mucosa (PubMed:24220297) |
| DSPP | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| DUX4 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human testis full-length protein is supported, while healthy muscle produces a shorter splice form and FSHD muscle produces full-length protein. Sequence/isoform-specific clearance is needed. |
| ENSG00000236543 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ENSG00000284779 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| ENSG00000285566 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| ENSG00000286224 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| ENSG00000287585 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| FAM240B | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| FAM246A | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| FAM246B | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A10 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A13 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A14 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A16/7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A18/9 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A22 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A23 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A26 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| FAM90A3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FAM90A9 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| FOXE3 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Transcription factor that controls lens epithelial cell growth through regulation of proliferation, apoptosis and cell cycle (PubMed:22527307, PubMed:25504734). During lens development, controls the ratio of the lens fiber cells to the cells of the anterior lens epithelium by regulating the rate of proliferation and differentiation (By similarity). Controls lens vesicle closure and subsequent separation of the lens vesicle from ectoderm (By similarity). Controls the expression of DNAJB1 in a pathway that is crucial for the development of the anterior segment of the eye (PubMed:27218149) Anterior segment dysgenesis 2: A form of anterior segment dysgenesis, a group of defects affecting anterior structures of the eye including cornea, iris, lens, trabecular meshwork, and Schlemm canal. Anterior segment dysgeneses result from abnormal migration or differentiation of the neural crest derived  |
| FOXI2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported enriched RNA in a defined adult human olfactory population. |
| FRMD7 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expressed in liver, kidney, pancreas and at low levels in brain and heart. Expressed in embryonic brain and developing neural retina |
| FUT5 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| GALR3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| GFRA4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| GFRAL | Other:Tumor-supported | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Normal human brainstem neurons express GFRAL. |
| GJE1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| GLYATL3 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| GOLGA6L26 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| GPR142 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. The cited physiological study is largely rodent; the human islet annotation requires a human-specific protein/localization assessment. |
| GPR148 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expression restricted to nervous system and testis. Is also detected in several tumors types, most notably prostate cancer |
| GPR152 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| GPR31 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Human intestinal cDC1s express GPR31; human cell models support receptor-dependent antigen-sampling function. |
| GPR33 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| GPX6 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| GSTT2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| GSX2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| H2BW2 | Other:Diffuse | Other:Diffuse | The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe. |
| HES3 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| HOXC12 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| HOXC13 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| HOXD12 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| HTR3C | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expressed in many tissues including adult brain, colon, intestine, lung, muscle and stomach as well as fetal colon and kidney |
| HTR3D | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expressed in liver, as well as fetal and adult colon and kidney |
| HUS1B | Other:Diffuse | Other:Diffuse | The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe. |
| IFNA1 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA10 | Other:Low | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA13 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA14 | Other:Inducible | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA16 | Other:Low | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA17 | Other:Inducible | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA2 | Other:Inducible | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA21 | Other:Weak | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA4 | Other:Tumor-supported | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA5 | Other:Inducible | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA6 | Other:Low | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA7 | Other:Inducible | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNA8 | Other:Low | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNB1 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNL2 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human plasma contains IFNL2 during immune activation; determine producing cell and assay specificity relative to IFNL3 before clearance. |
| IFNL3 | Other:Weak | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IFNL4 | Other:Weak | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Activated primary human hepatocytes from donors carrying the protein-producing allele express endogenous IFNL4 protein. |
| IFNW1 | Other:Tumor-supported | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Reported secreted normal immune-cell protein, increased after activation. |
| IL1F10 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| IL25 | Other:Weak | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Normal human bronchial epithelial cells contain IL-25 protein at rest, with release during activation. |
| IL36B | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| IL9 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Normal human activated T cells produce IL-9 protein. |
| ITIH6 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| KCNA10 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human evidence is in situ RNA; antibody localization is rat and must not be called human protein. |
| KCNJ18 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Gene-specific normal human muscle Northern blot and channel biology warrant a hold despite below-ceiling panel RNA; patient mutations are not therapeutic injury. |
| KCNK18 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human and mouse Treg findings must be separated; thymic permissiveness does not establish absence in circulating or tissue Tregs. |
| KDM4F | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KIR3DL3 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Specific antibody flow cytometry detects KIR3DL3 on rare normal lymphocytes, including peripheral blood and mucosal populations. |
| KRT28 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRT37 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRT38 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRT39 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRT82 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRT9 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP1-4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP10-1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-11 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-12 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-5 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-6 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-7 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-8 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP10-9 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP12-3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP12-4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| KRTAP13-1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP13-3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP13-4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP15-1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP19-1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP19-2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP19-4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP19-6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP19-7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP19-8 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP20-1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP20-2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| KRTAP20-3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP20-4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| KRTAP21-1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| KRTAP21-3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP22-1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP22-2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP23-1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP25-1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP29-1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP4-1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| KRTAP4-16 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP4-5 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-10 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-11 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-5 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-6 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-7 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP5-9 | Other:Diffuse | Other:Diffuse | The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe. |
| KRTAP6-1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP6-2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP6-3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP9-1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| KRTAP9-6 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| KRTAP9-7 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| LACTBL1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| LCE3B | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| LCE3C | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| LCE4A | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| LCE7A | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| LENEP | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human lens-cell expression/promoter evidence and animal localization require adult human protein adjudication. |
| MAGEB17 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| MAGEB6B | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| MARCOL | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| MAS1L | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| MDFIC2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| MLNR | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Normal-appearing human gut shows motilin-receptor protein in muscle and cholinergic enteric neurons, with orthogonal functional responses. |
| MRGPRG | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| MRGPRX4 | Other:Tumor-supported | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Human normal DRG neurons show localized RNA and protein; approximately 6–8% of neurons are positive. |
| MROH2A | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| MSLNL | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| MTCP1 | Other:Diffuse | Other:Diffuse | The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe. |
| MYMK | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Myoblast-specific protein that mediates myoblast fusion, an essential step for the formation of multi-nucleated muscle fibers (PubMed:28681861). Actively participates in the membrane fusion reaction by mediating the mixing of cell membrane lipids (hemifusion) upstream of MYMX. Acts independently of MYMX (By similarity). Involved in skeletal muscle regeneration in response to injury by mediating the fusion of satellite cells, a population of muscle stem cells, with injured myofibers (By similarity). Also involved in skeletal muscle hypertrophy, probably by mediating the fusion of satellite cells with myofibers (By similarity) Carey-Fineman-Ziter syndrome 1: An autosomal recessive multisystem disorder characterized by hypotonia, bilateral congenital facial palsy with impairment of ocular abduction (Moebius sequence), micrognathia, glossoptosis and high-arched or cleft palate (Pierre Robin  |
| NBPF4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| NCR2 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. NKp44 is a surface receptor of activated human normal lymphocytes. |
| NEU2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| NKX1-1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human CD34-positive progenitor transcripts include this historical homeobox sequence; exact contemporary identity and adult protein need resolution. |
| NKX2-6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| NLRP2B | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Widely detected RNA and activation in human monocytic cells warrant review; transfected protein function alone does not establish primary normal-cell protein. |
| NOBOX | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human RT-PCR detects pancreas as well as ovary/testis; oocyte expression does not resolve pancreatic safety. |
| NPBWR1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Found in cerebellum and frontal cortex. Detected at high levels in hippocampus, amygdala and trachea; at moderate levels in fetal brain, pituitary gland and prostate. Not in caudate, accumbens, kidney or liver. Also detected at high levels in lung carcinoma |
| NPS | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Ligand for NPSR1 (PubMed:15312648, PubMed:15947423, PubMed:16720571, PubMed:16790440). Binds to its receptor with nanomolar affinity and initiates a G(q)/GNAQ-dependent phospholipase C-activating signaling pathway. This results in Ca(2+) mobilization from intracellular stores and increased intracellular Ca(2+) levels (PubMed:15312648, PubMed:15947423, PubMed:16720571, PubMed:16790440, PubMed:25714705, PubMed:26865629). In addition to this pathway, NPS binding to its receptor activates cAMP/PKA signal transduction (PubMed:26865629). Finally, both pathways converge to activate ERK1/ERK2 phosphorylation and signaling cascade (PubMed:26865629). Modulates arousal and anxiety. May play an important anorexigenic role (By similarity)  |
| NPY4R2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OC90 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Major protein of the otoconia, a calcium carbonate structure in the saccule and utricle of the ear. Together with OTOL1, acts as a scaffold for otoconia biomineralization: sequesters calcium and forms interconnecting fibrils between otoconia that are incorporated into the calcium crystal structure. Together with OTOL1, modulates calcite crystal morphology and growth kinetics. It is unlikely that this protein has phospholipase A2 activity  |
| OCM2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OCSTAMP | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OOSP1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OOSP3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR10A2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR10A3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR10A4 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR10A5 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR10A6 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR10A7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10AC1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR10AG1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10C1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10D3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10G6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10G7 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR10G8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10G9 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10H2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10H3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10H4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10J3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR10J4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR10J5 | Other:Low | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Human vascular and endothelial OR10J5 expression has orthogonal functional evidence. |
| OR10K2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR10P1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR10Q1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10S1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10T2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10V1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10W1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR10X1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR11G2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR11H1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR11H12 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR11H2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR11H4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR11H6 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR11H7 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR11L1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR12D1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR12D2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR12D3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR13C2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR13C3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR13C4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR13C5 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR13C7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR13C8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR13C9 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR13D1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR13F1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR13G1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR13H1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR14A16 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR14C36 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR14K1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR1A2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR1B1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR1C1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR1D5 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR1E1 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR1E2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR1G1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR1J1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR1J4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR1K1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR1L1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR1L3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR1L4 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR1L6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR1N2 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR1P1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR1S1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR1S2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2A12 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR2A14 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2A2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2A5 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR2AE1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2AG1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2AJ1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR2B2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR2B3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2B6 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. A raw non-permitted normal observation exceeds 30 nTPM but does not meet the source support rule. |
| OR2B8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2C1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2D2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2D3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2F2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR2G2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2G3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2J1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2J2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2K2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2M2 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR2M5 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR2M7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2S2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR2T1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2T2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2T27 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2T29 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2T3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2T33 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2T34 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR2T35 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2T4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2T5 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR2T6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2T7 | Other:Reproductive | CTA | Independent cancer presence, reproductive presence, concentration, normal safety and required coverage pass. Positive cancer RNA below 2 TPM now establishes presence; it does not establish therapeutic abundance. |
| OR2V2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR2W1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2Y1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR2Z1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR3A3 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4A15 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4A16 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4A47 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4A5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4A8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4B1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C11 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C12 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C13 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C15 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C16 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C45 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C46 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4C6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4D10 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4D11 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4D2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4D5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4D6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4E1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR4E2 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR4F15 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4F16/29/3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR4F21 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR4F4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4F5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4F6 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4K1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4K13 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4K14 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4K15 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4K2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4K3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4K5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4L1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4M2B | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR4N5 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4P4 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR4Q2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4Q3 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR4S1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4S2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4X1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR4X2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51A2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51A4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51A7 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR51B2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR51B4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR51B6 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR51C1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR51D1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR51F1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51F2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51G1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51G2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR51I1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR51I2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51L1 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR51Q1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51S1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51T1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR51V1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52A5 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52B2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR52B4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR52B6 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR52D1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR52E1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR52E2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR52E5 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52E6 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR52E8 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR52H1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52I1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52I2 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52J3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR52L1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR52M1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR52N2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR52P1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR52R1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR52W1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR56A1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR56A5 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR56B1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR56B2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR56B4 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR5AC1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5AC2 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR5AK2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR5AK3P | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR5AL1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5AN1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5AP2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5AR1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5AU1 | Other:Diffuse | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5B12 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5B17 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5B2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5B21 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR5B3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5BS1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5C1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR5D13 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5D14 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5D16 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5D18 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR5D3P | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5F1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5H14 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR5H15 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5H2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5H6 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR5H8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5I1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5J2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5K1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5K2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR5K3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5K4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5L1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5L2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5M1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5M10 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5M11 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR5M3 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5M8 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5M9 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5T1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5T2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR5T3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR5W2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6A2 | Other:Weak | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Human OR6A2 participates in macrophage responses; the study separately tests mouse Olfr2. Reported RNA in normal human olfactory sensory neurons. |
| OR6B1 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR6C1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR6C2 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. A raw non-permitted normal observation exceeds 30 TPM but does not meet the source support rule. |
| OR6C3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6C4 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR6C6 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR6C65 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6C68 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR6C70 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR6C75 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR6C76 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR6K2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6K6 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| OR6M1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR6P1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6Q1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6S1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR6T1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR6V1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR6X1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR7A10 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR7C2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR7G1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR7G2 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR7G3 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR8B12 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR8B4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8D1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR8D2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8G2P | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8G3 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OR8H1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR8H2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR8H3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8I2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8J1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8J2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8J3 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR8K1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8K3 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR8K5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8S1 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR8U1 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR8U3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| OR9A1P | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| OR9G1 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR9G4 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR9I1 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR9K2 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OR9Q2 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported RNA in normal human olfactory sensory neurons. |
| OTOL1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Collagen-like protein specifically expressed in the inner ear, which provides an organic scaffold for otoconia, a calcium carbonate structure in the saccule and utricle of the ear. Acts as a scaffold for biomineralization: sequesters calcium and forms interconnecting fibrils between otoconia that are incorporated into the calcium crystal structure. Together with OC90, modulates calcite crystal morphology and growth kinetics  |
| OTOP1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| OTOR | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human fetal cochlear cloning and human/avian Northern evidence require stage-specific adult human localization. |
| OXGR1 | Other:Diffuse | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Normal human kidney RNA is reported; protein abundance and cellular localization remain unresolved in this record. |
| P2RY4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| PANX3 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| PCDHB1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Reported enriched RNA in a defined adult human olfactory population. |
| PERCC1 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human congenital intestinal disease plus mouse expression/function identifies a coverage concern; it does not establish adult human protein abundance. |
| PLA2G2E | Other:Disposable | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human normal-organ RNA contradicts a Disposable-only interpretation; activated macrophage localization in this paper is mouse. |
| POTEA | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| POTEB/2 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| POTED | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| POTEJ | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| POU4F3 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human retina supplied transcript sequences; the described protein-localization experiments are other vertebrates. Adult human sensory production needs direct adjudication. |
| PPIAL4A | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Highly expressed in brain, ovary and mammary gland. Moderately expressed in lung, salivary gland, kidney, skin, adipose tissue, intestine and spleen. Weakly expressed in skeletal muscle, liver and stomach. Expressed in pleiomorphic and undifferentiated liposarcomas, osteosarcomas and breast carcinomas |
| PPIAL4E/F | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| PRAMEF33 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| PRAMEF5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| PROKR1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| PRR23C | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| PWWP4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| RAB6D | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| RLN3 | Other:Normal | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. Postmortem human brainstem contains relaxin-3-positive neurons in the nucleus incertus. |
| SCGB1C1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| SCGB1C2 | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. A raw non-permitted normal observation exceeds 30 nTPM but does not meet the source support rule. |
| SCN10A | Other:Tumor-supported | Unsafe | Established non-permitted normal production supplies exclusion evidence; this does not assert a known clinical injury. SCN10A is a marker of multiple normal human sensory-neuron populations with spatial validation. |
| SCYGR1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| SCYGR10 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| SCYGR2 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| SCYGR3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SCYGR5 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| SCYGR6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SCYGR7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SCYGR8 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| SEBOX | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| SH2D4B | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| SLC2A7 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human intestinal Northern RNA is positive, but Western/IHC in this study is rat. Adult human protein needs resolution. |
| SLC35G4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SMIM28 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| SPATA31A3 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| SPATA31A5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SPATA31A6 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SPATA31A7 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| SPDYE13/4/5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| SRMS | Other:Diffuse | Other:Diffuse | The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe. |
| SSTR4 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Specifically expressed in fetal and adult brain, lung tissue, stomach, and in lesser quantities in the kidney, pituitary and adrenals |
| SUPT20HL1 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| TAAR2 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| TAAR5 | Other:Low | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Expressed almost exclusively in skeletal muscle and selected areas of the brain, such amygdala, hippocampus, caudate nucleus, thalamus and hypothalamus (PubMed:9464258). Weak expression is also find in substantia nigra (PubMed:9464258) |
| TAAR6 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human normal-brain RT-PCR requires compartment and protein assessment; association with schizophrenia is not itself a toxicity finding. |
| TAAR8 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Mixed human/rodent receptor-family study supports a normal-expression concern; no cross-species or paralog substitution is accepted. |
| TAAR9 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TAF11L2 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| TAF11L3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TAF11L4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TAS1R2 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TAS2R13 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| TAS2R16 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TAS2R41 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| TAS2R42 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| TAS2R43 | Other:Normal | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Normal human airway bitter-receptor localization raises a relevant safety concern; confirm gene-specific reagent and receptor identity. |
| TAS2R50 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| TAS2R60 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| TAS2R7 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TAS2R8 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| TAS2R9 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| TCAF2C | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| TCHHL1 | Other:Disposable | Other:Disposable | Appreciable expression is confined to permitted Disposable contexts under the measured screen; no qualifying CTA identity. |
| TECTB | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. One of the major non-collagenous components of the tectorial membrane (By similarity). The tectorial membrane is an extracellular matrix of the inner ear that covers the neuroepithelium of the cochlea and contacts the stereocilia bundles of specialized sensory hair cells. Sound induces movement of these hair cells relative to the tectorial membrane, deflects the stereocilia and leads to fluctuations in hair-cell membrane potential, transducing sound into electrical signals  |
| TEX13B | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| TGM6 | Other:Tumor-supported | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Catalyzes the cross-linking of proteins and the conjugation of polyamines to proteins Spinocerebellar ataxia 35: A form of spinocerebellar ataxia, a clinically and genetically heterogeneous group of cerebellar disorders. Patients show progressive incoordination of gait and often poor coordination of hands, speech and eye movements, due to degeneration of the cerebellum with variable involvement of the brainstem and spinal cord. SCA35 patients commonly show upper limb involvement and torticollis. There is no cognitive impairment. |
| TOMT | Other:Weak | Unassessed | Normal evidence is unresolved; review before direct targeting or shared-sequence clearance. Human deafness and alternative-reading-frame protein evidence warrant review of adult sensory cells and the exact selected protein. |
| TRIM51G | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TRIM64C | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| TTLL8 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| UCP1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| USP17L1 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| USP17L10 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L11 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| USP17L12 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L13 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| USP17L17 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L18 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L19 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L20 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| USP17L21 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L23 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L3 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L4 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| USP17L5 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| USP17L7 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| USP17L8 | Other:Tumor-supported | Other:Tumor-supported | Cancer evidence is substantial, but reproductive identity and tumor restriction are not established. Not a direct target. |
| VHLL | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| VN1R2 | Other:Diffuse | Other:Diffuse | The complete adult expression profile is too dispersed for CTA identity; normal expression does not establish Unsafe. |
| VN1R4 | Other:Low | Other:Low | Measured normal RNA stays below the descriptive floor, with no detected normal protein. Not a demonstrated absence or a direct target. |
| VN1R5 | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| WNT8B | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ZNF648 | Other:Normal | Other:Normal | Supported normal expression is below the tissue-specific safety ceilings; reproductive identity does not qualify. |
| ZNF705B | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |
| ZNF705D | Other:Weak | Other:Weak | No qualifying reproductive route; observations lack sufficient abundance or support. Retain raw values and assay gaps. |

## Reviewed source findings

### IFNA1, IFNA13 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA2 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA4 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA5 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA6 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA7 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA8 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA10 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA14 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA16 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA17 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNA21 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNB1 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNL1 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNL3 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### IFNW1 — exclude

Reported secreted normal immune-cell protein, increased after activation.

**Context:** blood; activated plasmacytoid dendritic-cell cultures; activated. **Method:** secretome mass spectrometry.

**Consequence:** Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background.

**Limitations:** Relative protein signal is not TPM or absolute concentration. Culture composition and highly homologous interferons limit cell/paralog attribution; IFNA1/13 is explicitly shared. No therapeutic injury or peptide presentation demonstrated.

[Proteomic and Single-Cell Transcriptomic Dissection of Human Plasmacytoid Dendritic Cell Response to Influenza Virus](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) — Table 1; protein identification and secretome results.

### GFRAL — exclude

Normal human brainstem neurons express GFRAL.

**Context:** brainstem; area postrema and nucleus tractus solitarius neurons; normal. **Method:** localized protein and expression analysis.

**Consequence:** Localized normal neuronal protein outweighs a low bulk-panel maximum.

**Limitations:** Species-specific findings are distinguished from animal functional experiments. No TCR toxicity demonstrated.

[GFRAL is the receptor for GDF15 and the ligand promotes weight loss in mice and nonhuman primates](https://pubmed.ncbi.nlm.nih.gov/28846097/) — Human/mouse neuronal localization; abstract and expression results.

### SCN10A — exclude

SCN10A is a marker of multiple normal human sensory-neuron populations with spatial validation.

**Context:** dorsal root ganglion; nociceptive sensory neurons; normal. **Method:** single-nucleus RNA and spatial RNA validation.

**Consequence:** A rare normal neural population is not represented adequately by the low routine-panel bulk maximum.

**Limitations:** Human cell localization and established sensory-channel biology support screening exclusion; native atlas measurements are not converted to TPM and no clinical cytotoxicity is asserted.

[Single-nucleus transcriptomic analysis of human dorsal root ganglion neurons](https://pmc.ncbi.nlm.nih.gov/articles/PMC8626086/) — Figures 1–5; human DRG population and in situ analyses.

### MRGPRX4 — exclude

Human normal DRG neurons show localized RNA and protein; approximately 6–8% of neurons are positive.

**Context:** dorsal root ganglion; sensory neurons; normal. **Method:** in situ RNA and immunofluorescence.

**Consequence:** Direct normal neuronal protein evidence blocks unrestricted targeting.

**Limitations:** Percentage of positive neurons is not TPM or a toxicity probability; functional inhibition differs from killing expressing cells.

[MRGPRX4 is a bile acid receptor for human cholestatic itch](https://elifesciences.org/articles/48431) — Figure 4 and normal human DRG expression results.

### NCR2 — exclude

NKp44 is a surface receptor of activated human normal lymphocytes.

**Context:** blood; activated NK cells and a subset of gamma-delta T cells; activated. **Method:** surface protein, transcript and molecular characterization.

**Consequence:** Normal activation does not exempt immune cells from safety.

**Limitations:** In vitro activation establishes a normal cell state, not resting expression or clinical injury.

[NKp44, a triggering receptor involved in tumor cell lysis by activated human natural killer cells, is a novel member of the immunoglobulin superfamily.](https://pubmed.ncbi.nlm.nih.gov/10049942/) — Abstract and surface-expression characterization.

### DCANP1 — exclude

Human DCNP1/DCANP1 protein localizes in dendritic cells; normal-tissue immunostaining supports production.

**Context:** blood and lymphoid tissues; dendritic cells; normal. **Method:** immunofluorescence, IHC and RNA.

**Consequence:** Normal antigen-presenting-cell production is relevant to shared-sequence safety.

**Limitations:** Primary antibody evidence predates modern HPA grades; nuclear protein detection does not demonstrate a presented epitope.

[Identification and immunocytochemical analysis of DCNP1, a dendritic cell-associated nuclear protein.](https://pubmed.ncbi.nlm.nih.gov/11798177/) — Abstract; immunofluorescence and tissue IHC results.

### OR10J5 — exclude

Human vascular and endothelial OR10J5 expression has orthogonal functional evidence.

**Context:** aorta and coronary artery; vascular endothelial cells; normal. **Method:** protein/expression localization and functional knockdown.

**Consequence:** Adult vascular expression is not made permissible by the olfactory receptor name.

**Limitations:** HUVEC experiments are an extraembryonic vascular model; adult aortic/coronary findings are distinct. No therapeutic toxicity is demonstrated.

[Expression of human olfactory receptor 10J5 in heart aorta, coronary artery, and endothelial cells and its functional role in angiogenesis.](https://pubmed.ncbi.nlm.nih.gov/25791473/) — Expression results and endothelial functional experiments.

### OR6A2 — exclude

Human OR6A2 participates in macrophage responses; the study separately tests mouse Olfr2.

**Context:** blood and vasculature; macrophages; activated. **Method:** human expression and receptor-dependent macrophage response.

**Consequence:** Normal macrophage production is safety evidence independent of the low resting panel.

**Limitations:** Mouse atherosclerosis experiments are not human toxicity evidence. The human finding is kept separate from mouse orthology.

[Olfactory receptor 2 in vascular macrophages drives atherosclerosis by NLRP3-dependent IL-1 production.](https://pubmed.ncbi.nlm.nih.gov/35025664/) — Human macrophage results; species-specific receptor analysis.

### BPY2, BPY2B, BPY2C — reproductive

VCY2/BPY2 protein is detected in normal human male germ cells.

**Context:** testis; spermatogonia, spermatocytes and round spermatids; normal. **Method:** Western blot and IHC.

**Consequence:** Reliable reproductive protein is an independent presence route; normal safety and cancer presence still apply.

**Limitations:** Antibody evidence supports the identical selected reference protein shared by BPY2/B/C; it does not distinguish locus transcription or establish clinical safety.

[Specific expression of VCY2 in human male germ cells and its involvement in the pathogenesis of male infertility.](https://pubmed.ncbi.nlm.nih.gov/12724276/) — Normal-spermatogenesis biopsy Western blot and cellular localization.

### NOBOX — review

Human RT-PCR detects pancreas as well as ovary/testis; oocyte expression does not resolve pancreatic safety.

**Context:** pancreas; unspecified pancreatic cells; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human RT-PCR detects pancreas as well as ovary/testis; oocyte expression does not resolve pancreatic safety. No native TPM or clinical injury inferred.

[cDNA cloning and expression of the human NOBOX gene in oocytes and ovarian follicles.](https://pubmed.ncbi.nlm.nih.gov/16597639/) — Abstract and cited tissue/function results.

### KCNA10 — review

Human evidence is in situ RNA; antibody localization is rat and must not be called human protein.

**Context:** kidney and vasculature; proximal tubular, endothelial and smooth-muscle cells; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human evidence is in situ RNA; antibody localization is rat and must not be called human protein. No native TPM or clinical injury inferred.

[Expression of KCNA10, a voltage-gated K channel, in glomerular endothelium and at the apical membrane of the renal proximal tubule.](https://pubmed.ncbi.nlm.nih.gov/12444201/) — Abstract and cited tissue/function results.

### LENEP — review

Human lens-cell expression/promoter evidence and animal localization require adult human protein adjudication.

**Context:** lens; lens epithelial cells; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human lens-cell expression/promoter evidence and animal localization require adult human protein adjudication. No native TPM or clinical injury inferred.

[Functional analysis of the promoter and chromosomal localization for human LEP503, a novel lens epithelium gene.](https://pubmed.ncbi.nlm.nih.gov/11376938/) — Abstract and cited tissue/function results.

### PLA2G2E — review

Human normal-organ RNA contradicts a Disposable-only interpretation; activated macrophage localization in this paper is mouse.

**Context:** brain, heart and lung; unresolved cell type; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human normal-organ RNA contradicts a Disposable-only interpretation; activated macrophage localization in this paper is mouse. No native TPM or clinical injury inferred.

[Structures, enzymatic properties, and expression of novel human and mouse secretory phospholipase A(2)s.](https://pubmed.ncbi.nlm.nih.gov/10681567/) — Abstract and cited tissue/function results.

### TAAR6 — review

Human normal-brain RT-PCR requires compartment and protein assessment; association with schizophrenia is not itself a toxicity finding.

**Context:** brain; unresolved neuronal population; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human normal-brain RT-PCR requires compartment and protein assessment; association with schizophrenia is not itself a toxicity finding. No native TPM or clinical injury inferred.

[Polymorphisms in the trace amine receptor 4 (TRAR4) gene on chromosome 6q23.2 are associated with susceptibility to schizophrenia.](https://pubmed.ncbi.nlm.nih.gov/15329799/) — Abstract and cited tissue/function results.

### TAAR8 — review

Mixed human/rodent receptor-family study supports a normal-expression concern; no cross-species or paralog substitution is accepted.

**Context:** kidney and amygdala; unresolved cell type; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Mixed human/rodent receptor-family study supports a normal-expression concern; no cross-species or paralog substitution is accepted. No native TPM or clinical injury inferred.

[Trace amines: identification of a family of mammalian G protein-coupled receptors.](https://pubmed.ncbi.nlm.nih.gov/11459929/) — Abstract and cited tissue/function results.

### POU4F3 — review

Human retina supplied transcript sequences; the described protein-localization experiments are other vertebrates. Adult human sensory production needs direct adjudication.

**Context:** retina and inner ear; sensory neurons and hair cells; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human retina supplied transcript sequences; the described protein-localization experiments are other vertebrates. Adult human sensory production needs direct adjudication. No native TPM or clinical injury inferred.

[The Brn-3 family of POU-domain factors: primary structure, binding specificity, and expression in subsets of retinal ganglion cells and somatosensory neurons.](https://pubmed.ncbi.nlm.nih.gov/7623109/) — Abstract and cited tissue/function results.

### KCNJ18 — review

Gene-specific normal human muscle Northern blot and channel biology warrant a hold despite below-ceiling panel RNA; patient mutations are not therapeutic injury.

**Context:** skeletal muscle; myocytes; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Gene-specific normal human muscle Northern blot and channel biology warrant a hold despite below-ceiling panel RNA; patient mutations are not therapeutic injury. No native TPM or clinical injury inferred.

[Mutations in potassium channel Kir2.6 cause susceptibility to thyrotoxic hypokalemic periodic paralysis.](https://pubmed.ncbi.nlm.nih.gov/20074522/) — Abstract and cited tissue/function results.

### PERCC1 — review

Human congenital intestinal disease plus mouse expression/function identifies a coverage concern; it does not establish adult human protein abundance.

**Context:** intestine; enteroendocrine lineage; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human congenital intestinal disease plus mouse expression/function identifies a coverage concern; it does not establish adult human protein abundance. No native TPM or clinical injury inferred.

[Noncoding deletions reveal a gene that is critical for intestinal function.](https://pubmed.ncbi.nlm.nih.gov/31217582/) — Abstract and cited tissue/function results.

### GPR142 — review

The cited physiological study is largely rodent; the human islet annotation requires a human-specific protein/localization assessment.

**Context:** pancreas and intestine; islet and enteroendocrine cells; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** The cited physiological study is largely rodent; the human islet annotation requires a human-specific protein/localization assessment. No native TPM or clinical injury inferred.

[The aromatic amino acid sensor GPR142 controls metabolism through balanced regulation of pancreatic and gut hormones.](https://pubmed.ncbi.nlm.nih.gov/30472415/) — Abstract and cited tissue/function results.

### A3GALT2 — review

A normal human dendritic-cell glycosphingolipid product is not direct measurement of the encoded enzyme; review pathway and gene attribution.

**Context:** blood; monocyte-derived dendritic cells; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** A normal human dendritic-cell glycosphingolipid product is not direct measurement of the encoded enzyme; review pathway and gene attribution. No native TPM or clinical injury inferred.

[Lack of iGb3 and Isoglobo-Series Glycosphingolipids in Pig Organs Used for Xenotransplantation: Implications for Natural Killer T-Cell Biology.](https://pubmed.ncbi.nlm.nih.gov/23378701/) — Abstract and cited tissue/function results.

### KCNK18 — review

Human and mouse Treg findings must be separated; thymic permissiveness does not establish absence in circulating or tissue Tregs.

**Context:** blood and sensory ganglia; regulatory T cells and sensory neurons; normal. **Method:** primary RNA/function evidence requiring resolution.

**Consequence:** Resolve relevant adult normal-cell production before target or shared-sequence clearance.

**Limitations:** Human and mouse Treg findings must be separated; thymic permissiveness does not establish absence in circulating or tissue Tregs. No native TPM or clinical injury inferred.

[K2P18.1 translates T cell receptor signals into thymic regulatory T cell development.](https://pubmed.ncbi.nlm.nih.gov/34702947/) — Abstract and cited tissue/function results.

### OR10A6 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR10A6.

### OR6C4 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6C4.

### OR5AU1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5AU1.

### OR5AN1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5AN1.

### OR8J3 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR8J3.

### OR3A3 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR3A3.

### OR4F6 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4F6.

### OR9K2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR9K2.

### OR11G2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR11G2.

### OR6B1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6B1.

### OR8D1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR8D1.

### OR6C76 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6C76.

### OR9I1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR9I1.

### OR6C1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6C1.

### OR56A1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR56A1.

### OR9G1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR9G1.

### OR10A2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR10A2.

### OR9G4 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR9G4.

### OR4K2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4K2.

### OR2V2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2V2.

### OR5K1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5K1.

### OR5D14 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5D14.

### OR4K13 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4K13.

### OR4N5 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4N5.

### OR5T1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5T1.

### OR52N2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR52N2.

### OR51L1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR51L1.

### OR2M2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2M2.

### OR51G2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR51G2.

### OR5M3 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5M3.

### OR6C75 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6C75.

### OR8H1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR8H1.

### OR52I1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR52I1.

### OR6A2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6A2.

### OR2F2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2F2.

### OR52I2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR52I2.

### OR4A16 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4A16.

### OR6C2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6C2.

### OR2Z1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2Z1.

### OR4K14 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4K14.

### OR52E5 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR52E5.

### OR5M10 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5M10.

### OR1C1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR1C1.

### OR8H2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR8H2.

### OR5D3P — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5D3P.

### OR7G2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR7G2.

### OR7A10 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR7A10.

### OR4D10 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4D10.

### OR2AJ1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2AJ1.

### OR2A5 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2A5.

### OR56B1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR56B1.

### OR4F15 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4F15.

### OR5B2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5B2.

### OR52H1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR52H1.

### OR6C6 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR6C6.

### OR10A5 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR10A5.

### OR1E1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR1E1.

### OR9Q2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR9Q2.

### OR51V1 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR51V1.

### OR10G7 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR10G7.

### OR5T2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR5T2.

### OR8K3 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR8K3.

### OR11H6 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR11H6.

### OR2A12 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR2A12.

### OR10K2 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR10K2.

### OR4P4 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR4P4.

### OR8B12 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR8B12.

### OR52A5 — review

Reported RNA in normal human olfactory sensory neurons.

**Context:** olfactory epithelium; olfactory sensory neurons; normal. **Method:** single-cell RNA.

**Consequence:** Rare normal-cell evidence must be resolved before interpreting low bulk RNA as shared-sequence clearance.

**Limitations:** Counts are positive neurons, not TPM, protein, donor replication or clinical toxicity. A sparse observation is a review hold rather than an automatic Unsafe finding.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 7; OR52A5.

### PCDHB1 — review

Reported enriched RNA in a defined adult human olfactory population.

**Context:** olfactory epithelium; mature olfactory neurons; normal. **Method:** single-cell RNA enrichment.

**Consequence:** A missing rare-cell panel cannot establish normal restriction.

**Limitations:** Fraction of detected cells is not TPM or protein. Confirm normal-cell role and abundance; no clinical toxicity is asserted.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 3; PCDHB1.

### CLDN25 — review

Reported enriched RNA in a defined adult human olfactory population.

**Context:** olfactory epithelium; olfactory microvillar cells; normal. **Method:** single-cell RNA enrichment.

**Consequence:** A missing rare-cell panel cannot establish normal restriction.

**Limitations:** Fraction of detected cells is not TPM or protein. Confirm normal-cell role and abundance; no clinical toxicity is asserted.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 3; CLDN25.

### FOXI2 — review

Reported enriched RNA in a defined adult human olfactory population.

**Context:** olfactory epithelium; olfactory microvillar cells; normal. **Method:** single-cell RNA enrichment.

**Consequence:** A missing rare-cell panel cannot establish normal restriction.

**Limitations:** Fraction of detected cells is not TPM or protein. Confirm normal-cell role and abundance; no clinical toxicity is asserted.

[Single cell analysis of olfactory neurogenesis and differentiation in adult humans](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) — Supplementary Table 3; FOXI2.

### TAAR5 — review

Expressed almost exclusively in skeletal muscle and selected areas of the brain, such amygdala, hippocampus, caudate nucleus, thalamus and hypothalamus (PubMed:9464258). Weak expression is also find in substantia nigra (PubMed:9464258)

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Trace amine-associated receptor 5 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/O14804/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### NPBWR1 — review

Found in cerebellum and frontal cortex. Detected at high levels in hippocampus, amygdala and trachea; at moderate levels in fetal brain, pituitary gland and prostate. Not in caudate, accumbens, kidney or liver. Also detected at high levels in lung carcinoma

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Neuropeptides B/W receptor type 1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P48145/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### SSTR4 — review

Specifically expressed in fetal and adult brain, lung tissue, stomach, and in lesser quantities in the kidney, pituitary and adrenals

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Somatostatin receptor type 4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P31391/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### FRMD7 — review

Expressed in liver, kidney, pancreas and at low levels in brain and heart. Expressed in embryonic brain and developing neural retina

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[FERM domain-containing protein 7 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6ZUT3/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### DSG4 — review

In scalp follicles, present in the inner root sheath (IRS) and all layers of the matrix and precortex (at protein level) (PubMed:12705872, PubMed:26173648). Expressed in the brain, muscle, kidney, pancreas, spleen, thymus and weakly expressed in the lung, liver, heart and placenta (PubMed:12705872). Highly expressed in skin, testis and prostate; less in salivary gland (PubMed:12648213). Expressed in the esophageal mucosa (PubMed:24220297)

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Desmoglein-4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q86SJ6/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### PPIAL4A — review

Highly expressed in brain, ovary and mammary gland. Moderately expressed in lung, salivary gland, kidney, skin, adipose tissue, intestine and spleen. Weakly expressed in skeletal muscle, liver and stomach. Expressed in pleiomorphic and undifferentiated liposarcomas, osteosarcomas and breast carcinomas

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Peptidyl-prolyl cis-trans isomerase A-like 4A — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9Y536/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### HTR3C — review

Expressed in many tissues including adult brain, colon, intestine, lung, muscle and stomach as well as fetal colon and kidney

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[5-hydroxytryptamine receptor 3C — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8WXA8/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### HTR3D — review

Expressed in liver, as well as fetal and adult colon and kidney

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[5-hydroxytryptamine receptor 3D — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q70Z44/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### CST9 — review

Expressed in heart, placenta, lung, liver, skeletal muscle and pancreas (PubMed:12535658). Not expressed in brain (PubMed:12535658). Expressed in epididymis, kidney, testis, spinal cord, and thymus with a strong expression in epididymis and kidney and a weak expression in the spinal cord and thymus (PubMed:20565543)

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Cystatin-9 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q5W186/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### GPR148 — review

Expression restricted to nervous system and testis. Is also detected in several tumors types, most notably prostate cancer

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Probable G protein-coupled receptor 148 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8TDV2/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### CHRNA9 — review

Expressed in cochlea, keratinocytes, pituitary gland, B-cells and T-cells

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Neuronal acetylcholine receptor subunit alpha-9 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9UGM1/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### FOXE3 — review

Transcription factor that controls lens epithelial cell growth through regulation of proliferation, apoptosis and cell cycle (PubMed:22527307, PubMed:25504734). During lens development, controls the ratio of the lens fiber cells to the cells of the anterior lens epithelium by regulating the rate of proliferation and differentiation (By similarity). Controls lens vesicle closure and subsequent separation of the lens vesicle from ectoderm (By similarity). Controls the expression of DNAJB1 in a pathway that is crucial for the development of the anterior segment of the eye (PubMed:27218149) Anterior segment dysgenesis 2: A form of anterior segment dysgenesis, a group of defects affecting anterior structures of the eye including cornea, iris, lens, trabecular meshwork, and Schlemm canal. Anterior segment dysgeneses result from abnormal migration or differentiation of the neural crest derived 

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Forkhead box protein E3 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q13461/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### CLRN2 — review

Plays a key role to hearing function. Required for normal organization and maintenance of the stereocilia bundle and for mechano-electrical transduction Deafness, autosomal recessive, 117: A form of non-syndromic deafness characterized by prelingual, moderate-to-profound sensorineural hearing loss. Sensorineural hearing loss results from damage to the neural receptors of the inner ear, the nerve pathways to the brain, or the area of the brain that receives sound information.

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Clarin-2 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/A0PK11/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### TGM6 — review

Catalyzes the cross-linking of proteins and the conjugation of polyamines to proteins Spinocerebellar ataxia 35: A form of spinocerebellar ataxia, a clinically and genetically heterogeneous group of cerebellar disorders. Patients show progressive incoordination of gait and often poor coordination of hands, speech and eye movements, due to degeneration of the cerebellum with variable involvement of the brainstem and spinal cord. SCA35 patients commonly show upper limb involvement and torticollis. There is no cognitive impairment.

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Protein-glutamine gamma-glutamyltransferase 6 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/O95932/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### MYMK — review

Myoblast-specific protein that mediates myoblast fusion, an essential step for the formation of multi-nucleated muscle fibers (PubMed:28681861). Actively participates in the membrane fusion reaction by mediating the mixing of cell membrane lipids (hemifusion) upstream of MYMX. Acts independently of MYMX (By similarity). Involved in skeletal muscle regeneration in response to injury by mediating the fusion of satellite cells, a population of muscle stem cells, with injured myofibers (By similarity). Also involved in skeletal muscle hypertrophy, probably by mediating the fusion of satellite cells with myofibers (By similarity) Carey-Fineman-Ziter syndrome 1: An autosomal recessive multisystem disorder characterized by hypotonia, bilateral congenital facial palsy with impairment of ocular abduction (Moebius sequence), micrognathia, glossoptosis and high-arched or cleft palate (Pierre Robin 

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Protein myomaker — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/A6NI61/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### OTOL1 — review

Collagen-like protein specifically expressed in the inner ear, which provides an organic scaffold for otoconia, a calcium carbonate structure in the saccule and utricle of the ear. Acts as a scaffold for biomineralization: sequesters calcium and forms interconnecting fibrils between otoconia that are incorporated into the calcium crystal structure. Together with OC90, modulates calcite crystal morphology and growth kinetics 

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Otolin-1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/A6NHN0/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### OC90 — review

Major protein of the otoconia, a calcium carbonate structure in the saccule and utricle of the ear. Together with OTOL1, acts as a scaffold for otoconia biomineralization: sequesters calcium and forms interconnecting fibrils between otoconia that are incorporated into the calcium crystal structure. Together with OTOL1, modulates calcite crystal morphology and growth kinetics. It is unlikely that this protein has phospholipase A2 activity 

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Otoconin-90 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q02509/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### TECTB — review

One of the major non-collagenous components of the tectorial membrane (By similarity). The tectorial membrane is an extracellular matrix of the inner ear that covers the neuroepithelium of the cochlea and contacts the stereocilia bundles of specialized sensory hair cells. Sound induces movement of these hair cells relative to the tectorial membrane, deflects the stereocilia and leads to fluctuations in hair-cell membrane potential, transducing sound into electrical signals 

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Beta-tectorin — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q96PL2/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### NPS — review

Ligand for NPSR1 (PubMed:15312648, PubMed:15947423, PubMed:16720571, PubMed:16790440). Binds to its receptor with nanomolar affinity and initiates a G(q)/GNAQ-dependent phospholipase C-activating signaling pathway. This results in Ca(2+) mobilization from intracellular stores and increased intracellular Ca(2+) levels (PubMed:15312648, PubMed:15947423, PubMed:16720571, PubMed:16790440, PubMed:25714705, PubMed:26865629). In addition to this pathway, NPS binding to its receptor activates cAMP/PKA signal transduction (PubMed:26865629). Finally, both pathways converge to activate ERK1/ERK2 phosphorylation and signaling cascade (PubMed:26865629). Modulates arousal and anxiety. May play an important anorexigenic role (By similarity) 

**Context:** source-specified normal tissue; requires primary localization; normal. **Method:** reviewed database annotation and linked literature.

**Consequence:** Relevant sensory, immune or organ-function evidence requires resolution before shared-sequence clearance.

**Limitations:** A database annotation or genetic phenotype is not a quantitative adult human protein measurement; inspect exact cell, stage, species and isoform.

[Neuropeptide S — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P0C0P6/entry) — Tissue specificity/function/disease; primary PMID links retained in annotation audit.

### KRTAP5-9 — context

Restricted to hair root, not detected in any other tissues. Expressed in cuticle layers of differentiating hair follicles

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-9 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P26371/entry) — Tissue specificity; primary references 15144888.

### KRTAP12-3 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 12-3 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60328/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP12-4 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 12-4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60329/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-2 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-2 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60368/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-5 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-5 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60370/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-6 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-6 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60371/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-4 — context

Restricted to hair root, not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60372/entry) — Tissue specificity; primary references .

### KRTAP10-7 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-7 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60409/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-8 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-8 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60410/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-9 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-9 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60411/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-11 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-11 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60412/entry) — Tissue specificity; primary references 14962103|15028290.

### KRT39 — context

Expressed in skin and scalp. In the hair follicle, it is present in the upper hair cuticle and the upper cortex. Also present in the in the upper portion of beard hairs (at protein level)

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin, type I cytoskeletal 39 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6A163/entry) — Tissue specificity; primary references 15617563|17301834.

### KRTAP5-11 — context

Restricted to hair root, not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-11 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6L8G4/entry) — Tissue specificity; primary references 15144888.

### KRTAP5-7 — context

Expressed in hair root but not in skin

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-7 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6L8G8/entry) — Tissue specificity; primary references 15144888.

### KRTAP5-6 — context

Expressed in hair root and not in skin. Expressed also in liver and skeletal muscle

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-6 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6L8G9/entry) — Tissue specificity; primary references 15144888.

### KRTAP5-3 — context

Restricted to hair root, not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-3 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6L8H2/entry) — Tissue specificity; primary references 15144888.

### KRT28 — context

Strongly expressed in skin and scalp, and weak expression observed in thymus. In the hair follicle, expressed in Henle layer, Huxley layer and in the irs cuticle. Expression extends from the bulb region up to the point of differentiation into the three layers. Also present in the medulla of beard hair (at protein level)

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin, type I cytoskeletal 28 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q7Z3Y7/entry) — Tissue specificity; primary references 15617563|16874310.

### KRTAP19-1 — context

Detected in the upper portion of the hair cortex

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 19-1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8IUB9/entry) — Tissue specificity; primary references .

### KRTAP4-1 — context

Expressed in the hair follicles

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 4-1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9BYQ7/entry) — Tissue specificity; primary references 15955084.

### KRTAP4-5 — context

Expressed in the hair follicles

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 4-5 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9BYR2/entry) — Tissue specificity; primary references 15955084.

### KRT82 — context

Expressed in the hair follicle, specifically in hair shaft cuticle during the anagen phase of the hair cycle

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin, type II cuticular Hb2 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NSB4/entry) — Tissue specificity; primary references 11445569|35145093.

### KRTAP10-1 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60331/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP10-12 — context

Restricted to a narrow region of the hair fiber cuticle, lying approximately 20 cell layers above the apex of the dermal papilla of the hair root; not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 10-12 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P60413/entry) — Tissue specificity; primary references 14962103|15028290.

### KRTAP5-5 — context

Restricted to hair root, not detected in any other tissues

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-5 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q701N2/entry) — Tissue specificity; primary references 15144888.

### KRTAP13-1 — context

Weak expression seen in the late matrix and entire cortex area of the hair follicle

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 13-1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8IUC0/entry) — Tissue specificity; primary references 12359730.

### KRTAP1-4 — context

Expressed in the middle/upper portions of the hair cortex, in the region termed the keratogenous zone

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 1-4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P0C5Y4/entry) — Tissue specificity; primary references 11279113.

### KRTAP5-10 — context

Expressed in hair root but not in skin. Expressed also in brain and skeletal muscle

**Context:** hair follicle; hair-forming epithelium; normal. **Method:** curated tissue annotation linked to primary expression studies.

**Consequence:** Hair-restricted production can support Other rather than CTA; all non-permitted normal safety gates remain active.

**Limitations:** Hair evidence does not override non-hair RNA/protein findings, and shared-family assays do not establish unique locus expression.

[Keratin-associated protein 5-10 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6L8G5/entry) — Tissue specificity; primary references 15144888.

### IL9 — exclude

Normal human activated T cells produce IL-9 protein.

**Context:** blood; activated CD4 and gamma-delta T cells; activated. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[Human Vδ2 T cells are a major source of interleukin-9](https://pubmed.ncbi.nlm.nih.gov/27791087/) — Human blood T-cell RNA, flow cytometry and supernatant protein results.

### KIR3DL3 — exclude

Specific antibody flow cytometry detects KIR3DL3 on rare normal lymphocytes, including peripheral blood and mucosal populations.

**Context:** blood, lung and intestine; CD8 and gamma-delta T lymphocytes; normal. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[Polymorphic KIR3DL3 expression modulates tissue-resident and innate-like T cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC10360443/) — Figure 2 and antibody-specificity supplement.

### IL25 — exclude

Normal human bronchial epithelial cells contain IL-25 protein at rest, with release during activation.

**Context:** airway; bronchial epithelial cells; resting and activated. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[Transcription of Interleukin-25 and Extracellular Release of the Protein Is Regulated by Allergen Proteases in Airway Epithelial Cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC5455604/) — Figures 1–2.

### IFNL4 — exclude

Activated primary human hepatocytes from donors carrying the protein-producing allele express endogenous IFNL4 protein.

**Context:** liver; primary hepatocytes; activated; genotype dependent. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[IFN-λ4 Attenuates Antiviral Responses by Enhancing Negative Regulation of IFN Signaling](https://pubmed.ncbi.nlm.nih.gov/29070670/) — Endogenous protein results in primary human hepatocytes.

### CNGA2 — exclude

Normal human arterial endothelium and smooth muscle show localized CNGA2 RNA and protein.

**Context:** coronary and cerebral arteries; endothelial and vascular smooth muscle cells; normal. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[Expression of olfactory-type cyclic nucleotide-gated channel (CNGA2) in vascular tissues](https://pubmed.ncbi.nlm.nih.gov/14618336/) — Human artery in situ hybridization and IHC results.

### MLNR — exclude

Normal-appearing human gut shows motilin-receptor protein in muscle and cholinergic enteric neurons, with orthogonal functional responses.

**Context:** stomach and intestine; smooth muscle and myenteric neurons; normal-appearing surgical tissue. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[Regional- and agonist-dependent facilitation of human neurogastrointestinal functions by motilin receptor agonists](https://pmc.ncbi.nlm.nih.gov/articles/PMC3575777/) — Figure 1; macroscopically normal tissue and functional results.

### RLN3 — exclude

Postmortem human brainstem contains relaxin-3-positive neurons in the nucleus incertus.

**Context:** pons; nucleus incertus neurons; postmortem adult. **Method:** human protein localization with orthogonal expression/function evidence.

**Consequence:** Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance.

**Limitations:** Evidence establishes normal production, not peptide presentation or clinical injury. No RNA concentration is inferred from protein detection. Patient-derived normal-appearing tissues and ex vivo activation do not describe all healthy resting donors.

[Unveiling a novel memory center in human brain: neurochemical identification of the nucleus incertus, a key pontine locus implicated in stress and neuropathology](https://pmc.ncbi.nlm.nih.gov/articles/PMC11253401/) — Human pontine immunohistochemistry and neuronal colocalization.

### GPR31 — exclude

Human intestinal cDC1s express GPR31; human cell models support receptor-dependent antigen-sampling function.

**Context:** intestine; type 1 conventional dendritic cells; normal. **Method:** human tissue RNA localization and receptor-dependent cell function.

**Consequence:** Normal antigen-presenting-cell production is relevant even when bulk panels are low.

**Limitations:** Human tissue findings and derived-cell functional experiments are separate. No therapeutic injury or quantitative TPM inferred.

[The pyruvate-GPR31 axis promotes transepithelial dendrite formation in human intestinal dendritic cells](https://pubmed.ncbi.nlm.nih.gov/39432783/) — Human intestinal localization and iPSC-derived cDC1 functional results.

### IFNL2 — review

Human plasma contains IFNL2 during immune activation; determine producing cell and assay specificity relative to IFNL3 before clearance.

**Context:** blood; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Human plasma contains IFNL2 during immune activation; determine producing cell and assay specificity relative to IFNL3 before clearance. This hold does not assert clinical injury or manufacture a TPM value.

[Reduced IFNL1 and/or IFNL2, but not IFNL3 is associated with worse outcome in patients with COVID-19.](https://pubmed.ncbi.nlm.nih.gov/38953458/) — Abstract and source-specific expression results.

### SLC2A7 — review

Human intestinal Northern RNA is positive, but Western/IHC in this study is rat. Adult human protein needs resolution.

**Context:** intestine; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Human intestinal Northern RNA is positive, but Western/IHC in this study is rat. Adult human protein needs resolution. This hold does not assert clinical injury or manufacture a TPM value.

[Cloning and functional characterization of the human GLUT7 isoform SLC2A7 from the small intestine.](https://pubmed.ncbi.nlm.nih.gov/15033637/) — Abstract and source-specific expression results.

### OTOR — review

Human fetal cochlear cloning and human/avian Northern evidence require stage-specific adult human localization.

**Context:** cochlea; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Human fetal cochlear cloning and human/avian Northern evidence require stage-specific adult human localization. This hold does not assert clinical injury or manufacture a TPM value.

[A novel conserved cochlear gene, OTOR: identification, expression analysis, and chromosomal mapping.](https://pubmed.ncbi.nlm.nih.gov/10873378/) — Abstract and source-specific expression results.

### TOMT — review

Human deafness and alternative-reading-frame protein evidence warrant review of adult sensory cells and the exact selected protein.

**Context:** inner ear; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Human deafness and alternative-reading-frame protein evidence warrant review of adult sensory cells and the exact selected protein. This hold does not assert clinical injury or manufacture a TPM value.

[Mutations of LRTOMT, a fusion gene with alternative reading frames, cause nonsyndromic deafness in humans.](https://pubmed.ncbi.nlm.nih.gov/18953341/) — Abstract and source-specific expression results.

### OXGR1 — review

Normal human kidney RNA is reported; protein abundance and cellular localization remain unresolved in this record.

**Context:** kidney; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Normal human kidney RNA is reported; protein abundance and cellular localization remain unresolved in this record. This hold does not assert clinical injury or manufacture a TPM value.

[GPR99, a new G protein-coupled receptor with homology to a new subgroup of nucleotide receptors.](https://pubmed.ncbi.nlm.nih.gov/12098360/) — Abstract and source-specific expression results.

### NKX1-1 — review

Human CD34-positive progenitor transcripts include this historical homeobox sequence; exact contemporary identity and adult protein need resolution.

**Context:** bone marrow; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Human CD34-positive progenitor transcripts include this historical homeobox sequence; exact contemporary identity and adult protein need resolution. This hold does not assert clinical injury or manufacture a TPM value.

[Identification of homeobox genes expressed in human haemopoietic progenitor cells.](https://pubmed.ncbi.nlm.nih.gov/7518789/) — Abstract and source-specific expression results.

### NLRP2B — review

Widely detected RNA and activation in human monocytic cells warrant review; transfected protein function alone does not establish primary normal-cell protein.

**Context:** blood; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Widely detected RNA and activation in human monocytic cells warrant review; transfected protein function alone does not establish primary normal-cell protein. This hold does not assert clinical injury or manufacture a TPM value.

[The CLRX.1/NOD24 (NLRP2P) pseudogene codes a functional negative regulator of NF-κB, pyrin-only protein 4.](https://pubmed.ncbi.nlm.nih.gov/24871464/) — Abstract and source-specific expression results.

### DUX4 — review

Human testis full-length protein is supported, while healthy muscle produces a shorter splice form and FSHD muscle produces full-length protein. Sequence/isoform-specific clearance is needed.

**Context:** skeletal muscle and testis; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Human testis full-length protein is supported, while healthy muscle produces a shorter splice form and FSHD muscle produces full-length protein. Sequence/isoform-specific clearance is needed. This hold does not assert clinical injury or manufacture a TPM value.

[Facioscapulohumeral dystrophy: incomplete suppression of a retrotransposed gene.](https://pubmed.ncbi.nlm.nih.gov/21060811/) — Abstract and source-specific expression results.

### TAS2R43 — review

Normal human airway bitter-receptor localization raises a relevant safety concern; confirm gene-specific reagent and receptor identity.

**Context:** airway; requires exact adult cell/isoform adjudication; normal. **Method:** human RNA, functional or isoform evidence.

**Consequence:** Resolve a substantive normal-cell or sequence concern before unrestricted targeting or shared-sequence clearance.

**Limitations:** Normal human airway bitter-receptor localization raises a relevant safety concern; confirm gene-specific reagent and receptor identity. This hold does not assert clinical injury or manufacture a TPM value.

[Motile cilia of human airway epithelia are chemosensory.](https://pubmed.ncbi.nlm.nih.gov/19628819/) — Abstract and source-specific expression results.

### MLNR — context

Expressed only in thyroid, stomach, and bone marrow

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Motilin receptor — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/O43193/entry) — Tissue specificity; primary references .

### UCP1 — context

Brown adipose tissue

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Mitochondrial brown fat uncoupling protein 1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P25874/entry) — Tissue specificity; primary references 3165741.

### KRT9 — context

Expressed in the terminally differentiated epidermis of palms and soles

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Keratin, type I cytoskeletal 9 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P35527/entry) — Tissue specificity; primary references 7507869.

### AWAT1 — context

Predominantly expressed in skin, where it is limited to the sebaceous gland. Expressed in more mature, centrally located cells just before their rupture and sebum release. Also expressed in all tissues except spleen. Expressed at higher level in thymus, prostate and testis

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Acyl-CoA wax alcohol acyltransferase 1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q58HT5/entry) — Tissue specificity; primary references 15671038.

### LCE3C — context

Skin-specific. Expression was readily detected in adult trunk skin, adult arm skin, fetal skin, penal skin, vulva, esophagus and tongue. Not expressed in the cervix, rectum, lung, colon, or placenta

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Late cornified envelope protein 3C — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q5T5A8/entry) — Tissue specificity; primary references 15854049.

### LCE3B — context

Skin-specific. Expression was readily detected in adult trunk skin, adult arm skin, fetal skin, penal skin, vulva, esophagus and tongue. Not expressed in the cervix, rectum, lung, colon, or placenta

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Late cornified envelope protein 3B — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q5TA77/entry) — Tissue specificity; primary references 15854049.

### LCE4A — context

Skin-specific. Expression was readily detected in adult trunk skin, adult arm skin, fetal skin, penal skin, vulva, esophagus and tongue. Not expressed in the cervix, rectum, lung, colon, or placenta

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Late cornified envelope protein 4A — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q5TA78/entry) — Tissue specificity; primary references 15854049.

### CYP26C1 — context

Detected in most tissues at very low level

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Cytochrome P450 26C1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6V0L0/entry) — Tissue specificity; primary references 14532297.

### ASCL4 — context

Expressed in skin. 7-fold higher expression in fetal skin than in adult skin. Weak expression also detected in fetal lung, aorta and brain, and in adult stomach, kidney, ovary and breast

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Achaete-scute homolog 4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6XD76/entry) — Tissue specificity; primary references 15475265.

### DGAT2L6 — context

Expressed in all tissues tested except pancreas

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Diacylglycerol O-acyltransferase 2-like protein 6 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6ZPD8/entry) — Tissue specificity; primary references 15671038.

### IL1F10 — context

Expressed in fetal skin, spleen and tonsil. Expressed mostly in the basal epithelia of skin and in proliferating B-cells of the tonsil

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Interleukin-1 family member 10 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8WWZ1/entry) — Tissue specificity; primary references .

### WNT8B — context

Expression is restricted to the brain, and more specifically to the forebrain

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Protein Wnt-8b — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q93098/entry) — Tissue specificity; primary references 9536085.

### GFRA4 — context

Predominantly expressed in the adult thyroid gland. Low levels also found in fetal adrenal and thyroid glands

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[GDNF family receptor alpha-4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9GZZ7/entry) — Tissue specificity; primary references 11116144.

### IL36B — context

Expression at low levels in tonsil, bone marrow, heart, placenta, lung, testis and colon but not in any hematopoietic cell lines. Not detected in adipose tissue. Expressed at higher levels in psoriatic plaques than in symptomless psoriatic skin or healthy control skin. Increased levels are not detected in inflamed joint tissue

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Interleukin-36 beta — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NZH7/entry) — Tissue specificity; primary references 20300079|21242515.

### DSPP — context

Expressed in teeth. DPP is synthesized by odontoblast and transiently expressed by pre-ameloblasts

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Dentin sialophosphoprotein — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NZW4/entry) — Tissue specificity; primary references .

### NEU2 — context

Expressed in skeletal muscle, fetal liver and embryonic carcinoma cell line NT2-D1

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Sialidase-2 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9Y3R4/entry) — Tissue specificity; primary references 10191093.

### GPX6 — context

Expressed in olfactory epithelium and embryos

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Glutathione peroxidase 6 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P59796/entry) — Tissue specificity; primary references 12775843.

### BPIFB3 — context

Detected in nasal septal epithelium

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[BPI fold-containing family B member 3 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P59826/entry) — Tissue specificity; primary references 11971875.

### GPR33 — context

Expressed in spleen, lung, heart, liver, kidney, pancreas, thymus, gonads and leukocytes

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Probable G protein-coupled receptor 33 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q49SQ1/entry) — Tissue specificity; primary references 15987686.

### VHLL — context

Abundantly expressed in the placenta

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[von Hippel-Lindau-like protein — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q6RSH7/entry) — Tissue specificity; primary references 14757845.

### PROKR1 — context

Localizes to glandular epithelium, stroma and vascular endothelial cells of first trimester decidua (at protein level). Up-regulated in first trimester decidua when compared with non-pregnant endometrium. Expressed in the stomach, throughout the small intestine, colon, rectum, thyroid gland, pituitary gland, salivary gland, adrenal gland, testis, ovary, brain, spleen, prostate and pancreas

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Prokineticin receptor 1 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8TCW9/entry) — Tissue specificity; primary references 18339712.

### APOL5 — context

Low level of expression; detected in uterus, testis, skeletal muscle and stomach

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Apolipoprotein L5 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9BWW9/entry) — Tissue specificity; primary references .

### HUS1B — context

Expressed strongly in testis, less in spleen, thymus, prostate, colon and leukocytes

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Checkpoint protein HUS1B — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q8NHY5/entry) — Tissue specificity; primary references 11944979.

### NBPF4 — context

Expressed in testis

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[NBPF family member NBPF4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q96M43/entry) — Tissue specificity; primary references 16079250.

### POTED — context

Expressed in prostate, ovary, testis, placenta and prostate cancer cell lines. Localizes to basal and terminal prostate epithelial cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[POTE ankyrin domain family member D — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q86YR6/entry) — Tissue specificity; primary references 12475935.

### TAF11L2 — context

Expressed in fetal brain and testis

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[TATA-box binding protein associated factor 11 like protein 2 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/A6NLC8/entry) — Tissue specificity; primary references 21078170.

### TAF11L3 — context

Expressed in fetal brain and testis

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[TATA-box-binding protein-associated factor 11-like protein 3 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/A0A1W2PRV1/entry) — Tissue specificity; primary references 21078170.

### TAF11L4 — context

Expressed in fetal brain and testis

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[TATA-box-binding protein-associated factor 11-like protein 4 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/A0A1W2PPE2/entry) — Tissue specificity; primary references 21078170.

### TEX13B — context

Testis specific

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Testis-expressed protein 13B — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9BXU2/entry) — Tissue specificity; primary references .

### TAS2R7 — context

Expressed in subsets of taste receptor cells of the tongue and palate epithelium and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 7 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NYW3/entry) — Tissue specificity; primary references .

### TAS2R8 — context

Expressed in subsets of taste receptor cells of the tongue and palate epithelium and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 8 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NYW2/entry) — Tissue specificity; primary references .

### TAS2R9 — context

Expressed in subsets of taste receptor cells of the tongue and palate epithelium and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 9 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NYW1/entry) — Tissue specificity; primary references .

### TAS2R13 — context

Expressed in subsets of taste receptor cells of the tongue and palate epithelium and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 13 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NYV9/entry) — Tissue specificity; primary references .

### TAS2R16 — context

Expressed in a subset of gustducin-positive taste receptor cells of the tongue. Expressed in circumvallate papillae and testis (PubMed:16720576)

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 16 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/Q9NYV7/entry) — Tissue specificity; primary references 16720576.

### TAS2R41 — context

Expressed in subsets of taste receptor cells of the tongue and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 41 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P59536/entry) — Tissue specificity; primary references .

### TAS2R50 — context

Expressed in subsets of taste receptor cells of the tongue and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 50 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P59544/entry) — Tissue specificity; primary references .

### TAS2R60 — context

Expressed in subsets of taste receptor cells of the tongue and exclusively in gustducin-positive cells

**Context:** source-specified tissue; source-specified cells; normal. **Method:** curated annotation linked to literature.

**Consequence:** Retain measured classification; annotation alone does not rescue reproductive admission or establish normal exclusion.

**Limitations:** Annotation can combine stages, species and assays. It is not proof of absent expression elsewhere or sufficient direct protein evidence for an override.

[Taste receptor type 2 member 60 — UniProtKB/Swiss-Prot 2026_03](https://www.uniprot.org/uniprotkb/P59551/entry) — Tissue specificity; primary references .

### OR2B6 — review

A raw non-permitted normal observation exceeds 30 nTPM but does not meet the source support rule.

**Context:** blood; PBMC platelet cluster; 5 reads; source-reported. **Method:** RNA; high raw observation without sufficient recurrence/read support.

**Consequence:** Hold unresolved high raw normal expression outside Other until support and identity are adjudicated.

**Limitations:** Five reads fail the existing 20-read support rule; the high normalized raw signal remains unresolved, not zero. No protein production or therapeutic injury inferred.

[HPA v23 tissue clusters](https://v23.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip) — PBMC platelet cluster; 5 reads.

### SCGB1C2 — review

A raw non-permitted normal observation exceeds 30 nTPM but does not meet the source support rule.

**Context:** blood; PBMC platelet cluster; 5 reads; source-reported. **Method:** RNA; high raw observation without sufficient recurrence/read support.

**Consequence:** Hold unresolved high raw normal expression outside Other until support and identity are adjudicated.

**Limitations:** Five reads fail the existing 20-read support rule; the high normalized raw signal remains unresolved, not zero. No protein production or therapeutic injury inferred.

[HPA v23 tissue clusters](https://v23.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip) — PBMC platelet cluster; 5 reads.

### OR6C2 — review

A raw non-permitted normal observation exceeds 30 TPM but does not meet the source support rule.

**Context:** blood; DICE blood M2; highest observed donor; source-reported. **Method:** RNA; high raw observation without sufficient recurrence/read support.

**Consequence:** Hold unresolved high raw normal expression outside Other until support and identity are adjudicated.

**Limitations:** The supported third-highest value is 0.803293 TPM. A high isolated donor observation needs sample and mapping review; this is not a confirmed recurrent signal. No protein production or therapeutic injury inferred.

[DICE human immune-cell RNA](https://dice-database.org/) — DICE blood M2; highest observed donor.
