# How OncoRef classifies human antigen candidates

**Research catalog · Policy `oncoref.target-selection.2026-11-19.4` · Default biological sets (unfiltered)**
[Searchable protein catalog](audits/target-selection-20261110/rebuild/report.html) · [Organ & cell injury hierarchy](tissue-hierarchy.md) · [Data & downloads](antigen-data.md)

## Executive summary

OncoRef screens the human coding proteome—**20,099 coding loci** collapsed into **19,893 protein-group records**—to separate proteins with selective cancer or reproductive/developmental expression from proteins expressed across healthy adult organs. Against one shared non-thymic, non-reproductive healthy baseline, **397** protein groups qualify as **Cancer Antigens (CA)** and **455** qualify as **Reproductive Antigens (RA)**. The **327** proteins that clear both routes independently form the **Cancer–Testis Antigen (CTA)** intersection (in this catalog's broad reproductive/developmental sense); the candidate union contains **525** proteins.

| Stage or set | Protein groups | Notes |
| --- | --- | --- |
| Human coding loci tracked | 20,099 | Collapsed by identical representative amino-acid sequence into 19,893 protein groups |
| Excluded by normal-tissue screens (`Unsafe`) | 18,250 | Bulk RNA, single-cell clusters, immune, cornea, rare-cell RNA, normal IHC, or curated normal findings |
| Held or non-qualifying (`Other`) | 1,118 | Coverage/translation/rare-cell review (734), sequence containment (48), expression bounds (248), unscreened home tissue (88) |
| Shared-eligible before enrichment | 841 | Native TPM compared for 921 groups (including counterfactuals for held proteins) |
| Cancer Antigens (`CA`) | 397 | 327 CTA + 70 CA only (≥10× non-thymic healthy reference) |
| Reproductive Antigens (`RA`) | 455 | 327 CTA + 128 RA only (≥10× on the same reference) |
| Candidate antigen union (`CA ∪ RA`) | 525 | 356 previously nominated in curated CTA papers + 169 newly identified |

> **What a catalog call means—and does not mean.** Every OncoRef membership is a reproducible expression-screening decision, not a clinical target certification. Antigen processing and MHC presentation, tumor-specific isoform production, subcellular localization, and in vivo tolerability of a specific therapeutic modality each require direct experimental evidence.

**Contents:** [1. Counting unit](#1-defining-the-counting-unit-proteins-not-transcripts) · [2. Tissues & host choices](#2-tissue-framework-scope-injury-tiers-and-host-protection) · [3. Evidence sources](#3-evidence-sources-and-anatomical-coverage) · [4. Classification rules & funnel](#4-the-three-stage-classification-funnel) · [5. Published literature & studied CTAs](#5-how-published-cta-oncofetal-and-placental-lists-compare) · [6. Worked edge cases](#6-worked-case-studies-why-edge-cases-are-held-instead-of-forced) · [7. Highest tumor peaks](#7-highest-tumor-expression-by-category) · [8. Reproduction & reference index](#8-reproduction-open-limits-and-reference-index)

## 1. Defining the counting unit: proteins, not transcripts

Source gene identifiers across RNA, single-cell, and protein atlases are mapped to human **Ensembl 115** loci. Each locus receives one qualified reference translation—preferring **MANE Select** on the primary assembly, then **Ensembl canonical**—checked for coding biotype, complete start and stop codons, uninterrupted open reading frame, and absence of unresolved translation conflicts.

Because paralogous loci can encode identical proteins (such as `SPANXA1/2`, `XAGE1A/B`, or `VCY/1B`), loci sharing an identical representative **amino-acid sequence** are grouped into one protein record even when their cDNA sequences differ. Wherever sample identities are available, member-locus expression is summed **within the same sample** before comparing samples across cohorts. **563** records lack a qualified representative sequence hash; they stay visible in the catalog, and missing translation can block admission even when RNA is measured.

*Supporting detail:* [Reference selection and transcript-level evidence](audits/cta-transcript-selection-20261108/index.md) · [Sequence API, alternative translations, and peptide subtraction](antigen-sequences.md)

## 2. Tissue framework: scope, injury tiers, and host protection

Target-selection lists frequently conflate three separate questions. OncoRef keeps them orthogonal:

1. **Anatomical scope** — whether a tissue enters the healthy adult reference, supplies a reproductive/developmental RA witness, or is tracked separately (`thymic`).
2. **Adult injury tier** — the expected consequence of sustained organ or cell-type loss. A compartment can supply an RA witness while still carrying a protective rating (e.g., adult testis supplies RA witnesses while retaining a `Protected` tier; breast biopsies are partitioned into lactating vs. resting mammary using a coordinated milk-program signature: `LALBA`, `CSN1S1`, `CSN2`, `CSN3`).
3. **Host protection** — user-declared patient anatomy and tissue-loss tolerance (`Protect`, `Accept loss`, `Absent`), which filters downstream candidates and expands the peptide-subtraction background without altering biological CA/RA/CTA membership.

| Adult injury tier | Policy interpretation | Examples |
|---|---|---|
| **Vital** | Loss threatens life, essential repair, or major sensory/neurological function | Heart, neurons, corneal endothelium |
| **Protected** | Persistent loss is consequential and clinical compensation is poorly established | Salivary gland, whole skin, adult testis |
| **Manageable** | Defined loss has an established precedent for clinical compensation or replacement | Thyroid, mature B cells |
| **Disposable** | A narrow compartment carries a specific loss allowance | Cutaneous melanocytes, gallbladder, adult whole thymus |

Bulk measurements inherit the organ tier; permissive cell-type exceptions require resolved cell identity, and mixed populations retain the more protective tier.

### Adult thymus policy

Adult whole thymus is rated **Disposable** under the organ-loss policy while defined immune-cell populations retain their own injury tiers. Because medullary thymic epithelial cells (mTECs) express tissue-restricted antigens during central tolerance, bulk thymus is tracked under the `thymic` scope and excluded from both enrichment denominators (**N** and **D**). Independent normal-expression and protein screens still apply, and OncoRef records the observational association between adult thymectomy and increased mortality/cancer incidence so users can weigh thymic exposure explicitly.

### Biological membership vs. host protection

Default biological sets apply **no reproductive protection**. Host profiles (`female-patient`, `male-patient`, lactation, pregnancy, or custom compartment overrides) apply downstream: protecting testis reduces 327 biological CTAs to 54 retained candidates; pregnancy leaves 0 retained (123 unresolved due to incomplete developmental coverage, 204 failing a host requirement).

![Host choices and CTA screening outcomes](audits/target-selection-20261110/rebuild/summary/cta-tiers-and-profiles.png)

*Figure 1.* Biological CTA membership (327 proteins) is host-agnostic. Explicit host anatomy and protection choices filter the downstream research set and expand the subtraction background without changing biological identity.

*Supporting detail:* [Organ and cell injury hierarchy](tissue-hierarchy.md) · [Tissue calls and cell exceptions](tissue-tiers.md) · [Host protection and evidence gaps](host-selection.md) · [Adult thymectomy study](https://pubmed.ncbi.nlm.nih.gov/37530823/) · [Human mTEC central-tolerance study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2211762/)

## 3. Evidence sources and anatomical coverage

OncoRef pins **323 SHA-256-verified files** across **1,221 source/context rows** and **130 assessed tumor cohorts**. Incompatible assay scales are never pooled; every source stays in its native units and enters only the screening layer it can support:

| Evidence pillar | Primary resources | Native statistic & role |
|---|---|---|
| **Healthy-tissue baseline (`N`, `D`)** | GTEx v11, HPA v23 biopsies & brain, HPA/DICE/Monaco sorted immune | Within-site donor or biopsy P95 in native TPM sets the non-thymic healthy reference (`N`) and Disposable max (`D`). |
| **Normal-tissue safety screens** | HPA consensus & scRNA clusters, HPA normal IHC, CELLxGENE Census + Siletti, Collin cornea, specialized anatomy | Tiered RNA ceilings (TPM/nTPM/CP10K), reliable positive IHC, matched negative-IHC credit, and donor-replicated rare-cell UMI screens. |
| **Tumor & reproductive/embryo expression** | 130 QC-passing tumor cohorts, GTEx/HPA reproductive biopsies, Petropoulos/Tyser/Molè embryos | Individual QC-passing tumor libraries qualify CA; individual reproductive samples or supported embryo lineages qualify RA. |
| **Curated literature & clinical controls** | 149 catalog-mapped papers (CTA, oncofetal/DUX4, placental, clinical trials, normal counterevidence) | Provenance metadata and curated normal-expression exclusions; a positive nomination alone never confers admission. |

![Native-TPM tissue coverage and analytical role across sources](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage.png)

*Figure 2.* Native-TPM tissue coverage grouped by anatomical scope. Colored cells mark qualifying coverage for at least one assessed protein group; gray marks no qualifying coverage in that source. **RA witness** indicates individual reproductive samples or supported embryo-lineage means capable of qualifying RA; **Reproductive context** indicates tissue summaries kept for context only.

<details markdown="1">
<summary>Detailed source-family rules and versioned file inventory (323 files across 1,221 contexts)</summary>

| Source family | Evidence provided | Role in classification |
|---|---|---|
| **Ensembl 115, MANE & sequence registry** | Locus IDs, transcript flags, CDS qualification, sequence digests, read-through maps | Define protein-group counting units and test translation, nesting, and Unsafe 9-mer containment. |
| **GTEx v11** | Individual site TPM, donor/library IDs, site medians, QC annotations | Non-reproductive sites contribute the within-site linear donor P95 (after within-donor replicate median). Individual reproductive samples can independently witness RA. |
| **HPA bulk tissue & brain** (incl. GTEx via HPA) | Individual tissue biopsies, tissue/brain-region means, consensus nTPM, anatomical mappings | Native TPM enters the healthy reference (within-tissue biopsy P95 where available) or reproductive RA witnesses. Consensus nTPM runs its own normal gate. Redistributed GTEx overlaps GTEx v11. |
| **HPA sorted blood, DICE & Monaco** | Sorted immune populations, donor TPM, ranked values, published means, immune nTPM | HPA donor P95 and DICE single-locus P95 enter native comparison; DICE multi-locus groups use the conservative sum of member maxima when donor pairing is unavailable; Monaco retains the published mean. Also evaluated in the immune screen. |
| **HPA single-cell clusters & normal IHC** | Cell-type nTPM, staining grades, antibody IDs, reliability scores | Screen normal cell clusters and reliable positive protein findings; matched negative IHC can raise RNA ceilings. Reproductive staining/nTPM aids interpretation but never supplies a native-TPM 10× witness. |
| **CELLxGENE Census & Siletti brain** | Versioned H5AD studies, donor IDs, cell types, molecule counts | Grade sparse normal-cell expression by donor penetration, study replication, locus mapping risk, and independent confirmation. |
| **Collin cornea & conjunctiva** | Adult, developmental, and disease ocular cell counts with donor metadata | Dedicated adult corneal normal gate in counts per 10,000 (CP10K); non-adult/disease material stays scoped separately. |
| **Specialized anatomy studies** | Olfactory epithelium (Durante, Olender, Verbeurgt, MSSM), fetal dental epithelium, developmental skeletal atlas, adult osteoblast markers | Document localized expression in tissues underrepresented in bulk atlases; age, scale, and marker-list scope constrain what each study can establish (marker lists never prove tissue-wide absence). |
| **Petropoulos, Tyser & Molè embryos** | Per-embryo lineage TPM (preimplantation, Carnegie stage 7) plus day-9/day-11 hypoblast UMI/staining | Supported native TPM can qualify RA; hypoblast UMI counts and staining document developmental expression on their native scales. |
| **Tumor cohorts** | Per-sample matrices, sample/patient IDs, cross-cohort sharing, scale and QC checks | **130 assessed cohorts** in the [tumor coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv). Comparable QC-passing libraries supply CA witnesses (germ-cell tumors eligible; coordinated reproductive RNA advisory). |
| **Curated papers, proteomics & trials** | Nominated gene lists, normal counterevidence, supersession history, clinical/epitope records, CT83 RNAscope | Separate positive nominations from normal counterevidence and literature benchmarks. Curated normal findings can exclude a protein; a positive nomination alone never grants admission. |

| Inventory role | Versioned files |
| --- | --- |
| assay observations | 4 |
| cell atlas | 4 |
| comparison or sample metadata | 9 |
| curated data | 8 |
| curation input | 124 |
| native assay or supplement | 42 |
| tumor expression | 132 |

</details>

*Supporting detail:* [Registered file inventory (CSV)](antigen-whitepaper/sources.csv) · [Source-by-tissue table (CSV)](audits/target-selection-20261110/rebuild/normal-rna/source-tissues.csv) · [Native RNA panels, QC & denominators](audits/target-selection-20261110/rebuild/normal-rna/index.md) · [Download API & atlas inventory](antigen-data.md)

## 4. The three-stage classification funnel

Every protein group moves through three sequential stages—first excluding normal-tissue expression, next verifying technical eligibility and sequence uniqueness, and finally testing tumor and reproductive enrichment against the same healthy baseline:

| Stage | Core question | Decision rule | Outcomes |
|---|---|---|---|
| **1. Healthy-tissue screen** | Is the protein expressed in vital, protected, or manageable adult tissues? | Tiered bulk/cluster/immune RNA ceilings (with negative-IHC credit), adult cornea CP10K, reliable normal IHC, curated normal literature, and donor-replicated rare-cell RNA. | Confirmed normal expression → **`Unsafe`**; strong rare-cell signal with unresolved paralog attribution → **`Other:Review`** |
| **2. Eligibility & coverage** | Can we make a fair, sequence-specific like-for-like comparison? | Complete coverage across **35** healthy contexts, qualified translation, <95% 9-mer containment in `Unsafe`, target peak > healthy max, Disposable ≤ **32** TPM/nTPM, and measured adult normal home. | Gaps/holds → **`Other:Review`**, **`Other:Sequence`**, **`Other:Expression`**, or **`Other:Unscreened tissue`** |
| **3. Symmetric 10× enrichment** | Does tumor or reproductive/embryo RNA exceed the healthy baseline by ≥10×? | Evaluated independently against **H₀ = max(N, D, 1 TPM)** (both non-thymic, non-reproductive). | Tumor passes → **`CA`**; reproductive/embryo passes → **`RA`**; both pass → **`CTA`** |

> **Qualifying expression ≥ 10 × max(N, D, 1 TPM)**

<details markdown="1">
<summary>Exact quantitative gate thresholds (RNA ceilings, IHC credit, rare-cell replication, and embryo support)</summary>

* **Bulk, cluster, and immune RNA ceilings:** Base limits in TPM and nTPM are **4 for Vital, 8 for Protected, 16 for Manageable**. Matched negative IHC across routine contexts of the tissue raises the ceiling by 2× (Approved), 4× (Supported), or 8× (Enhanced), alongside limited registry clinical adjustments. Exceeding a `Vital` or `Protected` ceiling establishes `Unsafe`; `Manageable` establishes `Unsafe` only above **128**.
* **Donor-resolved rare-cell RNA:** Evaluated in native UMI units. In a `Vital` population, a donor signal requires **≥10 cells** and **≥1%** of cells with ≥2 molecules (replicated across studies), or **≥5%** in a single study; `Protected` populations double those frequency cutoffs. Replicated signals without locus mapping risk (or with independent bulk/gene-specific confirmation) establish `Unsafe`; unresolved locus attribution holds the protein in `Other:Review`.
* **Structural & home-tissue gates:** Unsupported nested sequences, read-through loci, or ≥95% 9-mer containment in `Unsafe` proteins produce `Other:Sequence`. A curated adult normal home outside the native-TPM panels (e.g., `GFY` in olfactory neurons; see Section 6.2) blocks CA as `Other:Unscreened tissue` while still allowing independent RA qualification; fetal-only homes (such as `AMELX` in developing tooth epithelium) keep their developmental scope.
* **Embryo lineage support for RA:** Well-sampled lineages require ≥20 cells, ≥5 cells at ≥1 TPM, and ≥10% positive cells; smaller lineages require ≥5 cells, ≥50% positive cells, and replication in two independent embryos of the same stage and lineage. Within-site donor/biopsy **P95** dampens isolated healthy outliers but cannot diagnose carryover or prove absence from an unmeasured rare cell type.

</details>

### Proteome-wide attrition and enrichment boundary

Of **19,893** protein groups, **18,250** are excluded as `Unsafe` in Stage 1 (most at the bulk RNA ceiling) and **1,118** land in `Other`, leaving **525** retained candidate antigens (**841** pass all shared eligibility checks; **921** have native-TPM comparisons including counterfactuals for held proteins).

![Attrition across the entire protein catalog](audits/target-selection-20261110/rebuild/summary/funnel-all.png)

*Figure 3.* Sequential proteome-wide attrition from **19,893** protein groups to **525** retained candidates (`CA ∪ RA`), assigning each excluded protein once at its first failing gate. The broken x-axis above 3,000 keeps every post-bulk filter proportional and legible. [Download exact step counts (CSV)](audits/target-selection-20261110/rebuild/summary/funnel-all.csv).

![Shared 10x enrichment rule across cancer and reproductive compartments](audits/target-selection-20261110/rebuild/summary/tumor-normal.png)

*Figure 4.* All **921** protein groups with native-TPM comparisons plotted against the shared `10 × max(N, D, 1 TPM)` boundary. Gold points above the dotted line clear the fold rule but remain held in `Other:Review` or `Other:Unscreened tissue`.

| Exclusive category | Protein records | Member loci |
| --- | --- | --- |
| CTA | 327 | 358 |
| Cancer Antigen only | 70 | 74 |
| Other:Expression | 248 | 257 |
| Other:Review | 734 | 735 |
| Other:Sequence | 48 | 50 |
| Other:Unscreened tissue | 88 | 88 |
| Reproductive Antigen only | 128 | 130 |
| Unsafe | 18,250 | 18,407 |

*Supporting detail:* [Thresholds and assay rules](target-selection.md) · [All curation figures & distributions](audits/target-selection-20261110/rebuild/summary.md)

## 5. How published CTA, oncofetal, and placental lists compare

OncoRef tracks curated literature across four distinct antigen families—**CTA / cancer-germline catalogs**, **oncofetal / DUX4 / early-embryo programs**, **placental / cancer-placenta and reproductive catalogs**, and **tumor-vasculature or focused single-gene reviews**. A published nomination is provenance metadata and never overrides a normal-tissue exclusion.

Retention varies sharply with how strictly each paper screened normal tissues: focused modern CTA and DUX4 studies (`Seager et al. 2024`, `Chew et al. 2019`, `Carter et al. 2023`, `Loriot et al. 2025`) retain **34%–67% CTAs**, whereas broad germ-cell, microarray, ZGA, and cancer-placenta lists (`Bruggeman et al. 2018`, `Smith et al. 2023`, `Jamin et al. 2021`, `Bradley et al. 2020`) are **92%–100% `Unsafe`** under adult healthy-tissue screening.

![Curated CTA, oncofetal, and placental literature sources in OncoRef](audits/target-selection-20261110/rebuild/summary/literature-sources.png)

*Figure 5.* Curated literature sources grouped by antigen class, showing mapped protein-group counts (left) and disposition share (right). [Download full source-by-source table (CSV)](audits/target-selection-20261110/rebuild/summary/literature-sources.csv).

### Overlap with rule-derived sets and published-nomination attrition

Across **2,233** published CTA/oncofetal nominations, **356** survive into the candidate union—including **252** of our **327** CTAs—while **75** CTAs are newly identified by the expression rules and **1,877** published nominations fail normal-tissue, eligibility, or enrichment screens.

![Overlap of rule-derived CA and RA sets with published nominations](audits/target-selection-20261110/rebuild/summary/venn.png)

*Figure 6.* Overlap between rule-derived Cancer Antigens (**397**), Reproductive Antigens (**455**), and published CTA/oncofetal nominations (**2,233**).

![Attrition among published CTA or oncofetal nominations](audits/target-selection-20261110/rebuild/summary/funnel-published.png)

*Figure 7.* Sequential attrition of the **2,233** published CTA/oncofetal nominations through the same screening gates, leaving **356** retained candidates. [Download step counts (CSV)](audits/target-selection-20261110/rebuild/summary/funnel-published.csv).

### Where the most-nominated individual proteins land

Ranking by **distinct curated papers nominating the protein as a CTA** (ties ordered by protein label and entity ID) gives a reproducible registry benchmark across **2,238** protein groups (**2,233** with a CTA nomination paper plus **5** additional familiar or clinical reference entries). Below are the top ten by nomination count plus five familiar anchors (`NY-ESO-1`, `PRAME`, `CT83`, `CTAG2`, `CAGE1`):

| Protein | CTA papers | Catalog | N TPM · site | RA peak TPM · context | Tumor TPM · cohort |
| --- | --- | --- | --- | --- | --- |
| MAGEA3 | 8 | Unsafe | NR | NR | NR |
| MAGEA4 | 8 | CTA | 2.95 · lymph node | 327.97 · 4 day / unspecified lineage | 3,201.68 · SCLC |
| DCAF4L2 | 7 | CTA | 0.12 · Brain_Substantia_nigra | 59.96 · Testis | 310.77 · NEC_LUNG_LARGECELL |
| MAGEA1 | 7 | CTA | 0.21 · spleen | 31.36 · Testis | 346.43 · SCLC |
| MAGEA12 | 7 | Unsafe | NR | NR | NR |
| MAGEB2 | 7 | Unsafe | NR | NR | NR |
| CT83 | 6 | CTA | 8.6 · Minor_Salivary_Gland | 445.31 · Testis | 193.17 · CESC |
| CTCFL | 6 | Unsafe | NR | NR | NR |
| DPPA2 | 6 | CTA | 0.94 · B_CELL_NAIVE | 167.76 · 7 day / epiblast | 398.83 · RT |
| GAGE1 | 6 | CTA | 0.7 · retina | 165.18 · 5 day / trophectoderm | 1,181.9 · SKCM |
| PRAME | 6 | CTA | 10.93 · Adrenal_Gland | 310.91 · Testis | 3,692.46 · SCLC |
| NY-ESO-1 | 5 | CTA | 0.84 · spleen | 16.88 · testis.244 | 944.32 · SARC_SYN |
| CAGE1 | 3 | Other:Review | 0.94 · CD4_STIM | 32.98 · testis.349 | 17.4 · T_ALL |
| CTAG2 | 3 | Other:Review | 6.92 · Bladder | 101.92 · Testis | 604.33 · OV |

*Reading this table:* **N** and the reproductive peak can originate from different assays or donors ([full context table](antigen-whitepaper/studied-proteins.csv)). **NR** ("not run") marks proteins excluded as `Unsafe` in Stage 1 before native-TPM comparison (MAGEA3, MAGEA12, MAGEB2, CTCFL); tumor or testis libraries may still contain those transcripts. **—** indicates an unavailable value.

### Studied proteins held or excluded

Of **58** distinct protein groups in the familiar or clinical reference panels, **30** are `Other` or `Unsafe`. Neither label predicts the clinical toxicity of a specific therapeutic product; each records the strongest normal-tissue or rare-cell finding observed under policy `oncoref.target-selection.2026-11-19.4`:

| Studied protein | Category | Excluding or holding gate | Strongest measured normal signal |
| --- | --- | --- | --- |
| MAGEA3 | Unsafe | HPA tissue clusters | colon, undifferentiated cells · 53.1 nTPM (ceiling 4) |
| MAGEA12 | Unsafe | Established rare-cell RNA (detection confirmed by a published gene-specific assay) | brain, pericyte · 0.40% of cells (≥2 molecules) |
| CTCFL | Unsafe | Normal bulk | white matter, bulk · 15 nTPM (ceiling 4) |
| CTAG2 | Other:Review | Unresolved rare-cell RNA (replicated normal-cell signal with unresolved locus attribution) | lung, pulmonary alveolar type 2 cell · 4.29% of cells (≥2 molecules) |
| DDX43 | Unsafe | HPA normal protein  /  HPA tissue clusters  /  Monaco immune  /  DICE immune  /  Established rare normal-cell RNA | blood, NK · 21.09 TPM (ceiling 8) |
| SPA17 | Unsafe | HPA normal protein  /  Normal bulk  /  HPA tissue clusters  /  Adult cornea  /  Established rare normal-cell RNA | lung, ciliated cells · 345.8 nTPM (ceiling 8) |
| PLAC1 | Unsafe | Established rare-cell RNA (replicated signal at a locus without mapping risk) | brain, L2/3-6 intratelencephalic projecting glutamatergic neuron · 9.17% of cells (≥2 molecules) |
| CLDN6 | Unsafe | Normal bulk  /  HPA tissue clusters | pancreas, ductal cells · 64.8 nTPM (ceiling 8) |
| DAZL | Other:Review | Unresolved rare-cell RNA (frequent single-study normal-cell signal with unresolved locus attribution) | eye, astrocyte · 5.86% of cells (≥2 molecules) |
| POTEC | Other:Review | Unresolved rare-cell RNA (replicated normal-cell signal with unresolved locus attribution) | eye, corneal endothelial cell · 12.07% of cells (≥2 molecules) |
| CAGE1 | Other:Review | Unresolved rare-cell RNA (replicated normal-cell signal with unresolved locus attribution) | brain, dopaminergic neuron · 3.29% of cells (≥2 molecules) |

*Supporting detail:* [Full studied-protein table (CSV)](antigen-whitepaper/studied-proteins.csv) · [Literature adjudication](audits/target-selection-20261110/rebuild/literature.md)

## 6. Worked case studies: why edge cases are held instead of forced

Two concrete examples show why OncoRef separates unresolved evidence or missing assay coverage from confirmed normal-tissue exclusion.

### 6.1 Paralog ambiguity in single-cell RNA: `CTAG2` (LAGE-1) → `Other:Review`

In peripheral bulk native TPM, `CTAG2` clears both enrichment rules: its highest peripheral reference is **urinary bladder** at **6.92 TPM** (10× cutoff **69.25 TPM**), testis reaches **101.92 TPM** (**14.72×**), and ovarian cancer reaches **604.33 TPM** across **375** known patients at 10×.

However, two adult organs—**lung** (alveolar type 2 cells) and **brain** (cortical glutamatergic neurons; shared Allen donor `H18.30.002` deduplicated across collections, with intratelencephalic neurons replicated in `H18.30.002` and `H20.33.014`)—meet the donor-replication rare-cell rule:

| Organ | Strongest population | Atlas collections | Donor IDs with signal (linked where known) | Cells with ≥2 molecules in top donor | Grading |
| --- | --- | --- | --- | --- | --- |
| lung | pulmonary alveolar type 2 cell | 2 | 3 | 166/3,865 (4.29%) | unresolved |
| brain | L5 extratelencephalic projecting glutamatergic cortical neuron | 2 | 2 | 70/2,201 (3.18%) | unresolved |
| eye | ciliary muscle cell | 2 | 2 | 21/1,049 (2.00%) | dismissed |

**Why it is held rather than excluded:** **50.9%** of distinct `CTAG2` transcript 31-mers are shared with paralogs (`CTAG1A`, `CTAG1B`). Neither matched bulk RNA nor a 15-sample normal-cell qPCR panel (which assayed fibroblasts and astrocytes rather than alveolar type 2 or cortical neurons) resolves whether those counts come from `CTAG2` or a paralog. Meanwhile, a 25-patient myeloma trial of a shared NY-ESO-1/LAGE-1 receptor reported no fatal serious adverse events (evidence specific to that receptor and population), and a separate heart review found discordant RNA and cardiomyocyte staining. OncoRef therefore holds `CTAG2` in **`Other:Review`** until locus-specific data resolve the implicated cell types.

*Case sources:* [Donor findings](audits/normal-cell-atlas-20261110/rare-cell/findings.csv) · [Per-donor counts](audits/normal-cell-atlas-20261110/rare-cell/donor-populations.csv.gz) · [Mapping audit](audits/normal-cell-atlas-20261110/rare-cell/gene-mapping-risk.csv.gz) · [Heart assay review](audits/ctag2-heart-evidence.md) · [HLCA](https://www.nature.com/articles/s41591-023-02327-2) · [LungMAP CellRef](https://doi.org/10.1038/s41467-023-40173-5) · [Normal-cell qPCR panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC3624013/) · [Myeloma trial follow-up](https://pubmed.ncbi.nlm.nih.gov/31289029/)

### 6.2 Unmeasured adult normal home: `GFY` (*Goofy*) → `Other:Unscreened tissue`

`GFY` (*Goofy*) regulates olfactory cilia signaling and reaches **1,575.56 TPM** in rhabdoid tumor (`RT`) over a 2.65 TPM bulk peripheral reference. Yet adult human olfactory-biopsy single-cell data detect `GFY` RNA in **94.1%** of mature and **70.3%** of immature olfactory neurons—an organ missing from the native-TPM comparison panels—and mouse *Gfy* knockout disrupts `ADCY3` localization, olfactory cilia length, and odor sensitivity (with HEK293T coexpression enhancing odorant-receptor trafficking). Because no like-for-like native-TPM olfactory denominator exists, OncoRef retains the tumor peak on the record and classifies `GFY` as **`Other:Unscreened tissue`**.

*Case sources:* [Adult human olfactory atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) · [Mouse knockout study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6619734/) · [Receptor trafficking study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12677133/) · [Cell and isoform findings](curation.md#cell-and-isoform-evidence)

## 7. Highest tumor expression by category

Across the **525** candidate union, **162** proteins reach **≥100 TPM**, **30** reach **≥1,000 TPM**, and **132** have **≥10 known independent patients** crossing the 10× cutoff (median peak tumor expression: **73.45 TPM** for CTA; **54.07 TPM** for CA only).

*Reading these tables:* Each table ranks up to ten protein groups in its exclusive category by **maximum QC-passing tumor TPM**. Because `Tumor fold` divides by `max(N, D, 1 TPM)`, fold equals `Tumor TPM` whenever both `N` and `D` are below 1 TPM (`LALBA`, `SCGB1D4`, `VCY/1B`). `Other` covers only groups that reached the quantitative comparison (high peaks such as `TBATA`, `GFY`, and `CTAG2` reflect review or unscreened-home holds). Peaks measure single-library abundance, not the fraction of malignant cells expressing the protein.

### CTA

Tumor comparison available for 327 of 327 groups.

| Rank | Protein | Peak cohort | Tumor TPM | N TPM | Reproductive peak TPM | Tumor fold | Known patients at 10× |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | LALBA | BRCA | 79,700.44 | 0.78 | 45,625.2 | 79,700.44 | 25 |
| 2 | TGM4 | PRAD | 17,629.45 | 5.67 | 22,106.5 | 3,107.11 | 94 |
| 3 | XAGE1A/B | BRCA | 4,068.64 | 5.74 | 362.73 | 709.23 | 988 |
| 4 | SPANXA1/2 | ESCA | 3,700.06 | 3 | 77.07 | 1,233.35 | 61 |
| 5 | PRAME | SCLC | 3,692.46 | 10.93 | 310.91 | 337.97 | 1,206 |
| 6 | SCGB1D4 | UCEC | 3,576.17 | 0.42 | 1,383.54 | 3,576.17 | 49 |
| 7 | MAGEA4 | SCLC | 3,201.68 | 2.95 | 327.97 | 1,083.95 | 745 |
| 8 | RLN1 | PRAD | 3,039.25 | 5.85 | 268.82 | 519.47 | 127 |
| 9 | SEMG2 | PRAD | 3,033.78 | 9.52 | 26,737.9 | 318.52 | 7 |
| 10 | VCY/1B | SKCM | 2,858.58 | 0.1 | 269.66 | 2,858.58 | 125 |

### CA only

Tumor comparison available for 70 of 70 groups.

| Rank | Protein | Peak cohort | Tumor TPM | N TPM | Reproductive peak TPM | Tumor fold | Known patients at 10× |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | AMELX | RT | 1,696.17 | 1.9 | 2.29 | 892.72 | 30 |
| 2 | HOXC12 | WILMS | 1,559.28 | 11.6 | 5.51 | 134.47 | 9 |
| 3 | RPS4Y2 | B_ALL | 1,394.88 | 12.19 | 12.74 | 114.4 | 50 |
| 4 | GALP | SARC_DSRCT | 1,332.56 | 1.02 | 1.05 | 1,308.13 | 22 |
| 5 | SP7 | SARC_OS | 590.49 | 1.92 | 7.9 | 308.19 | 222 |
| 6 | ZACN | SCLC | 526.62 | 22.49 | 21.96 | 23.41 | 1 |
| 7 | IL36B | SARC_UPS | 498.45 | 20.17 | 9.83 | 24.71 | 1 |
| 8 | SRMS | NET_RECTAL | 463.74 | 12.05 | 39.04 | 38.47 | 0 |
| 9 | CSAG3 | LUSC | 336.48 | 4.27 | 22.46 | 78.72 | 401 |
| 10 | DBX1 | SCLC | 297.51 | 0.78 | 1.65 | 297.51 | 74 |

### RA only

Tumor comparison available for 128 of 128 groups.

| Rank | Protein | Peak cohort | Tumor TPM | N TPM | Reproductive peak TPM | Tumor fold | Known patients at 10× |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | STPG3 | NEC_MERKEL | 136.45 | 16.79 | 196.89 | 8.13 | 0 |
| 2 | OR1N1 | OV | 98.33 | 0.8 | 12.16 | 98.33 | 12 |
| 3 | PABPC3 | BL | 69.66 | 6.26 | 95.08 | 11.13 | 3 |
| 4 | PABPN1L | B_ALL | 58.46 | 11.77 | 585.75 | 4.97 | 0 |
| 5 | PSG7 | NEC_MERKEL | 55.42 | 6.12 | 70.23 | 9.05 | 0 |
| 6 | KRTAP20-4 | NEC_LUNG_LARGECELL | 50.41 | 0.3 | 10.05 | 50.41 | 2 |
| 7 | OR51C1 | PRAD | 50 | 2.37 | 43.5 | 21.08 | 8 |
| 8 | TMEM262 | NEC_LUNG_LARGECELL | 47.43 | 7.03 | 123.95 | 6.74 | 0 |
| 9 | TOMM20L | T_ALL | 46.57 | 6.16 | 73.76 | 7.57 | 0 |
| 10 | ZNF878 | B_ALL | 39.39 | 3.94 | 48.31 | 10 | 0 |

### Other

Tumor comparison available for 396 of 1,118 groups.

| Rank | Protein | Peak cohort | Tumor TPM | N TPM | Reproductive peak TPM | Tumor fold | Known patients at 10× |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | TBATA | THYM | 2,250.3 | 6.26 | 116.44 | 359.45 | 101 |
| 2 | KRT9 | PCPG | 1,865.83 | 8.65 | 7.03 | 215.65 | 16 |
| 3 | DSPP | HEPB | 1,695.58 | 0.6 | 0.21 | 1,695.58 | 1 |
| 4 | GFY | RT | 1,575.56 | 2.65 | 14.05 | 594.83 | 41 |
| 5 | PANX3 | SARC_OS | 1,387.97 | 2 | 2.32 | 692.34 | 155 |
| 6 | CTAG2 | OV | 604.33 | 6.92 | 101.92 | 87.27 | 375 |
| 7 | GLYATL3 | NET_MIDGUT | 479.55 | 1.58 | 0.32 | 302.67 | 29 |
| 8 | RIPPLY3 | RT | 450.15 | 2.36 | 10.02 | 190.41 | 29 |
| 9 | VCX3B | BRCA | 420.74 | 1.33 | 181.54 | 316.01 | 87 |
| 10 | PSG3 | LUAD | 326.38 | 17.41 | 1,155.62 | 18.74 | 1 |

### Unsafe

Tumor comparison available for 0 of 18,250 groups.

*Full rankings:* [All 525 candidates with inclusive CA/RA/CTA ranks, sample/cohort IDs, folds, and patient support (CSV)](antigen-whitepaper/tumor-ranks.csv) · [Retained tumor witnesses (CSV.gz)](audits/target-selection-20261110/rebuild/tumor-witnesses.csv.gz)

## 8. Reproduction, open limits, and reference index

### Reproducing this whitepaper and catalog

```sh
# Regenerate the whitepaper, tissue hierarchy, and whitepaper CSVs
python scripts/render_antigen_whitepaper.py

# Verify all generated whitepaper artifacts byte-for-byte
python scripts/render_antigen_whitepaper.py --check

# Verify the full curation receipt and all tracked outputs
python -m oncoref.curation check
```

The whitepaper generator reads only the versioned catalog, evidence registries, and committed audit tables. Rebuilding the **underlying catalog** from scratch additionally requires the pinned upstream reference and tumor matrices documented in [Build instructions](curation.md).

### Open scientific limitations

* **Cross-cohort tumor comparability:** Tumor cohorts come from distinct studies and processing pipelines; OncoRef enforces nonnegative comparable TPM and library QC rather than cross-study batch correction.
* **Read-level attribution in rare-cell atlases:** Single-cell count matrices cannot separate paralog cross-mapping from genuine sparse transcription without read-level realignment or locus-specific assays (e.g., `CTAG2`).
* **Normal-protein replay gaps:** Hair-follicle normal-protein evidence replay remains an open curation item.
* **Bulk peaks vs. clinical targetability:** Single-sample RNA peaks do not measure intratumoral heterogeneity, translated protein abundance, peptide-MHC presentation, or modality-specific safety.

### Source and documentation index

* **Core catalog & explorer:** [Searchable protein catalog](audits/target-selection-20261110/rebuild/report.html) · [Curation summary & all figures](audits/target-selection-20261110/rebuild/summary.md) · [Classification rules & thresholds](target-selection.md) · [Build & reproduction guide](curation.md)
* **Tissue & host policy:** [Organ & cell injury hierarchy](tissue-hierarchy.md) · [Tissue calls & cell exceptions](tissue-tiers.md) · [Host protection & evidence gaps](host-selection.md) · [Sequence API & peptide subtraction](antigen-sequences.md)
* **Downloadable tables:** [Registered file inventory (`sources.csv`)](antigen-whitepaper/sources.csv) · [Studied CTA dispositions (`studied-proteins.csv`)](antigen-whitepaper/studied-proteins.csv) · [Candidate tumor rankings (`tumor-ranks.csv`)](antigen-whitepaper/tumor-ranks.csv) · [Curated literature sources by class (`literature-sources.csv`)](audits/target-selection-20261110/rebuild/summary/literature-sources.csv)
* **Audit ledgers:** [Normal RNA source coverage & QC](audits/target-selection-20261110/rebuild/normal-rna/index.md) · [Tumor cohort coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv) · [Literature adjudication](audits/target-selection-20261110/rebuild/literature.md) · [Transcript reference selection](audits/cta-transcript-selection-20261108/index.md) · [Data download API & atlas inventory](antigen-data.md)
