# How OncoRef classifies human antigen candidates

[Organ and cell injury hierarchy with supporting sources](tissue-hierarchy.md)


**Research catalog · `oncoref.target-selection.2026-11-19.4` · default biological sets**

Of **20,099 human coding loci**, OncoRef tracks **19,893 protein-group records**. **397** qualify as Cancer Antigens (CA) and **455** as Reproductive Antigens (RA). The **327** that qualify independently on both sides are Cancer–Testis Antigens (CTA) in this catalog's broad reproductive/developmental sense. The union contains **525** candidates. A catalog call describes expression and screening evidence. Antigen presentation, isoform production, tumor-cell localization and tolerance of an intervention each require evidence of their own.

## What counts as one protein?

OncoRef maps source gene identifiers to human Ensembl 115 loci. MANE Select on the primary locus is the first choice for the reference translation, followed by Ensembl canonical. Qualification checks the coding biotype, complete start and stop, an uninterrupted protein sequence and unresolved conflicting evidence. Loci with the same representative **amino-acid sequence** share one protein record. Different cDNAs can encode that same sequence. **563 records have no representative sequence hash** and remain visible in the catalog. Missing translation can block admission even when RNA is measured. [Reference selection and transcript-level evidence](audits/cta-transcript-selection-20261108/index.md)

Expression for a group uses its member loci. Where sample identities are available, the analysis adds measured loci **within the same sample** before comparing samples. The sequence API also exposes alternative translations, including sequences used for peptide subtraction. Its `matches_representative_sequence` field marks an exact amino-acid match to a catalog representative. A custom selection that omits reference hashes yields a missing match value. Neither group-level RNA nor that match identifies which proteoform a tumor produced. [Sequence API and coverage](antigen-sequences.md)

## Which tissues count, and why?

The registry records two independent properties. **Anatomical scope** determines whether a measurement can support reproductive association or enters the healthy reference. **Injury tier** describes the expected consequences of sustained tissue loss. A testis can supply an RA witness while loss of the adult organ still carries a Protected rating. Breast and prostate can also supply RA witnesses. Breast expression can be further classified by a milk-expression signature. Its label is an inference from the biopsy.

| Adult injury tier | Policy interpretation | Example |
|---|---|---|
| Vital | Loss threatens life, essential repair or major sensory or neurological function | Heart, neurons |
| Protected | Persistent loss is consequential and compensation is poorly established | Salivary gland, whole skin |
| Manageable | Defined loss has a precedent for compensation or treatment | Thyroid, mature B cells |
| Disposable | A narrow compartment has a specific loss allowance | Cutaneous melanocytes, gallbladder |

Bulk measurements inherit the organ assessment. A more permissive cell rule requires resolved cell identity. Mixed populations retain the more protective applicable concern. Unresolved anatomy and missing coverage are recorded. The defaults for both male and female hosts add **no reproductive protection**. Users may explicitly protect a compartment, record its absence or accept its loss. Host choices affect the selected subset and subtraction background. The biological memberships below use the unfiltered catalog. [Tissue calls, references and cell exceptions](tissue-tiers.md) · [Host protection and evidence gaps](host-selection.md)

## Evidence sources and what each contributes

The [download inventory](antigen-whitepaper/sources.csv) gives the selector, version, description, source URL, local or upstream asset name, and SHA-256 for **323 registered files**. Several files can come from the same study, and some are derived observations from upstream files. The [source-by-tissue table](audits/target-selection-20261110/rebuild/normal-rna/source-tissues.csv) lists **1221 source/context rows**, including tissue, finer anatomical site or cell type, and analytical role.

| Source family | Available evidence | Use in classification |
|---|---|---|
| Ensembl, MANE and the sequence registry | Locus IDs, transcript designations, CDS qualification, sequence digests, read-through relationships | Define the counting unit and test translation, nesting and Unsafe peptide containment. |
| GTEx v11 | Individual site TPM, donor and library identifiers, tissue medians and QC annotations | For each non-reproductive site take the linear donor P95, after the within-donor median of replicate libraries. Individual reproductive samples can independently support RA. |
| HPA bulk tissue and brain, including GTEx redistributed through HPA | Individual tissue biopsies, tissue and brain-region means, consensus nTPM and anatomical mappings | Native TPM contributes to the healthy reference, or individual reproductive biopsies to RA. Use within-tissue biopsy P95 when available. Consensus nTPM has its own normal gate. Redistributed GTEx measurements overlap GTEx. |
| HPA sorted blood, DICE and Monaco | Sorted immune populations, donor TPM, ranked values or published means, plus HPA immune nTPM | HPA donor P95 and DICE single-locus P95 enter the native comparison. DICE multi-locus groups use a conservative sum of member maxima because donor pairing is unavailable. Monaco retains the published mean. The immune screen also uses these panels. |
| HPA single-cell clusters and normal-tissue staining | Cell-type nTPM, staining level, antibody identity and reliability | Screen normal cell populations and supported protein findings. Reproductive staining and cell-cluster nTPM can support interpretation. Neither supplies a native-TPM tenfold witness. |
| CELLxGENE Census, including brain and eye populations, and the Siletti brain evidence | Versioned study files, donors, cell types, molecule counts and detection summaries | Grade sparse normal-cell expression by penetration, study replication, locus attribution and independent confirmation. |
| Collin cornea and conjunctiva | Cell counts and annotations with donor, age and disease metadata | Adult corneal observations have a dedicated normal gate in counts per 10,000. Developmental and disease material retains its own scope. |
| Specialized anatomy studies | Durante, Olender, Verbeurgt and MSSM olfactory data, fetal dental epithelium, developmental skeletal atlas and adult osteoblast markers | Supply tissue-localization evidence in underrepresented cell types. Age, assay scale and marker-list scope determine what each result can establish. A marker list cannot establish a tissue-wide negative result. |
| Petropoulos and Tyser embryos, Molè hypoblast | Per-embryo lineage TPM in preimplantation and Carnegie-stage-7 datasets, with day-9/day-11 UMI counts and staining in the hypoblast supplement | Supported native TPM can qualify RA. Hypoblast counts and protein localization document developmental expression on their own scales. |
| Tumor cohorts | Per-sample matrices, sample IDs, sharing across cohorts, source scale, QC and patient identifiers where available | **130 assessed cohort rows** in the [coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv). Only comparable QC-passing observations supply CA witnesses. Germ-cell tumors are eligible. Coordinated reproductive RNA is advisory. The [witness table](audits/target-selection-20261110/rebuild/tumor-witnesses.csv.gz) is a selected subset of the matrices. |
| Published gene lists, proteomics studies, clinical records and specific assays | Source-reported genes, locus mappings, allowed evidence uses, normal findings and supersession history, protein/epitope evidence, RPS4 primers and tissue panel, CT83 tumor RNAscope | Keep nominations, localized expression, sequence-specific evidence and normal counterevidence distinct. A curated normal finding can exclude a protein. Literature controls test the rules. A positive nomination alone cannot confer membership. |

| Inventory role | Versioned files |
| --- | --- |
| assay observations | 4 |
| cell atlas | 4 |
| comparison or sample metadata | 9 |
| curated data | 8 |
| curation input | 124 |
| native assay or supplement | 42 |
| tumor expression | 132 |

[Download API and atlas version inventory](antigen-data.md) · [Native RNA panels, measured sites and denominators](audits/target-selection-20261110/rebuild/normal-rna/index.md) · [Literature adjudication](audits/target-selection-20261110/rebuild/literature.md)

![Normal and reproductive tissue coverage by source](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage.png)

## The decision funnel

**1. Screen healthy tissues.** Bulk, cell-cluster and sorted-cell RNA have tier- and assay-specific ceilings. In TPM and nTPM the base ceilings are 4 for Vital, 8 for Protected, 16 for Manageable. Matched negative IHC can raise the applicable ceiling by a factor of two, four or eight depending on reliability. The registry also specifies limited clinical adjustments. An above-ceiling Vital or Protected observation can establish Unsafe. The Manageable Unsafe threshold is above 128 on those scales. Corneal counts have their own scaling. Reliable positive normal protein, reviewed normal findings and adult corneal RNA have their own gates. Single-cell and single-nucleus counts retain their native units. For a Vital population the donor signal requires at least 10 cells and 1% of the population with at least two molecules. A single-study signal needs 5% for the frequent-signal rule. Protected populations double those frequencies. Replication and locus attribution distinguish established evidence (Unsafe) from a consequential unresolved signal (Other:Review). Independent confirmation can establish a signal. [Thresholds, replication and assay rules](target-selection.md)

**2. Check eligibility and coverage.** Missing required native-TPM coverage across 35 healthy contexts, unresolved evidence or an unqualified translation can produce Other:Review. Unsupported nested sequences, read-through loci or ≥95% containment of the protein's 9-mers in Unsafe proteins can produce Other:Sequence. The catalog also checks a target maximum above background and a Disposable ceiling of 32 independently on TPM and nTPM. An unmeasured established home tissue blocks CA. RA can still qualify independently. GFY illustrates the distinction. Its adult human olfactory-neuron RNA establishes a normal home outside the native-TPM panel. Its tumor peak remains recorded. The Other:Unscreened tissue call reflects the missing like-for-like comparison. Developmental tooth expression for AMELX retains its fetal scope. [Adult olfactory study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) · [Cell and isoform findings](curation.md#cell-and-isoform-evidence)

GFY, also called *Goofy*, is a Golgi-associated regulator of olfactory signaling. In mice, loss of Gfy altered the localization of adenylyl cyclase III (ADCY3), shortened olfactory cilia and reduced odor sensitivity. Adult human biopsy data detected GFY RNA in 94.1% of sampled mature and 70.3% of immature olfactory neurons. A separate HEK293T-cell study found that coexpressing Gfy with other trafficking factors increased surface expression of tested odorant receptors. The human biopsy establishes the cell lineage, while the mouse knockout directly tests the effect of gene loss on smell. [Mouse knockout study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6619734/) · [Receptor trafficking experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC12677133/) · [Adult human olfactory atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/)

The adult whole thymus is Disposable under the catalog’s loss allowance. Defined immune cell populations retain their own injury tiers. A matched observational study associated adult thymectomy with increased mortality and cancer incidence. [Adult thymectomy study](https://pubmed.ncbi.nlm.nih.gov/37530823/) Medullary thymic epithelial cells express tissue-restricted antigens during central tolerance. Bulk thymus is retained in the healthy ledger under the `thymic` scope and excluded from both enrichment denominators. Independent normal-expression and protein screens still apply. If removing a former Unsafe call exposes a gap in the donor-resolved cell panel, the group stays in `Other:Review`. [Human thymic epithelial-cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2211762/)

**3. Evaluate cancer and reproductive expression independently.** Let **N** be the highest applicable non-thymic, non-reproductive normal reference and **D** the non-thymic Disposable maximum in native TPM. Each route requires **≥10 × max(N, D, 1 TPM)**. One passing tumor sample can qualify CA. One measured reproductive sample or supported developmental observation can qualify RA. A well-sampled embryo lineage needs ≥20 cells, ≥5 cells at ≥1 TPM and ≥10% positive cells. A smaller lineage needs ≥5 cells, ≥50% positive cells and qualifying observations in two independent embryos of the same age and lineage. P95 reduces the influence of isolated healthy extremes. It cannot diagnose contamination or establish absence from a rare population.

The following waterfall assigns each excluded record once. Other’s labels express the primary disposition and additional failures remain on the record. The pipeline compared native expression for **921 groups**, including counterfactual measurements for held proteins. **841** passed the shared eligibility checks. These two populations answer different questions.

![Attrition across the entire protein catalog](audits/target-selection-20261110/rebuild/summary/funnel-all.png)

![Attrition among published CTA or oncofetal nominations](audits/target-selection-20261110/rebuild/summary/funnel-published.png)

The second figure starts with the published nomination subset. Each figure attributes a protein to its first failing gate.

| First attributed exclusion | Removed | Remaining |
| --- | --- | --- |
| All assessed proteins | 0 | 19,893 |
| Bulk tissue RNA above the tissue ceiling | 17,507 | 2,386 |
| Cell-type cluster RNA (HPA single cell) | 436 | 1,950 |
| Sorted immune-cell RNA | 96 | 1,854 |
| Rare-cell RNA established by donor | 107 | 1,747 |
| Normal-tissue protein (IHC) | 50 | 1,697 |
| Curated or published normal expression | 49 | 1,648 |
| Adult cornea RNA | 5 | 1,643 |
| Coverage, identity, translation or unresolved evidence | 734 | 909 |
| Nested, Unsafe-contained or read-through sequence | 48 | 861 |
| No qualifying expression or bounds fail | 248 | 613 |
| Home tissue no panel resolves | 88 | 525 |

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

![Overlap of the CA and RA memberships](audits/target-selection-20261110/rebuild/summary/venn.png)

## Where the studied CTAs land

“Most studied” here means most **distinct curated papers nominating the protein as a CTA**, with ties ordered by protein label and entity ID. This is a reproducible registry measure. It cannot establish the ranking across all publications. The table shows the first ten plus selected familiar reference proteins. The [studied-protein table](antigen-whitepaper/studied-proteins.csv) covers **2238 groups** with a CTA paper, familiar-panel entry or curated clinical record. It retains source IDs, clinical evidence scope, tissue sources, peak tumor sample and separate known-patient counts. A family-level trial entry retains family-level scope.

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

N and the reproductive peak can come from different assays or donors. The downloadable table retains peak sources and contexts as recorded on the catalog row. At most ten tied normal contexts are shown there. All observations are available in the full comparison table. **NR** means the shared screen excluded the group before the native-TPM comparison ran (MAGEA3, MAGEA12, MAGEB2, CTCFL). The three cells on each such row refer to one skipped comparison. A tumor library or testis biopsy can still contain the transcript. **—** means a value is unavailable. The exclusion evidence and its native unit appear below.

### Studied proteins held or excluded

Of **58** distinct groups in the familiar or clinical panels, **30** are Other or Unsafe. Other:Review retains a consequential uncertainty. Unsafe records have an accepted normal-expression or clinical exclusion under this policy. Neither establishes the injury caused by a particular therapeutic product.

Selected examples below include additional paper-nominated proteins with consequential normal-cell findings. The [full study table](antigen-whitepaper/studied-proteins.csv) gives the disposition for every group in both panels. Filter on `familiar_names` or `clinical_scope` to reproduce this subset.

| Studied protein | Category | Recorded evidence | Matched assay where available |
| --- | --- | --- | --- |
| MAGEA3 | Unsafe | HPA tissue clusters | colon, undifferentiated cells · 53.1 nTPM (ceiling 4) |
| MAGEA12 | Unsafe | brain, pericyte (detection confirmed by a published gene-specific assay) | — |
| CTCFL | Unsafe | Normal bulk | — |
| CTAG2 | Other:Review | Normal-expression evidence blocks admission — lung, pulmonary alveolar type 2 cell (replicated normal-cell signal with unresolved locus attribution) | — |
| DDX43 | Unsafe | HPA normal protein  /  HPA tissue clusters  /  Monaco immune  /  DICE immune  /  eye, endothelial cell of venule (replicated signal at a locus without mapping risk) | — |
| SPA17 | Unsafe | HPA normal protein  /  Normal bulk  /  HPA tissue clusters  /  Adult cornea  /  lung, multiciliated columnar cell of tracheobronchial tree (replicated signal at a locus without mapping risk) | — |
| PLAC1 | Unsafe | brain, L2/3-6 intratelencephalic projecting glutamatergic neuron (replicated signal at a locus without mapping risk) | — |
| CLDN6 | Unsafe | Normal bulk  /  HPA tissue clusters | — |
| DAZL | Other:Review | Normal-expression evidence blocks admission — eye, astrocyte (frequent single-study normal-cell signal with unresolved locus attribution) | — |
| POTEC | Other:Review | Normal-expression evidence blocks admission — eye, corneal endothelial cell (replicated normal-cell signal with unresolved locus attribution) | — |
| CAGE1 | Other:Review | Normal-expression evidence blocks admission — brain, dopaminergic neuron (replicated normal-cell signal with unresolved locus attribution) | — |

### Why CTAG2 (LAGE-1) is unresolved

Two adult organs meet the donor-replication rule. The table shows each organ’s strongest qualifying donor/population. Counts refer to atlas collections and donor identifiers. The grader links shared Allen donor IDs across collections. Both lung collections integrate primary datasets. Cross-collection independence requires a donor and specimen crosswalk. [Human Lung Cell Atlas](https://www.nature.com/articles/s41591-023-02327-2) · [LungMAP CellRef](https://doi.org/10.1038/s41467-023-40173-5) In lung, the counted signals agree on alveolar type 2 cells. In brain, H18.30.002 occurs in both collections. The three source/donor records resolve to two named donors. The strongest layer 5 observation occurs in H18.30.002 alone. The broader intratelencephalic glutamatergic population crosses the screen in both H18.30.002 and H20.33.014. [Per-donor counts](audits/normal-cell-atlas-20261110/rare-cell/donor-populations.csv.gz) The eye observations did not meet the replication or frequency criteria.

| Organ | Strongest population | Atlas collections | Donor IDs with signal (linked where known) | Cells with ≥2 molecules in top donor | Grading |
| --- | --- | --- | --- | --- | --- |
| lung | pulmonary alveolar type 2 cell | 2 | 3 | 166/3,865 (4.29%) | unresolved |
| brain | L5 extratelencephalic projecting glutamatergic cortical neuron | 2 | 2 | 70/2,201 (3.18%) | unresolved |
| eye | ciliary muscle cell | 2 | 2 | 21/1,049 (2.00%) | dismissed |

**The open question is which locus produced the RNA.** 50.9% of the distinct CTAG2 transcript 31-mers also occur at other loci, with CTAG1A and CTAG1B identified as sequence partners. That fraction describes the reference sequences. It does not estimate how many observed reads were misassigned. The matched bulk check supplied no independent confirmation in lung or brain. A locus-specific observation in the implicated normal cell populations would distinguish genuine sparse CTAG2 expression from counts attributable to a paralog. [Donor findings](audits/normal-cell-atlas-20261110/rare-cell/findings.csv) · [Mapping audit](audits/normal-cell-atlas-20261110/rare-cell/gene-mapping-risk.csv.gz)

An independent qPCR study reported one low LAGE-1 signal in a panel of 15 normal-cell samples. The panel included pulmonary fibroblasts and astrocytes. Those measurements do not resolve the alveolar type 2 or cortical neuronal observations above. [Normal-cell qPCR panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC3624013/)

**Peripheral enrichment.** The highest peripheral reference is urinary bladder at 6.92 TPM, giving a cutoff of 69.25 TPM. The supported testis peak of 101.92 TPM reaches 14.72×. The ovarian tumor peak is 604.33 TPM, and 375 known patients meet the cancer cutoff across cohorts. CTAG2 clears both fold requirements. The uncertain origin of the lung and brain RNA keeps the catalog classification at Other:Review. Their cell fractions describe the strongest sampled donor. They cannot quantify clinical risk or establish presentation of a specific peptide. A 25-patient myeloma study of a receptor recognizing a peptide shared by NY-ESO-1 and LAGE-1 reported no fatal serious adverse events. That observation applies to the studied receptor and population. [Myeloma trial follow-up](https://pubmed.ncbi.nlm.nih.gov/31289029/) The separate heart review found discordant RNA and cardiomyocyte staining. [Heart assay review](audits/ctag2-heart-evidence.md)

## Which proteins reach the highest tumor expression?

Each table starts again at rank 1 and shows up to ten protein groups in that category, ranked by **maximum observed QC-passing tumor TPM**. Ties are ordered by protein label and entity ID. Across the candidate union, **106 proteins** have peaks strictly between 136 and 1,300 TPM. A protein contributes one peak even if many samples fall within this interval. In the [retained tumor witnesses](audits/target-selection-20261110/rebuild/tumor-witnesses.csv.gz), the same interval contains **4,488 candidate observations across 2,750 distinct source/sample IDs**. These counts cover the exported witnesses. The full matrices also contain observations below each protein’s witness threshold. Known patients at tenfold are counted once per identified person across cohorts.

The CA-only and RA-only tables exclude CTAs. The Other table covers only groups that reached the quantitative comparison. A strong tumor peak can coexist with a coverage or attribution hold. RA-only groups can also exceed the cancer fold when their recorded normal home lacks coverage. The downloadable ranks include inclusive CA and RA.

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

[All candidates, inclusive CA/RA/CTA ranks, sample and cohort identifiers, folds and patient support](antigen-whitepaper/tumor-ranks.csv). A peak measures abundance in a single library. The witness counts describe how often enrichment was observed in the available cohorts. Neither quantity measures the fraction of malignant cells expressing the protein.

## Reproduction and limits

Run `python scripts/render_antigen_whitepaper.py` from a source checkout, or use `--check` to compare all four artifacts byte for byte. The generator reads the versioned catalog, evidence registries and committed audit tables. `python -m oncoref.curation check` verifies the curation receipt. Rebuilding the **underlying catalog** additionally requires its pinned upstream inputs. The open scientific questions include tumor-source comparability, read-level attribution of rare-cell findings and the hair-follicle normal-protein replay. [Build instructions](curation.md) · [Searchable protein evidence](audits/target-selection-20261110/rebuild/report.html)
