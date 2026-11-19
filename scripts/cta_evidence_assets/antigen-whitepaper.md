# How OncoRef classifies human antigen candidates

**Research catalog · `$policy` · default biological sets**

Of **$loci human coding loci**, OncoRef tracks **$n protein-group records**. **$ca** qualify as Cancer Antigens (CA) and **$ra** as Reproductive Antigens (RA). The **$cta** that qualify independently on both sides are Cancer–Testis Antigens (CTA) in this catalog's broad reproductive/developmental sense. The union contains **$union** candidates. A catalog call describes expression and screening evidence. Antigen presentation, isoform production, tumor-cell localization and tolerance of an intervention each require evidence of their own.

## What counts as one protein?

OncoRef maps source gene identifiers to human Ensembl 115 loci. MANE Select on the primary locus is the first choice for the reference translation, followed by Ensembl canonical. Qualification checks the coding biotype, complete start and stop, an uninterrupted protein sequence and unresolved conflicting evidence. Loci with the same representative **amino-acid sequence** share one protein record. Different cDNAs can encode that same sequence. **$missing_seq records have no representative sequence hash** and remain visible in the catalog. Missing translation can block admission even when RNA is measured. [Reference selection and transcript-level evidence](audits/cta-transcript-selection-20261108/index.md)

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

The [download inventory](antigen-whitepaper/sources.csv) gives the selector, version, description, source URL, local or upstream asset name, and SHA-256 for **$files registered files**. Several files can come from the same study, and some are derived observations from upstream files. The [source-by-tissue table](audits/target-selection-20261110/rebuild/normal-rna/source-tissues.csv) lists **$tissue_rows source/context rows**, including tissue, finer anatomical site or cell type, and analytical role.

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
| Tumor cohorts | Per-sample matrices, sample IDs, sharing across cohorts, source scale, QC and patient identifiers where available | **$cohorts assessed cohort rows** in the [coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv). Only comparable QC-passing observations supply CA witnesses. Germ-cell tumors are eligible. Coordinated reproductive RNA is advisory. The [witness table](audits/target-selection-20261110/rebuild/tumor-witnesses.csv.gz) is a selected subset of the matrices. |
| Published gene lists, proteomics studies, clinical records and specific assays | Source-reported genes, locus mappings, allowed evidence uses, normal findings and supersession history, protein/epitope evidence, RPS4 primers and tissue panel, CT83 tumor RNAscope | Keep nominations, localized expression, sequence-specific evidence and normal counterevidence distinct. A curated normal finding can exclude a protein. Literature controls test the rules. A positive nomination alone cannot confer membership. |

$sources

[Download API and atlas version inventory](antigen-data.md) · [Native RNA panels, measured sites and denominators](audits/target-selection-20261110/rebuild/normal-rna/index.md) · [Literature adjudication](audits/target-selection-20261110/rebuild/literature.md)

![Normal and reproductive tissue coverage by source](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage.png)

## The decision funnel

**1. Screen healthy tissues.** Bulk, cell-cluster and sorted-cell RNA have tier- and assay-specific ceilings. In TPM and nTPM the base ceilings are $rna_bases. Matched negative IHC can raise the applicable ceiling by a factor of two, four or eight depending on reliability. The registry also specifies limited clinical adjustments. An above-ceiling Vital or Protected observation can establish Unsafe. The Manageable Unsafe threshold is above $manageable_ceiling on those scales. Corneal counts have their own scaling. Reliable positive normal protein, reviewed normal findings and adult corneal RNA have their own gates. Single-cell and single-nucleus counts retain their native units. For a Vital population the donor signal requires at least $cells cells and $pct% of the population with at least two molecules. A single-study signal needs $frequent% for the frequent-signal rule. Protected populations double those frequencies. Replication and locus attribution distinguish established evidence (Unsafe) from a consequential unresolved signal (Other:Review). Independent confirmation can establish a signal. [Thresholds, replication and assay rules](target-selection.md)

**2. Check eligibility and coverage.** Missing required native-TPM coverage across $required healthy contexts, unresolved evidence or an unqualified translation can produce Other:Review. Unsupported nested sequences, read-through loci or ≥95% containment of the protein's 9-mers in Unsafe proteins can produce Other:Sequence. The catalog also checks a target maximum above background and a Disposable ceiling of $cap independently on TPM and nTPM. An unmeasured established home tissue blocks CA. RA can still qualify independently. GFY illustrates the distinction. Its adult human olfactory-neuron RNA establishes a normal home outside the native-TPM panel. Its tumor peak remains recorded. The Other:Unscreened tissue call reflects the missing like-for-like comparison. Developmental tooth expression for AMELX retains its fetal scope. [Adult olfactory study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) · [Cell and isoform findings](curation.md#cell-and-isoform-evidence)

**3. Evaluate cancer and reproductive expression independently.** Let **N** be the highest applicable non-reproductive normal reference and **D** the Disposable maximum in native TPM. Each route requires **≥$fold × max(N, D, $floor TPM)**. One passing tumor sample can qualify CA. One measured reproductive sample or supported developmental observation can qualify RA. A well-sampled embryo lineage needs ≥20 cells, ≥5 cells at ≥1 TPM and ≥10% positive cells. A smaller lineage needs ≥5 cells, ≥50% positive cells and qualifying observations in two independent embryos of the same age and lineage. P95 reduces the influence of isolated healthy extremes. It cannot diagnose contamination or establish absence from a rare population.

The following waterfall assigns each excluded record once. Other’s labels express the primary disposition and additional failures remain on the record. The pipeline compared native expression for **$assessed groups**, including counterfactual measurements for held proteins. **$shared** passed the shared eligibility checks. These two populations answer different questions.

$funnel

$partition

![Overlap of the CA and RA memberships](audits/target-selection-20261110/rebuild/summary/venn.png)

## Where the studied CTAs land

“Most studied” here means most **distinct curated papers nominating the protein as a CTA**, with ties ordered by protein label and entity ID. This is a reproducible registry measure. It cannot establish the ranking across all publications. The table shows the first ten plus selected familiar reference proteins. The [studied-protein table](antigen-whitepaper/studied-proteins.csv) covers **$studied groups** with a CTA paper, familiar-panel entry or curated clinical record. It retains source IDs, clinical evidence scope, tissue sources, peak tumor sample and separate known-patient counts. A family-level trial entry retains family-level scope.

$publications

N and the reproductive peak can come from different assays or donors. The downloadable table retains peak sources and contexts as recorded on the catalog row. At most ten tied normal contexts are shown there. All observations are available in the full comparison table. **—** means the catalog lacks that comparison. Early Unsafe exclusions commonly have no harmonized tumor peak. The finding that triggered exclusion can still be available in another unit.

### Studied proteins held or excluded

Of **$anchors** distinct groups in the familiar or clinical panels, **$held** are Other or Unsafe. Other:Review retains a consequential uncertainty. Unsafe records have an accepted normal-expression or clinical exclusion under this policy. Neither establishes the injury caused by a particular therapeutic product.

Selected examples below include additional paper-nominated proteins with consequential normal-cell findings. The [full study table](antigen-whitepaper/studied-proteins.csv) gives the disposition for every group in both panels. Filter on `familiar_names` or `clinical_scope` to reproduce this subset.

$excluded

## Which proteins reach the highest tumor expression?

Rank by the **maximum observed QC-passing tumor TPM**, with ties ordered by protein label and entity ID. The first four in inclusive **CA** are **$top_ca**. The first four in inclusive **RA** are **$top_ra**. The table splits the memberships so that high-expressing CA-only and RA-only groups remain visible. An RA-only peak can miss cancer enrichment or another CA-specific requirement. Known patients are counted at the tenfold threshold across cohorts, once per identified person. Anonymous samples have a separate count in the download.

$peaks

[All candidates, inclusive CA/RA/CTA ranks, sample and cohort identifiers, folds and patient support](antigen-whitepaper/tumor-ranks.csv). A peak measures abundance in a single library. The witness counts describe how often enrichment was observed in the available cohorts. Neither quantity measures the fraction of malignant cells expressing the protein.

## Reproduction and limits

Run `python scripts/render_antigen_whitepaper.py` from a source checkout, or use `--check` to compare all four artifacts byte for byte. The generator reads the versioned catalog, evidence registries and committed audit tables. `python -m oncoref.curation check` verifies the curation receipt. Rebuilding the **underlying catalog** additionally requires its pinned upstream inputs. The open scientific questions include tumor-source comparability, read-level attribution of rare-cell findings and the hair-follicle normal-protein replay. [Build instructions](curation.md) · [Searchable protein evidence](audits/target-selection-20261110/rebuild/report.html)
