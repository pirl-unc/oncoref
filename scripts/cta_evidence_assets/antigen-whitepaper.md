# Cancer-selective and cancer–testis antigens in OncoRef

*Policy `$policy` · [Searchable catalog](audits/target-selection-20261110/rebuild/report.html) · [Tissue injury hierarchy](tissue-hierarchy.md) · [Data & downloads](antigen-data.md)*

Targeted cancer immunotherapies need proteins that are active in tumors and absent from healthy adult organs. Published cancer–testis antigen (CTA) lists often mix historical microarray nominations, incompatible expression scales, and patient-specific clinical assumptions into one table.

OncoRef screens the human coding proteome—**$loci** Ensembl 115 coding loci collapsed into **$n** identical-sequence protein groups—with one rule set that separates normal-tissue safety, biological enrichment, and patient-specific organ protection.

## Categories

Every protein group lands in one of five mutually exclusive categories:

- <span class="cat cat-unsafe">Unsafe</span> (`$unsafe_n`): confirmed expression in vital, protected, or high-expression manageable adult somatic tissues. Normal-tissue expression always takes precedence over tumor expression.
- <span class="cat cat-other">Other</span> (`$other_n`): held out of the candidate sets because required normal coverage is incomplete, a strong single-cell signal has unresolved paralog attribution (`Other:Review`), the sequence overlaps <span class="cat cat-unsafe">Unsafe</span> proteins (`Other:Sequence`), an adult normal home tissue lacks bulk panel coverage (`Other:Unscreened tissue`), or tumor and reproductive expression fall below the enrichment bar (`Other:Expression`).
- <span class="cat cat-ca">Cancer Antigens (CA)</span> (`$ca`): passes the shared safety and eligibility screens and reaches **≥$fold×** the non-thymic healthy adult baseline in at least one QC-passing tumor sample.
- <span class="cat cat-ra">Reproductive Antigens (RA)</span> (`$ra`): passes the same shared screen and reaches **≥$fold×** the identical healthy baseline in an adult reproductive tissue sample or supported early-embryo lineage.
- <span class="cat cat-cta">Cancer–Testis Antigens (CTA)</span> (`$cta`): the intersection `CA ∩ RA`, clearing both cancer and reproductive/developmental enrichment while staying silent across healthy somatic organs. The union `CA ∪ RA` contains **$union** candidate targets (`$cta` <span class="cat cat-cta">CTA</span>, `$ca_only` <span class="cat cat-ca">CA only</span>, `$ra_only` <span class="cat cat-ra">RA only</span>).

![Antigen classification overview](audits/target-selection-20261110/rebuild/summary/classification-overview.png)

*Figure 1.* Proteome partition across the five categories. Biological membership (<span class="cat cat-ca">CA</span>, <span class="cat cat-ra">RA</span>, <span class="cat cat-cta">CTA</span>) is independent of patient anatomy; patient-specific organ protection applies downstream without relabeling biological categories.

$headline

> These classifications are expression-screening calls, not clinical target certifications. Translation, antigen processing, peptide–MHC presentation, localization, and modality-specific safety require direct experimental validation.

## Why multiple normal atlases, single-cell RNA, and protein staining

No single assay catches every normal-tissue risk. Each modality in OncoRef covers a blind spot left by the others:

**Multiple bulk RNA atlases and sorted immune panels.** Bulk RNA removes most ubiquitously expressed proteins, but individual atlases differ in organ coverage and donor sampling. Pairing GTEx v11 donor distributions with HPA v23 biopsies and brain subregions prevents single-atlas gaps from reading as silence. Sorted blood subsets (HPA, DICE, Monaco) catch immune-lineage transcripts diluted in whole blood or spleen.

**Single-cell clusters and rare-cell UMI atlases.** Bulk RNA averages across millions of cells, so a gene expressed strongly in 1–5% of a vital population—such as alveolar type 2 cells, cardiomyocytes, or cortical neurons—can sit below bulk cutoffs. Single-cell atlases (HPA clusters, Collin cornea, donor-resolved CELLxGENE Census and Siletti brain) catch these localized risks. Requiring ≥2 UMI per cell replicated across independent studies and donors avoids discarding real targets over ambient RNA, while holding paralog-ambiguous genes for review instead of auto-excluding them.

**Normal-tissue protein staining (IHC).** Steady-state mRNA does not always track protein abundance: stable proteins can persist at modest RNA, while low-level transcript noise does not always yield detectable protein. Reliable positive HPA immunohistochemistry in a vital or protected organ independently excludes a target; consistent negative IHC across routine contexts of a tissue gives calibrated headroom above low bulk RNA.

**Curated normal literature and specialized anatomy.** Standard bulk panels omit small adult organs such as olfactory epithelium and taste buds. Reviewed normal-expression literature and targeted anatomical studies either exclude vital/protected homes or hold unmeasured adult homes as `Other:Unscreened tissue`.

![Unsafe proteins caught by complementary assays beyond bulk tissue RNA](audits/target-selection-20261110/rebuild/summary/normal-evidence-layers.png)

*Figure 2.* `Unsafe` proteins passing bulk tissue RNA ceilings that are caught by single-cell clusters, sorted immune panels, normal IHC, donor-resolved rare-cell UMI counts, curated literature, or adult cornea scRNA.

## 1. Screening attrition and the 10× enrichment rule

Applying the normal-tissue safety layers, eligibility checks, and enrichment rules narrows **$n** protein groups to **$shared** shared-eligible proteins (**$assessed** with native-TPM comparisons, including counterfactuals for held proteins) and **$union** retained candidates (`CA ∪ RA`).

![Proteome-wide attrition across screening gates](audits/target-selection-20261110/rebuild/summary/funnel-all.png)

*Figure 3.* Sequential attrition from **$n** protein groups to **$union** retained candidates, assigning each excluded protein once at its first failing gate. The broken x-axis above 3,000 keeps post-bulk filters readable. [Step counts (CSV)](audits/target-selection-20261110/rebuild/summary/funnel-all.csv).

Both enrichment routes apply the same threshold—**≥$fold×** above **H₀ = max(N, D, $floor TPM)**—where `N` is the highest non-thymic, non-reproductive healthy adult reference and `D` is the non-thymic Disposable maximum. Gold points above the dotted line clear the fold rule numerically but remain held in `Other:Review` or `Other:Unscreened tissue` (such as `CTAG2` and `GFY`).

![Cancer enrichment relative to the healthy adult reference](audits/target-selection-20261110/rebuild/summary/cancer-enrichment.png)

*Figure 4.* Cancer enrichment (<span class="cat cat-ca">CA</span>) across all **$assessed** protein groups with native-TPM comparisons.

![Reproductive and developmental enrichment relative to the healthy adult reference](audits/target-selection-20261110/rebuild/summary/reproductive-enrichment.png)

*Figure 5.* Reproductive and developmental enrichment (<span class="cat cat-ra">RA</span>) against the identical non-thymic healthy adult baseline.

Exclusive category counts across the entire proteome:

$partition

## 2. Published CTA, oncofetal, and placental lists

OncoRef maps **$lit_sources** curated papers across four antigen families: classical CTA / cancer–germline catalogs, oncofetal / DUX4 / early-embryo programs, placental and reproductive catalogs, and focused single-gene or tumor-vasculature studies. Literature citations are provenance metadata; a published nomination never overrides normal-tissue expression.

Retention tracks how strictly each paper screened healthy tissues. Focused modern CTA and DUX4 studies (`Seager et al. 2024`, `Chew et al. 2019`, `Carter et al. 2023`, `Loriot et al. 2025`) retain **34%–67% CTAs**, whereas broad germ-cell, microarray, ZGA, and cancer–placenta lists (`Bruggeman et al. 2018`, `Smith et al. 2023`, `Jamin et al. 2021`, `Bradley et al. 2020`) are **92%–100% `Unsafe`** under modern adult somatic screening.

![Mapped protein counts across curated literature sources](audits/target-selection-20261110/rebuild/summary/literature-sources-counts.png)

*Figure 6.* Mapped protein-group counts per curated source, grouped by antigen family and colored by OncoRef disposition. [Full table (CSV)](audits/target-selection-20261110/rebuild/summary/literature-sources.csv).

![Retention versus normal-tissue exclusion across curated literature sources](audits/target-selection-20261110/rebuild/summary/literature-sources-retention.png)

*Figure 7.* Share of mapped proteins retained as <span class="cat cat-cta">CTA</span>, <span class="cat cat-ca">CA only</span>, or <span class="cat cat-ra">RA only</span> versus held or excluded as <span class="cat cat-unsafe">Unsafe</span>.

Of **$published_noms** published CTA or oncofetal nominations, **$published_surviving** survive into the candidate union—including **$cta_published** of our **$cta** CTAs—while **$cta_new** CTAs are newly identified by the expression rules and **$published_outside** nominations fail normal-tissue, eligibility, or enrichment screens.

![Overlap between rule-derived sets and published nominations](audits/target-selection-20261110/rebuild/summary/venn.png)

*Figure 8.* Overlap between rule-derived Cancer Antigens (**$ca**), Reproductive Antigens (**$ra**), and published CTA/oncofetal nominations (**$published_noms**).

![Attrition of published CTA and oncofetal nominations](audits/target-selection-20261110/rebuild/summary/funnel-published.png)

*Figure 9.* Sequential attrition of **$published_noms** published CTA/oncofetal nominations through the same gates, leaving **$published_surviving** candidates. [Step counts (CSV)](audits/target-selection-20261110/rebuild/summary/funnel-published.csv).

### Most-nominated individual proteins

Across **$studied** studied protein groups (**$published_noms** with a CTA nomination paper plus **$study_extra** clinical or familiar reference entries), below are the top ten by distinct curated CTA papers plus five familiar clinical anchors (`NY-ESO-1`, `PRAME`, `CT83`, `CTAG2`, `CAGE1`):

$publications

`N` and the reproductive peak may come from different assays or donors ([full context table](antigen-whitepaper/studied-proteins.csv)). `NR` ("not run") marks proteins excluded as `Unsafe` in Stage 1 before native-TPM comparison ($nr_proteins); tumor or testis libraries may still contain those transcripts. A dash (`—`) marks an unavailable value.

### Studied proteins held or excluded

Of **$anchors** protein groups in the familiar or clinical reference panels, **$held** are held in `Other` or excluded as `Unsafe`. Neither label predicts clinical toxicity for a specific therapeutic product; each records the strongest normal-tissue or rare-cell finding observed under policy `$policy`:

$excluded

## 3. Peak tumor expression and held edge cases

Across the **$union** candidate union, **$union_ge_100** protein groups reach **≥100 TPM** in at least one QC-passing tumor sample, **$union_ge_1000** reach **≥1,000 TPM**, and **$union_patients_ge_10** have **≥10 known independent patients** crossing the $fold× cutoff. Median peak tumor expression is **$cta_median_tpm TPM** for `CTA` and **$ca_only_median_tpm TPM** for `CA only`.

Each table ranks up to ten protein groups in its exclusive category by maximum QC-passing tumor TPM ([all $union candidates](antigen-whitepaper/tumor-ranks.csv)). Because `Tumor fold` divides by `max(N, D, $floor TPM)`, fold equals `Tumor TPM` whenever both `N` and `D` are below $floor TPM (`LALBA`, `SCGB1D4`, `VCY/1B`). `Other` covers only groups reaching quantitative comparison; high peaks such as `TBATA`, `GFY`, and `CTAG2` reflect review or unscreened-home holds. Peak TPM measures single-library abundance, not the fraction of malignant cells expressing the protein.

$peaks

### Worked edge cases

**`CTAG2` (LAGE-1) → `$ctag2_classification`.** In peripheral bulk TPM, `CTAG2` clears both enrichment rules: its highest peripheral reference is **$ctag2_site** at **$ctag2_reference TPM** ($fold× cutoff **$ctag2_bar TPM**), testis reaches **$ctag2_repro TPM** (**$ctag2_fold×**), and ovarian cancer reaches **$ctag2_tumor TPM** across **$ctag2_patients** known patients at $fold×. In donor-resolved single-cell atlases, however, **lung** (alveolar type 2 cells) and **brain** (cortical glutamatergic neurons) both pass the multi-study donor-replication rule even after consortium and integrated-atlas donor deduplication:

$ctag2_cells

Because **$ctag2_shared** of distinct `CTAG2` transcript 31-mers are shared with paralogs (`CTAG1A`, `CTAG1B`), neither matched bulk RNA nor a 15-sample normal-cell qPCR panel (which assayed fibroblasts and astrocytes rather than alveolar type 2 or cortical neurons) resolves locus attribution. A 25-patient myeloma trial of a shared NY-ESO-1/LAGE-1 receptor reported no fatal serious adverse events (specific to that receptor and cohort), and a separate heart review found discordant RNA and cardiomyocyte staining. OncoRef holds `CTAG2` in **`$ctag2_classification`** pending locus-specific resolution ([donor findings](audits/normal-cell-atlas-20261110/rare-cell/findings.csv), [mapping audit](audits/normal-cell-atlas-20261110/rare-cell/gene-mapping-risk.csv.gz), [heart review](audits/ctag2-heart-evidence.md), [HLCA](https://www.nature.com/articles/s41591-023-02327-2), [LungMAP CellRef](https://doi.org/10.1038/s41467-023-40173-5), [qPCR panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC3624013/), [myeloma follow-up](https://pubmed.ncbi.nlm.nih.gov/31289029/)).

**`GFY` (*Goofy*) → `Other:Unscreened tissue`.** `GFY` regulates olfactory cilia signaling and reaches **1,575.56 TPM** in rhabdoid tumor (`RT`) over a 2.65 TPM bulk peripheral reference. Yet adult human olfactory-biopsy single-cell data detect `GFY` in **94.1%** of mature and **70.3%** of immature olfactory neurons ([adult olfactory atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/))—an organ absent from bulk native-TPM panels—and mouse *Gfy* knockout impairs `ADCY3` localization, olfactory cilia length, and odor sensitivity ([knockout study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6619734/), [trafficking study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12677133/)). Without a comparable olfactory TPM denominator, OncoRef keeps the tumor measurement visible and holds `GFY` in **`Other:Unscreened tissue`**.

## 4. Tissue injury tiers, thymus, and patient-specific filtering

OncoRef separates three tissue-policy questions: **anatomical scope** (healthy baseline vs. <span class="cat cat-ra">RA</span> witness vs. `thymic`), **adult injury tier** (clinical consequence of sustained organ or cell-type loss), and **patient-specific organ protection** (downstream anatomy and tolerance via `HostContext` / `--profile`). A tissue can witness <span class="cat cat-ra">RA</span> while retaining a protective rating—adult testis witnesses <span class="cat cat-ra">RA</span> while carrying a `Protected` tier, and breast biopsies are split into lactating versus resting mammary via a coordinated milk-program signature (`LALBA`, `CSN1S1`, `CSN2`, `CSN3`).

| Adult injury tier | Policy interpretation | Examples |
|---|---|---|
| **Vital** | Loss threatens life, essential repair, or major sensory/neurological function | Heart, neurons, corneal endothelium |
| **Protected** | Persistent loss is consequential and clinical compensation is poorly established | Salivary gland, whole skin, adult testis |
| **Manageable** | Defined loss has an established precedent for clinical compensation or replacement | Thyroid, mature B cells |
| **Disposable** | A narrow compartment carries a specific loss allowance | Cutaneous melanocytes, hair follicles, gallbladder, adult whole thymus |

Bulk measurements inherit the organ tier; permissive cell-type exceptions require resolved cell identity, and mixed populations keep the more protective tier ([organ and cell hierarchy](tissue-hierarchy.md), [tissue tier policy](tissue-tiers.md)).

### Adult thymus

Adult whole thymus is rated `Disposable` under the organ-loss hierarchy, while defined immune subsets retain their own injury tiers. Because medullary thymic epithelial cells (mTECs) express tissue-restricted antigens during central tolerance ([human mTEC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2211762/)), bulk thymus sits in the `thymic` scope and is excluded from both enrichment denominators (`N` and `D`). Independent normal-expression and protein screens still apply, alongside the observational association between adult thymectomy and higher mortality and cancer incidence ([thymectomy cohort](https://pubmed.ncbi.nlm.nih.gov/37530823/)).

### Biological membership versus patient-specific filtering

Biological <span class="cat cat-cta">CTA</span>, <span class="cat cat-ca">CA</span>, and <span class="cat cat-ra">RA</span> labels record where a protein is expressed, not what a particular patient can safely lose. Default biological sets leave reproductive tissues unprotected until you opt in. Downstream patient profiles (`female-patient`, `male-patient`, lactation, pregnancy, or custom overrides) mark each compartment as <span class="cat cat-protect">Protect</span>, <span class="cat cat-loss">Accept loss</span>, or <span class="cat cat-absent">Absent</span> ([rules and limits](host-selection.md)). Protecting testis narrows the **$cta** biological CTAs to **$testis_protected_cta** retained candidates; pregnancy leaves **$pregnant_retained** retained (**$pregnant_unresolved** unresolved due to incomplete developmental coverage, **$pregnant_excluded** excluded by a patient-anatomy rule).

![Patient anatomy and organ-loss choices across profiles](audits/target-selection-20261110/rebuild/summary/host-profiles-matrix.png)

*Figure 10.* Compartment choices (<span class="cat cat-protect">Protect</span>, <span class="cat cat-loss">Accept loss</span>, <span class="cat cat-absent">Absent</span>, or unprotected) across standard patient profiles.

![Patient-profile screening outcomes for biological CTAs](audits/target-selection-20261110/rebuild/summary/host-profiles-outcomes.png)

*Figure 11.* Downstream screening outcomes for all **$cta** biological CTAs under each patient profile. Filtering a candidate out for one patient narrows the target list and expands the peptide-subtraction background without relabeling its biological <span class="cat cat-cta">CTA</span> identity.

### Sex-restricted reproductive antigens and cross-sex CTAs

Because male and female reproductive organs are anatomically distinct, an antigen restricted to one sex's reproductive tract has no normal target organ in non-pregnant patients of the opposite sex. Oncology has no standard shorthand for this cross-sex distinction, so OncoRef partitions Reproductive Antigens (<span class="cat cat-ra">RA</span>) and their cancer intersection (<span class="cat cat-cta">CTA</span>) into four explicit sex-restriction classes under strict protection of all opposite-sex adult organs (including resting and lactating breast):

- **Male-Restricted (`MRA`, `$mra_n` RAs · `CTA-M`, `$cta_m_n` CTAs — female-safe):** expressed in male reproductive organs (`testis`, `epididymis`, `prostate`, `seminal vesicle`) and cleanly clearing every female adult reproductive and mammary organ (`XAGE1A/B`, `SPANXA1/2`, `MAGEA4`, `NY-ESO-1`, `VCY/1B`, plus male-accessory/prostate targets such as `PRR23D1/2`, `SPAG11A`, `POTEH`, `EDDM3B`). Relaxing resting-breast protection widens `CTA-M` to **$cta_m_sex_only_n** CTAs, bringing in high-peak prostate and seminal-vesicle genes (`TGM4`, `SEMG2`) that show resting-breast RNA in GTEx.
- **Embryo / Placenta-Only (`ERA`, `$era_n` RAs · `CTA-E`, `$cta_e_n` CTAs — adult-absent in both sexes):** cleanly clear both all adult male and all non-pregnant adult female organs (`TBXT`, `MBD3L2`, `SSU72L2`, `TRIM48`, `TRIM43B`).
- **Female-Restricted (`FRA`, `$fra_n` RAs · `CTA-F`, `$cta_f_n` CTAs — male-safe):** confined to female reproductive compartments (`ovary`, tract, `placenta`/`lactating_breast`) while cleanly clearing all male adult organs (`ERVV-1`, `ZP4`).
- **Shared Male & Female (`SRA`, `$sra_n` RAs · `CTA-MF`, `$cta_mf_n` CTAs):** expressed in adult reproductive or mammary compartments of both sexes (`PRAME`, `LALBA`, `SCGB1D4`).

The asymmetry is striking: **$cta_female_safe_n** of **$cta** CTAs (`CTA-M + CTA-E`) clear every adult female organ in non-pregnant women, whereas only **$cta_male_safe_n** (`CTA-F + CTA-E`) clear every adult male organ, because transcriptional scanning in adult testis expresses nearly the entire ovarian, uterine, and mammary repertoire.

![Sex-restricted Reproductive Antigens and cross-sex CTA subsets](audits/target-selection-20261110/rebuild/summary/sex-restricted-ra-cta.png)

*Figure 12.* Partition of all **$ra** Reproductive Antigens, **$cta** CTAs, and **$ra_only** RA-only proteins into `MRA`/`CTA-M`, `ERA`/`CTA-E`, `FRA`/`CTA-F`, and `SRA`/`CTA-MF`. [Download full sex-restriction table (CSV)](audits/target-selection-20261110/rebuild/summary/sex-restricted-antigens.csv).

![Male versus female reproductive-organ exposure across all CTAs](audits/target-selection-20261110/rebuild/summary/sex-restricted-landscape.png)

*Figure 13.* Maximum male versus female sex-organ RNA across all **$cta** CTAs. Points in the lower-right (`CTA-M`, blue) are anatomically absent in non-pregnant women; points in the lower-left (`CTA-E`, green) are absent in both adult sexes; points in the upper-left (`CTA-F`, plum) are anatomically absent in men.

## 5. Methods: protein grouping, sources, and thresholds

### Protein-group counting unit

Source gene identifiers across RNA, single-cell, and protein atlases map to human **Ensembl 115** loci. Each locus gets one qualified reference translation—preferring **MANE Select** on the primary assembly, then **Ensembl canonical**—checked for coding biotype, complete start/stop codons, uninterrupted ORF, and absence of translation conflicts ([transcript selection audit](audits/cta-transcript-selection-20261108/index.md)).

Paralogous loci encoding identical proteins (`SPANXA1/2`, `XAGE1A/B`, `VCY/1B`) collapse into one protein group by representative amino-acid sequence even when cDNAs differ. Wherever sample identities exist, member-locus expression is summed within each sample before comparing across cohorts. The **$missing_seq** records lacking a qualified representative sequence stay visible; missing translation can block admission even when RNA is measured ([sequence and peptide-subtraction API](antigen-sequences.md)).

### Evidence sources and native-TPM coverage

OncoRef pins **$files** SHA-256-verified files across **$tissue_rows** source/context rows and **$cohorts** assessed tumor cohorts. Incompatible assay units are never pooled; every dataset stays in its native scale and enters only the layer it supports:

| Evidence pillar | Primary resources | Native statistic and role |
|---|---|---|
| **Healthy baseline (`N`, `D`)** | GTEx v11, HPA v23 biopsies and brain, HPA/DICE/Monaco sorted immune | Within-site donor or biopsy P95 in native TPM sets the non-thymic healthy reference (`N`) and Disposable maximum (`D`). |
| **Normal-tissue safety screens** | HPA consensus and scRNA clusters, HPA normal IHC, CELLxGENE Census + Siletti, Collin cornea, specialized anatomy | Tiered RNA ceilings (TPM/nTPM/CP10K), reliable positive IHC, matched negative-IHC credit, and consortium-deduplicated rare-cell UMI screens. |
| **Tumor and reproductive/embryo RNA** | $cohorts QC-passing tumor cohorts, GTEx/HPA reproductive biopsies, Petropoulos/Tyser/Molè embryos | Individual QC-passing tumor libraries qualify `CA`; individual reproductive samples or supported embryo lineages qualify `RA`. |
| **Curated literature and clinical controls** | $lit_sources catalog-mapped papers (CTA, oncofetal/DUX4, placental, clinical trials, normal counterevidence) | Provenance metadata and curated normal-expression exclusions; positive nominations alone never grant admission. |

![Native-TPM coverage across somatic, immune, disposable, and thymic tissues](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage-somatic.png)

*Figure 14.* Native-TPM tissue coverage and analytical role across somatic organs, immune subsets, disposable tissues, and thymus. Colored cells mark qualifying coverage for at least one assessed protein group; gray marks no qualifying coverage.

![Native-TPM coverage across neural, ocular, reproductive, and developmental tissues](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage-neural-repro.png)

*Figure 15.* Native-TPM coverage across brain, eye, nerve, reproductive, and embryonic tissues. **RA witness** marks individual reproductive samples or supported embryo-lineage means capable of qualifying `RA`; **Reproductive context** marks summaries kept for context only.

<details markdown="1">
<summary>Exact quantitative gate thresholds and rare-cell provenance rules</summary>

* **Bulk, cluster, and immune RNA ceilings:** Base limits in TPM and nTPM are **$rna_bases**. Matched negative IHC across routine contexts of the tissue raises the ceiling by 2× (Approved), 4× (Supported), or 8× (Enhanced), alongside limited registry clinical adjustments. Exceeding a `Vital` or `Protected` ceiling establishes `Unsafe`; `Manageable` establishes `Unsafe` only above **$manageable_ceiling**.
* **Donor-resolved rare-cell RNA and consortium provenance:** Evaluated in native UMI counts per adult donor and named cell type. In a `Vital` population, a donor signal requires **≥$cells cells** and **≥$pct%** of cells with ≥2 molecules (replicated across independent study groups), or **≥$frequent%** in a single study; `Protected` populations double those frequency cutoffs. Shared consortium collections (Baylor/Chen eye, Sanes eye, Allen reference brain, Michigan kidney, Teichmann heart) and donors embedded in integrated atlases (HLCA, Baylor/Chen eye, SEA-AD) are resolved to their original study group and biological donor ID before counting replication. Replicated signals without locus mapping risk (or with independent bulk or gene-specific confirmation) establish `Unsafe`; unresolved locus attribution or single-study signals remain held in `Other:Review`.
* **Structural and home-tissue gates:** Unsupported nested sequences, read-through loci, or ≥95% 9-mer containment in `Unsafe` proteins produce `Other:Sequence`. A curated adult normal home outside the native-TPM panels (`GFY` in olfactory neurons) blocks `CA` as `Other:Unscreened tissue` while still allowing independent `RA` qualification; fetal-only homes (`AMELX` in developing tooth epithelium) keep their developmental scope.
* **Eligibility and Disposable cap:** Shared eligibility requires complete coverage across **$required** healthy contexts, target peak above every measured non-reproductive and Disposable value, and Disposable expression ≤ **$cap** TPM and ≤ **$cap** nTPM.
* **Embryo lineage support for RA:** Well-sampled lineages require ≥20 cells, ≥5 cells at ≥1 TPM, and ≥10% positive cells; smaller lineages require ≥5 cells, ≥50% positive cells, and replication in two independent embryos of the same stage and lineage. Within-site donor/biopsy **P95** dampens isolated healthy outliers but cannot diagnose carryover or prove absence from an unmeasured rare cell type.

</details>

<details markdown="1">
<summary>Source-family rules and versioned file inventory ($files files across $tissue_rows contexts)</summary>

| Source family | Evidence provided | Role in classification |
|---|---|---|
| **Ensembl 115, MANE and sequence registry** | Locus IDs, transcript flags, CDS qualification, sequence digests, read-through maps | Define protein-group counting units and test translation, nesting, and `Unsafe` 9-mer containment. |
| **GTEx v11** | Individual site TPM, donor/library IDs, site medians, QC annotations | Non-reproductive sites contribute the within-site linear donor P95 (after within-donor replicate median). Individual reproductive samples can independently witness `RA`. |
| **HPA bulk tissue and brain** (incl. GTEx via HPA) | Individual tissue biopsies, tissue/brain-region means, consensus nTPM, anatomical mappings | Native TPM enters the healthy reference (within-tissue biopsy P95 where available) or reproductive `RA` witnesses. Consensus nTPM runs its own normal gate. Redistributed GTEx overlaps GTEx v11. |
| **HPA sorted blood, DICE and Monaco** | Sorted immune populations, donor TPM, ranked values, published means, immune nTPM | HPA donor P95 and DICE single-locus P95 enter native comparison; DICE multi-locus groups use the conservative sum of member maxima when donor pairing is unavailable; Monaco retains the published mean. Also evaluated in the immune screen. |
| **HPA single-cell clusters and normal IHC** | Cell-type nTPM, staining grades, antibody IDs, reliability scores | Screen normal cell clusters and reliable positive protein findings (including hair-follicle cell-pattern rules); matched negative IHC can raise RNA ceilings. Reproductive staining/nTPM aids interpretation but never supplies a native-TPM $fold× witness. |
| **CELLxGENE Census and Siletti brain** | Versioned H5AD studies, donor IDs, cell types, molecule counts | Grade sparse normal-cell expression by donor penetration, independent study replication, locus mapping risk, and independent confirmation. |
| **Collin cornea and conjunctiva** | Adult, developmental, and disease ocular cell counts with donor metadata | Dedicated adult corneal normal gate in counts per 10,000 (CP10K); non-adult and disease material stays scoped separately. |
| **Specialized anatomy studies** | Olfactory epithelium (Durante, Olender, Verbeurgt, MSSM), fetal dental epithelium, developmental skeletal atlas, adult osteoblast markers | Document localized expression in tissues underrepresented in bulk atlases; age, scale, and marker-list scope constrain what each study can establish (marker lists never prove tissue-wide absence). |
| **Petropoulos, Tyser and Molè embryos** | Per-embryo lineage TPM (preimplantation, Carnegie stage 7) plus day-9/day-11 hypoblast UMI/staining | Supported native TPM can qualify `RA`; hypoblast UMI counts and staining document developmental expression on their native scales. |
| **Tumor cohorts** | Per-sample matrices, sample/patient IDs, cross-cohort sharing, scale and QC checks | **$cohorts assessed cohorts** in the [tumor coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv). Comparable QC-passing libraries supply `CA` witnesses (germ-cell tumors eligible; coordinated reproductive RNA advisory). |
| **Curated papers, proteomics and trials** | Nominated gene lists, normal counterevidence, supersession history, clinical/epitope records, CT83 RNAscope | Separate positive nominations from normal counterevidence and literature benchmarks. Curated normal findings can exclude a protein; a positive nomination alone never grants admission. |

$sources

</details>

## 6. Reproducing the catalog and limitations

Every external reference file and per-sample tumor matrix required for a clean rebuild is SHA-256-pinned and fetched through OncoRef's uniform CLI (no out-of-band downloads or machine-local paths):

```sh
# 1. Fetch and verify all pinned reference datasets (add --include-tumors for all 132 tumor cohorts)
python -m oncoref.curation fetch --include-tumors

# 2. Rebuild the full curation graph (catalog, host exposure, figures, reports, whitepaper)
python -m oncoref.curation build

# 3. Verify every tracked implementation file, catalog table, and whitepaper artifact byte-for-byte
python -m oncoref.curation check
python scripts/render_antigen_whitepaper.py --check
```

### Limitations

Cross-cohort tumor matrices come from distinct studies and processing pipelines; OncoRef enforces nonnegative comparable TPM and library QC rather than cross-study batch correction. In rare-cell atlases, count matrices cannot separate paralog cross-mapping from genuine sparse transcription without read-level realignment or locus-specific assays (as in `CTAG2`), and while known consortium collections and integrated-atlas donor overlaps are deduplicated, undocumented donor reuse across unrelated public studies cannot be ruled out. Single-sample RNA peaks do not measure intratumoral heterogeneity, translated protein abundance, peptide–MHC presentation, or modality-specific clinical tolerability.

### Reference and download index

| Resource | Contents |
|---|---|
| [Searchable protein catalog](audits/target-selection-20261110/rebuild/report.html) | Interactive browser for all **$n** protein groups, per-gate reasons, witnesses, and literature links |
| [Summary figures and distributions](audits/target-selection-20261110/rebuild/summary.md) | Threshold heatmaps, scale distributions, rare-cell plots, and stage breakdowns |
| [Classification rules](target-selection.md) · [Build guide](curation.md) | Machine-readable policy contract, CLI commands, and build instructions |
| [Tissue injury hierarchy](tissue-hierarchy.md) · [Tissue tiers](tissue-tiers.md) · [Patient-specific filtering](host-selection.md) | Cited adult organ/cell tiers, hair-follicle and thymus revisions, and patient-profile filtering |
| [Sequence API and peptide subtraction](antigen-sequences.md) · [Data download API](antigen-data.md) | Ensembl 115 reference translations, k-mer subtraction, and programmatic dataset retrieval |
| [Registered file inventory (`sources.csv`)](antigen-whitepaper/sources.csv) | All **$files** pinned source files, catalog selectors, and SHA-256 digests |
| [Studied CTA dispositions (`studied-proteins.csv`)](antigen-whitepaper/studied-proteins.csv) | Full measurements and exclusion reasons for all **$studied** published or clinical anchor proteins |
| [Candidate tumor rankings (`tumor-ranks.csv`)](antigen-whitepaper/tumor-ranks.csv) | Peak tumor TPM, folds, cohort/sample IDs, and patient counts for all **$union** candidates |
| [Curated literature sources (`literature-sources.csv`)](audits/target-selection-20261110/rebuild/summary/literature-sources.csv) | Disposition counts across all **$lit_sources** catalog-mapped papers grouped by antigen class |
| [Sex-restricted antigens (`sex-restricted-antigens.csv`)](audits/target-selection-20261110/rebuild/summary/sex-restricted-antigens.csv) | Male-restricted (`MRA`/`CTA-M`), female-restricted (`FRA`/`CTA-F`), embryo/placenta-only (`ERA`/`CTA-E`), and shared (`SRA`/`CTA-MF`) calls with organ TPMs |
| [Normal RNA coverage](audits/target-selection-20261110/rebuild/normal-rna/index.md) · [Tumor coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv) · [Literature adjudication](audits/target-selection-20261110/rebuild/literature.md) | Audit ledgers for healthy panels, tumor cohorts, and literature reviews |
