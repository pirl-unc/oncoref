# Cancer-selective and cancer–testis antigens in OncoRef

*[Interactive protein catalog](audits/target-selection-20261110/rebuild/report.html) · [Multi-species whitepaper](species-antigen-whitepaper.md) · [Tissue injury hierarchy](tissue-hierarchy.md) · [Data & downloads](antigen-data.md)*

<p class="section-nav"><strong>Jump to:</strong> <a href="#categories">Categories</a> · <a href="#why-bulk-rna-alone-is-insufficient">Why multi-assay safety</a> · <a href="#clinical-toxicity-gold-standard-potency-proven-tcr-t-and-tcr-bispecific-trials">Clinical TCR gold standard</a> · <a href="#decision-tree-attrition-and-the-10-enrichment-rule">Decision tree &amp; 10× rule</a> · <a href="#published-cta-oncofetal-and-placental-lists">Published lists</a> · <a href="#peak-tumor-expression-and-edge-cases">Top tumor targets</a> · <a href="#tissue-injury-tiers-patient-profiles-and-sex-restricted-ctas">Tissue tiers &amp; sex-restricted CTAs</a> · <a href="#methods-protein-grouping-sources-and-thresholds">Methods</a> · <a href="#comparative-multi-species-framework-human-dog-and-mouse">Multi-species</a> · <a href="#reproducing-the-catalog-and-downloads">Reproducibility &amp; downloads</a></p>

Targeted cancer immunotherapies need proteins that are active in tumors and absent from healthy adult organs. Published cancer–testis antigen (CTA) lists often mix historical microarray nominations, incompatible expression scales, and patient-specific clinical assumptions into one table.

OncoRef screens the human coding proteome—**$loci** Ensembl 115 coding loci collapsed into **$n** identical-sequence protein groups—with one rule set that separates normal-tissue safety, biological enrichment, and patient-specific organ protection.

## Categories

Every protein group lands in one of five mutually exclusive categories:

- <span class="cat cat-unsafe">Unsafe</span> (`$unsafe_n`): confirmed expression in vital, protected, or high-expression manageable adult somatic tissues. Normal-tissue expression always takes precedence over tumor expression.
- <span class="cat cat-other">Other</span> (`$other_n`): kept out of the candidate sets without being confirmed <span class="cat cat-unsafe">Unsafe</span>—either because required evidence is incomplete or ambiguous (missing normal-tissue coverage, unresolved paralog attribution in single-cell data `Other:Review`, or sequence overlap with <span class="cat cat-unsafe">Unsafe</span> proteins `Other:Sequence`), or because tumor and reproductive expression do not qualify (below the $fold× threshold `Other:Expression`, or an adult normal organ missing from bulk panels `Other:Unscreened tissue`).
- <span class="cat cat-ca">Cancer Antigens (CA)</span> (`$ca`): passes the safety and eligibility screens and reaches **≥$fold×** the non-thymic healthy adult baseline in at least one QC-passing tumor sample.
- <span class="cat cat-ra">Reproductive Antigens (RA)</span> (`$ra`): passes the same safety and eligibility screens and reaches **≥$fold×** the identical healthy baseline in an adult reproductive tissue sample or supported early-embryo lineage.
- <span class="cat cat-cta">Cancer–Testis Antigens (CTA)</span> (`$cta`): the intersection `CA ∩ RA`, clearing both cancer and reproductive/developmental enrichment while staying silent across healthy somatic organs. The union `CA ∪ RA` contains **$union** candidate targets (`$cta` <span class="cat cat-cta">CTA</span>, `$ca_only` <span class="cat cat-ca">CA only</span>, `$ra_only` <span class="cat cat-ra">RA only</span>).

![Antigen classification overview](audits/target-selection-20261110/rebuild/summary/classification-overview.png)

*Figure 1.* Proteome partition across the five categories. In `H₀ = max(Normal N, Disposable D, 1 TPM)`, `N` is the protein's highest native TPM across vital, protected, and manageable healthy adult somatic organs and sorted immune subsets, and `D` is its maximum across non-thymic Disposable organs (appendix, gallbladder, tonsil). Biological membership (<span class="cat cat-ca">CA</span>, <span class="cat cat-ra">RA</span>, <span class="cat cat-cta">CTA</span>) is independent of patient anatomy; patient-specific organ protection applies downstream without relabeling biological categories.

> These classifications are expression-screening calls, not clinical target certifications. Translation, antigen processing, peptide–MHC presentation, localization, and modality-specific safety require direct experimental validation.

## Why bulk RNA alone is insufficient

No single assay catches every normal-tissue risk. Each modality in OncoRef covers a blind spot left by the others:

**Multiple bulk RNA atlases and sorted immune panels.** Bulk RNA removes most ubiquitously expressed proteins, but individual atlases differ in organ coverage and donor sampling. Pairing GTEx v11 donor distributions with HPA v23 biopsies and brain subregions prevents single-atlas gaps from reading as silence. Sorted blood subsets (HPA, DICE, Monaco) catch immune-lineage transcripts diluted in whole blood or spleen.

**Single-cell clusters and rare-cell UMI atlases.** Bulk RNA averages across millions of cells, so a gene expressed strongly in 1–5% of a vital population—such as alveolar type 2 cells, cardiomyocytes, or cortical neurons—can sit below bulk cutoffs. Single-cell atlases (HPA clusters, Collin cornea, donor-resolved CELLxGENE Census and Siletti brain) catch these localized risks. Requiring ≥2 UMI per cell replicated across independent studies and donors avoids discarding real targets over ambient RNA, while placing paralog-ambiguous genes in `Other:Review` instead of auto-excluding them.

**Normal-tissue protein staining (IHC).** Steady-state mRNA does not always track protein abundance: stable proteins can persist at modest RNA, while low-level transcript noise does not always yield detectable protein. Reliable positive HPA immunohistochemistry in a vital or protected organ independently excludes a target; consistent negative IHC across routine contexts of a tissue gives calibrated headroom above low bulk RNA.

**Curated normal literature and specialized anatomy.** Standard bulk panels omit small adult organs such as olfactory epithelium and taste buds. Reviewed normal-expression literature and targeted anatomical studies either exclude vital/protected tissues or assign unmeasured adult tissues to `Other:Unscreened tissue`.

![Unsafe proteins caught beyond bulk tissue RNA](audits/target-selection-20261110/rebuild/summary/normal-evidence-layers.png)

*Figure 2.* `Unsafe` proteins passing bulk tissue RNA ceilings that are caught by single-cell clusters, sorted immune subsets, normal IHC, donor-replicated rare-cell UMI counts, curated literature, or adult cornea scRNA. Bulk tissue RNA alone flags **$bulk_unsafe_n** `Unsafe` proteins; these six assays add **$non_bulk_added_n** more (bars overlap where multiple non-bulk assays agree).

## Clinical toxicity gold standard: potency-proven TCR-T and TCR-bispecific trials

Healthy-tissue bulk TPM, single-cell clusters, normal immunohistochemistry, donor-resolved rare-cell UMI counts, and `1-mismatch` peptide-overlap screens are **mechanistic proxies**—not direct measurements of human T-cell toxicity. The empirical gold standard for intracellular target safety is **human clinical experience with peptide–MHC-directed T-cell receptors (`TCR-T` and soluble `ImmTAC` bispecifics) targeting endogenous self-encoded proteins in trials that achieved objective anti-tumor responses (`ORR > 0%`)**. Requiring confirmed clinical tumor shrinkage (or fulminant *in vivo* target-cell cytolysis) is essential: when a trial sees `0% ORR` (`IMA202` against `MAGEA1`, `0/16` responses), the absence of healthy-organ injury cannot prove the target tissue is antigen-negative, because TCR avidity, cell dose, or persistence may simply have been insufficient to kill antigen-bearing cells (`oncoref.tcr_trials.gold_standard`, `oncoref tcr-trials --table benchmark`).

Across **$tcr_gold_n** curated endogenous-protein clinical `TCR-T` / `ImmTAC` benchmark entries (**$tcr_potent_n** potency-proven regimens plus **$tcr_nonpotent_n** segregated `0% ORR` or unreadout controls; somatic neoantigens such as `KRAS G12D` and viral antigens such as `HPV16 E6/E7` are excluded because they lack a germline normal-proteome counterpart), OncoRef's `1.0×` gate ceiling cleanly separates clinically safe CTAs from lethal or dose-limiting specificity toxicities ([download full benchmark CSV](antigen-whitepaper/clinical-tcr-gold-standard.csv)):

- **`100%` sensitivity to severe human clinical specificity toxicities (`$tcr_toxic_caught_n/$tcr_toxic_n` flagged `Unsafe`, `0` false-safe misses):**
  - *On-target / off-tumor healthy-organ destruction (`4/4` `Unsafe`):* **`CEACAM5`** (`33%` ORR → `3/3` dose-limiting Grade 3 pan-colitis; `563×` colon enterocyte ceiling), **`MLANA` / MART-1** (`DMF5`, `30%` ORR → `70%` rash/vitiligo, `55%` anterior ocular uveitis, `50%` stria vascularis sensorineural hearing loss; `451×` ocular melanocyte `Vital` ceiling), **`PMEL` / gp100** (`19%` TCR-T ORR → `94%` rash, `25%` uveitis, `31%` hearing loss; soluble `tebentafusp` FDA-approved with `>80%` rash but minimal uveitis/ototoxicity across the blood–ocular/labyrinthine barrier; `889×` ocular melanocyte ceiling), and **`TYR` / tyrosinase** (`TIL1383I`, `33%` ORR → diffuse vitiligo; `75×` ocular melanocyte ceiling). Together, `MLANA`, `PMEL`, and `TYR` directly validate why OncoRef separates `Disposable` cutaneous melanocytes from `Vital` ocular and inner-ear pigment cells.
  - *Fatal paralog and structural off-target cross-reactivity (`2/2` `Unsafe`):* **`MAGEA3` (`HLA-A*02:01 / KVAELVHFL`, `56%` ORR → `2/9` fatal necrotizing leukoencephalopathy)** is missed in `MAGEA12` bulk brain RNA (`0.5 nTPM`, `0.13×`) but caught on primary `MAGEA3` colon `scRNA-seq` (`13.3×`), qPCR-corroborated rare-cell brain RNA on `MAGEA12` (`4.75%` ≥1-UMI, `0.40%` ≥2-UMI), `EPS8L2` (`46.1×`), and `#769` `1-mismatch` peptide overlap (`KVAELVHFL ↔ KMAELVHFL`). **`MAGEA3` (`HLA-A*01:01 / EVDPIGHLY` affinity-enhanced `a3a` TCR → `2/2` fatal cardiogenic shock within 4–5 days)** cross-recognized cardiac **`TTN`** (`ESDPIVAQY`, `53,894 nTPM` cardiomyocytes = `13,473×` ceiling; primary `MAGEA3` itself is also `Unsafe` at `13.3×`); because `4` of `9` residues differ, `a3a` also marks the boundary where sequence k-mers end and experimental X-scan / iPSC-cardiomyocyte screening is required for affinity-matured receptors.
- **`100%` recall on potency-proven clinically safe CTAs (`$tcr_safe_cta_kept_n/$tcr_safe_cta_n` retained as `CTA`):** **`NY-ESO-1`** (`58%` ORR in Robbins 2015, `30%–80%` across MRCLS and myeloma, `0.06×`), **`MAGEA4`** (FDA-approved `afami-cel`, `37%` Phase 2 ORR, `0.01×`), **`PRAME`** (`IMA203`, `52.5%` u/cORR [`28.9%` overall / `54.5%` Ph1b cORR] across `40` patients without nephrotoxicity or adrenal insufficiency, `0.75×`), **`CTAG2` (`LAGE-1a`)** (`lete-cel`, `30%` pooled / `40%` standard-LDR MRCLS and `80%` myeloma ORR without pulmonary or urothelial toxicity, `0.78×`), and **`MAGEA10`** (`ADP-A2M10`, `9.1%` ORR, `0.30×`). Notably, **two of the five clinically proven safe CTAs (`PRAME` and `CTAG2`) require OncoRef's multi-modal calibration**—matched negative-IHC RNA headroom (`PRAME` shifts from `3.0×` naive RNA to `0.75×` under negative IHC) and isoform-scoped vascular perfusion clearance (`CTAG2` `ENST00000369585` `#766`)—to avoid false exclusion.
- **Conservative gene-level precision frontier (`$tcr_conservative_n` tolerated non-CTA regimens flagged `Unsafe`):** **`WT1`** (`TTCR-C4`, `100%` post-HCT AML relapse-free survival without podocyte nephrotic syndrome; `9.9×` mesothelial RNA + podocyte IHC `High`), **`MAGEA3` (`HLA-DPB1*04:01` `6F9` CD4 TCR-T, `24%` ORR without CNS/cardiac toxicity; `13.3×` colon `scRNA`)**, and **`AFP`** (`ADP-A2AFP`, `9.5%` HCC ORR without autoimmune hepatitis; `98×` renal proximal tubule `scRNA`) were clinically tolerated for their specific `(peptide, HLA)` complexes—reflecting immunoproteasome-dependent differential processing (`WT1` `RMFPNAPYL`), MHC-II restriction, or sub-threshold parenchymal pMHC density. Gene-level screening intentionally remains `Unsafe` on these proteins because epitope-specific processing escape cannot be generalized across other peptides or surface-targeted modalities.
- **Segregated non-potent or unreadout controls (`$tcr_nonpotent_n` entries excluded from safety calibration):** **`MAGEA1`** (`IMA202`, `0/16 = 0% ORR`) and **`CT83` / KK-LC-1** (ongoing Phase 1 `NCT05035407` / `NCT05483491`) are retained as `CTA-M` solely on independent multi-assay somatic silence (`0.00×` and `0.33×`), never credited as clinically cleared.

![Potency-proven clinical TCR-T and ImmTAC gold standard versus OncoRef normal-tissue safety proxies](antigen-whitepaper/clinical-tcr-gold-standard.png)

*Figure 3.* **(A)** Observed clinical objective response rate (or post-HCT AML relapse-free protection) versus strongest somatic or cross-reactive normal proxy signal (in multiples of OncoRef's `Vital`/`Protected` gate ceiling) across all **$tcr_gold_n** curated endogenous-protein `TCR-T` and `ImmTAC` regimens (`1.0×` vertical dashed line = OncoRef `Unsafe` boundary; gray band below `3.5%` = `0% ORR` / unreadout trials excluded from safety calibration). **(B)** Multi-assay proxy attribution matrix showing which screening layer catches each toxic target (`×`) or calibrates each clinically validated CTA (`✓`).

<details markdown="1">
<summary>Full 16-regimen clinical TCR-T &amp; ImmTAC gold-standard calibration table (download CSV)</summary>

Download the complete machine-readable benchmark: [`clinical-tcr-gold-standard.csv`](antigen-whitepaper/clinical-tcr-gold-standard.csv).

$tcr_gold_standard_table

</details>

## Decision tree, attrition, and the 10× enrichment rule

Every protein group follows a deterministic decision tree with strict precedence—normal-tissue safety first, then coverage and sequence eligibility, and finally symmetric **$fold×** enrichment in tumors and reproductive/embryonic compartments—narrowing **$n** protein groups to **$shared** safety- and coverage-eligible proteins and **$union** retained candidates (`CA ∪ RA`).

![OncoRef classification decision tree](audits/target-selection-20261110/rebuild/summary/decision-tree.png)

*Figure 4.* Exact classification decision tree showing every binary gate, precedence order, and terminal leaf count across all **$n** protein groups.

![Proteome-wide attrition across screening gates](audits/target-selection-20261110/rebuild/summary/funnel-all.png)

*Figure 5.* Sequential attrition from **$n** protein groups to **$union** retained candidates, assigning each excluded protein once at its first failing gate. The broken x-axis above 3,000 keeps post-bulk filters readable. [Step counts (CSV)](audits/target-selection-20261110/rebuild/summary/funnel-all.csv).

Both enrichment routes apply the same threshold—**≥$fold×** above **H₀ = max(Normal N, Disposable D, $floor TPM)**—across **$assessed** protein groups with complete healthy-panel TPM (**$shared** shared-eligible groups plus **$other_compared** ineligible `Other` groups with complete bulk panels). Points above the dotted line in `Other` clear the fold cutoff numerically while remaining held for locus-specific single-cell attribution or unmeasured anatomical sub-niches (`TBATA`, `DSPP`, `CAGE1`).

![Cancer enrichment relative to the healthy adult baseline](audits/target-selection-20261110/rebuild/summary/cancer-enrichment.png)

*Figure 6.* Cancer enrichment (<span class="cat cat-ca">CA</span>) across all **$assessed** compared protein groups (values below 0.05 TPM clipped to the axis floor).

![Reproductive and developmental enrichment relative to the healthy adult baseline](audits/target-selection-20261110/rebuild/summary/reproductive-enrichment.png)

*Figure 7.* Reproductive and developmental enrichment (<span class="cat cat-ra">RA</span>) against the identical healthy adult baseline (values below 0.05 TPM clipped to the axis floor).

Exclusive category and locus counts across the coding proteome:

$partition

## Published CTA, oncofetal, and placental lists

OncoRef maps **$lit_sources** curated papers across four antigen families: classical CTA / cancer–germline catalogs, oncofetal / DUX4 / early-embryo programs, placental and reproductive catalogs, and focused single-gene or tumor-vasculature studies. Literature citations are provenance metadata; a published nomination never overrides normal-tissue expression.

Retention tracks how strictly each paper screened healthy tissues. Focused modern CTA and DUX4 studies (`Seager et al. 2024`, `Chew et al. 2019`, `Carter et al. 2023`, `Loriot et al. 2025`) retain **34%–67% CTAs**, whereas broad germ-cell, microarray, ZGA, and cancer–placenta lists (`Bruggeman et al. 2018`, `Smith et al. 2023`, `Jamin et al. 2021`, `Bradley et al. 2020`) are **92%–100% `Unsafe`** under modern adult somatic screening.

![Published CTA, oncofetal, and placental lists under OncoRef screening](audits/target-selection-20261110/rebuild/summary/literature-sources-retention.png)

*Figure 8.* Share of mapped protein groups (`n` shown per study) retained as <span class="cat cat-cta">CTA</span>, <span class="cat cat-ca">CA only</span>, or <span class="cat cat-ra">RA only</span> versus assigned to <span class="cat cat-other">Other</span> or <span class="cat cat-unsafe">Unsafe</span>. [Full 149-source table (CSV)](audits/target-selection-20261110/rebuild/summary/literature-sources.csv).

Of **$published_noms** published CTA or oncofetal nominations, **$published_surviving** survive into the candidate union (`CA ∪ RA`)—including **$cta_published** of our **$cta** CTAs—while **$cta_new** CTAs are newly identified by the expression rules and **$published_outside** nominations fail normal-tissue, eligibility, or enrichment screens.

![Overlap between rule-derived sets and published nominations](audits/target-selection-20261110/rebuild/summary/venn.png)

*Figure 9.* Overlap between rule-derived Cancer Antigens (**$ca**), Reproductive Antigens (**$ra**), and published CTA/oncofetal nominations (**$published_noms**).

![Attrition of published CTA and oncofetal nominations](audits/target-selection-20261110/rebuild/summary/funnel-published.png)

*Figure 10.* Sequential attrition of **$published_noms** published CTA/oncofetal nominations through the same gates, leaving **$published_surviving** candidates. [Step counts (CSV)](audits/target-selection-20261110/rebuild/summary/funnel-published.csv).

### Most-nominated individual proteins

Across **$studied** studied protein groups (**$published_noms** with a CTA nomination paper plus **$study_extra** clinical or familiar reference entries), below are the top ten by distinct curated CTA papers plus five familiar clinical anchors (`NY-ESO-1`, `PRAME`, `CT83`, `CTAG2`, `CAGE1`):

$publications

`N` and the reproductive peak may come from different assays or donors ([full context table](antigen-whitepaper/studied-proteins.csv)). `NR` ("not run") marks proteins excluded as `Unsafe` in Stage 1 before native-TPM comparison ($nr_proteins); tumor or testis libraries may still contain those transcripts.

### Studied proteins in Other or Unsafe

Of **$anchors** protein groups in the familiar or clinical reference panels, **$held** land in `Other` or `Unsafe`. Neither label predicts clinical toxicity for a specific therapeutic product; each records the strongest normal-tissue or rare-cell finding observed across the healthy adult panels:

$excluded

## Peak tumor expression and edge cases

Across the **$union** candidate union, **$union_ge_100** protein groups reach **≥100 TPM** in at least one QC-passing tumor sample, **$union_ge_1000** reach **≥1,000 TPM**, and **$union_patients_ge_10** have **≥10 known independent patients** crossing the $fold× cutoff. Median peak tumor expression is **$cta_median_tpm TPM** for `CTA` and **$ca_only_median_tpm TPM** for `CA only`.

Each table below ranks the top ten protein groups in its category by peak QC-passing tumor TPM ([all $union candidates](antigen-whitepaper/tumor-ranks.csv)); `Unsafe` proteins (`$unsafe_n`) exit in Stage 1 before tumor comparison. Because `Tumor fold` divides by `max(N, D, $floor TPM)`, fold equals `Tumor TPM` whenever both `N` and `D` are below $floor TPM (`LALBA`, `SCGB1D4`, `VCY/1B`). In `Other`, high tumor peaks (`TBATA`, `DSPP`, `CAGE1`) reflect unresolved single-cell or unmeasured sub-niche reviews rather than low expression. Peak TPM measures single-library abundance, not the fraction of malignant cells expressing the protein.

$peaks

### Worked edge cases

**`CTAG2` (LAGE-1) → `$ctag2_classification` (isoform-scoped clinical pMHC clearance, `#766`).** In peripheral bulk TPM, `CTAG2` clears both enrichment rules: its highest peripheral reference is **$ctag2_site** at **$ctag2_reference TPM** ($fold× cutoff **$ctag2_bar TPM**), testis reaches **$ctag2_repro TPM** (**$ctag2_fold×**), and ovarian cancer reaches **$ctag2_tumor TPM** across **$ctag2_patients** known patients at $fold×. In donor-resolved single-cell atlases, **lung** (alveolar type 2 cells) and **brain** (cortical glutamatergic neurons) pass the multi-study donor-replication screen while sharing **$ctag2_shared** of transcript 31-mers across the closed `Xq28` `CTAG1A`/`CTAG1B`/`CTAG2` cluster (`0.0%` external `Unsafe` 31-mer overlap):

$ctag2_cells

Under `oncoref.clinical_organ_clearance` ([audit `#766`](audits/clinical-organ-clearance-20261123/index.md)), adoptive `lete-cel` (`c259`, `HLA-A*02:01 / SLLMWITQC`) TCR-T infusions across synovial sarcoma, MRCLS, NSCLC, and myeloma cohorts without on-target pulmonary or urothelial injury clear `lung` and `urinary bladder` specifically for `ENST00000369585` (`LAGE-1a`, `180 aa`), while `#711` `Xq28` attribution (`0.15 TPM` bulk brain) resolves the cortical nuclear signal—restoring `CTAG2` to **`$ctag2_classification`**. Because `ENST00000247306` (`LAGE-1b`, `210 aa`) lacks `SLLMWITQC` and a separate heart review recorded a confounded day-161 cardiac event ([heart review](audits/ctag2-heart-evidence.md)), `CTAG2` carries `clinical_clearance_scope = "epitope"` and `clinical_unverified_isoforms = ["ENST00000247306"]` for downstream epitope selection.

**`PLAC1` & `DUX4` → `CTA` (nuclear `snRNA-seq` spillover `#765` + closed `D4Z4` repeat equivalence `#767`).** Intron-inclusive cortical `snRNA-seq` flags `PLAC1` across a `197.7 kb` non-coding 5'-UTR intron `5.2 kb` downstream of somatic `PABIR2` (`54.1 TPM`), yet cytoplasmic polyA bulk brain RNA is `0.0051 TPM`, somatic IHC is `0/4`, and MS immunopeptidomics is strictly `CANCER_ONLY`—clearing `PLAC1` (`152.6 TPM` Wilms tumor, `187.0 TPM` placenta) to `CTA` (`CTA-M`) with `rare_cell_status = "nuclear_only"`, whereas negative-control `MAGEA12` (`0.30 TPM` bulk brain, fatal neurotoxicity precedent) stays `Unsafe` ([audit `#765`](audits/nuclear-rna-spillover-20261123/index.md)). Symmetrically summing `DUX4` with all `20` primary-assembly `D4Z4` `DUX4L*` repeat units confirms `0.00 TPM` donor P95 in GTEx Skeletal Muscle (`FSHD` locus) and Skin, advancing cleavage-stage `DUX4` jointly under `#765 + #767` to `CTA` (`CTA-P`, [audit `#767`](audits/repeat-family-equivalence-20261123/index.md)).

**`GFY` (*Goofy*) & `SP7` → `Unsafe` (adult unmeasured-home panels `#735` & `2.5%` OR gDNA tumor QC).** `GFY` (`1,575.56 TPM` in `RT`) and `SP7` (`590.49 TPM` in `SARC_OS`) previously sat in `Other:Unscreened tissue` because olfactory epithelium and mineralized bone are absent from GTEx/HPA bulk panels. Integrating SHA-256-pinned adult human olfactory (`GFY`: **94.1%** of mature olfactory sensory neurons, `>300 TPM`) and mandible bone (`SP7`: **193.0 TPM** donor P95 vs. `4.0 TPM` Vital ceiling) panels—paired with hair-bearing GTEx skin `P99 ≤ 0.50 TPM` / `max < 4.0 TPM` bounds across `1,405` leg/suprapubic biopsies and a tightened **`≤ 2.5%` olfactory-receptor `≥ 5 TPM` genomic-DNA ceiling** across baseline-admitted tumor libraries (`92` gDNA-contaminated tumor samples newly excluded across `NEC_MERKEL`, `SARC_PEC`, `BL`, `MCL`, `HL`, `B_ALL`, `LAML`, and `T_ALL`)—reclassifies `5` ceiling-violating proteins (`GFY`, `SP7`, `PANX3`, `BPIFB3`, `GLYATL3`) to **`Unsafe`**, retains `33` home-gated proteins (`TAS1R2`, `KRTAP20-3`, `KRTAP21-2`, `KRTAP21-3`, `KRTAP6-2`) in `CA only` plus **`OR10J1` and `KRTAP20-4` promoted to `CTA` (`CTA-M`)** (alongside `OR1N1` and `OR4M2`), demotes `22` single-cohort gDNA-background `CA only` proteins (`17` `OR*`s, `KRTAP9-7`, `C19orf85`, `SCYGR8`, `SEBOX`, `VHLL` → `Other:Expression`) plus `10` single-cohort male-reproductive proteins (`CAPZA3`, `DNAJB8`, `IQCF2`, `IQCF6`, `MAGEB18`, `RNASE9`, `TMEM225`, `SPATA31A6`, `SMIM42`, `ZCCHC13` → `RA-M-only`), and keeps unmeasured sub-niches (`DSPP` odontoblast process, `CALHM3` type II taste cells) in `Other:Unscreened tissue`.

## Tissue injury tiers, patient profiles, and sex-restricted CTAs

OncoRef separates three tissue-policy questions: **anatomical role** (healthy baseline vs. <span class="cat cat-ra">RA</span>-qualifying reproductive/embryonic tissue vs. `thymic`), **adult injury tier** (clinical consequence of sustained organ or cell-type loss), and **patient-specific organ protection** (downstream anatomy and tolerance via `--profile`). A tissue can qualify <span class="cat cat-ra">RA</span> while retaining a protective rating—adult testis qualifies <span class="cat cat-ra">RA</span> while carrying a `Protected` injury tier, and breast biopsies are split into lactating versus resting mammary via a coordinated milk-program signature (`LALBA`, `CSN1S1`, `CSN2`, `CSN3`).

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

Biological <span class="cat cat-cta">CTA</span>, <span class="cat cat-ca">CA</span>, and <span class="cat cat-ra">RA</span> labels record where a protein is expressed, not what a particular patient can safely lose. Default biological sets leave reproductive tissues unprotected until you opt in. Downstream patient profiles (`female-patient`, `male-patient`, lactation, pregnancy, or custom overrides) mark each compartment as <span class="cat cat-protect">Protect</span>, <span class="cat cat-loss">Accept loss</span>, or <span class="cat cat-absent">Absent</span> ([rules and limits](host-selection.md)). Protecting testis narrows the **$cta** biological CTAs to **$testis_protected_cta** retained candidates; pregnancy leaves **$pregnant_retained** retained (**$pregnant_unresolved** unresolved due to incomplete developmental coverage, **$pregnant_excluded** excluded by expression in a protected organ).

![Patient anatomy and organ-loss choices across profiles](audits/target-selection-20261110/rebuild/summary/host-profiles-matrix.png)

*Figure 11.* Compartment choices (<span class="cat cat-protect">Protect</span>, <span class="cat cat-loss">Accept loss</span>, <span class="cat cat-absent">Absent</span>, or unprotected) across standard patient profiles.

![Patient-profile screening outcomes for biological CTAs](audits/target-selection-20261110/rebuild/summary/host-profiles-outcomes.png)

*Figure 12.* Downstream screening outcomes across all **$cta** biological CTAs under each patient profile. Filtering a candidate out for one patient narrows the target list and expands the peptide-subtraction background without relabeling its biological <span class="cat cat-cta">CTA</span> identity.

### Sex-restricted reproductive antigens and cross-sex CTAs

Because male and female reproductive organs are anatomically distinct, an antigen restricted to one sex's reproductive tract has no normal target organ in non-pregnant, non-lactating patients of the opposite sex. Crucially, `HostContext` (`oncoref.host_context`) uses the **same `Protected` + negative-IHC calibration as Stage 1 somatic safety (`8.0 TPM / 8.0 nTPM` base ceiling, raised `2×` to `16.0 TPM / 16.0 nTPM` in human compartments where routine IHC is `not_detected`)**—empirically anchored by the Clinical TCR Gold Standard (`IMA203` on `PRAME` at `10.93 TPM` adrenal / `12.0 nTPM` heart fibroblasts with negative IHC and `0` organ toxicity at `52.5%` ORR). Every Reproductive Antigen (**$ra** RAs) and Cancer–Testis Antigen (**$cta** CTAs) belongs to **one of four mutually exclusive host-physiology classes** (**`RA-M` / `CTA-M`**, **`RA-P` / `CTA-P`**, **`RA-F` / `CTA-F`**, **`RA-MF` / `CTA-MF`**, with `RA-{M,P,F,MF}-only` denoting `RA \ CA` rows and legacy aliases `MRA`, `PRA`/`RA-E`/`CTA-E`/`ERA`, `FRA`, `SRA`):

$sex_partition_table

- **Male-Restricted (`RA-M`, `$ra_m_n` RAs · `CTA-M`, `$cta_m_n` CTAs · `RA-M-only`, `$ra_m_only_n` — non-pregnant female-safe):** expressed in male reproductive organs (`testis`, `epididymis`, `prostate`, `seminal vesicle`) and cleanly clearing every non-pregnant, non-lactating female adult reproductive and mammary organ (`XAGE1A/B`, `SPANXA1/2`, `MAGEA4`, `NY-ESO-1`, `VCY/1B`, `PRM1`, `SPATA19`, `CT83`, `ODF1`, `ADAM2`, `DKKL1`, `SYCP1`, `GAGE2A`, `PAGE1`, `KRTAP20-4`, plus male-accessory/prostate targets such as `SEMG2`, `PRR23D1/2`, `SPAG11A`, `POTEH`, `EDDM3B`). Because **336 of 514** (`65.4%`) GTEx `Breast_Mammary_Tissue` biopsies come from genetically `XY` male cadavers (`#776`, `#779`) and **~60% of `XX` GTEx donors are postmenopausal (`ages 50–79`)**, OncoRef stratifies female `HostContext` exposures across **reproductive-age `XX` women (`ages 20–49`, `gtex.donor_age()`)**, integrates stage-replicated adult human folliculogenesis/oocyte RNA-seq (`GSE107746`, `128` QC-passing stages across primordial-to-antral oocytes and granulosa cells, `oncoref.reproductive_age_female`), and backfills GTEx-v11-missing `vagina`/`uterus` measurements (`TSSK2` → `RA-M-only`, leaving `0` human CTAs or RAs unresolved for missing female-tract coverage).
- **Pregnancy / Peripartum-Only (`RA-P`, `$ra_p_n` RAs · `CTA-P`, `$cta_p_n` CTAs · `RA-P-only`, `$ra_p_only_n` — adult-absent across non-pregnant, non-lactating men & women; legacy aliases `CTA-E` / `RA-E` / `ERA` / `PRA`):** expressed in `embryo`, `placenta`, or `lactating_breast` (`lactating_mammary`) and silent (`< 8` TPM/nTPM; `< 16` with `IHC−`) across **every** adult male and non-pregnant, non-lactating adult female reproductive and resting-mammary organ (`DUX4`, `ERVV-1`, `ERVV-2`, `LGALS16`, `BTN1A1`, `TBXT`, `MBD3L2`, `SSU72L2`, `TRIM48`, `TRIM43B`, `CER1`, `LEUTX`). Embryogenesis, placentation, and postpartum lactation represent a single physiological peripartum continuum absent from non-pregnant, non-lactating adult hosts; for the `$peri_cta_n` CTAs (`$peri_ra_n` RAs) active in `placenta` or `lactating_breast` that also express in adult testis (`PLAC1`, `XAGE3`, `MAGEA4/8/10`) or resting female breast (`LALBA`, `CSN2`), `HostContext` (`pregnant=True`, `lactating=True`) and `peripartum_antigen_ids()` enforce organ-level peripartum exclusions without breaking disjointness of the 4-way partition.
- **Female-Restricted (`RA-F`, `$ra_f_n` RAs · `CTA-F`, `$cta_f_n` CTAs · `RA-F-only`, `$ra_f_only_n` — male-safe):** expressed in non-pregnant adult female reproductive or resting-mammary compartments (`ovary`/adult oocytes, female tract, `XX resting_breast`, alongside lactation/embryo peaks) while cleanly clearing all male adult reproductive organs and non-secretory `XY resting_breast` under the calibrated `Protected` + `IHC−` gate (`CSN2`, `LALBA`, `ZP4`, `DPPA3`, `CPSF4L`, `IZUMO1R`, `ACCSL`, `OR1N1`, `KHDC3L`, `ZAR1L`, `VRTN`, `PRR23E`).
- **Shared Male & Female (`RA-MF`, `$ra_mf_n` RAs · `CTA-MF`, `$cta_mf_n` CTAs · `RA-MF-only`, `$ra_mf_only_n`):** exceed the calibrated `HostContext` ceiling (`≥8 TPM/nTPM`, `≥16` with `not_detected` IHC, or positive IHC) in adult reproductive or resting-mammary organs of both sexes (`PRAME`, `SCGB1D4`, `CTAG2`, `MOS`, `ACTL8`, `NLRP4`, `OOSP2`, `DMRTC2`, `SPO11`, `PIWIL3`, `FOXR1`).

The cross-sex asymmetry remains stark even after rescuing adult female oocyte expression: **$cta_female_safe_n** of **$cta** CTAs (`CTA-M + CTA-P`, **$cta_female_safe_pct%**) clear every adult female organ in non-pregnant, non-lactating women, whereas **$cta_male_safe_n** (`CTA-F + CTA-P`, **$cta_male_safe_pct%**) belong to the male-safe `CTA-F + CTA-P` sets, with **`CTA-M` outnumbering `CTA-F` by $cta_m_f_ratio to 1 (`$cta_m_n` vs. `$cta_f_n`)** and **`RA-M` outnumbering `RA-F` by $ra_m_f_ratio to 1 (`$ra_m_n` vs. `$ra_f_n`)** ([download full sex-restricted table](antigen-whitepaper/sex-restricted-antigens.csv)).

![Why female-lineage genes rarely become CTA-F and live cross-sex CTA exposure landscape](antigen-whitepaper/sex-restriction-funnel-and-landscape.png)

*Figure 13.* **(A)** Fate of **$female_audit_n** canonical mammary, oocyte/follicular, meiotic, Müllerian/endometrial, and placental protein groups across OncoRef's safety, enrichment, and sex-restriction gates ([full `$female_audit_n`-gene audit CSV](antigen-whitepaper/female-lineage-audit.csv)). **(B)** Live deconvolved adult male versus non-pregnant, non-lactating adult female exposure across all **$cta** human CTAs (`8 TPM / 8 nTPM` `Protected` base dashed lines; `16 TPM / 16 nTPM` `IHC−` ceiling dotted lines).

#### Why so few `CTA-F` and `RA-F` antigens exist—even across `BRCA_TNBC`, `OV`, `UCEC`, `UCS`, and `CESC`

Across the entire human coding proteome (`$n` protein groups), **10,000** protein groups reach `≥10 nTPM` in at least one non-placental adult female reproductive or mammary organ (`ovary`, `endometrium`, `fallopian tube`, `cervix`, `vagina`, `breast`), yet **9,993 (99.93%) fail Stage 1 normal-tissue safety (`Unsafe`)** (`9,980` by bulk somatic RNA, `9,924` by HPA somatic `scRNA-seq` clusters, and `6,805` by HPA somatic immunohistochemistry), leaving just **27** protein groups in the entire genome with non-placental female consensus `≥10 nTPM` and somatic consensus `<4 nTPM`. Auditing **$female_audit_n** canonical female-program genes across five biological lineages reveals four distinct quantitative mechanisms driving this attrition:

$female_audit_summary_table

1. **Evolutionary exocrine and ciliated mucosal homology (`Stage 1 Unsafe` eliminates 63%–94% of canonical female-organ markers).** Unlike seminiferous tubules—which evolved a specialized blood–testis barrier and dedicated chromatin/sperm-tail machinery—adult female reproductive and mammary epithelia share developmental origin and secretory machinery with vital or protected somatic organs:
   - *Mammary gland as a modified cutaneous apocrine gland:* Luminal breast markers massively expressed in `BRCA` and `BRCA_TNBC` remain active in **salivary serous/mucus acini, eccrine/apocrine skin glands, bronchus, and cornea**—excluding `PIP` (`GCDFP-15`, `141,212 TPM` `BRCA` / `15,008 TPM` `TNBC` → `16,000 nTPM` salivary gland), `SCGB2A2` (mammaglobin-A, `112,378 TPM` `BRCA` / `56,831 TPM` `TNBC` → `690 nTPM` skin), `MUCL1` (`85,652 TPM` `TNBC` → `740 nTPM` skin + `630 nTPM` salivary), `CSN3` (κ-casein, `14,566 TPM` `TNBC` → `260 nTPM` salivary mucus + `66 nTPM` bulk), `ELF5` (`1,044 TPM` `TNBC` → `230 nTPM` salivary), `ROPN1`/`ROPN1B` (`129–329 TPM` `TNBC` → `32–33 nTPM` salivary serous + corneal melanocytes), `VTCN1`/`B7-H4` (`3,716 TPM` `TNBC` → `52 nTPM` bronchial club + `17 nTPM` kidney), and `ANKRD30A`/`NY-BR-1` (`273 TPM` `BRCA` / `80.8 TPM` `TNBC` → `4.9 nTPM` cortical excitatory neurons).
   - *Müllerian tract (`endometrium`, `fallopian tube`, `endocervix`) as intermediate mesoderm and ciliated mucosa:* Gynecologic lineage markers highly expressed in `OV`, `UCEC`, `UCS`, and `CESC` share intermediate-mesoderm identity with **renal tubules** (`PAX8`: `2,535 TPM` `OV`, `681 TPM` `UCEC` → `350 TPM` kidney), secretory transport with **choroid plexus, alveolar type 2 cells, bronchial club cells, and pancreatic ducts** (`FOLR1`: `4,869 TPM` `OV`, `1,598 TPM` `UCEC` → `2,900 nTPM` choroid plexus + `800 nTPM` renal proximal tubule; `SLC34A2`/`NaPi2b`: `2,966 TPM` `OV`, `927 TPM` `UCEC` → `2,200 nTPM` lung AT2; `WFDC2`/`HE4`: `41,779 TPM` `UCEC`, `20,031 TPM` `OV` → `18,000 nTPM` bronchus; `MUC16`/`CA-125` & `MSLN` → corneal epithelium + bronchus; `PAEP`/glycodelin & `CLDN6` → `19–65 nTPM` pancreatic ducts; `OVGP1`: `3,501 TPM` `UCEC` → `65 TPM` naive CD4 T cells + `20 TPM` cerebellum), while oviduct-enriched structural genes (`CFAP141`, `ANKRD66`, `DNAAF6`) are motile-cilia axonemal subunits shared with respiratory bronchial ciliated cells and ependyma.
   - *Placental hormones and syncytial retroelements:* Trophoblast hormones reactivated in `CESC`, `UCS`, `UCEC`, and basal `BRCA` are shared with or paralogous to **pituitary gonadotrophs/somatotrophs and endocrine pancreas** (`CGA`: `8,938 TPM` `BRCA_Basal`, `3,118 TPM` `CESC` → `4,300 nTPM` pituitary; `CGB3/5/8`: `318 TPM` `CESC`, `182 TPM` `UCS` → `28 nTPM` pancreas + `13 nTPM` pituitary; `GH2` → `100 nTPM` pituitary; `ERVW-1`/syncytin-1: `227 TPM` `UCEC` → `24 nTPM` brain white matter; `HLA-G`: `718 TPM` `BRCA` → `320 TPM` memory Tregs; `ALPP` → `72 nTPM` lung AT1).

2. **Adult testis pachytene transcriptional scanning captures shared germline, endometrial, and placental antigens into `CTA-MF` and `CTA-M`.** Human pachytene spermatocytes and round spermatids open chromatin genome-wide and transcribe `>80%` of protein-coding loci, whereas adult oocytes arrest in dictyate prophase I without analogous scanning:
   - *Endometrial/ovarian tumor antigens active in both tracts (`CTA-MF`):* The strongest safe antigens in `UCEC`, `UCS`, and `OV`—including `SCGB1D4` (`3,576 TPM` `UCEC`, `1,384 TPM` fallopian tube, `383 TPM` uterus → `51.2 TPM` `epididymis` + `20.2 nTPM` `testis`), `PRAME` (`505 TPM` `UCEC`, `470 TPM` `UCS`, `59.3 TPM` `ovary` → `310.9 TPM` `testis`), and `CTAG2` (`604 TPM` `OV`, `77.0 TPM` `vagina`, `35.7 TPM` `ovary` → `101.9 TPM` `testis`)—strongly exceed the calibrated ceiling in both sexes (`CTA-MF`).
   - *Oocyte/SCMC and meiotic factors co-expressed in spermatogenesis (`CTA-MF`):* Even when postmenopausal bulk `ovary` looks low, HPA ovarian oocyte `scRNA-seq` and adult `testis` confirm substantial dual-sex germline expression across `FOXR1` (`660 nTPM` oocytes, `23.9 TPM` testis), `OOSP2` (`540 nTPM` oocytes, `103.9 TPM` testis), `ACTL8` (`401 nTPM` oocytes, `38.5 TPM` testis), `NLRP4` (`32.5 TPM` `UCEC`, `56.6 nTPM` oocytes, `18.5 TPM` testis), `PIWIL3` (`30.9 nTPM` oocytes, `40.8 TPM` testis), `DMRTB1`, `ADAD1`, and `MPC1L` (`CTA-MF`), whereas oocyte/maternal factors with only low IHC-negative testis background (`DPPA3` at `5.6 TPM`, `VRTN` at `6.9 TPM`, `PRR23E` at `4.6 TPM` testis) cleanly classify as **`CTA-F`**, and `SYCP1`/`DKKL1`/`GAGE2A`/`PAGE1` (`<16` IHC-negative female background) classify as **`CTA-M`**.
   - *Placental `CT-X` antigens silent in non-pregnant women (`CTA-M` & `CTA-P`):* Trophoblast-expressed X-linked CTAs frequently derepressed in `BRCA_TNBC`, `UCEC`, `UCS`, `OV`, and `CESC`—such as `XAGE3` (`234 TPM` `TNBC`, `206 TPM` `UCEC`), `PLAC1` (`88.5 TPM` `BRCA`, `52.5 TPM` `UCS`, `47.8 TPM` `OV`), `CT83` (`193 TPM` `CESC`, `90.4 TPM` `TNBC`), and `MAGEA10` (`78.6 TPM` `CESC`, `53.6 TPM` `TNBC`)—are silent across non-pregnant adult female organs (`<2–3.3 TPM`) while expressed at `4.1–445 TPM` in adult `testis`, making them **female-safe `CTA-M` targets** for `TNBC` and gynecologic cancers (with placental members simultaneously indexed in `CTA-P`)!

3. **Postmenopausal bulk ovary dilution rescued by reproductive-age GTEx (`20–49`) and adult human oocyte profiling (`GSE107746`) — yet autosomal oocyte genes stay epigenetically silent in most epithelial female tumors (`RA-F-only` / `CTA-F`).** GTEx (`n = 192`, ~60% postmenopausal ages `50–79`) and HPA bulk `ovary` sample follicularly depleted ovarian cortex where oocytes represent `<0.01%` of tissue mass, diluting pure oocyte mRNAs `100–1,000×` (`BMP15 = 1.36 TPM`, `WEE2 = 1.44 TPM`, `PABPN1L = 1.54 TPM`, `ZAR1L = 1.24 TPM`, `KHDC3L = 3.85 TPM`). Stratifying GTEx female organs to reproductive-age `XX` women (`ages 20–49`) and integrating stage-replicated adult human folliculogenesis RNA-seq (`GSE107746`, `128` QC-passing stages across primordial, primary, secondary, and antral oocytes and granulosa cells) directly quantifies adult female germ-cell expression—promoting **`CPSF4L`** (`265.7 TPM` adult oocytes, male `≤3.69 TPM`) from `CA only → CTA-F` and **`MOS`** (`734.1 TPM` adult antral oocytes + `10.8 nTPM` spermatocytes/spermatids) from `CA only → CTA-MF`, reclassifying **`IZUMO1R` (`JUNO`, `33.6 TPM` adult oocytes)**, **`ACCSL` (`28.7 TPM`)**, and **`OR1N1` (`13.5 TPM`)** into **`CTA-F`** (`12` total `CTA-F`s), moving **`BMP15` (`369.5 TPM` adult oocytes)** and **`OR52D1` (`11.2 TPM`)** into **`RA-F-only`**, and unmasking adult oocyte transcription (`8.6–156.1 TPM`) across `11` dual-sex meiotic/germline CTAs (`DMRTC2`, `EBLN1`, `LYPD4`, `DCAF4L2`, `MAGEC1`, `FTHL17`, `MAGEB16`, `CCDC42`, `ANKRD60`, `SPO11`, `TGM4` → `CTA-MF`). Even with adult oocyte expression fully resolved, autosomal maternal-effect genes require oocyte-specific transcription factors (`FIGLA`, `NOBOX`, `LHX8`) rather than passive DNA hypomethylation and stay `<4–11 TPM` across epithelial `OV`, `UCEC`, `UCS`, and `CESC` (`PABPN1L`, `WEE2`, `BMP15`, `TCL1B`, `OR52D1` → `RA-F-only`).

4. **Single-assay somatic gates on remaining canonical oocyte factors.** The remaining textbook oocyte/zona genes fail specific somatic screens: `ZP1` (`18 TPM` pituitary), `ZP2` (`100.5 TPM` `TNBC` → `59 TPM` cerebellum + IHC), `ZP3` (`143.8 TPM` `UCEC` → `150 nTPM` bronchial basal cells), `NLRP5`/`MATER` (`50.8 TPM` `BRCA` → `22 nTPM` parathyroid), `OOEP` (`49.7%` alveolar type 2 rare-cell UMI), `FIGLA` and `PADI6` (`14 TPM` in `DICE` M2 macrophages / NK cells vs. `8 TPM` immune ceiling, and `<10 TPM` in all female tumors), and `TUBB8` (`10.3 TPM` `UCEC`, `121.8 TPM` `DSRCT` → `Other:Review` due to β-tubulin paralog cross-mapping in cortical `snRNA-seq`).

All **$ra_f_n** human `RA-F` antigens (**$cta_f_n** `CTA-F` and **$ra_f_only_n** `RA-F-only`):

$ra_f_all_table

#### Case study: why both `CSN2` (β-casein) and `LALBA` (α-lactalbumin) are `CTA-F` (and `BTN1A1` is `CTA-P`)

Milk-droplet proteins are the canonical expectation for female-specific breast cancer targets (`LALBA` reaches `79,700 TPM` in `BRCA`, `2,877 TPM` in `BRCA_TNBC`, and `45,625 TPM` in lactating breast and is an active clinical `TNBC` vaccine target, e.g. `NCT04674306`; `CSN2` reaches `99,762 TPM` in `BRCA` and `988 TPM` in `BRCA_TNBC`; `BTN1A1` reaches `427.5 TPM` in `BRCA`/`BRCA_TNBC`). Interrogating all **514** deconvolved GTEx `Breast_Mammary_Tissue` biopsies (`178 XX` female, `336 XY` male) alongside **424** GTEx/HPA `Testis` and **282** `Prostate` libraries—and calibrating `HostContext` to OncoRef's `Protected` (`8 TPM / 8 nTPM`) + negative-IHC (`2×` = `16`) gate—explains why **`CSN2` and `LALBA` both classify as `CTA-F`**, whereas **`BTN1A1` classifies as `CTA-P`**:

![Donor-level interrogation of LALBA, CSN2, and BTN1A1 across XX/XY mammary tissue, male testis, and female cancers](antigen-whitepaper/lalba-csn2-mammary-testis-interrogation.png)

*Figure 14.* **(A)** Donor-level `LALBA`, `CSN2`, and `BTN1A1` expression across `514` GTEx `Breast_Mammary_Tissue` biopsies split by genetic sex (`XX` vs. `XY`) and coordinated 4-gene milk-secretory state, isolating `8/336` (`2.4%`) `XY` cadaver galactorrhea donors (`XIST ≤ 0.09 TPM`, `RPS4Y1 = 178–280 TPM`) from `328` non-secretory resting male breasts. **(B)** Adult male reproductive-organ distributions (`414` GTEx `Testis`, `282` GTEx `Prostate`, and HPA `testis` germ-cell `scRNA-seq`), showing complete male silence of `CSN2` and `BTN1A1` and why `LALBA`'s `4.87 TPM` P95 / `7.76 TPM` P99 spermatid transcription with negative `Enhanced` testis IHC clears the calibrated `8 TPM` `Protected` / `16 TPM` `IHC−` gate (`LALBA → CTA-F`). **(C)** Peak tumor TPM across `BRCA`, `BRCA_TNBC`, `OV`, `UCEC`, `UCS`, and `CESC` for retained candidates versus excluded somatic-sharing clinical markers.

- **Isolating ectopic male galactorrhea in `XY` GTEx `Breast_Mammary_Tissue` clears resting male breast (`≤2.34 TPM`) across all three milk proteins.** Across all `336` genetically confirmed `XY` male breast donors (`XIST = 0.00–0.09 TPM`, `RPS4Y1 = 178–280 TPM`, `DDX3Y = 54–107 TPM`—ruling out sample swaps or female splash), **97.6% (`328/336`) have near-zero milk-gene RNA** (`XY` resting `P50 = 0.00 TPM`; `P95 = 0.10 TPM` for `CSN2` and `0.38 TPM` for `LALBA`; purity-weighted `host_exposure` `P99 = 2.14 TPM` for `CSN2`, `2.34 TPM` for `LALBA`, and `0.67 TPM` for `BTN1A1`). Only **8 hospitalized/ICU male cadaver donors** (`GTEX-11PRG` at `LALBA = 1,797 / CSN2 = 1,572 / CSN1S1 = 675 / CSN3 = 377 TPM`, `GTEX-139T8` at `327 / 62 / 195 / 62 TPM`, `GTEX-XMD1` at `161 / 156 / 287 / 180 TPM`, `GTEX-R55F` at `63 / 202 / 81 / 8 TPM`, plus `GTEX-1CB4J`, `GTEX-1P4AB`, `GTEX-1N2DW`, `GTEX-T5JC`) exhibit **coordinated 4-gene prolactin-driven secretory induction** (male galactorrhea from critical-illness hyperprolactinemia, dopamine-antagonist antiemetics/antipsychotics, cirrhosis, or spironolactone). Separating coordinated milk-operon activation (`≥3` of `{LALBA, CSN1S1, CSN2, CSN3}` `≥5 TPM` with peak `≥15 TPM`) from physiological `resting_breast`—and unifying peripartum `embryo` + `placenta` + `lactating_breast` into `CTA-P`—cleanly promotes **`CSN2` (`male_sex_organ_max_tpm = 1.22 TPM`, `XX` breast glandular `Low` IHC) to `CTA-F`** and **`BTN1A1` (`male_sex_organ_max_tpm = 2.11 TPM`, `XX resting_breast = 1.02 TPM`, `86.2 TPM` exclusively in `lactating_breast`) to `CTA-P`**.
- **Why `LALBA` (`45,625 TPM` lactating breast vs. `7.76 TPM` IHC-negative testis) belongs in `CTA-F` rather than `CTA-MF`.** Unlike `CSN2` (`0/414` GTEx testes `≥4 TPM`), `LALBA` (`ENSG00000167531`, `0.0%` 31-mer paralog overlap, `ENST00000301046` MANE Select) exhibits genuine low-level spermatogenic transcription reflecting its evolutionary origin as a duplicate of **c-type lysozyme (`LYZ` / sperm-acrosomal `LYZL1/2/4/6`)**: across **414** deconvolved GTEx `Testis` biopsies, `LALBA` measures **median `2.57 TPM`, P95 `4.87 TPM`, `host_exposure` P99 `7.76 TPM`, max `9.24 TPM`**, with **4.3 nTPM** in HPA late spermatids (`seminal vesicle` `3.58 TPM`, `prostate` `0.22 TPM`, `epididymis` `0.00 TPM`), while **HPA `Enhanced`-reliability immunohistochemistry is `Not detected` across `testis` (`seminiferous ducts` and `Leydig cells`), `epididymis`, `prostate`, and `seminal vesicle`**. Applying an unadjusted `<4 TPM` (`Vital`-organ) gate to `testis` would penalize `~5 TPM` untranslated spermatid chromatin leakage (`6,000×` below lactating mammary RNA!) more strictly than Stage 1 screens somatic `Protected` organs (`8 TPM` base, `16–64 TPM` with negative IHC—validated clinically by `IMA203` on `PRAME` at `10.93 TPM` adrenal / `12.0 nTPM` heart fibroblasts with `0` organ toxicity). Under OncoRef's harmonized `Protected` (`<8 TPM / <8 nTPM`) + negative-IHC (`<16 TPM / <16 nTPM`) `HostContext`, `LALBA` clears every adult male compartment (`7.76 < 8.0` base and `< 16.0` `IHC−`) while crossing in `XX resting_breast` (`160.2 nTPM` luminal cells, `12.8 TPM` bulk) and `lactating_breast` (`45,625 TPM`), classifying **`LALBA` as `CTA-F`** (and excluded under `HostContext(lactating=True)`); callers seeking an ultra-strict `<4 TPM` unadjusted RNA screen can still pass `HostContext(protect=("testis",), max_tpm=4.0, negative_ihc_multiplier=1.0)`.

<details markdown="1">
<summary>Representative gene-by-gene disposition table across mammary, oocyte, meiotic, Müllerian, and placental lineages (download all $female_audit_n genes as CSV)</summary>

Download the complete machine-readable audit: [`female-lineage-audit.csv`](antigen-whitepaper/female-lineage-audit.csv).

$female_audit_highlights_table

</details>

## Methods: protein grouping, sources, and thresholds

### Protein-group counting unit

Source gene identifiers across RNA, single-cell, and protein atlases map to human **Ensembl 115** loci. Each locus gets one qualified reference translation—preferring **MANE Select** on the primary assembly, then **Ensembl canonical**—checked for coding biotype, complete start/stop codons, uninterrupted ORF, and absence of translation conflicts ([transcript selection audit](audits/cta-transcript-selection-20261108/index.md)).

Paralogous loci encoding identical proteins (`SPANXA1/2`, `XAGE1A/B`, `VCY/1B`) collapse into one protein group by representative amino-acid sequence even when cDNAs differ. Wherever sample identities exist, member-locus expression is summed within each sample before comparing across cohorts. The **$missing_seq** records lacking a qualified representative sequence stay visible; missing translation blocks admission even when RNA is measured ([sequence and peptide-subtraction API](antigen-sequences.md)).

### Evidence sources and native-TPM coverage

OncoRef pins **$files** SHA-256-verified files across **$tissue_rows** source/context rows and **$cohorts** assessed tumor cohorts. Incompatible assay units are never pooled; every dataset stays in its native scale and enters only the layer it supports:

| Evidence pillar | Primary resources | Native statistic and role |
|---|---|---|
| **Healthy baseline (`N`, `D`)** | GTEx v11, HPA v23 biopsies and brain, HPA/DICE/Monaco sorted immune | Within-site donor or biopsy P95 in native TPM (optional `mean_plus_2std` sensitivity mode) sets the non-thymic healthy reference (`N`) and Disposable maximum (`D`). |
| **Normal-tissue safety screens** | HPA consensus and scRNA clusters, HPA normal IHC, CELLxGENE Census + Siletti, Collin cornea, specialized anatomy | Tiered RNA ceilings (TPM/nTPM/CP10K), reliable positive IHC, matched negative-IHC credit, and consortium-deduplicated rare-cell UMI screens. |
| **Tumor and reproductive/embryo RNA** | $cohorts QC-passing tumor cohorts, GTEx/HPA reproductive biopsies, Petropoulos/Tyser/Molè embryos | Individual QC-passing tumor libraries qualify `CA`; individual reproductive samples or supported embryo lineages qualify `RA`. |
| **Curated literature and clinical controls** | $lit_sources catalog-mapped papers (CTA, oncofetal/DUX4, placental, clinical trials, normal counterevidence) | Provenance metadata and curated normal-expression exclusions; positive nominations alone never grant admission. |

![Native-TPM coverage across somatic, immune, disposable, and thymic tissues](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage-somatic.png)

*Figure 15.* Native-TPM tissue coverage across somatic organs, immune subsets, disposable tissues, and thymus (empty out-of-panel columns omitted). Colored cells mark qualifying coverage for at least one assessed protein group; gray marks no qualifying coverage.

![Native-TPM coverage across neural, ocular, reproductive, and developmental tissues](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage-neural-repro.png)

*Figure 16.* Native-TPM coverage across brain, eye, nerve, reproductive, and embryonic tissues. **Qualifies RA** marks individual reproductive samples or supported embryo-lineage means capable of qualifying `RA`; **Context only** marks tissue summaries retained for display.

<details markdown="1">
<summary>Exact quantitative gate thresholds and rare-cell provenance rules</summary>

* **Bulk, cluster, and immune RNA ceilings:** Base limits in TPM and nTPM are **$rna_bases**. Matched negative IHC across routine contexts of the tissue raises the ceiling by 2× (Approved), 4× (Supported), or 8× (Enhanced), alongside limited registry clinical adjustments. Exceeding a `Vital` or `Protected` ceiling establishes `Unsafe`; `Manageable` establishes `Unsafe` only above **$manageable_ceiling**.
* **Donor-resolved rare-cell RNA and consortium provenance:** Evaluated in native UMI counts per adult donor and named cell type. In a `Vital` population, a donor signal requires **≥$cells cells** and **≥$pct%** of cells with ≥2 molecules (replicated across independent study groups), or **≥$frequent%** in a single study; `Protected` populations double those frequency cutoffs. Shared consortium collections (Baylor/Chen eye, Sanes eye, Allen reference brain, Michigan kidney, Teichmann heart) and donors embedded in integrated atlases (HLCA, Baylor/Chen eye, SEA-AD) are resolved to their original study group and biological donor ID before counting replication. Replicated signals without locus mapping risk (or with independent bulk or gene-specific confirmation) establish `Unsafe`; unresolved locus attribution or single-study signals go to `Other:Review`.
* **Structural and home-tissue gates:** Unsupported nested sequences, read-through loci, or ≥95% 9-mer containment in `Unsafe` proteins produce `Other:Sequence`. A curated adult normal tissue outside the native-TPM panels (`GFY` in olfactory neurons) assigns `Other:Unscreened tissue` for `CA` while still allowing independent `RA` qualification; fetal-only tissues (`AMELX` in developing tooth epithelium) keep their developmental scope.
* **Eligibility and Disposable cap:** Eligibility requires complete coverage across **$required** healthy contexts, target peak above every measured non-reproductive and Disposable value, and Disposable expression ≤ **$cap** TPM and ≤ **$cap** nTPM.
* **Embryo lineage support for RA:** Well-sampled lineages require ≥20 cells, ≥5 cells at ≥1 TPM, and ≥10% positive cells; smaller lineages require ≥5 cells, ≥50% positive cells, and replication in two independent embryos of the same stage and lineage. Within-site donor/biopsy **P95** dampens isolated healthy outliers but cannot diagnose carryover or prove absence from an unmeasured rare cell type.

</details>

<details markdown="1">
<summary>Source-family rules and versioned file inventory ($files files across $tissue_rows contexts)</summary>

| Source family | Evidence provided | Role in classification |
|---|---|---|
| **Ensembl 115, MANE and sequence registry** | Locus IDs, transcript flags, CDS qualification, sequence digests, read-through maps | Define protein-group counting units and test translation, nesting, and `Unsafe` 9-mer containment. |
| **GTEx v11** | Individual site TPM, donor/library IDs, site medians, QC annotations | Non-reproductive sites contribute the within-site linear donor P95 (after within-donor replicate median). Individual reproductive samples can independently qualify `RA`. |
| **HPA bulk tissue and brain** (incl. GTEx via HPA) | Individual tissue biopsies, tissue/brain-region means, consensus nTPM, anatomical mappings | Native TPM enters the healthy reference (within-tissue biopsy P95 where available) or `RA`-qualifying reproductive samples. Consensus nTPM runs its own normal gate. Redistributed GTEx overlaps GTEx v11. |
| **HPA sorted blood, DICE and Monaco** | Sorted immune populations, donor TPM, ranked values, published means, immune nTPM | HPA donor P95 and DICE single-locus P95 enter native comparison; DICE multi-locus groups use the conservative sum of member maxima when donor pairing is unavailable; Monaco retains the published mean. Also evaluated in the immune screen. |
| **HPA single-cell clusters and normal IHC** | Cell-type nTPM, staining grades, antibody IDs, reliability scores | Screen normal cell clusters and reliable positive protein findings (including resolved cell-pattern rules); matched negative IHC can raise RNA ceilings. Reproductive staining/nTPM aids interpretation but never qualifies $fold× `RA` on its own. |
| **CELLxGENE Census and Siletti brain** | Versioned H5AD studies, donor IDs, cell types, molecule counts | Grade sparse normal-cell expression by donor penetration, independent study replication, locus mapping risk, and independent confirmation. |
| **Collin cornea and conjunctiva** | Adult, developmental, and disease ocular cell counts with donor metadata | Dedicated adult corneal normal gate in counts per 10,000 (CP10K); non-adult and disease material stays scoped separately. |
| **Specialized anatomy studies** | Olfactory epithelium (Durante, Olender, Verbeurgt, MSSM), fetal dental epithelium, developmental skeletal atlas, adult osteoblast markers | Document localized expression in tissues underrepresented in bulk atlases; age, scale, and marker-list scope constrain what each study can establish (marker lists never prove tissue-wide absence). |
| **Petropoulos, Tyser and Molè embryos** | Per-embryo lineage TPM (preimplantation, Carnegie stage 7) plus day-9/day-11 hypoblast UMI/staining | Supported native TPM can qualify `RA`; hypoblast UMI counts and staining document developmental expression on their native scales. |
| **Tumor cohorts** | Per-sample matrices, sample/patient IDs, cross-cohort sharing, scale and QC checks | **$cohorts assessed cohorts** in the [tumor coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv). Comparable QC-passing libraries qualify `CA` (germ-cell tumors eligible; coordinated reproductive RNA advisory). |
| **Curated papers, proteomics and trials** | Nominated gene lists, normal counterevidence, clinical/epitope records, CT83 RNAscope | Separate positive nominations from reviewed normal counterevidence and literature benchmarks. Curated normal findings can exclude a protein; a positive nomination alone never grants admission. |

$sources

</details>

## Comparative multi-species framework: human, dog, and mouse

OncoRef applies one mathematical definition of <span class="cat cat-ca">CA</span>, <span class="cat cat-ra">RA</span>, <span class="cat cat-cta">CTA</span>, and the 4-way physiological sex-restriction partition (`RA-M`/`CTA-M`, `RA-F`/`CTA-F`, `RA-P`/`CTA-P`, `RA-MF`/`CTA-MF`; legacy aliases `MRA`/`FRA`/`ERA`/`PRA`/`SRA`/`CTA-E`/`RA-E`) across **human (`9606`)**, **domestic dog (`9615`)**, and **laboratory mouse (`10090`)**, while parameterizing each species with its own verified genome reference, comparative mammalian anatomy, native expression datasets, and unit-specific cutoffs (`oncoref.species_policies`, `oncoref.species_target_selection`). No expression value, safety exclusion, or threshold credit is ever transferred across species by orthology ([full multi-species whitepaper](species-antigen-whitepaper.md), [dataset bindings CSV](antigen-whitepaper/species-datasets.csv), [cross-species ortholog concordance CSV](antigen-whitepaper/cross-species-orthologs.csv)).

$species_summary_table

### *Homo sapiens* (taxon 9606)

All eight screening pillars are active across **$n** Ensembl 115 protein groups (**$loci** coding loci), combining GTEx v11 (reproductive-age `XX` `20–49` female stratification) and HPA v23 donor P95 baselines, HPA/DICE/Monaco immune panels, consortium-deduplicated CELLxGENE and Siletti rare-cell UMI atlases, HPA immunohistochemistry, `GSE107746` adult ovarian follicles, Petropoulos/Tyser/Molè embryos, and **$cohorts** QC-passing tumor cohorts (`$cta_m_n` `CTA-M`, `$cta_p_n` `CTA-P`, `$cta_f_n` `CTA-F`, `$cta_mf_n` `CTA-MF`).

### *Canis lupus familiaris* (taxon 9615)

Canine comparative oncology is evaluated genome-wide across **$dog_n** Ensembl 98 `CanFam3.1` (`GCA_000002285.2`, cross-referenced to `RefSeq 106` `UU_Cfam_GSD_1.0`) protein groups (**$dog_loci** coding loci, **7/8** active pillars): **$dog_cta** `CTA` (`$dog_cta_m_n` `CTA-M`, `$dog_cta_p_n` `CTA-P`, `$dog_cta_f_n` `CTA-F`, `$dog_cta_mf_n` `CTA-MF`), **$dog_ca_only** `CA only`, and **$dog_ra_only** `RA only` across **8** spontaneous companion-dog tumor cohorts (`OSA_NCI`, `OSA_COTC030`, `OSA_CCOGC`, `CMT`, `HSA`, `LYM`, `DLBCL`, `GLIOMA`), `Bgee 15.2` + `GSE219045` + `GSE119810` normal/reproductive linear `TPM` panels, and `Pfizer` 38-tissue quantitative TMT proteomics (`MSV000099744` with within-plex bleed subtraction). Across the **14** literature-curated canine RT-PCR nominations (`hindriks2025`, `nemec2019`), only **`PRAME` (`$dog_prame_call`, `CTA-M`)** survives quantitative somatic safety and tumor enrichment, while `MAGEB10` and `SSX6` classify as `Reproductive Antigen only` and somatic-expressed nominees (`AKR1E2`, `IGF2BP3`, `NT5C1B`, `PAGE3-like`, `SPECC1`, `TPX2`) are excluded as `Unsafe` ([canine whitepaper section](species-antigen-whitepaper.md#canis-lupus-familiaris-domestic-dog-taxon-9615)).

### *Mus musculus* (taxon 10090)

Murine preclinical screening is anchored on the complete Ensembl 115 GRCm39 (`GCA_000001635.9`, C57BL/6J) proteome (**$mouse_n** protein groups across all **8/8** active pillars under uniform `≤ 2.5%` ectopic `Olfr` gDNA library QC: **$mouse_cta** `CTA` [`$mouse_cta_m_n` `CTA-M`, `$mouse_cta_p_n` `CTA-P`, `$mouse_cta_f_n` `CTA-F`, `$mouse_cta_mf_n` `CTA-MF`], **$mouse_ca_only** `CA only`, **$mouse_ra_only** `RA only`). Murine anatomy includes `coagulating_gland` (anterior prostate) and `seminal vesicle` on the male axis and `uterine_horn` on the female axis while omitting `resting_mammary` from adult male protection (fetal androgen ablation at E14.5) and omitting `appendix` and palatine `tonsil` from `Disposable` (`D`). Across all **$ortholog_total_keys_n** distinct mammalian `CTA` ortholog families spanning human, dog, and mouse, only **$conserved_cta_pairs_n** (`CSN2`, `PRAME`, `PSMA8`, `ANTXRL`, `LUZP4`, `TRIM42`, `PASD1`, `PLAC1`/`Plac1`, `TEX19`/`Tex19.1`) independently qualify as `CTA` in two or more species—demonstrating why CTA catalogs must be constructed natively per species rather than projected by orthology ([murine & cross-species orthology section](species-antigen-whitepaper.md#cross-species-cta-orthology-across-human-dog-and-mouse)).

## Reproducing the catalog and downloads

Every external reference file and per-sample tumor matrix required for a clean rebuild (policy `$policy`) is SHA-256-pinned and fetched through OncoRef's uniform CLI (no out-of-band downloads or machine-local paths):

```sh
# 1. Fetch and verify all pinned reference datasets (add --include-tumors for all 132 tumor cohorts)
python -m oncoref.curation fetch --include-tumors

# 2. Rebuild the full curation graph (catalog, host exposure, figures, reports, whitepaper;
#    pass --normal-stat mean_plus_2std to gtex_sample_evidence.py / bulk_normal_percentiles.py for sensitivity runs)
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
| [Interactive protein catalog](audits/target-selection-20261110/rebuild/report.html) | Searchable browser for all **$n** protein groups, per-gate reasons, top samples, and literature links |
| [Summary figures and distributions](audits/target-selection-20261110/rebuild/summary.md) | Threshold heatmaps, scale distributions, rare-cell plots, and stage breakdowns |
| [Classification rules](target-selection.md) · [Build guide](curation.md) | Policy contract, CLI commands, and reproducible build instructions |
| [Tissue injury hierarchy](tissue-hierarchy.md) · [Tissue tiers](tissue-tiers.md) · [Patient-specific filtering](host-selection.md) | Cited adult organ and cell-type injury assessments and patient-profile filtering |
| [Sequence API and peptide subtraction](antigen-sequences.md) · [Data download API](antigen-data.md) | Ensembl 115 reference translations, k-mer subtraction, and programmatic dataset retrieval |
| [Registered file inventory (`sources.csv`)](antigen-whitepaper/sources.csv) | All **$files** pinned source files, catalog selectors, and SHA-256 digests |
| [Studied CTA dispositions (`studied-proteins.csv`)](antigen-whitepaper/studied-proteins.csv) | Full measurements and exclusion reasons for all **$studied** published or clinical anchor proteins |
| [Candidate tumor rankings (`tumor-ranks.csv`)](antigen-whitepaper/tumor-ranks.csv) | Peak tumor TPM, folds, cohort/sample IDs, and patient counts for all **$union** candidates |
| [Curated literature sources (`literature-sources.csv`)](audits/target-selection-20261110/rebuild/summary/literature-sources.csv) | Disposition counts across all **$lit_sources** catalog-mapped papers grouped by antigen class |
| [Clinical TCR gold-standard benchmark (`clinical-tcr-gold-standard.csv`)](antigen-whitepaper/clinical-tcr-gold-standard.csv) | Potency-proven human `TCR-T` and `ImmTAC` target, response, and clinical toxicity curation (**$tcr_gold_n** regimens) evaluated against OncoRef's safety proxies |
| [Sex-restricted antigens (`sex-restricted-antigens.csv`)](antigen-whitepaper/sex-restricted-antigens.csv) | Disjoint 4-way `RA-M`/`CTA-M`, `RA-F`/`CTA-F`, `RA-P`/`CTA-P`, and `RA-MF`/`CTA-MF` calls across all **$ra** human RAs with organ TPMs |
| [Female-lineage CTA audit (`female-lineage-audit.csv`)](antigen-whitepaper/female-lineage-audit.csv) | Gene-by-gene disposition, `BRCA`/`TNBC`/`OV`/`UCEC`/`UCS`/`CESC` tumor peaks, and gating normal/male measurements across **$female_audit_n** canonical mammary, oocyte, meiotic, Müllerian, and placental genes |
| [Cross-species CTA orthologs (`cross-species-orthologs.csv`)](antigen-whitepaper/cross-species-orthologs.csv) | Ortholog mapping and native catalog calls across **$ortholog_total_keys_n** mammalian CTA families in human (`9606`), dog (`9615`), and mouse (`10090`) |
| [Normal RNA coverage](audits/target-selection-20261110/rebuild/normal-rna/index.md) · [Tumor coverage ledger](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv) · [Literature adjudication](audits/target-selection-20261110/rebuild/literature.md) | Audit ledgers for healthy panels, tumor cohorts, and literature reviews |
