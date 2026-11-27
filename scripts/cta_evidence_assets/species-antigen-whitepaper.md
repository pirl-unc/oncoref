# Multi-species antigen whitepaper: human, canine, and murine parameterization

*[Interactive antigen explorer](antigen-explorer.html) · [Main antigen whitepaper](antigen-whitepaper.md) · [Species & data groups API](species.md) · [Dataset bindings (CSV)](antigen-whitepaper/species-datasets.csv) · [Non-human nominations (CSV)](antigen-whitepaper/species-nominations.csv) · [Cross-species orthologs (CSV)](antigen-whitepaper/cross-species-orthologs.csv)*

<p class="section-nav"><strong>Jump to:</strong> <a href="#goals-and-species-invariant-architecture">Goals &amp; architecture</a> · <a href="#ectopic-activation-in-vivo-purity-and-chemosensory-gpcr-guards">Ectopic, TILs &amp; GPCRs</a> · <a href="#comparative-anatomy-and-sex-restricted-subsets">Anatomy &amp; sex subsets</a> · <a href="#eight-pillar-cross-species-parameterization-matrix">8-pillar matrix</a> · <a href="#homo-sapiens-human-taxon-9606">Human (9606)</a> · <a href="#canis-lupus-familiaris-domestic-dog-taxon-9615">Dog (9615)</a> · <a href="#mus-musculus-laboratory-mouse-taxon-10090">Mouse (10090)</a> · <a href="#cross-species-cta-orthology-across-human-dog-and-mouse">Cross-species orthology</a> · <a href="#programmatic-and-cli-access">API &amp; CLI</a></p>

## Abstract

Translational cancer immunotherapy moves vaccine constructs, T-cell receptor designs, and preclinical safety benchmarks across humans, outbred companion dogs with spontaneous cancers, and syngeneic or genetically engineered mouse models. Treating "cancer–testis antigen (`CTA`)" as a human gene list projected onto animals by symbol orthology fails in both directions: `CT-X` and zygotic-genome-activation (`ZGA`) gene arrays evolve rapidly across mammals (`NY-ESO-1`, `SPANX`, `VCX/Y`, `GAGE`, `XAGE`, `DUX4`, `PRAMEF` vs. murine `Mageb`, `Pramel`, `Sbp`, `Pbsn`), accessory reproductive and disposable lymphoid organs differ anatomically across species, and heterogeneous assay scales (`TPM`, `nTPM`, `TMM-CPM`, droplet `UMI`, `IHC`, pooled `TMT-iBAQ`) cannot be merged or converted by rule of thumb.

OncoRef solves comparative antigen discovery by separating **species-invariant mathematical definitions** from **per-species genome, anatomical, dataset, and threshold parameterizations** (`oncoref.species_policies`, `oncoref.species_target_selection`). Across **human (`9606`, `$n` PEGs, `$cta` CTAs)**, **domestic dog (`9615`, `$dog_n` PEGs, `$dog_cta` CTAs)**, and **laboratory mouse (`10090`, `$mouse_n` PEGs, `$mouse_cta` CTAs)**, **all `8/8` canonical screening pillars are active in every species** with zero cross-species label borrowing.

## Goals and species-invariant architecture

Within any mammalian species *S* with SHA-256-verified proteome reference `R_S`, comparative anatomy profile `A_S`, and declared single-panel comparison unit `u_S`, every coding locus traverses the identical four-step pipeline:

1. **Stage 0 — Identical-sequence protein-equivalence grouping.** Qualified coding translations within `R_S` collapse into identical-sequence **protein-equivalence groups**. Within each expression library, complete member-locus abundance is summed before summarizing across donors, while every annotated isoform translation is retained for **≥95% 9-mer containment** against that species' <span class="cat cat-unsafe">Unsafe</span> proteins and downstream host self-proteome subtraction.
2. **Stage 1 — Multi-modal normal somatic & activated-immune safety (`Unsafe` precedence).** Each normal-tissue dataset is evaluated in its own native unit and assay rule—bulk somatic and sorted/activated immune RNA ceilings (`Vital ≤ 4`, `Protected ≤ 8`, `Manageable ≤ 16` hold / `> 128` <span class="cat cat-unsafe">Unsafe</span>, after uniform `≤ 2.5%` ectopic olfactory-receptor `≥ 5 TPM` genomic-DNA library QC), donor-resolved rare-cell `UMI` screens (`≥10` cells and `≥1%` multi-molecule cells across `≥2` independent studies or `≥5%` in one study; `2×` in `Protected`), normal protein staining or quantitative mass-spectrometry proteomics, and curated adult normal-organ literature. Confirmed Vital/Protected adult somatic expression permanently assigns <span class="cat cat-unsafe">Unsafe</span>.
3. **Stage 2 — Sequence, unscreened-niche & monogenic chemosensory GPCR eligibility (`Other`).** Surviving proteins with incomplete required normal coverage or unresolved paralog/host attribution enter `Other:Review`; proteins nested `≥95%` in `Unsafe` 9-mers or readthrough loci enter `Other:Sequence`; and proteins residing in unscreened adult micro-niches—including **all unvalidated monogenic chemosensory GPCRs (`OR*`/`COR*`/`Olfr*`, `TAS1R*`, `TAS2R*`, `VN1R*`/`Vmn*r*`)**—carry `home_tissue_unmeasured = True`, blocking <span class="cat cat-ca">CA</span>/<span class="cat cat-cta">CTA</span> admission (`Other:Unscreened tissue` or `Reproductive Antigen only`).
4. **Stage 3 — Symmetric `10 × H₀` *in vivo* enrichment, `tumor_activation_mode`, and *in vitro* cell-line references.** Across non-thymic, non-reproductive somatic baseline `N` (within-organ donor P95), species-applicable non-thymic `Disposable` maximum `D` (capped at `≤ 32`), and detection floor `F_S = 1` in unit `u_S`, let **H₀ = max(N, D, F_S)**:
   - <span class="cat cat-ca">Cancer Antigen (CA)</span>: clears shared eligibility, has `home_tissue_unmeasured == False`, and reaches **≥ 10 × H₀** in at least one QC-passing ***in vivo*** patient, syngeneic/orthotopic graft, `GEMM`, or carcinogen-induced tumor sample of species *S*. Pure *in vitro* cancer cell lines (`comparison_cell_lines`) never grant `CA` or `CTA` on their own; instead they populate `cell_line_max_tpm`, `cell_line_top_model`, and `qualifying_cell_lines`.
   - <span class="cat cat-ra">Reproductive Antigen (RA)</span>: clears shared eligibility and reaches **≥ 10 × H₀** in an adult reproductive sample or cell-supported embryonic lineage of species *S*.
   - <span class="cat cat-cta">Cancer–Testis Antigen (CTA)</span>: the intersection **CTA = CA ∩ RA**, annotated with `tumor_activation_mode` (`"ectopic"`, `"ectopic_and_lineage_retained"`, or `"lineage_retained"`).

![Conceptual three-stage antigen screening architecture and three-species proteome attrition](antigen-whitepaper/conceptual-pipeline-architecture.png)

*Figure 1.* **(A)** Species-invariant three-stage screening and downstream `HostContext` architecture. **(B)** Empirical proteome-wide attrition across human (`$n` PEGs), domestic dog (`$dog_n` PEGs), and laboratory mouse (`$mouse_n` PEGs), all evaluated across **8/8** active evidence pillars.

> [!NOTE]
> Negative-IHC RNA ceiling multipliers (`2×` Approved, `4×` Supported, `8×` Enhanced) apply **only** where a species has routine-context, cell-type-resolved antibody validation (`human` HPA v23). Pooled mass-spectrometry proteomics (`dog` Pfizer `TMT-iBAQ`) and literature-curated developmental/adult staining (`mouse` MGI `GXD` `TS28`) supply affirmative <span class="cat cat-unsafe">Unsafe</span> / review-hold gates when positive, and never inflate RNA ceilings when absent.

## Ectopic activation, *in vivo* purity, and chemosensory GPCR guards

Harmonizing target discovery across human, dog, and mouse uncovered three cross-cutting methodological requirements that apply identically across species:

1. **Tracking `tumor_activation_mode` without discarding homologous-organ targets:** Across all three species, **75%–97%** of biological CTAs exhibit true **ectopic** activation in non-homologous somatic cancers (**$cta_ectopic_any_n/$cta = $cta_ectopic_any_pct%** in human; **$dog_cta_ectopic_any_n/$dog_cta = $dog_cta_ectopic_any_pct%** in dog; **$mouse_cta_ectopic_any_n/$mouse_cta = $mouse_cta_ectopic_any_pct%** in mouse). Retaining homologous-organ witnesses (`SEMG2` in human `PRAD`, `CSN2`/`CSN3` in canine `CMT`, `Pbsn`/`Sbp` in murine `BMPC_PRAD`) inside `CTA` with `tumor_activation_mode = "lineage_retained"` preserves valid clinical/preclinical targets for hosts where the homologous organ is absent or expendable (`HostContext`), while letting pan-cancer callers filter directly on `tumor_activation_mode != "lineage_retained"`.
2. **Separating *in vivo* tumors from pure *in vitro* cancer cell lines diagnoses host `TIL` contamination:** *In vivo* tumors mix malignant cells with host leukocytes, whereas pure *in vitro* cell lines lack host stroma but acquire culture-specific drift. Requiring *in vivo* tumor witnesses for `CA`/`CTA` while evaluating **$dog_cell_lines_n** canine and **$mouse_cell_lines_n** murine pure cancer cell lines separately as an experimental model reference immediately unmasks host `LAK`/`NK` cytolytic granzymes (`Gzmd`, `Gzme`, `Gzmf`, `Gzmg`), which reach `443–2,138 TPM` inside *in vivo* `MC38`/`CT26`/`B16F10` grafts (via infiltrating host NK cells) yet collapse to `0.05–2.26 TPM` in pure cultured tumor cells—confirming their exclusion as activated-NK `Unsafe` host transcripts while simultaneously identifying which genuine CTAs (`Hmga2`, `Tex19.1`, `Wfdc13`, `Pramel13`, canine `HORMAD1`) have ready-to-use syngeneic or canine *in vitro* target lines (`qualifying_cell_lines`).
3. **Uniform monogenic chemosensory GPCR hold (`0 CA / 0 CTA`):** Singular-receptor choice in `~0.1%` of olfactory sensory neurons and taste/solitary chemosensory cells dilutes `OR*`/`COR*`/`Olfr*`, `TAS1R*`, `TAS2R*`, and `VN1R*`/`Vmn*r*` `~1,000×` in whole-mucosa RNA. Enforcing `home_tissue_unmeasured = True` across all three species prevents `29` human `CA only` GPCRs and `3` human `CTA` ORs (as well as canine `COR*` and murine `Olfr*`/`Vmn*r*` loci) from entering `CA` or `CTA`.

![Ectopic versus homologous lineage activation, in-vivo versus pure cell-line TIL separation, and chemosensory GPCR guard](antigen-whitepaper/ectopic-vs-lineage-and-cell-lines.png)

*Figure 2.* **(A)** `CTA` and `CA only` breakdown by `tumor_activation_mode` across human, dog, and mouse. **(B)** *In vivo* murine tumor TPM versus pure *in vitro* cell-line TPM separating cancer-intrinsic CTAs from excluded host `TIL` NK granzymes (`Gzmd–g`). **(C)** Monogenic chemosensory GPCR dispositions across human (`$human_gpcr_n`), dog (`$dog_gpcr_n`), and mouse (`$mouse_gpcr_n`), yielding `0 CA / 0 CTA` in all three species.

## Comparative anatomy and sex-restricted subsets

Every <span class="cat cat-ra">RA</span> and <span class="cat cat-cta">CTA</span> in every species belongs to **one of four mutually exclusive host-physiology classes**, defined by the **same `HostContext` protection algebra** parameterized by that species' actual reproductive, mammary, thymic, and disposable organs in `normal-tissue-taxonomy.yaml`:

- **`RA-M` (`CTA-M` / `RA-M-only` · Male-restricted, non-pregnant female-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every non-pregnant, non-lactating female adult reproductive and resting-mammary compartment of species *S* (`< 8` native RNA units—`< 16` in `human` where routine IHC is `not_detected`—with measured coverage and resolved protein status), while expressed in male organs.
- **`RA-F` (`CTA-F` / `RA-F-only` · Female-restricted, male-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every male adult reproductive compartment (plus resting mammary tissue in species where adult males retain mammary ducts, `human` and `dog`; in `mouse`, fetal androgen regression at E14.5 ablates male mammary buds), while expressed in non-pregnant female reproductive or resting-mammary compartments.
- **`RA-P` (`CTA-P` / `RA-P-only` · Pregnancy / peripartum-only, adult-absent across non-pregnant, non-lactating men & women):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible under full non-pregnant, non-lactating adult organ protection in **both** sexes—expressed in `embryo`, `placenta`, or `lactating_mammary` (`lactating_breast`) and silent across every adult male and non-pregnant, non-lactating female reproductive and resting-mammary compartment. For male/female CTAs that are additionally active in `placenta` or `lactating_mammary` (`PLAC1`, `XAGE3`, `LALBA`, `CSN2`), `HostContext(pregnant=True)` / `HostContext(lactating=True)` and `peripartum_antigen_ids()` enforce organ-level exclusions.
- **`RA-MF` (`CTA-MF` / `RA-MF-only` · Shared male & female):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) expressed in adult reproductive/resting-mammary organs of both sexes (or lacking measured opposite-sex clearance, failing closed).

$comparative_anatomy_table

Comparative anatomy directly shapes both enrichment denominators and host protection:

- **Accessory male glands:** Dogs anatomically lack seminal vesicles and bulbourethral glands (`prostate` is the sole major accessory gland), whereas mice carry a distinct `coagulating_gland` (anterior prostate) alongside `seminal vesicle` and `prostate`. Requesting `seminal vesicle` in a canine `HostContext` raises an immediate `ValueError` rather than silently dropping the constraint.
- **Bicornuate female tract & gonadectomy:** Both dogs and mice add `uterine_horn` (normalized automatically under `uterus`), and companion-animal oncology frequently treats spayed (ovariohysterectomized) or castrated hosts where `HostContext(absent=("ovary", "uterus", "uterine_horn", "endometrium"), species="dog")` explicitly records surgically absent organs distinct from accepted tissue injury.
- **Disposable somatic organs (`D`):** Neither dogs nor mice possess a vermiform `appendix`, and mice also lack palatine `tonsil`, narrowing non-thymic `Disposable` (`D`) from five compartments in human to four in dog and three in mouse.

![Comparative canine and murine in-vivo cohorts, cell lines, and four-way sex restriction across species](antigen-whitepaper/comparative-species-funnel-and-cohorts.png)

*Figure 3.* **(A)** Spontaneous companion-dog patient tumor cohorts (**$dog_cohorts_n** cohorts, **$dog_tumor_samples_n** primary patient samples, plus **$dog_cell_lines_n** canine cell lines). **(B)** Murine *in vivo* cohorts (**$mouse_cohorts_n** cohorts, **$mouse_tumor_samples_n** *in vivo* tumors across `H-2b`, `H-2d`, and `H-2q`, plus **$mouse_cell_lines_n** murine cell lines). **(C)** 4-way physiological sex-restriction partition (`CTA-M`, `CTA-P`, `CTA-F`, `CTA-MF`) across human, dog, and mouse.

## Eight-pillar cross-species parameterization matrix

All **8/8** canonical screening pillars are active across human (`9606`), domestic dog (`9615`), and laboratory mouse (`10090`) ([download full binding table](antigen-whitepaper/species-datasets.csv)):

$eight_pillar_matrix

## *Homo sapiens* (human · taxon 9606)

Human target selection is evaluated across **$n** protein-equivalence groups (**$loci** Ensembl 115 coding loci) under policy `$policy`: **$cta** <span class="cat cat-cta">CTA</span> (`$cta_ectopic_any_n/$cta = $cta_ectopic_any_pct%` with ectopic activation, `$cta_lineage_only_n` homologous lineage-retained), **$ca_only** <span class="cat cat-ca">CA only</span>, **$ra_only** <span class="cat cat-ra">RA only</span> (**$union** total candidates), **$other_n** <span class="cat cat-other">Other</span>, and **$unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Stratifying GTEx v11 female reproductive and mammary organs to reproductive-age `XX` women (`ages 20–49`, `gtex.donor_age()`) and integrating stage-replicated adult human folliculogenesis/oocyte RNA-seq (`GSE107746`, `128` QC-passing libraries across `9` cell-type × stage strata) partitions the **$ra** Reproductive Antigens (**$cta** CTAs) into **$ra_m_n** `RA-M` (**$cta_m_n** `CTA-M`), **$ra_p_n** `RA-P` (**$cta_p_n** `CTA-P`), **$ra_f_n** `RA-F` (**$cta_f_n** `CTA-F`), and **$ra_mf_n** `RA-MF` (**$cta_mf_n** `CTA-MF`). See the [main antigen whitepaper](antigen-whitepaper.md) for clinical TCR gold-standard benchmarks, human figures, and edge-case audits.

$human_datasets_table

## *Canis lupus familiaris* (domestic dog · taxon 9615)

Spontaneous canine cancers—appendicular osteosarcoma, multicentric B- and T-cell lymphoma, visceral hemangiosarcoma, mammary carcinoma, soft-tissue sarcoma, invasive urothelial carcinoma, oral squamous cell carcinoma, acanthomatous ameloblastoma, histiocytic sarcoma, glioma, meningioma, adrenocortical & thyroid carcinoma, and acute leukemia—develop in outbred, immunocompetent companion dogs sharing human environmental exposures. Canine target selection is evaluated genome-wide across **$dog_n** Ensembl 98 `CanFam3.1` (`GCA_000002285.2`, cross-referenced to NCBI `RefSeq 106` `UU_Cfam_GSD_1.0`) protein-equivalence groups (**$dog_loci** coding loci, **8/8** evidence pillars active): **$dog_cta** <span class="cat cat-cta">CTA</span> (`$dog_cta_ectopic_any_n/$dog_cta = $dog_cta_ectopic_any_pct%` with ectopic activation; `$dog_cta_in_cell_lines_n` also ≥10× in pure canine cell lines), **$dog_ca_only** <span class="cat cat-ca">CA only</span>, **$dog_ra_only** <span class="cat cat-ra">RA only</span> (**$dog_union** total candidates), **$dog_other_n** <span class="cat cat-other">Other</span>, and **$dog_unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Across the **$dog_ra** canine Reproductive Antigens (**$dog_cta** CTAs), comparative canine `HostContext` partitioning yields **$dog_mra_n** `RA-M` (**$dog_cta_m_n** `CTA-M`), **$dog_era_n** `RA-P` (**$dog_cta_p_n** `CTA-P`), **$dog_fra_n** `RA-F` (**$dog_cta_f_n** `CTA-F`), and **$dog_sra_n** `RA-MF` (**$dog_cta_mf_n** `CTA-MF`).

### Canine dataset and threshold parameterization

$dog_datasets_table

Key canine methodological rules enforced in code (`oncoref.canine`, `oncoref.species_target_selection`, `scripts/build_dog_target_catalog.py`):

- **Single-scale linear `TPM` somatic & macrophage panels + Cufflinks paralog guard:** `Bgee 15.2` (`25` healthy canine RNA-seq experiments deduplicated by SRA Run ID, routing BarkBase gestational `d33–d44` libraries in `SRP114662` to `embryo` and excluding neonatal/puppy stages), the **300-sample** sex-stratified Preclinical Species Gene Expression Database (`GSE219045`, normalized within-sample from `FPKM` to exact `TPM`), **64** normal female mammary biopsies (`GSE119810`, split by coordinated 4-gene milk-operon induction into `lactating_mammary` vs. `resting_mammary`, keeping unmeasured male canine `resting_mammary` fail-closed), primary companion-dog adrenal, blood, bone marrow, gingiva, brain, liver, lung, and spleen controls (`GSE294265`, `GSE294391`, `GSE272909`, `GSE252794`, `GSE288068`), and `GSE283373` primary BAL alveolar macrophages cover all **24** required (`36` evaluated) somatic contexts plus `13` RNA-measured reproductive/embryonic compartments (`+1` `proteome_only` `fallopian tube`) under uniform `≤ 2.5%` ectopic olfactory-receptor (`COR*`/`OR*`) gDNA QC. BarkBase `TMM-CPM` (`46` shared SRA runs) is never pooled with `TPM`, Cufflinks multi-mapping `FPKM` cohorts mask phase-invariant `≥5%` 31-mer paralog-ambiguous loci, and unannotated single-exon intronless loci sharing `≥50%` 9-mers with `Unsafe` somatic paralogs (`CanFam3.1` processed-retrocopy risk) are held in `Other:Review` across all cohorts.
- **Healthy companion-dog circulating leukocyte 10x scRNA-seq (`GSE225599`, closing the 8th pillar) & activated-TIL gates:** Donor-resolved 10x Chromium UMI matrices across **7 healthy companion dogs (`H_1..H_7`)** (`Ammons et al. 2023`, `CanFam3.1`) resolve `T/NK lymphocyte`, `monocyte`, and `neutrophil` (`Protected`) plus `platelet` (`Vital`) per donor via `assess_rare_cell_rna` (`150,155` donor × cell-type evaluations, excluding `OS_1..OS_10` tumor-bearing dogs and `Manageable` B cells), paired with curated canine activated-TIL checkpoint/cytokine exclusions (`PDCD1`, `TNFRSF4`, `TNFSF4`, `IFNG`, `IFNB1`, `IL17A`, `IL17F`, `CD70`, `TNFRSF8`, `SIGLEC15`) so tumor-infiltrating T-cell transcripts in lymphoma/melanoma biopsies never masquerade as tumor antigens.
- **Pooled 38-tissue proteomics with within-plex TMT bleed subtraction (`pfizer-dog-protein`):** MassIVE `MSV000099744` (`9,929` canine protein groups converted from `rIBAQ` to within-channel molar `ppm`, subtracting `10%` same-plex reproductive reporter-ion bleed) supplies affirmative <span class="cat cat-unsafe">Unsafe</span> (`Vital > 25 ppm`, `Protected > 50 ppm`) and review-hold (`Vital 10–25 ppm`, `Protected 20–50 ppm`) gates with **zero** negative-IHC RNA ceiling multipliers. `fallopian tube` (`Pool3` `126`) carries `coverage = "proteome_only"` and fails closed (`unresolved`) when directly protected.

### Spontaneous canine patient tumor cohorts ($dog_cohorts_n *in vivo* cohorts, $dog_tumor_samples_n samples)

$dog_cohorts_table

### Pure *in vitro* canine cancer cell-line purity & model-selection panel ($dog_cell_lines_n lines)

Pure *in vitro* canine lines (`canine.cell_lines()`) are strictly segregated from *in vivo* `CA`/`CTA` admission and populate `cell_line_max_tpm`, `cell_line_top_model`, and `qualifying_cell_lines`:

$dog_cell_lines_table

### Top canine CTAs by peak *in vivo* tumor expression

$dog_top_cta_table

### Packaged canine nomination ledger (14 literature RT-PCR candidates)

$dog_nominations_table

Evaluating all **14** published canine RT-PCR nominations (`hindriks2025`, `nemec2019`) against genome-wide canine somatic RNA, `GSE225599` healthy PBMC 10x scRNA-seq, and 38-tissue proteomics shows why quantitative somatic gates are essential: **`PRAME` (`$dog_prame_call`, `CTA-M`)** and **`SSX6` (`CTA-M`)** are validated as genuine canine CTAs across spontaneous patient cohorts, `MAGEB10` qualifies as `Reproductive Antigen only`, `MAGEA13P` (`Other:Expression`) and ambiguous `MAGE` (`Other:Expression / Reproductive Antigen only`) stay below `CA`, `CAGE1`/`CTCFL`/`DDX53` remain held in `Other:Review`, and somatic-expressed nominees (`AKR1E2`, `IGF2BP3`, `NT5C1B`, `PAGE3-like`, `SPECC1`, `TPX2`) are excluded as `Unsafe` due to healthy canine blood, bone marrow, lung, brain, or liver expression.

## *Mus musculus* (laboratory mouse · taxon 10090)

Preclinical vaccine, TCR, and bispecific studies in mice require host-proteome peptide subtraction and endogenous normal-tissue screening on mouse-native coordinates rather than human-to-mouse symbol projection—especially across X-linked CTA families (`Magea`, `Mageb`, `Xlr`, `Trap1a`, `Pramef`) where copy number and tissue regulation diverge sharply from primates.

Murine target selection is evaluated genome-wide across **$mouse_n** Ensembl 115 `GRCm39` protein-equivalence groups (**$mouse_loci** primary-assembly coding loci, **8/8** evidence pillars active under harmonized `≤ 2.5%` ectopic `Olfr` gDNA library QC and monogenic chemosensory GPCR holds): **$mouse_cta** <span class="cat cat-cta">CTA</span> (`$mouse_cta_ectopic_any_n/$mouse_cta = $mouse_cta_ectopic_any_pct%` with ectopic activation; `$mouse_cta_in_cell_lines_n` also ≥10× in pure murine cancer cell lines), **$mouse_ca_only** <span class="cat cat-ca">CA only</span>, **$mouse_ra_only** <span class="cat cat-ra">RA only</span> (**$mouse_union** total candidates), **$mouse_other_n** <span class="cat cat-other">Other</span>, and **$mouse_unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Across the **$mouse_ra** murine Reproductive Antigens (**$mouse_cta** CTAs), 16-compartment `HostContext` partitioning yields **$mouse_mra_n** `RA-M` (**$mouse_cta_m_n** `CTA-M`), **$mouse_era_n** `RA-P` (**$mouse_cta_p_n** `CTA-P`), **$mouse_fra_n** `RA-F` (**$mouse_cta_f_n** `CTA-F`), and **$mouse_sra_n** `RA-MF` (**$mouse_cta_mf_n** `CTA-MF`).

### Murine dataset and threshold parameterization

$mouse_datasets_table

Key murine methodological rules enforced in code (`oncoref.mouse`, `oncoref.species_target_selection`, `scripts/build_mouse_target_catalog.py`):

- **Complete pinned Ensembl 115 GRCm39 host proteome:** `mouse.references.ensure_reference("ensembl115")` verifies SHA-256-pinned GTF, cDNA, ncRNA, and protein FASTAs (**66,037** qualified primary-assembly coding translations across **21,433** `PEG`s and **21,760** coding loci), enabling full murine `PEG` grouping, 9-mer <span class="cat cat-unsafe">Unsafe</span> containment, and offline `PeptideIndex` self-proteome subtraction (`mouse-proteome-sequences`).
- ***In vivo* tumor cohorts vs. pure *in vitro* murine cancer cell lines:** Reference coordinates belong to **C57BL/6J** (`GRCm39`), while **$mouse_cohorts_n** *in vivo* cohorts (`$mouse_tumor_samples_n` samples) span **C57BL/6 (`H-2b`)** (`MC38`, `B16F10`, `B78_MEL`, `LLC`, `RM1`, `KP_LUAD`, `GL261_TUMOR`, `CT2A_TUMOR`, `SB28_TUMOR`, `MADR_GLIOMA`, `DEN_HCC`), **BALB/c (`H-2d`)** (`CT26`, `EMT6`), and **FVB (`H-2q`)** (`BMPC_PRAD`), paired with **$mouse_cell_lines_n** pure *in vitro* cancer cell lines (`4T1`, `P815`, `EO771`, `Hepa1-6`, `ID8`, `GL261_LINE`, `MB49`, `MOC1`, `MOC22` in `mouse.cell_lines()`). Importantly, whereas expanding outbred spontaneous companion-dog cohorts (`7 → 19` cohorts, `935` patients) more than doubled canine CTAs (`54 → $dog_cta`), expanding *in vivo* murine models (`6 → 14` cohorts) while stripping out `4` host `TIL` NK granzymes (`Gzmd–g`) and evaluating pure *in vitro* `4T1`/`P815` separately in `mouse.cell_lines()` yields **$mouse_cta** *in vivo* murine CTAs—reflecting both the absence of primate `CT-X` gene arrays in rodents and the rapid single-oncogene kinetics of inbred `GEMM`/syngeneic models compared with multi-year epigenetic hypomethylation in spontaneous human and canine cancers.
- **Multi-modal normal & cytokine-activated NK/Treg/macrophage safety without negative-IHC inflation:** Adult somatic bulk RNA (`Bgee 15.2` + `GSE219045` + `GSE322844` untreated WT liver), sorted & activated immune populations (`ImmGen` ULI `GSE109125`, `ImmGen` MNP `GSE122108`, WT `LPS`/`IL-4` primary BMDMs `GSE181223`, splenic/IL-15 NK `GSE287306`, IL-2/ionomycin Treg `GSE336796`, splenic CD4/B `GSE302070`, and curated `IL-2`/`IL-12`/`IL-18` `LAK`/`NK` cytolytic granzymes `Gzmd–g`), `Tabula Muris` 10x droplet UMI (`assess_rare_cell_rna`), and `MGI GXD 6.27` adult `TS28` immunohistochemistry/immunoblot detections gate safety in their native units with zero negative-IHC RNA ceiling multipliers, while mapped-universe zero-count completion resolves all 16 murine `HostContext` compartments.

### Murine *in vivo* syngeneic, GEMM, and carcinogen tumor cohorts ($mouse_cohorts_n cohorts, $mouse_tumor_samples_n samples)

$mouse_cohorts_table

### Pure *in vitro* murine cancer cell-line purity & model-selection panel ($mouse_cell_lines_n lines)

$mouse_cell_lines_table

### Top murine CTAs by peak *in vivo* tumor expression

$mouse_top_cta_table

### Historical murine nomination check (`P1A` / `Trap1a`)

$mouse_nominations_table

`P1A` (`Trap1a`, `ENSMUSG00000051257` / `ENSMUST00000054534.4` / `ENSMUSP00000059448.4`, 224 aa) illustrates why OncoRef separates nomination provenance and pure cell-line expression from *in vivo* safety admission: the DBA/2 `P815` rejection antigen (`M36387.1` / `AAA40475.1`, Van den Eynde et al. 1991) encodes an identical 224-aa protein to C57BL/6J `ENSMUSP00000059448.4` and reaches high expression in pure `P815`, `EO771`, `Hepa1-6`, and `MOC22` cancer lines (`cell_line_max_tpm`), while endogenous somatic counterevidence in BALB/c lung, spleen, liver, and thymus (Sarma et al. 1999) plus quantitative multi-organ screening classify `Trap1a` as `$mouse_p1a_call` rather than admitting it to `CTA`.

## Cross-species CTA orthology across human, dog, and mouse

With independent genome-wide catalogs evaluated natively across **human (`$cta` CTAs)**, **dog (`$dog_cta` CTAs)**, and **mouse (`$mouse_cta` CTAs)**, OncoRef quantifies directly how often a `CTA` in one mammal maps to a `CTA` in another ([download full `$ortholog_total_keys_n`-family ortholog concordance table](antigen-whitepaper/cross-species-orthologs.csv)):

$ortholog_summary_table

Why do only **$conserved_cta_pairs_n** of **$ortholog_total_keys_n** mammalian `CTA` ortholog families independently qualify as `CTA` in two or more species (and **$conserved_cta_or_ra_n** as conserved `CTA` or `RA-only`)? Three evolutionary and experimental mechanisms govern cross-species mapping:

1. **Rapid mammalian lineage amplification of CT-X and ZGA arrays:** Nearly half of human CTAs—including `NY-ESO-1` (`CTAG1A/B`), `SPANXA/B/C/D/N`, `VCX/VCY`, `GAGE*`, `XAGE1A/B`, `CT45A*`, `CT47A*`, `DUX4` (`D4Z4`), `PRAMEF1–26`, `MBD3L2–5`, and `SSU72L1–5`—originated or amplified in primates after boreoeutherian divergence and have **no 1:1 coding ortholog** in dog or mouse (`human_cta_with_dog_orth_n = $human_cta_with_dog_orth_n / $cta`, `human_cta_with_mouse_orth_n = $human_cta_with_mouse_orth_n / $cta`), while rodents independently expanded `Pramel*`, `Mageb*`, and `Sbp`/`Pbsn` (`CTA-M` in mouse).
2. **Conserved boreoeutherian germline/lactation core (`$conserved_cta_or_ra_n` families):** Single-copy autosomal and ancient X-linked germline, placental, and milk genes conserve reproductive restriction across 95 million years of mammalian evolution—yielding **conserved multi-species CTAs** ($conserved_cta_examples), alongside conserved `RA`s (`DPPA5`, `CIMIP4`/`TEX33`, `SPMIP9`/`TEX37`, `MEIOSIN`/`BHMG1`, `TBXT`/`T`, `DEFB123`/`Defb23`, `MAGEB18`, `CYLC1`, `ACTL8`, `ADAM2`, `AKAP4`) whose tumor reactivation tracks histology breadth (`135` human cohorts vs. `$dog_cohorts_n` canine cohorts vs. `$mouse_cohorts_n` murine models).
3. **Divergent somatic regulation and single-cell resolution:** Several conserved 1:1 orthologs (`HORMAD1`, `DDX43`, `CLDN6`, `IGF2BP1`) qualify in canine somatic panels while failing human ultra-deep single-cell tissue cluster, adult cornea, or sorted-immune gates—underscoring why target safety must be evaluated natively in every species and never inferred across species boundaries.

<details markdown="1">
<summary>Conserved multi-species CTA and RA-ortholog concordance table ($conserved_cta_or_ra_n families; download full CSV)</summary>

Download all `$ortholog_total_keys_n` ortholog families: [`cross-species-orthologs.csv`](antigen-whitepaper/cross-species-orthologs.csv).

$conserved_ortholog_ctas_table

</details>

## Programmatic and CLI access

Every species exposes the identical plain-language `oncoref.antigens` front door alongside per-species dataset, cohort, cell-line, and policy modules:

```python
from oncoref import antigens, canine, human, mouse
from oncoref.species_policies import role_comparison_table

# 1. Cross-species summary and plain-keyword queries across human, dog, and mouse
antigens.summary()
female_safe_ctas = antigens.table(
    "CTA", sex="female", ectopic_only=True
)  # ectopic subset of $cta_m_n CTA-M + $cta_p_n CTA-P
pregnancy_only = antigens.table("CTA-P")  # $cta_p_n adult-absent pregnancy/peripartum CTAs
dog_ctas = antigens.table("CTA", species="dog", ectopic_only=True, replicated_only=True)
mouse_ctas = antigens.table("CTA", species="mouse")

# 2. Spayed/castrated or organ-protected host queries and self-proteome peptide subtraction
spayed_female_dog_ctas = antigens.table(
    "CTA",
    species="dog",
    sex="female",
    absent=("ovary", "uterus", "uterine_horn", "endometrium"),
)
antigens.explain("PRAME", species="dog", sex="female")
canine_unique_9mers = antigens.subtract("PRAME", species="dog", sex="female", lengths=[9])

# 3. Inspect declarative per-species thresholds, 8/8 dataset bindings, cohorts, and cell lines
human.screening_policy()
canine.tumor_cohorts()  # $dog_cohorts_n spontaneous companion-dog tumor cohorts
canine.cell_lines()  # $dog_cell_lines_n pure canine cancer cell lines
mouse.tumor_cohorts()  # $mouse_cohorts_n in-vivo syngeneic / GEMM / carcinogen cohorts
mouse.cell_lines()  # $mouse_cell_lines_n pure murine cancer cell lines
role_comparison_table()
```

```sh
# CLI queries across human, dog, and mouse, including --sex, --ectopic-only, and --replicated
oncoref antigens summary
oncoref antigens list --sex female --ectopic-only
oncoref antigens list --set CTA-P --format json
oncoref antigens list --set CTA --species dog --replicated
oncoref antigens explain PRAME --species dog --sex female
oncoref antigens explain P1A --species mouse
```
