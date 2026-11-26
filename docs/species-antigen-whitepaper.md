# Multi-species antigen whitepaper: human, canine, and murine parameterization

*[Main antigen whitepaper](antigen-whitepaper.md) · [Species & data groups API](species.md) · [Dataset bindings (CSV)](antigen-whitepaper/species-datasets.csv) · [Non-human nominations (CSV)](antigen-whitepaper/species-nominations.csv) · [Cross-species orthologs (CSV)](antigen-whitepaper/cross-species-orthologs.csv)*

<p class="section-nav"><strong>Jump to:</strong> <a href="#one-definition-three-species-parameterizations">Universal definitions</a> · <a href="#comparative-anatomy-and-sex-restricted-subsets">Anatomy &amp; sex subsets</a> · <a href="#eight-pillar-cross-species-parameterization-matrix">8-pillar matrix</a> · <a href="#homo-sapiens-human-taxon-9606">Human (9606)</a> · <a href="#canis-lupus-familiaris-domestic-dog-taxon-9615">Dog (9615)</a> · <a href="#mus-musculus-laboratory-mouse-taxon-10090">Mouse (10090)</a> · <a href="#cross-species-cta-orthology-across-human-dog-and-mouse">Cross-species orthology</a> · <a href="#programmatic-and-cli-access">API &amp; CLI</a></p>

Translational immunotherapy moves targets, vaccine designs, and preclinical benchmarks across humans, companion dogs with spontaneous cancers, and syngeneic or genetically engineered mouse models. Treating "CTA" as a human gene list copied by orthology fails in both directions: paralog families diverge rapidly across mammals (`MAGE`, `CT45`, `PRAME`, `SPANX`, `XAGE`), accessory reproductive and lymphoid organs differ between species, and assay scales (`TPM`, `nTPM`, `TMM-CPM`, UMI counts, IHC, pooled TMT `iBAQ`) cannot be pooled or converted by rule of thumb.

OncoRef solves this by separating **species-invariant set definitions** from **per-species anatomical, dataset, and threshold parameterizations** (`oncoref.species_policies`, `oncoref.species_target_selection`). Every species runs through the same three-stage pipeline on its own genome reference and native measurements; where an individual assay layer is absent for a species (such as a normal single-cell UMI atlas in dogs), the policy records the gap explicitly rather than borrowing calls from another species.

## One definition, three species parameterizations

Within any species *S* with verified proteome reference `R_S`, comparative anatomy profile `A_S`, and declared single-panel comparison unit `u_S`:

1. **Protein-equivalence counting unit (`PEG`).** Qualified coding translations within `R_S` collapse into identical-sequence **protein-equivalence groups** (`PEG:<sha256[:20]>`). Within each expression library, complete member-locus abundance is summed before summarizing across donors, while every annotated isoform translation is retained for **≥95% 9-mer containment** against that species' <span class="cat cat-unsafe">Unsafe</span> proteins and downstream `PeptideIndex` subtraction.
2. **Stage 1 — Multi-modal normal safety and shared eligibility.** Each normal-tissue dataset is evaluated in its own native unit and assay rule—bulk somatic and immune RNA ceilings (`Vital ≤ 4`, `Protected ≤ 8`, `Manageable ≤ 16` hold / `> 128` <span class="cat cat-unsafe">Unsafe</span>, after uniform `≤ 2.5%` ectopic olfactory-receptor `≥ 5 TPM` genomic-DNA library QC), donor-resolved rare-cell UMI screens (`≥10` cells and `≥1%` multi-molecule cells across `≥2` independent studies or `≥5%` in one study; `2×` in `Protected`), normal protein staining or proteomics, and curated adult normal-organ literature. Confirmed Vital/Protected adult somatic expression always takes precedence (<span class="cat cat-unsafe">Unsafe</span>); incomplete required normal coverage, unresolved paralog/host attribution, or unscreened adult home organs stay in <span class="cat cat-other">Other</span> (`Other:Review`, `Other:Sequence`, `Other:Unscreened tissue`).
3. **Stage 2 — Symmetric enrichment (`CA`, `RA`, `CTA`).** Across non-thymic, non-reproductive somatic baseline `N` (within-organ donor P95), species-applicable non-thymic `Disposable` maximum `D` (capped at `≤ 32`), and detection floor `F_S = 1` in unit `u_S`, let **H₀ = max(N, D, F_S)**:
   - <span class="cat cat-ca">Cancer Antigen (CA)</span>: clears shared eligibility, has no unscreened adult home organ, and reaches **≥ 10 × H₀** in at least one QC-passing tumor sample of species *S*.
   - <span class="cat cat-ra">Reproductive Antigen (RA)</span>: clears shared eligibility and reaches **≥ 10 × H₀** in an adult reproductive sample or cell-supported embryonic lineage of species *S*.
   - <span class="cat cat-cta">Cancer–Testis Antigen (CTA)</span>: the intersection **CTA = CA ∩ RA**, partitioning the candidate union **CA ∪ RA** into `CTA`, `Cancer Antigen only`, and `Reproductive Antigen only`.

> Negative-IHC RNA ceiling multipliers (`2×` Approved, `4×` Supported, `8×` Enhanced) apply **only** where a species has routine-context, cell-type-resolved antibody validation (`human` HPA v23). Pooled mass-spectrometry proteomics (`dog` Pfizer TMT-iBAQ) and literature-curated staining (`mouse` GXD) supply affirmative <span class="cat cat-unsafe">Unsafe</span> / hold evidence when positive, and never inflate RNA ceilings when absent.

## Comparative anatomy and sex-restricted subsets

Every <span class="cat cat-ra">RA</span> and <span class="cat cat-cta">CTA</span> in every species belongs to **one of four mutually exclusive host-physiology classes**, defined by the **same `HostContext` protection algebra** parameterized by that species' actual reproductive, mammary, thymic, and disposable organs in `normal-tissue-taxonomy.yaml`:

- **`RA-M` (`CTA-M` / `RA-M-only`, legacy alias `MRA` · Male-restricted, non-pregnant female-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every non-pregnant, non-lactating female adult reproductive and resting-mammary compartment of species *S* (`< 8` native RNA units—`< 16` in `human` where routine IHC is `not_detected`—with measured coverage and resolved protein status), while expressed in male organs.
- **`RA-F` (`CTA-F` / `RA-F-only`, legacy alias `FRA` · Female-restricted, male-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every male adult reproductive compartment (plus resting mammary tissue in species where adult males retain mammary ducts, `human` and `dog`; in `mouse`, fetal androgen regression at E14.5 ablates male mammary buds), while expressed in non-pregnant female reproductive or resting-mammary compartments.
- **`RA-P` (`CTA-P` / `RA-P-only`, legacy aliases `PRA`, `RA-E`, `CTA-E`, `ERA` · Pregnancy / peripartum-only, adult-absent across non-pregnant, non-lactating men & women):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible under full non-pregnant, non-lactating adult organ protection in **both** sexes—expressed in `embryo`, `placenta`, or `lactating_mammary` (`lactating_breast`) and silent across every adult male and non-pregnant, non-lactating female reproductive and resting-mammary compartment. For male/female CTAs that are additionally active in `placenta` or `lactating_mammary` (`PLAC1`, `XAGE3`, `LALBA`, `CSN2`), `HostContext(pregnant=True)` / `HostContext(lactating=True)` and `peripartum_antigen_ids()` enforce organ-level exclusions.
- **`RA-MF` (`CTA-MF` / `RA-MF-only`, legacy alias `SRA` · Shared male & female):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) expressed in adult reproductive/resting-mammary organs of both sexes (or lacking measured opposite-sex clearance, failing closed).

| Species | Core reproductive / embryonic | Male accessory | Female tract | Mammary | Non-thymic Disposable (D) |
| --- | --- | --- | --- | --- | --- |
| Homo sapiens (`9606`) | embryo, placenta, testis, ovary | epididymis, prostate, seminal vesicle | fallopian tube, uterus, endometrium, cervix, vagina | lactating_breast, resting_breast | appendix, gallbladder, tonsil, hair follicle, cutaneous_melanocyte |
| Canis lupus familiaris (`9615`) | embryo, placenta, testis, ovary | epididymis, prostate | fallopian tube, uterus, uterine_horn, endometrium, cervix, vagina | lactating_mammary, resting_mammary | gallbladder, tonsil, hair follicle, cutaneous_melanocyte |
| Mus musculus (`10090`) | embryo, placenta, testis, ovary | epididymis, prostate, coagulating_gland, seminal vesicle | fallopian tube, uterus, uterine_horn, endometrium, cervix, vagina | lactating_mammary, resting_mammary | gallbladder, hair follicle, cutaneous_melanocyte |

Comparative anatomy changes both enrichment denominators and host protection:

- **Accessory male glands:** Dogs anatomically lack seminal vesicles and bulbourethral glands (`prostate` is the sole major accessory gland), whereas mice carry a distinct `coagulating_gland` (anterior prostate) alongside `seminal vesicle` and `prostate`. Requesting `seminal vesicle` in a canine `HostContext` raises an immediate `ValueError` rather than silently dropping the constraint.
- **Bicornuate female tract & gonadectomy:** Both dogs and mice add `uterine_horn` (normalized automatically under `uterus`), and companion-animal oncology frequently treats spayed (ovariohysterectomized) or castrated hosts where `HostContext(absent=("ovary", "uterus", "uterine_horn", "endometrium"), species="dog")` explicitly records surgically absent organs distinct from accepted tissue injury.
- **Disposable somatic organs (`D`):** Neither dogs nor mice possess a vermiform `appendix`, and mice also lack palatine `tonsil`, narrowing non-thymic `Disposable` (`D`) from five compartments in human to four in dog and three in mouse.

## Eight-pillar cross-species parameterization matrix

The table below aligns all eight canonical screening roles across human, dog, and mouse ([download full binding table](antigen-whitepaper/species-datasets.csv)):

| Screening pillar | Human (9606) | Domestic dog (9615) | Laboratory mouse (10090) |
| --- | --- | --- | --- |
| **Reference proteome & identical-sequence protein-equivalence groups (PEG)** | ensembl115 (amino-acid sequence) — *active_in_catalog* | ensembl98 + refseq106 (amino-acid sequence) — *active_in_catalog* | ensembl115 (amino-acid sequence) — *active_in_catalog* |
| **Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D)** | gtex_v11 + hpa_v23_biopsies + sorted_immune (TPM) — *active_in_catalog* | bgee_v15_2_dog_25exps + gse219045_dog_atlas (TPM) — *active_in_catalog* | bgee_v15_2_mouse_12exps + gse219045_mouse_atlas (TPM) — *active_in_catalog* |
| **Normal somatic bulk & sorted immune RNA safety ceilings** | hpa_rna_consensus + gtex_tissue_rna + hpa_brain + dice + monaco (TPM; nTPM) — *active_in_catalog* | bgee_v15_2_dog_somatic + gse219045_dog_atlas (TPM) — *active_in_catalog* | bgee_v15_2_mouse_somatic + gse219045 + immgen_gse109125_gse122108 + gse181223_bmdm + sarma1999 (TPM; qualitative RT-PCR/Southern) — *active_in_catalog* |
| **Single-cell cluster & donor-resolved rare-cell UMI safety screens** | cellxgene_census + siletti_brain + hpa_single_cell + collin_cornea (UMI counts; nTPM; CP10K) — *active_in_catalog* | canine_normal_scrna_gap (UMI counts) — *identified_analog* | tabula_muris_10x (UMI counts) — *active_in_catalog* |
| **Normal-tissue protein staining (IHC) or mass-spectrometry proteomics** | hpa_normal_tissue + external_ihc (categorical IHC grade + reliability) — *active_in_catalog* | pfizer-dog-protein (within-plex bleed-corrected molar ppm (from rIBAQ)) — *active_in_catalog* | mgi_gxd_6_27 (curated adult IHC/Immunoblot detections) — *active_in_catalog* |
| **Adult reproductive & embryonic/placental enrichment (RA) and host exposure** | gtex_repro + hpa_repro + petropoulos + tyser + mole (TPM; UMI/staining context) — *active_in_catalog* | bgee_repro_embryo + gse219045_dog_repro + gse119810_normal_mammary + pfizer_fallopian_endometrium (TPM; molar ppm (Fallopian_Tube)) — *active_in_catalog* | bgee_repro_embryo + gse219045 + accessory_tract_lactation_7gse (TPM) — *active_in_catalog* |
| **QC-passing tumor cohorts for cancer enrichment (CA)** | oncoref_tumor_union_132_cohorts (TPM) — *active_in_catalog* | dog_8_tumor_cohorts (TPM) — *active_in_catalog* | mouse_8_cohorts_syngeneic_gemm (TPM) — *active_in_catalog* |
| **Curated normal counterevidence, unscreened adult homes & published nominations** | normal_expression_literature + rare_cell_attribution_reviews + cta_literature_149 (curated verdicts & citations) — *active_in_catalog* | hindriks2025 + nemec2019 (curated RT-qPCR & cross-reference index) — *active_in_catalog* | vandeneynde1991 + sarma1999 (curated primary-source XMLs & GenBank records) — *active_in_catalog* |

## *Homo sapiens* (human · taxon 9606)

Human target selection is fully evaluated across **19,893** protein-equivalence groups (**20,099** Ensembl 115 coding loci) under policy `oncoref.target-selection.2026-11-26.1`: **325** <span class="cat cat-cta">CTA</span>, **95** <span class="cat cat-ca">CA only</span>, **143** <span class="cat cat-ra">RA only</span> (**563** total candidates), **1,081** <span class="cat cat-other">Other</span>, and **18,249** <span class="cat cat-unsafe">Unsafe</span>. Stratifying GTEx v11 female reproductive and mammary organs to reproductive-age `XX` women (`ages 20–49`, `gtex.donor_age()`) and integrating stage-replicated adult human folliculogenesis/oocyte RNA-seq (`GSE107746`, `128` QC-passing libraries across `9` cell-type × stage strata) partitions the **468** Reproductive Antigens (**325** CTAs) into **342** `RA-M` (**224** `CTA-M`), **60** `RA-P` (**48** `CTA-P`), **18** `RA-F` (**13** `CTA-F`), and **48** `RA-MF` (**40** `CTA-MF`). See the [main antigen whitepaper](antigen-whitepaper.md) for full human figures, clinical anchor tables, and edge-case audits.

| Canonical role | Dataset(s) | Native unit | Status | Threshold / screening rule |
| --- | --- | --- | --- | --- |
| Reference proteome & identical-sequence protein-equivalence groups (PEG) | `ensembl115` — Human Ensembl 115 GRCh38 + MANE Select | amino-acid sequence | `active_in_catalog` | Prefer MANE Select then Ensembl canonical; collapse identical translations into PEG; >=95% 9-mer Unsafe containment -> Other:Sequence |
| Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D) | `gtex_v11 + hpa_v23_biopsies + sorted_immune` — GTEx v11, HPA v23 tissues/brain & HPA/DICE/Monaco immune | TPM | `active_in_catalog` | Within-site linear donor P95 sets non-thymic somatic N and Disposable D (<=32 TPM); H0 = max(N, D, 1 TPM) |
| Normal somatic bulk & sorted immune RNA safety ceilings | `hpa_rna_consensus + gtex_tissue_rna + hpa_brain + dice + monaco` — HPA v23 consensus/brain/immune, GTEx v11 medians & DICE/Monaco | TPM; nTPM | `active_in_catalog` | Base ceilings Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) in native unit; matched negative IHC raises ceiling 2x/4x/8x |
| Single-cell cluster & donor-resolved rare-cell UMI safety screens | `cellxgene_census + siletti_brain + hpa_single_cell + collin_cornea` — CELLxGENE Census, Siletti brain, HPA v23 scRNA clusters & Collin cornea | UMI counts; nTPM; CP10K | `active_in_catalog` | UMI: >=10 cells and >=1% >=2-UMI cells across >=2 consortium-deduplicated study groups (or >=5% in 1 study; 2x in Protected); cornea CP10K /100 nominal ceiling |
| Normal-tissue protein staining (IHC) or mass-spectrometry proteomics | `hpa_normal_tissue + external_ihc` — HPA v23 normal immunohistochemistry & curated external staining | categorical IHC grade + reliability | `active_in_catalog` | Approved/Supported/Enhanced positive staining in Vital/Protected tissues -> Unsafe; routine-context negative IHC grants RNA headroom |
| Adult reproductive & embryonic/placental enrichment (RA) and host exposure | `gtex_repro + hpa_repro + petropoulos + tyser + mole` — GTEx/HPA adult reproductive samples & Petropoulos/Tyser/Molè embryos | TPM; UMI/staining context | `active_in_catalog` | RA qualifies at >=10 x H0 in adult reproductive samples or cell-supported embryo lineages; per-compartment maxima feed HostContext (<8 TPM / <8 nTPM, 2x with negative IHC) |
| QC-passing tumor cohorts for cancer enrichment (CA) | `oncoref_tumor_union_132_cohorts` — 132 QC-passing per-sample human tumor cohorts (TCGA, TARGET, BEAT-AML, Treehouse, GEO) | TPM | `active_in_catalog` | CA qualifies when >=1 QC-passing library reaches >=10 x H0; known-patient counts deduplicate samples per donor |
| Curated normal counterevidence, unscreened adult homes & published nominations | `normal_expression_literature + rare_cell_attribution_reviews + cta_literature_149` — Reviewed normal-tissue literature, #711 Ensembl 115 attribution ledger & 149 CTA/oncofetal papers | curated verdicts & citations | `active_in_catalog` | Verified adult Vital/Protected expression -> Unsafe; unscreened adult home (olfactory, taste bud, bone, dental pulp, hair follicle) -> Other:Unscreened tissue |

## *Canis lupus familiaris* (domestic dog · taxon 9615)

Spontaneous canine cancers—appendicular osteosarcoma, multicentric and diffuse large B-cell lymphoma, visceral hemangiosarcoma, mammary carcinoma, and high-grade oligodendroglioma—develop in outbred, immunocompetent companion dogs sharing human environmental exposures. Canine target selection is evaluated genome-wide across **20,082** Ensembl 98 `CanFam3.1` (`GCA_000002285.2`, cross-referenced to NCBI `RefSeq 106` `UU_Cfam_GSD_1.0`) protein-equivalence groups (**20,199** coding loci, **7/8** evidence pillars active): **54** <span class="cat cat-cta">CTA</span>, **59** <span class="cat cat-ca">CA only</span>, **421** <span class="cat cat-ra">RA only</span> (**534** total candidates), **2,573** <span class="cat cat-other">Other</span>, and **16,975** <span class="cat cat-unsafe">Unsafe</span>. Across the **475** canine Reproductive Antigens (**54** CTAs), comparative canine `HostContext` partitioning yields **439** `RA-M` (**42** `CTA-M`), **8** `RA-P` (**4** `CTA-P`), **6** `RA-F` (**0** `CTA-F`), and **22** `RA-MF` (**8** `CTA-MF`).

### Canine dataset and threshold parameterization

| Canonical role | Dataset(s) | Native unit | Status | Threshold / screening rule |
| --- | --- | --- | --- | --- |
| Reference proteome & identical-sequence protein-equivalence groups (PEG) | `ensembl98 + refseq106` — Canine Ensembl 98 CanFam3.1 (GCA_000002285.2, primary) & NCBI RefSeq 106 (UU_Cfam_GSD_1.0) | amino-acid sequence | `active_in_catalog` | Complete GTF/cDNA/ncRNA/peptide FASTA pins; identical translations collapse into canine PEG:<sha256>; >=95% canine 9-mer Unsafe containment |
| Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D) | `bgee_v15_2_dog_25exps + gse219045_dog_atlas` — Bgee 15.2 adult canine multi-organ RNA-seq (SRA-deduplicated) & GSE219045 300-sample sex-stratified organ atlas | TPM | `active_in_catalog` | Within-organ donor P95 across 24 required (36 evaluated) non-thymic adult somatic contexts sets N and canine Disposable D (gallbladder, tonsil, hair follicle); <=2.5% ectopic OR gDNA QC |
| Normal somatic bulk & sorted immune RNA safety ceilings | `bgee_v15_2_dog_somatic + gse219045_dog_atlas` — Bgee 15.2 adult canine somatic/leukocyte RNA & GSE219045 multi-organ whole-blood/marrow/lymph-node atlas | TPM | `active_in_catalog` | Tiered somatic ceilings Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) TPM without negative-IHC multipliers |
| Single-cell cluster & donor-resolved rare-cell UMI safety screens | `canine_normal_scrna_gap` — Donor-resolved canine normal single-cell UMI atlases | UMI counts | `identified_analog` | Same assess_rare_cell_rna contract (>=10 cells, >=1% >=2-UMI across >=2 canine studies, 2x in Protected) when multi-organ single-cell atlases are ingested |
| Normal-tissue protein staining (IHC) or mass-spectrometry proteomics | `pfizer-dog-protein` — Pfizer 38-tissue canine quantitative proteomics (MSV000099744, 9,929 protein groups) | within-plex bleed-corrected molar ppm (from rIBAQ) | `active_in_catalog` | Sequence-reconciled Vital (>25 ppm Unsafe, 10-25 ppm hold) & Protected (>50 ppm Unsafe, 20-50 ppm hold) tissue-pool protein gates after 10% same-plex reproductive TMT bleed subtraction; never grants negative-IHC RNA multipliers |
| Adult reproductive & embryonic/placental enrichment (RA) and host exposure | `bgee_repro_embryo + gse219045_dog_repro + gse119810_normal_mammary + pfizer_fallopian_endometrium` — Bgee 15.2 testis/ovary/placenta/d33-d44 embryos, GSE219045 male/female reproductive tract, GSE119810 normal resting/lactating mammary & Pfizer Fallopian_Tube/Endometrium proteomics | TPM; molar ppm (Fallopian_Tube) | `active_in_catalog` | RA qualifies at >=10 x H0 across canine compartments (uterine_horn, absent seminal vesicle, milk-operon secretory separation) and populates canine HostContext |
| QC-passing tumor cohorts for cancer enrichment (CA) | `dog_8_tumor_cohorts` — 8 canine tumor cohorts (7 spontaneous patient cohorts + 1 DLBCL line): OSA_NCI (GSE238110), OSA_COTC030 (GSE304066), OSA_CCOGC (GSE239948), CMT (GSE119810), HSA (GSE95183), LYM (GSE95184), DLBCL (GSE287971 CLBL-1), GLIOMA (GSE284559) | TPM | `active_in_catalog` | CA qualifies when >=1 QC-passing (<=2.5% OR gDNA) canine tumor reaches >=10 x H0 on CanFam3.1 linear TPM; replicate tumor samples deduplicate per dog |
| Curated normal counterevidence, unscreened adult homes & published nominations | `hindriks2025 + nemec2019` — Curated canine CTA nominations and normal/tumor RT-qPCR observations (14 candidates) | curated RT-qPCR & cross-reference index | `active_in_catalog` | Reconciles all 14 published canine nominations onto CanFam3.1 while enforcing somatic RT-PCR counterevidence holds and zero cross-species CTA transfer |

Key canine methodological rules enforced in code (`oncoref.canine`, `oncoref.species_target_selection`, `scripts/build_dog_target_catalog.py`):

- **Single-scale linear `TPM` panels & Cufflinks paralog guard:** `Bgee 15.2` (`25` healthy canine RNA-seq experiments deduplicated by SRA Run ID, routing BarkBase gestational `d33–d44` libraries in `SRP114662` to `embryo` and excluding neonatal/puppy stages), the **300-sample** sex-stratified Preclinical Species Gene Expression Database (`GSE219045`, normalized within-sample from `FPKM` to exact `TPM`), and **64** normal female mammary biopsies (`GSE119810`, split by coordinated 4-gene milk-operon induction into `lactating_mammary` vs. `resting_mammary`, keeping unmeasured male canine `resting_mammary` fail-closed) cover all **24** required (`36` evaluated) somatic contexts plus `13` RNA-measured reproductive/embryonic compartments (`+1` `proteome_only` `fallopian tube`) under uniform `≤ 2.5%` ectopic olfactory-receptor (`COR*`/`OR*`) gDNA QC. BarkBase `TMM-CPM` (`46` shared SRA runs) is never pooled with `TPM`, and Cufflinks multi-mapping `FPKM` cohorts (`GSE119810`, `GSE95183`, `GSE95184`) mask phase-invariant `≥5%` 31-mer paralog-ambiguous loci and single-exon intronless retrogenes sharing `≥50%` `Unsafe` 9-mers so multi-mapped housekeeping/ribosomal reads never inflate `CA` or `RA`.
- **Pooled 38-tissue proteomics with within-plex TMT bleed subtraction (`pfizer-dog-protein`):** MassIVE `MSV000099744` (`9,929` canine protein groups converted from `rIBAQ` to within-channel molar `ppm`, subtracting `10%` same-plex reproductive reporter-ion bleed) supplies affirmative <span class="cat cat-unsafe">Unsafe</span> (`Vital > 25 ppm`, `Protected > 50 ppm`) and review-hold (`Vital 10–25 ppm`, `Protected 20–50 ppm`) gates with **zero** negative-IHC RNA ceiling multipliers. Because no public canine oviduct RNA-seq dataset exists, `fallopian tube` (`Pool3` `126` in `MSV000099744`) is recorded with `coverage = "proteome_only"`, remaining fail-closed (`unresolved`) when directly protected in `HostContext(protect=("fallopian tube",), species="dog")`.
- **Explicit single-cell UMI gap (`canine_normal_scrna_gap`):** No multi-organ adult canine single-cell UMI atlas comparable to human CELLxGENE or murine `Tabula Muris` is yet integrated (`7/8` pillars active).

### Canine tumor cohorts (7 spontaneous patient cohorts + 1 DLBCL line)

| Cohort | Model | Breeds / DLA | Accession | Samples | Quantification |
| --- | --- | --- | --- | --- | --- |
| `OSA_NCI` — Spontaneous primary & metastatic appendicular osteosarcoma (NCI COTC021/022) | spontaneous_patient_tumor | multi-breed companion dogs (`diverse DLA`) | `GSE238110` | 196 | linear TPM from raw integer counts and CanFam3.1 exonic_bp |
| `OSA_COTC030` — Spontaneous primary appendicular osteosarcoma (NCI COTC030 trial) | spontaneous_patient_tumor | multi-breed companion dogs (`diverse DLA`) | `GSE304066` | 26 | linear TPM from verified raw integer counts and CanFam3.1 exonic_bp |
| `OSA_CCOGC` — Canine Comparative Oncology & Genomics Consortium (CCOGC) osteosarcoma | spontaneous_patient_tumor | Labrador, Greyhound, Rottweiler, Golden Retriever & mixed (`diverse DLA`) | `GSE239948` | 43 | linear TPM from expected counts and CanFam3.1 exonic_bp |
| `CMT` — Spontaneous canine mammary gland carcinoma | spontaneous_patient_tumor | multi-breed female companion dogs (`diverse DLA`) | `GSE119810` | 158 | within-sample linear TPM from FPKM |
| `HSA` — Spontaneous canine visceral hemangiosarcoma | spontaneous_patient_tumor | multi-breed companion dogs (`diverse DLA`) | `GSE95183` | 45 | within-sample linear TPM from FPKM |
| `LYM` — Spontaneous canine multicentric B- and T-cell lymphoma | spontaneous_patient_tumor | multi-breed companion dogs (`diverse DLA`) | `GSE95184` | 76 | within-sample linear TPM from FPKM |
| `DLBCL` — Canine diffuse large B-cell lymphoma (CLBL-1 vehicle controls) | canine_tumor_line | German Shepherd (`DLA-88`) | `GSE287971` | 9 | linear TPM from raw counts and CanFam3.1 exonic_bp |
| `GLIOMA` — Spontaneous prosencephalic high-grade oligodendroglioma | spontaneous_patient_tumor | French Bulldog, Boxer & Boston Terrier (`diverse DLA`) | `GSE284559` | 18 | linear TPM from raw counts and CanFam3.1 exonic_bp |

### Top canine CTAs by peak tumor expression

| Rank | Protein | Sex subset | Peak cohort | Tumor TPM | N TPM | Reproductive peak TPM | Qualifying reproductive compartments |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `CSN3` | `CTA-MF` | `CMT` | 65,068.53 | 5 | 60,012.11 | lactating_mammary  /  resting_mammary |
| 2 | `CSN2` | `CTA-MF` | `CMT` | 30,394.7 | 5.43 | 28,450.54 | lactating_mammary  /  resting_mammary |
| 3 | `CSN1S1` | `CTA-MF` | `CMT` | 24,030.79 | 2.51 | 16,753.21 | testis  /  lactating_mammary  /  resting_mammary |
| 4 | `ENSCAFG00000012690` | `CTA-M` | `LYM` | 174.94 | 2.85 | 522.03 | testis |
| 5 | `ENSCAFG00000013406` | `CTA-M` | `GLIOMA` | 173.04 | 0.68 | 199.5 | testis |
| 6 | `HORMAD1` | `CTA-MF` | `OSA_CCOGC` | 163.11 | 3.44 | 351.1 | testis |
| 7 | `NAA11` | `CTA-M` | `OSA_COTC030` | 158.54 | 6.57 | 100.7 | testis |
| 8 | `GPC2` | `CTA-M` | `HSA` | 152.12 | 10.25 | 126.44 | embryo  /  testis |
| 9 | `PLD6` | `CTA-M` | `OSA_NCI` | 143.72 | 1.54 | 96.42 | testis |
| 10 | `IGF2BP1` | `CTA-P` | `OSA_COTC030` | 137.67 | 1.42 | 46.65 | embryo |
| 11 | `GPR50` | `CTA-P` | `OSA_CCOGC` | 135.72 | 0.73 | 12.26 | embryo |
| 12 | `PRAME` | `CTA-M` | `LYM` | 110.71 | 0.97 | 100.82 | testis |

### Packaged canine nomination ledger (14 literature RT-PCR candidates)

| Candidate | Reference gene ID(s) | Mapping | Source(s) | Native obs. | Catalog call | Missing for quantitative CA/RA |
| --- | --- | --- | --- | --- | --- | --- |
| `AKR1E2` | `608809` | unique | nemec2019 | 2 | `Unsafe` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; quantification_reference |
| `CAGE1` | `609230` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `CTCFL` | `477271` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `DDX53` | `119868264` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `IGF2BP3` | `482357` | unique | nemec2019 | 1 | `Unsafe` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; tumor_observations; quantification_reference |
| `MAGE` | `403710, 492173` | ambiguous | hindriks2025 | 47 | `Other:Expression / Reproductive Antigen only` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference; unique_complete_gene_mapping |
| `MAGEA13P` | `492173` | unique | hindriks2025 | 2 | `Other:Expression` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; somatic_observations; quantification_reference |
| `MAGEB10` | `403710` | unique | hindriks2025 | 1 | `Reproductive Antigen only` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; somatic_observations; tumor_observations; quantification_reference |
| `NT5C1B` | `100856179` | unique | nemec2019, hindriks2025 | 24 | `Unsafe` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `PAGE3-like` | `119868444` | unique | hindriks2025 | 22 | `Unsafe` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `PRAME` | `477552` | unique | hindriks2025 | 22 | `CTA` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `SPECC1` | `489524` | unique | nemec2019 | 2 | `Unsafe` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; quantification_reference |
| `SSX6` | `612524` | unique | hindriks2025 | 22 | `Reproductive Antigen only` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `TPX2` | `477186` | unique | nemec2019 | 2 | `Unsafe` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; quantification_reference |

Evaluating all **14** published canine RT-PCR nominations (`hindriks2025`, `nemec2019`) against genome-wide canine somatic RNA and 38-tissue proteomics shows why quantitative somatic gates are essential: **`PRAME` (`CTA`, `CTA-M`)** is validated as a genuine canine CTA (`110.7 TPM` in spontaneous canine lymphoma, `100.8 TPM` in testis, `0.97 TPM` somatic baseline), `MAGEB10` and `SSX6` qualify as `Reproductive Antigen only`, `MAGEA13P` (`Other:Expression`) and ambiguous `MAGE` (`Other:Expression / Reproductive Antigen only`) stay below `CA`, `CAGE1`, `CTCFL`, and `DDX53` remain held in `Other:Review` on published somatic RT-PCR counterevidence, and **6 of 14 nominations (`AKR1E2`, `IGF2BP3`, `NT5C1B` [mapped to chimeric Ensembl 98 `NT5C1B`+`RDH14` locus `ENSCAFG00000003760`], `PAGE3-like`, `SPECC1`, `TPX2`) are excluded as `Unsafe`** (`1 + 2 + 1 + 1 + 3 + 6 = 14`) due to substantial healthy canine blood, bone marrow, lung, brain, or liver expression.

## *Mus musculus* (laboratory mouse · taxon 10090)

Preclinical vaccine, TCR, and bispecific studies in mice require host-proteome peptide subtraction and endogenous normal-tissue screening on mouse-native coordinates rather than human-to-mouse symbol projection—especially across X-linked CTA families (`Magea`, `Mageb`, `Xlr`, `Trap1a`, `Pramef`) where copy number and tissue regulation diverge sharply from primates.

Murine target selection is evaluated genome-wide across **21,433** Ensembl 115 `GRCm39` protein-equivalence groups (**21,760** primary-assembly coding loci, **8/8** evidence pillars active under harmonized `≤ 2.5%` ectopic `Olfr` gDNA library QC): **12** <span class="cat cat-cta">CTA</span>, **24** <span class="cat cat-ca">CA only</span>, **305** <span class="cat cat-ra">RA only</span> (**341** total candidates), **2,541** <span class="cat cat-other">Other</span>, and **18,551** <span class="cat cat-unsafe">Unsafe</span>. Across the **317** murine Reproductive Antigens (**12** CTAs), 16-compartment `HostContext` partitioning yields **160** `RA-M` (**5** `CTA-M`), **77** `RA-P` (**2** `CTA-P`), **29** `RA-F` (**3** `CTA-F`), and **51** `RA-MF` (**2** `CTA-MF`).

### Murine dataset and threshold parameterization

| Canonical role | Dataset(s) | Native unit | Status | Threshold / screening rule |
| --- | --- | --- | --- | --- |
| Reference proteome & identical-sequence protein-equivalence groups (PEG) | `ensembl115` — Mouse Ensembl 115 GRCm39 (GCA_000001635.9, C57BL/6J, 66,668 proteins) | amino-acid sequence | `active_in_catalog` | Complete GTF/cDNA/ncRNA/peptide FASTA pins; identical translations collapse into murine PEG:<sha256>; >=95% murine 9-mer Unsafe containment |
| Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D) | `bgee_v15_2_mouse_12exps + gse219045_mouse_atlas` — Bgee 15.2 adult murine multi-organ RNA-seq (12 experiments) & GSE219045 sex-stratified preclinical organ atlas | TPM | `active_in_catalog` | Within-organ donor P95 across 30 non-thymic adult somatic contexts sets N; H0 = max(N, D, 1 TPM); thymus stays in thymic scope |
| Normal somatic bulk & sorted immune RNA safety ceilings | `bgee_v15_2_mouse_somatic + gse219045 + immgen_gse109125_gse122108 + gse181223_bmdm + sarma1999` — Bgee 15.2 somatic RNA, GSE219045 organ atlas, ImmGen sorted/activated immune (GSE109125, GSE122108), LPS-stimulated BMDMs (GSE181223) & Sarma 1999 RT-PCR counterevidence | TPM; qualitative RT-PCR/Southern | `active_in_catalog` | Quantitative somatic & sorted/activated immune panels apply Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) TPM without negative-IHC multipliers; Sarma 1999 somatic counterevidence holds Trap1a |
| Single-cell cluster & donor-resolved rare-cell UMI safety screens | `tabula_muris_10x` — Tabula Muris 10x droplet single-cell UMI atlas across 10 adult Vital/Protected somatic organs | UMI counts | `active_in_catalog` | Enters assess_rare_cell_rna (>=10 cells, >=5% single-study or bulk/GXD-confirmed >=1% >=2-UMI in Vital, 2x in Protected) with Ensembl 115 31-mer/host mapping-risk flags |
| Normal-tissue protein staining (IHC) or mass-spectrometry proteomics | `mgi_gxd_6_27` — MGI Gene Expression Database (GXD 6.27) TS28 adult Immunohistochemistry (MMO:0000498) & Immunoblot (MMO:0000672) | curated adult IHC/Immunoblot detections | `active_in_catalog` | >=2 independent publications in Vital/Protected adult somatic organs -> Unsafe; 1 publication -> Other:Review hold; never grants negative-IHC RNA ceiling multipliers |
| Adult reproductive & embryonic/placental enrichment (RA) and host exposure | `bgee_repro_embryo + gse219045 + accessory_tract_lactation_7gse` — Bgee 15.2 & GSE219045 gonads/placenta/embryo plus WT anterior prostate (coagulating_gland), ventral prostate, seminal vesicle, oviduct, uterine horn, endometrium, cervix, vagina & lactating mammary | TPM | `active_in_catalog` | RA qualifies at >=10 x H0 across all 16 murine compartments (including coagulating_gland, seminal vesicle, uterine_horn, lactating_mammary) and populates 16-organ HostContext |
| QC-passing tumor cohorts for cancer enrichment (CA) | `mouse_8_cohorts_syngeneic_gemm` — 8 strain-explicit murine tumor cohorts: MC38, B16F10, LLC, KP_LUAD (C57BL/6); CT26, 4T1 (BALB/c); BMPC_PRAD (FVB); P815 (DBA/2) | TPM | `active_in_catalog` | CA qualifies when >=1 QC-passing control tumor sample reaches >=10 x H0 on Ensembl 115 TPM; cell-line replicates deduplicate donor counts per line |
| Curated normal counterevidence, unscreened adult homes & published nominations | `vandeneynde1991 + sarma1999` — Van den Eynde 1991 (P1A DBA/2 P815 clone 211) & Sarma 1999 endogenous normal-tissue follow-up | curated primary-source XMLs & GenBank records | `active_in_catalog` | Preserves exact 224-aa AAA40475.1 / ENSMUSP00000059448.4 protein equivalence alongside BALB/c lung/spleen/liver/thymus counterevidence |

Key murine methodological rules enforced in code (`oncoref.mouse`, `oncoref.species_target_selection`):

- **Complete pinned Ensembl 115 GRCm39 host proteome:** `mouse.references.ensure_reference("ensembl115")` verifies SHA-256-pinned GTF, cDNA, ncRNA, and protein FASTAs (**66,037** qualified primary-assembly coding translations across **21,433** `PEG`s and **21,760** coding loci), enabling full murine `PEG` grouping, 9-mer <span class="cat cat-unsafe">Unsafe</span> containment, and offline `PeptideIndex` self-proteome subtraction (`mouse-proteome-sequences`).
- **Strain-explicit tumor and normal provenance:** Reference coordinates belong to **C57BL/6J** (`GRCm39`), whereas syngeneic grafts and GEMMs span **C57BL/6** (`MC38`, `B16F10`, `LLC`, `KP_LUAD`), **BALB/c** (`CT26`, `4T1`), **FVB** (`BMPC_PRAD`), and **DBA/2** (`P815`). Cell-line replicates deduplicate donor counts per line rather than inflating independent animal counts.
- **Multi-modal normal & sorted/activated immune safety without negative-IHC inflation:** Adult somatic bulk RNA (`Bgee 15.2` + `GSE219045` across 24 required and additional somatic contexts after `≤ 2.5%` `Olfr` gDNA library QC), sorted & activated immune populations (`ImmGen` ULI `GSE109125`, `ImmGen` MNP `GSE122108`, and WT `LPS`/`IL-4` primary BMDMs `GSE181223`, summarized via finite-sample replicated `min(P95, X_{(n-1)})` across `≥2`-replicate populations), `Tabula Muris` 10x droplet UMI (`assess_rare_cell_rna`), and `MGI GXD 6.27` adult `TS28` immunohistochemistry/immunoblot detections gate safety in their native units with zero negative-IHC RNA ceiling multipliers. Mapped-universe zero-count completion on pre-filtered GEO reproductive count matrices (`GSE287166` `coagulating_gland`, `GSE324827` `fallopian tube`, `GSE60450` mammary) resolves all 16 murine `HostContext` compartments (`Sbp`, `Pbsn`, `Tex19.1`, `Plac1`, `Zfp1010` → `CTA-M`; `Serpinb9e`, `Pramel13` → `CTA-P`; `Gzmg`, `Gzmd`, `Gzmf` → `CTA-F`; `Hmga2`, `Gzme` → `CTA-MF`).

### Murine syngeneic and GEMM tumor cohorts

| Cohort | Model | Host strain (MHC) | Accession | Samples | Quantification |
| --- | --- | --- | --- | --- | --- |
| `MC38` — MC38 subcutaneous colon adenocarcinoma | syngeneic_graft | C57BL/6 (`H-2b`) | `GSE355313` | 4 | linear TPM from raw counts and gene_length |
| `B16F10` — B16-F10 subcutaneous melanoma | syngeneic_graft | C57BL/6 (`H-2b`) | `GSE334005` | 3 | author linear TPM |
| `CT26` — CT26 subcutaneous colon carcinoma | syngeneic_graft | BALB/c (`H-2d`) | `GSE334005` | 4 | author linear TPM |
| `LLC` — Lewis lung carcinoma (LLC) subcutaneous graft | syngeneic_graft | C57BL/6 (`H-2b`) | `GSE334005` | 3 | author linear TPM |
| `4T1` — 4T1 triple-negative mammary carcinoma | syngeneic_line | BALB/c (`H-2d`) | `GSE332912` | 3 | author linear TPM |
| `KP_LUAD` — Primary Kras;Trp53 lung adenocarcinoma GEMM | primary_gemm | C57BL/6 (`H-2b`) | `GSE277076` | 10 | linear TPM from raw counts and Ensembl 115 exonic_bp |
| `BMPC_PRAD` — BMPC MYC;Pten primary & metastatic prostate adenocarcinoma GEMM | primary_gemm | FVB (`H-2q`) | `GSE287166` | 6 | linear TPM from normalized counts and Ensembl 115 exonic_bp |
| `P815` — P815 mastocytoma (native DBA/2 source of P1A/Trap1a) | syngeneic_line | DBA/2 (`H-2d`) | `GSE181773` | 3 | within-sample linear TPM from FPKM |

### Top murine CTAs by peak tumor expression

| Rank | Protein | Sex subset | Peak cohort | Tumor TPM | N TPM | Reproductive peak TPM | Qualifying reproductive compartments |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `Sbp` | `CTA-M` | `BMPC_PRAD` | 729.12 | 3.94 | 320,740.64 | prostate  /  coagulating_gland  /  seminal vesicle |
| 2 | `Gzmg` | `CTA-F` | `MC38` | 600.88 | 1.03 | 187.72 | placenta  /  uterine_horn |
| 3 | `Gzmd` | `CTA-F` | `MC38` | 578.9 | 2.68 | 183.42 | placenta  /  uterine_horn |
| 4 | `Gzmf` | `CTA-F` | `MC38` | 334.69 | 2.49 | 215.57 | placenta  /  uterine_horn |
| 5 | `Pbsn` | `CTA-M` | `BMPC_PRAD` | 312.65 | 2.67 | 60,870.76 | prostate  /  coagulating_gland  /  seminal vesicle |
| 6 | `Hmga2` | `CTA-MF` | `CT26` | 259.46 | 6.93 | 350.18 | embryo |
| 7 | `Gzme` | `CTA-MF` | `MC38` | 223.64 | 1.59 | 94.95 | placenta |
| 8 | `Tex19.1` | `CTA-M` | `P815` | 96.78 | 6 | 98.1 | placenta |
| 9 | `Serpinb9e` | `CTA-P` | `BMPC_PRAD` | 30.81 | 1.21 | 1,088.84 | placenta |
| 10 | `Pramel13` | `CTA-P` | `MC38` | 29.86 | 2.78 | 183.89 | placenta |

### Historical murine nomination check (`P1A` / `Trap1a`)

| Candidate | Reference gene ID(s) | Mapping | Source(s) | Native obs. | Catalog call | Missing for quantitative CA/RA |
| --- | --- | --- | --- | --- | --- | --- |
| `P1A (Trap1a)` | `ENSMUSG00000051257` | unique | vandeneynde1991, sarma1999 | 6 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; quantification_reference |

`P1A` (`Trap1a`, `ENSMUSG00000051257` / `ENSMUST00000054534.4` / `ENSMUSP00000059448.4`, 224 aa) illustrates why OncoRef separates nomination provenance from quantitative admission: the DBA/2 `P815` rejection antigen (`M36387.1` / `AAA40475.1`, Van den Eynde et al. 1991) encodes an identical 224-aa protein to C57BL/6J `ENSMUSP00000059448.4` and reaches high expression in `P815` mastocytoma, while endogenous somatic counterevidence in BALB/c lung, spleen, liver, and thymus (Sarma et al. 1999) plus quantitative multi-organ screening classify `Trap1a` as `Other:Review` rather than admitting it to `CTA`.

## Cross-species CTA orthology across human, dog, and mouse

With independent genome-wide catalogs evaluated natively across **human (`325` CTAs)**, **dog (`54` CTAs)**, and **mouse (`12` CTAs)**, OncoRef quantifies directly how often a `CTA` in one mammal maps to a `CTA` in another ([download full `382`-family ortholog concordance table](antigen-whitepaper/cross-species-orthologs.csv)):

| Cross-species orthology regime (union of human, dog & mouse CTAs) | Ortholog families | Biological mechanism |
| --- | --- | --- |
| **Conserved multi-species CTA** | 9 | Independently passes Stage 1 somatic safety AND 10× cancer + reproductive enrichment in ≥2 species (`ANTXRL`, `CSN2`, `LUZP4`, `PASD1`, `PLAC1`/`Plac1`, `PRAME`, `PSMA8`, `TEX19`/`Tex19.1`, `TRIM42`) |
| **CTA + conserved RA-only ortholog** | 75 | Conserved germline/reproductive restriction (`RA`) in ≥2 species, with ≥10× tumor reactivation observed in the broader cohort panel of one species (`DPPA5`, `MAGEB18`, `CYLC1`, `ACTL8`, `ADAM2`, `AKAP4`) |
| **Divergent somatic / review gate** | 143 | 1:1 ortholog exists across species, but somatic expression diverges evolutionary or is caught by deeper human single-cell/IHC screens (`HORMAD1`, `DDX43`, `CLDN6`, `IGF2BP1`) |
| **Lineage-specific / paralog-diverged CTA** | 155 | Primate-amplified CT-X/subtelomeric families (`NY-ESO-1`, `SPANX*`, `VCX/Y`, `GAGE*`, `XAGE*`, `CT45A*`, `CT47A*`, `PRAMEF*`, `DUX4`), canine-specific loci, or murine `Gm*`/`Gzmd-g`/`Pbsn`/`Sbp` expansions |

Why do only **9** of **382** mammalian `CTA` ortholog families independently qualify as `CTA` in two or more species (and **84** as conserved `CTA` or `RA-only`)? Three evolutionary and experimental mechanisms govern cross-species mapping:

1. **Rapid mammalian lineage amplification of CT-X and ZGA arrays:** Nearly half of human CTAs—including `NY-ESO-1` (`CTAG1A/B`), `SPANXA/B/C/D/N`, `VCX/VCY`, `GAGE*`, `XAGE1A/B`, `CT45A*`, `CT47A*`, `DUX4` (`D4Z4`), `PRAMEF1–26`, `MBD3L2–5`, and `SSU72L1–5`—originated or amplified in primates after boreoeutherian divergence and have **no 1:1 coding ortholog** in dog or mouse (`human_cta_with_dog_orth_n = 131 / 325`, `human_cta_with_mouse_orth_n = 161 / 325`), while rodents independently expanded `Pramel*`, `Mageb*`, uterine NK `Gzmd/e/f/g` (`CTA-F` in mouse), and `Sbp`/`Pbsn` (`CTA-M` in mouse).
2. **Conserved boreoeutherian germline/lactation core (`84` families):** Single-copy autosomal germline, placental, and milk genes conserve reproductive restriction across 95 million years of mammalian evolution—yielding **conserved dual-species CTAs** (`ANTXRL`, `CSN2`, `LUZP4`, `PASD1`, `PLAC1`/`Plac1`, `PRAME`, `PSMA8`, `TEX19`/`Tex19.1`, `TRIM42`), alongside conserved `RA`s (`DPPA5`, `CIMIP4`/`TEX33`, `SPMIP9`/`TEX37`, `MEIOSIN`/`BHMG1`, `TBXT`/`T`, `DEFB123`/`Defb23`, `MAGEB18`, `CYLC1`, `ACTL8`, `ADAM2`, `AKAP4`) whose tumor reactivation depends on cohort breadth (`132` human cohorts vs. `8` canine cohorts vs. `8` murine models).
3. **Divergent somatic regulation and single-cell resolution:** Several conserved 1:1 orthologs (`HORMAD1`, `DDX43`, `CLDN6`, `IGF2BP1`) qualify in canine bulk RNA + 38-tissue proteomics while failing human single-cell cluster, adult cornea, or sorted-immune gates—underscoring why canine single-cell UMI coverage (`canine_normal_scrna_gap`) is tracked explicitly and why target safety is never inferred across species.

<details markdown="1">
<summary>Conserved multi-species CTA and RA-ortholog concordance table (84 families; download full CSV)</summary>

Download all `382` ortholog families: [`cross-species-orthologs.csv`](antigen-whitepaper/cross-species-orthologs.csv).

| Ortholog family | Human call (subset) | Dog call (subset) | Mouse call (subset) | Regime |
| --- | --- | --- | --- | --- |
| `ACTL8` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `ACTL9` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ACTRT2` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ADAD1` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-MF-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ADAM2` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `AKAP4` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ANTXRL` | `CTA` (`CTA-M`) | `CTA` (`CTA-M`) | `Unsafe` | Conserved multi-species CTA |
| `C10ORF120` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `C12ORF50` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `C16ORF90` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `C2ORF78` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `C3ORF22` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `CCDC116` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CCDC166` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CCDC179` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CCDC182` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CCDC27` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CCDC63` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CCT8L2` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `CER1` | `CTA` (`CTA-P`) | `Other:Expression` | `Reproductive Antigen only` (`RA-P-only`) | CTA + conserved RA-only ortholog |
| `CHRNG` | `Unsafe` | `CTA` (`CTA-P`) | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `CIMIP4` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `CSN2` | `CTA` (`CTA-F`) | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-MF-only`) | Conserved multi-species CTA |
| `CYLC1` | `Reproductive Antigen only` (`RA-MF-only`) | `CTA` (`CTA-M`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `DEFB112` | `CTA` (`CTA-M`) | `Unsafe` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `DEFB123` | `CTA` (`CTA-M`) | `Other:Review` | `Reproductive Antigen only` (`RA-MF-only`) | CTA + conserved RA-only ortholog |
| `DEFB126` | `CTA` (`CTA-M`) | `Other:Review` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `DEFB129` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `DMRT1` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `DMRTC2` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `DPPA5` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `DPRX` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-F-only`) | `absent` | CTA + conserved RA-only ortholog |
| `DUSP21` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ESX1` | `CTA` (`CTA-M`) | `Other:Expression` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `FAM187B` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `FMR1NB` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `FOXL3` | `CTA` (`CTA-M`) | `absent` | `Reproductive Antigen only` (`RA-MF-only`) | CTA + conserved RA-only ortholog |
| `FSCB` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `FSCN3` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `GK2` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `IQCA1L` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `LCN9` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `LIPI` | `Unsafe` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `LRRC52` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `LUZP4` | `CTA` (`CTA-M`) | `CTA` (`CTA-M`) | `Cancer Antigen only` | Conserved multi-species CTA |
| `LYPD4` | `CTA` (`CTA-MF`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `MAGEB18` | `Reproductive Antigen only` (`RA-M-only`) | `CTA` (`CTA-M`) | `Other:Expression` | CTA + conserved RA-only ortholog |
| `MEIOSIN` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Other:Expression` | CTA + conserved RA-only ortholog |
| `MROH2B` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ODF1` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `ODF4` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `OLIG3` | `Other:Review` | `CTA` (`CTA-P`) | `Reproductive Antigen only` (`RA-F-only`) | CTA + conserved RA-only ortholog |
| `PASD1` | `CTA` (`CTA-M`) | `CTA` (`CTA-M`) | `Unsafe` | Conserved multi-species CTA |
| `PATE1` | `CTA` (`CTA-M`) | `Other:Review` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `PATE3` | `CTA` (`CTA-M`) | `Other:Review` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `PLAC1` | `CTA` (`CTA-M`) | `absent` | `CTA` (`CTA-M`) | Conserved multi-species CTA |
| `PRAME` | `CTA` (`CTA-MF`) | `CTA` (`CTA-M`) | `absent` | Conserved multi-species CTA |
| `PRDM14` | `CTA` (`CTA-P`) | `Reproductive Antigen only` (`RA-M-only`) | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `PRSS37` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `PRSS58` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `PSMA8` | `CTA` (`CTA-M`) | `CTA` (`CTA-M`) | `Unsafe` | Conserved multi-species CTA |
| `RNASE11` | `CTA` (`CTA-M`) | `absent` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `SMCP` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `SPAG11A` | `CTA` (`CTA-M`) | `absent` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `SPAG11B` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `SPATA19` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `SPEM2` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `SPMIP9` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `SYCP1` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TBC1D21` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TBXT` | `CTA` (`CTA-F`) | `Other:Expression` | `Reproductive Antigen only` (`RA-M-only`) | CTA + conserved RA-only ortholog |
| `TEX19` | `CTA` (`CTA-M`) | `absent` | `CTA` (`CTA-M`) | Conserved multi-species CTA |
| `TEX38` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TEX44` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TEX55` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TGIF2LX` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `TMEM270` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TPD52L3` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |
| `TRIM42` | `CTA` (`CTA-MF`) | `CTA` (`CTA-M`) | `Unsafe` | Conserved multi-species CTA |
| `TRIML1` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `TRIML2` | `Unsafe` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-P-only`) | CTA + conserved RA-only ortholog |
| `UBQLN3` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `Unsafe` | CTA + conserved RA-only ortholog |
| `VRTN` | `CTA` (`CTA-F`) | `Other:Expression` | `Reproductive Antigen only` (`RA-P-only`) | CTA + conserved RA-only ortholog |
| `WFDC10A` | `CTA` (`CTA-M`) | `Reproductive Antigen only` (`RA-M-only`) | `absent` | CTA + conserved RA-only ortholog |

</details>

## Programmatic and CLI access

Every species exposes the identical policy, dataset-binding, catalog, sex-subset, and evaluation surface:

```python
from oncoref import antigen_sets, canine, host_context, human, mouse
from oncoref.species_policies import role_comparison_table

# 1. Inspect declarative per-species thresholds, anatomy, and 8-pillar dataset bindings
human.screening_policy()
canine.screening_datasets()
mouse.screening_datasets()
role_comparison_table()

# 2. Query biological sets and the 4-way physiological sex-restriction partition (RA-M/F/P/MF, CTA-M/F/P/MF)
cta_female_safe = antigen_sets.select("CTA-M")  # 224 male-restricted female-safe human CTAs
ra_female_only = antigen_sets.select("RA-F")  # 18 female-restricted male-safe human RAs
pregnancy_only = antigen_sets.select("CTA-P")  # 48 adult-absent pregnancy/peripartum CTAs
dog_ctas = antigen_sets.select("CTA", species="dog")  # 54 spontaneous canine tumor CTAs
mouse_ctas = antigen_sets.select("CTA", species="mouse")  # 12 murine GRCm39 CTAs

# 3. Parameterize canine or murine HostContext and subtract host self-proteome peptides
spayed_bitch = host_context.HostContext(
    absent=(
        "testis",
        "epididymis",
        "prostate",
        "ovary",
        "uterus",
        "uterine_horn",
        "endometrium",
        "lactating_mammary",
    ),
    protect=("cervix", "vagina", "resting_mammary"),
    species="dog",
)
canine_selection = antigen_sets.select("CTA", species="dog", host=spayed_bitch)
canine_bg = canine_selection.background_index("Unsafe|Other")
```

```sh
# CLI queries across human, dog, and mouse, including sex-restricted subsets
oncoref antigens list --set CTA-M --format table
oncoref antigens list --set CTA-P --format json
oncoref antigens list --set CTA --species dog --format table
oncoref antigens explain PRAME --species dog --format table
oncoref antigens explain P1A --species mouse --format table
```
