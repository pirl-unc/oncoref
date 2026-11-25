# Multi-species antigen whitepaper: human, canine, and murine parameterization

*[Main antigen whitepaper](antigen-whitepaper.md) · [Species & data groups API](species.md) · [Dataset bindings (CSV)](antigen-whitepaper/species-datasets.csv) · [Non-human nominations (CSV)](antigen-whitepaper/species-nominations.csv)*

<p class="section-nav"><strong>Jump to:</strong> <a href="#one-definition-three-species-parameterizations">Universal definitions</a> · <a href="#comparative-anatomy-and-sex-restricted-subsets">Anatomy &amp; sex subsets</a> · <a href="#eight-pillar-cross-species-parameterization-matrix">8-pillar matrix</a> · <a href="#homo-sapiens-human-taxon-9606">Human (9606)</a> · <a href="#canis-lupus-familiaris-domestic-dog-taxon-9615">Dog (9615)</a> · <a href="#mus-musculus-laboratory-mouse-taxon-10090">Mouse (10090)</a> · <a href="#programmatic-and-cli-access">API &amp; CLI</a></p>

Translational immunotherapy moves targets, vaccine designs, and preclinical benchmarks across humans, companion dogs with spontaneous cancers, and syngeneic or genetically engineered mouse models. Treating "CTA" as a human gene list copied by orthology fails in both directions: paralog families diverge rapidly across mammals (`MAGE`, `CT45`, `PRAME`, `SPANX`, `XAGE`), accessory reproductive and lymphoid organs differ between species, and assay scales (`TPM`, `nTPM`, `TMM-CPM`, UMI counts, IHC, pooled TMT `iBAQ`) cannot be pooled or converted by rule of thumb.

OncoRef solves this by separating **species-invariant set definitions** from **per-species anatomical, dataset, and threshold parameterizations** (`oncoref.species_policies`, `oncoref.species_target_selection`). Every species runs through the same three-stage pipeline on its own genome reference and native measurements; where analogous datasets are not yet reconciled for a species, the catalog reports explicit `Other:Review` diagnostics rather than borrowing calls from another species.

## One definition, three species parameterizations

Within any species *S* with verified proteome reference `R_S`, comparative anatomy profile `A_S`, and declared single-panel comparison unit `u_S`:

1. **Protein-equivalence counting unit (`PEG`).** Qualified coding translations within `R_S` collapse into identical-sequence **protein-equivalence groups** (`PEG:<sha256[:20]>`). Within each expression library, complete member-locus abundance is summed before summarizing across donors, while every annotated isoform translation is retained for **≥95% 9-mer containment** against that species' <span class="cat cat-unsafe">Unsafe</span> proteins and downstream `PeptideIndex` subtraction.
2. **Stage 1 — Multi-modal normal safety and shared eligibility.** Each normal-tissue dataset is evaluated in its own native unit and assay rule—bulk somatic and immune RNA ceilings (`Vital ≤ 4`, `Protected ≤ 8`, `Manageable ≤ 16` hold / `> 128` <span class="cat cat-unsafe">Unsafe</span>), donor-resolved rare-cell UMI screens (`≥10` cells and `≥1%` multi-molecule cells across `≥2` independent studies or `≥5%` in one study; `2×` in `Protected`), normal protein staining or proteomics, and curated adult normal-organ literature. Confirmed Vital/Protected adult somatic expression always takes precedence (<span class="cat cat-unsafe">Unsafe</span>); incomplete required normal coverage, unresolved paralog/host attribution, or unscreened adult home organs stay in <span class="cat cat-other">Other</span> (`Other:Review`, `Other:Sequence`, `Other:Unscreened tissue`).
3. **Stage 2 — Symmetric enrichment (`CA`, `RA`, `CTA`).** Across non-thymic, non-reproductive somatic baseline `N` (within-organ donor P95), species-applicable non-thymic `Disposable` maximum `D` (capped at `≤ 32`), and detection floor `F_S = 1` in unit `u_S`, let **H₀ = max(N, D, F_S)**:
   - <span class="cat cat-ca">Cancer Antigen (CA)</span>: clears shared eligibility, has no unscreened adult home organ, and reaches **≥ 10 × H₀** in at least one QC-passing tumor sample of species *S*.
   - <span class="cat cat-ra">Reproductive Antigen (RA)</span>: clears shared eligibility and reaches **≥ 10 × H₀** in an adult reproductive sample or cell-supported embryonic lineage of species *S*.
   - <span class="cat cat-cta">Cancer–Testis Antigen (CTA)</span>: the intersection **CTA = CA ∩ RA**, partitioning the candidate union **CA ∪ RA** into `CTA`, `Cancer Antigen only`, and `Reproductive Antigen only`.

> Negative-IHC RNA ceiling multipliers (`2×` Approved, `4×` Supported, `8×` Enhanced) apply **only** where a species has routine-context, cell-type-resolved antibody validation (`human` HPA v23). Pooled mass-spectrometry proteomics (`dog` Pfizer TMT-iBAQ) and literature-curated staining (`mouse` GXD) supply affirmative <span class="cat cat-unsafe">Unsafe</span> / hold evidence when positive, and never inflate RNA ceilings when absent.

## Comparative anatomy and sex-restricted subsets

Male- and female-restricted subsets are defined by the **same `HostContext` protection algebra** in every species, parameterized by that species' actual reproductive, mammary, thymic, and disposable organs in `normal-tissue-taxonomy.yaml`:

- **`MRA` (`CTA-M` / `MRA-only`, Male-restricted · female-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every female adult reproductive and mammary compartment of species *S* (`< 4` native RNA units with measured coverage and resolved protein status), while expressed in male organs.
- **`FRA` (`CTA-F` / `FRA-only`, Female-restricted · male-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every male adult reproductive compartment plus resting mammary tissue of species *S*, while expressed in female or placental compartments.
- **`ERA` (`CTA-E` / `ERA-only`, Embryo / placenta-only · adult-absent in both sexes):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible under full adult organ protection in **both** sexes—expressed in embryo or placenta and silent across every adult male and non-pregnant female reproductive and mammary compartment.
- **`SRA` (`CTA-MF` / `SRA-only`, Shared male & female):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) expressed in adult reproductive/mammary organs of both sexes (or lacking measured opposite-sex clearance, failing closed).

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
| **Reference proteome & identical-sequence protein-equivalence groups (PEG)** | ensembl115 (amino-acid sequence) — *active_in_catalog* | refseq106 + ensembl98 (amino-acid sequence) — *active_in_catalog* | ensembl115 (amino-acid sequence) — *active_in_catalog* |
| **Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D)** | gtex_v11 + hpa_v23_biopsies + sorted_immune (TPM) — *active_in_catalog* | bgee-srp009687 + bgee-srp181218  /  barkbase (TPM or TMM-CPM (single-scale panel)) — *ingested_pending_crosswalk* | bgee_v15_2_mouse_12exps + gse219045_mouse_atlas (TPM) — *active_in_catalog* |
| **Normal somatic bulk & sorted immune RNA safety ceilings** | hpa_rna_consensus + gtex_tissue_rna + hpa_brain + dice + monaco (TPM; nTPM) — *active_in_catalog* | bgee-srp009687 + bgee-srp181218 + barkbase + doga + pfizer-dog-rna (TPM; TMM-CPM; STRT2 ERCC (metadata)) — *ingested_pending_crosswalk* | bgee_v15_2_mouse_somatic + gse219045 + sarma1999 (TPM; qualitative RT-PCR/Southern) — *active_in_catalog* |
| **Single-cell cluster & donor-resolved rare-cell UMI safety screens** | cellxgene_census + siletti_brain + hpa_single_cell + collin_cornea (UMI counts; nTPM; CP10K) — *active_in_catalog* | canine_normal_scrna_gap (UMI counts) — *identified_analog* | tabula_muris_10x (UMI counts) — *active_in_catalog* |
| **Normal-tissue protein staining (IHC) or mass-spectrometry proteomics** | hpa_normal_tissue + external_ihc (categorical IHC grade + reliability) — *active_in_catalog* | pfizer-dog-protein + pfizer-dog-peptide-evidence (pooled iBAQ + peptide assignments) — *ingested_pending_crosswalk* | mgi_gxd_6_27 (curated adult IHC/Immunoblot detections) — *active_in_catalog* |
| **Adult reproductive & embryonic/placental enrichment (RA) and host exposure** | gtex_repro + hpa_repro + petropoulos + tyser + mole (TPM; UMI/staining context) — *active_in_catalog* | bgee_repro + barkbase_repro_embryo (TPM or TMM-CPM) — *ingested_pending_crosswalk* | bgee_repro_embryo + gse219045 + accessory_tract_lactation_7gse (TPM) — *active_in_catalog* |
| **QC-passing tumor cohorts for cancer enrichment (CA)** | oncoref_tumor_union_132_cohorts (TPM) — *active_in_catalog* | gse304066-osteosarcoma + hindriks2025 + nemec2019 (blocked raw/FPKM; qualitative study observations) — *metadata_or_blocked* | mouse_8_cohorts_syngeneic_gemm (TPM) — *active_in_catalog* |
| **Curated normal counterevidence, unscreened adult homes & published nominations** | normal_expression_literature + rare_cell_attribution_reviews + cta_literature_149 (curated verdicts & citations) — *active_in_catalog* | hindriks2025 + nemec2019 (source-native qualitative/semi-quantitative) — *native_seed_only* | vandeneynde1991 + sarma1999 (curated primary-source XMLs & GenBank records) — *active_in_catalog* |

## *Homo sapiens* (human · taxon 9606)

Human target selection is fully evaluated across **19,893** protein-equivalence groups (**20,099** Ensembl 115 coding loci) under policy `oncoref.target-selection.2026-11-23.1`: **332** <span class="cat cat-cta">CTA</span>, **117** <span class="cat cat-ca">CA only</span>, **133** <span class="cat cat-ra">RA only</span> (**582** total candidates), **1,062** <span class="cat cat-other">Other</span>, and **18,249** <span class="cat cat-unsafe">Unsafe</span>. Across the **465** Reproductive Antigens (**332** CTAs), sex-restriction partitioning yields **340** `MRA` (**232** `CTA-M`), **52** `ERA` (**44** `CTA-E`), **10** `FRA` (**5** `CTA-F`), and **63** `SRA` (**51** `CTA-MF`). See the [main antigen whitepaper](antigen-whitepaper.md) for full human figures, clinical anchor tables, and edge-case audits.

| Canonical role | Dataset(s) | Native unit | Status | Threshold / screening rule |
| --- | --- | --- | --- | --- |
| Reference proteome & identical-sequence protein-equivalence groups (PEG) | `ensembl115` — Human Ensembl 115 GRCh38 + MANE Select | amino-acid sequence | `active_in_catalog` | Prefer MANE Select then Ensembl canonical; collapse identical translations into PEG; >=95% 9-mer Unsafe containment -> Other:Sequence |
| Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D) | `gtex_v11 + hpa_v23_biopsies + sorted_immune` — GTEx v11, HPA v23 tissues/brain & HPA/DICE/Monaco immune | TPM | `active_in_catalog` | Within-site linear donor P95 sets non-thymic somatic N and Disposable D (<=32 TPM); H0 = max(N, D, 1 TPM) |
| Normal somatic bulk & sorted immune RNA safety ceilings | `hpa_rna_consensus + gtex_tissue_rna + hpa_brain + dice + monaco` — HPA v23 consensus/brain/immune, GTEx v11 medians & DICE/Monaco | TPM; nTPM | `active_in_catalog` | Base ceilings Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) in native unit; matched negative IHC raises ceiling 2x/4x/8x |
| Single-cell cluster & donor-resolved rare-cell UMI safety screens | `cellxgene_census + siletti_brain + hpa_single_cell + collin_cornea` — CELLxGENE Census, Siletti brain, HPA v23 scRNA clusters & Collin cornea | UMI counts; nTPM; CP10K | `active_in_catalog` | UMI: >=10 cells and >=1% >=2-UMI cells across >=2 consortium-deduplicated study groups (or >=5% in 1 study; 2x in Protected); cornea CP10K /100 nominal ceiling |
| Normal-tissue protein staining (IHC) or mass-spectrometry proteomics | `hpa_normal_tissue + external_ihc` — HPA v23 normal immunohistochemistry & curated external staining | categorical IHC grade + reliability | `active_in_catalog` | Approved/Supported/Enhanced positive staining in Vital/Protected tissues -> Unsafe; routine-context negative IHC grants RNA headroom |
| Adult reproductive & embryonic/placental enrichment (RA) and host exposure | `gtex_repro + hpa_repro + petropoulos + tyser + mole` — GTEx/HPA adult reproductive samples & Petropoulos/Tyser/Molè embryos | TPM; UMI/staining context | `active_in_catalog` | RA qualifies at >=10 x H0 in adult reproductive samples or cell-supported embryo lineages; per-compartment maxima feed HostContext (<4 TPM / <4 nTPM) |
| QC-passing tumor cohorts for cancer enrichment (CA) | `oncoref_tumor_union_132_cohorts` — 132 QC-passing per-sample human tumor cohorts (TCGA, TARGET, BEAT-AML, Treehouse, GEO) | TPM | `active_in_catalog` | CA qualifies when >=1 QC-passing library reaches >=10 x H0; known-patient counts deduplicate samples per donor |
| Curated normal counterevidence, unscreened adult homes & published nominations | `normal_expression_literature + rare_cell_attribution_reviews + cta_literature_149` — Reviewed normal-tissue literature, #711 Ensembl 115 attribution ledger & 149 CTA/oncofetal papers | curated verdicts & citations | `active_in_catalog` | Verified adult Vital/Protected expression -> Unsafe; unscreened adult home (olfactory, taste bud, bone, dental pulp, hair follicle) -> Other:Unscreened tissue |

## *Canis lupus familiaris* (domestic dog · taxon 9615)

Spontaneous canine cancers—osteosarcoma, diffuse large B-cell lymphoma, urothelial carcinoma, mammary carcinoma, oral melanoma, and glioma—develop in outbred, immunocompetent hosts sharing human environmental exposures. OncoRef ships pinned canine genome references, native normal RNA and proteomics readers, comparative canine `HostContext` profiles, and **14** literature-curated canine CTA nominations (`hindriks2025`, `nemec2019`).

### Canine dataset and threshold parameterization

| Canonical role | Dataset(s) | Native unit | Status | Threshold / screening rule |
| --- | --- | --- | --- | --- |
| Reference proteome & identical-sequence protein-equivalence groups (PEG) | `refseq106 + ensembl98` — Canine NCBI RefSeq 106 (UU_Cfam_GSD_1.0) & Ensembl 98 (CanFam3.1) | amino-acid sequence | `active_in_catalog` | Complete GTF/FASTA/feature_table pins; identical translations collapse into canine PEG:<sha256>; >=95% canine 9-mer Unsafe containment |
| Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D) | `bgee-srp009687 + bgee-srp181218  /  barkbase` — Bgee 15.2 canine normal RNA-seq (13 + 75 libraries, TPM) & BarkBase (150 libraries, TMM-CPM) | TPM or TMM-CPM (single-scale panel) | `ingested_pending_crosswalk` | Within-organ donor P95 sets somatic N and canine Disposable D (gallbladder, tonsil, hair follicle, cutaneous melanocyte; no appendix); never pool TPM with TMM-CPM |
| Normal somatic bulk & sorted immune RNA safety ceilings | `bgee-srp009687 + bgee-srp181218 + barkbase + doga + pfizer-dog-rna` — Bgee TPM, BarkBase TMM-CPM & DoGA / Pfizer canine RNA metadata | TPM; TMM-CPM; STRT2 ERCC (metadata) | `ingested_pending_crosswalk` | Tiered somatic ceilings Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) in declared panel unit without negative-IHC inflation |
| Single-cell cluster & donor-resolved rare-cell UMI safety screens | `canine_normal_scrna_gap` — Donor-resolved canine normal single-cell UMI atlases | UMI counts | `identified_analog` | Same assess_rare_cell_rna contract (>=10 cells, >=1% >=2-UMI across >=2 canine studies, 2x in Protected) when single-cell atlases are ingested |
| Normal-tissue protein staining (IHC) or mass-spectrometry proteomics | `pfizer-dog-protein + pfizer-dog-peptide-evidence` — Pfizer 38-tissue canine quantitative proteomics (MSV000099744, 9,929 groups) | pooled iBAQ + peptide assignments | `ingested_pending_crosswalk` | Affirmative Vital/Protected tissue-pool protein detection supplies Unsafe/hold review once searched-FASTA groups reconcile to canine genes; never grants negative-IHC RNA credit |
| Adult reproductive & embryonic/placental enrichment (RA) and host exposure | `bgee_repro + barkbase_repro_embryo` — Bgee & BarkBase canine testis, ovary, oviduct, uterine horn, prostate, mammary & embryonic libraries | TPM or TMM-CPM | `ingested_pending_crosswalk` | RA qualifies at >=10 x H0 across canine compartments (includes uterine_horn; omits absent seminal vesicle and bulbourethral gland); feeds canine HostContext |
| QC-passing tumor cohorts for cancer enrichment (CA) | `gse304066-osteosarcoma + hindriks2025 + nemec2019` — GSE304066 canine osteosarcoma & published canine tumor RT-qPCR/RNA observations | blocked raw/FPKM; qualitative study observations | `metadata_or_blocked` | CA requires >=1 QC-passing tumor library at >=10 x H0 on the same reconciled reference and linear unit as the canine normal panel |
| Curated normal counterevidence, unscreened adult homes & published nominations | `hindriks2025 + nemec2019` — Curated canine CTA nominations and normal/tumor study observations (14 candidates) | source-native qualitative/semi-quantitative | `native_seed_only` | Retains native canine gene mappings, paralog assay ambiguity, and positive normal somatic observations; packaged rows remain Other:Review until panels reconcile |

Key canine methodological rules enforced in code (`oncoref.canine`, `oncoref.species_target_selection`):

- **Single-scale enrichment panels (`TPM` vs. `TMM-CPM`):** Bgee (`bgee-srp009687`, `bgee-srp181218`) quantifies canine libraries in `TPM`, whereas BarkBase (`barkbase`) reports StringTie genes in `TMM-CPM` and shares **46 SRA runs** with the global Bgee inventory (`canine.expression.lineage_overlaps()`). `SpeciesScreeningSpec` accepts either a pure `TPM` panel or a pure `TMM-CPM` panel, forbidding cross-unit pooling or duplicate donor counting across shared runs.
- **Pooled proteomics (`pfizer-dog-protein`):** MassIVE `MSV000099744` provides quantitative `iBAQ` across **38 pooled canine tissues** (**9,929** protein groups) plus optional peptide-to-protein assignments (`pfizer-dog-peptide-evidence`). Because every channel pools multiple animals without cell-type resolution, `assess_pooled_protein_gates` uses high Vital/Protected `iBAQ` only as affirmative <span class="cat cat-unsafe">Unsafe</span> / hold evidence once searched-FASTA accessions map to canine genes, and never awards negative-IHC RNA ceiling multipliers.
- **Explicit data gaps blocking packaged genome-wide calls:** Source quantification GTF crosswalks for Bgee/BarkBase and raw-count versus FPKM provenance in `gse304066-osteosarcoma` remain unresolved, so packaged canine nominations stay `Other:Review` rather than inheriting human CTA labels.

### Packaged canine nomination ledger (14 native candidates)

| Candidate | Reference gene ID(s) | Mapping | Source(s) | Native obs. | Catalog call | Missing for quantitative CA/RA |
| --- | --- | --- | --- | --- | --- | --- |
| `AKR1E2` | `608809` | unique | nemec2019 | 2 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; quantification_reference |
| `CAGE1` | `609230` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `CTCFL` | `477271` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `DDX53` | `119868264` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `IGF2BP3` | `482357` | unique | nemec2019 | 1 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; tumor_observations; quantification_reference |
| `MAGE` | `403710, 492173` | ambiguous | hindriks2025 | 47 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference; unique_complete_gene_mapping |
| `MAGEA13P` | `492173` | unique | hindriks2025 | 2 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; somatic_observations; quantification_reference |
| `MAGEB10` | `403710` | unique | hindriks2025 | 1 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; somatic_observations; tumor_observations; quantification_reference |
| `NT5C1B` | `100856179` | unique | nemec2019, hindriks2025 | 24 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `PAGE3-like` | `119868444` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `PRAME` | `477552` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `SPECC1` | `489524` | unique | nemec2019 | 2 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; quantification_reference |
| `SSX6` | `612524` | unique | hindriks2025 | 22 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; quantification_reference |
| `TPX2` | `477186` | unique | nemec2019 | 2 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; reproductive_observations; quantification_reference |

## *Mus musculus* (laboratory mouse · taxon 10090)

Preclinical vaccine, TCR, and bispecific studies in mice require host-proteome peptide subtraction and endogenous normal-tissue screening on mouse-native coordinates rather than human-to-mouse symbol projection—especially across X-linked CTA families (`Magea`, `Mageb`, `Xlr`, `Trap1a`, `Pramef`) where copy number and tissue regulation diverge sharply from primates.

Murine target selection is evaluated genome-wide across **21,433** Ensembl 115 `GRCm39` protein-equivalence groups (**21,760** primary-assembly coding loci, **8/8** evidence pillars active): **20** <span class="cat cat-cta">CTA</span>, **29** <span class="cat cat-ca">CA only</span>, **296** <span class="cat cat-ra">RA only</span> (**345** total candidates), **2,593** <span class="cat cat-other">Other</span>, and **18,495** <span class="cat cat-unsafe">Unsafe</span>. Across the **316** murine Reproductive Antigens (**20** CTAs), 16-compartment `HostContext` sex-restriction partitioning yields **69** `MRA` (**3** `CTA-M`), **8** `ERA` (**0** `CTA-E`), **17** `FRA` (**2** `CTA-F`), and **222** `SRA` (**15** `CTA-MF`).

### Murine dataset and threshold parameterization

| Canonical role | Dataset(s) | Native unit | Status | Threshold / screening rule |
| --- | --- | --- | --- | --- |
| Reference proteome & identical-sequence protein-equivalence groups (PEG) | `ensembl115` — Mouse Ensembl 115 GRCm39 (GCA_000001635.9, C57BL/6J, 66,668 proteins) | amino-acid sequence | `active_in_catalog` | Complete GTF/cDNA/ncRNA/peptide FASTA pins; identical translations collapse into murine PEG:<sha256>; >=95% murine 9-mer Unsafe containment |
| Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D) | `bgee_v15_2_mouse_12exps + gse219045_mouse_atlas` — Bgee 15.2 adult murine multi-organ RNA-seq (12 experiments) & GSE219045 sex-stratified preclinical organ atlas | TPM | `active_in_catalog` | Within-organ donor P95 across 30 non-thymic adult somatic contexts sets N; H0 = max(N, D, 1 TPM); thymus stays in thymic scope |
| Normal somatic bulk & sorted immune RNA safety ceilings | `bgee_v15_2_mouse_somatic + gse219045 + sarma1999` — Bgee 15.2 somatic/leukocyte/marrow RNA, GSE219045 whole-blood/lymph-node/CNS & Sarma 1999 RT-PCR counterevidence | TPM; qualitative RT-PCR/Southern | `active_in_catalog` | Quantitative panels apply Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) TPM without negative-IHC multipliers; Sarma 1999 somatic counterevidence holds Trap1a |
| Single-cell cluster & donor-resolved rare-cell UMI safety screens | `tabula_muris_10x` — Tabula Muris 10x droplet single-cell UMI atlas across 10 adult Vital/Protected somatic organs | UMI counts | `active_in_catalog` | Enters assess_rare_cell_rna (>=10 cells, >=5% single-study or bulk/GXD-confirmed >=1% >=2-UMI in Vital, 2x in Protected) with Ensembl 115 31-mer/host mapping-risk flags |
| Normal-tissue protein staining (IHC) or mass-spectrometry proteomics | `mgi_gxd_6_27` — MGI Gene Expression Database (GXD 6.27) TS28 adult Immunohistochemistry (MMO:0000498) & Immunoblot (MMO:0000672) | curated adult IHC/Immunoblot detections | `active_in_catalog` | >=2 independent publications in Vital/Protected adult somatic organs -> Unsafe; 1 publication -> Other:Review hold; never grants negative-IHC RNA ceiling multipliers |
| Adult reproductive & embryonic/placental enrichment (RA) and host exposure | `bgee_repro_embryo + gse219045 + accessory_tract_lactation_7gse` — Bgee 15.2 & GSE219045 gonads/placenta/embryo plus WT anterior prostate (coagulating_gland), ventral prostate, seminal vesicle, oviduct, uterine horn, endometrium, cervix, vagina & lactating mammary | TPM | `active_in_catalog` | RA qualifies at >=10 x H0 across all 16 murine compartments (including coagulating_gland, seminal vesicle, uterine_horn, lactating_mammary) and populates 16-organ HostContext |
| QC-passing tumor cohorts for cancer enrichment (CA) | `mouse_8_cohorts_syngeneic_gemm` — 8 strain-explicit murine tumor cohorts: MC38, B16F10, LLC, KP_LUAD (C57BL/6); CT26, 4T1 (BALB/c); BMPC_PRAD (FVB); P815 (DBA/2) | TPM | `active_in_catalog` | CA qualifies when >=1 QC-passing control tumor sample reaches >=10 x H0 on Ensembl 115 TPM; cell-line replicates deduplicate donor counts per line |
| Curated normal counterevidence, unscreened adult homes & published nominations | `vandeneynde1991 + sarma1999` — Van den Eynde 1991 (P1A DBA/2 P815 clone 211) & Sarma 1999 endogenous normal-tissue follow-up | curated primary-source XMLs & GenBank records | `active_in_catalog` | Preserves exact 224-aa AAA40475.1 / ENSMUSP00000059448.4 protein equivalence alongside BALB/c lung/spleen/liver/thymus counterevidence |

Key murine methodological rules enforced in code (`oncoref.mouse`, `oncoref.species_target_selection`):

- **Complete pinned Ensembl 115 GRCm39 host proteome:** `mouse.references.ensure_reference("ensembl115")` verifies SHA-256-pinned GTF, cDNA, ncRNA, and protein FASTAs (**66,037** qualified primary-assembly coding translations across **21,433** `PEG`s and **21,760** coding loci), enabling full murine `PEG` grouping, 9-mer <span class="cat cat-unsafe">Unsafe</span> containment, and offline `PeptideIndex` self-proteome subtraction (`mouse-proteome-sequences`).
- **Strain-explicit tumor and normal provenance:** Reference coordinates belong to **C57BL/6J** (`GRCm39`), whereas syngeneic grafts and GEMMs span **C57BL/6** (`MC38`, `B16F10`, `LLC`, `KP_LUAD`), **BALB/c** (`CT26`, `4T1`), **FVB** (`BMPC_PRAD`), and **DBA/2** (`P815`). Cell-line replicates deduplicate donor counts per line rather than inflating independent animal counts.
- **Multi-modal normal safety without negative-IHC inflation:** Adult somatic bulk RNA (`Bgee 15.2` + `GSE219045` across 24 required and additional somatic contexts), `Tabula Muris` 10x droplet UMI (`assess_rare_cell_rna`), and `MGI GXD 6.27` adult `TS28` immunohistochemistry/immunoblot detections gate safety in their native units with zero negative-IHC RNA ceiling multipliers.

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
| 1 | `Sbp` | `CTA-M` | `BMPC_PRAD` | 729.12 | 4.53 | 320,740.64 | prostate  /  coagulating_gland  /  seminal vesicle |
| 2 | `Gzmg` | `CTA-MF` | `MC38` | 600.88 | 1.29 | 187.72 | placenta  /  uterine_horn |
| 3 | `Gzmd` | `CTA-MF` | `MC38` | 578.9 | 2.76 | 183.42 | placenta  /  uterine_horn |
| 4 | `Gzmf` | `CTA-F` | `MC38` | 334.69 | 2.67 | 215.57 | placenta  /  uterine_horn |
| 5 | `Pbsn` | `CTA-MF` | `BMPC_PRAD` | 312.65 | 3.15 | 60,870.76 | prostate  /  coagulating_gland  /  seminal vesicle |
| 6 | `Hmga2` | `CTA-MF` | `CT26` | 259.46 | 6.93 | 350.18 | embryo |
| 7 | `Gzme` | `CTA-MF` | `MC38` | 223.64 | 1.74 | 94.95 | placenta |
| 8 | `Utf1` | `CTA-MF` | `BMPC_PRAD` | 130.78 | 1.07 | 61.07 | embryo  /  placenta |
| 9 | `Csf3` | `CTA-MF` | `4T1` | 114.79 | 1.1 | 1,031.21 | ovary  /  fallopian tube  /  uterine_horn  /  endometrium  /  lactating_mammary  /  resting_mammary |
| 10 | `Tex19.1` | `CTA-M` | `P815` | 96.78 | 6 | 98.1 | placenta |

### Historical murine nomination check (`P1A` / `Trap1a`)

| Candidate | Reference gene ID(s) | Mapping | Source(s) | Native obs. | Catalog call | Missing for quantitative CA/RA |
| --- | --- | --- | --- | --- | --- | --- |
| `P1A (Trap1a)` | `ENSMUSG00000051257` | unique | vandeneynde1991, sarma1999 | 6 | `Other:Review` | complete_tissue_panel; independent_sample_total; independent_donor_total; measured_negative_or_censored_observations; quantification_reference |

`P1A` (`Trap1a`, `ENSMUSG00000051257` / `ENSMUST00000054534.4` / `ENSMUSP00000059448.4`, 224 aa) illustrates why OncoRef separates nomination provenance from quantitative admission: the DBA/2 `P815` rejection antigen (`M36387.1` / `AAA40475.1`, Van den Eynde et al. 1991) encodes an identical 224-aa protein to C57BL/6J `ENSMUSP00000059448.4` and reaches high expression in `P815` mastocytoma, while endogenous somatic counterevidence in BALB/c lung, spleen, liver, and thymus (Sarma et al. 1999) plus quantitative multi-organ screening classify `Trap1a` as `Other:Review` rather than admitting it to `CTA`.

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

# 2. Query biological sets and sex-restricted subsets (MRA/FRA/ERA/SRA, CTA-M/F/E/MF)
cta_female_safe = antigen_sets.select("CTA-M")          # 232 male-restricted female-safe human CTAs
era_targets = antigen_sets.select("CTA-E")              # 44 adult-absent embryo/placenta CTAs
human.sex_restricted_antigens()

# 3. Parameterize canine or murine HostContext and run full catalog evaluation on native panels
spayed_bitch = host_context.HostContext(
    absent=("ovary", "uterus", "uterine_horn", "endometrium"),
    protect=("resting_mammary",),
    species="dog",
)
spec = canine.screening_policy().to_screening_spec(
    comparison_id="canine-tmmcpm-v1",
    comparability_reference="canine-panel-receipt",
    comparison_unit="TMM-CPM",
    required_normal_contexts=frozenset({"liver", "heart", "lung", "kidney", "brain"}),
)
# result = canine.evaluate_target_catalog(spec, proteins, sequences=seqs, ...)
# canine_cta_m = result.select("CTA-M")
```

```sh
# CLI queries across human, dog, and mouse, including sex-restricted subsets
oncoref antigens list --set CTA-M --format table
oncoref antigens list --set CTA-E --format json
oncoref antigens funnel --species dog --format table
oncoref antigens explain P1A --species mouse --format table
```
