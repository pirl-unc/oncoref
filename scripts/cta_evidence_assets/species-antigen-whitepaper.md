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

$comparative_anatomy_table

Comparative anatomy changes both enrichment denominators and host protection:

- **Accessory male glands:** Dogs anatomically lack seminal vesicles and bulbourethral glands (`prostate` is the sole major accessory gland), whereas mice carry a distinct `coagulating_gland` (anterior prostate) alongside `seminal vesicle` and `prostate`. Requesting `seminal vesicle` in a canine `HostContext` raises an immediate `ValueError` rather than silently dropping the constraint.
- **Bicornuate female tract & gonadectomy:** Both dogs and mice add `uterine_horn` (normalized automatically under `uterus`), and companion-animal oncology frequently treats spayed (ovariohysterectomized) or castrated hosts where `HostContext(absent=("ovary", "uterus", "uterine_horn", "endometrium"), species="dog")` explicitly records surgically absent organs distinct from accepted tissue injury.
- **Disposable somatic organs (`D`):** Neither dogs nor mice possess a vermiform `appendix`, and mice also lack palatine `tonsil`, narrowing non-thymic `Disposable` (`D`) from five compartments in human to four in dog and three in mouse.

## Eight-pillar cross-species parameterization matrix

The table below aligns all eight canonical screening roles across human, dog, and mouse ([download full binding table](antigen-whitepaper/species-datasets.csv)):

$eight_pillar_matrix

## *Homo sapiens* (human · taxon 9606)

Human target selection is fully evaluated across **$n** protein-equivalence groups (**$loci** Ensembl 115 coding loci) under policy `$policy`: **$cta** <span class="cat cat-cta">CTA</span>, **$ca_only** <span class="cat cat-ca">CA only</span>, **$ra_only** <span class="cat cat-ra">RA only</span> (**$union** total candidates), **$other_n** <span class="cat cat-other">Other</span>, and **$unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Stratifying GTEx v11 female reproductive and mammary organs to reproductive-age `XX` women (`ages 20–49`, `gtex.donor_age()`) and integrating stage-replicated adult human folliculogenesis/oocyte RNA-seq (`GSE107746`, `128` QC-passing stages) partitions the **$ra** Reproductive Antigens (**$cta** CTAs) into **$ra_m_n** `RA-M` (**$cta_m_n** `CTA-M`), **$ra_p_n** `RA-P` (**$cta_p_n** `CTA-P`), **$ra_f_n** `RA-F` (**$cta_f_n** `CTA-F`), and **$ra_mf_n** `RA-MF` (**$cta_mf_n** `CTA-MF`). See the [main antigen whitepaper](antigen-whitepaper.md) for full human figures, clinical anchor tables, and edge-case audits.

$human_datasets_table

## *Canis lupus familiaris* (domestic dog · taxon 9615)

Spontaneous canine cancers—appendicular osteosarcoma, multicentric and diffuse large B-cell lymphoma, visceral hemangiosarcoma, mammary carcinoma, and high-grade oligodendroglioma—develop in outbred, immunocompetent companion dogs sharing human environmental exposures. Canine target selection is evaluated genome-wide across **$dog_n** Ensembl 98 `CanFam3.1` (`GCA_000002285.2`, cross-referenced to NCBI `RefSeq 106` `UU_Cfam_GSD_1.0`) protein-equivalence groups (**$dog_loci** coding loci, **7/8** evidence pillars active): **$dog_cta** <span class="cat cat-cta">CTA</span>, **$dog_ca_only** <span class="cat cat-ca">CA only</span>, **$dog_ra_only** <span class="cat cat-ra">RA only</span> (**$dog_union** total candidates), **$dog_other_n** <span class="cat cat-other">Other</span>, and **$dog_unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Across the **$dog_ra** canine Reproductive Antigens (**$dog_cta** CTAs), comparative canine `HostContext` partitioning yields **$dog_mra_n** `RA-M` (**$dog_cta_m_n** `CTA-M`), **$dog_era_n** `RA-P` (**$dog_cta_p_n** `CTA-P`), **$dog_fra_n** `RA-F` (**$dog_cta_f_n** `CTA-F`), and **$dog_sra_n** `RA-MF` (**$dog_cta_mf_n** `CTA-MF`).

### Canine dataset and threshold parameterization

$dog_datasets_table

Key canine methodological rules enforced in code (`oncoref.canine`, `oncoref.species_target_selection`, `scripts/build_dog_target_catalog.py`):

- **Single-scale linear `TPM` panels & Cufflinks paralog guard:** `Bgee 15.2` (`25` healthy canine RNA-seq experiments deduplicated by SRA Run ID, routing BarkBase gestational `d33–d44` libraries in `SRP114662` to `embryo`), the **300-sample** sex-stratified Preclinical Species Gene Expression Database (`GSE219045`, normalized within-sample from `FPKM` to exact `TPM`), and **64** normal mammary biopsies (`GSE119810`, split by coordinated 4-gene milk-operon induction into `lactating_mammary` vs. `resting_mammary`) cover all **22** required somatic contexts plus `14` measured reproductive/embryonic compartments under uniform `≤ 2.5%` ectopic olfactory-receptor (`COR*`/`OR*`) gDNA QC. BarkBase `TMM-CPM` (`46` shared SRA runs) is never pooled with `TPM`, and Cufflinks multi-mapping `FPKM` cohorts (`GSE119810`, `GSE95183`, `GSE95184`) mask `≥5%` 31-mer paralog-ambiguous retrocopies so multi-mapped ribosomal reads never inflate `CA` or `RA`.
- **Pooled 38-tissue proteomics with within-plex TMT bleed subtraction (`pfizer-dog-protein`):** MassIVE `MSV000099744` (`9,929` canine protein groups converted from `rIBAQ` to within-channel molar `ppm`, subtracting `10%` same-plex reproductive reporter-ion bleed) supplies affirmative <span class="cat cat-unsafe">Unsafe</span> (`Vital > 25 ppm`, `Protected > 50 ppm`) and review-hold (`Vital 10–25 ppm`, `Protected 20–50 ppm`) gates with **zero** negative-IHC RNA ceiling multipliers. Because no public canine oviduct RNA-seq dataset exists, `fallopian tube` (`Pool3` `126` in `MSV000099744`) is recorded with `coverage = "proteome_only"`, remaining fail-closed (`unresolved`) when directly protected in `HostContext(protect=("fallopian tube",), species="dog")`.
- **Explicit single-cell UMI gap (`canine_normal_scrna_gap`):** No multi-organ adult canine single-cell UMI atlas comparable to human CELLxGENE or murine `Tabula Muris` is yet integrated (`7/8` pillars active).

### Spontaneous canine patient tumor cohorts

$dog_cohorts_table

### Top canine CTAs by peak tumor expression

$dog_top_cta_table

### Packaged canine nomination ledger (14 literature RT-PCR candidates)

$dog_nominations_table

Evaluating the **14** published canine RT-PCR nominations (`hindriks2025`, `nemec2019`) against genome-wide canine somatic RNA and 38-tissue proteomics shows why quantitative somatic gates are essential: **`PRAME` (`$dog_prame_call`, `CTA-M`)** is validated as a genuine canine CTA (`110.7 TPM` in spontaneous canine lymphoma, `100.8 TPM` in testis, `0.97 TPM` somatic baseline), `MAGEB10` and `SSX6` qualify as `Reproductive Antigen only`, `CAGE1`, `CTCFL`, and `DDX53` remain held in `Other:Review` on published somatic RT-PCR counterevidence, and **6 of 14 nominations (`AKR1E2`, `IGF2BP3`, `NT5C1B`, `PAGE3-like`, `SPECC1`, `TPX2`) are excluded as `Unsafe`** due to substantial healthy canine blood, bone marrow, lung, brain, or liver expression.

## *Mus musculus* (laboratory mouse · taxon 10090)

Preclinical vaccine, TCR, and bispecific studies in mice require host-proteome peptide subtraction and endogenous normal-tissue screening on mouse-native coordinates rather than human-to-mouse symbol projection—especially across X-linked CTA families (`Magea`, `Mageb`, `Xlr`, `Trap1a`, `Pramef`) where copy number and tissue regulation diverge sharply from primates.

Murine target selection is evaluated genome-wide across **$mouse_n** Ensembl 115 `GRCm39` protein-equivalence groups (**$mouse_loci** primary-assembly coding loci, **8/8** evidence pillars active under harmonized `≤ 2.5%` ectopic `Olfr` gDNA library QC): **$mouse_cta** <span class="cat cat-cta">CTA</span>, **$mouse_ca_only** <span class="cat cat-ca">CA only</span>, **$mouse_ra_only** <span class="cat cat-ra">RA only</span> (**$mouse_union** total candidates), **$mouse_other_n** <span class="cat cat-other">Other</span>, and **$mouse_unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Across the **$mouse_ra** murine Reproductive Antigens (**$mouse_cta** CTAs), 16-compartment `HostContext` partitioning yields **$mouse_mra_n** `RA-M` (**$mouse_cta_m_n** `CTA-M`), **$mouse_era_n** `RA-P` (**$mouse_cta_p_n** `CTA-P`), **$mouse_fra_n** `RA-F` (**$mouse_cta_f_n** `CTA-F`), and **$mouse_sra_n** `RA-MF` (**$mouse_cta_mf_n** `CTA-MF`).

### Murine dataset and threshold parameterization

$mouse_datasets_table

Key murine methodological rules enforced in code (`oncoref.mouse`, `oncoref.species_target_selection`):

- **Complete pinned Ensembl 115 GRCm39 host proteome:** `mouse.references.ensure_reference("ensembl115")` verifies SHA-256-pinned GTF, cDNA, ncRNA, and protein FASTAs (**66,037** qualified primary-assembly coding translations across **21,433** `PEG`s and **21,760** coding loci), enabling full murine `PEG` grouping, 9-mer <span class="cat cat-unsafe">Unsafe</span> containment, and offline `PeptideIndex` self-proteome subtraction (`mouse-proteome-sequences`).
- **Strain-explicit tumor and normal provenance:** Reference coordinates belong to **C57BL/6J** (`GRCm39`), whereas syngeneic grafts and GEMMs span **C57BL/6** (`MC38`, `B16F10`, `LLC`, `KP_LUAD`), **BALB/c** (`CT26`, `4T1`), **FVB** (`BMPC_PRAD`), and **DBA/2** (`P815`). Cell-line replicates deduplicate donor counts per line rather than inflating independent animal counts.
- **Multi-modal normal & sorted/activated immune safety without negative-IHC inflation:** Adult somatic bulk RNA (`Bgee 15.2` + `GSE219045` across 24 required and additional somatic contexts after `≤ 2.5%` `Olfr` gDNA library QC), sorted & activated immune populations (`ImmGen` ULI `GSE109125`, `ImmGen` MNP `GSE122108`, and WT `LPS`/`IL-4` primary BMDMs `GSE181223`, summarized via finite-sample replicated `min(P95, X_{(n-1)})` across `≥2`-replicate populations), `Tabula Muris` 10x droplet UMI (`assess_rare_cell_rna`), and `MGI GXD 6.27` adult `TS28` immunohistochemistry/immunoblot detections gate safety in their native units with zero negative-IHC RNA ceiling multipliers. Mapped-universe zero-count completion on pre-filtered GEO reproductive count matrices (`GSE287166` `coagulating_gland`, `GSE324827` `fallopian tube`, `GSE60450` mammary) resolves all 16 murine `HostContext` compartments (`Sbp`, `Pbsn`, `Tex19.1`, `Plac1`, `Zfp1010` → `CTA-M`; `Serpinb9e`, `Pramel13` → `CTA-P`; `Gzmg`, `Gzmd`, `Gzmf` → `CTA-F`; `Hmga2`, `Gzme` → `CTA-MF`).

### Murine syngeneic and GEMM tumor cohorts

$mouse_cohorts_table

### Top murine CTAs by peak tumor expression

$mouse_top_cta_table

### Historical murine nomination check (`P1A` / `Trap1a`)

$mouse_nominations_table

`P1A` (`Trap1a`, `ENSMUSG00000051257` / `ENSMUST00000054534.4` / `ENSMUSP00000059448.4`, 224 aa) illustrates why OncoRef separates nomination provenance from quantitative admission: the DBA/2 `P815` rejection antigen (`M36387.1` / `AAA40475.1`, Van den Eynde et al. 1991) encodes an identical 224-aa protein to C57BL/6J `ENSMUSP00000059448.4` and reaches high expression in `P815` mastocytoma, while endogenous somatic counterevidence in BALB/c lung, spleen, liver, and thymus (Sarma et al. 1999) plus quantitative multi-organ screening classify `Trap1a` as `$mouse_p1a_call` rather than admitting it to `CTA`.

## Cross-species CTA orthology across human, dog, and mouse

With independent genome-wide catalogs evaluated natively across **human (`$cta` CTAs)**, **dog (`$dog_cta` CTAs)**, and **mouse (`$mouse_cta` CTAs)**, OncoRef quantifies directly how often a `CTA` in one mammal maps to a `CTA` in another ([download full `$ortholog_total_keys_n`-family ortholog concordance table](antigen-whitepaper/cross-species-orthologs.csv)):

$ortholog_summary_table

Why do only **$conserved_cta_pairs_n** of **$ortholog_total_keys_n** mammalian `CTA` ortholog families independently qualify as `CTA` in two or more species (and **$conserved_cta_or_ra_n** as conserved `CTA` or `RA-only`)? Three evolutionary and experimental mechanisms govern cross-species mapping:

1. **Rapid mammalian lineage amplification of CT-X and ZGA arrays:** Nearly half of human CTAs—including `NY-ESO-1` (`CTAG1A/B`), `SPANXA/B/C/D/N`, `VCX/VCY`, `GAGE*`, `XAGE1A/B`, `CT45A*`, `CT47A*`, `DUX4` (`D4Z4`), `PRAMEF1–26`, `MBD3L2–5`, and `SSU72L1–5`—originated or amplified in primates after boreoeutherian divergence and have **no 1:1 coding ortholog** in dog or mouse (`human_cta_with_dog_orth_n = $human_cta_with_dog_orth_n / $cta`, `human_cta_with_mouse_orth_n = $human_cta_with_mouse_orth_n / $cta`), while rodents independently expanded `Pramel*`, `Mageb*`, uterine NK `Gzmd/e/f/g` (`CTA-F` in mouse), and `Sbp`/`Pbsn` (`CTA-M` in mouse).
2. **Conserved boreoeutherian germline/lactation core (`$conserved_cta_or_ra_n` families):** Single-copy autosomal germline, placental, and milk genes conserve reproductive restriction across 95 million years of mammalian evolution—yielding **conserved dual-species CTAs** (`CSN2` [`CTA-F` in human `BRCA` & canine `CMT`], `PRAME` [`CTA` in human & canine `LYM`], `PSMA8`, `ANTXRL`, `LUZP4`, `TRIM42`, `PASD1` [`CTA-M`/`CTA-MF` in human & dog], plus `PLAC1`/`Plac1` and `TEX19`/`Tex19.1` [`CTA-M` in both human & mouse]), alongside conserved `RA`s (`DPPA5`, `MAGEB18`, `CYLC1`, `ACTL8`, `ADAM2`, `AKAP4`) whose tumor reactivation depends on cohort breadth (`132` human cohorts vs. `8` canine patient cohorts vs. `8` murine models).
3. **Divergent somatic regulation and single-cell resolution:** Several conserved 1:1 orthologs (`HORMAD1`, `DDX43`, `CLDN6`, `IGF2BP1`) qualify in canine bulk RNA + 38-tissue proteomics while failing human single-cell cluster, adult cornea, or sorted-immune gates—underscoring why canine single-cell UMI coverage (`canine_normal_scrna_gap`) is tracked explicitly and why target safety is never inferred across species.

<details markdown="1">
<summary>Conserved multi-species CTA and RA-ortholog concordance table ($conserved_cta_or_ra_n families; download full CSV)</summary>

Download all `$ortholog_total_keys_n` ortholog families: [`cross-species-orthologs.csv`](antigen-whitepaper/cross-species-orthologs.csv).

$conserved_ortholog_ctas_table

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
cta_female_safe = antigen_sets.select("CTA-M")  # $cta_m_n male-restricted female-safe human CTAs
ra_female_only = antigen_sets.select("RA-F")  # $ra_f_n female-restricted male-safe human RAs
pregnancy_only = antigen_sets.select("CTA-P")  # $cta_p_n adult-absent pregnancy/peripartum CTAs
dog_ctas = antigen_sets.select("CTA", species="dog")  # $dog_cta spontaneous canine tumor CTAs
mouse_ctas = antigen_sets.select("CTA", species="mouse")  # $mouse_cta murine GRCm39 CTAs

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
