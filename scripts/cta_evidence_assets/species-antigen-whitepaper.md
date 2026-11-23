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

- **`MRA` (`CTA-M` / `MRA-only`, Male-restricted · female-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every female adult reproductive and mammary compartment of species *S* (`< 4` native RNA units with measured coverage, maintaining ≥10× tumor margin over protected host organs for `CTA-M`), while expressed in male organs.
- **`FRA` (`CTA-F` / `FRA-only`, Female-restricted · male-safe):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible when protecting every male adult reproductive compartment plus resting mammary tissue of species *S*, while expressed in female or placental compartments.
- **`ERA` (`CTA-E` / `ERA-only`, Embryo / placenta-only · adult-absent in both sexes):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) eligible under full adult organ protection in **both** sexes—expressed in embryo or placenta and silent across every adult male and non-pregnant female reproductive and mammary compartment.
- **`SRA` (`CTA-MF` / `SRA-only`, Shared male & female):** <span class="cat cat-ra">RA</span> (or <span class="cat cat-cta">CTA</span>) expressed in adult reproductive/mammary organs of both sexes (or lacking measured opposite-sex clearance, failing closed).

$comparative_anatomy_table

Comparative anatomy changes both enrichment denominators and host protection:

- **Accessory male glands:** Dogs anatomically lack seminal vesicles and bulbourethral glands (`prostate` is the sole major accessory gland), whereas mice carry a distinct `coagulating_gland` (anterior prostate) alongside `seminal vesicle` and `prostate`. Requesting `seminal vesicle` in a canine `HostContext` raises an immediate `ValueError` rather than silently dropping the constraint.
- **Bicornuate female tract & gonadectomy:** Both dogs and mice add `uterine_horn` (normalized automatically under `uterus`), and companion-animal oncology frequently treats spayed (ovariohysterectomized) or castrated hosts where `HostContext(absent=("ovary", "uterus", "uterine_horn", "endometrium"), species="dog")` explicitly records surgically absent organs distinct from accepted tissue injury.
- **Disposable somatic organs (`D`):** Neither dogs nor mice possess a vermiform `appendix`, and mice also lack palatine `tonsil`, narrowing non-thymic `Disposable` (`D`) from five compartments in human to four in dog and three in mouse.

## Eight-pillar cross-species parameterization matrix

The table below aligns all eight canonical screening roles across human, dog, and mouse ([download full binding table](antigen-whitepaper/species-datasets.csv)):

$eight_pillar_matrix

## *Homo sapiens* (human · taxon 9606)

Human target selection is fully evaluated across **$n** protein-equivalence groups (**$loci** Ensembl 115 coding loci) under policy `$policy`: **$cta** <span class="cat cat-cta">CTA</span>, **$ca_only** <span class="cat cat-ca">CA only</span>, **$ra_only** <span class="cat cat-ra">RA only</span> (**$union** total candidates), **$other_n** <span class="cat cat-other">Other</span>, and **$unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Across the **$ra** Reproductive Antigens (**$cta** CTAs), sex-restriction partitioning yields **$mra_n** `MRA` (**$cta_m_n** `CTA-M`), **$era_n** `ERA` (**$cta_e_n** `CTA-E`), **$fra_n** `FRA` (**$cta_f_n** `CTA-F`), and **$sra_n** `SRA` (**$cta_mf_n** `CTA-MF`). See the [main antigen whitepaper](antigen-whitepaper.md) for full human figures, clinical anchor tables, and edge-case audits.

$human_datasets_table

## *Canis lupus familiaris* (domestic dog · taxon 9615)

Spontaneous canine cancers—osteosarcoma, diffuse large B-cell lymphoma, urothelial carcinoma, mammary carcinoma, oral melanoma, and glioma—develop in outbred, immunocompetent hosts sharing human environmental exposures. OncoRef ships pinned canine genome references, native normal RNA and proteomics readers, comparative canine `HostContext` profiles, and **14** literature-curated canine CTA nominations (`hindriks2025`, `nemec2019`).

### Canine dataset and threshold parameterization

$dog_datasets_table

Key canine methodological rules enforced in code (`oncoref.canine`, `oncoref.species_target_selection`):

- **Single-scale enrichment panels (`TPM` vs. `TMM-CPM`):** Bgee (`bgee-srp009687`, `bgee-srp181218`) quantifies canine libraries in `TPM`, whereas BarkBase (`barkbase`) reports StringTie genes in `TMM-CPM` and shares **46 SRA runs** with the global Bgee inventory (`canine.expression.lineage_overlaps()`). `SpeciesScreeningSpec` accepts either a pure `TPM` panel or a pure `TMM-CPM` panel, forbidding cross-unit pooling or duplicate donor counting across shared runs.
- **Pooled proteomics (`pfizer-dog-protein`):** MassIVE `MSV000099744` provides quantitative `iBAQ` across **38 pooled canine tissues** (**9,929** protein groups) plus optional peptide-to-protein assignments (`pfizer-dog-peptide-evidence`). Because every channel pools multiple animals without cell-type resolution, `assess_pooled_protein_gates` uses high Vital/Protected `iBAQ` only as affirmative <span class="cat cat-unsafe">Unsafe</span> / hold evidence once searched-FASTA accessions map to canine genes, and never awards negative-IHC RNA ceiling multipliers.
- **Explicit data gaps blocking packaged genome-wide calls:** Source quantification GTF crosswalks for Bgee/BarkBase and raw-count versus FPKM provenance in `gse304066-osteosarcoma` remain unresolved, so packaged canine nominations stay `Other:Review` rather than inheriting human CTA labels.

### Packaged canine nomination ledger (14 native candidates)

$dog_nominations_table

## *Mus musculus* (laboratory mouse · taxon 10090)

Preclinical vaccine, TCR, and bispecific studies in mice require host-proteome peptide subtraction and endogenous normal-tissue screening on mouse-native coordinates rather than human-to-mouse symbol projection—especially across X-linked CTA families (`Magea`, `Mageb`, `Xlr`, `Trap1a`, `Pramef`) where copy number and tissue regulation diverge sharply from primates.

Murine target selection is evaluated genome-wide across **$mouse_n** Ensembl 115 `GRCm39` protein-equivalence groups (**$mouse_loci** primary-assembly coding loci, **8/8** evidence pillars active): **$mouse_cta** <span class="cat cat-cta">CTA</span>, **$mouse_ca_only** <span class="cat cat-ca">CA only</span>, **$mouse_ra_only** <span class="cat cat-ra">RA only</span> (**$mouse_union** total candidates), **$mouse_other_n** <span class="cat cat-other">Other</span>, and **$mouse_unsafe_n** <span class="cat cat-unsafe">Unsafe</span>. Across the **$mouse_ra** murine Reproductive Antigens (**$mouse_cta** CTAs), 16-compartment `HostContext` sex-restriction partitioning yields **$mouse_mra_n** `MRA` (**$mouse_cta_m_n** `CTA-M`), **$mouse_era_n** `ERA` (**$mouse_cta_e_n** `CTA-E`), **$mouse_fra_n** `FRA` (**$mouse_cta_f_n** `CTA-F`), and **$mouse_sra_n** `SRA` (**$mouse_cta_mf_n** `CTA-MF`).

### Murine dataset and threshold parameterization

$mouse_datasets_table

Key murine methodological rules enforced in code (`oncoref.mouse`, `oncoref.species_target_selection`):

- **Complete pinned Ensembl 115 GRCm39 host proteome:** `mouse.references.ensure_reference("ensembl115")` verifies SHA-256-pinned GTF, cDNA, ncRNA, and protein FASTAs (**66,037** qualified primary-assembly coding translations across **21,433** `PEG`s and **21,760** coding loci), enabling full murine `PEG` grouping, 9-mer <span class="cat cat-unsafe">Unsafe</span> containment, and offline `PeptideIndex` self-proteome subtraction (`mouse-proteome-sequences`).
- **Strain-explicit tumor and normal provenance:** Reference coordinates belong to **C57BL/6J** (`GRCm39`), whereas syngeneic grafts and GEMMs span **C57BL/6** (`MC38`, `B16F10`, `LLC`, `KP_LUAD`), **BALB/c** (`CT26`, `4T1`), **FVB** (`BMPC_PRAD`), and **DBA/2** (`P815`). Cell-line replicates deduplicate donor counts per line rather than inflating independent animal counts.
- **Multi-modal normal safety without negative-IHC inflation:** Adult somatic bulk RNA (`Bgee 15.2` + `GSE219045` across 24 required and additional somatic contexts), `Tabula Muris` 10x droplet UMI (`assess_rare_cell_rna`), and `MGI GXD 6.27` adult `TS28` immunohistochemistry/immunoblot detections gate safety in their native units with zero negative-IHC RNA ceiling multipliers.

### Murine syngeneic and GEMM tumor cohorts

$mouse_cohorts_table

### Top murine CTAs by peak tumor expression

$mouse_top_cta_table

### Historical murine nomination check (`P1A` / `Trap1a`)

$mouse_nominations_table

`P1A` (`Trap1a`, `ENSMUSG00000051257` / `ENSMUST00000054534.4` / `ENSMUSP00000059448.4`, 224 aa) illustrates why OncoRef separates nomination provenance from quantitative admission: the DBA/2 `P815` rejection antigen (`M36387.1` / `AAA40475.1`, Van den Eynde et al. 1991) encodes an identical 224-aa protein to C57BL/6J `ENSMUSP00000059448.4` and reaches high expression in `P815` mastocytoma, while endogenous somatic counterevidence in BALB/c lung, spleen, liver, and thymus (Sarma et al. 1999) plus quantitative multi-organ screening classify `Trap1a` as `$mouse_p1a_call` rather than admitting it to `CTA`.

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
cta_female_safe = antigen_sets.select("CTA-M")          # 201 female-safe human CTAs
era_targets = antigen_sets.select("CTA-E")              # 20 adult-absent embryo/placenta CTAs
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
