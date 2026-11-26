# Cancer and Reproductive Antigens

**CTA = CA ∩ RA. Candidate antigen proteins = CA ∪ RA.** Cancer-only proteins remain
Cancer Antigens. These are research screening memberships, not demonstrated antigen
presentation or clinical safety.

Catalog `$policy_id`: $categories.

[Antigen whitepaper](antigen-whitepaper.md) ·
[Interactive protein catalog](audits/target-selection-20261110/rebuild/report.html) ·
[Summary figures](audits/target-selection-20261110/rebuild/summary.md) ·
[Patient-anatomy rules](host-selection.md) · [Patient profiles & filters](target-filters.md)

## Three separate questions

1. **Safety & eligibility:** Does the protein pass the normal-expression, coverage,
   identity, translation, sequence, and expression-bounds checks?
2. **Biological enrichment:** Does tumor expression qualify (`CA`), reproductive or
   developmental expression qualify (`RA`), or both (`CTA`)?
3. **Patient-specific filtering:** Which reproductive or physiological compartments must be
   protected for a given patient scenario, and is coverage sufficient to assess them?

![Catalog and host selection](audits/target-selection-20261110/rebuild/summary/source-to-tiers-architecture.png)

## One threshold, two compartments

Both CA and RA independently require **≥10 × H₀**, where
**H₀ = max(non-thymic, non-reproductive normal N, Disposable D, 1 TPM)**, after shared eligibility.

| Set | Positive RNA evidence | Proteins |
|---|---|---:|
| CA | At least one QC-passing tumor sample ≥10 × H₀ | $ca |
| RA | Reproductive sample or supported embryo-lineage mean ≥10 × H₀ | $ra |
| CTA | Both routes | $cta |
| CA∪RA | Either route | $union |

The adult whole thymus is Disposable under the catalog’s loss allowance. Defined immune cell populations retain their own injury tiers. A matched observational study associated adult thymectomy with increased mortality and cancer incidence. [Adult thymectomy study](https://pubmed.ncbi.nlm.nih.gov/37530823/) Medullary thymic epithelial cells express tissue-restricted antigens during central tolerance. Bulk thymus is retained in the healthy ledger under the `thymic` scope and excluded from both enrichment denominators. The independent normal-expression and protein screens still apply. [Human thymic epithelial-cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2211762/)

GTEx v11 contributes the **linearly interpolated donor P95 within each non-reproductive tissue site**. Complete member-locus TPMs are summed within a library first. Reproductive sites, including breast and prostate, retain individual samples as RA witnesses. HPA v23 contributes the within-tissue biopsy P95 and the within-population sorted immune donor P95. DICE contributes a P95 for single loci. Its ranks lack specimen IDs, so identical-protein groups with multiple loci use the sum of measured member maxima as an upper bound. GTEx medians, HPA tissue and brain means, and Monaco cell means retain their published statistics. The highest applicable value across these sources sets N or D. Single-cell and single-nucleus counts, nTPM and protein staining have separate normal-expression gates.

[All RNA sources, tissue coverage and set definitions](audits/target-selection-20261110/rebuild/normal-rna/index.md)

![Normal and reproductive TPM sources by tissue](audits/target-selection-20261110/rebuild/normal-rna/tpm-coverage.png)

A small cohort can have a P95 close to its maximum. The underlying measurements and sample QC remain available. [Download the source data and literature controls](antigen-data.md).

Missing required normal coverage does not become zero.
[Source receipt](audits/target-selection-20261110/rebuild/comparison.json).

Developmental means are calculated within embryo and lineage. A well-sampled
population requires ≥20 measured cells, ≥5 cells with ≥1 TPM and ≥10% positive cells.
Smaller populations require ≥5 measured cells, ≥50% positive cells and two independent
embryos of the same age and lineage reaching the tenfold threshold. Cells or lineages
from one embryo are not independent donors. These positive-support requirements do
**not** determine whether an observation can inform normal host exposure.

## Shared eligibility and Other

| Evidence or gap | Classification |
|---|---|
| Accepted normal RNA, normal protein or applicable clinical exclusion | Unsafe |
| Missing coverage, uncertain identity/translation or unresolved normal finding | Other:Review |
| Unsupported nested protein, ≥95% Unsafe 9-mer containment or read-through locus | Other:Sequence |
| Unmeasured home tissue without qualifying reproductive evidence | Other:Unscreened tissue |
| Neither enrichment route qualifies, or expression bounds fail | Other:Expression |

Insufficient coverage of a documented normal home blocks the CA route even when RA qualifies. Unsafe takes precedence, then Other:Review. Every failed gate remains on the row.
**Other describes unresolved or nonqualifying proteins, not absent normal expression.**
Disposable expression has separate 32 TPM and 32 nTPM bounds in the catalog; these
are not the limits for host-protected organs.

## Rare normal-cell RNA

Single-cell and single-nucleus measurements can identify populations diluted in bulk
RNA. Donor-resolved UMI evidence is graded independently of cluster nTPM and bulk TPM.

$rare_table

A mapping-risk flag is not proof of read misassignment. Strong unresolved signals
remain Other:Review; low bulk expression alone cannot erase them.
[Donor records and rationale](audits/normal-cell-atlas-20261110/rare-cell/index.md).

### Nuclear snRNA-seq compartment & transcriptional spillover (`#765`)

Single-nucleus RNA-seq (`snRNA-seq`, `assay_compartment == "nucleus"`) captures unspliced pre-mRNA across long non-coding 5'-UTR introns (`PLAC1`: `197.7 kb` 5'-UTR span, `0.32%` CDS), same-strand termination readthrough downstream of active somatic genes (`PABIR2 → PLAC1`, `EPO → ZAN`, `ERCC5 → METTL21EP`), and divergent promoter-upstream transcripts (`NOP10 ↔ NUTM1` across a `135 bp` head-to-head CpG island at `4,501 TPM` in neutrophils). Because pMHC presentation requires cytoplasmic translation of a coding reading frame, a replicated rare-cell flag clears to `rare_cell_status = "nuclear_only"` (`nuclear_rna_caution = "snrnaseq_intronic_spillover_cleared"`) **only** when all four guardrails hold simultaneously:

1. **Zero cytoplasmic single-cell support:** `nucleus_only_signal == True` across all adult whole-cell `scRNA-seq` collections.
2. **Strict bulk-dilution ceiling (`MAGEA12` regression anchor):** `finding_bulk_tpm < 0.05 TPM` in the flagged organ. A genuine `30–50 TPM` rare somatic population in `<1%` of cells dilutes to `0.1–0.5 TPM` in whole-organ polyA bulk RNA (`MAGEA12` in brain pericytes measures `0.30 TPM` with `confirm_detection` literature and remains strictly `Unsafe`).
3. **Negative somatic IHC & literature:** zero positive somatic antibody staining and no `exclude` or `confirm_detection` normal-expression literature record.
4. **Verified Ensembl 115 spillover geometry:** `< 1 kb` bidirectional promoter (`≥ 10 TPM` partner), `< 10 kb` same-strand upstream partner (`≥ 5 TPM`), or `≥ 50 kb` non-coding 5'-UTR span (`CDS/span ≤ 5%`). Single-CDS-exon / intronless coding genes (`DUX4`, `PLAC1`, `PRM1`) record `is_intronless_cds = True` and accept `cds_exonic_reads` rather than requiring exon–exon splice junctions.

Cleared loci that lack unresolved paralog ambiguity enter normal native-TPM evaluation directly (`PLAC1 → CTA`, `ZAN → CTA`, `TENT5D → RA only`). [Full nuclear spillover audit](audits/nuclear-rna-spillover-20261123/index.md).

## Stoichiometric multiplex index-hopping vs. boundary carryover (`#764`, `#776`)

On patterned Illumina flowcells, hyper-abundant reproductive transcripts (`≥ 1,000 TPM`, such as `PRM1` at `58,751 TPM` and `PRM2` at `50,821 TPM` in testis or `KLK3`/`MSMB`/`ACP3` in prostate) hop into co-sequenced somatic and opposite-sex GTEx libraries at `10⁻⁵` to `5 × 10⁻³` of source abundance. `oncoref.normal_library_qc` flags a non-source GTEx biopsy for cross-sample bleed **only** when `≥ 3` source sentinels co-occur at `≥ 0.25 TPM` in stoichiometric proportion (`IQR(log₁₀ α) ≤ 0.5`) **and** orthogonal HPA v23 single-cell clusters, HPA normal IHC, and non-GTEx HPA bulk biopsies confirm absence in that organ. Flagged biopsies are **dropped from the `gtex_donor_p95` denominator** rather than zero-filled (cleaning `PRM1` Thyroid bleed from `8.05` to `5.71 TPM`), whereas anatomic boundary carryover in adjacent pelvic/urogenital organs (`SEMG2` and `TGM4` in `Bladder`) is distinguished and never masked in somatic `N`. In downstream `host-exposure` (`#776`, `#779`), GTEx `Breast_Mammary_Tissue` (`336/514` XY male cadaver biopsies vs. `178` XX female biopsies) is stratified by biological donor sex, reproductive age (`GTEx AGE 20–49` across female reproductive & mammary organs), and secretory state (`female_max_tpm` vs. `male_max_tpm`), adult human folliculogenesis/oocyte expression (`GSE107746`, `128` QC-passing stages) is integrated into `ovary`, GTEx-v11-missing `vagina` and `uterus` measurements (`TSSK2`) are backfilled from pinned `rna_tissue_gtex:v23` and `hpa_tissue_samples:v23`, and `classify_sex_restriction` partitions all `RA`s and `CTA`s into four mutually exclusive physiological classes—`RA-M` / `CTA-M` (`226 CTA-M`), `RA-F` / `CTA-F` (`12 CTA-F`), `RA-MF` / `CTA-MF` (`38 CTA-MF`), and `RA-P` / `CTA-P` (`49 CTA-P`, pregnancy/peripartum-only across `embryo`, `placenta`, and `lactating_breast`; legacy aliases `CTA-E`/`RA-E`/`ERA`/`PRA`). [Multiplex bleed QC audit](audits/multiplex-bleed-qc-20261123/index.md).

## Isoform- & epitope-scoped clinical pMHC organ clearance (`#766`)

Human adoptive TCR-T trials (`lete-cel` targeting `HLA-A*02:01 / SLLMWITQC` on `NY-ESO-1` and `CTAG2/LAGE-1a`; `afami-cel` on `MAGEA4`; `IMA203` on `PRAME`) directly test whether vascularly perfused organs present lethal levels of a specific pMHC complex *in vivo*, while fatal clinical toxicities (`MAGEA3` titin cardiotoxicity, `MAGEA3/A12` necrotizing leukoencephalopathy) remain permanent `Unsafe` negative controls. Clinical clearance (`oncoref.clinical_organ_clearance`) is strictly **epitope-, HLA-, and Ensembl-115-isoform-scoped**:

- **Isoform specificity (`CTAG2`):** `ENST00000369585` (`LAGE-1a`, `180 aa`) encodes `SLLMWITQC` (`aa 157–165`), whereas `ENST00000247306` (`LAGE-1b`, `210 aa`) retains intron-2 sequence and **lacks `SLLMWITQC`**. `CTAG2` advances from `Other:Review → CTA` (`604.3 TPM` in `OV` vs. `6.92 TPM` GTEx Urinary Bladder baseline, `87.3×`), while `clinical_clearance_scope = "epitope"`, `clinical_cleared_transcripts = ["ENST00000369585"]`, and `clinical_unverified_isoforms = ["ENST00000247306"]` warn downstream epitope selectors (`tsarina`) not to generalize lung/bladder safety to unverified `LAGE-1b` peptides or other HLA alleles.
- **Barrier & modality guardrails:** Systemic IV TCR-T never clears blood–brain/retina/testis barrier-protected organs (`brain` on `CTAG2` is resolved separately by `#711` Xq28 31-mer attribution) and never clears `heart` when confounded cardiac events exist; peptide/protein vaccines cannot grant organ clearance. [Clinical organ-clearance audit](audits/clinical-organ-clearance-20261123/index.md).

## Closed homologous repeat-family sequence equivalence (`#767`)

Tandem repeat arrays and recently duplicated cancer–germline gene clusters (`DUX4` across `D4Z4` on `4q35`/`10q26`, `CT47A*`, `TSPY*`, `SPANX*`, `CSAG*`, `CTAG*`) trigger `paralog_ambiguous = True` in single-gene aligners even when every sequence-sharing copy in the human genome is equally silent across healthy adult somatic organs. `oncoref.repeat_family_equivalence` partitions sequence partners of each paralog-flagged locus `g` into:

- **`P_self(g)` (closed self-family homologs):** primary-assembly `protein_coding` and `unprocessed_pseudogene` copies sharing `≥ 85%` CDS 31-mers (or `≥ 90%` amino-acid identity) with zero `Unsafe` somatic expression (`unsafe_sequence_overlap_pct < 5%` and family-summed somatic donor `P95 < 5.0 TPM`).
- **`P_ext(g)` (external / somatic paralogs):** `Unsafe` coding paralogs (`PSG8 ↔ PSG1`, `POTEC ↔ POTED/E`, `GAGE2E ↔ GAGE12F/G`), `processed_pseudogene` retrocopies, or chimeric somatic hosts, which remain blocked in `Other:Review`.

When `P_ext(g)` is empty, symmetric family-summed expression across `{g} ∪ P_self(g)` tests somatic silence honestly (`DUX4` + `20` `DUX4L*` repeat units = `0.00 TPM` in GTEx Skeletal Muscle and Skin, `0.30 TPM` somatic max, advancing `DUX4 → CTA` [`CTA-P`] jointly with `#765`). Catalog rows also expose `closed_homolog_family_id`, `co_silenced_coding_sibling_ids`, and `closed_self_paralog_ids` so downstream peptide subtraction (`tsarina`) does not self-subtract shared epitopes against co-silenced germline siblings (`CT47C1`, `TSPY9`, `SPANXN2/5`, `CSAG2/3`). [Closed repeat-family audit](audits/repeat-family-equivalence-20261123/index.md).

## Protein staining & assay-boundary reconciliations (`#768`)

Identity-resolved, reliable normal protein evidence can independently mark a protein
Unsafe under the tissue-tier policy. Antibody reliability, specificity, cell type and
source scope remain part of each finding. Four explicit assay-boundary reconciliations (`oncoref.assay_boundary_adjudications`) prevent single-assay artifacts from misclassifying targets:

1. **Non-admixed tumor witness requirement (`#713` Note `R` `TCGA-TGCT`):** `11` post-meiotic spermatid markers (`PRM2`, `CETN1`, `PGK2`, `BOD1L2`, `CCDC70`, `TSSK2`, `PRR30`, `ZNRF4`, `IQCF5`, `BPIFA3`, `CYLC2`) whose `≥ 10×` cancer signal rested solely on `library_qc` Note `R` seminiferous-parenchyma-admixed `TCGA-TGCT` biopsies reclassify to `Reproductive Antigen only`, whereas `PRM1` (`58.97 TPM` in female `BRCA`) and `SPATA19` (`13.28 TPM` in `#737` `cSCC`) retain independent somatic tumor witnesses and remain `CTA`.
2. **RNA-zero single-antibody HPA IHC (`cta-ihc-unreliable.csv`):** Single-antibody staining contradicted by `100%` complete required somatic RNA silence (`< 0.5 TPM` and `< 1.0 nTPM`, e.g. `MAGEB6` at `0.012 TPM` somatic max) clears to `CTA`, whereas loci with somatic RNA (`NLRP9` at `2.32 TPM`) remain `Unsafe`.
3. **Adult cornea (`GSE155683`) multi-UMI support floor:** Requires `≥ 10` adult corneal cells with `≥ 2` UMIs, matching the rare-cell screen's minimum support (`DKKL1` with `8 / 21,343` cells advances `Unsafe → CTA`; `INSL4`, `GPHB5`, `QRFP`, `C20orf202` with `24–64` cells stay `Unsafe`).
4. **Protamine-diluted germline transcription factors:** Nuclear germ-cell regulators compressed below `10 TPM` in whole-testis bulk RNA by `> 100,000 TPM` late-spermatid protamines qualify on `RA` when somatic normal expression is silent (`≤ 0.10 TPM`, `0.0 nTPM` somatic scRNA), germline bulk reaches `≥ 5.0 TPM` (`≥ 50×` over somatic) across `≥ 2` GTEx Testis donors `≥ 3.0 TPM`, and HPA germ-cell `scRNA-seq` confirms germline lineage expression (`≥ 1.0 nTPM`), advancing `FOXR2`, `SPATA31A3`, and `PPP4R3C` from `CA only → CTA` (`CTA-M`). [Assay-boundary & 1-mm near-self audit](audits/assay-boundary-near-self-20261123/index.md).

## Published normal expression

Reviewed gene-specific normal expression
retains sources, adjudication and supersession history.
[Literature rubric](audits/target-selection-20261110/rebuild/literature.md).

## Home tissue the panels do not resolve

Olfactory epithelium, nasal mucosa, taste buds, hair follicle, dental pulp and mineralized
bone are not resolved by the standard GTEx/HPA bulk native-TPM panels. `oncoref.unmeasured_home_panels` integrates SHA-256-pinned adult human bulk TPM and donor-resolved single-cell UMI panels across these six tissues (`#735`): proteins exceeding their home-tissue injury-tier ceiling (`SP7` and `PANX3` in bone, `GFY` and `BPIFB3` in olfactory epithelium, `GLYATL3`) become `Unsafe`, proteins clearing adult home panels or Disposable follicle bounds (`33` proteins in `CA only` under clean tumor QC, incl. `TAS1R2`, `KRTAP20-3`, `KRTAP21-2`, `KRTAP21-3`, `KRTAP6-2`, plus `OR10J1` and `KRTAP20-4` promoted to `CTA`) enter the candidate sets, and proteins whose sub-niche remains unmeasured (`DSPP` odontoblast process, `CALHM3` type II taste cells) stay in `Other:Unscreened tissue`. A protein with a documented home in an unmeasured niche can still qualify as `RA`.

## Tumor quality and replication

Library QC enforces the olfactory-receptor genomic-DNA background ceiling (`≤ 2.5%` of measured olfactory receptors at `≥ 5 TPM`, newly excluding `92` baseline-admitted gDNA-contaminated tumor libraries across `NEC_MERKEL`, `SARC_PEC`, `BL`, `MCL`, `HL`, `B_ALL`, `LAML`, and `T_ALL`) and unrelated-tissue screens.
Coordinated late-spermatid RNA is an **advisory note only**: normal reproductive
specificity does not establish absence from cancer or prove normal-cell admixture.
These notes never change sample QC status or remove tumor or reproductive evidence,
including when selecting only `sample_qc="pass"`. Independent QC exclusions still apply.
The [library audit](audits/target-selection-20261110/rebuild/library-qc.csv.gz) separates
blocking `reasons` from advisory `notes`; the public expression QC API exposes
`sample_qc_reasons` and `sample_qc_notes` separately.
$library_note_legend
Cohort availability and exclusions remain explicit
in the [source coverage table](audits/target-selection-20261110/rebuild/tumor-source-coverage.csv).
Catalog rows expose samples, independent patients, anonymous samples and cancer types
at tenfold. `independent_patient_tumor_folds` holds descending patient maxima; multiple
samples from one known patient count once. Unknown donors do not count as known patients.

## Biological origin and patient-specific protection

`reproductive_qualifying_compartments` records the reproductive or embryonic compartments reaching ≥10×. Protect, Accept loss, and Absent are independent downstream choices: anatomy presets protect no reproductive tissues by default, and filtering a candidate out for a patient scenario leaves its biological `CA` / `RA` / `CTA` membership unchanged while adding it to that selection's peptide-subtraction background.
[Patient-anatomy rules and limits](host-selection.md).

![Host choices and outcomes](audits/target-selection-20261110/rebuild/summary/cta-tiers-and-profiles.png)

## One selection API and downstream `tsarina` columns

```python
from oncoref import antigen_sets, target_selection

selection = antigen_sets.select("CTA", mode="standard")
proteins = selection.proteins()
genes = selection.gene_ids()
sequences = selection.sequences()  # offline sequence_bundle or local Ensembl 115
background = selection.background_index("Unsafe|Other")

# Inspect structured clinical clearance, closed-family sibling, and 1-mm near-self columns
cols = [
    "protein_group_label",
    "classification",
    "sex_restriction_code",
    "clinical_clearance_scope",
    "clinical_cleared_transcripts",
    "closed_homolog_family_id",
    "co_silenced_coding_sibling_ids",
    "exact_unique_9mers",
    "safe_2mm_9mers",
    "unsafe_near_self_1mm_pct",
]
target_selection.results()[cols].query(
    "protein_group_label in ['CTAG2', 'DUX4', 'PLAC1', 'DKKL1', 'MAGEB6', 'MAGEA6']"
)
```

[Sequence and peptide API](antigen-sequences.md) · [CLI and rebuild](curation.md).
The [rebuild](audits/target-selection-20261110/rebuild/index.md) records every input,
output and decision. Outstanding scientific adjudications remain explicit holds;
reproductive-content notes do not resolve malignant-cell localization.
