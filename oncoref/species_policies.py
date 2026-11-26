# Licensed under the Apache License, Version 2.0.
"""Declarative per-species parameterization of the OncoRef antigen-selection framework.

Human (9606), dog (9615) and mouse (10090) share one mathematical and anatomical
definition of Cancer Antigens (CA), Reproductive Antigens (RA), Cancer-Testis
Antigens (CTA = CA & RA), and the 4-way physiological sex-restriction partition
(RA-M/CTA-M, RA-F/CTA-F, RA-P/CTA-P, RA-MF/CTA-MF; legacy aliases MRA/FRA/ERA/PRA/SRA).
Each species parameterizes that invariant framework with:

1. Its own verified genome/proteome ReferenceIdentity registrations.
2. Its own comparative mammalian anatomy from ``normal-tissue-taxonomy.yaml``.
3. Explicit per-unit quantitative thresholds (never converting or pooling scales).
4. Species-native datasets bound to eight canonical evidence roles, with honest
   integration status where analogous atlases are incomplete or pending crosswalks.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import cache

import pandas as pd

from . import cta_selection_policy, target_selection
from .reference_identity import ReferenceIdentity
from .species import Species, available_species, resolve_species
from .species_assets import registered_references
from .species_target_selection import (
    SpeciesAnatomyProfile,
    SpeciesScreeningSpec,
    anatomy_for_species,
)

POLICY_SCHEMA_ID = "oncoref.species-screening-policy.2026-11-21.1"

EVIDENCE_ROLES: tuple[tuple[str, str], ...] = (
    (
        "reference_proteome",
        "Reference proteome & identical-sequence protein-equivalence groups (PEG)",
    ),
    (
        "healthy_baseline_rna",
        "Healthy non-thymic, non-reproductive baseline (N) & Disposable maximum (D)",
    ),
    (
        "normal_bulk_and_immune_safety",
        "Normal somatic bulk & sorted immune RNA safety ceilings",
    ),
    (
        "normal_cell_and_rare_umi",
        "Single-cell cluster & donor-resolved rare-cell UMI safety screens",
    ),
    (
        "normal_protein",
        "Normal-tissue protein staining (IHC) or mass-spectrometry proteomics",
    ),
    (
        "reproductive_and_developmental_rna",
        "Adult reproductive & embryonic/placental enrichment (RA) and host exposure",
    ),
    (
        "tumor_rna",
        "QC-passing tumor cohorts for cancer enrichment (CA)",
    ),
    (
        "curated_literature_and_clinical",
        "Curated normal counterevidence, unscreened adult homes & published nominations",
    ),
)
ROLE_IDS = tuple(role_id for role_id, _ in EVIDENCE_ROLES)
ROLE_LABELS = dict(EVIDENCE_ROLES)

DATASET_STATUSES = (
    "active_in_catalog",
    "ingested_pending_crosswalk",
    "native_seed_only",
    "metadata_or_blocked",
    "identified_analog",
)


@dataclass(frozen=True)
class SpeciesThresholdConfig:
    """Explicit quantitative parameters in a species' native measurement units."""

    comparison_units: tuple[str, ...]
    normal_summary_statistic: str = "p95"
    enrichment_fold: float = target_selection.ENRICHMENT_FOLD
    detection_floor: float = target_selection.DETECTION_FLOOR
    disposable_cap: float = target_selection.DISPOSABLE_CAP
    normal_rna_tier_ceilings: dict[str, float] | None = None
    negative_ihc_multipliers: dict[str, float] | None = None
    cornea_cp10k_scale: float | None = None
    rare_cell_min_cells: int = target_selection.RARE_CELL_MIN_CELLS
    rare_cell_multi_molecules: int = target_selection.RARE_CELL_MULTI_MOLECULES
    rare_cell_vital_replicated_pct: float = target_selection.RARE_CELL_DETECTION_PCT
    rare_cell_vital_single_study_pct: float = target_selection.RARE_CELL_SINGLE_STUDY_PCT
    rare_cell_protected_factor: float = target_selection.RARE_CELL_PROTECTED_FACTOR
    rare_cell_replicate_studies: int = target_selection.RARE_CELL_REPLICATE_STUDIES
    rare_cell_replicate_donors: int = target_selection.RARE_CELL_REPLICATE_DONORS
    unsafe_ninemer_containment_pct: float = 95.0
    host_protected_max_tpm: float = 8.0
    host_protected_max_ntpm: float = 8.0
    host_negative_ihc_multiplier: float | None = None

    def __post_init__(self):
        if self.normal_rna_tier_ceilings is None:
            object.__setattr__(
                self,
                "normal_rna_tier_ceilings",
                {
                    "vital": cta_selection_policy.TIER_BASES["vital"],
                    "protected": cta_selection_policy.TIER_BASES["protected"],
                    "manageable_hold": cta_selection_policy.TIER_BASES["manageable"],
                    "manageable_unsafe": cta_selection_policy.MANAGEABLE_UNSAFE_FLOOR,
                },
            )
        if self.negative_ihc_multipliers is None:
            object.__setattr__(self, "negative_ihc_multipliers", {})
        if self.host_negative_ihc_multiplier is None:
            object.__setattr__(
                self,
                "host_negative_ihc_multiplier",
                min(self.negative_ihc_multipliers.values())
                if self.negative_ihc_multipliers
                else 1.0,
            )


@dataclass(frozen=True)
class SpeciesDatasetBinding:
    """One species-native dataset bound to one canonical screening role."""

    role: str
    source_id: str
    label: str
    native_unit: str
    status: str
    threshold_or_use: str
    limitations: str

    def __post_init__(self):
        if self.role not in ROLE_LABELS:
            raise ValueError(f"Unknown evidence role {self.role!r}; expected one of {ROLE_IDS}")
        if self.status not in DATASET_STATUSES:
            raise ValueError(
                f"Unknown dataset status {self.status!r}; expected one of {DATASET_STATUSES}"
            )


@dataclass(frozen=True)
class SpeciesPolicy:
    """Complete per-species parameterization of the OncoRef antigen-selection rules."""

    species: Species
    policy_id: str
    primary_reference_alias: str
    registered_reference_aliases: tuple[str, ...]
    anatomy: SpeciesAnatomyProfile
    thresholds: SpeciesThresholdConfig
    datasets: tuple[SpeciesDatasetBinding, ...]
    overview: str

    @property
    def primary_reference(self) -> ReferenceIdentity:
        from .species_antigens import cta_snapshot

        if self.species.taxon_id == 9606:
            return ReferenceIdentity.from_dict(
                cta_snapshot("human").as_dict()["references"][0]["identity"]
            )
        reg = next(
            r
            for r in registered_references(self.species)
            if r["id"] == self.primary_reference_alias
        )
        return ReferenceIdentity.from_dict(reg["identity"])

    def sex_subset_compartments(self) -> dict[str, dict[str, tuple[str, ...]]]:
        """Exact species-applicable compartment sets defining RA-{M,F,P,MF} and CTA-{M,F,P,MF}."""
        a = self.anatomy
        dev = (
            *(c for c in ("embryo", "placenta") if c in a.core_reproductive_compartments),
            a.lactating_mammary_compartment,
        )
        return {
            "RA-M / CTA-M": {
                "expressed_in": a.male_sex_compartments,
                "must_clear_protected": a.female_adult_protected_compartments,
            },
            "RA-F / CTA-F": {
                "expressed_in": (*a.female_sex_compartments, a.resting_mammary_compartment),
                "must_clear_protected": a.male_adult_protected_compartments,
            },
            "RA-P / CTA-P": {
                "expressed_in": dev,
                "must_clear_protected": tuple(
                    dict.fromkeys(
                        (
                            *a.female_adult_protected_compartments,
                            *a.male_adult_protected_compartments,
                        )
                    )
                ),
            },
            "RA-MF / CTA-MF": {
                "expressed_in": a.reproductive_compartments,
                "must_clear_protected": (),
            },
        }

    def dataset_table(self, *, role: str | None = None, status: str | None = None) -> pd.DataFrame:
        if role is not None and role not in ROLE_LABELS:
            raise ValueError(f"Unknown role {role!r}; expected one of {ROLE_IDS}")
        if status is not None and status not in DATASET_STATUSES:
            raise ValueError(f"Unknown status {status!r}; expected one of {DATASET_STATUSES}")
        rows = [
            {
                "taxon_id": self.species.taxon_id,
                "species": self.species.name,
                "role": b.role,
                "role_label": ROLE_LABELS[b.role],
                "source_id": b.source_id,
                "label": b.label,
                "native_unit": b.native_unit,
                "status": b.status,
                "threshold_or_use": b.threshold_or_use,
                "limitations": b.limitations,
            }
            for b in self.datasets
            if (role is None or b.role == role) and (status is None or b.status == status)
        ]
        return pd.DataFrame(
            rows,
            columns=[
                "taxon_id",
                "species",
                "role",
                "role_label",
                "source_id",
                "label",
                "native_unit",
                "status",
                "threshold_or_use",
                "limitations",
            ],
        )

    def to_screening_spec(
        self,
        *,
        comparison_id: str = "",
        comparability_reference: str = "",
        comparison_unit: str | None = None,
        required_normal_contexts: frozenset[str] = frozenset(),
        reference: ReferenceIdentity | None = None,
    ) -> SpeciesScreeningSpec:
        chosen_unit = comparison_unit or self.thresholds.comparison_units[0]
        if chosen_unit not in self.thresholds.comparison_units:
            raise ValueError(
                f"Unit {chosen_unit!r} is not declared in {self.species.name} policy comparison_units "
                f"{self.thresholds.comparison_units}"
            )
        return SpeciesScreeningSpec(
            species=self.species,
            reference=reference or self.primary_reference,
            policy_id=self.policy_id,
            comparison_id=comparison_id,
            comparability_reference=comparability_reference,
            comparison_unit=chosen_unit,
            required_normal_contexts=frozenset(required_normal_contexts),
            anatomy=self.anatomy,
            disposable_cap=self.thresholds.disposable_cap,
            enrichment_fold=self.thresholds.enrichment_fold,
            detection_floor=self.thresholds.detection_floor,
        )

    def as_dict(self) -> dict:
        return {
            "schema": POLICY_SCHEMA_ID,
            "taxon_id": self.species.taxon_id,
            "species": self.species.name,
            "scientific_name": self.species.scientific_name,
            "policy_id": self.policy_id,
            "primary_reference_alias": self.primary_reference_alias,
            "registered_reference_aliases": list(self.registered_reference_aliases),
            "overview": self.overview,
            "anatomy": self.anatomy.as_dict(),
            "sex_subset_compartments": {
                k: {field_k: list(field_v) for field_k, field_v in v.items()}
                for k, v in self.sex_subset_compartments().items()
            },
            "thresholds": asdict(self.thresholds),
            "datasets": [asdict(d) for d in self.datasets],
        }


@cache
def _policies() -> dict[int, SpeciesPolicy]:
    human_sp = resolve_species("human")
    dog_sp = resolve_species("dog")
    mouse_sp = resolve_species("mouse")

    human_policy = SpeciesPolicy(
        species=human_sp,
        policy_id=target_selection.POLICY_ID,
        primary_reference_alias="ensembl115",
        registered_reference_aliases=("ensembl115",),
        anatomy=anatomy_for_species(human_sp),
        thresholds=SpeciesThresholdConfig(
            comparison_units=("TPM",),
            negative_ihc_multipliers={"Approved": 2.0, "Supported": 4.0, "Enhanced": 8.0},
            cornea_cp10k_scale=0.01,
        ),
        overview=(
            "Complete proteome-wide evaluation across 20,099 Ensembl 115 coding loci "
            "(19,893 protein-equivalence groups) with all eight evidence pillars populated."
        ),
        datasets=(
            SpeciesDatasetBinding(
                "reference_proteome",
                "ensembl115",
                "Human Ensembl 115 GRCh38 + MANE Select",
                "amino-acid sequence",
                "active_in_catalog",
                "Prefer MANE Select then Ensembl canonical; collapse identical translations into PEG; >=95% 9-mer Unsafe containment -> Other:Sequence",
                "Representative-protein equivalence does not assert equal abundance of alternative isoforms.",
            ),
            SpeciesDatasetBinding(
                "healthy_baseline_rna",
                "gtex_v11 + hpa_v23_biopsies + sorted_immune",
                "GTEx v11, HPA v23 tissues/brain & HPA/DICE/Monaco immune",
                "TPM",
                "active_in_catalog",
                "Within-site linear donor P95 sets non-thymic somatic N and Disposable D (<=32 TPM); H0 = max(N, D, 1 TPM)",
                "56 required healthy contexts; bulk thymus stays in thymic scope outside N and D.",
            ),
            SpeciesDatasetBinding(
                "normal_bulk_and_immune_safety",
                "hpa_rna_consensus + gtex_tissue_rna + hpa_brain + dice + monaco",
                "HPA v23 consensus/brain/immune, GTEx v11 medians & DICE/Monaco",
                "TPM; nTPM",
                "active_in_catalog",
                "Base ceilings Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) in native unit; matched negative IHC raises ceiling 2x/4x/8x",
                "TPM and nTPM run separate native gates without cross-unit conversion.",
            ),
            SpeciesDatasetBinding(
                "normal_cell_and_rare_umi",
                "cellxgene_census + siletti_brain + hpa_single_cell + collin_cornea",
                "CELLxGENE Census, Siletti brain, HPA v23 scRNA clusters & Collin cornea",
                "UMI counts; nTPM; CP10K",
                "active_in_catalog",
                "UMI: >=10 cells and >=1% >=2-UMI cells across >=2 consortium-deduplicated study groups (or >=5% in 1 study; 2x in Protected); cornea CP10K /100 nominal ceiling",
                "Paralog/3'-host mapping ambiguity without bulk or qPCR confirmation stays Other:Review.",
            ),
            SpeciesDatasetBinding(
                "normal_protein",
                "hpa_normal_tissue + external_ihc",
                "HPA v23 normal immunohistochemistry & curated external staining",
                "categorical IHC grade + reliability",
                "active_in_catalog",
                "Approved/Supported/Enhanced positive staining in Vital/Protected tissues -> Unsafe; routine-context negative IHC grants RNA headroom",
                "Antibody reliability is categorical support, not a calibrated false-negative rate.",
            ),
            SpeciesDatasetBinding(
                "reproductive_and_developmental_rna",
                "gtex_repro + hpa_repro + petropoulos + tyser + mole",
                "GTEx/HPA adult reproductive samples & Petropoulos/Tyser/Molè embryos",
                "TPM; UMI/staining context",
                "active_in_catalog",
                "RA qualifies at >=10 x H0 in adult reproductive samples or cell-supported embryo lineages; per-compartment maxima feed HostContext (<8 TPM / <8 nTPM, 2x with negative IHC)",
                "Early embryo atlases do not cover mid/late fetal stages; pregnancy protection remains conservative.",
            ),
            SpeciesDatasetBinding(
                "tumor_rna",
                "oncoref_tumor_union_135_cohorts",
                "135 QC-passing per-sample in-vivo human tumor cohorts (TCGA, TARGET, BEAT-AML, Treehouse, GEO)",
                "TPM",
                "active_in_catalog",
                "CA qualifies when >=1 QC-passing in-vivo tumor library reaches >=10 x H0; tags ectopic vs lineage_retained activation",
                "Cross-study cohorts use library gDNA QC and nonnegative linear TPM checks rather than batch correction.",
            ),
            SpeciesDatasetBinding(
                "curated_literature_and_clinical",
                "normal_expression_literature + rare_cell_attribution_reviews + cta_literature_149",
                "Reviewed normal-tissue literature, #711 Ensembl 115 attribution ledger & 149 CTA/oncofetal papers",
                "curated verdicts & citations",
                "active_in_catalog",
                "Verified adult Vital/Protected expression -> Unsafe; unscreened adult home (olfactory/chemosensory GPCRs, taste bud, bone, dental pulp, hair follicle) -> Other:Unscreened tissue",
                "Published CTA nominations are provenance only and never override normal-tissue counterevidence.",
            ),
        ),
    )

    dog_policy = SpeciesPolicy(
        species=dog_sp,
        policy_id=POLICY_SCHEMA_ID,
        primary_reference_alias="ensembl98",
        registered_reference_aliases=("ensembl98", "refseq106"),
        anatomy=anatomy_for_species(dog_sp),
        thresholds=SpeciesThresholdConfig(
            comparison_units=("TPM", "TMM-CPM"),
            negative_ihc_multipliers={},
            cornea_cp10k_scale=None,
        ),
        overview=(
            "Complete 8-pillar genome-wide quantitative canine oncology evaluation across the "
            "SHA-256-pinned Ensembl 98 CanFam3.1 host proteome (20,199 coding genes), 24 required "
            "(36 evaluated) adult somatic normal contexts, GSE225599 healthy companion-dog "
            "circulating leukocyte 10x scRNA-seq (7 donors), Pfizer 38-tissue quantitative TMT "
            "proteomics (MSV000099744), comparative canine reproductive/embryonic compartments, "
            "19 spontaneous canine patient tumor cohorts across 13 histologies, and 12 canine "
            "cancer cell lines."
        ),
        datasets=(
            SpeciesDatasetBinding(
                "reference_proteome",
                "ensembl98 + refseq106",
                "Canine Ensembl 98 CanFam3.1 (GCA_000002285.2, primary) & NCBI RefSeq 106 (UU_Cfam_GSD_1.0)",
                "amino-acid sequence",
                "active_in_catalog",
                "Complete GTF/cDNA/ncRNA/peptide FASTA pins; identical translations collapse into canine PEG:<sha256>; >=95% canine 9-mer Unsafe containment",
                "CanFam3.1 is a female Boxer assembly; unplaced scaffolds are retained and Y-linked single-copy loci remain unrepresented.",
            ),
            SpeciesDatasetBinding(
                "healthy_baseline_rna",
                "bgee_v15_2_dog_25exps + gse219045_dog_atlas + companion_normal_5gse",
                "Bgee 15.2 adult canine multi-organ RNA-seq, GSE219045 300-sample organ atlas & healthy companion-dog adrenal/blood/marrow/gingiva/brain/liver/lung/spleen controls",
                "TPM",
                "active_in_catalog",
                "Within-organ donor P95 across 24 required (36 evaluated) non-thymic adult somatic contexts sets N and canine Disposable D (gallbladder, tonsil, hair follicle); <=2.5% ectopic OR gDNA QC",
                "Unvalidated monogenic chemosensory GPCRs, dental-pulp and taste-bud home loci remain Other:Unscreened tissue.",
            ),
            SpeciesDatasetBinding(
                "normal_bulk_and_immune_safety",
                "bgee_v15_2_dog_somatic + gse219045_dog_atlas + gse283373_alveolar_mac",
                "Bgee 15.2 adult canine somatic/leukocyte RNA, GSE219045 whole-blood/marrow/lymph-node atlas & GSE283373 primary canine alveolar macrophages",
                "TPM",
                "active_in_catalog",
                "Tiered somatic & macrophage ceilings Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) TPM without negative-IHC multipliers",
                "Sorted activated canine T/NK cytokine subsets remain narrower than human DICE/Monaco and murine ImmGen.",
            ),
            SpeciesDatasetBinding(
                "normal_cell_and_rare_umi",
                "gse225599_healthy_dog_pbmc_10x",
                "GSE225599 healthy companion-dog circulating leukocyte 10x single-cell UMI atlas (7 healthy donors H_1..H_7)",
                "UMI counts",
                "active_in_catalog",
                "Evaluates donor-resolved circulating T/NK lymphocytes, monocytes, neutrophils (Protected) and platelets (Vital) on native CanFam3.1 ENSCAFG IDs via assess_rare_cell_rna",
                "Covers circulating hematopoietic & platelet populations; solid-organ adult canine normal scRNA atlases remain future additions.",
            ),
            SpeciesDatasetBinding(
                "normal_protein",
                "pfizer-dog-protein",
                "Pfizer 38-tissue canine quantitative proteomics (MSV000099744, 9,929 protein groups)",
                "within-plex bleed-corrected molar ppm (from rIBAQ)",
                "active_in_catalog",
                "Sequence-reconciled Vital (>25 ppm Unsafe, 10-25 ppm hold) & Protected (>50 ppm Unsafe, 20-50 ppm hold) tissue-pool protein gates after 10% same-plex reproductive TMT bleed subtraction; never grants negative-IHC RNA multipliers",
                "All 38 tissue samples are pooled TMT channels without individual-donor counts or cell-type resolution.",
            ),
            SpeciesDatasetBinding(
                "reproductive_and_developmental_rna",
                "bgee_repro_embryo + gse219045_dog_repro + gse119810_normal_mammary + pfizer_fallopian_endometrium",
                "Bgee 15.2 testis/ovary/placenta/d33-d44 embryos, GSE219045 male/female reproductive tract, GSE119810 normal resting/lactating mammary & Pfizer Fallopian_Tube/Endometrium proteomics",
                "TPM; molar ppm (Fallopian_Tube)",
                "active_in_catalog",
                "RA qualifies at >=10 x H0 across canine compartments (uterine_horn, absent seminal vesicle, milk-operon secretory separation) and populates canine HostContext",
                "Canine fallopian tube lacks a public RNA-seq atlas and is screened via Pfizer quantitative proteomics (coverage='proteome_only').",
            ),
            SpeciesDatasetBinding(
                "tumor_rna",
                "dog_19_spontaneous_cohorts_plus_12_cell_lines",
                "19 spontaneous companion-dog patient tumor cohorts across 13 histologies (OSA, CMT, HSA, LYM, GLIOMA, URO, STS, COSCC, CAA, HS, MNG, ACC, THCA, LEUK) + 12 canine cancer cell lines (CLBL-1, 8 MGT, 3 HS)",
                "TPM",
                "active_in_catalog",
                "CA qualifies strictly in >=1 QC-passing (<=2.5% OR gDNA) spontaneous in-vivo patient tumor at >=10 x H0; cell lines populate cell_line_* model-selection sidecar only",
                "Spontaneous oral melanoma linear-TPM patient cohorts remain a priority future histology.",
            ),
            SpeciesDatasetBinding(
                "curated_literature_and_clinical",
                "hindriks2025 + nemec2019 + canine_til_effector_gates",
                "Curated canine CTA nominations (14 candidates), somatic RT-qPCR holds & essential canine activated-TIL checkpoint/cytokine exclusions",
                "curated RT-qPCR & cross-reference index",
                "active_in_catalog",
                "Reconciles all 14 published canine nominations onto CanFam3.1, holds somatic RT-PCR positives, excludes essential TIL checkpoints/cytokines, and enforces zero cross-species CTA transfer",
                "Human CTA labels are never transferred to canine orthologs.",
            ),
        ),
    )

    mouse_policy = SpeciesPolicy(
        species=mouse_sp,
        policy_id=POLICY_SCHEMA_ID,
        primary_reference_alias="ensembl115",
        registered_reference_aliases=("ensembl115",),
        anatomy=anatomy_for_species(mouse_sp),
        thresholds=SpeciesThresholdConfig(
            comparison_units=("TPM",),
            negative_ihc_multipliers={},
            cornea_cp10k_scale=None,
        ),
        overview=(
            "Complete 8-pillar genome-wide murine evaluation across the SHA-256-pinned Ensembl 115 "
            "GRCm39 host proteome (66,668 proteins), 30 adult somatic normal contexts, ImmGen + "
            "IL-15/IL-2/ionomycin NK & Treg sorted immune gates, Tabula Muris 10x single-cell UMI, "
            "MGI GXD 6.27 adult protein staining, all 16 murine reproductive/embryonic compartments, "
            "14 in-vivo syngeneic/GEMM/carcinogen tumor cohorts across 9 histologies, and 9 pure "
            "in-vitro cancer cell lines."
        ),
        datasets=(
            SpeciesDatasetBinding(
                "reference_proteome",
                "ensembl115",
                "Mouse Ensembl 115 GRCm39 (GCA_000001635.9, C57BL/6J, 66,668 proteins)",
                "amino-acid sequence",
                "active_in_catalog",
                "Complete GTF/cDNA/ncRNA/peptide FASTA pins; identical translations collapse into murine PEG:<sha256>; >=95% murine 9-mer Unsafe containment",
                "Reference coordinates are C57BL/6J; DBA/2, BALB/c and FVB polymorphisms require explicit sequence verification.",
            ),
            SpeciesDatasetBinding(
                "healthy_baseline_rna",
                "bgee_v15_2_mouse_12exps + gse219045_mouse_atlas + gse322844_wt_liver",
                "Bgee 15.2 adult murine multi-organ RNA-seq (12 experiments), GSE219045 sex-stratified preclinical organ atlas & GSE322844 WT liver",
                "TPM",
                "active_in_catalog",
                "Within-organ donor P95 across 30 non-thymic adult somatic contexts sets N; H0 = max(N, D, 1 TPM); thymus stays in thymic scope",
                "Unvalidated monogenic chemosensory GPCRs, hair-follicle, dental-pulp and taste-bud home loci remain Other:Unscreened tissue.",
            ),
            SpeciesDatasetBinding(
                "normal_bulk_and_immune_safety",
                "bgee_v15_2_mouse_somatic + gse219045 + immgen_gse109125_gse122108 + gse181223_bmdm + gse287306_gse336796_gse302070_nk_treg + sarma1999",
                "Bgee 15.2 somatic RNA, GSE219045 organ atlas, ImmGen sorted/activated immune, LPS BMDMs, splenic/IL-15/IL-2/ionomycin NK & Treg panels & LAK granzyme/Sarma 1999 gates",
                "TPM; qualitative RT-PCR/Southern",
                "active_in_catalog",
                "Quantitative somatic & sorted/activated immune panels apply Vital <=4, Protected <=8, Manageable <=16 (Unsafe >128) TPM without negative-IHC multipliers; excludes LAK cytolytic granzymes Gzmd/e/f/g and holds Trap1a",
                "Sorted ULI FACS populations require >=2-replicate support (min(P95, second_highest)) to reject single-tube carryover.",
            ),
            SpeciesDatasetBinding(
                "normal_cell_and_rare_umi",
                "tabula_muris_10x",
                "Tabula Muris 10x droplet single-cell UMI atlas across 10 adult Vital/Protected somatic organs",
                "UMI counts",
                "active_in_catalog",
                "Enters assess_rare_cell_rna (>=10 cells, >=5% single-study or bulk/GXD-confirmed >=1% >=2-UMI in Vital, 2x in Protected) with Ensembl 115 31-mer/host mapping-risk flags",
                "Single-study 10x UMI signals between 1% and 5% without bulk or GXD confirmation await a second independent murine UMI atlas.",
            ),
            SpeciesDatasetBinding(
                "normal_protein",
                "mgi_gxd_6_27",
                "MGI Gene Expression Database (GXD 6.27) TS28 adult Immunohistochemistry (MMO:0000498) & Immunoblot (MMO:0000672)",
                "curated adult IHC/Immunoblot detections",
                "active_in_catalog",
                ">=2 independent publications in Vital/Protected adult somatic organs -> Unsafe; 1 publication -> Other:Review hold; never grants negative-IHC RNA ceiling multipliers",
                "Literature-curated assays cover studied proteins unevenly; unassayed organs never imply protein absence.",
            ),
            SpeciesDatasetBinding(
                "reproductive_and_developmental_rna",
                "bgee_repro_embryo + gse219045 + accessory_tract_lactation_7gse",
                "Bgee 15.2 & GSE219045 gonads/placenta/embryo plus WT anterior prostate (coagulating_gland), ventral prostate, seminal vesicle, oviduct, uterine horn, endometrium, cervix, vagina & lactating mammary",
                "TPM",
                "active_in_catalog",
                "RA qualifies at >=10 x H0 across all 16 murine compartments (including coagulating_gland, seminal vesicle, uterine_horn, lactating_mammary) and populates 16-organ HostContext",
                "Embryo libraries remain bulk/tissue developmental stages (scope_complete=False under pregnancy protection).",
            ),
            SpeciesDatasetBinding(
                "tumor_rna",
                "mouse_14_invivo_cohorts_plus_9_cell_lines",
                "14 strain-explicit in-vivo murine tumor cohorts (MC38, CT26, B16F10, B78_MEL, LLC, EMT6, RM1, GL261_TUMOR, CT2A_TUMOR, SB28_TUMOR, MADR_GLIOMA, KP_LUAD, BMPC_PRAD, DEN_HCC) + 9 pure cancer cell lines (4T1, P815, EO771, Hepa1-6, ID8, GL261_LINE, MB49, MOC1, MOC22)",
                "TPM",
                "active_in_catalog",
                "CA qualifies strictly in >=1 QC-passing in-vivo tumor sample at >=10 x H0 on Ensembl 115 TPM; pure cell lines populate cell_line_* purity/model-selection sidecar only",
                "Covers 14 widely used syngeneic, orthotopic, primary GEMM, and carcinogen-induced models across C57BL/6, BALB/c, and FVB.",
            ),
            SpeciesDatasetBinding(
                "curated_literature_and_clinical",
                "vandeneynde1991 + sarma1999 + lak_granzyme_counterevidence",
                "Van den Eynde 1991 (P1A DBA/2 P815 clone 211), Sarma 1999 normal-tissue follow-up & IL-2/IL-12/IL-18 LAK granzyme counterevidence",
                "curated primary-source XMLs & GenBank records",
                "active_in_catalog",
                "Preserves exact 224-aa AAA40475.1 / ENSMUSP00000059448.4 protein equivalence alongside BALB/c somatic counterevidence and host TIL granzyme exclusions",
                "Transgenic tolerance experiments are excluded from endogenous normal-expression claims.",
            ),
        ),
    )

    return {
        human_sp.taxon_id: human_policy,
        dog_sp.taxon_id: dog_policy,
        mouse_sp.taxon_id: mouse_policy,
    }


def species_policy(species: Species | str | int = "human") -> SpeciesPolicy:
    """Return the declarative screening policy and dataset bindings for ``species``."""
    resolved = resolve_species(species)
    return _policies()[resolved.taxon_id]


def all_species_policies() -> tuple[SpeciesPolicy, ...]:
    """Policies for all registered species (human, dog, mouse) in canonical order."""
    return tuple(species_policy(sp) for sp in available_species())


def role_comparison_table() -> pd.DataFrame:
    """Side-by-side comparison of datasets, units and status across human, dog and mouse."""
    policies = {sp.species.name: sp for sp in all_species_policies()}
    rows = []
    for role_id, role_label in EVIDENCE_ROLES:
        entry: dict[str, str] = {"role": role_id, "role_label": role_label}
        for name in ("human", "dog", "mouse"):
            bindings = [b for b in policies[name].datasets if b.role == role_id]
            entry[f"{name}_sources"] = "; ".join(
                f"{b.source_id} ({b.native_unit})" for b in bindings
            )
            entry[f"{name}_status"] = ", ".join(dict.fromkeys(b.status for b in bindings))
            entry[f"{name}_rule"] = "; ".join(b.threshold_or_use for b in bindings)
        rows.append(entry)
    return pd.DataFrame(rows)
