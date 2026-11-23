"""Cancer and Reproductive Antigen sets, with a shared eligibility screen.

``results()`` exposes the rebuilt set-based research catalog. The earlier
``cta_selection.results()`` catalog remains a versioned compatibility snapshot.
Evaluator inputs must follow this contract; legacy >0 tumor presence and 5-unit
reproductive flags are not substitutes for the new evidence checks.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from copy import deepcopy
from functools import lru_cache
from importlib.resources import files

import numpy as np
import pandas as pd

from .load_dataset import _register_derived_cache

POLICY_ID = "oncoref.target-selection.2026-11-20.2"
DETECTION_FLOOR = 1.0
#: Both positive compartments use the same fold over the non-reproductive reference.
ENRICHMENT_FOLD = 10.0
TUMOR_FOLD = ENRICHMENT_FOLD  # Public compatibility name.
REPRODUCTIVE_FOLD = ENRICHMENT_FOLD
#: No histology is excluded wholesale; library-level QC also applies to TGCT.
GERM_CELL_TUMORS = ()
#: Separate nominal caps for native TPM and HPA nTPM; never a unit conversion.
DISPOSABLE_CAP = 32.0
#: Rare normal-cell RNA, per adult donor and named cell population.
RARE_CELL_DETECTION_PCT = 1.0
RARE_CELL_MIN_CELLS = 10
RARE_CELL_MULTI_MOLECULES = 2
RARE_CELL_REPLICATE_STUDIES = 2
#: Two studies replicate a signal when they agree on the cell type, or with this many donors.
RARE_CELL_REPLICATE_DONORS = 3
#: Without replication, only a frequent multi-molecule signal is consequential.
RARE_CELL_SINGLE_STUDY_PCT = 5.0
#: Protected populations need twice the frequency, as their bulk ceiling is twice Vital's.
RARE_CELL_PROTECTED_FACTOR = 2.0
RARE_CELL_STATUSES = ("established", "unresolved", "dismissed", "below_screen", "unassessed")
EVIDENCE_STATUSES = ("qualifies", "does_not_qualify", "unassessed")
CATEGORIES = ("CTA", "Cancer Antigen only", "Reproductive Antigen only", "Other", "Unsafe")
SHARED_CHECKS = (
    "normal_screen_pass",
    "required_coverage_complete",
    "identity_translation_pass",
    "evidence_resolved",
    "expression_bounds_pass",
    "sequence_pass",
)
#: Failing one of these is Other:Review; bounds and sequence have their own subtypes.
REVIEW_CHECKS = SHARED_CHECKS[:4]
#: Normal tissues that no expression panel resolves. Bulk skin and tongue dilute hair
#: follicle and taste buds; nose, tooth and mineralized bone (adult osteoblasts/osteocytes;
#: bone marrow aspirate is hematopoietic) are absent. A tumor value was never compared
#: with them, so it cannot on its own admit a protein whose recorded home is one of them.
UNMEASURED_HOME_TISSUES = (
    "olfactory epithelium",
    "nasal mucosa",
    "taste bud",
    "hair follicle",
    "dental pulp",
    "bone",
)
#: Germline, gonadal, placental and embryonic compartments defining Core CTAs.
CORE_REPRODUCTIVE_COMPARTMENTS = ("embryo", "placenta", "testis", "ovary")
#: Adult accessory reproductive, tract and mammary compartments (sex- or state-conditional).
MALE_ACCESSORY_COMPARTMENTS = ("epididymis", "prostate", "seminal vesicle")
FEMALE_TRACT_COMPARTMENTS = ("fallopian tube", "uterus", "endometrium", "cervix", "vagina")
MAMMARY_COMPARTMENTS = ("lactating_breast", "resting_breast")
EXTENDED_REPRODUCTIVE_COMPARTMENTS = (
    *MALE_ACCESSORY_COMPARTMENTS,
    *FEMALE_TRACT_COMPARTMENTS,
    *MAMMARY_COMPARTMENTS,
)
REPRODUCTIVE_COMPARTMENTS = (
    *CORE_REPRODUCTIVE_COMPARTMENTS,
    *EXTENDED_REPRODUCTIVE_COMPARTMENTS,
)
CTA_TIERS = ("core", "core_with_extended_exposure", "extended_only")


#: Supported within-tissue normal-reference aggregation methods.
NORMAL_TISSUE_STATISTICS = ("p95", "mean_plus_2std")


def summarize_normal_tissue_tpm(
    values: pd.DataFrame | pd.Series,
    *,
    method: str = "p95",
    quantile: float = 0.95,
    n_std: float = 2.0,
) -> pd.Series | float:
    """Aggregate donor or biopsy TPM within a normal tissue.

    ``values`` has observations (donors/biopsies) on the index and entities on
    columns (or is a 1-D Series of observations). Supported ``method`` values:

    - ``"p95"`` (default): linearly interpolated quantile (default 0.95).
    - ``"mean_plus_2std"``: sample mean plus ``n_std`` sample standard
      deviations (``ddof=1``; single-observation tissues contribute ``std=0``).
    """
    if method == "p95":
        return values.quantile(quantile, interpolation="linear")
    if method == "mean_plus_2std":
        mean = values.mean()
        if isinstance(values, pd.DataFrame):
            std = values.std(ddof=1).fillna(0.0).where(values.count().gt(0))
        else:
            count = int(values.notna().sum())
            if count == 0:
                return np.nan
            std = float(values.std(ddof=1)) if count > 1 else 0.0
        return mean + n_std * std
    raise ValueError(
        f"Unknown normal tissue statistic {method!r}; choose from {NORMAL_TISSUE_STATISTICS}"
    )


def definition() -> dict:
    """Return the versioned, machine-readable target-selection contract."""
    return {
        "policy_id": POLICY_ID,
        "status": "research_catalog_rebuilt",
        "species": "human",
        "categories": list(CATEGORIES),
        "shared_checks": list(SHARED_CHECKS),
        "sets": {
            "Cancer Antigens": "eligible AND cancer evidence qualifies AND recorded normal home is covered",
            "Reproductive Antigens": "eligible AND reproductive/developmental evidence qualifies",
            "CTA": "Cancer Antigens INTERSECTION Reproductive Antigens",
            "Candidate antigen proteins": "Cancer Antigens UNION Reproductive Antigens",
        },
        "precedence": [
            "Unsafe",
            "Other:Review",
            "Other:Sequence",
            "Other:Unscreened tissue",
            "set memberships",
        ],
        "other_subtypes": {
            "Review": "Required evidence is missing, unresolved, or fails a shared admission check.",
            "Sequence": (
                "Nested protein without independent support, >=95% Unsafe 9-mer containment, or "
                "a read-through locus whose every coding transcript joins two separate genes."
            ),
            "Expression": (
                "Neither positive expression route qualifies, or the target maximum or "
                "Disposable cap fails."
            ),
            "Unscreened tissue": (
                "The protein's recorded normal home is a tissue no expression panel resolves, "
                "and no reproductive evidence qualifies for RA. Its tumor value was "
                "never compared with that tissue, so cancer evidence alone does not admit it."
            ),
        },
        "home_tissue": {
            "unmeasured_tissues": list(UNMEASURED_HOME_TISSUES),
            "recorded_by": (
                "A family rule, a recorded literature finding, a scoped adult normal home in the publication registry, or the nose and skin reference."
            ),
            "effect": (
                "Insufficient home-tissue coverage blocks cancer admission. "
                "Reproductive enrichment can independently qualify RA."
            ),
        },
        "published_normal_expression": {
            "source": "oncoref/data/normal-expression-literature.json",
            "scope": "Every Target candidate of the previous catalog and every rare-cell dependent.",
            "exclude": (
                "A gene-specific report in normal adult tissue that the tissue-tier registry or its "
                "updates tier Vital or Protected, with a quote retrieved from the source. Unsafe."
            ),
            "recorded": (
                "Confirmed findings that do not meet the rubric, and family rules by protein "
                "name, keep their tissue and tier on the row. The only category effect is "
                "the home-tissue rule."
            ),
        },
        "expression_bounds": {
            "target_maximum": (
                "max(tumor, adult reproductive, developmental) >= 1 and above every measured "
                "nonreproductive normal and Disposable value on the comparison scale."
            ),
            "disposable_cap_tpm": DISPOSABLE_CAP,
            "disposable_cap_ntpm": DISPOSABLE_CAP,
            "cap_scope": "Each source scale separately; a screening convention, not a tolerance.",
            "failure": "Other:Expression",
        },
        "rare_cell_rna": {
            "unit": "One adult donor x one named Vital or Protected cell population, UMI counts",
            "signal": (
                f">= {RARE_CELL_DETECTION_PCT:g}% of a donor's cells and >= {RARE_CELL_MIN_CELLS} "
                f"cells carry >= {RARE_CELL_MULTI_MOLECULES} molecules."
            ),
            "replication": (
                f"The signal in >= {RARE_CELL_REPLICATE_STUDIES} independent study groups of the same "
                f"organ, which agree on a cell type or include >= {RARE_CELL_REPLICATE_DONORS} "
                "distinct biological donors. Cell-type concordance also requires distinct biological donors. "
                "Shared consortium collections (Baylor/Chen eye, Sanes eye, Allen reference brain, Michigan kidney, "
                "Teichmann heart) and donors embedded in integrated atlases (HLCA, Baylor/Chen eye, SEA-AD) are "
                "attributed to their original study and biological donor ID before counting replication. "
                "Donors of one study share its libraries and reference; one donor in each of two studies, "
                "in different cell types, agrees on nothing but the organ."
            ),
            "single_study": (
                f"Without replication, a signal counts only when >= {RARE_CELL_SINGLE_STUDY_PCT:g}% "
                "of a donor's cells carry multiple molecules."
            ),
            "protected": (
                f"Protected populations need {RARE_CELL_PROTECTED_FACTOR:g} times these frequencies."
            ),
            "single_molecules": (
                "Single-molecule detection alone does not establish expression. Ambient RNA "
                "and read assignment are possible confounders, not demonstrated explanations."
            ),
            "mapping_risk": (
                "Shared transcript sequence with another locus, or an exon inside another "
                "same-strand gene. Attribution needs read-level review; a flag does not prove artifact."
            ),
            "confirmation": (
                f"Bulk RNA >= {DETECTION_FLOOR:g} native TPM in the same organ in both the HPA and "
                "GTEx sample sets, a different assay with its own read assignment."
            ),
            "established": (
                "A replicated signal at a locus without mapping risk, or any signal with bulk "
                "confirmation. Unsafe."
            ),
            "dismissed": (
                "A weak or inconsistent signal below replication and frequency requirements. "
                "Recorded with its reason; no category effect."
            ),
            "unresolved": (
                "A replicated signal with unresolved locus attribution, or a frequent single-study "
                "signal without independent confirmation. Other:Review pending adjudication."
            ),
            "literature": (
                "Where a reviewed publication reports the gene in the same organ by a "
                "gene-specific assay, 1% detection in one donor is established, single "
                "molecules included: two assay types then agree. Published protein in a Vital "
                "or Protected tissue excludes through the knowledge layer without any counts."
            ),
            "consequences": {
                "established": "Unsafe",
                "unresolved": "Other:Review; strong evidence needs attribution or confirmation.",
                "dismissed": "No finding; the reason and counts stay in the record.",
                "below_screen": "No finding. Not a demonstration of absence or safety.",
                "unassessed": "Recorded as unmeasured; never treated as negative.",
            },
            "reported_sensitivity_cutoffs_pct": [1, 2, 5],
            "interpretation": (
                "Detection frequency measures captured RNA, not the fraction of cells producing "
                "protein. No cutoff is a calibrated toxicity threshold."
            ),
        },
        "cancer_rna": {
            "detection_floor": DETECTION_FLOOR,
            "tumor_fold": TUMOR_FOLD,
            "rule": (
                "T >= 10 * max(N, D, 1) in at least one qualified tumor sample of any comparable "
                "TPM cohort. N is the highest non-thymic, non-reproductive normal value, D the highest "
                "non-thymic Disposable value. This is the only cancer route."
            ),
            "normal_reference": (
                "N: maximum QC-passing, identity-resolved measurement in the declared compatible "
                "peripheral panel, excluding thymus, Disposable and reproductive compartments. D is taken "
                "from the non-thymic Disposable tissues. Reproductive and developmental observations only supply the positive route."
            ),
            "descriptive": (
                "Counts at the actual tenfold threshold distinguish samples, known patients, "
                "anonymous samples and cancer types. Above-normal counts are legacy descriptors."
            ),
            "excluded_tumors": {
                "codes": list(GERM_CELL_TUMORS),
                "reason": (
                    "No histology is excluded wholesale. Reproductive programs are advisory "
                    "sample notes, not evidence of normal-cell contamination, and never exclude "
                    "TGCT or other libraries. Independent QC exclusions still apply."
                ),
            },
            "comparability": (
                "One documented comparison scale and complete required normal contexts. Equal unit "
                "names alone do not establish comparability. No TPM/nTPM/pTPM/CP10K pooling."
            ),
            "missing_reference": "unassessed, never zero",
            "protein": "Independent annotation; qualitative staining cannot satisfy a numerical RNA fold.",
        },
        "reproductive_evidence": {
            "branches": ["adult reproductive", "placental", "extraembryonic", "embryonic", "fetal"],
            "rna_detection_floor": DETECTION_FLOOR,
            "reproductive_fold": REPRODUCTIVE_FOLD,
            "rule": "R >= 10 * max(N, D, 1), using the identical normal panel and comparison scale as CA.",
            "bulk": "At least one QC-passing, identity-resolved reproductive tissue sample on the declared comparison scale.",
            "cell": "Cell-type nTPM and UMI counts remain separate annotations and normal checks; no unit pooling.",
            "developmental": (
                "Mean native TPM within one embryo and lineage. One well-sampled population: "
                ">=20 cells, >=5 cells and >=10% cells at >=1 TPM. Smaller populations: "
                ">=5 measured cells, >=50% at >=1 TPM, and >=2 independent embryos at the same "
                "age and lineage reaching the fold. The mean must reach the tenfold threshold; "
                "higher fold filters retain replication. "
                "Embryos are biological donors; cells are not independent witnesses."
            ),
            "protein": "Reliable reproductive IHC and reviewed protein findings annotate association; they cannot bypass the RNA fold.",
            "concentration": "The same tenfold enrichment requirement as CA, assessed independently.",
            "compartments": {
                "core": list(CORE_REPRODUCTIVE_COMPARTMENTS),
                "male_accessory": list(MALE_ACCESSORY_COMPARTMENTS),
                "female_tract": list(FEMALE_TRACT_COMPARTMENTS),
                "mammary": list(MAMMARY_COMPARTMENTS),
                "lactating_breast_rule": (
                    "Individual breast biopsy samples with coordinated milk-protein RNA "
                    "(>= 3 of LALBA, CSN1S1, CSN2, CSN3 at >= 100 native TPM) are assigned to "
                    "lactating_breast (inferred secretory signature). Complete panels without that signature "
                    "retain the legacy resting_breast label; incomplete marker coverage stays unknown. "
                    "These are not provider-confirmed physiological states."
                ),
                "cta_tier_interpretation": "Legacy origin/exposure annotation retained for compatibility, not a toxicity tier or host eligibility rule. Use reproductive_qualifying_compartments for origins and HostContext for host choices.",
                "cta_tiers": {
                    "core": (
                        "Qualifies in embryo, placenta, testis or ovary, with all extended "
                        "adult compartments <= 32 TPM and <= tumor_max_tpm / 10."
                    ),
                    "core_with_extended_exposure": (
                        "Qualifies in embryo, placenta, testis or ovary, and also exceeds 32 TPM "
                        "or lacks a 10x tumor margin in an extended adult compartment."
                    ),
                    "extended_only": (
                        "Qualifies only via epididymis, prostate, seminal vesicle, fallopian tube, "
                        "uterus, endometrium, cervix, vagina, lactating_breast or resting_breast."
                    ),
                },
            },
        },
        "shared_normal_policy": {
            "policy_id": "oncoref.cta-selection.2026-11-08.1",
            "scope": "Normal assessment only; legacy cancer, reproductive and concentration gates excluded.",
            "allowances": (
                "Same resolved reproductive and Disposable compartment allowances for both sets. "
                "Immune, endothelial and other protected cells do not inherit an organ's allowance."
            ),
            "threshold_interpretation": (
                "RNA/IHC cutoffs are screening conventions, not measured injury thresholds. "
                "Negative IHC does not establish a calibrated safety margin."
            ),
        },
        "translation": (
            "Qualifying coding reference or supported noncanonical ORF plus complete approved "
            "MANE-first/canonical representative; resolve conflicting coding evidence. "
            "Experimental protein detection is retained separately."
        ),
        "published_target_support": (
            "Positive evidence annotation for a curated, exact-scope target nomination. Preserve "
            "every protein-paper link and its attested claim; biology or background-table coverage "
            "alone is not a target nomination. Publication cannot override an exclusion or hold."
        ),
        "only_meaning": "Membership in one qualified set, not proven absence from the other context.",
    }


def _require_columns(frame: pd.DataFrame, columns: set[str]) -> None:
    if missing := columns - set(frame):
        raise ValueError(f"Missing fields: {sorted(missing)}")


def metadata() -> dict:
    """Counts, comparison scope and the hash of the packaged rebuilt catalog."""
    return json.loads(files("oncoref").joinpath("data/target-selection.json").read_text())


@lru_cache(maxsize=1)
def _cached_results_bundle() -> tuple[pd.DataFrame, dict[str, list[int]]]:
    from . import clinical_organ_clearance

    with files("oncoref").joinpath("data/target-selection.csv.gz").open("rb") as handle:
        base = pd.read_csv(handle, compression="gzip", low_memory=False)
    base = clinical_organ_clearance.annotate_catalog(base, enabled=True)
    gene_index: dict[str, list[int]] = {}
    for idx, (lbl, raw_ids) in enumerate(zip(base["protein_group_label"], base["member_gene_ids"])):
        gene_index.setdefault(str(lbl), []).append(idx)
        for gid in json.loads(raw_ids):
            if gid != lbl:
                gene_index.setdefault(str(gid), []).append(idx)
    return base, gene_index


_register_derived_cache(_cached_results_bundle.cache_clear)


def results(
    gene: str | None = None, *, category: str | None = None, changed_only: bool = False
) -> pd.DataFrame:
    """Rebuilt human protein decisions, exact transcripts, values and paper links.

    ``gene`` is an exact protein label or a member Ensembl gene ID. Category uses
    the exclusive display labels in ``CATEGORIES``. ``is_cancer_antigen`` and
    ``is_reproductive_antigen`` include their intersection, CTA. Missing values
    remain missing; the released-category columns retain the migration baseline.
    """
    if category is not None and category not in CATEGORIES:
        raise ValueError(f"Unknown category {category!r}; expected one of {CATEGORIES}")
    base, gene_index = _cached_results_bundle()
    if gene is not None:
        idxs = gene_index.get(gene, [])
        result = base.iloc[idxs].copy()
    else:
        result = base.copy()
    if category is not None:
        result = result[result.category.eq(category)]
    if changed_only:
        result = result[result.category.ne(result.released_category)]
    return result.reset_index(drop=True)


def _require_booleans(frame: pd.DataFrame, columns) -> None:
    for column in columns:
        if not frame[column].map(lambda v: isinstance(v, (bool, np.bool_))).all():
            raise ValueError(f"{column} requires explicit booleans, without missing values")


def classify(facts: pd.DataFrame) -> pd.DataFrame:
    """Assign sets to assessed proteins, retaining every input evidence column.

    Required columns are ``unsafe``, the six ``SHARED_CHECKS``, and
    ``cancer_evidence``/``reproductive_evidence`` (``EVIDENCE_STATUSES``).
    Protein IDs, source links, findings and publication annotations pass through.
    An unassessed positive association does not erase the other set's membership;
    missing required *normal* coverage blocks both sets.

    Optional ``home_tissue_unmeasured`` marks a protein whose recorded normal home no
    panel resolves. Cancer admission requires adequate home-tissue coverage.
    Reproductive enrichment can still establish RA.
    """
    flags = ("unsafe", *SHARED_CHECKS)
    _require_columns(facts, {*flags, "cancer_evidence", "reproductive_evidence"})
    _require_booleans(facts, flags)
    for field in ("cancer_evidence", "reproductive_evidence"):
        if not facts[field].isin(EVIDENCE_STATUSES).all():
            raise ValueError(f"{field} must be one of {EVIDENCE_STATUSES}")
    d = facts.copy(deep=True)
    eligible = ~d.unsafe & d[list(SHARED_CHECKS)].all(axis=1)
    reproductive = eligible & d.reproductive_evidence.eq("qualifies")
    if "home_tissue_unmeasured" in d:
        _require_booleans(d, ["home_tissue_unmeasured"])
        home = d.home_tissue_unmeasured.astype(bool)
    else:
        home = pd.Series(False, index=d.index)
    uncompared = eligible & d.cancer_evidence.eq("qualifies") & home
    cancer = eligible & d.cancer_evidence.eq("qualifies") & ~uncompared
    d["eligible"] = eligible
    d["is_cancer_antigen"] = cancer
    d["is_reproductive_antigen"] = reproductive
    d["is_cta"] = cancer & reproductive
    d["is_target_candidate"] = cancer | reproductive
    d["category"] = "Other"
    d.loc[cancer, "category"] = "Cancer Antigen only"
    d.loc[reproductive, "category"] = "Reproductive Antigen only"
    d.loc[cancer & reproductive, "category"] = "CTA"
    d.loc[d.unsafe, "category"] = "Unsafe"
    d["subtype"] = ""
    other = d.category.eq("Other")
    d.loc[other, "subtype"] = "Expression"
    d.loc[other & uncompared, "subtype"] = "Unscreened tissue"
    d.loc[other & ~d.sequence_pass, "subtype"] = "Sequence"
    review = ~d[list(REVIEW_CHECKS)].all(axis=1) | (
        eligible & d[["cancer_evidence", "reproductive_evidence"]].eq("unassessed").any(axis=1)
    )
    d.loc[other & review, "subtype"] = "Review"
    d["classification"] = d.category
    d.loc[other, "classification"] = "Other:" + d.loc[other, "subtype"]
    d["target_policy_id"] = POLICY_ID
    return d


def _assess_enrichment_rna(
    tumors: pd.DataFrame,
    normals: pd.DataFrame,
    *,
    required_normal_contexts: set[str],
    comparison_id: str,
    comparability_reference: str,
    unit: str,
    compartment: str,
) -> dict:
    """Evaluate one protein in one reviewed tumor/normal comparison panel.

    Both frames require ``entity_id, source_id, comparison_id, unit, value,
    qc_pass, identity_resolved``. Tumors also require ``sample_id, participant_id``;
    normals require ``context_id, scope``. Scope is ``normal``, ``thymic``, ``disposable`` or
    ``reproductive`` and must come from reviewed compartment mappings, never from
    the organ name alone. Required contexts are the nonexempt normal contexts.

    Cancer evidence qualifies when one tumor sample reaches ``TUMOR_FOLD`` times the
    highest peripheral non-reproductive value: max(normal, Disposable, ``DETECTION_FLOOR``).
    The receipt retains observed maxima, supporting samples, and missing coverage.
    Missing/failed normal observations cannot lower the comparison to zero. Counts of
    independent patients above the normal maximum are descriptive.
    """
    common = {
        "entity_id",
        "source_id",
        "comparison_id",
        "unit",
        "value",
        "qc_pass",
        "identity_resolved",
    }
    _require_columns(tumors, common | {"sample_id", "participant_id"})
    _require_columns(normals, common | {"context_id", "scope"})
    if unit not in {"TPM", "nTPM", "pTPM", "TMM-CPM"}:
        raise ValueError(
            "This rule requires a declared TPM/nTPM/pTPM/TMM-CPM comparison, not CP10K or counts"
        )
    if not comparison_id.strip() or not comparability_reference.strip():
        raise ValueError("A comparison ID and documented comparability reference are required")
    if not required_normal_contexts or any(
        not isinstance(c, str) or not c.strip() for c in required_normal_contexts
    ):
        raise ValueError("Declare the nonempty required nonexempt normal context panel")
    entity_ids = set()
    for frame in (tumors, normals):
        _require_booleans(frame, ["qc_pass", "identity_resolved"])
        if not frame.unit.eq(unit).all() or not frame.comparison_id.eq(comparison_id).all():
            raise ValueError("Do not combine different units or comparison panels")
        for field in ("entity_id", "source_id"):
            if not frame[field].map(lambda v: isinstance(v, str) and bool(v.strip())).all():
                raise ValueError(f"Every observation requires {field}")
        entity_ids.update(frame.entity_id)
        values = pd.to_numeric(frame.value, errors="raise")
        if ((values < 0) | (~np.isfinite(values) & values.notna())).any():
            raise ValueError("Expression must be finite and nonnegative, or missing")
    if len(entity_ids) != 1:
        raise ValueError("Assess one resolved protein identity at a time")
    if not normals.scope.isin(["normal", "disposable", "reproductive", "thymic"]).all():
        raise ValueError("Unknown normal compartment scope")
    if not normals.context_id.map(lambda v: isinstance(v, str) and bool(v.strip())).all():
        raise ValueError("Every normal observation requires a context ID")
    if normals.context_id.eq("thymus").ne(normals.scope.eq("thymic")).any():
        raise ValueError("Thymus observations require a separate thymic scope")
    if "thymus" in required_normal_contexts:
        raise ValueError("Thymus cannot be a required peripheral comparison context")
    if normals.groupby("context_id").scope.nunique().gt(1).any():
        raise ValueError("A normal context has conflicting compartment allowances")
    if not tumors.sample_id.map(lambda v: isinstance(v, str) and bool(v.strip())).all():
        raise ValueError("Every tumor observation requires a sample ID")
    for field in ("participant_id", "value"):
        if tumors.groupby(["source_id", "sample_id"])[field].nunique(dropna=False).gt(1).any():
            raise ValueError(f"Conflicting {field} for the same tumor observation")
    n = normals.copy()
    t = tumors.copy()
    n["value"] = pd.to_numeric(n.value)
    t["value"] = pd.to_numeric(t.value)
    valid_normal = n.scope.eq("normal") & n.qc_pass & n.identity_resolved & n.value.notna()
    reference = n[valid_normal]
    missing = sorted(required_normal_contexts - set(reference.context_id))
    normal_maximum = None if reference.empty else float(reference.value.max())
    disposable = n[n.scope.eq("disposable") & n.qc_pass & n.identity_resolved & n.value.notna()]
    disposable_maximum = float(disposable.value.max()) if not disposable.empty else 0.0
    receipt = {
        "policy_id": POLICY_ID,
        "entity_id": next(iter(entity_ids)),
        "comparison_id": comparison_id,
        "comparability_reference": comparability_reference,
        "unit": unit,
        "positive_compartment": compartment,
        "required_normal_contexts": sorted(required_normal_contexts),
        "missing_normal_contexts": missing,
        "normal_maximum": normal_maximum,
        "normal_maximum_observations": reference[reference.value.eq(normal_maximum)].to_dict(
            "records"
        ),
        "disposable_maximum": disposable_maximum,
        "single_sample_threshold": None,
        "independent_participants_above_normal": 0,
        "samples_above_normal": 0,
        "unresolved_participant_samples_above_normal": 0,
        "samples_at_tenfold": 0,
        "independent_participants_at_tenfold": 0,
        "anonymous_samples_at_tenfold": 0,
        "participant_or_sample_witnesses_at_tenfold": 0,
        "independent_participant_folds": [],
        "supported_maximum": None,
        "supported_maximum_observations": [],
        "supporting_tumors": [],
        "routes": [],
        "status": "unassessed",
    }
    if missing or normal_maximum is None:
        receipt["reason"] = "Required normal comparison is incomplete"
        return receipt
    threshold = ENRICHMENT_FOLD * max(normal_maximum, disposable_maximum, DETECTION_FLOOR)
    receipt["single_sample_threshold"] = threshold
    valid_tumor = t.qc_pass & t.identity_resolved & t.value.notna()
    t = t[valid_tumor].drop_duplicates(["source_id", "sample_id"]).copy()
    if t.empty:
        receipt["reason"] = f"No qualified {compartment} observations"
        return receipt
    t["above_normal"] = t.value.ge(DETECTION_FLOOR) & t.value.gt(normal_maximum)
    # A replicated developmental observation is limited by the required donor's
    # expression, not the strongest cell pool from one embryo.
    support = t.support_ceiling if compartment == "reproductive" else t.value
    t["single_sample_pass"] = t.value.ge(threshold) & support.ge(threshold)
    receipt["supported_maximum"] = float(support.max()) if support.notna().any() else None
    receipt["supported_maximum_observations"] = t[support.eq(receipt["supported_maximum"])].to_dict(
        "records"
    )
    known_participant = t.participant_id.map(lambda v: isinstance(v, str) and bool(v.strip()))
    participants = t.loc[t.above_normal & known_participant, "participant_id"].nunique()
    unresolved = int((t.above_normal & ~known_participant).sum())
    receipt["independent_participants_above_normal"] = int(participants)
    receipt["samples_above_normal"] = int(t.above_normal.sum())
    receipt["unresolved_participant_samples_above_normal"] = unresolved
    qualifying = t.single_sample_pass
    receipt["independent_participant_folds"] = (
        (
            t.loc[qualifying & known_participant].groupby("participant_id").value.max()
            / (threshold / ENRICHMENT_FOLD)
        )
        .sort_values(ascending=False)
        .tolist()
    )
    receipt["samples_at_tenfold"] = int(qualifying.sum())
    receipt["independent_participants_at_tenfold"] = int(
        t.loc[qualifying & known_participant, "participant_id"].nunique()
    )
    receipt["anonymous_samples_at_tenfold"] = int((qualifying & ~known_participant).sum())
    receipt["participant_or_sample_witnesses_at_tenfold"] = (
        receipt["independent_participants_at_tenfold"] + receipt["anonymous_samples_at_tenfold"]
    )
    receipt["supporting_tumors"] = t[t.above_normal | t.single_sample_pass].to_dict("records")
    if t.single_sample_pass.any():
        receipt["routes"].append(f"{compartment}_tenfold_above_nonreproductive_normal")
    receipt["status"] = "qualifies" if receipt["routes"] else "does_not_qualify"
    receipt["reason"] = (
        f"A {compartment} sample is at least 10 x every non-reproductive normal value"
        if receipt["routes"]
        else f"No {compartment} sample reaches 10 x the non-reproductive normal maximum"
    )
    return deepcopy(receipt)


def assess_cancer_rna(tumors: pd.DataFrame, normals: pd.DataFrame, **comparison) -> dict:
    """Assess cancer samples using the shared tenfold RNA comparison contract.

    See :func:`assess_reproductive_rna` for the symmetric positive compartment.
    Tumors require source, sample and participant IDs (empty participant ID means
    unknown), value, unit, comparison ID, QC and identity flags. Normal observations
    require context IDs and reviewed normal/disposable/reproductive scopes.
    """
    return _assess_enrichment_rna(tumors, normals, compartment="tumor", **comparison)


def developmental_donor_requirement(cells: int, positive_cells: int) -> int:
    """Donors required for a native-TPM cell mean; zero means inadequate cell support."""
    if (
        any(
            isinstance(v, (bool, np.bool_)) or not np.isfinite(v) or v < 0 or int(v) != v
            for v in (cells, positive_cells)
        )
        or positive_cells > cells
    ):
        raise ValueError("Require integer counts with 0 <= positive_cells <= cells")
    if cells >= 20 and positive_cells >= 5 and positive_cells / cells >= 0.1:
        return 1
    if cells >= 5 and positive_cells / cells >= 0.5:
        return 2
    return 0


def assess_reproductive_rna(
    reproductive: pd.DataFrame, normals: pd.DataFrame, **comparison
) -> dict:
    """Assess reproductive/developmental samples against the same normal reference.

    In addition to the cancer observation schema, every positive row requires a
    reviewed ``scope`` of reproductive or developmental. nTPM, UMI counts and
    qualitative protein findings cannot be mixed into a native-TPM comparison.
    Receipt counts distinguish observed samples, known donors and anonymous samples.
    """
    _require_columns(reproductive, {"scope"})
    if not reproductive.scope.isin(["reproductive", "developmental"]).all():
        raise ValueError("RA witnesses must be reproductive or developmental observations")
    reproductive = reproductive.reset_index(drop=True).copy()
    _require_booleans(reproductive, ["qc_pass", "identity_resolved"])
    reproductive["support_ceiling"] = reproductive.value.where(
        reproductive.qc_pass & reproductive.identity_resolved
    )
    developmental = reproductive.scope.eq("developmental")
    if developmental.any():
        _require_columns(
            reproductive, {"source_context", "developmental_stage", "cells", "positive_cells"}
        )
        for field in ("source_context", "developmental_stage"):
            if (
                not reproductive.loc[developmental, field]
                .map(lambda v: isinstance(v, str) and bool(v.strip()))
                .all()
            ):
                raise ValueError(
                    "Developmental observations require a known stage and lineage context"
                )
        required = pd.Series(0, index=reproductive.index)
        required.loc[developmental] = [
            developmental_donor_requirement(row.cells, row.positive_cells)
            for row in reproductive[developmental].itertuples()
        ]
        supported = reproductive.participant_id.map(
            lambda v: isinstance(v, str) and bool(v.strip())
        )
        supported &= required.gt(0)
        reproductive.loc[developmental & ~supported, "support_ceiling"] = np.nan
        for _, group in reproductive[developmental].groupby(
            ["source_id", "developmental_stage", "source_context"]
        ):
            known = group.participant_id.map(lambda v: isinstance(v, str) and bool(v.strip()))
            donors = group[known].groupby("participant_id").support_ceiling.max().dropna()
            ceiling = donors.nlargest(2).iloc[-1] if len(donors) >= 2 else np.nan
            small = group[required.loc[group.index].eq(2)]
            reproductive.loc[small.index, "support_ceiling"] = (
                small.support_ceiling.clip(upper=ceiling) if pd.notna(ceiling) else np.nan
            )
    receipt = _assess_enrichment_rna(
        reproductive, normals, compartment="reproductive", **comparison
    )
    receipt["supporting_reproductive_samples"] = receipt.pop("supporting_tumors")
    return receipt


#: Verified consortium groupings and integrated-atlas donor provenance rules for
#: CELLxGENE Census rare-cell replication (#715). Separate collections from the
#: same consortium or integrated atlases embedding primary-study donors cannot
#: fabricate independent multi-study replication.
_BAYLOR_EYE_COLLECTIONS = frozenset(
    {
        "05e3d0fc-c9dd-4f14-9163-2b242b3bb5c2",
        "0cbf8ef8-87bd-4b51-8521-ec3f68183e11",
        "0f7d022a-46c7-4e64-be4c-e34adbb78089",
        "48c15b0c-6039-4e0f-9668-b6b3c0b830ad",
        "4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e",
    }
)
_SANES_EYE_COLLECTIONS = frozenset(
    {
        "63d03351-06be-478e-a0db-f7a653b6b19b",
        "7855daa4-3f34-4aae-b65a-87720c02e7cd",
    }
)
_ALLEN_REFERENCE_COLLECTIONS = frozenset(
    {
        "283d65eb-dd53-496d-adb7-7570c7caa443",
        "d17249d2-0e6e-4500-abb8-e6c93fa1ac6f",
    }
)
_ALLEN_LINKED_DONOR_COLLECTIONS = _ALLEN_REFERENCE_COLLECTIONS | {
    "1ca90a2d-2943-483d-b678-b809bf464c30"
}
_MICHIGAN_KIDNEY_COLLECTIONS = frozenset(
    {
        "4cbb929b-b03b-4aa8-a943-00f61dc22641",
        "a98b828a-622a-483a-80e0-15703678befd",
    }
)
_TEICHMANN_HEART_COLLECTIONS = frozenset(
    {
        "3116d060-0a8e-4767-99bb-e866badea1ed",
        "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    }
)
_HLCA_COLLECTION = "6f6d381a-7701-4781-935c-db10d30de293"
_ILD_LUNG_COLLECTION = "07e12576-b41b-4350-adfc-059bfc4328ea"
_MADISSOON_LUNG_COLLECTION = "4d74781b-8186-4c9a-b659-ff4dc4601d91"
_MADISSOON_DONORS = frozenset(
    {"296C", "325C", "328C", "343B", "356C", "362C", "367C", "368C", "390C"}
)
_SANES_DONOR_RE = re.compile(r"^(?:Hu\d{1,4}|Pt\d+)$")
_ALLEN_DONOR_RE = re.compile(r"^(?:H\d{2}\.\d{2}\.\d{3}|H200\.\d{4})$")
_ILD_DONOR_RE = re.compile(r"^(?:THD|VUHD)\d+")


def resolve_rare_cell_provenance(source_id: str, donor_id: str) -> tuple[str, str]:
    """Map a Census (collection_id, donor_id) pair to (study_group, linked_donor).

    Collections from the same consortium share one ``study_group``. Donors embedded
    in integrated atlases (such as Sanes eye donors in Chen's retinal atlas or ILD
    and Madissoon donors in HLCA) are attributed back to their original study group
    and biological donor identity so they cannot count as independent replication.
    """
    stripped_sanes = re.sub(r"^sanes(?:2022)?_", "", donor_id)
    if source_id in _SANES_EYE_COLLECTIONS or (
        source_id in _BAYLOR_EYE_COLLECTIONS and _SANES_DONOR_RE.fullmatch(stripped_sanes)
    ):
        return "sanes_eye_atlas", f"sanes_eye_atlas:{stripped_sanes}"
    if source_id in _BAYLOR_EYE_COLLECTIONS:
        norm = (
            "Atlanta_" + donor_id[: -len("_Atlanta")] if donor_id.endswith("_Atlanta") else donor_id
        )
        return "baylor_eye_atlas", f"baylor_eye_atlas:{norm}"
    if source_id in _ALLEN_LINKED_DONOR_COLLECTIONS:
        group = "allen_reference_brain" if source_id in _ALLEN_REFERENCE_COLLECTIONS else source_id
        if _ALLEN_DONOR_RE.fullmatch(donor_id):
            return group, f"Allen:{donor_id}"
        return group, f"{group}:{donor_id}"
    if source_id == _HLCA_COLLECTION and "483747" in donor_id:
        suffix = donor_id.split("483747", 1)[1]
        if _ILD_DONOR_RE.match(suffix):
            return _ILD_LUNG_COLLECTION, f"{_ILD_LUNG_COLLECTION}:{suffix}"
        if suffix in _MADISSOON_DONORS:
            return _MADISSOON_LUNG_COLLECTION, f"{_MADISSOON_LUNG_COLLECTION}:{suffix}"
    if source_id in _MICHIGAN_KIDNEY_COLLECTIONS:
        return "michigan_kidney", f"michigan_kidney:{donor_id}"
    if source_id in _TEICHMANN_HEART_COLLECTIONS:
        return "teichmann_heart", f"teichmann_heart:{donor_id}"
    return source_id, f"{source_id}:{donor_id}"


def assess_rare_cell_rna(
    observations: pd.DataFrame,
    *,
    mapping_risk: bool,
    confirmed_tissues: set[str],
    literature_tissues: frozenset[str] | set[str] = frozenset(),
    provenance_resolver: Callable[[str, str], tuple[str, str]] | None = None,
) -> dict:
    """Grade sparse normal-cell RNA for one gene from donor-resolved molecule counts.

    Each row is one adult donor in one named cell population of primary normal
    tissue: ``gene_id, source_id, donor_id, tissue, cell_type, tier, n_cells,
    n_positive, n_multi``. ``source_id`` identifies an independent study and
    ``tissue`` the organ. ``n_multi`` counts cells with at least
    ``RARE_CELL_MULTI_MOLECULES`` deduplicated molecules; supply UMI counts only.

    A signal needs multi-molecule cells in at least 1% of a donor population.
    Single-molecule detection is
    reported and never graded: one mis-assigned read per cell is what a paralog
    or overlapping gene produces. ``mapping_risk`` says the reference makes such
    mis-assignment possible for this gene; ``confirmed_tissues`` lists organs with
    orthogonal bulk-RNA support. A signal is a multi-molecule frequency seen in
    two independent studies of the organ, or in at least 5% of one donor's cells.
    Two studies replicate it when both find it in the same cell type, or when at
    least ``RARE_CELL_REPLICATE_DONORS`` donors carry it.
    It is ``established`` when replicated at a locus without mapping risk, or
    when bulk RNA confirms it. Strong unresolved signals require review; weak or
    inconsistent signals are dismissed with their reason.
    ``literature_tissues`` lists organs where a reviewed publication reports the
    gene by a gene-specific assay. There, a second assay already vouches for the
    locus, so 1% detection in one donor is enough, single molecules included.
    A reference mapping flag is not evidence that the observed reads are artifacts.
    """
    columns = {
        "gene_id",
        "source_id",
        "donor_id",
        "tissue",
        "cell_type",
        "tier",
        "n_cells",
        "n_positive",
        "n_multi",
    }
    _require_columns(observations, columns)
    if not isinstance(mapping_risk, (bool, np.bool_)):
        raise ValueError("mapping_risk requires an explicit boolean")
    d = observations.copy()
    if d.gene_id.nunique() > 1:
        raise ValueError("Assess one gene at a time; never add frequencies across paralogs")
    if not d.tier.isin(["vital", "protected"]).all():
        raise ValueError("Rare-cell tiers must be vital or protected")
    for field in ("source_id", "donor_id", "tissue", "cell_type"):
        if not d[field].map(lambda v: isinstance(v, str) and bool(v.strip())).all():
            raise ValueError(f"Every observation requires {field}")
    if d.duplicated(["source_id", "donor_id", "tissue", "cell_type"]).any():
        raise ValueError("One row per donor and cell population")
    counts = d[["n_cells", "n_positive", "n_multi"]].apply(pd.to_numeric, errors="raise")
    if (
        counts.isna().any().any()
        or (counts < 0).any().any()
        or (counts.n_multi > counts.n_positive).any()
        or (counts.n_positive > counts.n_cells).any()
    ):
        raise ValueError("Counts must satisfy 0 <= n_multi <= n_positive <= n_cells")
    d[["n_cells", "n_positive", "n_multi"]] = counts
    receipt = {
        "policy_id": POLICY_ID,
        "gene_id": d.gene_id.iloc[0] if len(d) else None,
        "status": "unassessed",
        "tier": "",
        "action": "none",
        "mapping_risk": bool(mapping_risk),
        "assessable_donor_populations": 0,
        "max_detection_pct": None,
        "max_multi_molecule_pct": None,
        "findings": [],
    }
    # A donor population too small to reach ten cells at 1% cannot be assessed.
    floor = RARE_CELL_MIN_CELLS * 100 / RARE_CELL_DETECTION_PCT
    d = d[d.n_cells.ge(floor)].copy()
    if d.empty:
        return receipt
    d["detection_pct"] = 100 * d.n_positive / d.n_cells
    d["multi_molecule_pct"] = 100 * d.n_multi / d.n_cells
    d["detected"] = d.n_positive.ge(RARE_CELL_MIN_CELLS) & d.detection_pct.ge(
        RARE_CELL_DETECTION_PCT
    )
    factor = d.tier.map({"vital": 1.0, "protected": RARE_CELL_PROTECTED_FACTOR})
    enough = d.n_multi.ge(RARE_CELL_MIN_CELLS)
    d["multi_molecule"] = enough & d.multi_molecule_pct.ge(RARE_CELL_DETECTION_PCT * factor)
    d["frequent"] = enough & d.multi_molecule_pct.ge(RARE_CELL_SINGLE_STUDY_PCT * factor)
    resolver = resolve_rare_cell_provenance if provenance_resolver is None else provenance_resolver
    resolved = [resolver(s, did) for s, did in zip(d.source_id, d.donor_id)]
    d["study_group"] = [r[0] for r in resolved]
    d["donor"] = [r[1] for r in resolved]
    findings = []
    for (tier, tissue), group in d.groupby(["tier", "tissue"], sort=True):
        multi = group[group.multi_molecule]
        two_studies = multi.study_group.nunique() >= RARE_CELL_REPLICATE_STUDIES
        # Two lone donors in different cell types share only the organ label.
        concordance = multi.groupby("cell_type")[["study_group", "donor"]].nunique()
        same_cell = concordance.study_group.ge(RARE_CELL_REPLICATE_STUDIES) & concordance.donor.ge(
            2
        )
        replicated = two_studies and (
            multi.donor.nunique() >= RARE_CELL_REPLICATE_DONORS or same_cell.any()
        )
        published = tissue in literature_tissues and group.detected.any()
        if not two_studies and not group.frequent.any() and not published:
            continue
        confirmed = tissue in confirmed_tissues
        if published:
            status, basis = "established", "detection confirmed by a published gene-specific assay"
        elif replicated and not mapping_risk:
            status, basis = "established", "replicated signal at a locus without mapping risk"
        elif confirmed:
            status, basis = "established", "signal confirmed by bulk RNA"
        elif replicated or group.frequent.any():
            status = "unresolved"
            basis = (
                (
                    "replicated normal-cell signal with unresolved locus attribution"
                    if replicated
                    else "frequent single-study normal-cell signal with unresolved locus attribution"
                )
                if mapping_risk
                else "frequent single-study normal-cell signal needs independent confirmation"
            )
        elif mapping_risk:
            status, basis = (
                "dismissed",
                "insufficient concordant signal; locus attribution uncertain",
            )
        elif two_studies:
            status, basis = (
                "dismissed",
                "insufficient independent donor and cell-type concordance; no bulk RNA",
            )
        else:
            status, basis = "dismissed", "one study only; no bulk RNA"
        ranked = multi if len(multi) else group[group.detected]
        top = ranked.sort_values(["multi_molecule_pct", "detection_pct"], ascending=False).iloc[0]
        findings.append(
            {
                "status": status,
                "basis": basis,
                "tier": tier,
                "tissue": tissue,
                "replicated": replicated,
                "orthogonal_confirmation": confirmed,
                "published_confirmation": published,
                "studies_multi_molecule": int(multi.study_group.nunique()),
                "collections_multi_molecule": int(multi.source_id.nunique()),
                "donors_multi_molecule": int(multi.donor.nunique()),
                "studies_detected": int(group.loc[group.detected, "study_group"].nunique()),
                "donors_detected": int(group.loc[group.detected, "donor"].nunique()),
                "donors_assessed": int(group.donor.nunique()),
                "cell_types_multi_molecule": int(multi.cell_type.nunique()),
                "top_cell_type": top.cell_type,
                "top_donor": top.donor,
                "top_n_cells": int(top.n_cells),
                "top_n_positive": int(top.n_positive),
                "top_n_multi": int(top.n_multi),
                "top_detection_pct": float(top.detection_pct),
                "top_multi_molecule_pct": float(top.multi_molecule_pct),
            }
        )
    rank = {"established": 0, "unresolved": 1, "dismissed": 2}
    findings.sort(
        key=lambda f: (rank[f["status"]], not f["replicated"], -f["top_multi_molecule_pct"])
    )
    receipt["assessable_donor_populations"] = len(d)
    receipt["max_detection_pct"] = float(d.detection_pct.max())
    receipt["max_multi_molecule_pct"] = float(d.multi_molecule_pct.max())
    receipt["findings"] = findings
    receipt["status"] = findings[0]["status"] if findings else "below_screen"
    receipt["tier"] = findings[0]["tier"] if findings else ""
    receipt["action"] = {"established": "unsafe", "unresolved": "review"}.get(
        receipt["status"], "none"
    )
    return deepcopy(receipt)
