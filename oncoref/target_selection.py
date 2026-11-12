"""Cancer and Reproductive Antigen sets, with a shared eligibility screen.

``results()`` exposes the rebuilt set-based research catalog. The earlier
``cta_selection.results()`` catalog remains a versioned compatibility snapshot.
Evaluator inputs must follow this contract; legacy >0 tumor presence and 5-unit
reproductive flags are not substitutes for the new evidence checks.
"""

from __future__ import annotations

import json
from copy import deepcopy
from importlib.resources import files

import numpy as np
import pandas as pd

POLICY_ID = "oncoref.target-selection.2026-11-10.1"
DETECTION_FLOOR = 1.0
SINGLE_TUMOR_FOLD = 10.0
REPEATED_TUMORS = 2
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
RARE_CELL_STATUSES = ("established", "dismissed", "below_screen", "unassessed")
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


def definition() -> dict:
    """Return the versioned, machine-readable target-selection contract."""
    return {
        "policy_id": POLICY_ID,
        "status": "research_catalog_rebuilt",
        "species": "human",
        "categories": list(CATEGORIES),
        "shared_checks": list(SHARED_CHECKS),
        "sets": {
            "Cancer Antigens": "eligible AND cancer evidence qualifies",
            "Reproductive Antigens": "eligible AND reproductive/developmental evidence qualifies",
            "CTA": "Cancer Antigens INTERSECTION Reproductive Antigens",
            "Target candidates": "Cancer Antigens UNION Reproductive Antigens",
        },
        "precedence": ["Unsafe", "Other:Review", "Other:Sequence", "set memberships"],
        "other_subtypes": {
            "Review": "Required evidence is missing, unresolved, or fails a shared admission check.",
            "Sequence": "Nested protein without independent support or >=95% Unsafe 9-mer containment.",
            "Expression": (
                "Neither positive expression route qualifies, or the target maximum or "
                "Disposable cap fails."
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
                f"The signal in >= {RARE_CELL_REPLICATE_STUDIES} independent studies of the same "
                f"organ, which agree on a cell type or include >= {RARE_CELL_REPLICATE_DONORS} "
                "donors. Donors of one study share its libraries and reference; one donor in "
                "each of two studies, in different cell types, agrees on nothing but the organ."
            ),
            "single_study": (
                f"Without replication, a signal counts only when >= {RARE_CELL_SINGLE_STUDY_PCT:g}% "
                "of a donor's cells carry multiple molecules."
            ),
            "protected": (
                f"Protected populations need {RARE_CELL_PROTECTED_FACTOR:g} times these frequencies."
            ),
            "single_molecules": (
                "Cells with one molecule never change a category. A mis-assigned read from a "
                "paralog, pseudogene or overlapping gene produces exactly this pattern."
            ),
            "mapping_risk": (
                "Shared transcript sequence with another locus, or an exon inside another "
                "same-strand gene. Counts for a flagged gene cannot be attributed to it."
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
                "A signal that is neither: counts at a flagged locus, or from one study, that no "
                "other assay confirms. Recorded with its reason; no category effect."
            ),
            "literature": (
                "Where a reviewed publication reports the gene in the same organ by a "
                "gene-specific assay, 1% detection in one donor is established, single "
                "molecules included: two assay types then agree. Published protein in a Vital "
                "or Protected tissue excludes through the knowledge layer without any counts."
            ),
            "consequences": {
                "established": "Unsafe",
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
            "independent_participants_required": REPEATED_TUMORS,
            "repeated_rule": "T >= 1 AND T > N in at least two independent participants",
            "single_tumor_fold": SINGLE_TUMOR_FOLD,
            "single_rule": "T >= 10 * max(N, 1) in at least one qualified tumor sample",
            "normal_reference": (
                "Maximum QC-passing, identity-resolved measurement in the declared compatible "
                "normal panel, excluding specifically allowed Disposable and reproductive compartments."
            ),
            "comparability": (
                "One documented comparison scale and complete required normal contexts. Equal unit "
                "names alone do not establish comparability. No TPM/nTPM/pTPM/CP10K pooling."
            ),
            "missing_reference": "unassessed, never zero",
            "alternative": "Reviewed reliable tumor protein with exact identity and source scope.",
        },
        "reproductive_evidence": {
            "branches": ["adult reproductive", "placental", "extraembryonic", "embryonic", "fetal"],
            "rna_detection_floor": DETECTION_FLOOR,
            "bulk": "Mapped human reproductive compartment >=1 native TPM/nTPM.",
            "cell": "Mapped reproductive cells >=1 native TPM/nTPM and >=20 reads.",
            "developmental": (
                "Human stage/lineage mean >=1 TPM, >=20 cells, >=5 cells at >=1 TPM, "
                "and >=10% cells at >=1 TPM. Retain embryos, donors and stage."
            ),
            "protein": "Reliable mapped reproductive IHC or reviewed human reproductive protein.",
            "concentration": "Descriptive measurement; not an admission gate for either set.",
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
    with files("oncoref").joinpath("data/target-selection.csv.gz").open("rb") as handle:
        result = pd.read_csv(handle, compression="gzip", low_memory=False)
    if gene is not None:
        result = result[
            result.protein_group_label.eq(gene)
            | result.member_gene_ids.map(lambda ids: gene in json.loads(ids))
        ]
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
    """
    flags = ("unsafe", *SHARED_CHECKS)
    _require_columns(facts, {*flags, "cancer_evidence", "reproductive_evidence"})
    _require_booleans(facts, flags)
    for field in ("cancer_evidence", "reproductive_evidence"):
        if not facts[field].isin(EVIDENCE_STATUSES).all():
            raise ValueError(f"{field} must be one of {EVIDENCE_STATUSES}")
    d = facts.copy(deep=True)
    eligible = ~d.unsafe & d[list(SHARED_CHECKS)].all(axis=1)
    cancer = eligible & d.cancer_evidence.eq("qualifies")
    reproductive = eligible & d.reproductive_evidence.eq("qualifies")
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
    d.loc[other & ~d.sequence_pass, "subtype"] = "Sequence"
    review = ~d[list(REVIEW_CHECKS)].all(axis=1) | (
        eligible & d[["cancer_evidence", "reproductive_evidence"]].eq("unassessed").any(axis=1)
    )
    d.loc[other & review, "subtype"] = "Review"
    d["classification"] = d.category
    d.loc[other, "classification"] = "Other:" + d.loc[other, "subtype"]
    d["target_policy_id"] = POLICY_ID
    return d


def assess_cancer_rna(
    tumors: pd.DataFrame,
    normals: pd.DataFrame,
    *,
    required_normal_contexts: set[str],
    comparison_id: str,
    comparability_reference: str,
    unit: str,
) -> dict:
    """Evaluate one protein in one reviewed tumor/normal comparison panel.

    Both frames require ``entity_id, source_id, comparison_id, unit, value,
    qc_pass, identity_resolved``. Tumors also require ``sample_id, participant_id``;
    normals require ``context_id, scope``. Scope is ``normal``, ``disposable`` or
    ``reproductive`` and must come from reviewed compartment mappings, never from
    the organ name alone. Required contexts are the nonexempt normal contexts.

    The receipt retains observed maxima, supporting samples, and missing coverage.
    Missing/failed normal observations cannot lower the comparison to zero.
    Participant IDs must be resolved across source/cohort aliases. An unknown
    participant can support the single-sample route, but cannot count as replication.
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
    if unit not in {"TPM", "nTPM", "pTPM"}:
        raise ValueError(
            "This rule requires a declared TPM/nTPM/pTPM comparison, not CP10K or counts"
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
    if not normals.scope.isin(["normal", "disposable", "reproductive"]).all():
        raise ValueError("Unknown normal compartment scope")
    if not normals.context_id.map(lambda v: isinstance(v, str) and bool(v.strip())).all():
        raise ValueError("Every normal observation requires a context ID")
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
    receipt = {
        "policy_id": POLICY_ID,
        "entity_id": next(iter(entity_ids)),
        "comparison_id": comparison_id,
        "comparability_reference": comparability_reference,
        "unit": unit,
        "required_normal_contexts": sorted(required_normal_contexts),
        "missing_normal_contexts": missing,
        "normal_maximum": normal_maximum,
        "normal_maximum_observations": reference[reference.value.eq(normal_maximum)].to_dict(
            "records"
        ),
        "single_sample_threshold": None,
        "independent_participants_above_normal": 0,
        "samples_above_normal": 0,
        "unresolved_participant_samples_above_normal": 0,
        "supporting_tumors": [],
        "routes": [],
        "status": "unassessed",
    }
    if missing or normal_maximum is None:
        receipt["reason"] = "Required normal comparison is incomplete"
        return receipt
    threshold = SINGLE_TUMOR_FOLD * max(normal_maximum, DETECTION_FLOOR)
    receipt["single_sample_threshold"] = threshold
    valid_tumor = t.qc_pass & t.identity_resolved & t.value.notna()
    t = t[valid_tumor].drop_duplicates(["source_id", "sample_id"]).copy()
    if t.empty:
        receipt["reason"] = "No qualified tumor observations"
        return receipt
    t["above_normal"] = t.value.ge(DETECTION_FLOOR) & t.value.gt(normal_maximum)
    t["single_sample_pass"] = t.value.ge(threshold)
    known_participant = t.participant_id.map(lambda v: isinstance(v, str) and bool(v.strip()))
    participants = t.loc[t.above_normal & known_participant, "participant_id"].nunique()
    unresolved = int((t.above_normal & ~known_participant).sum())
    receipt["independent_participants_above_normal"] = int(participants)
    receipt["samples_above_normal"] = int(t.above_normal.sum())
    receipt["unresolved_participant_samples_above_normal"] = unresolved
    receipt["supporting_tumors"] = t[t.above_normal | t.single_sample_pass].to_dict("records")
    if participants >= REPEATED_TUMORS:
        receipt["routes"].append("two_independent_tumors_above_normal")
    if t.single_sample_pass.any():
        receipt["routes"].append("one_tumor_tenfold_above_normal")
    receipt["status"] = "qualifies" if receipt["routes"] else "does_not_qualify"
    receipt["reason"] = (
        "Cancer RNA route qualifies" if receipt["routes"] else "No cancer RNA route qualifies"
    )
    if not receipt["routes"] and unresolved and participants + unresolved >= REPEATED_TUMORS:
        receipt["status"] = "unassessed"
        receipt["reason"] = (
            "Patient identity is unresolved for potentially qualifying repeated tumor samples"
        )
    return deepcopy(receipt)


def assess_rare_cell_rna(
    observations: pd.DataFrame,
    *,
    mapping_risk: bool,
    confirmed_tissues: set[str],
    literature_tissues: frozenset[str] | set[str] = frozenset(),
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
    when bulk RNA confirms it; otherwise it is ``dismissed`` with its reason.
    ``literature_tissues`` lists organs where a reviewed publication reports the
    gene by a gene-specific assay. There, a second assay already vouches for the
    locus, so 1% detection in one donor is enough, single molecules included.
    Every outcome is terminal: there is no review state.
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
    d["donor"] = d.source_id + ":" + d.donor_id
    findings = []
    for (tier, tissue), group in d.groupby(["tier", "tissue"], sort=True):
        multi = group[group.multi_molecule]
        two_studies = multi.source_id.nunique() >= RARE_CELL_REPLICATE_STUDIES
        # Two lone donors in different cell types share only the organ label.
        replicated = two_studies and (
            multi.donor.nunique() >= RARE_CELL_REPLICATE_DONORS
            or multi.groupby("cell_type").source_id.nunique().max() >= RARE_CELL_REPLICATE_STUDIES
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
        elif mapping_risk:
            status, basis = "dismissed", "counts cannot be attributed to this locus; no bulk RNA"
        elif two_studies:
            status, basis = "dismissed", "two donors in different cell types; no bulk RNA"
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
                "studies_multi_molecule": int(multi.source_id.nunique()),
                "donors_multi_molecule": int(multi.donor.nunique()),
                "studies_detected": int(group.loc[group.detected, "source_id"].nunique()),
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
    findings.sort(key=lambda f: (f["status"] != "established", -f["top_multi_molecule_pct"]))
    receipt["assessable_donor_populations"] = len(d)
    receipt["max_detection_pct"] = float(d.detection_pct.max())
    receipt["max_multi_molecule_pct"] = float(d.multi_molecule_pct.max())
    receipt["findings"] = findings
    receipt["status"] = findings[0]["status"] if findings else "below_screen"
    receipt["tier"] = findings[0]["tier"] if findings else ""
    receipt["action"] = "unsafe" if receipt["status"] == "established" else "none"
    return deepcopy(receipt)
