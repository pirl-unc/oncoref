"""Current CTA/oncofetal candidate decisions with explicit species coverage.

``results``, ``coverage`` and ``sources`` use the same contract for human, dog
and mouse. ``audit`` adds source-list accounting and measured decision contexts.
Human rows describe proteins: identical selected reference sequences collapsed. Native
nonhuman nominations remain unassigned/unassessed until their reference identities,
measurement scales and required panels support a species-specific policy.
Human expression, antibody grades and clinical outcomes never transfer through
orthology. The historical nomination/evidence APIs remain separately accessible.

Safety and CTA identity are separate assessments. Only normal-tissue or applicable
clinical safety evidence supplies Unsafe; failing cancer/reproductive admission
does not. Cancer presence requires positive (>0) native TPM/pTPM in one measured tumor sample or reviewed reliable tumor protein.
Reproductive is the umbrella for reproductive, placental, extraembryonic, embryonic
and fetal contexts; its measured branches have positive-expression rules. Required
coverage includes the GTEx v11 52-context bulk core for every member locus.
Neither presence check imposes a normal-expression ratio. Whole-profile
diffuseness is an active adult/placental identity gate. A supported diffuse
profile can establish Other:Diffuse, never Unsafe; developmental concentration remains
explicitly unassessed without compatible measurements.

Reviewed TCR-T evidence can raise a protein-profile RNA ceiling by one ladder
step, capped at Enhanced-negative. The IHC grade is unchanged. This inference
retains exact peptide/HLA/product limitations and does not override protein or
clinical-injury exclusions. See :mod:`oncoref.tcr_trials` for cohorts and events.
"""

from __future__ import annotations

import hashlib
import json
import warnings
from importlib.resources import files
from typing import Literal

import pandas as pd

from .cta_selection_audit import audit as audit
from .species import Species, resolve_species

Category = Literal["CTA", "Unsafe", "Other", "Neither"]
OtherSubtype = Literal[
    "Low", "Disposable", "Diffuse", "Reproductive", "Tumor-supported", "Inducible", "Normal", "Weak"
]
SpeciesLike = Species | str | int
CATEGORIES = ("CTA", "Unsafe", "Other")
OTHER_SUBTYPES = (
    "Low",
    "Disposable",
    "Diffuse",
    "Reproductive",
    "Tumor-supported",
    "Inducible",
    "Normal",
    "Weak",
)
CANDIDATE_CATEGORIES = ("CTA",)
EVIDENCE_CATEGORIES = (
    "nomination",
    "reproductive_expression",
    "developmental_expression",
    "tumor_expression",
    "normal_bulk_rna",
    "normal_cell_rna",
    "immune_rna",
    "normal_protein",
    "clinical_safety",
)
RESULT_COLUMNS = (
    "taxon_id",
    "species",
    "entity_id",
    "entity_kind",
    "label",
    "gene_ids",
    "gene_identities",
    "identity_status",
    "reference_id",
    "category",
    "subtype",
    "classification",
    "assessment_status",
    "policy_id",
    "reasons",
    "evidence_flags",
    "source_snapshot_id",
)
JSON_COLUMNS = ("gene_ids", "gene_identities", "reasons", "evidence_flags")


def _metadata() -> dict:
    """Read the packaged release receipt without network or mutable defaults."""
    metadata = json.loads(files("oncoref").joinpath("data/cta-selection.json").read_text())
    if metadata.get("schema") != "oncoref.cta-selection" or metadata.get("schema_version") != 1:
        raise ValueError("Unsupported CTA selection metadata schema")
    return metadata


def _checked_table(name: str) -> pd.DataFrame:
    path = files("oncoref").joinpath("data", name)
    expected = _metadata()["data_sha256"][name]
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"CTA selection artifact checksum mismatch: {name}")
    with path.open("rb") as handle:
        return pd.read_csv(
            handle, compression="gzip" if name.endswith(".gz") else None, keep_default_na=False
        )


def _native_nominations(species: Species) -> pd.DataFrame:
    """Expose native candidates and unknowns without synthesizing protein groups."""
    from .species_antigens import cta_snapshot

    rows = []
    for candidate in cta_snapshot(species).candidates():
        mapping = candidate["mapping"]
        identities = tuple(
            {
                "gene_id": target["id"],
                "symbol": target["symbol"],
                "namespace": mapping["gene_namespace"],
                "historical_gene_ids": None,
                "alternate_gene_ids": None,
                "synonyms": None,
                "identity_status": mapping["status"],
                "history_status": "not_curated",
            }
            for target in mapping["targets"]
        )
        rows.append(
            {
                "taxon_id": species.taxon_id,
                "species": species.name,
                "entity_id": f"nomination:{species.taxon_id}:{candidate['id']}",
                "entity_kind": "native_nomination",
                "label": candidate["label"],
                "gene_ids": tuple(candidate["gene_ids"]),
                "gene_identities": identities,
                "identity_status": "resolved"
                if mapping["status"] == "unique" and mapping["complete"]
                else mapping["status"],
                "reference_id": mapping["reference_id"],
                "category": None,
                "subtype": None,
                "classification": "Unassessed",
                "assessment_status": "unassessed",
                "policy_id": None,
                "reasons": (
                    "species_selection_policy_unavailable",
                    *candidate["evidence_coverage"]["missing"],
                ),
                "evidence_flags": (f"native_evidence_{candidate['evidence_status']}",),
                "source_snapshot_id": candidate["snapshot_id"],
            }
        )
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def results(
    *,
    species: SpeciesLike = "human",
    category: Category | None = None,
    subtype: OtherSubtype | None = None,
) -> pd.DataFrame:
    """Return detached current decisions in the common species-scoped schema.

    CTA means every gate of the named research policy passes. Unsafe denotes
    exclusion evidence under that policy, not demonstrated toxicity. Other
    is the adequately assessed normal-pass noncandidate family. Every Other row
    has a subtype; classification is the single categorical label (e.g. Other:Low).
    Tumor-supported remains an Other category, not admitted target clearance. A null
    category denotes unassigned evidence; ``assessment_status`` explains why.
    ``entity_kind`` prevents comparing human protein groups with unresolved
    nonhuman assay/family nominations. Tuple-valued identity fields retain every
    supported alternative. No download, orthology inference or zero filling occurs.

    Cancer presence and reproductive/developmental association are separate from
    safety. All must pass for CTA admission. Neither presence check compares with
    a normal-tissue maximum. The complete HPA adult/placental concentration gate
    blocks that route above 10%; supported lineage cells and localized indirect
    cells provide positive association without overriding it. Normal induction
    counterevidence can invalidate indirect-only support. Reviewed normal production separately supplies safety exclusion or a hold. Reviewed reproductive protein supplies presence only. Supported developmental evidence remains an
    independent route with concentration explicitly unassessed.
    """
    selected = resolve_species(species)
    if category == "Neither":
        warnings.warn(
            "category='Neither' is deprecated; use 'Other'. Returned rows use Other.",
            DeprecationWarning,
            stacklevel=2,
        )
        category = "Other"
    if category is not None and category not in CATEGORIES:
        raise ValueError(f"Unknown category {category!r}; choose from {CATEGORIES}")
    if subtype is not None and (subtype not in OTHER_SUBTYPES or category not in {None, "Other"}):
        raise ValueError(f"Unknown or incompatible Other subtype {subtype!r}")
    if selected.taxon_id == 9606:
        frame = _checked_table("cta-selection-human.csv.gz")
        frame["category"] = frame.category.replace("", pd.NA)
        frame["subtype"] = frame.subtype.replace("", pd.NA)
        for column in JSON_COLUMNS:
            frame[column] = frame[column].map(lambda value: tuple(json.loads(value)))
    else:
        frame = _native_nominations(selected)
    if tuple(frame.columns) != RESULT_COLUMNS or not frame.taxon_id.eq(selected.taxon_id).all():
        raise ValueError("CTA result schema or species mismatch")
    if not (frame.category.isna() | frame.category.isin(CATEGORIES)).all():
        raise ValueError("Unknown packaged CTA decision")
    other = frame.category.eq("Other").fillna(False)
    if (
        not frame.loc[other, "subtype"].isin(OTHER_SUBTYPES).all()
        or frame.loc[~other, "subtype"].notna().any()
    ):
        raise ValueError("Other requires a known subtype; other categories cannot have one")
    expected = frame.category.fillna("Unassessed").astype(str)
    expected.loc[other] = "Other:" + frame.loc[other, "subtype"]
    if not frame.classification.eq(expected).all():
        raise ValueError("Categorical classification does not match category and subtype")
    frame["taxon_id"] = frame.taxon_id.astype("int64")
    for column in RESULT_COLUMNS:
        if column not in (*JSON_COLUMNS, "taxon_id"):
            frame[column] = frame[column].astype("string")
    if category is not None:
        frame = frame.loc[frame.category.eq(category)]
    if subtype is not None:
        frame = frame.loc[frame.category.eq("Other") & frame.subtype.eq(subtype)]
    return frame.reset_index(drop=True)


def candidates(*, species: SpeciesLike = "human") -> pd.DataFrame:
    """Return CTA candidates; Other:Tumor-supported remains a reconsideration lead.

    This is a research screening set, not a clinical safety assertion or an
    exact-peptide clearance. Existing ``results(category='CTA')`` keeps its scope.
    """
    frame = results(species=species)
    return frame.loc[frame.category.isin(CANDIDATE_CATEGORIES)].reset_index(drop=True)


def sequence_dispositions(*, species: SpeciesLike = "human") -> pd.DataFrame:
    """Return current screening roles for direct targets and sequence overlap.

    CTA is a direct-target candidate. Unsafe contributes to an exclusion
    background; unresolved rows require review and must not be silently treated
    as allowed. Other contributes no exclusion under this named screening policy.
    These are protein-level research roles, not epitope/HLA clearance. They do
    not modify the independently versioned legacy exact-protein backgrounds or
    cover every alternative isoform. IDs, reasons and the policy receipt remain
    attached so downstream code can choose an explicit review policy.
    """
    frame = results(species=species)
    frame["direct_target_candidate"] = frame.category.eq("CTA").fillna(False)
    frame["shared_sequence"] = (
        frame.category.map({"CTA": "allow", "Other": "allow", "Unsafe": "exclude"})
        .fillna("review")
        .astype("string")
    )
    return frame


def sources(*, species: SpeciesLike = "human", category: str | None = None) -> pd.DataFrame:
    """Common evidence categories with native units and actual integration status.

    A source may support multiple categories. ``used_in_selection`` differs from
    available matrices, metadata-only sources and literature-identified resources.
    Listing an atlas never establishes that its measurements have been ingested.
    """
    selected = resolve_species(species)
    if category is not None and category not in EVIDENCE_CATEGORIES:
        raise ValueError(f"Unknown evidence category {category!r}")
    frame = _checked_table("cta-selection-sources.csv")
    frame = frame.loc[frame.taxon_id.eq(selected.taxon_id)]
    if category is not None:
        frame = frame.loc[frame.evidence_category.eq(category)]
    return frame.reset_index(drop=True)


def coverage(*, species: SpeciesLike = "human") -> dict:
    """Report the selection denominator, source categories and missing capabilities.

    Dog and mouse counts describe curated nominations only, not their genomes.
    Unavailable selection means no quantitative species policy. Source
    counts deduplicate source IDs within each evidence category, not donors.
    """
    selected = resolve_species(species)
    metadata = _metadata()
    policy = metadata["species"][selected.name]
    frame = results(species=selected)
    inventory = sources(species=selected)
    categories = {}
    for name in EVIDENCE_CATEGORIES:
        part = inventory[inventory.evidence_category.eq(name)]
        categories[name] = {
            "used_sources": int(
                part.loc[part.status.eq("used_in_selection"), "source_id"].nunique()
            ),
            "available_sources": int(
                part.loc[
                    part.status.isin(
                        ["used_in_selection", "native_observations", "native_matrix", "annotation"]
                    ),
                    "source_id",
                ].nunique()
            ),
            "identified_sources": int(part.source_id.nunique()),
        }
    return {
        "taxon_id": selected.taxon_id,
        "species": selected.name,
        **policy,
        "entities": len(frame),
        "counts": {
            **{category: int(frame.category.eq(category).sum()) for category in CATEGORIES},
            "Unassigned": int(frame.category.isna().sum()),
        },
        "other_subtypes": frame.loc[frame.category.eq("Other"), "subtype"].value_counts().to_dict(),
        "evidence_categories": categories,
    }
