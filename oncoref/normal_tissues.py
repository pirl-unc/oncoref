# Licensed under the Apache License, Version 2.0.
"""Cited, species-aware normal-tissue scopes and source vocabulary mappings.

Scopes describe permitted anatomical/developmental contexts, not clinical safety.
Restriction is a subset test over *all* detections, after taking tissue unions.
Neither missing observations nor a source's bulk tissue label proves cell-type
specificity. Frozen CTA facts remain unchanged; these are versioned annotations.
"""

from __future__ import annotations

import csv
import hashlib
import math
from copy import deepcopy
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import yaml

from .species import resolve_species

_DATA = Path(__file__).parent / "data"
SCREENED_NORMAL_TISSUE_TIERS = ("vital", "protected", "manageable")


def _key(value):
    return " ".join(str(value or "").strip().casefold().split())


@cache
def _registry():
    taxonomy = (_DATA / "normal-tissue-taxonomy.yaml").read_bytes()
    mappings = (_DATA / "normal-tissue-mappings.csv").read_bytes()
    doc = yaml.safe_load(taxonomy)
    if doc.get("schema") != "oncoref.normal-tissue-taxonomy" or doc.get("schema_version") != 1:
        raise ValueError("Unsupported normal-tissue taxonomy schema")
    rows = list(csv.DictReader(mappings.decode().splitlines()))
    concepts, scopes = doc["concepts"], doc["scopes"]
    aliases = {}
    for cid, concept in concepts.items():
        if not set(concept["parents"]) <= concepts.keys():
            raise ValueError(f"Unknown anatomy parent: {cid}")
        if not set(concept["references"]) <= doc["references"].keys():
            raise ValueError(f"Unknown anatomy citation: {cid}")
        for name in [cid, cid.replace("_", " "), concept["label"], *concept["aliases"]]:
            key = _key(name)
            if key in aliases and aliases[key] != cid:
                raise ValueError(f"Ambiguous anatomy alias: {name}")
            aliases[key] = cid

    def applicable_taxa(cid, path=()):
        if cid in path:
            raise ValueError(f"Cyclic normal-tissue anatomy: {cid}")
        taxa = set(concepts[cid].get("taxon_ids", doc["taxon_ids"]))
        if not taxa or not taxa <= set(doc["taxon_ids"]):
            raise ValueError(f"Invalid anatomy species: {cid}")
        for parent in concepts[cid]["parents"]:
            if not taxa <= applicable_taxa(parent, (*path, cid)):
                raise ValueError(f"Anatomy species exceed parent applicability: {cid}")
        return taxa

    for cid in concepts:
        applicable_taxa(cid)
    scope_aliases = {}
    for name, scope in scopes.items():
        if not set(scope["roots"]) <= concepts.keys() or not set(scope["union"]) <= scopes.keys():
            raise ValueError(f"Unknown scope member: {name}")
        for alias in [name, *scope["aliases"]]:
            if alias in scope_aliases:
                raise ValueError(f"Duplicate scope alias: {alias}")
            scope_aliases[alias] = name
    for name, policy in doc.get("screening_policies", {}).items():
        if policy["permitted_scope"] not in scopes:
            raise ValueError(f"Unknown screening policy scope: {name}")
        if policy["required_positive_scope"] not in scopes:
            raise ValueError(f"Unknown required positive scope: {name}")
        floor = policy["minimum_reproductive_rna"]
        if not isinstance(floor, (int, float)) or not math.isfinite(floor) or floor < 0:
            raise ValueError(f"Invalid reproductive expression floor: {name}")
        if policy["reproductive_rna_comparison"] != "greater_than":
            raise ValueError(f"Unknown reproductive RNA comparison: {name}")
        roots = set(policy["unresolved_roots"])
        for tier in SCREENED_NORMAL_TISSUE_TIERS:
            roots.update(policy[tier]["roots"])
            limit = policy[tier]["maximum_rna"]
            if not isinstance(limit, (int, float)) or not math.isfinite(limit) or limit <= 0:
                raise ValueError(f"Invalid baseline RNA ceiling: {name}/{tier}")
            selectivity = policy[tier]["minimum_reproductive_selectivity"]
            if not math.isfinite(selectivity) or not 0 < selectivity < 1:
                raise ValueError(f"Invalid baseline selectivity: {name}/{tier}")
        ladder = policy.get("negative_ihc_rna_limits")
        if ladder is not None:
            evidence = ("Missing", "Uncertain", "Approved", "Supported", "Enhanced")
            if set(ladder) != set(evidence):
                raise ValueError(f"Incomplete IHC evidence ladder: {name}")
            previous = [policy[tier]["maximum_rna"] for tier in SCREENED_NORMAL_TISSUE_TIERS]
            for reliability in evidence:
                row = ladder[reliability]
                if set(row) != set(SCREENED_NORMAL_TISSUE_TIERS):
                    raise ValueError(f"Incomplete tissue-category grid: {name}/{reliability}")
                values = [row[tier] for tier in SCREENED_NORMAL_TISSUE_TIERS]
                if any(not math.isfinite(v) or v <= 0 for v in values):
                    raise ValueError(f"Invalid adaptive RNA ceiling: {name}/{reliability}")
                if values != sorted(values) or any(v < p for v, p in zip(values, previous)):
                    raise ValueError(f"Non-monotonic adaptive RNA grid: {name}/{reliability}")
                if reliability in {"Missing", "Uncertain"} and values != previous:
                    raise ValueError(f"Missing/uncertain IHC must use baseline: {name}")
                previous = values
        fractions = policy.get("negative_ihc_minimum_reproductive_selectivity")
        if fractions is not None:
            if set(fractions) != set(ladder):
                raise ValueError(f"Incomplete selectivity evidence ladder: {name}")
            previous = [
                policy[tier]["minimum_reproductive_selectivity"]
                for tier in SCREENED_NORMAL_TISSUE_TIERS
            ]
            for reliability in ("Missing", "Uncertain", "Approved", "Supported", "Enhanced"):
                row = fractions[reliability]
                if set(row) != set(SCREENED_NORMAL_TISSUE_TIERS):
                    raise ValueError(f"Incomplete selectivity category grid: {name}/{reliability}")
                values = [row[tier] for tier in SCREENED_NORMAL_TISSUE_TIERS]
                if any(not math.isfinite(v) or not 0 < v < 1 for v in values):
                    raise ValueError(f"Invalid selectivity threshold: {name}/{reliability}")
                if values != sorted(values, reverse=True) or any(
                    v > p for v, p in zip(values, previous)
                ):
                    raise ValueError(f"Non-monotonic selectivity grid: {name}/{reliability}")
                if reliability in {"Missing", "Uncertain"} and values != previous:
                    raise ValueError(f"Missing/uncertain IHC must use baseline selectivity: {name}")
                previous = values
        if (
            not roots <= concepts.keys()
            or not set(policy["references"]) <= doc["references"].keys()
        ):
            raise ValueError(f"Unknown screening policy anatomy or citation: {name}")
        for review in policy.get("context_reviews", {}).values():
            if (
                not set(review["roots"]) <= concepts.keys()
                or not set(review["references"]) <= doc["references"].keys()
                or not review["flags"]
                or not review["rationale"]
            ):
                raise ValueError(f"Invalid screening context review: {name}")
    index = {}
    for row in rows:
        row["taxon_id"] = int(row["taxon_id"])
        if row["concept_id"] and row["concept_id"] not in concepts:
            raise ValueError("Unknown source anatomy mapping")
        if row["concept_id"] and row["taxon_id"] not in concepts[row["concept_id"]].get(
            "taxon_ids", doc["taxon_ids"]
        ):
            raise ValueError("Source anatomy mapping is not applicable to its species")
        if row["mapping_kind"] not in {
            "exact",
            "reviewed_equivalent",
            "source_subcompartment",
            "organ_context",
            "unresolved",
        } or not row["source_reference"].startswith("https://"):
            raise ValueError("Source mappings require a kind and citation")
        key = (
            row["taxon_id"],
            row["source_id"],
            row["source_version"],
            _key(row["source_tissue"]),
            _key(row["source_cell_type"]),
        )
        if key in index:
            raise ValueError(f"Duplicate source tissue mapping: {key}")
        index[key] = row
    doc["taxonomy_id"] = (
        "oncoref:normal-tissues:v1:sha256:"
        + hashlib.sha256(taxonomy + b"\0" + mappings).hexdigest()
    )
    return doc, aliases, scope_aliases, index


def normal_tissue_taxonomy() -> dict:
    """Detached definitions, scope unions, species applicability and bibliography."""
    return deepcopy(_registry()[0])


def normal_tissue_screening_policy(policy="tiered_cta") -> dict:
    """Detached screening parameters; RNA units and species are explicit."""
    try:
        return deepcopy(_registry()[0]["screening_policies"][policy])
    except KeyError:
        raise ValueError(f"Unknown normal-tissue screening policy {policy!r}") from None


def normal_tissue_mappings(*, species, source_id=None, source_version=None) -> tuple[dict, ...]:
    """Reviewed native labels; an absent mapping never implies absent expression."""
    taxon = resolve_species(species).taxon_id
    return tuple(
        deepcopy(row)
        for row in _registry()[3].values()
        if row["taxon_id"] == taxon
        and (source_id is None or row["source_id"] == source_id)
        and (source_version is None or row["source_version"] == source_version)
    )


def _ancestors(cid, path=()):
    if cid in path:
        raise ValueError("Cyclic normal-tissue anatomy")
    result = {cid}
    for parent in _registry()[0]["concepts"][cid]["parents"]:
        result.update(_ancestors(parent, (*path, cid)))
    return result


def _scope_roots(scope, path=()):
    if scope in path:
        raise ValueError("Cyclic normal-tissue scope union")
    record = _registry()[0]["scopes"][scope]
    roots = set(record["roots"])
    for child in record["union"]:
        roots.update(_scope_roots(child, (*path, scope)))
    return roots


def _scope_name(scope):
    try:
        return _registry()[2][scope]
    except KeyError:
        raise ValueError(f"Unknown normal-tissue scope {scope!r}") from None


@cache
def _members(scope, taxon_id):
    roots = _scope_roots(scope)
    doc = _registry()[0]
    return frozenset(
        cid
        for cid, concept in doc["concepts"].items()
        if taxon_id in concept.get("taxon_ids", doc["taxon_ids"]) and roots & _ancestors(cid)
    )


def normal_tissue_scope(scope, *, species) -> frozenset[str]:
    """Canonical anatomy members after union and species filtering.

    Embryonic/fetal scopes additionally admit explicitly staged organ samples.
    Use resolve_normal_tissue() to evaluate stage-aware source contexts.
    """
    return _members(_scope_name(scope), resolve_species(species).taxon_id)


def normal_tissue_scopes(*, species) -> dict[str, frozenset[str]]:
    """All canonical scopes for a species, including their nested unions."""
    return {name: normal_tissue_scope(name, species=species) for name in _registry()[0]["scopes"]}


@dataclass(frozen=True)
class NormalTissueAnnotation:
    taxon_id: int
    taxonomy_id: str
    policy_version: str
    source_id: str
    source_version: str | None
    source_tissue: str | None
    source_cell_type: str | None
    source_stage: str | None
    stage: str
    physiological_state: str | None
    concept_id: str | None
    mapping_status: str
    mapping_kind: str
    scopes: tuple[str, ...]
    references: tuple[str, ...]
    tolerability_evidence: str
    source_sex: str | None = None
    sex: str = "unspecified"

    def as_dict(self):
        return asdict(self)


def resolve_normal_tissue(
    tissue=None,
    *,
    species,
    cell_type=None,
    stage=None,
    source_id="anatomy",
    source_version=None,
    context="normal",
    sex=None,
    physiological_state=None,
) -> NormalTissueAnnotation:
    """Resolve an exact native label without guessing a narrower source compartment.

    ``anatomy`` accepts only the documented canonical names and aliases. Other
    sources require an exact species/source/version mapping. No human-source
    fallback is used for animals. Preserve all original source labels in output.
    A normal scope is never assigned to tumor/mixed/unknown-context observations.
    """
    selected = resolve_species(species)
    doc, aliases, _, mappings = _registry()
    status, kind, cid, state = "unmapped", "unresolved", None, None
    references = {doc["references"]["project_policy"]["url"]}
    native_stage = _key(stage)
    sex_name = {
        "": "unspecified",
        "unknown": "unspecified",
        "unspecified": "unspecified",
        "male": "male",
        "m": "male",
        "female": "female",
        "f": "female",
        "mixed": "mixed",
    }.get(_key(sex), "unresolved")
    stage_name = {
        "": "unspecified",
        "unknown": "unspecified",
        "life cycle": "unspecified",
        "embryonic": "embryonic",
        "embryo": "embryonic",
        "embryonic stage": "embryonic",
        "fetal": "fetal",
        "foetal": "fetal",
        "fetus": "fetal",
        "fetal stage": "fetal",
        "adult": "postnatal",
        "prime adult stage": "postnatal",
        "postnatal": "postnatal",
        "juvenile": "postnatal",
        "neonatal": "postnatal",
    }.get(native_stage, "unresolved")
    known_source = source_id == "anatomy" or any(
        key[:3] == (selected.taxon_id, source_id, source_version) for key in mappings
    )
    if source_id == "anatomy":
        cid = aliases.get(_key(tissue))
        cell_cid = aliases.get(_key(cell_type))
        if _key(tissue) in {"skin", "skin 1", "skin 2", "zone of skin"} and _key(cell_type) in {
            "melanocyte",
            "melanocytes",
        }:
            cell_cid = "cutaneous_melanocyte"
        if cell_cid and ((cid is None and tissue is None) or cid in _ancestors(cell_cid)):
            cid = cell_cid
        elif cell_cid:
            cid, status = None, "conflicting_anatomy"
        kind = "canonical" if cid else "unresolved"
        if _key(tissue) == "lactating breast":
            state = "lactating"
    else:
        key = (selected.taxon_id, source_id, source_version, _key(tissue), _key(cell_type))
        row = mappings.get(key)
        if row is None and cell_type:
            # Broader organ mapping is explicit; cell specificity is never inferred.
            row = mappings.get((*key[:4], ""))
            kind = "organ_context"
        if row:
            cid = row["concept_id"] or None
            if kind != "organ_context":
                kind = row["mapping_kind"]
            references.add(row["source_reference"])
            state = row["physiological_state"] or None
        if not known_source:
            status = "source_unavailable"
    if cid:
        concept = doc["concepts"][cid]
        references.update(doc["references"][r]["url"] for r in concept["references"])
        if selected.taxon_id not in concept.get("taxon_ids", doc["taxon_ids"]):
            status = "not_applicable_to_species"
        else:
            status = "mapped"
    if cid in {"embryo", "fetus"}:
        inferred_stage = "embryonic" if cid == "embryo" else "fetal"
        if stage_name in {"unspecified", inferred_stage}:
            stage_name = inferred_stage
        else:
            status = "conflicting_stage"
    if stage_name == "unresolved":
        status = "unresolved_stage"
    if sex_name == "unresolved":
        status = "unresolved_sex"
    if physiological_state is not None:
        supplied_state = _key(physiological_state)
        if state and supplied_state != state:
            status = "conflicting_physiological_state"
        else:
            state = supplied_state or None
    memberships = []
    if context != "normal":
        status = "not_normal_context"
    elif status == "mapped" or (
        known_source and status == "unmapped" and stage_name in {"embryonic", "fetal"}
    ):
        for name in doc["scopes"]:
            members = _members(name, selected.taxon_id)
            developmental = {"embryonic": "embryo", "fetal": "fetus"}.get(stage_name)
            if cid in members or developmental in members:
                memberships.append(name)
        if status != "mapped":
            status = "stage_only"
    if "disposable_normal" in memberships:
        references.update(
            doc["references"][r]["url"]
            for r in ("blagosklonny2003", "johnson2009", "nci_alopecia", "purba2019")
        )
    if "permitted_normal" in memberships:
        references.add(doc["references"]["project_policy"]["url"])
    return NormalTissueAnnotation(
        selected.taxon_id,
        doc["taxonomy_id"],
        doc["policy_version"],
        source_id,
        source_version,
        tissue,
        cell_type,
        stage,
        stage_name,
        state,
        cid,
        status,
        kind,
        tuple(memberships),
        tuple(sorted(references)),
        "human_context_dependent_rationale"
        if selected.taxon_id == 9606
        else "anatomical_transfer_only_not_validated_tolerability",
        sex,
        sex_name,
    )


def normal_tissue_screening_tier(annotation, policy="tiered_cta") -> str:
    """Return permitted, vital, protected, manageable, or unresolved for current anatomy.

    This assigns anatomy only. Human HPA nTPM cutoffs must not be transferred to
    another species or assay. Explicit developmental permission takes precedence.
    """
    if not isinstance(annotation, NormalTissueAnnotation):
        raise TypeError("Expected a NormalTissueAnnotation")
    if annotation.taxonomy_id != _registry()[0]["taxonomy_id"]:
        raise ValueError("Resolve anatomy against the current taxonomy before assigning tiers")
    config = normal_tissue_screening_policy(policy)
    if annotation.mapping_status not in {"mapped", "stage_only"}:
        return "unresolved"
    if config["permitted_scope"] in annotation.scopes:
        return "permitted"
    if not annotation.concept_id:
        return "unresolved"
    ancestors = _ancestors(annotation.concept_id)
    tiers = [
        tier for tier in SCREENED_NORMAL_TISSUE_TIERS if ancestors & set(config[tier]["roots"])
    ]
    if len(tiers) > 1:
        raise ValueError("Overlapping screened tissue-category definitions")
    return tiers[0] if tiers else "unresolved"


def normal_tissue_screening_details(annotation, policy="tiered_cta") -> dict:
    """Tier, readable label and context-specific review conditions.

    ``permitted`` is the stable identifier for conditionally acceptable contexts;
    ``vital`` includes critical sensory function. Review flags are annotations,
    not additional numerical tiers or evidence that a therapy is acceptable.
    Developmental permission does not inherit adult organ-loss review conditions.
    """
    tier = normal_tissue_screening_tier(annotation, policy)
    config = normal_tissue_screening_policy(policy)
    ancestors = (
        _ancestors(annotation.concept_id)
        if annotation.mapping_status == "mapped" and annotation.stage not in {"embryonic", "fetal"}
        else set()
    )
    reviews = {
        name: review
        for name, review in config.get("context_reviews", {}).items()
        if ancestors & set(review["roots"])
    }
    doc = _registry()[0]
    category = normal_tissue_category(annotation, policy)
    return {
        "tissue_category": category,
        "tissue_category_label": {
            "reproductive": "Reproductive",
            "disposable": "Disposable",
        }.get(category, config["tier_labels"][tier]),
        "screening_policy_id": config["policy_id"],
        "screening_tier": tier,
        "screening_tier_label": config["tier_labels"][tier],
        "review_flags": sorted({flag for review in reviews.values() for flag in review["flags"]}),
        "context_reviews": deepcopy(reviews),
        "references": sorted(
            set(annotation.references)
            | {
                doc["references"][reference]["url"]
                for review in reviews.values()
                for reference in review["references"]
            }
        ),
        "tolerability_evidence": annotation.tolerability_evidence,
        "policy_references": [doc["references"][key]["url"] for key in config["references"]],
    }


def normal_tissue_category(annotation, policy="tiered_cta") -> str:
    """Category ID; unresolved is a QC status outside the five tissue categories.

    Readable names are returned by normal_tissue_screening_details; each name
    is the capitalized category ID.
    Disposable names the conditional project allowance. It does not
    assert complete organ loss is harmless or inflammation is steroid-responsive.
    The older, narrower ``disposable_normal`` anatomical scope remains distinct.
    """
    tier = normal_tissue_screening_tier(annotation, policy)
    if tier == "permitted":
        return "reproductive" if "reproductive" in annotation.scopes else "disposable"
    return tier


def normal_tissue_inventory(*, species, policy="tiered_cta") -> tuple[dict, ...]:
    """Every canonical concept, including explicitly inapplicable species anatomy.

    This is an anatomy inventory, not an expression-coverage assertion. Dataset
    mapping availability is reported separately by normal_tissue_mapping_inventory.
    """
    doc = _registry()[0]
    result = []
    for cid, concept in doc["concepts"].items():
        annotation = resolve_normal_tissue(cid, species=species)
        details = normal_tissue_screening_details(annotation, policy)
        result.append(
            {
                "concept_id": cid,
                "label": concept["label"],
                "parents": tuple(concept["parents"]),
                "applicable_taxon_ids": tuple(concept.get("taxon_ids", doc["taxon_ids"])),
                "applicability_notes": concept.get("notes", ""),
                "taxon_id": annotation.taxon_id,
                "mapping_status": annotation.mapping_status,
                "tissue_category": details["tissue_category"],
                "tissue_category_label": details["tissue_category_label"],
                "category_policy_id": details["screening_policy_id"],
                "category_policy_references": details["policy_references"],
                "review_flags": details["review_flags"],
                "scopes": annotation.scopes,
                "references": annotation.references,
                "expression_coverage": "not_assessed",
            }
        )
    return tuple(result)


def normal_tissue_mapping_inventory(*, species, source_id=None, source_version=None):
    """Native source/version labels with canonical anatomy and five-way category.

    An empty result means no mapping inventory is bundled, never no expression.
    Exact rows retain compartment versus organ-context mapping resolution.
    """
    result = []
    for row in normal_tissue_mappings(
        species=species, source_id=source_id, source_version=source_version
    ):
        annotation = resolve_normal_tissue(
            row["source_tissue"] or None,
            species=species,
            cell_type=row["source_cell_type"] or None,
            source_id=row["source_id"],
            source_version=row["source_version"],
        )
        details = normal_tissue_screening_details(annotation)
        result.append(
            {
                **row,
                "mapping_status": annotation.mapping_status,
                "tissue_category": details["tissue_category"],
                "tissue_category_label": details["tissue_category_label"],
                "category_policy_id": details["screening_policy_id"],
                "category_policy_references": details["policy_references"],
                "review_flags": details["review_flags"],
                "taxonomy_id": annotation.taxonomy_id,
                "expression_coverage": "not_assessed",
            }
        )
    return tuple(result)


def classify_normal_tissue_pattern(detections) -> dict:
    """Find scopes containing EVERY supplied detection, including mixed branches.

    Callers decide detection using each assay's native semantics. This reports
    compatibility of positive observations, not proof of absence outside a scope.
    Unknown mappings block a claim; empty detections never qualify vacuously.
    """
    detections = tuple(detections)
    if any(not isinstance(row, NormalTissueAnnotation) for row in detections):
        raise TypeError("Expected NormalTissueAnnotation detections")
    if len({row.taxon_id for row in detections}) > 1:
        raise ValueError("Cannot combine detections from different species")
    if len({row.taxonomy_id for row in detections}) > 1:
        raise ValueError("Cannot combine different tissue taxonomy versions")
    if any(row.taxonomy_id != _registry()[0]["taxonomy_id"] for row in detections):
        raise ValueError(
            "Resolve detections against the current tissue taxonomy before classifying"
        )
    unresolved = tuple(
        i for i, row in enumerate(detections) if row.mapping_status not in {"mapped", "stage_only"}
    )
    matching = (
        set.intersection(*(set(row.scopes) for row in detections))
        if detections and not unresolved
        else set()
    )
    narrowest = set(matching)
    if matching:
        scopes = normal_tissue_scopes(species=detections[0].taxon_id)
        for name in matching:
            if any(scopes[other] < scopes[name] for other in matching):
                narrowest.discard(name)
    return {
        "taxonomy_id": _registry()[0]["taxonomy_id"],
        "status": "no_detections"
        if not detections
        else "unresolved"
        if unresolved
        else "compatible"
        if matching
        else "outside_permitted_scopes",
        "matching_scopes": tuple(sorted(matching)),
        "narrowest_scopes": tuple(sorted(narrowest)),
        "unresolved_detection_indices": unresolved,
        "absence_evidence": "not_assessed",
    }
