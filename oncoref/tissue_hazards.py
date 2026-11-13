# Licensed under the Apache License, Version 2.0.
"""Central, evidence-linked adult human tissue and cell-compartment hazard policy.

Functional consequence is independent of reproductive membership, expression and
receptor recognition. No RNA cutoff or clinical probability is inferred here.
"""

from __future__ import annotations

import hashlib
import re
from copy import deepcopy
from functools import cache
from importlib.resources import files

import pandas as pd
import yaml

TIERS = ("vital", "protected", "manageable", "disposable")
CELL_SOURCES = frozenset(
    {
        "hpa_single_cell_tissue",
        "hpa_normal_tissue",
        "hpa_immune_cell",
        "hpa_immune_monaco",
        "cornea_collin_adult_counts",
        *[
            "dice_" + x
            for x in (
                "b_cell_naive",
                "monocytes",
                "m2",
                "nk",
                "treg_mem",
                "cd4_naive",
                "cd4_stim",
                "treg_naive",
                "tfh",
                "th1",
                "thstar",
                "th17",
                "th2",
                "cd8_naive",
                "cd8_stim",
            )
        ],
    }
)


@cache
def _load():
    payload = files("oncoref").joinpath("data/normal-tissue-hazards.yaml").read_bytes()
    doc = yaml.safe_load(payload)
    validate_registry(doc)
    doc["registry_id"] = "oncoref:tissue-hazards:v1:sha256:" + hashlib.sha256(payload).hexdigest()
    return doc


def validate_registry(doc: dict) -> None:
    """Reject incomplete, ambiguous or uncited policy records before use."""
    if doc.get("schema") != "oncoref.normal-tissue-hazards" or doc.get("schema_version") != 1:
        raise ValueError("Unsupported tissue-hazard registry schema")
    if set(doc["tier_definitions"]) != set(TIERS):
        raise ValueError("Incomplete hazard tier definitions")
    required = {
        "id",
        "label",
        "level",
        "tier",
        "system",
        "tissues",
        "mechanisms",
        "acute_loss",
        "persistent_loss",
        "repair",
        "management",
        "rationale",
        "references",
        "inference",
        "evidence_scope",
        "recognition_evidence",
        "taxon_ids",
        "stage",
        "priority",
    }
    ids, bulk = set(), set()
    for r in doc["assessments"]:
        if not required <= r.keys() or any(not r[k] for k in required):
            raise ValueError(f"Incomplete hazard assessment: {r.get('id')}")
        if r["id"] in ids or r["tier"] not in TIERS or r["level"] not in {"bulk", "cell"}:
            raise ValueError(f"Invalid hazard identity or tier: {r['id']}")
        ids.add(r["id"])
        if r["taxon_ids"] != [9606] or r["stage"] != "adult":
            raise ValueError("Current hazard evidence policy is adult human only")
        if not set(r["references"]) <= doc["references"].keys():
            raise ValueError(f"Unknown hazard reference: {r['id']}")
        if not set(r["mechanisms"]) <= doc["mechanism_definitions"].keys():
            raise ValueError(f"Unknown hazard mechanism: {r['id']}")
        if r["level"] == "cell":
            re.compile(r["cell_pattern"])
        else:
            if bulk & set(r["tissues"]):
                raise ValueError("Overlapping bulk hazard roots")
            bulk.update(r["tissues"])
    for key, reference in doc["references"].items():
        if not all(
            reference.get(k)
            for k in ("title", "model", "design", "finding", "limitation", "checked")
        ):
            raise ValueError(f"Incomplete hazard reference: {key}")
        if not reference.get("url", "").startswith("https://"):
            raise ValueError(f"Unlinked hazard reference: {key}")


def registry() -> dict:
    """Detached policy, assessments, bibliography and content fingerprint."""
    return deepcopy(_load())


def assessments() -> pd.DataFrame:
    """One row per reviewed bulk or cell rule; nested lists are preserved."""
    return pd.DataFrame(deepcopy(_load()["assessments"]))


def references() -> pd.DataFrame:
    """Evidence design, biological model, supported finding and limitations."""
    return (
        pd.DataFrame.from_dict(deepcopy(_load()["references"]), orient="index")
        .rename_axis("reference_id")
        .reset_index()
    )


def screening_roots() -> dict[str, list[str]]:
    """Current somatic bulk tiers, excluding separate permitted-scope membership."""
    return {
        tier: sorted(
            t
            for r in _load()["assessments"]
            if r["level"] == "bulk" and r["tier"] == tier and r["system"] != "Reproductive"
            for t in r["tissues"]
        )
        for tier in TIERS
    }


@cache
def _anatomy(tissue: str):
    from .normal_tissues import _ancestors, resolve_normal_tissue

    a = resolve_normal_tissue(tissue, species="human", stage="adult")
    return a, (_ancestors(a.concept_id) if a.concept_id else set())


@cache
def _resolve(tissue: str, cell_type: str, taxon_id: int, stage: str, resolution: str):
    result = {
        "tier": "unresolved",
        "assessment_id": None,
        "matched_level": None,
        "reason": "unresolved_anatomy",
        "registry_id": _load()["registry_id"],
        "references": [],
        "anatomical_membership": [],
        "cell_type": cell_type,
        "tissue": tissue,
        "taxon_id": taxon_id,
        "stage": stage,
        "resolution": resolution,
    }
    if taxon_id != 9606 or stage != "adult":
        result["reason"] = "outside_adult_human_policy"
        return result
    annotation, ancestors = _anatomy(tissue)
    result["anatomical_membership"] = list(annotation.scopes)
    if not ancestors:
        return result
    bulk = [
        r
        for r in _load()["assessments"]
        if r["level"] == "bulk" and ancestors.intersection(r["tissues"])
    ]
    # Explicit anatomical child overrides its parent (e.g. hair follicle / skin).
    exact = [r for r in bulk if annotation.concept_id in r["tissues"]]
    if exact:
        bulk = exact
    if len(bulk) != 1:
        result["reason"] = "unreviewed_anatomy" if not bulk else "ambiguous_anatomy"
        return result
    chosen = bulk[0]
    reason = "bulk_rule" if resolution == "bulk" or not cell_type else "organ_fallback"
    # Mixed population labels cannot establish restricted exposure. Keep the
    # broad organ's concern instead of picking one matching substring.
    mixed = bool(re.search(r"\b(and|mixed|mixture)\b|[+/;]", cell_type))
    if resolution == "cell" and cell_type:
        candidates = [
            r
            for r in _load()["assessments"]
            if r["level"] == "cell"
            and ancestors.intersection(r["tissues"])
            and re.search(r["cell_pattern"], cell_type)
        ]
        if mixed and candidates:
            chosen = min([chosen, *candidates], key=lambda r: TIERS.index(r["tier"]))
            reason = "mixed_cell_label_conservative_fallback"
        elif candidates:
            priority = max(r["priority"] for r in candidates)
            best = [r for r in candidates if r["priority"] == priority]
            if len(best) > 1:
                # Never choose a lower tier by incidental registry ordering.
                result["reason"] = "ambiguous_cell_rule"
                result["candidate_assessments"] = sorted(r["id"] for r in best)
                result["fallback_tier"] = min(
                    [chosen["tier"], *(r["tier"] for r in best)], key=TIERS.index
                )
                return result
            chosen, reason = best[0], "cell_rule"
    if mixed and reason != "mixed_cell_label_conservative_fallback":
        reason = "mixed_cell_label_organ_fallback"
    result.update(
        tier=chosen["tier"],
        assessment_id=chosen["id"],
        matched_level=chosen["level"],
        reason=reason,
        references=chosen["references"],
    )
    return result


def assess(
    tissue: str,
    cell_type: str | None = None,
    *,
    taxon_id: int = 9606,
    stage: str = "adult",
    resolution: str = "bulk",
) -> dict:
    """Resolve a hazard with provenance; cell-specific calls require resolution='cell'.

    Unknown cell labels inherit their organ. Unknown anatomy, conflicting rules,
    children and nonhuman species are unresolved rather than implicitly permitted.
    Reproductive membership does not overwrite the independently returned hazard.
    """
    if resolution not in {"bulk", "cell"}:
        raise ValueError("Resolution must be bulk or cell")

    def normalize(x):
        return " ".join(str(x or "").casefold().strip().split())

    cell = re.sub(r"[_-]+", " ", normalize(cell_type))
    cell = _load().get("cell_aliases", {}).get(cell, cell)
    return deepcopy(_resolve(normalize(tissue), cell, taxon_id, normalize(stage), resolution))


def context_tiers(contexts: pd.DataFrame) -> pd.Series:
    """Apply the central policy once per unique measured tissue/cell/source context.

    Existing reproductive/developmental membership stays independent. Unknown
    anatomy retains the incoming tier and cannot earn an undocumented downgrade.
    """
    tiers = contexts.category.replace({"renewing_mucosa": "protected"}).copy()
    eligible = contexts.mapping_status.eq("mapped") & contexts.category.isin(
        ["vital", "protected", "manageable", "disposable", "renewing_mucosa"]
    )
    if "taxonomy_id" in contexts:
        eligible &= contexts.taxonomy_id.eq(9606)
    if "stage" in contexts:
        eligible &= contexts.stage.eq("adult")
    keys = pd.DataFrame(
        {
            "tissue": contexts.tissue.fillna(""),
            "cell_type": contexts.cell_type.fillna(""),
            "cellular": contexts.source.isin(CELL_SOURCES),
        }
    )
    mapping = {}
    for t, c, cellular in keys.loc[eligible].drop_duplicates().itertuples(index=False, name=None):
        finding = assess(t, c, resolution="cell" if cellular else "bulk")
        mapping[(t, c, cellular)] = finding.get("fallback_tier", finding["tier"])
    assigned = pd.Series(
        [mapping.get(key, "unresolved") for key in keys.itertuples(index=False, name=None)],
        index=contexts.index,
    )
    # Disposable overrides require resolved narrow anatomy/cell identity above.
    valid = eligible & assigned.isin(TIERS)
    tiers.loc[valid] = assigned.loc[valid]
    return tiers


def anatomy_inventory() -> pd.DataFrame:
    """Human anatomy with separate selection scope and adult injury assessment.

    Reproductive includes placental, extraembryonic, embryonic and fetal
    contexts. Broad location labels require localization; they are neither a
    fifth hazard tier nor a Disposable allowance. This inventory describes
    canonical concepts, not completeness of measured expression coverage.
    """
    from .normal_tissues import normal_tissue_taxonomy, resolve_normal_tissue

    rows = []
    for cid, concept in normal_tissue_taxonomy()["concepts"].items():
        if 9606 not in concept.get("taxon_ids", [9606, 9615, 10090]):
            continue
        stage = {"embryo": "embryonic", "fetus": "fetal"}.get(cid, "adult")
        annotation = resolve_normal_tissue(cid, species="human", stage=stage)
        r = assess(cid, stage=stage)
        reproductive = "reproductive" in annotation.scopes
        broad = cid in {"soft_tissue", "head", "foot"}
        rows.append(
            {
                "concept_id": cid,
                "label": concept["label"],
                **r,
                "anatomical_membership": list(annotation.scopes),
                "selection_category": "Reproductive"
                if reproductive
                else r["tier"].title()
                if r["tier"] in TIERS
                else None,
                "context_kind": "broad_location" if broad else "specific_context",
                "requires_localization": broad,
                "requires_hazard_review": not reproductive and not broad and r["tier"] not in TIERS,
            }
        )
    return pd.DataFrame(rows)
