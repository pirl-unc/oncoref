"""Reviewed TCR-T experience and its limited protein-level RNA screening credit.

Clinical evidence belongs to a peptide, HLA, receptor, regimen and population.
The credit is a research-policy inference to a protein's aggregated RNA profile,
not clinical clearance of every isoform, peptide, receptor or patient. Sentinel
events are not an exhaustive adverse-event table. Unreviewed is never negative.
"""

from __future__ import annotations

import json
from copy import deepcopy
from importlib.resources import files

import pandas as pd

EVIDENCE_FILE = "tcr-trial-evidence.json"
SPECIFICITY_FINDINGS = {"none_identified", "suspected_or_confirmed", "unassessed"}
SPECIFICITY_ATTRIBUTIONS = {
    "suspected",
    "confirmed",
    "not_supported",
    "not_established",
    "unresolved",
    "unassessed",
}


def validate(record: dict) -> None:
    """Reject ambiguous identities, denominators and contradictory safety reviews."""
    if (record.get("schema"), record.get("schema_version")) != ("oncoref.tcr-trial-evidence", 1):
        raise ValueError("Unsupported TCR trial evidence schema")
    products = {p["product_id"]: p for p in record["products"]}
    cohorts = {c["cohort_id"]: c for c in record["cohorts"]}
    if len(products) != len(record["products"]) or len(cohorts) != len(record["cohorts"]):
        raise ValueError("Duplicate product or cohort identity")
    if len({e["event_id"] for e in record["events"]}) != len(record["events"]):
        raise ValueError("Duplicate event identity")
    sources = record["sources"]
    for p in products.values():
        if p["modality"] != "TCR-T" or p["taxon_id"] != 9606:
            raise ValueError("This curation is scoped to human TCR-T")
        if not p["hla"] or not p["targets"] or not p["identity_sources"]:
            raise ValueError("Missing peptide/HLA/product identity")
        if not set(p["identity_sources"]) <= sources.keys():
            raise ValueError("Unknown identity source")
        for target in p["targets"]:
            if not all(target.get(k) for k in ["gene_id", "symbol", "peptide", "isoform_scope"]):
                raise ValueError("Incomplete target identity")
    for c in cohorts.values():
        if c["product_id"] not in products or c["source"] not in sources:
            raise ValueError("Unknown cohort product/source")
        n = c["patients_treated"]
        if n is not None and (type(n) is not int or n <= 0):
            raise ValueError("Treated denominator must be positive or unknown")
        if c["grade45_specificity"] not in SPECIFICITY_FINDINGS:
            raise ValueError("Unknown specificity finding")
        if c["review_status"] not in {"reviewed", "pending"}:
            raise ValueError("Unknown review status")
        if c["grade45_specificity"] == "none_identified" and (
            c["review_status"] != "reviewed" or n is None or not c["source_sections"]
        ):
            raise ValueError("Favorable evidence requires reviewed results and denominator")
    for e in record["events"]:
        if e["product_id"] not in products or e["source"] not in sources:
            raise ValueError("Unknown event product/source")
        if e["specificity_attribution"] not in SPECIFICITY_ATTRIBUTIONS:
            raise ValueError("Unknown specificity attribution")
        if e["grade"] is not None and (type(e["grade"]) is not int or not 1 <= e["grade"] <= 5):
            raise ValueError("Grade must be 1-5 or unreported")
        cohort = e["cohort_id"]
        if cohort is not None:
            if cohort not in cohorts or cohorts[cohort]["product_id"] != e["product_id"]:
                raise ValueError("Event/cohort product mismatch")
            if (
                _specificity_block(e)
                and cohorts[cohort]["grade45_specificity"] == "none_identified"
            ):
                raise ValueError("Favorable cohort contradicts a specificity event")


def _specificity_block(event: dict) -> bool:
    # An ungraded suspected specificity injury cannot masquerade as below grade 4.
    return event["specificity_attribution"] in {"suspected", "confirmed"} and (
        event["grade"] is None or event["grade"] >= 4
    )


def evidence() -> dict:
    """Return a validated, detached snapshot with linked primary sources, offline."""
    record = json.loads(files("oncoref").joinpath("data", EVIDENCE_FILE).read_text())
    validate(record)
    return record


def products(*, record: dict | None = None) -> pd.DataFrame:
    """One row per reviewed receptor/product, never one observation per gene alias.

    Favorable credit needs a reviewed treated cohort with no identified grade 4-5
    specificity signal, and no such signal in any reviewed cohort/event for that
    receptor. A pending cohort cannot supply credit; an explicitly unassessed
    specificity event blocks credit until reviewed. Treatment-related CRS,
    conditioning toxicity and events with no established specificity mechanism
    remain visible without being renamed antigen-specific tissue destruction.
    """
    record = evidence() if record is None else deepcopy(record)
    validate(record)
    rows = []
    for p in record["products"]:
        cohorts = [c for c in record["cohorts"] if c["product_id"] == p["product_id"]]
        events = [e for e in record["events"] if e["product_id"] == p["product_id"]]
        positive = any(c["grade45_specificity"] == "suspected_or_confirmed" for c in cohorts)
        blocked = positive or any(_specificity_block(e) for e in events)
        unresolved = any(
            e["specificity_attribution"] in {"unassessed", "unresolved"} for e in events
        )
        favorable = [
            c["cohort_id"]
            for c in cohorts
            if c["review_status"] == "reviewed" and c["grade45_specificity"] == "none_identified"
        ]
        credit = bool(favorable) and not blocked and not unresolved
        rows.append(
            p
            | {
                "rna_credit": credit,
                "credit_cohorts": tuple(sorted(favorable)) if credit else (),
                "clinical_finding": "specificity_signal"
                if blocked
                else "unassessed"
                if unresolved or not favorable
                else "no_identified_severe_specificity_signal",
                "event_ids": tuple(e["event_id"] for e in events),
            }
        )
    return pd.DataFrame(rows)


def cohorts() -> pd.DataFrame:
    """Reported cohort denominators; do not add repeated publications or cases."""
    return pd.DataFrame(evidence()["cohorts"])


def events() -> pd.DataFrame:
    """Sentinel clinical events, including ungraded and unresolved attribution."""
    return pd.DataFrame(evidence()["events"])


def protein_credit(gene_ids, *, taxon_id: int = 9606, record: dict | None = None) -> dict:
    """Assess modest protein-profile credit; every member locus must be supported.

    No symbol-prefix/family/orthology transfer. Multiple trials never stack.
    Returning credit does not clear other peptides, untested isoforms or a new
    receptor. Use ``products`` and its exact target mappings for peptide work.
    """
    ids = set(gene_ids)
    if not ids or isinstance(gene_ids, str):
        raise ValueError("Supply a nonempty collection of gene IDs")
    if taxon_id != 9606:
        return {
            "rna_credit": False,
            "product_ids": (),
            "cohort_ids": (),
            "scope": "no_species_transfer",
        }
    p = products(record=record)
    selected = []
    covered = set()
    for row in p.to_dict("records"):
        matches = ids & {t["gene_id"] for t in row["targets"]}
        if row["rna_credit"] and matches:
            selected.append(row)
            covered |= matches
    credit = ids <= covered
    return {
        "rna_credit": credit,
        "product_ids": tuple(sorted(row["product_id"] for row in selected)) if credit else (),
        "cohort_ids": tuple(sorted({c for row in selected for c in row["credit_cohorts"]}))
        if credit
        else (),
        "scope": "protein_profile_only_tested_epitope_limits_apply"
        if credit
        else "unassessed_or_not_qualified",
    }


def annotate_proteins(facts: pd.DataFrame) -> pd.DataFrame:
    """Attach credit and provenance to the replay's reference-bound member sets."""
    record = evidence()
    supported = {t["gene_id"] for p in record["products"] for t in p["targets"]}
    d = facts.copy()
    credits = []
    for members in d.member_gene_ids:
        ids = json.loads(members) if isinstance(members, str) else members
        credits.append(
            protein_credit(ids, record=record)
            if supported.intersection(ids)
            else {"rna_credit": False, "product_ids": (), "cohort_ids": (), "scope": "unreviewed"}
        )
    d["clinical_rna_credit"] = [c["rna_credit"] for c in credits]
    d["clinical_product_ids"] = [";".join(c["product_ids"]) for c in credits]
    d["clinical_cohort_ids"] = [";".join(c["cohort_ids"]) for c in credits]
    d["clinical_credit_scope"] = [c["scope"] for c in credits]
    return d
