"""Reviewed gene findings, with tissue/state and source-specific limitations.

The reference separates observations from their screening consequences. A
literature search hit, homology, a developmental phenotype or an unquantified
annotation alone cannot establish adult normal production. These are research
adjudications, not estimates of clinical toxicity or peptide presentation.
"""

from __future__ import annotations

import json
import re
from importlib.resources import files

import pandas as pd

EFFECTS = {"exclude", "review", "reproductive", "context"}


def validate(records: list[dict]) -> None:
    """Reject incomplete provenance, identity and incompatible conclusions."""
    seen = set()
    for row in records:
        if row["id"] in seen:
            raise ValueError("Duplicate knowledge finding")
        seen.add(row["id"])
        if not row["gene_ids"] or any(not re.fullmatch(r"ENSG\d{11}", g) for g in row["gene_ids"]):
            raise ValueError("Knowledge requires exact Ensembl gene identities")
        if row["effect"] not in EFFECTS or row["review_status"] != "accepted":
            raise ValueError("Unreviewed leads cannot enter the knowledge reference")
        if not all(
            row.get(k)
            for k in [
                "finding",
                "tissue",
                "cell_type",
                "state",
                "assay",
                "source_title",
                "source_location",
                "limitations",
                "decision_reason",
            ]
        ):
            raise ValueError("Knowledge requires context, method, finding and limitations")
        if not row["source_url"].startswith("https://"):
            raise ValueError("Knowledge requires a source link")
        if row["effect"] == "exclude" and (
            row["taxonomy_id"] != 9606
            or row["stage"] != "adult"
            or row["tissue_tier"] not in {"vital", "protected"}
            or row["support"] not in {"normal_protein", "localized_rna_and_function"}
        ):
            raise ValueError("Exclusion requires established adult human normal production")
        if row["effect"] == "reproductive" and (
            row["taxonomy_id"] != 9606
            or row["tissue_tier"] != "reproductive"
            or row["support"] not in {"normal_protein", "localized_rna_and_function"}
        ):
            raise ValueError("Reproductive credit requires direct human production evidence")
        # pandas 2/3 differ in whether mixed columns retain None or use NaN.
        # Neither representation may turn an unquantified finding into a value.
        if bool(pd.isna(row.get("value"))) != bool(pd.isna(row.get("unit"))):
            raise ValueError("Measurements require native value and unit together")


def findings(gene: str | None = None, *, taxonomy_id: int = 9606) -> pd.DataFrame:
    """Return detached reviewed findings by exact gene ID or case-insensitive symbol.

    No findings means no curated evidence in this reference, not normal absence.
    Species are filtered explicitly; human findings never transfer by orthology.
    ``effect`` records the reviewed consequence: exclude, review, reproductive
    presence, or context only. Native units (including fold changes) are retained.
    """
    data = json.loads(files("oncoref").joinpath("data/gene-knowledge.json").read_text())
    if data["schema"] != "oncoref.gene-knowledge" or data["version"] != 1:
        raise ValueError("Unsupported gene knowledge schema")
    validate(data["findings"])
    rows = [r for r in data["findings"] if r["taxonomy_id"] == taxonomy_id]
    if gene is not None:
        rows = [
            r
            for r in rows
            if gene in r["gene_ids"] or gene.casefold() in {s.casefold() for s in r["gene_symbols"]}
        ]
    return pd.DataFrame(rows, columns=data["columns"])


def annotate(facts: pd.DataFrame, evidence: pd.DataFrame | None = None) -> pd.DataFrame:
    """Apply exact-member findings without copying evidence from a paralog.

    A normal exclusion/hold on any member applies to the identical reference
    protein. Reproductive presence requires evidence covering every member;
    explicitly reviewed protein-equivalent members can share a finding.
    """
    if evidence is None:
        evidence = findings()
    validate(evidence.to_dict("records"))
    d = facts.copy()
    members = d.member_gene_ids.map(lambda x: set(json.loads(x)) if isinstance(x, str) else set(x))
    if members.map(lambda group: not group).any():
        raise ValueError("Protein findings require at least one explicit member identity")
    human = evidence[evidence.taxonomy_id.eq(9606)]
    for effect, column in [
        ("exclude", "knowledge_normal_exclusion"),
        ("review", "knowledge_normal_hold"),
        ("reproductive", "knowledge_reproductive_support"),
    ]:
        ids = {g for genes in human.loc[human.effect.eq(effect), "gene_ids"] for g in genes}
        d[column] = members.map(
            lambda group, ids=ids, effect=effect: (
                group <= ids if effect == "reproductive" else bool(group & ids)
            )
        )
    d["knowledge_finding_ids"] = members.map(
        lambda group: json.dumps(
            sorted(r.id for r in human.itertuples() if group & set(r.gene_ids))
        )
    )
    return d
