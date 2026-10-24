"""Pinned, licensed UniProt annotations with explicit Ensembl gene mappings.

Annotations are biological source facts, not clinical outcome assessments or
screen overrides. Shared entries remain shared; symbols never substitute for
absent Ensembl cross-references. The source retains evidence codes, by-similarity
qualifiers, disease details and complete bibliographic references.
"""

from __future__ import annotations

import gzip
import json
from collections import defaultdict

import pandas as pd

from . import reference_data


def reviewed_human() -> list[dict]:
    """Fresh native Swiss-Prot entries from the shared, checksum-verified cache."""
    path = reference_data.ensure("uniprot_reviewed_human")
    with gzip.open(path, "rt") as handle:
        return json.load(handle)["results"]


def _pubmed_ids(value):
    if isinstance(value, dict):
        if value.get("source") == "PubMed" and value.get("id"):
            yield value["id"]
        for item in value.values():
            yield from _pubmed_ids(item)
    elif isinstance(value, list):
        for item in value:
            yield from _pubmed_ids(item)


def _gene_rows(entries):
    rows = []
    for entry in entries:
        gene_ids = set()
        for crossref in entry.get("uniProtKBCrossReferences", []):
            if crossref.get("database") == "Ensembl":
                for prop in crossref.get("properties", []):
                    if prop.get("key") == "GeneId":
                        gene_id = prop["value"].split(".")[0]
                        if gene_id.startswith("ENSG"):
                            gene_ids.add(gene_id)
        texts = defaultdict(list)
        evidence = defaultdict(set)
        for comment in entry.get("comments", []):
            kind = comment["commentType"]
            texts[kind].extend(t["value"] for t in comment.get("texts", []))
            if kind == "DISEASE":
                disease = comment.get("disease", {})
                texts[kind].append(
                    ": ".join(
                        v for v in [disease.get("diseaseId"), disease.get("description")] if v
                    )
                )
            if kind == "SUBCELLULAR LOCATION":
                texts[kind].extend(
                    loc["location"]["value"]
                    for loc in comment.get("subcellularLocations", [])
                    if "location" in loc
                )
            evidence[kind].update(_pubmed_ids(comment))
        description = entry.get("proteinDescription", {})
        names = description.get("recommendedName", description.get("submissionNames", [{}])[0])
        for gene_id in sorted(gene_ids):
            row = {
                "gene_id": gene_id,
                "accession": entry["primaryAccession"],
                "protein_name": names.get("fullName", {}).get("value", ""),
                "entry_gene_ids": "|".join(sorted(gene_ids)),
                "shared_entry": len(gene_ids) > 1,
                "protein_existence": entry.get("proteinExistence", ""),
                "url": "https://www.uniprot.org/uniprotkb/" + entry["primaryAccession"] + "/entry",
            }
            for kind, column in [
                ("FUNCTION", "function"),
                ("TISSUE SPECIFICITY", "tissue_specificity"),
                ("DISEASE", "disease"),
                ("PATHWAY", "pathway"),
                ("SUBCELLULAR LOCATION", "location"),
                ("CAUTION", "caution"),
            ]:
                row[column] = " ".join(dict.fromkeys(t for t in texts[kind] if t))
                row[column + "_pubmed_ids"] = "|".join(sorted(evidence[kind]))
            rows.append(row)
    return rows


def gene_annotations(gene_ids=None) -> pd.DataFrame:
    """One row per explicit gene/entry pair, preserving shared-entry ambiguity.

    Ensembl version suffixes are stripped, but genes are never remapped by symbol
    or paralog. A missing row means no exact cross-reference in this snapshot.
    """
    result = pd.DataFrame(_gene_rows(reviewed_human()))
    if gene_ids is not None:
        if isinstance(gene_ids, str):
            gene_ids = [gene_ids]
        result = result.loc[result.gene_id.isin(gene_ids)]
    return result.reset_index(drop=True)
