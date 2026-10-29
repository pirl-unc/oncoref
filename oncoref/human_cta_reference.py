# Licensed under the Apache License, Version 2.0.
"""Frozen human annotation mappings; historical assay references are separate."""

import gzip
import json
from importlib.resources import files

from .reference_identity import ReferenceIdentity, _identity
from .species_assets import registered_references


def human_cta_reference():
    """Read the verified sequence-backed mapping without downloads or indexes."""
    path = files("oncoref").joinpath("data/species/human/cta-reference-mappings.json.gz")
    record = json.loads(gzip.decompress(path.read_bytes()))
    if set(record) != {
        "schema",
        "schema_version",
        "reference",
        "mappings",
        "sequence_audit",
        "mapping_id",
    } or (record["schema"], record["schema_version"]) != ("oncoref.cta-reference-mappings", 1):
        raise ValueError("Unsupported human CTA reference mapping")
    identity = record.pop("mapping_id")
    if identity != _identity("cta-reference-mappings", 1, record):
        raise ValueError("Human CTA reference mapping digest mismatch")
    registration = next(r for r in registered_references("human") if r["id"] == "ensembl115")
    reference = ReferenceIdentity.from_dict(record["reference"])
    if reference != ReferenceIdentity.from_dict(registration["identity"]):
        raise ValueError("Human CTA mapping and registered reference disagree")
    for gene_id, mapping in record["mappings"].items():
        if mapping["reference_id"] != reference.reference_id or any(
            target["id"] != gene_id for target in mapping["targets"]
        ):
            raise ValueError("Human CTA mapping endpoint mismatch")
    return reference, record["mappings"]
