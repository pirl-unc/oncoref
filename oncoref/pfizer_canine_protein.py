# Licensed under the Apache License, Version 2.0.
"""Pfizer MSV000099744 canine bulk protein abundance.

Use `source()` to inspect provenance and acquisition status without a download.
"""

from __future__ import annotations

from . import expression_registry as _registry
from . import species_expression as _expression
from .expression_registry import ExpressionSource
from .species_assets import DEFAULT_MAX_BYTES
from .species_expression import ExpressionDataset

__all__ = ["read", "source"]


def source() -> ExpressionSource:
    """Return the bundled source record, including its measurement contract."""
    record = _registry.expression_source("pfizer-dog-protein", species=9615)
    if record is None:
        raise KeyError("Missing bundled canine source: pfizer-dog-protein")
    return record


def read(
    *,
    offline: bool = False,
    cache_root=None,
    local_overrides=None,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> ExpressionDataset:
    """Native protein-group iBAQ across tissue pools, not independent donors.

    Uses the pinned ORFdata dataset and shared managed cache. Requires
    oncoref[quantitative] on Python 3.11+. An override names
    ``9615:pfizer-atlas-dog-protein`` and the released dataset.sqlite file.
    """
    return _expression.read_species_expression(
        "pfizer-dog-protein",
        species=9615,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )


def assignment_contexts(expression, dataset, identifiers):
    """Resolve summary representatives through the producer's reconciled index.

    Summary group IDs and native numeric group IDs occupy different namespaces.
    This bridge retains that distinction and does not allocate group abundance.
    """
    from orfdata.sources.pfizer_attribution import CanineEvidence

    summary = expression.provenance["manifest"]
    source = next(row for row in summary["sources"] if row["asset_id"] == "protein.tsv")
    expected = dict(dataset.manifest.processing["protein_summary"])
    if (
        expected["dataset_id"] != summary["dataset_id"]
        or expected["version"] != summary["version"]
        or expected["source_sha256"] != source["sha256"]
        or dict(expected["source_counts"]) != summary["coverage"]["source_counts"]
    ):
        raise ValueError("Peptide evidence identifies a different protein summary")
    native = CanineEvidence(dataset.path)
    if native.search_reference.reference != dataset.manifest.references[0]:
        raise ValueError("Native searched reference disagrees with evidence manifest")
    contexts, targets = [], {}
    for identifier in identifiers:
        context = native.summary_group(identifier)
        if (
            context is None
            or context["baseline"] != expected
            or context["representative"] != expression._features.at[identifier, "Protein"]
        ):
            raise ValueError("Native summary context disagrees with expression feature")
        protein = native.search_protein(context["representative"])
        if protein is None:
            raise ValueError("Summary representative is absent from the pinned search reference")
        contexts.append(_portable_native(context))
        targets[protein.target.key] = protein.target
    return contexts, tuple(targets[key] for key in sorted(targets))


def _portable_native(value):
    """Lossless JSON for native BLOBs, nonfinite doubles and large integer IDs."""
    import base64
    import math
    import struct

    if isinstance(value, bytes):
        return {"__oncoref_native__": "bytes", "base64": base64.b64encode(value).decode("ascii")}
    if type(value) is float and not math.isfinite(value):
        return {"__oncoref_native__": "float64", "hex": struct.pack("<d", value).hex()}
    if type(value) is int and abs(value) > 2**53 - 1:
        return {"__oncoref_native__": "integer", "decimal": str(value)}
    if isinstance(value, dict):
        return {key: _portable_native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_portable_native(item) for item in value]
    return value


def assignment_records(dataset, targets):
    """Use public indexed queries, retaining distinct reported and compatible sets."""
    from dataclasses import replace

    from orfdata import ScopedTarget
    from orfdata.sources.pfizer_attribution import CanineEvidence

    native = CanineEvidence(dataset.path)
    groups, sequences, seen = set(), set(), set()

    def complete(mapping):
        if mapping.mapping_id in seen:
            return
        seen.add(mapping.mapping_id)
        provenance = dict(mapping.provenance)
        if "native_table" in provenance:
            table, key = provenance["native_table"], tuple(int(x) for x in provenance["native_key"])
            (record,) = native.records(table, key=key)
            links = [
                {"link": row["source_link"], "target": row["right"], "status": row["status"]}
                for row in native.linked(table, key, via=provenance["link_table"])
            ]
            provenance["oncoref_native_context"] = _portable_native(
                {"record": record, "reported_links": links}
            )
            mapping = replace(mapping, provenance=provenance)
        yield "mappings", mapping

    for target in targets:
        for mapping in dataset.mappings(target=target):
            yield from complete(mapping)
            if mapping.source_namespace == "amino-acid-sequence":
                sequences.add(mapping.source_identifier)
            elif mapping.provenance.get("native_table") == "TargetProteinGroups":
                groups.add((mapping.source_identifier, mapping.source_namespace))
    # A sequence-compatible representative need not occur in the reported set.
    # Query all records of each selected peptide to keep that disagreement visible.
    for sequence in sorted(sequences):
        for mapping in dataset.mappings(
            source_identifier=sequence, source_namespace="amino-acid-sequence"
        ):
            yield from complete(mapping)
    for identifier, namespace in sorted(groups):
        group = ScopedTarget(dataset.manifest.references[0], "protein_group", identifier, namespace)
        yield from (("relations", row) for row in dataset.relations(left_target=group))
