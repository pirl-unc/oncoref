# Licensed under the Apache License, Version 2.0.
"""Adapt existing human CTA facts to the shared species evidence API.

This is a lossless wrapper of the historical summary, not a new human curation.
The original flat APIs and their default-panel decisions remain unchanged.
"""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache

from .expression_contract import ExpressionContract, ExpressionMeasurement
from .human_cta_reference import human_cta_reference
from .load_dataset import _register_derived_cache
from .reference_identity import _canonical_json, _identity
from .species_antigens import CTAEvidenceSnapshot
from .version import __version__


@lru_cache(maxsize=1)
def human_cta_snapshot():
    from .cta import cta_evidence

    # Hash logical records, including the historical human-specific curation.
    # This also honors the existing dataset cache invalidation/override hooks.
    rows = json.loads(cta_evidence().to_json(orient="records", double_precision=15))
    reference, mappings = human_cta_reference()
    source_hash = hashlib.sha256(_canonical_json(rows).encode("ascii")).hexdigest()
    source_id = "oncoref-human-cta"
    version = f"oncoref-{__version__}-human-cta-summary"
    contract = ExpressionContract(
        9606,
        source_id,
        version,
        "mixed",
        None,
        "legacy_human",
        "Ensembl gene",
        ExpressionMeasurement("gene", None, "mixed native summary fields", None),
        {"source": source_hash},
        "Historical human CTA summary. RNA *_ntpm fields retain nTPM; protein fields retain categorical IHC calls; fractions and curation flags retain their original meaning. A summary is not one expression measurement. Quantification-reference provenance is legacy_human.",
    )
    source = {
        "id": source_id,
        "version": version,
        "url": "https://github.com/pirl-unc/oncoref/blob/main/oncoref/cta.py",
        "doi": None,
        "pmid": None,
        "sha256": source_hash,
        "scope": "SHA-256 of canonical JSON logical records returned by the historical cta_evidence() accessor. All original fields are retained in evidence.source_fields. This is a derived summary receipt, not an upstream publication checksum.",
        "study_id": "oncoref-human-cta-summary",
        "lineage": ["oncoref-human-cta-summary", "HPA-normal-tissue-summary", "HPA:v23"],
    }
    payload = {
        "schema": "oncoref.cta-evidence",
        "schema_version": 1,
        "taxon_id": 9606,
        "source_version": version,
        "coverage": {
            "status": "partial",
            "limitations": [
                "Lossless adapter of the existing human CTA summary; patient-level admission is not supplied.",
                "Summary fields mix native nTPM, categorical IHC, ratios and curation decisions. They are not interchangeable measurements.",
                "Sample/donor counts, study independence and complete upstream annotation provenance are not reconstructed from the summary.",
                "Ensembl 115 mapping enumerates pinned annotation RNA/protein alternatives. Noncoding transcripts retain null protein IDs. Sequence mapping does not establish equivalence to historical assay references or measured isoform abundance.",
                "Nomination provenance remains legacy_human; source_databases and specificity decisions are retained verbatim in source_fields.",
            ],
        },
        "references": [{"id": "ensembl115", "identity": reference.as_dict()}],
        "sources": [source],
        "contracts": [contract.as_dict()],
        "candidates": [],
        "evidence": [],
    }
    for row in rows:
        gid, label = row["Ensembl_Gene_ID"], row["Symbol"]
        eid = f"human-summary:{gid}"
        mapping = mappings.get(
            gid,
            {
                "reference_id": reference.reference_id,
                "gene_namespace": "Ensembl gene",
                "status": "unresolved",
                "method": "Exact lookup in the bundled sequence-backed CTA annotation inventory",
                "complete": False,
                "targets": [],
                "notes": "Candidate is outside the verified annotation inventory; no symbol or historical canonical-transcript fallback is applied.",
            },
        )
        payload["candidates"].append(
            {
                "id": gid,
                "label": label,
                "nomination_kind": "legacy_human",
                "evidence_status": "limited",
                "derivation": [
                    {
                        "kind": "legacy_human",
                        "source_id": source_id,
                        "locator": f"cta_evidence(): {gid}",
                        "note": "Historical candidate provenance retained in source_fields.source_databases; not reclassified by this adapter.",
                    }
                ],
                "mappings": [mapping],
                "evidence_ids": [eid],
                "notes": "Historical human curation is preserved; these rules are never applied to other species.",
            }
        )
        payload["evidence"].append(
            {
                "id": eid,
                "candidate_ids": [gid],
                "source_id": source_id,
                "study_id": source["study_id"],
                "contract_id": contract.contract_id,
                "locator": f"cta_evidence(): {gid}",
                "assay": "Historical HPA RNA/protein and CTA curation summary",
                "tissue": None,
                "tissue_role": "mixed",
                "sample_ids": [],
                "sample_count": None,
                "donor_count": None,
                "technical_replicates": None,
                "observation": {
                    "status": "qualitative",
                    "value": None,
                    "reported": "Native summary fields preserved in source_fields.",
                    "positive_count": None,
                },
                "source_fields": row,
                "notes": "No missing measurements or donor counts are inferred. See the expression contract for field units and source_fields for the original curation decisions.",
            }
        )
    return CTAEvidenceSnapshot.from_dict(
        {**payload, "snapshot_id": _identity("cta-evidence", 1, payload)}
    )


_register_derived_cache(human_cta_snapshot.cache_clear)
