#!/usr/bin/env python3
"""Rebuild the sparse P1A seed from checked primary-source bytes, offline.

Curated statements refer only to endogenous RNA in Sarma 1999 Figure 1 and
its accompanying results. Transgenic experiments and presentation are excluded.
See docs/species-cta-evidence.md for source acquisition and identity scope.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from oncoref.cta_evidence_adapters import normalize_cta_snapshot
from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.reference_identity import ReferenceIdentity
from oncoref.species_antigens import CTAEvidenceSnapshot

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "oncoref/data/species/mouse"
REFERENCE_ALIAS = "p1a-refseq-rs2024-02"
VERSION = "2026-10-29.1"


def build(source_dir: Path) -> CTAEvidenceSnapshot:
    """Verify every curated input; preserve assay and target reference separately."""
    receipt = json.loads((DATA / "cta-curation.json").read_text())
    for source in receipt["inputs"]:
        path = source_dir / source["name"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Source checksum mismatch: {path}")
    pins = {row["name"]: row["sha256"] for row in receipt["inputs"]}
    reference = ReferenceIdentity(
        taxon_id=10090,
        assembly_accession="GCF_000001635.27",
        annotation_source="NCBI RefSeq",
        annotation_release="RS_2024_02",
        source_version="NCBI Gene 22037 single-gene subset retrieved 2026-10-29",
        identifier_namespaces={"gene": "NCBI Gene", "transcript": "RefSeq", "protein": "RefSeq"},
        content_provenance={
            "gene_record_subset": pins["gene-22037.json"],
            "product_report_subset": pins["gene-22037-products.json"],
            "rna": pins["NM_011635.2.gb"],
            "protein": pins["NP_035765.1.gb"],
        },
    )
    sources = [
        {
            "id": "vandeneynde1991",
            "version": "Europe PMC full-text XML retrieved 2026-10-29",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2190831/",
            "doi": "10.1084/jem.173.6.1373",
            "pmid": "1903428",
            "sha256": pins["vandeneynde1991.xml"],
            "scope": "P1A nomination and native M36387.1/AAA40475.1 identity from DBA/2 P815 clone 211; expression figures are not digitized.",
            "study_id": "PMID:1903428",
            "lineage": ["PMID:1903428", "P1A:P815-clone-211"],
        },
        {
            "id": "sarma1999",
            "version": "Europe PMC full-text XML retrieved 2026-10-29",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2192944/",
            "doi": "10.1084/jem.189.5.811",
            "pmid": "10049945",
            "sha256": pins["sarma1999.xml"],
            "scope": "Figure 1 and accompanying normal-tissue results: endogenous qualitative RNA only; no transgenic, tolerance or presentation results imported.",
            "study_id": "PMID:10049945",
            "lineage": ["PMID:10049945", "P1A:normal-tissue-RNA"],
        },
    ]
    payload = {
        "schema": "oncoref.cta-evidence",
        "schema_version": 1,
        "taxon_id": 10090,
        "source_version": VERSION,
        "references": [{"id": REFERENCE_ALIAS, "identity": reference.as_dict()}],
        "sources": sources,
        "contracts": [],
        "candidates": [
            {
                "id": "P1A",
                "label": "P1A (Trap1a)",
                "nomination_kind": "published_nomination",
                "derivation": [
                    {
                        "kind": "published_nomination",
                        "source_id": "vandeneynde1991",
                        "locator": "Abstract; P1A cDNA clone 211 and accession M36387",
                        "note": "Native DBA/2 P815 rejection antigen; no mouse-wide restriction claim.",
                    },
                    {
                        "kind": "assay_followup",
                        "source_id": "sarma1999",
                        "locator": "Figure 1 and accompanying results",
                        "note": "Normal somatic RNA counterevidence and J558 tumor-cell expression.",
                    },
                ],
                "evidence_status": "conflicting",
                "mappings": [
                    {
                        "reference_id": reference.reference_id,
                        "gene_namespace": "NCBI Gene",
                        "status": "unique",
                        "method": "Native AAA40475.1 protein exactly matches NP_035765.1; NCBI Gene 22037 product report binds NM_011635.2/NP_035765.1 to GRCm39 NC_000086.8. All products in this pinned single-gene subset retained.",
                        "complete": True,
                        "targets": [
                            {
                                "id": "22037",
                                "symbol": "Trap1a",
                                "sequences": [
                                    {"transcript_id": "NM_011635.2", "protein_id": "NP_035765.1"}
                                ],
                            }
                        ],
                        "notes": "Reference assembly is C57BL/6J; RefSeq RNA source strain is C57BL/6. Native P815 source is DBA/2 with unspecified substrain. Proteins match (224 aa), but RNA/CDS differ. This is an annotation association, not measured reference-strain or isoform expression. See cta-curation.json for checked asset scope.",
                    }
                ],
                "evidence_ids": [],
                "notes": "A sparse nomination with positive normal-tissue RNA counterevidence. No tumor restriction, safety or vaccine admission is inferred.",
            }
        ],
        "evidence": [],
        "coverage": {
            "status": "partial",
            "limitations": [
                "One P1A nomination; not a mouse CTA panel, genome-wide screen or quantitative admission policy.",
                "Only Sarma 1999 Figure 1 and accompanying endogenous RNA results are imported. Protein abundance, placental expression, presentation and transgenic experiments are outside this seed.",
                "Five normal organs from BALB/c mice do not establish a complete normal panel; strain, stage, cell type and donor coverage remain limited.",
                "Donor/sample counts, biological independence and technical replicate counts are unreported for these imported observations; no counts are inferred from figure lanes.",
                "Qualitative assay signals and source-relative comparisons are not TPM, protein abundance or isoform-specific measurements. Less sensitive Northern results do not negate sensitive RT-PCR signals.",
                "Reference assets are an explicitly scoped NCBI Gene 22037 subset, not a full annotation release. Historical assay-to-reference sequence equivalence remains unresolved.",
            ],
        },
    }
    source = sources[1]
    observations = [
        (
            "testis",
            "reproductive",
            "normal",
            "BALB/c",
            None,
            "Detected; highest RNA signal among the normal tissues examined.",
            "Detected reproducibly by Northern; other normal organs lacked reproducible significant Northern signal.",
        ),
        (
            "lung",
            "somatic",
            "normal",
            "BALB/c",
            None,
            "Detected; RNA signal greater than liver in the reported comparison.",
            None,
        ),
        (
            "spleen",
            "somatic",
            "normal",
            "BALB/c",
            None,
            "Detected; RNA signal greater than liver and lower than J558 tumor cells.",
            None,
        ),
        (
            "liver",
            "somatic",
            "normal",
            "BALB/c",
            None,
            "Detected at a lower RNA signal than lung and spleen.",
            None,
        ),
        (
            "thymus",
            "somatic",
            "normal",
            "BALB/c",
            None,
            "RNA detected in some experiments; intermittent detection is not a measured zero.",
            None,
        ),
        (
            "J558 tumor cells",
            "tumor",
            "cell_line",
            None,
            "J558",
            "Detected; RNA signal greater than normal spleen in the reported comparison.",
            None,
        ),
    ]
    for tissue, role, context, strain, cell_line, reported, northern in observations:
        contract = ExpressionContract(
            10090,
            source["id"],
            source["version"],
            context,
            None,
            "unresolved",
            "study assay target P1A",
            ExpressionMeasurement("gene", None, "source_native", False),
            {
                "source": source["sha256"],
                "nomination_sequence": pins["M36387.1.gb"],
                "curation_receipt": hashlib.sha256(
                    (DATA / "cta-curation.json").read_bytes()
                ).hexdigest(),
            },
            "Historical assay reference and measured isoform are unresolved. Reference mapping does not transfer expression between strains.",
        )
        if contract.as_dict() not in payload["contracts"]:
            payload["contracts"].append(contract.as_dict())
        eid = f"sarma1999:{len(payload['evidence']) + 1:03d}"
        payload["evidence"].append(
            {
                "id": eid,
                "candidate_ids": ["P1A"],
                "source_id": source["id"],
                "study_id": source["study_id"],
                "contract_id": contract.contract_id,
                "locator": "Figure 1 and Results: Tumor Antigen P1A Is Expressed at Low Levels in Normal Tissue",
                "assay": "RT-PCR with Southern detection",
                "tissue": tissue,
                "tissue_role": role,
                "sample_ids": [],
                "sample_count": None,
                "donor_count": None,
                "technical_replicates": None,
                "observation": {
                    "status": "qualitative",
                    "value": None,
                    "reported": reported,
                    "positive_count": None,
                },
                "source_fields": {
                    "host_taxon_id": 10090,
                    "strain": strain,
                    "substrain": None,
                    "cell_line": cell_line,
                    "expression_origin": "endogenous",
                    "measured_isoform": None,
                    "nomination_identity": {
                        "role": "Nomination sequence comparison, not the measured assay isoform",
                        **receipt["identity_audit"],
                    },
                    "sample_ids_reported": False,
                    "northern_context": northern,
                    "relative_comparison": {
                        "reference": "J558 tumor cells",
                        "reported_fold_lower_range": [100, 1000],
                    }
                    if tissue == "spleen"
                    else {"reference": "spleen", "reported_fold_lower_minimum": 10}
                    if tissue == "thymus"
                    else None,
                },
                "notes": "Same study/assay lineage; figure comparisons are not independent donor replication. Strain and substrain are not inferred for J558. Quantitative titration bounds remain source-relative fields, not absolute measurements.",
            }
        )
        payload["candidates"][0]["evidence_ids"].append(eid)
    return normalize_cta_snapshot(CTAEvidenceSnapshot.create(payload))


def refresh_tables():
    """Refresh packaged native views without replaying historical human decisions."""
    from build_cta_selection_audit import TABLE_COLUMNS, membership_tables, native_assessment
    from build_cta_selection_release import mouse_selection_metadata, source_inventory

    package = DATA.parent.parent
    inventory = package / "cta-selection-sources.csv"
    source_inventory().to_csv(inventory, index=False)
    path = package / "cta-selection.json"
    selection = json.loads(path.read_text())
    selection["data_sha256"][inventory.name] = hashlib.sha256(inventory.read_bytes()).hexdigest()
    selection["species"]["mouse"] = mouse_selection_metadata()
    selection["inputs"]["oncoref/data/species/mouse/cta-evidence.json"] = hashlib.sha256(
        (DATA / "cta-evidence.json").read_bytes()
    ).hexdigest()
    path.write_text(json.dumps(selection, indent=2) + "\n")
    audit = {
        "taxon_id": 10090,
        "status": "unassessed",
        "policy_id": None,
        **membership_tables("mouse"),
        **native_assessment("mouse"),
    }
    audit["criteria"] = [
        {"stage": row["stage"], "status": row["status"], "rule": row["note"], "source_ids": []}
        for row in audit["funnel"]
    ]
    assert all(
        set(row) == set(columns) for key, columns in TABLE_COLUMNS.items() for row in audit[key]
    )
    name = "cta-selection-audit-mouse.json.gz"
    content = gzip.compress(
        json.dumps(audit, allow_nan=False, separators=(",", ":")).encode(), mtime=0
    )
    (package / name).write_bytes(content)
    path = package / "cta-selection-audit.json"
    receipt = json.loads(path.read_text())
    receipt["files"][name] = hashlib.sha256(content).hexdigest()
    receipt["native_updates"] = {
        "mouse": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), DATA / "cta-evidence.json", DATA / "cta-curation.json")
        }
    }
    path.write_text(json.dumps(receipt, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA / "cta-evidence.json")
    parser.add_argument(
        "--refresh-tables",
        action="store_true",
        help="Refresh packaged selection inventory and mouse audit",
    )
    args = parser.parse_args()
    if args.refresh_tables and args.output.resolve() != (DATA / "cta-evidence.json").resolve():
        parser.error("--refresh-tables requires the packaged snapshot output")
    snapshot = build(args.source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(snapshot.as_dict(), indent=2, sort_keys=True) + "\n")
    if args.refresh_tables:
        refresh_tables()
    print(snapshot.snapshot_id)


if __name__ == "__main__":
    main()
