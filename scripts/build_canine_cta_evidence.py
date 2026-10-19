#!/usr/bin/env python3
"""Rebuild the curated canine CTA snapshot from pinned primary-source bytes.

No network access. See docs/species-cta-evidence.md for inputs and curation scope.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement
from oncoref.reference_identity import ReferenceIdentity, _identity
from oncoref.species_antigens import CTAEvidenceSnapshot
from oncoref.species_assets import registered_references

HINDRIKS = {
    "id": "hindriks2025",
    "version": "2025-05-06-publisher-version",
    "url": "https://research-portal.uu.nl/ws/portalfiles/portal/264124617/1-s2.0-S0165242725000662-main.pdf",
    "doi": "10.1016/j.vetimm.2025.110946",
    "pmid": "40359694",
    "sha256": "825450315ad1433b7c98d49a2e8c715f3e3760ca2889c5bb1a1041c2bfcdc391",
    "scope": "Publisher PDF: Tables 1-3 and sections 2.1-3.4. Supplementary files and plotted values are not digitized.",
    "study_id": "PMID:40359694",
    "lineage": ["PMID:40359694"],
}
NEMEC = {
    "id": "nemec2019",
    "version": "2019-article-abstract-europepmc-2026-10-13",
    "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30854786%20AND%20SRC:MED&resultType=core&format=json",
    "doi": "10.1111/vco.12475",
    "pmid": "30854786",
    "sha256": "6f27557f46c3607dde7bb32c158e7cc0a26adc6a75a230dbdf9dfc748bd11b73",
    "scope": "Article abstract indexed by Europe PMC only. Full-text peptide sequences, assay units and sample-level values are unavailable in this import.",
    "study_id": "PMID:30854786",
    "lineage": ["PMID:30854786"],
}

# Table 1 IDs, retaining LOC symbols in the target annotation. No human resolver.
GENE_IDS = {
    "CAGE1": "609230",
    "CTCFL": "477271",
    "DDX53": "119868264",
    "NT5C1B": "100856179",
    "PAGE3-like": "119868444",
    "PRAME": "477552",
    "SSX6": "612524",
    "MAGEA13P": "492173",
}
PCR_TARGETS = ["CAGE1", "CTCFL", "DDX53", "MAGE", "NT5C1B", "PAGE3-like", "PRAME", "SSX6"]

# Table 2 as printed (including inconsistent sign/value pairs). Empty cells mean
# no amplification measured, NOT zero, and '<' retains an upper bound.
# Columns follow PCR_TARGETS. Cell line spelling is deliberately source-native.
PCR_ROWS = [
    ("melanoma", "COMM.2", ["", "", "", "+:0.15", "", "-:0.05", "", ""]),
    ("melanoma", "COMM.O", ["-:0.01", "", "", "++:5.3", "", "+++:25.0", "-:0.04", "+:0.6"]),
    ("melanoma", "COMM.K", ["", "", "", "++:2.6", "", "+:0.58", "-:0.01", ""]),
    ("melanoma", "LMCK", ["-:0.08", "", "", "-:0.9", "", "", "", ""]),
    ("melanoma", "CMM10", ["-:0.04", "", "", "+:0.11", "", "", "", "-:0.04"]),
    ("melanoma", "CMM12", ["-:0.05", "-:0.05", "", "+:0.10", "", "", "", "-:0.07"]),
    ("melanoma", "UCKD9M3", ["-:0.02", "", "", "+++:24.9", "", "-:<0.01", "", ""]),
    ("melanoma", "UCKD9M4", ["+:0.19", "", "", "+:0.14", "", "+++:74.2", "++:20.2", ""]),
    ("metastasis", "UCKD9M2", ["+:0.25", "", "", "+:0.13", "", "+++:97.8", "+++:31.2", ""]),
    ("metastasis", "UCKD9M5", ["-:<0.01", "", "", "++:1.4", "", "", "", ""]),
    ("metastasis", "OLGA", ["-:0.01", "", "", "+++:13.2", "", "-:<0.01", "++:3.7", "-:0.05"]),
    ("metastasis", "NM.G", ["-:<0.01", "", "", "+:0.21", "", "+:0.24", "", ""]),
    ("other", "GL-1", ["", "+:0.11", "", "", "", "", "", ""]),
    ("other", "CLBL-1", ["", "", "", "", "", "", "", ""]),
    ("other", "3132", ["-:0.03", "", "", "++:3.9", "-:<0.01", "-:0.04", "", "+:0.5"]),
    ("other", "CL-1", ["-:<0.01", "-:0.02", "", "++:3.0", "-:<0.01", "", "", ""]),
    ("healthy", "Placenta", ["", "", "", "", "", "+++:15.2", "+:0.49", ""]),
    ("healthy", "Kidney", ["", "++:2.2", "++:4.9", "", "-:<0.01", "", "", ""]),
    ("healthy", "MSCEK", ["", "", "", "", "", "", "", ""]),
    ("healthy", "CPEK", ["-:<0.01", "", "", "", "-:0.01", "", "", ""]),
    ("healthy", "RPE", ["", "", "", "", "", "", "+:0.25", ""]),
]


def verify(path, digest):
    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
        raise ValueError(f"Source checksum mismatch: {path}")


def reference_targets(path):
    """Preserve every annotation sequence pair, including one-to-many mappings."""
    genes = {}
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gid = row["GeneID"]
            if not gid:
                continue
            target = genes.setdefault(
                gid, {"id": gid, "symbol": row["symbol"] or gid, "sequences": []}
            )
            transcript = protein = None
            if row["# feature"] == "mRNA":
                transcript, protein = row["product_accession"], row["related_accession"]
            elif row["# feature"] == "CDS":
                protein, transcript = row["product_accession"], row["related_accession"]
            if transcript or protein:
                sequence = {"transcript_id": transcript or None, "protein_id": protein or None}
                if sequence not in target["sequences"]:
                    target["sequences"].append(sequence)
    for target in genes.values():
        target["sequences"].sort(key=lambda s: (s["transcript_id"] or "", s["protein_id"] or ""))
    return genes


def build(feature_table, hindriks_pdf, nemec_record):
    registration = next(r for r in registered_references("dog") if r["id"] == "refseq106")
    reference = ReferenceIdentity.from_dict(registration["identity"])
    verify(feature_table, reference.content_provenance["feature_table"])
    verify(hindriks_pdf, HINDRIKS["sha256"])
    verify(nemec_record, NEMEC["sha256"])
    targets = reference_targets(feature_table)
    payload = {
        "schema": "oncoref.cta-evidence",
        "schema_version": 1,
        "taxon_id": 9615,
        "source_version": "2026-10-13.1",
        "references": [{"id": "refseq106", "identity": reference.as_dict()}],
        "sources": [HINDRIKS, NEMEC],
        "contracts": [],
        "candidates": [],
        "evidence": [],
        "coverage": {
            "status": "partial",
            "limitations": [
                "Published nominations and assay follow-ups only; no expression-derived panel or patient admission policy is supplied.",
                "Only RefSeq annotation 106 mappings are curated; no cross-annotation equivalence is inferred.",
                "Nemec 2019 is abstract-only: peptide sequences, per-gene tissue details, native numeric units and sample/donor counts are not imported.",
                "Hindriks 2025 Table 2 footnote thresholds conflict with section 3.1; LMCK MAGE (-, 0.9%) and UCKD9M4 PRAME (++, 20.2%) signs conflict with the stated thresholds. Values and labels remain as reported.",
                "Blank PCR cells mean no amplification measured; unreported tissues and donor counts remain unknown, never zero.",
                "Hindriks section 2.1 describes healthy tissues from four dogs or client-owned dogs; assay-specific donor identities and cross-assay overlap are not given.",
                "MAGE family assays do not identify individual gene abundance. Named MAGEA13P/MAGEB10 follow-ups do not enumerate all primer or antibody targets.",
                "Study lineage groups repeated assays; unknown shared donors or cell-line ancestry must not be counted as independent biological replication.",
            ],
        },
    }
    candidates = {}
    sources = {s["id"]: s for s in payload["sources"]}

    def candidate(
        label, source, *, kind="published_nomination", locator="Abstract", ids=None, ambiguous=False
    ):
        if label in candidates:
            candidates[label]["derivation"].append(
                {
                    "kind": kind,
                    "source_id": source,
                    "locator": locator,
                    "note": "Additional published nomination; orthology alone does not establish canine tissue restriction.",
                }
            )
            return
        if ids is None:
            if label in GENE_IDS:
                ids = [GENE_IDS[label]]
                method = "Table 1 NCBI Gene ID exact match to pinned feature table"
            else:
                ids = [gid for gid, t in targets.items() if t["symbol"] == label]
                method = "Source symbol exact match to pinned canine feature table; assay sequence resolution remains unknown"
        else:
            method = "Reported MAGE subtype follow-ups (section 3.2), not exhaustive family assay mapping"
        mapped = [targets[gid] for gid in sorted(ids) if gid in targets]
        mapping = {
            "reference_id": reference.reference_id,
            "gene_namespace": "NCBI Gene",
            "status": "unresolved"
            if not mapped
            else "ambiguous"
            if ambiguous or len(mapped) > 1
            else "unique",
            "method": method,
            "complete": bool(mapped) and not ambiguous and len(mapped) == len(ids),
            "targets": mapped,
            "notes": "Sequence pairs are annotation associations, not measured isoforms. All supported alternatives are retained; family assay membership is incomplete."
            if ambiguous
            else "Annotation sequence associations do not allocate measured abundance to transcripts or proteins.",
        }
        candidates[label] = {
            "id": label,
            "label": label,
            "nomination_kind": "published_nomination"
            if kind != "assay_followup"
            else "assay_followup",
            "derivation": [
                {
                    "kind": kind,
                    "source_id": source,
                    "locator": locator,
                    "note": "Published canine nomination or assay follow-up; no human CTA policy is transferred.",
                }
            ],
            "evidence_status": "conflicting"
            if label
            in {
                "AKR1E2",
                "SPECC1",
                "TPX2",
                "IGF2BP3",
                "NT5C1B",
                "CTCFL",
                "DDX53",
                "PAGE3-like",
                "PRAME",
                "MAGE",
            }
            else "limited",
            "mappings": [mapping],
            "evidence_ids": [],
            "notes": "Evidence is limited to the imported assay and tissue scope; nomination is not proof of tumor restriction or eligibility.",
        }

    for label in ("IGF2BP3", "AKR1E2", "SPECC1", "TPX2", "NT5C1B"):
        candidate(label, "nemec2019")
    for label in PCR_TARGETS:
        candidate(
            label,
            "hindriks2025",
            kind="orthology_nomination",
            locator="Tables 1-2; section 3.1",
            ids=["492173", "403710"] if label == "MAGE" else None,
            ambiguous=label == "MAGE",
        )
    candidate("MAGEA13P", "hindriks2025", kind="assay_followup", locator="Table 1; section 3.2")
    candidate("MAGEB10", "hindriks2025", kind="assay_followup", locator="Section 3.2")

    def evidence(
        labels,
        source,
        locator,
        assay,
        tissue,
        role,
        context,
        feature,
        unit,
        reported,
        *,
        value=None,
        status="qualitative",
        samples=None,
        n=None,
        donors=None,
        repeats=None,
        positives=None,
        notes="Sample/donor coverage beyond the reported assay is unknown.",
    ):
        src = sources[source]
        contract = ExpressionContract(
            9615,
            source,
            src["version"],
            context,
            None,
            "unresolved",
            "study assay target",
            ExpressionMeasurement(feature, unit, "source_native", False),
            {"source": src["sha256"]},
            "Historical assay annotation is unresolved. Query annotation mappings are separate and do not establish measurement/reference equivalence.",
        )
        if not any(c["contract_id"] == contract.contract_id for c in payload["contracts"]):
            payload["contracts"].append(contract.as_dict())
        eid = f"{source}:{len(payload['evidence']) + 1:03d}"
        row = {
            "id": eid,
            "candidate_ids": labels,
            "source_id": source,
            "study_id": src["study_id"],
            "contract_id": contract.contract_id,
            "locator": locator,
            "assay": assay,
            "tissue": tissue,
            "tissue_role": role,
            "sample_ids": samples or [],
            "sample_count": n,
            "donor_count": donors,
            "technical_replicates": repeats,
            "source_fields": {},
            "observation": {
                "status": status,
                "value": value,
                "reported": reported,
                "positive_count": positives,
            },
            "notes": notes,
        }
        payload["evidence"].append(row)
        for label in labels:
            candidates[label]["evidence_ids"].append(eid)

    evidence(
        ["IGF2BP3"],
        "nemec2019",
        "Abstract",
        "quantitative/end-point RT-PCR",
        "peripheral blood mononuclear cells",
        "somatic",
        "normal",
        "gene",
        None,
        "High PBMC expression reported.",
    )
    evidence(
        ["AKR1E2", "SPECC1", "TPX2"],
        "nemec2019",
        "Abstract",
        "quantitative/end-point RT-PCR",
        None,
        "somatic",
        "normal",
        "gene",
        None,
        "The three-gene group showed high expression in critical normal tissues, including kidney, brain and marrow.",
        notes="This is a group-level abstract statement: which individual gene was high in which tissue is not resolved. Do not expand into gene-by-tissue values.",
    )
    evidence(
        ["AKR1E2", "SPECC1", "TPX2"],
        "nemec2019",
        "Abstract",
        "quantitative/end-point RT-PCR",
        "histiocytic sarcoma and T-cell lymphoma biopsies",
        "tumor",
        "tumor",
        "gene",
        None,
        "Variable tumor expression reported for the three-gene group.",
    )
    evidence(
        ["NT5C1B"],
        "nemec2019",
        "Abstract",
        "quantitative/end-point RT-PCR",
        "T-cell lymphoma",
        "tumor",
        "tumor",
        "gene",
        None,
        "Detected in T-cell lymphomas.",
    )
    evidence(
        ["NT5C1B"],
        "nemec2019",
        "Abstract",
        "quantitative/end-point RT-PCR",
        None,
        "somatic",
        "normal",
        "gene",
        None,
        "Low expression in some normal dog tissues; tissue identities not given in the abstract.",
    )

    for group, sample, cells in PCR_ROWS:
        tissue = {
            "Placenta": "placenta",
            "Kidney": "kidney",
            "RPE": "retinal pigment epithelium",
            "MSCEK": "epidermal keratinocytes",
            "CPEK": "epidermal keratinocytes",
        }.get(sample, group)
        role = (
            "tumor" if group != "healthy" else "reproductive" if sample == "Placenta" else "somatic"
        )
        for label, cell in zip(PCR_TARGETS, cells):
            value = None
            status = "no_amplification"
            if cell:
                number = cell.split(":")[1]
                status = "below_limit" if number.startswith("<") else "reported"
                value = float(number.lstrip("<"))
            evidence(
                [label],
                "hindriks2025",
                f"Table 2: {sample}, {label}",
                "RT-PCR relative to RPS5 and testis",
                tissue,
                role,
                "normal" if sample in {"Placenta", "Kidney"} else "cell_line",
                "unknown" if label == "MAGE" else "gene",
                "percent of testis-normalized expression",
                cell.replace(":", " ") if cell else "blank: no amplification measured",
                value=value,
                status=status,
                samples=[sample],
                n=None if sample in {"Placenta", "Kidney"} else 1,
                repeats=3,
                notes="Technical triplicates are not donors. Table signs and numeric values are preserved despite conflicting threshold descriptions; no positivity cutoff is recomputed. MAGE PCR is a multi-gene assay.",
            )
    for label in PCR_TARGETS:
        evidence(
            [label],
            "hindriks2025",
            "Sections 2.5, 3.1; Figure 1",
            "RT-PCR relative to RPS5 and testis",
            "testis",
            "reproductive",
            "normal",
            "unknown" if label == "MAGE" else "gene",
            "percent of testis-normalized expression",
            "Testis calibrator; all eight assays amplified.",
            value=100,
            status="reported",
            repeats=3,
            notes="100% is assigned by calibration, not an absolute expression measurement. Individual Ct values and donor counts are not tabulated.",
        )
    evidence(
        ["MAGEA13P"],
        "hindriks2025",
        "Section 3.2",
        "MAGEA13P-specific RT-PCR",
        "oral melanoma",
        "tumor",
        "cell_line",
        "gene",
        None,
        "Similar MAGE and MAGEA13P Ct values in CMM10 supported subtype attribution.",
        samples=["CMM10"],
        n=1,
    )
    evidence(
        ["MAGEA13P"],
        "hindriks2025",
        "Section 3.2",
        "MAGEA13P-specific RT-PCR",
        "testis",
        "reproductive",
        "normal",
        "gene",
        "percent of MAGE expression",
        "Estimated MAGEA13P contribution using CMM10 calibration.",
        value=5.2,
        status="reported",
        notes="Subtype contribution estimate; not a transcript or protein abundance measurement.",
    )
    evidence(
        ["MAGEB10"],
        "hindriks2025",
        "Section 3.2",
        "Sanger sequencing of family RT-PCR product",
        "testis",
        "reproductive",
        "normal",
        "unknown",
        None,
        "PCR product resembled MAGEB10; a qualitative family-product attribution.",
        notes="No exact product sequence or accession is provided here; sequence equivalence and isoform abundance remain unresolved.",
    )

    for tissue, positive, total in [
        ("cutaneous melanoma", 0, 12),
        ("uveal melanoma", 0, 2),
        ("digital melanoma", 2, 8),
        ("oral melanoma", 2, 11),
        ("cutaneous melanocytoma", 0, 5),
        ("uveal melanocytoma", 0, 4),
        ("oral melanocytoma", 0, 1),
    ]:
        evidence(
            ["MAGE"],
            "hindriks2025",
            "Table 3",
            "IHC anti-MAGE clone 6C1",
            tissue,
            "tumor",
            "tumor",
            "protein_group",
            "positive samples / examined samples",
            f"{positive}/{total} positive",
            status="reported",
            positives=positive,
            n=total,
            notes="Family-reactive antibody; subtype/isoform unresolved. Disjoint anatomical rows reproduce 4/33 melanomas and 0/10 melanocytomas without duplicating totals. Donor uniqueness is not asserted.",
        )
    for tissue, count in [
        (t, 3)
        for t in [
            "cerebellum",
            "cerebrum",
            "duodenum",
            "heart",
            "kidney",
            "liver",
            "lung",
            "lymph node",
            "nail bed",
            "buccal oral mucosa",
            "pancreas",
            "skin",
            "spleen",
            "thyroid",
            "uterus",
        ]
    ] + [("ovary", 2), ("testis", 1), ("eye", 1)]:
        evidence(
            ["MAGE"],
            "hindriks2025",
            "Sections 3.3-3.4",
            "IHC anti-MAGE clone 6C1",
            tissue,
            "reproductive" if tissue in {"uterus", "ovary", "testis"} else "somatic",
            "normal",
            "protein_group",
            "positive samples / examined samples",
            "Positive spermatogonia and primary spermatocytes; Sertoli and Leydig cells negative."
            if tissue == "testis"
            else "No immunopositivity reported.",
            status="reported",
            positives=1 if tissue == "testis" else 0,
            n=count,
            notes="Assay-negative in the sampled tissue is not zero gene abundance or a complete normal-tissue survey. Family specificity and donor overlap remain unresolved.",
        )

    payload["candidates"] = sorted(candidates.values(), key=lambda c: c["id"])
    payload["contracts"].sort(key=lambda c: c["contract_id"])
    return CTAEvidenceSnapshot.from_dict(
        {**payload, "snapshot_id": _identity("cta-evidence", 1, payload)}
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-table", type=Path, required=True)
    parser.add_argument("--hindriks-pdf", type=Path, required=True)
    parser.add_argument("--nemec-record", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    snapshot = build(args.feature_table, args.hindriks_pdf, args.nemec_record)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(snapshot.as_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(snapshot.snapshot_id)


if __name__ == "__main__":
    main()
