"""Build evidence facets without altering CTA classifications or paper coverage."""

from __future__ import annotations

from collections import Counter, defaultdict
from statistics import median

LOCATION_LABELS = {
    "secreted": "Secreted / extracellular matrix",
    "extracellular_detection": "Detected outside cells; route unresolved",
    "cell_surface": "Cell surface / integral plasma membrane",
    "cell_membrane_unresolved": "Cell membrane; exposed side unresolved",
    "intracellular_membrane": "Membrane inside the cell",
    "membrane_unspecified": "Membrane; compartment unspecified",
    "intracellular": "Intracellular / cell structure",
    "unknown": "No resolved localization",
}
EVIDENCE_LABELS = {
    "ECO:0000269": "Experimental evidence",
    "ECO:0000250": "By similarity",
    "ECO:0000255": "Sequence model / prediction",
    "ECO:0000305": "Curator inference",
    "ECO:0000303": "Author statement",
    "ECO:0000312": "Imported annotation",
}


def location_facet(location):
    """Classify the stated compartment; do not infer a location from a name."""
    name = location.get("location", {}).get("value", "").lower()
    topology = location.get("topology", {}).get("value", "").lower()
    orientation = location.get("orientation", {}).get("value", "").lower()
    if name.startswith("secreted") or name == "zona pellucida":
        return "secreted"
    if name == "cell surface":
        return "cell_surface"
    plasma = "cell membrane" in name or name in {
        "cell projection, cilium, flagellum membrane",
        "cell projection, microvillus membrane",
    }
    if plasma:
        if orientation == "cytoplasmic side":
            return "cell_membrane_unresolved"
        if "pass" in topology or "gpi-anchor" in topology or orientation == "extracellular side":
            return "cell_surface"
        return "cell_membrane_unresolved"
    if name == "membrane":
        return "membrane_unspecified"
    if "membrane" in name:
        return "intracellular_membrane"
    return "intracellular" if name else "unknown"


def evidence_labels(value):
    return sorted(
        {
            EVIDENCE_LABELS.get(e["evidenceCode"], e["evidenceCode"])
            for e in value.get("evidences", [])
        }
    ) or ["Source annotation; no evidence code"]


def build_localization(coverage, snapshot, notes):
    by_gene = defaultdict(list)
    for entry in snapshot["entries"]:
        for gid in entry["gene_ids"]:
            by_gene[gid].append(entry)
    manual = defaultdict(list)
    for note in notes:
        manual[note["protein"]].append(note)
    result = []
    for row in coverage:
        ids = set(row["gene_ids"])
        entries = {e["accession"]: e for gid in row["gene_ids"] for e in by_gene.get(gid, [])}
        annotations, entry_records = [], []
        supported, shared = defaultdict(set), set()
        for accession, entry in sorted(entries.items()):
            outside = sorted(set(entry["gene_ids"]) - ids)
            matched = sorted(set(entry["gene_ids"]) & ids)
            entry_records.append(
                {
                    "accession": accession,
                    "url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
                    "protein_name": entry["protein_name"],
                    "gene_ids": entry["gene_ids"],
                    "matched_gene_ids": matched,
                    "outside_gene_ids": outside,
                    "protein_existence": entry["protein_existence"],
                    "cautions": [c for c in entry["comments"] if c["commentType"] == "CAUTION"],
                }
            )
            for comment in entry["comments"]:
                for location in comment.get("subcellularLocations", []):
                    facet = location_facet(location)
                    annotations.append(
                        {
                            "accession": accession,
                            "gene_ids": entry["gene_ids"],
                            "outside_gene_ids": outside,
                            "mapping": "shared_with_other_loci" if outside else "within_record",
                            "isoform": comment.get("molecule"),
                            "note": comment.get("note"),
                            "facet": facet,
                            "location": location,
                            "location_evidence": evidence_labels(location.get("location", {})),
                            "topology_evidence": evidence_labels(location.get("topology", {})),
                        }
                    )
                    if outside:
                        shared.add(facet)
                    else:
                        supported[facet].update(matched)
        for note in manual[row["protein"]]:
            for facet in note["apply_facets"]:
                supported[facet].update(ids)
        facets = sorted(k for k, genes in supported.items() if genes >= ids)
        partial = sorted(k for k, genes in supported.items() if genes < ids)
        result.append(
            {
                "entity_id": row["entity_id"],
                "protein": row["protein"],
                "gene_ids": row["gene_ids"],
                "paper_gap": not row["exact_paper_sources"],
                "facets": facets or ["unknown"],
                "facet_gene_ids": {k: sorted(v) for k, v in sorted(supported.items())},
                "partial_facets": partial,
                "shared_facets": sorted(shared),
                "entries": entry_records,
                "annotations": annotations,
                "paper_notes": manual[row["protein"]],
                "sequence_localization_verified": False,
            }
        )
    return result


def build_followup(coverage, reviews, rows, notes, snapshot):
    reviews = {r["protein"]: r for r in reviews}
    coverage_map = {r["protein"]: r for r in coverage}
    groups = {g["id"]: g for g in notes["concern_groups"]}
    concerns = []
    for note in notes["concerns"]:
        review = reviews[note["protein"]]
        group = groups[note["concern_group"]]
        if review["outcome"] != group["review_outcome"]:
            raise ValueError(f"Concern outcome mismatch: {note['protein']}")
        concerns.append(
            {
                **note,
                "category": group["label"],
                "outcome": review["outcome"],
                "entity_id": review["entity_id"],
                "gene_ids": review["gene_ids"],
                "citations": review["citations"],
                "paper_gap": not coverage_map[note["protein"]]["exact_paper_sources"],
                "classification_effect": "Review finding; frozen CTA assignment unchanged.",
            }
        )
    expected = {
        p for p, r in reviews.items() if r["outcome"] in {"normal_review", "identity_review"}
    }
    if {r["protein"] for r in concerns} != expected:
        raise ValueError("The concern summary must cover every flagged individual review")
    gap_groups = {p: g for g in notes["gap_groups"] for p in g["proteins"]}
    gap_names = {r["protein"] for r in coverage if not r["exact_paper_sources"]}
    if set(gap_groups) != gap_names or sum(len(g["proteins"]) for g in notes["gap_groups"]) != len(
        gap_names
    ):
        raise ValueError("Each paper gap needs exactly one descriptive biological group")
    normal_questions = {"CALHM3", "PRAMEF18", "TAF11L11", "TAF11L13", "TAF11L14"}
    gaps = []
    for row in coverage:
        if row["protein"] not in gap_names:
            continue
        protein = row["protein"]
        measurement = rows[row["entity_id"]]
        review = reviews.get(protein)
        group = gap_groups[protein]
        recommendation = (
            "Resolve the normal-cell question before presenting this as a cleared target. The current finding does not prove exact human protein expression."
            if protein in normal_questions
            else "Retain as a discovery CTA under the current admission rules; prioritize human protein and tumor validation."
        )
        if protein in {"ATF7-NPFF", "ENSG00000273155"}:
            recommendation += " Resolve the exact product before transferring parent-gene or shared-entry experiments."
        gaps.append(
            {
                **row,
                "biological_group": group["label"],
                "group_rationale": group["rationale"],
                "evidence_summary": review["summary"]
                if review
                else " ".join(f["finding"] for f in measurement["literature_review"]),
                "identity_note": review["identity"]
                if review
                else "Related-species experiment; exact human localization remains unresolved.",
                "normal_question": protein in normal_questions,
                "exact_product_question": protein in {"ATF7-NPFF", "ENSG00000273155"},
                "recommendation": recommendation,
                "translation_qualified": measurement["translation_qualified"],
                "translation_basis": measurement["translation_support_reasons"],
                "tumor_max_tpm": measurement["tumor_max_tpm"],
                "normal_max_ntpm": measurement["normal_max_ntpm"],
                "normal_max_tpm": measurement["normal_max_tpm"],
                "developmental_lineage_max_tpm": measurement["developmental_lineage_max_tpm"],
                "reproductive_max_ntpm": measurement["reproductive_max_ntpm"],
            }
        )
    comparison = {}
    for name, is_gap in [("paper_gap", True), ("paper_covered", False)]:
        subset = [
            rows[r["entity_id"]] for r in coverage if (not r["exact_paper_sources"]) == is_gap
        ]
        values = [r["tumor_max_tpm"] for r in subset if r["tumor_max_tpm"] is not None]
        comparison[name] = {
            "proteins": len(subset),
            "routes": dict(sorted(Counter(r["route"] for r in subset).items())),
            "measured_tumor_tpm": len(values),
            "median_highest_tumor_tpm": median(values),
            "highest_tumor_below_10_tpm": sum(v < 10 for v in values),
            "translation_qualified": sum(r["translation_qualified"] for r in subset),
        }
    localization = build_localization(coverage, snapshot, notes["localization_notes"])
    return {
        "gap_assessments": gaps,
        "gap_groups": notes["gap_groups"],
        "gap_comparison": comparison,
        "concerns": concerns,
        "concern_groups": [
            {**g, "proteins": [r["protein"] for r in concerns if r["concern_group"] == g["id"]]}
            for g in notes["concern_groups"]
        ],
        "localization": localization,
        "localization_labels": LOCATION_LABELS,
        "localization_source": snapshot["source"],
        "localization_counts": {
            scope: dict(
                sorted(
                    Counter(
                        f
                        for r in localization
                        if scope == "all_ctas" or r["paper_gap"]
                        for f in r["facets"]
                    ).items()
                )
            )
            for scope in ["all_ctas", "paper_gaps"]
        },
        "localization_rules": [
            "Locations are gene-mapped annotations, not validation of the selected OncoRef sequence or tumor-cell accessibility.",
            "A summary facet requires mapped evidence for every member of the protein record; partial and shared-entry annotations remain separate.",
            "Cell-surface facets require an explicit surface location or an integral/GPI-anchored plasma-membrane annotation. Peripheral association without orientation remains unresolved.",
            "Secretory vesicles and organelle membranes remain intracellular. An unspecified membrane is not assumed to be the cell surface.",
            "UniProt location evidence and topology evidence are retained separately. Protein-existence level does not establish experimental localization.",
            "Primary-paper additions retain their species, cell model and proteoform scope. They do not change classification or paper-coverage counts.",
        ],
    }
