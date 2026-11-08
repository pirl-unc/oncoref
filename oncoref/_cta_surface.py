"""Sequence-scoped extracellular features and separately measured cancer RNA/IHC."""

from __future__ import annotations

import hashlib


def interval(feature):
    """Return exact inclusive UniProt coordinates; an uncertain endpoint is not zero."""
    location = feature["location"]
    start, end = location["start"], location["end"]
    if start.get("modifier") != "EXACT" or end.get("modifier") != "EXACT":
        return None
    return [start["value"], end["value"]]


def dimensions(intervals):
    if not intervals:
        return None, None
    sizes = [end - start + 1 for start, end in intervals]
    return sum(sizes), max(sizes)


def build_surface(localization, rows, snapshot):
    entries = {r["primaryAccession"]: r for r in snapshot["uniprot_entries"]}
    sequences = snapshot["selected_sequences"]
    result, rna, ihc = [], [], []
    for loc in localization:
        if "cell_surface" not in loc["facets"]:
            continue
        row = rows[loc["entity_id"]]
        candidates = [entries[e["accession"]] for e in loc["entries"] if e["accession"] in entries]
        if len(candidates) != 1:
            raise ValueError(f"Ambiguous topology source: {loc['protein']}")
        entry = candidates[0]
        accession = entry["primaryAccession"]
        sequence = entry["sequence"]["value"]
        sha = hashlib.sha256(sequence.encode()).hexdigest()
        selected = row["representative_sequence_sha256s"]
        for key in selected:
            if hashlib.sha256(sequences[key].encode()).hexdigest() != key:
                raise ValueError("Selected sequence hash mismatch")
        exact = len(selected) == 1 and selected[0] == sha
        features = entry["features"]
        external = [
            f
            for f in features
            if f["type"] == "Topological domain" and f.get("description") == "Extracellular"
        ]
        basis = "UniProt extracellular topological domains"
        if accession in {"A6ND01", "Q6UWN0", "P38567"}:
            # GPI maturation removes the signal peptide and the C-terminal propeptide.
            external = [f for f in features if f["type"] == "Chain" and interval(f)]
            if len(external) != 1:
                raise ValueError("Expected one exact mature chain for GPI protein")
            basis = "Mature GPI-anchored chain; derived from Chain and GPI annotations"
        intervals = [interval(f) for f in external]
        known = bool(intervals) and all(i is not None for i in intervals)
        intervals = intervals if known else []
        total, longest = dimensions(intervals)
        match = "exact_sequence" if exact else "different_sequence"
        projected = intervals if exact else []
        notes = []
        if accession == "P38567":
            if len(selected) != 1 or not all(
                sequence[a - 1 : b] == sequences[selected[0]][a - 1 : b] for a, b in intervals
            ):
                raise ValueError("SPAM1 mature chain no longer matches the selected sequence")
            match = "mature_chain_identical_anchor_tail_differs"
            projected = intervals
            notes.append(
                "The selected 511-aa sequence matches the canonical 36–490 mature chain, "
                "but differs after residue 495. The 455-aa chain is a coordinate projection; "
                "GPI processing of the selected tail is not established."
            )
        if accession == "F5H094":
            notes.append(
                "The 687-aa UniProt readthrough differs from the selected 748-aa sequence. "
                "Do not transfer its domain coordinates. Human liver experiments localize "
                "LST-3TM12 mainly to ER; plasma-membrane localization was shown after "
                "heterologous expression. ER-luminal loops are not exposed on intact cells."
            )
        if accession == "Q12836":
            notes.append(
                "487 aa describes the extracellular precursor region (19–505). "
                "The exact annotated mature chain is 19–462 (444 aa); another processing "
                "annotation has an unknown endpoint. Mature zona-pellucida matrix and "
                "uncleaved membrane-tethered precursor are distinct states."
            )
        if not known:
            notes.append(
                "Transmembrane segments are annotated, but outside-facing boundaries are "
                "not. Extracellular size is unknown, not zero."
            )
        if accession == "E9PQX1":
            notes.append(
                "Released label TMEM262; current UniProt gene name CATSPERH. The sequence matches."
            )
        if accession == "Q5QGT7":
            notes.append("The only annotated external segment is the six-residue C-terminal tail.")
        if loc["protein"].startswith("OR"):
            notes.append(
                "Seven-pass olfactory receptor: total external residues span several "
                "separate loops/N-terminus, not one continuous domain."
            )
        ident = {k: loc[k] for k in ["entity_id", "protein", "gene_ids"]}
        rr = [{**r, **ident} for r in snapshot["rna"] if r["gene_id"] in loc["gene_ids"]]
        ii = [{**r, **ident} for r in snapshot["ihc"] if r["gene_id"] in loc["gene_ids"]]
        rna.extend(rr)
        ihc.extend(ii)
        tops = {}
        for cohort in ["TCGA", "validation"]:
            cohort_rows = [r for r in rr if r["cohort"] == cohort and r["samples"]]
            tops[cohort] = sorted(
                cohort_rows,
                key=lambda r: (
                    -r["expressed_samples_ptpm_ge_1"] / r["samples"],
                    -r["expressed_samples_ptpm_ge_5"] / r["samples"],
                    r["cancer"],
                ),
            )[:3]
        result.append(
            {
                **ident,
                "accession": accession,
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
                "uniprot_sequence_sha256": sha,
                "uniprot_length": len(sequence),
                "selected_sequence_sha256s": selected,
                "selected_lengths": [len(sequences[h]) for h in selected],
                "sequence_match": match,
                "basis": basis,
                "extracellular_intervals": intervals,
                "extracellular_total_aa": total,
                "extracellular_longest_aa": longest,
                "selected_extracellular_intervals": projected,
                "selected_extracellular_total_aa": dimensions(projected)[0],
                "selected_extracellular_longest_aa": dimensions(projected)[1],
                "topology_features": [
                    f
                    for f in features
                    if f["type"]
                    in {
                        "Topological domain",
                        "Transmembrane",
                        "Signal",
                        "Chain",
                        "Propeptide",
                        "Lipidation",
                    }
                ],
                "extracellular_evidence_codes": sorted(
                    {e["evidenceCode"] for f in external for e in f.get("evidences", [])}
                ),
                "transmembrane_count": sum(f["type"] == "Transmembrane" for f in features),
                "notes": notes,
                "top_rna": tops,
                "rna_tcga_cancers": len([r for r in rr if r["cohort"] == "TCGA"]),
                "ihc_measured_cancers": sum(r["measurement_status"] == "measured" for r in ii),
                "rna_ihc_scope": "Gene-level measurements; neither resolves selected proteoform or exposed epitope.",
            }
        )
    return {
        "surface_domains": result,
        "surface_cancer_rna": rna,
        "surface_cancer_ihc": ihc,
        "surface_sources": {
            k: snapshot[k]
            for k in ["uniprot_receipt", "hpa_source", "sequence_bundle_id", "reference"]
        },
    }
