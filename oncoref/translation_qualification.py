"""Translation quality is an admission gate, independent of normal-tissue hazard."""

from __future__ import annotations

import hashlib
import json
import re
from importlib.resources import files

import pandas as pd

POLICY_ID = "oncoref.translation-qualification.2026-11-04.1"
PRIMARY_CONTIGS = frozenset([*(str(n) for n in range(1, 23)), "X", "Y", "MT"])
CODING_BIOTYPES = frozenset(
    [
        "protein_coding",
        "IG_V_gene",
        "IG_C_gene",
        "IG_J_gene",
        "IG_D_gene",
        "TR_V_gene",
        "TR_C_gene",
        "TR_J_gene",
        "TR_D_gene",
    ]
)


def valid_segments(sequence: str) -> tuple[dict, ...]:
    """AA20 runs with original zero-based half-open coordinates; never bridge X/*/U."""
    if not isinstance(sequence, str) or not sequence.isascii():
        raise ValueError("Require an ASCII source sequence")
    return tuple(
        {
            "start": m.start(),
            "end": m.end(),
            "sequence": m.group(),
            "sha256": hashlib.sha256(m.group().encode()).hexdigest(),
        }
        for m in re.finditer(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence)
    )


def qualify(sequence: str, annotation: dict) -> dict:
    """Require the exact primary-locus annotated CDS, complete and without LoF/NMD."""
    reasons = []
    if annotation.get("contig") not in PRIMARY_CONTIGS:
        reasons.append("alternate_or_unplaced_locus")
    if annotation.get("transcript_biotype") not in CODING_BIOTYPES:
        reasons.append("nonqualifying_transcript_biotype")
    if annotation.get("gene_biotype") == "protein_coding_LoF":
        reasons.append("gene_loss_of_function_annotation")
    if not annotation.get("cds_bases", 0):
        reasons.append("missing_annotated_cds")
    if annotation.get("start_codon_bases") != 3:
        reasons.append("incomplete_start")
    if annotation.get("stop_codon_bases") != 3:
        reasons.append("incomplete_stop")
    if set(annotation.get("tags", [])) & {
        "cds_start_NF",
        "cds_end_NF",
        "mRNA_start_NF",
        "mRNA_end_NF",
    }:
        reasons.append("incomplete_annotation_tag")
    # FASTA terminal stop is optional; internal stops and non-AA20 are retained.
    residues = sequence.removesuffix("*")
    if not residues or re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", residues) is None:
        reasons.append("ambiguous_or_interrupted_sequence")
    if not residues.startswith("M"):
        reasons.append("noncanonical_start_residue")
    return {
        "target_qualified": not reasons,
        "qualification_reasons": reasons,
        "qualification_policy_id": POLICY_ID,
    }


def records() -> pd.DataFrame:
    """Offline source-bound transcript/protein quality table; returns a fresh frame."""
    data = files("oncoref").joinpath("data")
    meta = json.loads(data.joinpath("translation-qualification.json").read_text())
    path = data.joinpath("translation-qualification.csv.gz")
    if hashlib.sha256(path.read_bytes()).hexdigest() != meta["records_sha256"]:
        raise ValueError("Translation qualification checksum mismatch")
    with path.open("rb") as handle:
        frame = pd.read_csv(handle, compression="gzip", keep_default_na=False)
    for column in ["tags", "qualification_reasons"]:
        frame[column] = frame[column].map(json.loads)
    return frame


def annotate_legacy_longest(facts: pd.DataFrame, table: pd.DataFrame | None = None) -> pd.DataFrame:
    """Revalidate common longest primary-locus sequence under the pinned reference.

    An alternate-locus alias or a shorter complete isoform cannot silently rescue
    the selected longest protein. A separate ncORF nomination cannot establish
    translation of an NMD reference peptide. Missing mappings require review.
    """
    table = records() if table is None else table.copy()
    primary = table[table.contig.isin(PRIMARY_CONTIGS)]
    longest = primary[primary.length.eq(primary.groupby("gene_id").length.transform("max"))]
    by_gene = dict(iter(longest.groupby("gene_id")))
    qualified, reasons, hashes = [], [], []
    for ids in facts.member_gene_ids.map(json.loads):
        selections = [by_gene.get(g) for g in ids]
        if any(t is None or t.empty for t in selections):
            qualified.append(False)
            reasons.append("missing_primary_translation")
            hashes.append("[]")
            continue
        common = set.intersection(*(set(t.sequence_sha256) for t in selections))
        passing = set.intersection(
            *(set(t.loc[t.target_qualified, "sequence_sha256"]) for t in selections)
        )
        qualified.append(bool(passing))
        hashes.append(json.dumps(sorted(common)))
        failures = sorted({r for t in selections for rr in t.qualification_reasons for r in rr})
        reasons.append(
            "" if passing else " | ".join(failures or ["representative_sequence_mismatch"])
        )
    result = facts.copy()
    result["translation_qualified"] = qualified
    result["translation_qualification_reasons"] = reasons
    result["representative_sequence_sha256s"] = hashes
    result["translation_qualification_policy"] = POLICY_ID
    return result


def annotate(
    facts: pd.DataFrame, table: pd.DataFrame | None = None, *, ledger: dict | None = None
) -> pd.DataFrame:
    """Attach current representative decisions without conflating support and quality."""
    from . import transcript_selection

    selected = (
        transcript_selection.records()
        if table is None
        else transcript_selection.select(
            table, transcript_selection.evidence() if ledger is None else ledger
        )
    )
    by_gene = selected.set_index("gene_id").to_dict("index")
    rows = []
    for ids in facts.member_gene_ids.map(json.loads):
        candidates = [by_gene.get(g) for g in ids]
        if not candidates or any(r is None for r in candidates):
            rows.append(
                {
                    "translation_qualified": False,
                    "translation_qualification_reasons": "missing_primary_translation",
                    "representative_selection_approved": False,
                    "representative_selection_reasons": "missing_primary_translation",
                    "representative_sequence_sha256s": "[]",
                    "representative_candidate_sha256s": "[]",
                    "representative_evidence_ids": "[]",
                    "representative_transcripts": "[]",
                }
            )
            continue
        common = set.intersection(*(set(r["candidate_sequence_sha256s"]) for r in candidates))
        qualified = bool(common) and all(r["qualification_passed"] for r in candidates)
        approved_hashes = {r["selected_sequence_sha256"] for r in candidates}
        approved = (
            all(r["selection_status"] == "selected" for r in candidates)
            and len(approved_hashes) == 1
        )
        reasons = sorted({s for r in candidates for s in r["selection_reasons"]})
        if not common:
            reasons.append("representative_sequence_mismatch")
        rows.append(
            {
                "translation_qualified": qualified,
                "translation_qualification_reasons": ""
                if qualified
                else " | ".join(reasons or ["designated_translation_incomplete"]),
                "representative_selection_approved": approved,
                "representative_selection_reasons": " | ".join(sorted(set(reasons))),
                "representative_sequence_sha256s": json.dumps(sorted(approved_hashes))
                if approved
                else "[]",
                "representative_candidate_sha256s": json.dumps(sorted(common)),
                "representative_evidence_ids": json.dumps(
                    sorted({v for r in candidates for v in r["evidence_ids"]})
                ),
                "representative_transcripts": json.dumps(
                    sorted({v for r in candidates for v in r["selected_transcripts"]})
                ),
            }
        )
    result = facts.copy()
    for column in rows[0] if rows else []:
        result[column] = [r[column] for r in rows]
    result["translation_qualification_policy"] = POLICY_ID
    result["representative_selection_policy"] = transcript_selection.POLICY_ID
    return result
