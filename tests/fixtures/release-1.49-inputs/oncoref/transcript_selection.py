"""Evidence-ranked reference proteins, separate from tissue and tumor evidence."""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files

import pandas as pd

POLICY_ID = "oncoref.transcript-selection.2026-11-08.1"
PRIMARY_CONTIGS = frozenset([*(str(n) for n in range(1, 23)), "X", "Y", "MT"])
DISPOSITIONS = frozenset({"support", "hold", "resolved", "prefer"})


def evidence() -> dict:
    """Reviewed, source-linked claims and explicit decisions; no live lookups."""
    return json.loads(files("oncoref").joinpath("data/transcript-evidence.json").read_text())


def validate_evidence(ledger: dict, table: pd.DataFrame) -> list[dict]:
    """Reject stale references, guessed mappings and unreviewed overrides."""
    if ledger.get("schema") != "oncoref.transcript-evidence" or ledger.get("schema_version") != 1:
        raise ValueError("Unsupported transcript evidence schema")
    entries = ledger["findings"]
    ids = [r["id"] for r in entries]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate transcript evidence ID")
    for r in entries:
        if r.get("scope") not in {"gene", "exact_sequence"}:
            raise ValueError("Transcript finding needs gene or exact_sequence scope")
        if r["scope"] == "exact_sequence" and not r.get("sequence_sha256s"):
            raise ValueError("Exact-sequence finding requires sequence hashes")
        if r["disposition"] not in DISPOSITIONS or not r.get("rationale"):
            raise ValueError("Transcript finding requires a disposition and rationale")
        if not r.get("sources") or any(
            not s["url"].startswith("https://") or not s.get("title") for s in r["sources"]
        ):
            raise ValueError("Transcript finding requires named, linked sources")
        if not r.get("reviewed_on") or not r.get("scope") or not r.get("limitations"):
            raise ValueError("Transcript finding requires review date, scope and limitations")
        pd.to_datetime(r["reviewed_on"], format="%Y-%m-%d", errors="raise")
        locus = table.loc[table.gene_id.eq(r["gene_id"])]
        if locus.empty or set(locus.reference_id) != {r["reference_id"]}:
            raise ValueError(f"Transcript evidence reference mismatch: {r['id']}")
        for seq in r.get("sequence_sha256s", []):
            if seq not in set(locus.sequence_sha256):
                raise ValueError(f"Transcript evidence sequence absent from locus: {r['id']}")
        if r["disposition"] == "prefer":
            if r["scope"] != "exact_sequence" or len(r.get("sequence_sha256s", [])) != 1:
                raise ValueError("A preference requires one exact reference sequence")
            if r.get("evidence_kind") not in {
                "isoform_resolved_protein",
                "reviewed_annotation_reconciliation",
            }:
                raise ValueError(
                    "Gene RNA, peptide fragments and nomination cannot select an isoform"
                )
        for previous in r.get("resolves", []):
            if previous not in ids or previous == r["id"]:
                raise ValueError("Unknown or self-resolved transcript finding")
            old = next(x for x in entries if x["id"] == previous)
            if old["gene_id"] != r["gene_id"] or old["disposition"] != "hold":
                raise ValueError("A resolution must address a hold at the same locus")
            if r["disposition"] not in {"resolved", "prefer"}:
                raise ValueError("Only an explicit resolution can dispose of a hold")
            if old["scope"] == "gene" and r["scope"] != "gene":
                raise ValueError("A sequence-specific decision cannot resolve a gene-wide hold")
            if r["scope"] != "gene" and not set(old.get("sequence_sha256s", [])) <= set(
                r.get("sequence_sha256s", [])
            ):
                raise ValueError("Resolution does not cover the challenged sequence")
    return entries


def select(table: pd.DataFrame, ledger: dict | None = None) -> pd.DataFrame:
    """One decision per exact locus; never select on length or silently fall back.

    Designated but incomplete proteins remain candidates with no approved
    default. Additional evidence can hold a candidate or explicitly prefer a
    qualified alternative. Annotation support never asserts cancer production.
    """
    required = {
        "gene_id",
        "reference_id",
        "transcript_accession",
        "protein_id",
        "sequence_sha256",
        "length",
        "contig",
        "tags",
        "target_qualified",
        "qualification_reasons",
    }
    if not required <= set(table):
        raise ValueError(f"Missing transcript columns: {sorted(required - set(table))}")
    if (
        table.transcript_accession.duplicated().any()
        or not table.target_qualified.isin([True, False]).all()
    ):
        raise ValueError("Require unique transcript versions and complete qualification booleans")
    findings = validate_evidence(ledger, table) if ledger is not None else []
    by_gene = {}
    for finding in findings:
        by_gene.setdefault(finding["gene_id"], []).append(finding)
    d = table.copy()
    d["selection_rank"] = d.tags.map(
        lambda tags: 0 if "MANE_Select" in tags else 1 if "Ensembl_canonical" in tags else 2
    )
    outputs = []
    for gid, all_rows in d.groupby("gene_id", sort=True):
        if all_rows.reference_id.nunique() != 1:
            raise ValueError(f"Mixed references at locus {gid}")
        primary = all_rows.loc[all_rows.contig.isin(PRIMARY_CONTIGS)]
        designated = primary.loc[
            primary.selection_rank.eq(primary.selection_rank.min()) & primary.selection_rank.lt(2)
        ]
        basis = (
            "MANE Select"
            if not designated.empty and designated.selection_rank.min() == 0
            else "Ensembl canonical"
            if not designated.empty
            else "No designated representative"
        )
        claims = by_gene.get(gid, [])
        resolutions = {v for r in claims for v in r.get("resolves", [])}
        preferences = {
            seq for r in claims if r["disposition"] == "prefer" for seq in r["sequence_sha256s"]
        }
        reasons = []
        if len(preferences) == 1:
            designated = primary.loc[primary.sequence_sha256.isin(preferences)]
            basis = "Reviewed sequence preference"
        elif len(preferences) > 1:
            reasons.append("conflicting_reviewed_preferences")
        designated = designated.sort_values("transcript_accession")
        if designated.empty:
            reasons.append(
                "missing_primary_translation" if primary.empty else "no_designated_representative"
            )
        hashes = sorted(set(designated.sequence_sha256))
        if len(hashes) > 1:
            reasons.append("conflicting_designated_sequences")
        qualified = designated.loc[designated.target_qualified]
        passing = sorted(set(qualified.sequence_sha256))
        if not designated.empty and not passing:
            reasons.append("designated_translation_incomplete")
            reasons.extend(
                sorted({v for values in designated.qualification_reasons for v in values})
            )
        active_holds = [
            r
            for r in claims
            if r["disposition"] == "hold"
            and r["id"] not in resolutions
            and (r["scope"] == "gene" or set(r.get("sequence_sha256s", [])).intersection(hashes))
        ]
        if active_holds:
            reasons.append("unresolved_translation_evidence")
        approved = len(passing) == 1 and not reasons
        support = []
        flags = set()
        for r in designated.to_dict("records"):
            tsl = r.get("transcript_support_level", "")
            support.append(
                {
                    "transcript": r["transcript_accession"],
                    "tags": r["tags"],
                    "transcript_support_level": tsl,
                    "ccds_id": r.get("ccds_id", ""),
                    "transcript_source": r.get("transcript_source", ""),
                }
            )
            if not tsl:
                flags.add("transcript_support_not_supplied")
            elif tsl.startswith("NA"):
                flags.add("transcript_support_not_assessed")
            elif tsl[0] in "45":
                flags.add("limited_full_transcript_alignment_support")
            if "previous version" in tsl:
                flags.add("support_assigned_to_previous_version")
        if basis == "Ensembl canonical":
            flags.add("no_primary_MANE_Select")
        outputs.append(
            {
                "gene_id": gid,
                "reference_id": all_rows.reference_id.iloc[0],
                "selection_policy_id": POLICY_ID,
                "selection_basis": basis,
                "selection_status": "selected" if approved else "review",
                "selection_reasons": sorted(set(reasons)),
                "candidate_sequence_sha256s": hashes,
                "selected_sequence_sha256": passing[0] if approved else None,
                "candidate_transcripts": sorted(designated.transcript_accession),
                "selected_transcripts": sorted(qualified.transcript_accession) if approved else [],
                "candidate_protein_ids": sorted(designated.protein_id),
                "candidate_lengths": sorted(set(map(int, designated.length))),
                "qualification_passed": bool(len(passing) == 1 and len(hashes) == 1),
                "evidence_ids": sorted(r["id"] for r in claims),
                "unresolved_evidence_ids": sorted(r["id"] for r in active_holds),
                "alternative_count": int(primary.sequence_sha256.nunique() - len(hashes)),
                "candidate_support": support,
                "support_flags": sorted(flags),
                "protein_production": "not_established_by_reference_selection",
                "tumor_proteoform_expression": "not_established_by_reference_selection",
            }
        )
    return pd.DataFrame(outputs)


def records() -> pd.DataFrame:
    """Read the pinned current selection, verifying policy and evidence identities."""
    data = files("oncoref").joinpath("data")
    meta = json.loads(data.joinpath("transcript-selection.json").read_text())
    if meta["policy_id"] != POLICY_ID:
        raise ValueError("Transcript selection policy mismatch")
    for name, expected in meta["dependencies"].items():
        if hashlib.sha256(data.joinpath(name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Stale transcript selection dependency: {name}")
    path = data.joinpath("transcript-selection.csv.gz")
    if hashlib.sha256(path.read_bytes()).hexdigest() != meta["records_sha256"]:
        raise ValueError("Transcript selection checksum mismatch")
    with path.open("rb") as handle:
        frame = pd.read_csv(handle, compression="gzip", keep_default_na=False)
    for col in meta["json_columns"]:
        frame[col] = frame[col].map(json.loads)
    return frame
