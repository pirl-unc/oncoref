"""Current antigen memberships and all annotated translations, with explicit scope.

These APIs use target_selection, never the historical target_sets snapshot.
Sequences inherit gene-level evidence; membership does not establish translation
of every isoform, antigen presentation, or therapeutic safety.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from . import target_selection

SETS = {
    "CTA": "is_cta",
    "CA": "is_cancer_antigen",
    "RA": "is_reproductive_antigen",
    "CA|RA": "is_target_candidate",
    "Unsafe": None,
}
ENSEMBL_RELEASE = 115


def proteins(group="CTA") -> pd.DataFrame:
    """Current protein groups: CTA, CA, RA, CA|RA, or Unsafe.

    CA and RA include CTA. CA|RA is the union of candidate antigen proteins;
    RA-only membership does not imply tumor expression.
    """
    if group not in SETS:
        raise ValueError(f"Unknown antigen set {group!r}; choose from {tuple(SETS)}")
    rows = target_selection.results()
    mask = rows.category.eq("Unsafe") if group == "Unsafe" else rows[SETS[group]]
    return rows[mask].reset_index(drop=True)


def gene_ids(group="CTA") -> frozenset[str]:
    """All member gene IDs, preserving protein-equivalent loci."""
    return frozenset(g for ids in proteins(group).member_gene_ids for g in json.loads(ids))


def sequences(group="CTA", *, genome=None) -> pd.DataFrame:
    """All annotated protein translations, one row per transcript occurrence.

    Uses the local human Ensembl 115 reference (oncoref[genome]); does not download
    it. Missing genes/translations are recorded in ``attrs['missing_gene_ids']``.
    Stops/ambiguous residues remain in originals; peptide matching splits at them.
    """
    if genome is None:
        from pyensembl import EnsemblRelease

        genome = EnsemblRelease(ENSEMBL_RELEASE)
    if genome.release != ENSEMBL_RELEASE or genome.species.latin_name != "homo_sapiens":
        raise ValueError("Current antigen sets require human Ensembl 115")
    rows, missing = [], []
    for protein in proteins(group).itertuples():
        for gene_id in json.loads(protein.member_gene_ids):
            found = False
            try:
                transcripts = genome.gene_by_id(gene_id).transcripts
            except ValueError:
                transcripts = []
            for transcript in transcripts:
                # Decay-associated and protein_coding_LoF transcripts can have
                # annotated translations too. Do not omit those from subtraction.
                sequence = transcript.protein_sequence
                if not sequence:
                    continue
                found = True
                rows.append(
                    {
                        "entity_id": protein.entity_id,
                        "gene_id": gene_id,
                        "transcript_id": transcript.id,
                        "transcript_biotype": transcript.biotype,
                        "sequence": sequence,
                        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                        "set": group,
                    }
                )
            if not found:
                missing.append(gene_id)
    result = pd.DataFrame(
        rows,
        columns=[
            "entity_id",
            "gene_id",
            "transcript_id",
            "transcript_biotype",
            "sequence",
            "sequence_sha256",
            "set",
        ],
    )
    result.attrs.update(
        policy_id=target_selection.POLICY_ID,
        reference="homo_sapiens:GRCh38:Ensembl:115",
        missing_gene_ids=sorted(missing),
        sequence_scope="all annotated protein translations",
    )
    return result
