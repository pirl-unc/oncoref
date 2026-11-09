"""Current human CTA representatives without downloading the sequence bundle."""

from __future__ import annotations

import json
from importlib.resources import files

import pandas as pd

from . import cta_selection, transcript_selection, translation_qualification

COLUMNS = [
    "gene_id",
    "symbol",
    "entity_id",
    "protein",
    "transcript_id",
    "protein_id",
    "length",
    "sequence_sha256",
    "selection_basis",
    "selection_policy_id",
    "reference_id",
    "annotation_source",
    "annotation_release",
    "assembly_accession",
]


def records(
    gene: str | None = None, *, annotation_release: str | int | None = None
) -> pd.DataFrame:
    """Approved CTA representatives, one row per exact gene/transcript occurrence.

    ``gene`` is an exact current symbol or unversioned Ensembl gene ID. Unknown
    genes and genes outside CTA return an empty table with the same columns.
    Transcript and protein accessions include their versions. Shared proteins
    retain each source gene's own transcript. No network or sequence download.

    The returned release is pinned. If a caller supplies ``annotation_release``,
    a mismatch raises ValueError; this function never guesses a transcript in a
    different Ensembl release. Amino-acid bytes are available in the verified
    sequence bundle and must match ``sequence_sha256``.
    """
    if gene is not None and (not isinstance(gene, str) or not gene):
        raise ValueError("gene must be a nonempty exact symbol or Ensembl gene ID")
    meta = json.loads(files("oncoref").joinpath("data/translation-qualification.json").read_text())
    reference = meta["reference"]
    if (
        annotation_release is not None
        and str(annotation_release) != reference["annotation_release"]
    ):
        raise ValueError(
            f"Representatives use Ensembl {reference['annotation_release']}; "
            f"cannot resolve against requested release {annotation_release}"
        )
    profiles = cta_selection.profiles()
    members = []
    for row in profiles[profiles.category.eq("CTA")].to_dict("records"):
        for identity in row["gene_identities"]:
            gid, symbol = identity["gene_id"], identity["symbol"]
            if gene is None or gene in {gid, symbol}:
                members.append((gid, symbol, row))
    if not members:
        return pd.DataFrame(columns=COLUMNS)
    choices = transcript_selection.records().set_index("gene_id")
    translations = translation_qualification.records().set_index("transcript_accession")
    output = []
    for gid, symbol, profile in members:
        chosen = choices.loc[gid]
        if chosen.selection_status != "selected" or not chosen.selected_transcripts:
            raise ValueError(f"Admitted CTA lacks an approved representative: {gid}")
        for tid in chosen.selected_transcripts:
            record = translations.loc[tid]
            if (
                record.gene_id != gid
                or tid not in profile["representative_transcripts"]
                or not record.target_qualified
                or record.sequence_sha256 != chosen.selected_sequence_sha256
                or record.sequence_sha256 not in profile["representative_sequence_sha256s"]
                or record.reference_id != reference["reference_id"]
            ):
                raise ValueError(f"CTA representative provenance mismatch: {gid}, {tid}")
            output.append(
                {
                    "gene_id": gid,
                    "symbol": symbol,
                    "entity_id": profile["entity_id"],
                    "protein": profile["label"],
                    "transcript_id": tid,
                    "protein_id": record.protein_id,
                    "length": int(record.length),
                    "sequence_sha256": record.sequence_sha256,
                    "selection_basis": chosen.selection_basis,
                    "selection_policy_id": chosen.selection_policy_id,
                    **{
                        key: reference[key]
                        for key in [
                            "reference_id",
                            "annotation_source",
                            "annotation_release",
                            "assembly_accession",
                        ]
                    },
                }
            )
    return (
        pd.DataFrame(output, columns=COLUMNS)
        .sort_values(["gene_id", "transcript_id"])
        .reset_index(drop=True)
    )
