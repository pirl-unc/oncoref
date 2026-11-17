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
    "Other": None,
    "Unsafe|Other": None,
}
EXCLUSION_CATEGORIES = {
    "Unsafe": ("Unsafe",),
    "Other": ("Other",),
    "Unsafe|Other": ("Unsafe", "Other"),
}
ENSEMBL_RELEASE = 115


def proteins(
    group: str = "CTA",
    *,
    profile: str = "all",
    lactating: bool = False,
    pregnant: bool = False,
    allow_tissues: tuple[str, ...] = (),
    protect_tissues: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Current antigen groups or Unsafe, Other, and Unsafe|Other backgrounds.

    CA and RA include CTA. CA|RA is the union of candidate antigen proteins;
    RA-only membership does not imply tumor expression.
    Other denotes unresolved or nonqualifying proteins, not established hazards.
    Optional ``profile`` ("all", "core", "male-patient", "female-patient") and
    tissue-tolerance overrides narrow positive sets by host compartment rules.
    """
    if group not in SETS:
        raise ValueError(f"Unknown antigen set {group!r}; choose from {tuple(SETS)}")
    has_host_filter = (
        profile != "all" or lactating or pregnant or bool(allow_tissues) or bool(protect_tissues)
    )
    if group in EXCLUSION_CATEGORIES:
        if has_host_filter:
            raise ValueError(
                "Host profiles apply to candidate antigen sets, not Unsafe/Other backgrounds"
            )
        rows = target_selection.results()
        return rows[rows.category.isin(EXCLUSION_CATEGORIES[group])].reset_index(drop=True)
    rows = target_selection.results()
    if not has_host_filter:
        return rows[rows[SETS[group]]].reset_index(drop=True)
    from . import target_filters

    filt = target_filters.TargetFilter(
        host_profile=profile,
        lactating=lactating,
        pregnant=pregnant,
        allow_tissues=tuple(allow_tissues),
        protect_tissues=tuple(protect_tissues),
    )
    evaluated = target_filters.apply(filt, rows)
    ca, ra = evaluated.filtered_cancer_antigen, evaluated.filtered_reproductive_antigen
    member = {"CTA": ca & ra, "CA": ca, "RA": ra, "CA|RA": ca | ra}[group]
    selected = evaluated[evaluated.keep & member].copy()
    selected["category"] = selected.filtered_category
    return selected.reset_index(drop=True)


def gene_ids(
    group: str = "CTA",
    *,
    profile: str = "all",
    lactating: bool = False,
    pregnant: bool = False,
    allow_tissues: tuple[str, ...] = (),
    protect_tissues: tuple[str, ...] = (),
) -> frozenset[str]:
    """All member gene IDs, preserving protein-equivalent loci."""
    has_host_filter = (
        profile != "all" or lactating or pregnant or bool(allow_tissues) or bool(protect_tissues)
    )
    frame = (
        proteins(
            group,
            profile=profile,
            lactating=lactating,
            pregnant=pregnant,
            allow_tissues=allow_tissues,
            protect_tissues=protect_tissues,
        )
        if has_host_filter
        else proteins(group)
    )
    return frozenset(g for ids in frame.member_gene_ids for g in json.loads(ids))


def loss_consequences(gene: str) -> list[dict]:
    """Reviewed tissue-loss consequences beside observed reproductive contexts.

    Includes measured native RNA >=1 TPM, even below RA's tenfold threshold;
    context QC is explicit. This is not a complete exposure or injury prediction.
    Bulk data use organ assessments and never imply selective germ-cell loss.
    No context rows means unavailable evidence, not absence or tolerability.
    """
    from . import tissue_tiers

    rows = target_selection.results(gene)
    if rows.empty:
        raise ValueError(f"Unknown protein label or member gene ID: {gene}")
    result = []
    for protein in rows.itertuples():
        encoded = protein.reproductive_expression_contexts
        if not isinstance(encoded, str):
            continue
        for context in json.loads(encoded):
            tissue, cell = context["tissue"], context["cell_type"]
            assessed = tissue_tiers.consequences(
                tissue,
                None if cell == "bulk" else cell,
                stage=context["stage"],
                resolution="bulk" if cell == "bulk" else "cell",
            )
            result.append(
                {"entity_id": protein.entity_id, "expression_context": context, **assessed}
            )
    return result


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
