"""Current antigen memberships and all annotated translations, with explicit scope.

These APIs use target_selection, never the historical target_sets snapshot.
Sequences inherit gene-level evidence; membership does not establish translation
of every isoform, antigen presentation, or therapeutic safety.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, replace

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


@dataclass(frozen=True)
class AntigenSelection:
    """One evaluated selection, shared by protein, gene, sequence and peptide views."""

    group: str
    _evaluated: pd.DataFrame

    @property
    def metadata(self):
        from copy import deepcopy
        return deepcopy({"catalog_policy_id": target_selection.POLICY_ID, "set": self.group, **self._evaluated.attrs})

    def decisions(self):
        """All catalog rows with original membership, host status and filter reasons."""
        return self._evaluated.copy(deep=True)

    def proteins(self):
        d = self._evaluated
        ca, ra = d.filtered_cancer_antigen, d.filtered_reproductive_antigen
        member = {"CTA": ca & ra, "CA": ca, "RA": ra, "CA|RA": ca | ra}[self.group]
        result = d[d.keep & member].copy()
        for name in ("category", "classification", "is_cta", "is_cancer_antigen", "is_reproductive_antigen", "is_target_candidate"):
            result["catalog_" + name] = result[name]
        result["category"] = result.filtered_category
        result["classification"] = result.category
        result["is_cancer_antigen"] = result.filtered_cancer_antigen
        result["is_reproductive_antigen"] = result.filtered_reproductive_antigen
        result["is_cta"] = result.is_cancer_antigen & result.is_reproductive_antigen
        result["is_target_candidate"] = result.is_cancer_antigen | result.is_reproductive_antigen
        result.attrs.update(self.metadata)
        return result.reset_index(drop=True)

    def gene_ids(self):
        return frozenset(g for members in self.proteins().member_gene_ids for g in json.loads(members))

    def background_proteins(self, background="Unsafe|Other"):
        """Global exclusions plus host-excluded/unresolved candidates, never a CTA exemption."""
        if background not in EXCLUSION_CATEGORIES:
            raise ValueError("Background must be Unsafe, Other or Unsafe|Other")
        d = self._evaluated
        base = d.category.isin(EXCLUSION_CATEGORIES[background])
        host = d.is_target_candidate & ~d.pass_host
        result = d[base | host].copy()
        result["background_reason"] = result.category.where(base, "Host: " + result.host_status)
        result.attrs.update(self.metadata, background=background, background_scope="catalog exclusions plus host-excluded or unresolved candidates")
        return result.reset_index(drop=True)

    def sequences(self, *, genome=None):
        return _sequences(self.proteins(), self.group, genome=genome)

    def background_sequences(self, background="Unsafe|Other", *, genome=None):
        return _sequences(self.background_proteins(background), background, genome=genome)

    def background_index(self, background="Unsafe|Other", *, genome=None, policy=None, lengths=range(8, 16)):
        from .peptide_subtraction import PeptideIndex
        records = self.background_sequences(background, genome=genome)
        if records.attrs["missing_gene_ids"]:
            raise ValueError("Background has missing gene translations; subtraction coverage is incomplete")
        return PeptideIndex(records, policy=policy, lengths=lengths)


def select(group="CTA", *, mode="all", host=None, table=None, host_evidence=None, **overrides):
    """Select once, then obtain consistent protein, gene, sequence and background views.

    ``host`` is a HostContext. An omitted host means catalog screening only.
    Tissue-loss choices never change the broad biological catalog or global Unsafe.
    """
    from . import target_filters
    if group not in SETS or group in EXCLUSION_CATEGORIES:
        raise ValueError("Select a candidate set: CTA, CA, RA or CA|RA")
    base = target_filters.FILTER_MODES[mode] if isinstance(mode, str) else mode
    if host is not None:
        overrides["host"] = host
    filt = replace(base, **overrides)
    return AntigenSelection(group, target_filters.apply(filt, table, host_evidence=host_evidence))


def proteins(group="CTA", *, profile="all", lactating=False, pregnant=False, allow_tissues=(), protect_tissues=(), absent_tissues=(), host=None, mode="all"):
    """Compatibility convenience view; prefer select() when using sequences/backgrounds."""
    if group not in SETS:
        raise ValueError(f"Unknown antigen set {group!r}; choose from {tuple(SETS)}")
    active = any((profile != "all", lactating, pregnant, allow_tissues, protect_tissues, absent_tissues, host is not None, mode != "all"))
    if group in EXCLUSION_CATEGORIES:
        if active:
            raise ValueError("Host profiles and candidate filters require select(...).background_proteins()")
        rows = target_selection.results()
        return rows[rows.category.isin(EXCLUSION_CATEGORIES[group])].reset_index(drop=True)
    if not active:
        rows = target_selection.results()
        return rows[rows[SETS[group]]].reset_index(drop=True)
    return select(group, mode=mode, host=host, host_profile=profile, lactating=lactating, pregnant=pregnant, allow_tissues=allow_tissues, protect_tissues=protect_tissues, absent_tissues=absent_tissues).proteins()


def gene_ids(group="CTA", **kwargs):
    """All member gene IDs from the same selection contract."""
    return frozenset(g for members in proteins(group, **kwargs).member_gene_ids for g in json.loads(members))


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


def sequences(group="CTA", *, genome=None, **selection) -> pd.DataFrame:
    """Annotated translations of the exact same profiled protein selection."""
    return _sequences(proteins(group, **selection), group, genome=genome)


def _sequences(protein_rows, group, *, genome=None) -> pd.DataFrame:
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
    for protein in protein_rows.itertuples():
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
    result.attrs.update(protein_rows.attrs)
    result.attrs.update(
        policy_id=target_selection.POLICY_ID,
        reference="homo_sapiens:GRCh38:Ensembl:115",
        missing_gene_ids=sorted(missing),
        sequence_scope="all annotated protein translations",
    )
    return result
