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
SEX_SETS = ("MRA", "FRA", "ERA", "SRA", "CTA-M", "CTA-F", "CTA-E", "CTA-MF")
CANDIDATE_SETS = ("CTA", "CA", "RA", "CA|RA", *SEX_SETS)
ALL_SETS = (*SETS, *SEX_SETS)
EXCLUSION_CATEGORIES = {
    "Unsafe": ("Unsafe",),
    "Other": ("Other",),
    "Unsafe|Other": ("Unsafe", "Other"),
}
ENSEMBL_RELEASE = 115
_UNSET = object()


@dataclass(frozen=True)
class AntigenSelection:
    """One evaluated selection, shared by protein, gene, sequence and peptide views."""

    group: str
    _evaluated: pd.DataFrame
    _sequences_table: pd.DataFrame | None = None
    _species: str = "human"
    _sex_eligible_ids: frozenset[str] | None = None

    @property
    def metadata(self):
        from copy import deepcopy

        return deepcopy(
            {
                "catalog_policy_id": target_selection.POLICY_ID,
                "set": self.group,
                **self._evaluated.attrs,
            }
        )

    def decisions(self):
        """All catalog rows with original membership, host status and filter reasons."""
        return self._evaluated.copy(deep=True)

    def proteins(self):
        d = self._evaluated
        ca, ra = d.filtered_cancer_antigen, d.filtered_reproductive_antigen
        if self.group in SEX_SETS:
            from .species_target_selection import SEX_SETS as _SEX_SPEC

            _sex_class, requires_cta = _SEX_SPEC[self.group]
            base_member = (ca & ra) if requires_cta else ra
            member = base_member & d.entity_id.isin(self._sex_eligible_ids or frozenset())
        else:
            member = {"CTA": ca & ra, "CA": ca, "RA": ra, "CA|RA": ca | ra}[self.group]
        result = d[d.keep & member].copy()
        for name in (
            "category",
            "classification",
            "is_cta",
            "is_cancer_antigen",
            "is_reproductive_antigen",
            "is_target_candidate",
        ):
            result["catalog_" + name] = result[name]
        result["category"] = result.filtered_category
        result["classification"] = result.category
        result["is_cancer_antigen"] = result.filtered_cancer_antigen
        result["is_reproductive_antigen"] = result.filtered_reproductive_antigen
        result["is_cta"] = result.is_cancer_antigen & result.is_reproductive_antigen
        result["is_target_candidate"] = result.is_cancer_antigen | result.is_reproductive_antigen
        if "cta_tier" in result:
            result["catalog_cta_tier"] = result.cta_tier
            result["cta_tier"] = result.cta_tier.where(result.is_cta, "")
        result.attrs.update(self.metadata)
        return result.reset_index(drop=True)

    def gene_ids(self):
        return frozenset(
            g for members in self.proteins().member_gene_ids for g in json.loads(members)
        )

    def background_proteins(self, background="Unsafe|Other"):
        """Global exclusions plus host-excluded/unresolved candidates, never a CTA exemption."""
        if background not in EXCLUSION_CATEGORIES:
            raise ValueError("Background must be Unsafe, Other or Unsafe|Other")
        d = self._evaluated
        base = d.category.isin(EXCLUSION_CATEGORIES[background])
        host = d.is_target_candidate & ~d.pass_host
        result = d[base | host].copy()
        result["background_reason"] = result.category.where(base, "Host: " + result.host_status)
        result.attrs.update(
            self.metadata,
            background=background,
            background_scope="catalog exclusions plus host-excluded or unresolved candidates",
        )
        return result.reset_index(drop=True)

    def sequences(self, *, genome=None):
        return _sequences(
            self.proteins(),
            self.group,
            genome=genome,
            sequences_table=self._sequences_table,
            species=self._species,
        )

    def background_sequences(self, background="Unsafe|Other", *, genome=None):
        return _sequences(
            self.background_proteins(background),
            background,
            genome=genome,
            sequences_table=self._sequences_table,
            species=self._species,
        )

    def background_index(
        self, background="Unsafe|Other", *, genome=None, policy=None, lengths=range(8, 16)
    ):
        from .peptide_subtraction import PeptideIndex

        records = self.background_sequences(background, genome=genome)
        if records.attrs["missing_gene_ids"]:
            raise ValueError(
                "Background has missing gene translations; subtraction coverage is incomplete"
            )
        return PeptideIndex(records, policy=policy, lengths=lengths)


def select(
    group="CTA",
    *,
    mode="all",
    host=None,
    table=None,
    host_evidence=None,
    species="human",
    sequences_table=None,
    **overrides,
):
    """Select once, then obtain consistent protein, gene, sequence and background views.

    ``host`` is a HostContext. An omitted host means catalog screening only.
    Tissue-loss choices never change the broad biological catalog or global Unsafe.
    """
    from . import target_filters
    from .species import resolve_species

    resolved = resolve_species(species)
    if group not in CANDIDATE_SETS:
        raise ValueError(f"Select a candidate set: {', '.join(CANDIDATE_SETS)}")
    if table is None and resolved.taxon_id != 9606:
        from .species_target_selection import packaged_species_results

        table = packaged_species_results(resolved)
    sex_ids: frozenset[str] | None = None
    if group in SEX_SETS:
        from .species_target_selection import SEX_SETS as _SEX_SPEC
        from .species_target_selection import classify_sex_restriction

        sex_class, _ = _SEX_SPEC[group]
        sex_df = classify_sex_restriction(table, species=resolved, host_evidence=host_evidence)
        sex_ids = frozenset(sex_df.loc[sex_df.sex_restriction_class.eq(sex_class), "entity_id"])
    base = target_filters.FILTER_MODES[mode] if isinstance(mode, str) else mode
    if host is not None:
        overrides["host"] = host
    filt = replace(base, **overrides)
    return AntigenSelection(
        group,
        target_filters.apply(filt, table, host_evidence=host_evidence),
        _sequences_table=sequences_table,
        _species=resolved.name,
        _sex_eligible_ids=sex_ids,
    )


def proteins(
    group="CTA",
    *,
    profile=_UNSET,
    lactating=_UNSET,
    pregnant=_UNSET,
    allow_tissues=_UNSET,
    protect_tissues=_UNSET,
    absent_tissues=_UNSET,
    host=None,
    mode="all",
):
    """Compatibility convenience view; prefer select() when using sequences/backgrounds."""
    if group not in ALL_SETS:
        raise ValueError(f"Unknown antigen set {group!r}; choose from {ALL_SETS}")
    # Omitted convenience arguments must preserve a caller's TargetFilter.
    # Explicit False/()/"all" remain valid overrides of that filter.
    overrides = {
        key: value
        for key, value in {
            "host_profile": profile,
            "lactating": lactating,
            "pregnant": pregnant,
            "allow_tissues": allow_tissues,
            "protect_tissues": protect_tissues,
            "absent_tissues": absent_tissues,
        }.items()
        if value is not _UNSET
    }
    active = (
        any(
            value != "all" if key == "host_profile" else bool(value)
            for key, value in overrides.items()
        )
        or host is not None
        or mode != "all"
    )
    if group in EXCLUSION_CATEGORIES:
        if active:
            raise ValueError(
                "Host profiles and candidate filters require select(...).background_proteins()"
            )
        rows = target_selection.results()
        return rows[rows.category.isin(EXCLUSION_CATEGORIES[group])].reset_index(drop=True)
    if not active and group in SETS:
        rows = target_selection.results()
        return rows[rows[SETS[group]]].reset_index(drop=True)
    return select(group, mode=mode, host=host, **overrides).proteins()


def gene_ids(group="CTA", **kwargs):
    """All member gene IDs from the same selection contract."""
    return frozenset(
        g for members in proteins(group, **kwargs).member_gene_ids for g in json.loads(members)
    )


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
    """Translations from selected loci, with an exact catalog sequence match flag."""
    return _sequences(proteins(group, **selection), group, genome=genome)


def _sequences(
    protein_rows, group, *, genome=None, sequences_table=None, species="human"
) -> pd.DataFrame:
    """All annotated protein translations, one row per transcript occurrence.

    Uses an explicit ``sequences_table``, an ``orfdata.Proteome``, or a local
    Ensembl 115 reference (oncoref[genome]); never downloads assets. Missing
    genes/translations are recorded in ``attrs['missing_gene_ids']``.
    Stops/ambiguous residues remain in originals; peptide matching splits at them.
    """
    from .species import resolve_species

    resolved = resolve_species(species)
    if sequences_table is not None:
        by_gene = (
            dict(iter(sequences_table.groupby("gene_id")))
            if not sequences_table.empty and "gene_id" in sequences_table
            else {}
        )
        rows, missing = [], []
        for protein in protein_rows.itertuples():
            raw_hashes = getattr(protein, "representative_sequence_sha256s", None)
            representatives = None if pd.isna(raw_hashes) else set(json.loads(raw_hashes))
            for gene_id in json.loads(protein.member_gene_ids):
                sub = by_gene.get(gene_id)
                if sub is None or sub.empty:
                    missing.append(gene_id)
                    continue
                for rec in sub.itertuples(index=False):
                    sequence = str(rec.sequence)
                    if not sequence:
                        continue
                    digest = hashlib.sha256(sequence.encode()).hexdigest()
                    rows.append(
                        {
                            "entity_id": protein.entity_id,
                            "gene_id": gene_id,
                            "transcript_id": str(getattr(rec, "transcript_id", gene_id)),
                            "transcript_biotype": str(
                                getattr(rec, "transcript_biotype", "protein_coding")
                            ),
                            "sequence": sequence,
                            "sequence_sha256": digest,
                            "matches_representative_sequence": (
                                pd.NA if representatives is None else digest in representatives
                            ),
                            "set": group,
                            "catalog_category": getattr(
                                protein, "catalog_category", getattr(protein, "category", "")
                            ),
                            "selection_category": getattr(protein, "category", ""),
                            "background_reason": getattr(protein, "background_reason", ""),
                        }
                    )
        result = pd.DataFrame(
            rows,
            columns=[
                "entity_id",
                "gene_id",
                "transcript_id",
                "transcript_biotype",
                "sequence",
                "sequence_sha256",
                "matches_representative_sequence",
                "set",
                "catalog_category",
                "selection_category",
                "background_reason",
            ],
        )
        result["matches_representative_sequence"] = result.matches_representative_sequence.astype(
            "boolean"
        )
        result.attrs.update(protein_rows.attrs)
        result.attrs.update(
            policy_id=target_selection.POLICY_ID,
            reference=sequences_table.attrs.get(
                "reference", sequences_table.attrs.get("reference_id", f"{resolved.name}:explicit")
            ),
            missing_gene_ids=sorted(missing),
            sequence_scope="all annotated protein translations",
        )
        return result

    expected_latin = {9606: "homo_sapiens", 10090: "mus_musculus"}.get(resolved.taxon_id)
    expected_ref = {
        9606: "homo_sapiens:GRCh38:Ensembl:115",
        10090: "mus_musculus:GRCm39:Ensembl:115",
    }.get(resolved.taxon_id)
    if genome is None:
        if expected_latin is None:
            raise ValueError(
                f"Species {resolved.name!r} requires an explicit sequences_table or genome"
            )
        from pyensembl import EnsemblRelease

        py_species = "human" if resolved.taxon_id == 9606 else "mouse"
        genome = EnsemblRelease(ENSEMBL_RELEASE, species=py_species)
    if (
        expected_latin is None
        or genome.release != ENSEMBL_RELEASE
        or genome.species.latin_name != expected_latin
    ):
        raise ValueError(f"Current antigen sets require {resolved.name} Ensembl 115")
    rows, missing = [], []
    for protein in protein_rows.itertuples():
        raw_hashes = getattr(protein, "representative_sequence_sha256s", None)
        representatives = None if pd.isna(raw_hashes) else set(json.loads(raw_hashes))
        for gene_id in json.loads(protein.member_gene_ids):
            found = False
            try:
                transcripts = genome.gene_by_id(gene_id).transcripts
            except ValueError:
                transcripts = []
            for transcript in transcripts:
                sequence = transcript.protein_sequence
                if not sequence:
                    continue
                found = True
                digest = hashlib.sha256(sequence.encode()).hexdigest()
                rows.append(
                    {
                        "entity_id": protein.entity_id,
                        "gene_id": gene_id,
                        "transcript_id": transcript.id,
                        "transcript_biotype": transcript.biotype,
                        "sequence": sequence,
                        "sequence_sha256": digest,
                        "matches_representative_sequence": (
                            pd.NA if representatives is None else digest in representatives
                        ),
                        "set": group,
                        "catalog_category": getattr(
                            protein, "catalog_category", getattr(protein, "category", "")
                        ),
                        "selection_category": getattr(protein, "category", ""),
                        "background_reason": getattr(protein, "background_reason", ""),
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
            "matches_representative_sequence",
            "set",
            "catalog_category",
            "selection_category",
            "background_reason",
        ],
    )
    result["matches_representative_sequence"] = result.matches_representative_sequence.astype(
        "boolean"
    )
    result.attrs.update(protein_rows.attrs)
    result.attrs.update(
        policy_id=target_selection.POLICY_ID,
        reference=expected_ref,
        missing_gene_ids=sorted(missing),
        sequence_scope="all annotated protein translations",
    )
    return result
