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
SEX_SETS = (
    "RA-M",
    "RA-F",
    "RA-E",
    "RA-MF",
    "CTA-M",
    "CTA-F",
    "CTA-E",
    "CTA-MF",
    "MRA",
    "FRA",
    "ERA",
    "SRA",
)
CANDIDATE_SETS = ("CTA", "CA", "RA", "CA|RA", *SEX_SETS)
ALL_SETS = (*SETS, *SEX_SETS)
EXCLUSION_CATEGORIES = {
    "Unsafe": ("Unsafe",),
    "Other": ("Other",),
    "Unsafe|Other": ("Unsafe", "Other"),
}
SEQUENCE_SOURCES = ("auto", "sequence_bundle", "pyensembl")
NEAR_SELF_CONTRACT_ID = "oncoref.near-self-subtraction.2026-11-21.1"
ISOFORM_SELECTIVITY_CONTRACT_ID = "oncoref.isoform-selectivity.2026-11-22.1"
ENSEMBL_RELEASE = 115
_UNSET = object()
_BUNDLE_TABLE_CACHE: dict[str, pd.DataFrame] = {}


def _load_bundle_sequences_table() -> pd.DataFrame:
    """Return a cached occurrence table for the pinned current sequence bundle."""
    from . import sequence_bundle

    pin_id = sequence_bundle.descriptor()["bundle_id"]
    cached = _BUNDLE_TABLE_CACHE.get(pin_id)
    if cached is not None:
        return cached
    loaded = sequence_bundle.load(sequence_bundle.ensure(download=False)).to_sequences_table()
    _BUNDLE_TABLE_CACHE.clear()
    _BUNDLE_TABLE_CACHE[pin_id] = loaded
    return loaded


@dataclass(frozen=True)
class AntigenSelection:
    """One evaluated selection, shared by protein, gene, sequence and peptide views."""

    group: str
    _evaluated: pd.DataFrame
    _sequences_table: pd.DataFrame | None = None
    _species: str = "human"
    _sex_eligible_ids: frozenset[str] | None = None
    _sequence_source: str = "auto"

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

    def sequences(self, *, genome=None, sequence_source: str | None = None):
        return _sequences(
            self.proteins(),
            self.group,
            genome=genome,
            sequences_table=self._sequences_table,
            species=self._species,
            sequence_source=self._sequence_source if sequence_source is None else sequence_source,
        )

    def background_sequences(
        self, background="Unsafe|Other", *, genome=None, sequence_source: str | None = None
    ):
        return _sequences(
            self.background_proteins(background),
            background,
            genome=genome,
            sequences_table=self._sequences_table,
            species=self._species,
            sequence_source=self._sequence_source if sequence_source is None else sequence_source,
        )

    def background_index(
        self,
        background="Unsafe|Other",
        *,
        genome=None,
        sequence_source: str | None = None,
        policy=None,
        lengths=range(8, 16),
    ):
        from .peptide_subtraction import PeptideIndex

        records = self.background_sequences(
            background, genome=genome, sequence_source=sequence_source
        )
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
    sequence_source: str = "auto",
    **overrides,
):
    """Select once, then obtain consistent protein, gene, sequence and background views.

    ``host`` is a HostContext. An omitted host means catalog screening only.
    Tissue-loss choices never change the broad biological catalog or global Unsafe.
    """
    from . import target_filters
    from .species import resolve_species

    if sequence_source not in SEQUENCE_SOURCES:
        raise ValueError(
            f"Unknown sequence_source {sequence_source!r}; choose from {SEQUENCE_SOURCES}"
        )
    resolved = resolve_species(species)
    if group not in CANDIDATE_SETS:
        raise ValueError(f"Select a candidate set: {', '.join(CANDIDATE_SETS)}")
    if table is None or host_evidence is None:
        from .species_target_selection import (
            packaged_species_host_exposure,
            packaged_species_results,
        )

        if table is None:
            table = packaged_species_results(resolved)
        if host_evidence is None:
            host_evidence = packaged_species_host_exposure(resolved)
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
        _sequence_source=sequence_source,
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
    species="human",
):
    """Compatibility convenience view; prefer select() when using sequences/backgrounds."""
    from .species import resolve_species

    if group not in ALL_SETS:
        raise ValueError(f"Unknown antigen set {group!r}; choose from {ALL_SETS}")
    resolved = resolve_species(species)
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
        from .species_target_selection import packaged_species_results

        rows = packaged_species_results(resolved)
        result = rows[rows.category.isin(EXCLUSION_CATEGORIES[group])].reset_index(drop=True)
        result.attrs.update(
            catalog_policy_id=rows.attrs.get("policy_id", target_selection.POLICY_ID),
            set=group,
            species=resolved.name,
            taxon_id=resolved.taxon_id,
        )
        return result
    return select(group, mode=mode, host=host, species=resolved, **overrides).proteins()


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


def sequences(
    group="CTA",
    *,
    genome=None,
    species="human",
    sequence_source: str = "auto",
    sequences_table=None,
    **selection,
) -> pd.DataFrame:
    """Translations from selected loci, with an exact catalog sequence match flag."""
    protein_rows = (
        proteins(group)
        if species == "human" and not selection
        else proteins(group, species=species, **selection)
    )
    return _sequences(
        protein_rows,
        group,
        genome=genome,
        sequences_table=sequences_table,
        species=species,
        sequence_source=sequence_source,
    )


def forbidden_background(
    background="Unsafe|Other",
    *,
    host=None,
    mode="all",
    species="human",
    sequence_source: str = "auto",
    genome=None,
    sequences_table=None,
    **overrides,
) -> pd.DataFrame:
    """Return the version-pinned downstream forbidden/near-self background table (#640).

    Explicitly separates ``Unsafe`` from ``Other`` unless ``background='Unsafe|Other'``
    (the conservative default), never applies a blanket CTA exemption, incorporates any
    ``HostContext`` exclusions, fails closed on cross-species mismatch, and attaches
    explicit uncertainty/provenance fields in ``result.attrs``.
    """
    from .species import resolve_species

    if background not in EXCLUSION_CATEGORIES:
        raise ValueError("Background must be Unsafe, Other or Unsafe|Other")
    resolved = resolve_species(species)
    selected = select(
        "CA|RA",
        mode=mode,
        host=host,
        species=resolved,
        sequences_table=sequences_table,
        sequence_source=sequence_source,
        **overrides,
    )
    bg_proteins = selected.background_proteins(background)
    records = _sequences(
        bg_proteins,
        background,
        genome=genome,
        sequences_table=selected._sequences_table,
        species=resolved.name,
        sequence_source=sequence_source,
    )
    classifications = (
        bg_proteins["classification"].astype(str)
        if "classification" in bg_proteins
        else pd.Series("", index=bg_proteins.index, dtype="string")
    )
    reasons = (
        bg_proteins["background_reason"].astype(str)
        if "background_reason" in bg_proteins
        else pd.Series("", index=bg_proteins.index, dtype="string")
    )
    records.attrs.update(
        near_self_contract_id=NEAR_SELF_CONTRACT_ID,
        catalog_policy_id=target_selection.POLICY_ID,
        species=resolved.name,
        taxon_id=resolved.taxon_id,
        background=background,
        included_categories=list(EXCLUSION_CATEGORIES[background]),
        omitted_categories=[
            cat for cat in ("Unsafe", "Other") if cat not in EXCLUSION_CATEGORIES[background]
        ],
        cta_exemption_applied=False,
        unmeasured_home_entity_ids=sorted(
            bg_proteins.loc[classifications.eq("Other:Unscreened tissue"), "entity_id"].astype(str)
        ),
        review_hold_entity_ids=sorted(
            bg_proteins.loc[
                classifications.isin(["Other:Review", "Other:Sequence"]), "entity_id"
            ].astype(str)
        ),
        host_excluded_entity_ids=sorted(
            bg_proteins.loc[reasons.str.startswith("Host:"), "entity_id"].astype(str)
        ),
        downstream_policy_boundary=(
            "oncoref supplies version-pinned reference proteome sequences, catalog exclusions, "
            "host-filtered candidates, and explicit uncertainty metadata; MHC presentation, "
            "TCR cross-reactivity weighting, and clinical risk thresholds belong to downstream callers."
        ),
    )
    return records


def _sequences(
    protein_rows,
    group,
    *,
    genome=None,
    sequences_table=None,
    species="human",
    sequence_source: str = "auto",
) -> pd.DataFrame:
    """All annotated protein translations, one row per transcript occurrence.

    Uses an explicit ``sequences_table``, the verified offline ``sequence_bundle``
    (when ``sequence_source='sequence_bundle'`` or ``'auto'`` with a staged local
    human bundle), or a local Ensembl 115 reference (``oncoref[genome]``); never
    downloads assets implicitly. Missing genes/translations are recorded in
    ``attrs['missing_gene_ids']``. Stops/ambiguous residues remain in originals;
    peptide matching splits at them.
    """
    from .species import resolve_species

    if sequence_source not in SEQUENCE_SOURCES:
        raise ValueError(
            f"Unknown sequence_source {sequence_source!r}; choose from {SEQUENCE_SOURCES}"
        )
    resolved = resolve_species(species)
    if sequences_table is not None:
        table_taxon = sequences_table.attrs.get("taxon_id")
        if table_taxon is not None and int(table_taxon) != resolved.taxon_id:
            raise ValueError(
                f"Sequence reference taxon mismatch: sequences_table taxon_id={table_taxon} "
                f"cannot be used for species {resolved.name!r} (taxon_id={resolved.taxon_id})"
            )
    elif genome is None:
        if sequence_source == "sequence_bundle":
            from .species_target_selection import packaged_species_metadata

            meta = packaged_species_metadata(resolved)
            if meta is None or meta.get("sequence_provider") != "sequence_bundle":
                raise ValueError(
                    f"sequence_bundle currently pins human GRCh38/Ensembl 115 (taxon 9606), "
                    f"not {resolved.name!r}"
                )
            sequences_table = _load_bundle_sequences_table()
        elif sequence_source == "auto":
            from .species_target_selection import (
                packaged_species_metadata,
                packaged_species_sequences,
            )

            if packaged_species_metadata(resolved) is not None:
                try:
                    sequences_table = packaged_species_sequences(resolved)
                except FileNotFoundError:
                    sequences_table = None
    elif sequence_source == "sequence_bundle":
        raise ValueError("Cannot pass both genome=... and sequence_source='sequence_bundle'")
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
            sequence_source=sequences_table.attrs.get("sequence_source", "explicit_table"),
            missing_gene_ids=sorted(missing),
            sequence_scope="all annotated protein translations",
            **(
                {"bundle_id": sequences_table.attrs["bundle_id"]}
                if "bundle_id" in sequences_table.attrs
                else {}
            ),
            **(
                {"role_policy_id": sequences_table.attrs["role_policy_id"]}
                if "role_policy_id" in sequences_table.attrs
                else {}
            ),
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
        sequence_source="pyensembl",
        missing_gene_ids=sorted(missing),
        sequence_scope="all annotated protein translations",
    )
    return result


def isoform_selectivity_contract() -> dict:
    """Return the versioned isoform-selective target & locus-Unsafe-precedence contract (#736)."""
    return {
        "contract_id": ISOFORM_SELECTIVITY_CONTRACT_ID,
        "catalog_policy_id": target_selection.POLICY_ID,
        "default_max_normal_isoform_tpm": 1.0,
        "default_min_tumor_isoform_tpm": 10.0,
        "invariants": [
            "Gene-level Unsafe precedence NEVER yields to an isoform nomination without transcript/junction-level healthy-normal quantification",
            "Even when healthy-normal transcript expression is absent, only peptides subtracted against the full Unsafe proteome (including canonical Unsafe isoforms of the same locus) may qualify",
            "Qualifying an isoform-unique junction/exon peptide NEVER reclassifies the parent locus or canonical Unsafe protein",
        ],
    }


def _extract_isoform_max_tpm(
    table: pd.DataFrame | None,
    *,
    value_col: str,
) -> dict[str, float]:
    if table is None or table.empty:
        return {}
    if "transcript_id" not in table.columns:
        raise ValueError("Isoform expression table requires a 'transcript_id' column")
    work = table.copy()
    tids = work["transcript_id"].astype(str).str.strip().str.split(".", n=1).str[0]
    if value_col in work.columns:
        vals = pd.to_numeric(work[value_col], errors="coerce")
    else:
        num_cols = [c for c in work.columns if c not in {"transcript_id", "gene_id", "Symbol"}]
        if not num_cols:
            raise ValueError(
                f"Isoform expression table requires {value_col!r} or numeric sample/tissue columns"
            )
        vals = work[num_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
    grouped = pd.DataFrame({"tid": tids, "val": vals}).groupby("tid", sort=False)["val"].max()
    return {str(k): float(v) for k, v in grouped.items() if pd.notna(v)}


def evaluate_isoform_selectivity(
    candidate_sequences: pd.DataFrame,
    *,
    unsafe_index,
    normal_isoform_tpm: pd.DataFrame | None = None,
    tumor_isoform_tpm: pd.DataFrame | None = None,
    max_normal_tpm: float = 1.0,
    min_tumor_tpm: float = 10.0,
    lengths=range(8, 16),
) -> pd.DataFrame:
    """Evaluate isoform-selective peptide targetability without weakening locus Unsafe rules (#736).

    Enforces three conjunctive gates on every candidate transcript row:
    1. Transcript/junction-resolved healthy normal RNA quantification must be supplied
       for ``transcript_id`` and satisfy ``observed_max_normal_tpm <= max_normal_tpm``.
    2. Subtracting ``sequence`` against ``unsafe_index`` (which must include any
       normal-expressed canonical isoforms of the same gene) must leave >=1 unique k-mer.
    3. Transcript-resolved tumor expression must satisfy ``observed_max_tumor_tpm >= min_tumor_tpm``.
    """
    required = {"transcript_id", "sequence"}
    missing_cols = sorted(required - set(candidate_sequences.columns))
    if missing_cols:
        raise ValueError(f"candidate_sequences lacks required columns: {missing_cols}")

    normal_max = _extract_isoform_max_tpm(normal_isoform_tpm, value_col="max_normal_tpm")
    tumor_max = _extract_isoform_max_tpm(tumor_isoform_tpm, value_col="max_tumor_tpm")

    rows = []
    for rec in candidate_sequences.to_dict("records"):
        raw_tid = str(rec["transcript_id"]).strip()
        tid = raw_tid.split(".", 1)[0]
        gene_id = str(rec.get("gene_id") or "")
        seq = str(rec["sequence"])
        subtracted = unsafe_index.subtract(seq, lengths=lengths)
        unique_peptides = sorted(set(subtracted["peptide"].astype(str))) if len(subtracted) else []
        n_retained = len(unique_peptides)

        obs_norm = normal_max.get(tid)
        obs_tum = tumor_max.get(tid)

        if obs_norm is None:
            status = "held_missing_normal_isoform_quantification"
            qualifies = False
        elif obs_norm > max_normal_tpm:
            status = "excluded_normal_isoform_expression"
            qualifies = False
        elif n_retained == 0:
            status = "excluded_no_unsafe_unique_peptides"
            qualifies = False
        elif obs_tum is None or obs_tum < min_tumor_tpm:
            status = "held_insufficient_tumor_isoform_expression"
            qualifies = False
        else:
            status = "qualifies_isoform_selective_peptides_only"
            qualifies = True

        rows.append(
            {
                "gene_id": gene_id,
                "transcript_id": raw_tid,
                "observed_max_normal_isoform_tpm": obs_norm if obs_norm is not None else pd.NA,
                "observed_max_tumor_isoform_tpm": obs_tum if obs_tum is not None else pd.NA,
                "n_retained_unique_kmers": n_retained,
                "retained_unique_peptides": ";".join(unique_peptides),
                "isoform_status": status,
                "qualifies_isoform_selective": qualifies,
                "locus_unsafe_precedence_preserved": True,
            }
        )
    out = pd.DataFrame(rows)
    out.attrs.update(
        contract_id=ISOFORM_SELECTIVITY_CONTRACT_ID,
        catalog_policy_id=target_selection.POLICY_ID,
        max_normal_tpm=max_normal_tpm,
        min_tumor_tpm=min_tumor_tpm,
    )
    return out
