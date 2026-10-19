"""Auditable, annotation-scoped Ensembl gene identities (contract version 1).

The legacy alias resolver harmonizes tables without an annotation. This module
verifies those explicit ID mappings against one supplied annotation before they
can affect CTA admission or negative-reference exclusions. Symbols and identical
proteins are evidence checks, never a source of new aliases.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from functools import lru_cache

from .gene_ids import unversioned
from .load_dataset import _register_derived_cache, get_data
from .version import __version__

GENE_IDENTITY_CONTRACT_VERSION = 1
_PRIMARY_CONTIGS = frozenset([*(str(i) for i in range(1, 23)), "X", "Y", "MT"])


def _digest(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@lru_cache(maxsize=1)
def _reference():
    aliases = get_data("ensembl-id-aliases", copy=False).fillna("")
    canonical = get_data("canonical-gene-space", copy=False).fillna("")
    alias_rows = sorted(
        (unversioned(a), unversioned(p), str(s), str(m))
        for a, p, s, m in zip(
            aliases.alt_haplotype_id, aliases.primary_contig_id, aliases.symbol, aliases.source
        )
    )
    canonical_rows = sorted(
        (unversioned(g), str(s), str(c), str(r))
        for g, s, c, r in zip(
            canonical.ensembl_gene_id,
            canonical.symbol,
            canonical.seqname,
            canonical.ensembl_release,
        )
    )
    by_alias = {}
    for a, p, s, m in alias_rows:
        by_alias.setdefault(a, set()).add((p, s, m))
    by_id = {g: (s, c, r) for g, s, c, r in canonical_rows}
    primary_by_symbol = {}
    for g, s, c, _ in canonical_rows:
        if s and c in _PRIMARY_CONTIGS:
            primary_by_symbol.setdefault(s, set()).add(g)
    return (
        by_alias,
        by_id,
        primary_by_symbol,
        _digest(alias_rows),
        _digest(canonical_rows),
        tuple(sorted({r for _, _, _, r in canonical_rows})),
    )


_register_derived_cache(_reference.cache_clear)


@dataclass(frozen=True)
class GeneIdentity:
    """One gene-ID decision with original source identity and mapping evidence.

    Only ``canonical`` and ``alias`` are verified. Other outcomes have no resolved
    canonical ID; ``candidate_gene_ids`` records competing/unverified targets.
    ``as_dict`` returns a JSON-serializable audit record. Reference hashes cover
    sorted logical rows of the alias/canonical tables, not their compression.
    """

    input_gene_id: str
    source_gene_id: str
    canonical_gene_id: str | None
    candidate_gene_ids: tuple[str, ...]
    status: str
    mapping_method: str
    source_gene_name: str
    source_contig: str
    source_start: int | None
    source_end: int | None
    source_strand: str
    annotation_species: str
    annotation_assembly: str
    annotation_release: int
    alias_rows: tuple[tuple[str, str, str], ...]
    aliases_sha256: str
    canonical_gene_space_sha256: str
    canonical_reference_releases: tuple[str, ...]
    oncoref_version: str = __version__
    identity_contract_version: int = GENE_IDENTITY_CONTRACT_VERSION

    @property
    def verified(self) -> bool:
        return self.status in {"canonical", "alias"}

    def as_dict(self) -> dict:
        return asdict(self)


def _annotation(genome):
    species = getattr(getattr(genome, "species", None), "latin_name", None)
    assembly = getattr(genome, "reference_name", None)
    release = getattr(genome, "release", None)
    if species != "homo_sapiens" or assembly not in {"GRCh37", "GRCh38"}:
        raise ValueError("gene identity requires an explicit human GRCh37/GRCh38 annotation")
    if not isinstance(release, int) or release < 1:
        raise ValueError("gene identity requires an explicit positive Ensembl release")
    return species, assembly, release


def resolve_gene_identity(gene_id: str, *, genome) -> GeneIdentity:
    """Resolve an explicit Ensembl ID against one human Ensembl annotation.

    A shipped alternate-contig mapping needs its recorded symbol, a non-primary
    source contig, and a unique canonical primary-locus target. A historical ID
    additionally needs the same primary chromosome and no coexisting successor.
    Unknown methods, competing targets, inconsistent annotation evidence, and
    unknown same-symbol loci remain unresolved. Lookup errors other than the
    annotation's missing-ID ``ValueError`` propagate to the caller.
    """
    species, assembly, release = _annotation(genome)
    source_id = unversioned(gene_id).upper()
    by_alias, by_id, primary_by_symbol, aliases_hash, canonical_hash, releases = _reference()
    mappings = tuple(sorted(by_alias.get(source_id, ())))
    candidates = tuple(sorted({row[0] for row in mappings}))
    status, method, resolved = "unmapped", "", None
    try:
        gene = genome.gene_by_id(source_id)
    except ValueError:
        gene = None
    name = str(getattr(gene, "gene_name", "") or "")
    contig = str(getattr(gene, "contig", "") or "")
    if gene is None:
        status = "absent"
    elif len(candidates) > 1 or (source_id in by_id and mappings):
        status = "ambiguous"
    elif source_id in by_id:
        if contig == by_id[source_id][1]:
            status, method, resolved = "canonical", "stable_id", source_id
        else:
            status = "annotation_conflict"
    elif mappings:
        target = candidates[0]
        methods = {row[2] for row in mappings}
        symbols = {row[1] for row in mappings}
        method = ",".join(sorted(methods))
        if len(methods) != 1 or len(symbols) != 1:
            status = "ambiguous"
        elif target not in by_id or method not in {"ensembl_alt_contig", "ensembl_id_history_name"}:
            status = "unverified"
        else:
            symbol = next(iter(symbols))
            target_symbol, target_contig, _ = by_id[target]
            valid = (
                bool(symbol)
                and name == symbol == target_symbol
                and primary_by_symbol.get(symbol) == {target}
                and target_contig in _PRIMARY_CONTIGS
            )
            try:
                target_gene = genome.gene_by_id(target)
            except ValueError:
                target_gene = None
            if target_gene is not None:
                valid = valid and (
                    str(target_gene.contig) == target_contig and target_gene.gene_name == symbol
                )
            if method == "ensembl_alt_contig":
                valid = valid and bool(contig) and contig not in _PRIMARY_CONTIGS
            else:
                valid = valid and contig == target_contig and target_gene is None
            if valid:
                status, resolved = "alias", target
            else:
                status = "annotation_conflict"
    return GeneIdentity(
        input_gene_id=str(gene_id),
        source_gene_id=source_id,
        canonical_gene_id=resolved,
        candidate_gene_ids=candidates,
        status=status,
        mapping_method=method,
        source_gene_name=name,
        source_contig=contig,
        source_start=getattr(gene, "start", None),
        source_end=getattr(gene, "end", None),
        source_strand=str(getattr(gene, "strand", "") or ""),
        annotation_species=species,
        annotation_assembly=assembly,
        annotation_release=release,
        alias_rows=mappings,
        aliases_sha256=aliases_hash,
        canonical_gene_space_sha256=canonical_hash,
        canonical_reference_releases=releases,
    )


def cta_annotation_gene_identities(genome) -> tuple[GeneIdentity, ...]:
    """Audit annotated IDs potentially belonging to the full CTA candidate set.

    Includes unresolved mappings with a CTA among their candidate IDs, so rejected
    aliases remain visible. Unrelated genes are not enumerated in this CTA audit.
    Original transcript/protein/occurrence IDs must still be retained by consumers.
    """
    from .cta import cta_unfiltered_gene_ids

    _annotation(genome)
    candidates = {unversioned(g) for g in cta_unfiltered_gene_ids()}
    aliases = _reference()[0]
    relevant = candidates | {
        a for a, mappings in aliases.items() if any(p in candidates for p, _, _ in mappings)
    }
    return tuple(
        resolve_gene_identity(g, genome=genome)
        for g in sorted(set(genome.gene_ids()))
        if unversioned(g).upper() in relevant
    )


def cta_annotation_gene_ids(genome, *, unfiltered: bool = False) -> set[str]:
    """Verified annotation source IDs for default CTA admission or self exclusion.

    ``unfiltered=True`` selects the full candidate universe for negative self
    references; it does not change which targets the default panel admits.
    Ambiguous/unverified sources remain in the negative reference.
    """
    from .cta import cta_gene_ids, cta_unfiltered_gene_ids

    selected = cta_unfiltered_gene_ids() if unfiltered else cta_gene_ids()
    selected = {unversioned(g) for g in selected}
    return {
        record.source_gene_id
        for record in cta_annotation_gene_identities(genome)
        if record.verified and record.canonical_gene_id in selected
    }
