"""Exact and bounded weighted edit-distance subtraction against protein sets.

Distances are screening parameters, not probabilities of TCR cross-reactivity.
Coordinates are zero-based, half-open in the original translation. Invalid
residues and stops break segments; matches never cross them or protein records.
"""

from __future__ import annotations

import bisect
import math
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from types import MappingProxyType

import pandas as pd

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
SEGMENTS = re.compile("[ACDEFGHIKLMNPQRSTVWY]+")
NEAR_SELF_CONTRACT_ID = "oncoref.near-self-subtraction.2026-11-21.1"


def near_self_contract() -> dict:
    """Return the versioned downstream near-self / forbidden-background contract (#640)."""
    from . import sequence_bundle, target_selection

    return {
        "contract_id": NEAR_SELF_CONTRACT_ID,
        "catalog_policy_id": target_selection.POLICY_ID,
        "sequence_role_policy_id": sequence_bundle.ROLE_POLICY_ID,
        "default_background": "Unsafe|Other",
        "supported_backgrounds": ["Unsafe", "Other", "Unsafe|Other"],
        "sequence_bundle_roles": {
            "Unsafe": ["exclude"],
            "Unsafe|Other": ["exclude", "review"],
        },
        "invariants": [
            "Other is never silently relabeled Unsafe; callers choose 'Unsafe|Other' (conservative default) or 'Unsafe' explicitly",
            "No blanket CTA or candidate-antigen exemption is ever applied to background subtraction",
            "Host-excluded or host-unresolved candidates join the background when a HostContext is supplied",
            "Missing gene translations fail PeptideIndex construction in background_index() unless explicitly inspected",
            "Cross-species sequence tables and catalog selections fail closed",
        ],
        "downstream_responsibilities": [
            "MHC binding and antigen processing",
            "Position/HLA/TCR-aware cross-reactivity scoring",
            "Sample-specific somatic/germline genotype personalization beyond population AF filters",
            "Therapeutic risk-benefit adjudication",
        ],
        "population_polymorphism_support": {
            "helper": "PeptideIndex.with_coding_variants / apply_coding_variants",
            "default_min_af": 0.001,
            "coordinate_system": "1-based reference amino-acid position, validated fail-closed",
        },
    }


@dataclass(frozen=True)
class CodingVariant:
    """One germline missense polymorphism on a reference transcript or protein record (#736).

    ``position`` is 1-based in the reference amino-acid translation and is verified
    fail-closed against ``ref_aa`` before generating flanking background peptides.
    """

    target_id: str
    position: int
    ref_aa: str
    alt_aa: str
    af: float
    source_id: str
    variant_id: str = ""

    def __post_init__(self) -> None:
        tid = str(self.target_id).strip()
        if not tid:
            raise ValueError("CodingVariant target_id must be nonempty")
        if type(self.position) is not int or self.position < 1:
            raise ValueError("CodingVariant position must be a 1-based positive integer")
        if len(self.ref_aa) != 1 or self.ref_aa not in AMINO_ACIDS:
            raise ValueError(f"CodingVariant ref_aa must be a canonical AA; got {self.ref_aa!r}")
        if len(self.alt_aa) != 1 or self.alt_aa not in AMINO_ACIDS:
            raise ValueError(f"CodingVariant alt_aa must be a canonical AA; got {self.alt_aa!r}")
        if self.ref_aa == self.alt_aa:
            raise ValueError("CodingVariant ref_aa and alt_aa must differ")
        if not math.isfinite(self.af) or not (0.0 <= self.af <= 1.0):
            raise ValueError("CodingVariant af must be finite and within [0.0, 1.0]")
        src = str(self.source_id).strip()
        if not src:
            raise ValueError("CodingVariant source_id must be nonempty")
        vid = str(self.variant_id).strip() or f"{tid}:p.{self.ref_aa}{self.position}{self.alt_aa}"
        object.__setattr__(self, "target_id", tid)
        object.__setattr__(self, "source_id", src)
        object.__setattr__(self, "variant_id", vid)


def apply_coding_variants(
    records: pd.DataFrame,
    variants: Iterable[CodingVariant],
    *,
    min_af: float = 0.001,
    flank: int = 14,
    strict_targets: bool = True,
) -> pd.DataFrame:
    """Augment a background sequence table with common population missense windows (#736).

    Retains all original reference rows (`is_coding_variant=False`) and appends
    flanking mutated sequence windows (`[pos - 1 - flank : pos + flank]`) for every
    variant with ``af >= min_af``. Fails closed if any variant coordinate is out of
    bounds or disagrees with ``ref_aa``, or (when ``strict_targets=True``) if a
    qualifying ``target_id`` matches no row in ``records``.
    """
    if "sequence" not in records or not records.sequence.map(lambda s: isinstance(s, str)).all():
        raise ValueError("Every record requires a sequence string")
    if not math.isfinite(min_af) or not (0.0 <= min_af <= 1.0):
        raise ValueError("min_af must be finite and within [0.0, 1.0]")
    if type(flank) is not int or flank < 1:
        raise ValueError("flank must be a positive integer")

    qualified = [v for v in variants if v.af >= min_af]
    base = records.copy()
    if "is_coding_variant" not in base.columns:
        base["is_coding_variant"] = False
    for col, default in (
        ("variant_id", ""),
        ("variant_af", pd.NA),
        ("variant_source_id", ""),
        ("variant_position", pd.NA),
    ):
        if col not in base.columns:
            base[col] = default

    if not qualified:
        base.attrs.update(
            records.attrs,
            coding_variant_min_af=min_af,
            coding_variants_applied=0,
        )
        return base

    key_to_indices: dict[str, list[int]] = {}
    for idx, row in enumerate(base.to_dict("records")):
        keys: set[str] = set()
        for col in ("transcript_id", "transcript_accession", "record_id", "gene_id"):
            val = row.get(col)
            if val is not None and pd.notna(val) and str(val).strip():
                raw = str(val).strip()
                keys.add(raw)
                keys.add(raw.split(".", 1)[0])
        for key in keys:
            key_to_indices.setdefault(key, []).append(idx)

    base_records = base.to_dict("records")
    extra_rows: list[dict] = []
    applied = 0
    for variant in qualified:
        matched_indices = key_to_indices.get(variant.target_id) or key_to_indices.get(
            variant.target_id.split(".", 1)[0]
        )
        if not matched_indices:
            if strict_targets:
                raise ValueError(
                    f"CodingVariant {variant.variant_id!r} target_id={variant.target_id!r} "
                    "matched no background sequence record"
                )
            continue
        applied += 1
        pos0 = variant.position - 1
        for row_idx in matched_indices:
            orig = base_records[row_idx]
            seq = orig["sequence"]
            if pos0 >= len(seq):
                raise ValueError(
                    f"CodingVariant {variant.variant_id!r} position {variant.position} "
                    f"exceeds sequence length {len(seq)} on {variant.target_id!r}"
                )
            observed_aa = seq[pos0]
            if observed_aa != variant.ref_aa:
                raise ValueError(
                    f"CodingVariant {variant.variant_id!r} reference mismatch at 1-based "
                    f"position {variant.position} on {variant.target_id!r}: "
                    f"expected {variant.ref_aa!r}, observed {observed_aa!r}"
                )
            mutated = seq[:pos0] + variant.alt_aa + seq[pos0 + 1 :]
            win_start = max(0, pos0 - flank)
            win_end = min(len(mutated), pos0 + flank + 1)
            new_row = dict(orig)
            new_row["sequence"] = mutated[win_start:win_end]
            new_row["is_coding_variant"] = True
            new_row["variant_id"] = variant.variant_id
            new_row["variant_af"] = float(variant.af)
            new_row["variant_source_id"] = variant.source_id
            new_row["variant_position"] = int(variant.position)
            extra_rows.append(new_row)

    out = (
        pd.concat([base, pd.DataFrame(extra_rows, columns=base.columns)], ignore_index=True)
        if extra_rows
        else base
    )
    out.attrs.update(
        records.attrs,
        coding_variant_min_af=min_af,
        coding_variants_applied=applied,
    )
    return out


@dataclass(frozen=True)
class MatchPolicy:
    """Query-to-background Levenshtein costs; exact subtraction by default.

    Set max_distance=1 for ordinary edit distance, or supply positive substitution
    costs for weighted distance. Unspecified substitutions cost 1; insertions and
    deletions have separate positive costs. Equivalence groups (e.g. ('ILV',))
    alone have zero cost. Custom costs are directed (query residue, background
    residue); provide both directions for symmetric behavior. No default reduced
    alphabet or biologically calibrated weighting is implied.
    """

    max_distance: float = 0.0
    equivalence_groups: tuple[str, ...] = ()
    substitution_costs: dict[tuple[str, str], float] = field(default_factory=dict)
    insertion_cost: float = 1.0
    deletion_cost: float = 1.0

    def __post_init__(self):
        if not math.isfinite(self.max_distance) or self.max_distance < 0:
            raise ValueError("max_distance must be finite and nonnegative")
        seen = set()
        for group in self.equivalence_groups:
            if len(group) < 2 or len(set(group)) != len(group) or not set(group) <= AMINO_ACIDS:
                raise ValueError("Equivalence groups require distinct canonical amino acids")
            if seen & set(group):
                raise ValueError("Equivalence groups must be disjoint")
            seen.update(group)
        for pair in self.substitution_costs:
            if len(pair) != 2 or not set(pair) <= AMINO_ACIDS or pair[0] == pair[1]:
                raise ValueError("Substitution keys must be distinct canonical residue pairs")
        for cost in [self.insertion_cost, self.deletion_cost, *self.substitution_costs.values()]:
            if not math.isfinite(cost) or cost <= 0:
                raise ValueError(
                    "Edit costs must be finite and positive; use groups for equivalence"
                )
        object.__setattr__(self, "equivalence_groups", tuple(self.equivalence_groups))
        object.__setattr__(
            self, "substitution_costs", MappingProxyType(dict(self.substitution_costs))
        )

    def normalize(self, sequence):
        return sequence.translate(
            str.maketrans({c: g[0] for g in self.equivalence_groups for c in g})
        )

    def as_dict(self):
        """Detached JSON-serializable parameters for result provenance."""
        return {
            "max_distance": self.max_distance,
            "equivalence_groups": list(self.equivalence_groups),
            "insertion_cost": self.insertion_cost,
            "deletion_cost": self.deletion_cost,
            "substitution_costs": [
                {"query": a, "background": b, "cost": c}
                for (a, b), c in sorted(self.substitution_costs.items())
            ],
        }

    @property
    def edit_bound(self):
        return math.floor(
            (self.max_distance + 1e-12)
            / min(1.0, self.insertion_cost, self.deletion_cost, *self.substitution_costs.values())
        )


def edit_distance(query: str, background: str, *, policy: MatchPolicy | None = None) -> float:
    """Weighted global Levenshtein distance; all residues must be canonical."""
    policy = policy or MatchPolicy()
    if not set(query + background) <= AMINO_ACIDS:
        raise ValueError("Peptides must contain only canonical uppercase amino acids")
    nq, nb = policy.normalize(query), policy.normalize(background)
    previous = [j * policy.insertion_cost for j in range(len(background) + 1)]
    for i, aa in enumerate(query, 1):
        row = [i * policy.deletion_cost]
        for j, bb in enumerate(background, 1):
            substitution = (
                0 if nq[i - 1] == nb[j - 1] else policy.substitution_costs.get((aa, bb), 1.0)
            )
            row.append(
                min(
                    previous[j] + policy.deletion_cost,
                    row[-1] + policy.insertion_cost,
                    previous[j - 1] + substitution,
                )
            )
        previous = row
    return previous[-1]


def _lengths(lengths):
    result = tuple(sorted(set(lengths)))
    if not result or any(type(k) is not int or k < 1 for k in result):
        raise ValueError("lengths must contain positive integers")
    return result


def kmers(sequence: str, lengths=range(8, 16)):
    """Yield (peptide, start, end), including repeat occurrences, within valid runs."""
    lengths = _lengths(lengths)
    for match in SEGMENTS.finditer(sequence):
        for k in lengths:
            for start in range(match.start(), match.end() - k + 1):
                yield sequence[start : start + k], start, start + k


class PeptideIndex:
    """Reusable sequence search retaining every transcript/gene occurrence.

    Accept a sequence table (e.g. antigen_sets.sequences('Unsafe')). Uses a compact
    concatenated sequence search and exact seeds for bounded edit distance rather
    than materializing every 8–15mer. High distance/low costs can require exhaustive
    windows; restrict parameters and batch by protein for large screens.
    """

    def __init__(self, records: pd.DataFrame, *, policy=None, lengths=range(8, 16)):
        self.policy = policy or MatchPolicy()
        self.lengths = _lengths(lengths)
        if (
            "sequence" not in records
            or not records.sequence.map(lambda s: isinstance(s, str)).all()
        ):
            raise ValueError("Every record requires a sequence string")
        self._records = records.copy()
        self.metadata = dict(records.attrs)
        self._non_seq_columns = [column for column in records.columns if column != "sequence"]
        self._col_values = {col: records[col].tolist() for col in self._non_seq_columns}
        self._result_columns = list(
            dict.fromkeys(
                [*self._non_seq_columns, "query", "matched_peptide", "start", "end", "distance"]
            )
        )
        self._segments, self._starts, chunks = [], [], []
        offset = 0
        # Group identical sequences in a single pass while retaining every occurrence (#789).
        grouped_indices: dict[str, list[int]] = {}
        for row_idx, sequence in enumerate(records["sequence"].tolist()):
            grouped_indices.setdefault(sequence, []).append(row_idx)
        for sequence, row_indices in grouped_indices.items():
            for segment in SEGMENTS.finditer(sequence):
                self._starts.append(offset)
                self._segments.append((segment.group(), segment.start(), row_indices))
                chunk = self.policy.normalize(segment.group()) + "|"
                chunks.append(chunk)
                offset += len(chunk)
        self._text = "".join(chunks)

    @classmethod
    def from_sequence_bundle(
        cls,
        bundle=None,
        *,
        species: str | int | None = None,
        background: str = "Unsafe|Other",
        host=None,
        roles: tuple[str, ...] | None = None,
        policy: MatchPolicy | None = None,
        lengths=range(8, 16),
        path=None,
        cache_dir=None,
        download: bool = False,
    ) -> PeptideIndex:
        """Construct a ``PeptideIndex`` directly from offline sequence assets without ``pyensembl``.

        When ``roles`` is omitted, builds the catalog-aligned ``antigen_sets.forbidden_background``
        (``Unsafe``, ``Other``, or ``Unsafe|Other``, plus any ``HostContext`` exclusions) over the
        verified human ``sequence_bundle`` or packaged species proteome occurrences (``dog``,
        ``mouse``) and verifies zero missing gene translations. When ``roles`` is supplied (e.g.
        ``roles=('exclude', 'review')``), indexes every ``SequenceBundle`` record matching those
        roles directly, including alternate-locus Ensembl records outside the primary gene catalog.
        """
        from . import antigen_sets, sequence_bundle
        from .species import resolve_species
        from .species_target_selection import packaged_species_metadata, packaged_species_sequences

        if bundle is not None or path is not None:
            loaded_bundle = bundle if bundle is not None else sequence_bundle.load(path)
            if roles is not None:
                if host is not None:
                    raise ValueError(
                        "HostContext filtering requires catalog background mode (roles=None)"
                    )
                records = loaded_bundle.to_sequences_table(roles=roles)
                return cls(records, policy=policy, lengths=lengths)
            table = loaded_bundle.to_sequences_table()
            inferred_species = (
                resolve_species(species)
                if species is not None
                else resolve_species(int(table.attrs.get("taxon_id", 9606)))
            )
            seq_source = "sequence_bundle" if inferred_species.taxon_id == 9606 else "auto"
        else:
            inferred_species = resolve_species("human" if species is None else species)
            meta = packaged_species_metadata(inferred_species)
            if meta is not None and meta.get("sequence_provider") != "sequence_bundle":
                if roles is not None:
                    raise ValueError(
                        f"roles={roles!r} requires a SequenceBundle (human or explicit bundle=/path=), "
                        f"not {inferred_species.name!r} packaged proteome"
                    )
                table = packaged_species_sequences(inferred_species, allow_download=download)
                seq_source = "auto"
            else:
                loaded_bundle = sequence_bundle.load(
                    sequence_bundle.ensure(cache_dir=cache_dir, download=download)
                )
                if roles is not None:
                    if host is not None:
                        raise ValueError(
                            "HostContext filtering requires catalog background mode (roles=None)"
                        )
                    records = loaded_bundle.to_sequences_table(roles=roles)
                    return cls(records, policy=policy, lengths=lengths)
                table = loaded_bundle.to_sequences_table()
                seq_source = "sequence_bundle"

        records = antigen_sets.forbidden_background(
            background,
            host=host,
            species=inferred_species.name,
            sequence_source=seq_source,
            sequences_table=table,
        )
        if records.attrs.get("missing_gene_ids"):
            raise ValueError(
                "Background has missing gene translations; subtraction coverage is incomplete"
            )
        return cls(records, policy=policy, lengths=lengths)

    def matches(self, peptide: str) -> pd.DataFrame:
        """All background occurrences within the inclusive configured distance.

        Background windows use this index's lengths; indels can match a different
        length. Query may be any nonempty canonical peptide. Empty results retain
        the result schema and reference/policy metadata.
        """
        if not peptide or not set(peptide) <= AMINO_ACIDS:
            raise ValueError("Query must be a nonempty canonical peptide")
        normalized = self.policy.normalize(peptide)
        d = self.policy.edit_bound
        candidates = set()
        if d < len(peptide):
            # Any alignment with <=d edits leaves at least one of d+1 disjoint
            # query blocks unchanged, including when edits are insertions/deletions.
            blocks = [
                (i * len(peptide) // (d + 1), (i + 1) * len(peptide) // (d + 1))
                for i in range(d + 1)
            ]
            for left, right in blocks:
                seed = normalized[left:right]
                pos = self._text.find(seed)
                while pos >= 0:
                    segment_id = bisect.bisect_right(self._starts, pos) - 1
                    local = pos - self._starts[segment_id]
                    for start in range(max(0, local - left - d), local - left + d + 1):
                        candidates.add((segment_id, start))
                    pos = self._text.find(seed, pos + 1)
        else:
            candidates = {
                (i, start) for i, (s, _, _) in enumerate(self._segments) for start in range(len(s))
            }
        rows = []
        for i, start in sorted(candidates):
            sequence, offset, row_indices = self._segments[i]
            for k in self.lengths:
                lower_bound = (
                    (k - len(peptide)) * self.policy.insertion_cost
                    if k >= len(peptide)
                    else (len(peptide) - k) * self.policy.deletion_cost
                )
                if start + k > len(sequence) or lower_bound > self.policy.max_distance + 1e-12:
                    continue
                window = sequence[start : start + k]
                distance = edit_distance(peptide, window, policy=self.policy)
                if distance <= self.policy.max_distance + 1e-12:
                    for row_idx in row_indices:
                        rows.append(
                            {
                                **{
                                    col: self._col_values[col][row_idx]
                                    for col in self._non_seq_columns
                                },
                                "query": peptide,
                                "matched_peptide": window,
                                "start": offset + start,
                                "end": offset + start + k,
                                "distance": distance,
                            }
                        )
        result = pd.DataFrame(rows, columns=self._result_columns)
        result.attrs.update(
            self.metadata, match_policy=self.policy.as_dict(), background_lengths=self.lengths
        )
        return result

    def subtract(self, sequence: str, *, lengths=range(8, 16)) -> pd.DataFrame:
        """Every query occurrence with no configured background match.

        A retained peptide only passes this sequence screen. It is not certified
        targetable. Query coordinates and repeated occurrences are preserved.
        """
        cache, rows = {}, []
        for peptide, start, end in kmers(sequence, lengths):
            if peptide not in cache:
                cache[peptide] = self.matches(peptide).empty
            if cache[peptide]:
                rows.append({"peptide": peptide, "start": start, "end": end})
        result = pd.DataFrame(rows, columns=["peptide", "start", "end"])
        result.attrs.update(
            self.metadata, match_policy=self.policy.as_dict(), background_lengths=self.lengths
        )
        return result

    def with_coding_variants(
        self,
        variants: Iterable[CodingVariant],
        *,
        min_af: float = 0.001,
        flank: int = 14,
        strict_targets: bool = True,
    ) -> PeptideIndex:
        """Return a new ``PeptideIndex`` augmented with population missense variants (#736)."""
        augmented = apply_coding_variants(
            self._records,
            variants,
            min_af=min_af,
            flank=flank,
            strict_targets=strict_targets,
        )
        return PeptideIndex(augmented, policy=self.policy, lengths=self.lengths)
