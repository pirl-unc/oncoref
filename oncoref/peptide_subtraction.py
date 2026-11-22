"""Exact and bounded weighted edit-distance subtraction against protein sets.

Distances are screening parameters, not probabilities of TCR cross-reactivity.
Coordinates are zero-based, half-open in the original translation. Invalid
residues and stops break segments; matches never cross them or protein records.
"""

from __future__ import annotations

import bisect
import math
import re
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
            "Germline variant and sample-specific isoform personalization",
            "Therapeutic risk-benefit adjudication",
        ],
    }


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
        self.metadata = dict(records.attrs)
        self._result_columns = list(
            dict.fromkeys(
                [column for column in records.columns if column != "sequence"]
                + ["query", "matched_peptide", "start", "end", "distance"]
            )
        )
        self._segments, self._starts, chunks = [], [], []
        offset = 0
        # Group identical sequences once while retaining every occurrence.
        for sequence, occurrences in records.groupby("sequence", sort=False):
            for segment in SEGMENTS.finditer(sequence):
                self._starts.append(offset)
                self._segments.append(
                    (segment.group(), segment.start(), occurrences.to_dict("records"))
                )
                chunk = self.policy.normalize(segment.group()) + "|"
                chunks.append(chunk)
                offset += len(chunk)
        self._text = "".join(chunks)

    @classmethod
    def from_sequence_bundle(
        cls,
        bundle=None,
        *,
        background: str = "Unsafe|Other",
        host=None,
        roles: tuple[str, ...] | None = None,
        policy: MatchPolicy | None = None,
        lengths=range(8, 16),
        path=None,
        cache_dir=None,
        download: bool = False,
    ) -> PeptideIndex:
        """Construct a ``PeptideIndex`` directly from ``sequence_bundle`` without ``pyensembl``.

        When ``roles`` is omitted, builds the catalog-aligned ``antigen_sets.forbidden_background``
        (``Unsafe``, ``Other``, or ``Unsafe|Other``, plus any ``HostContext`` exclusions) over the
        verified bundle occurrences and verifies zero missing gene translations. When ``roles`` is
        supplied (e.g. ``roles=('exclude', 'review')``), indexes every bundle record matching those
        roles directly, including alternate-locus Ensembl records outside the primary gene catalog.
        """
        from . import antigen_sets, sequence_bundle

        loaded_bundle = (
            bundle
            if bundle is not None
            else sequence_bundle.load(path)
            if path is not None
            else sequence_bundle.load(
                sequence_bundle.ensure(cache_dir=cache_dir, download=download)
            )
        )
        if roles is not None:
            if host is not None:
                raise ValueError(
                    "HostContext filtering requires catalog background mode (roles=None)"
                )
            records = loaded_bundle.to_sequences_table(roles=roles)
            return cls(records, policy=policy, lengths=lengths)

        table = loaded_bundle.to_sequences_table()
        records = antigen_sets.forbidden_background(
            background,
            host=host,
            species="human",
            sequence_source="sequence_bundle",
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
            sequence, offset, occurrences = self._segments[i]
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
                    for occurrence in occurrences:
                        rows.append(
                            {
                                **{
                                    key: value
                                    for key, value in occurrence.items()
                                    if key != "sequence"
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
