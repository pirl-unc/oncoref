"""Three-way research classification with independent evidence at every layer.

This proposal does not mutate released expression labels or sequence references.
Unsafe denotes a policy exclusion; it does not assert demonstrated toxicity.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
import pandas as pd


class EvidenceLayer(NamedTuple):
    """One independent decision layer; masks align exactly to the facts index."""

    name: str
    candidate_hold: pd.Series
    unsafe_evidence: pd.Series


class FunnelResult(NamedTuple):
    """Entity decisions, sequential counts and overlapping independent effects."""

    classifications: pd.DataFrame
    funnel: pd.DataFrame
    independent_effects: pd.DataFrame


def classify_layers(
    facts: pd.DataFrame,
    layers: Sequence[EvidenceLayer | tuple[str, pd.Series, pd.Series]],
    *,
    entity_column: str = "gene_id",
) -> FunnelResult:
    """Apply (name, candidate_hold, unsafe_evidence) to the entire gene universe.

    Later evidence may promote an earlier Neither to Unsafe. Missing evidence
    must be supplied as a hold, never an Unsafe flag. At intermediate stages
    CTA means eligible so far, not a completed assessment. Input masks must
    be explicit, complete booleans aligned to facts. ``entity_column`` identifies
    either gene loci or protein equivalence groups; the two counting units must
    never be mixed. Returns copied facts with decisions and two count tables.
    """
    if (
        entity_column not in facts
        or facts[entity_column].isna().any()
        or not facts[entity_column].is_unique
    ):
        raise ValueError("Require unique entity identities")
    if not layers:
        raise ValueError("At least one evidence layer is required")
    names = [name for name, _, _ in layers]
    if len(set(names)) != len(names):
        raise ValueError("Layer names must be unique")
    holds, unsafe = {}, {}
    for name, hold, strong in layers:
        for mask in [hold, strong]:
            if (
                not mask.index.equals(facts.index)
                or mask.isna().any()
                or not pd.api.types.is_bool_dtype(mask)
            ):
                raise ValueError("Complete aligned boolean evidence required")
        holds[name], unsafe[name] = hold | strong, strong
    hold_table, unsafe_table = pd.DataFrame(holds), pd.DataFrame(unsafe)
    held = pd.Series(False, index=facts.index)
    excluded = held.copy()
    sequential, independent = [], []
    all_holds, all_unsafe = hold_table.any(axis=1), unsafe_table.any(axis=1)
    for name in names:
        before = ~(held | excluded)
        old_unsafe = excluded.copy()
        held |= hold_table[name]
        excluded |= unsafe_table[name]
        sequential.append(
            {
                "layer": name,
                "new_candidate_removals": int((before & (held | excluded)).sum()),
                "new_Unsafe": int((excluded & ~old_unsafe).sum()),
                "CTA_eligible_so_far": int((~(held | excluded)).sum()),
                "Unsafe": int(excluded.sum()),
                "Neither": int((held & ~excluded).sum()),
            }
        )
        other_holds = hold_table.drop(columns=name).any(axis=1)
        other_unsafe = unsafe_table.drop(columns=name).any(axis=1)
        independent.append(
            {
                "layer": name,
                "candidate_exclusions_anywhere": int(hold_table[name].sum()),
                "Unsafe_evidence_anywhere": int(unsafe_table[name].sum()),
                "CTA_added_if_layer_removed": int((all_holds & ~other_holds).sum()),
                "Unsafe_reclassified_if_layer_removed": int((all_unsafe & ~other_unsafe).sum()),
            }
        )
    # Replays may receive previous decision columns. Only the supplied layers
    # define this invocation's holds and exclusions.
    result = facts.drop(columns=[c for c in facts if c.startswith(("hold:", "unsafe:"))]).copy()
    result["category"] = np.select([all_unsafe, all_holds], ["Unsafe", "Neither"], default="CTA")
    result["exclusion_reasons"] = hold_table.apply(lambda row: " | ".join(row.index[row]), axis=1)
    result["Unsafe_reasons"] = unsafe_table.apply(lambda row: " | ".join(row.index[row]), axis=1)
    for name in names:
        result[f"hold:{name}"] = hold_table[name]
        result[f"unsafe:{name}"] = unsafe_table[name]
    return FunnelResult(result, pd.DataFrame(sequential), pd.DataFrame(independent))
