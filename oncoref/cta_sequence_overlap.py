# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Sequence overlap between candidate proteins and the Unsafe blacklist.

Unsafe proteins blacklist any matching subsequence from CTAs. A would-be CTA
whose 9-mers lie almost entirely inside Unsafe proteins therefore keeps no
usable peptide. Readthrough transcripts and young paralogs are typical: their
RNA and sequence belong to a neighbouring gene that is Unsafe.

The 9-mers of a protein group are the union over its member genes of every
protein-coding transcript's translated sequence in the pinned Ensembl release.
A group with no 9-mer (all sequences shorter than nine residues) has an
undefined overlap and is never contained.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd

K = 9
ENSEMBL_RELEASE = 115


def gene_ninemers(gene_ids: Sequence[str], release: int = ENSEMBL_RELEASE) -> dict[str, frozenset]:
    """9-mer set per gene across all protein-coding transcripts of ``release``."""
    from pyensembl import EnsemblRelease

    genome = EnsemblRelease(release)
    result = {}
    for gene_id in dict.fromkeys(gene_ids):
        sequences = set()
        try:
            for transcript in genome.gene_by_id(gene_id).transcripts:
                if transcript.is_protein_coding and transcript.protein_sequence:
                    sequences.add(transcript.protein_sequence)
        except ValueError:  # gene id absent from this release
            pass
        result[gene_id] = frozenset(
            seq[i : i + K] for seq in sequences for i in range(len(seq) - K + 1)
        )
    return result


def group_ninemers(member_gene_ids: pd.Series, release: int = ENSEMBL_RELEASE) -> list[frozenset]:
    """9-mer set for each protein group, given JSON lists of member gene IDs."""
    members = [json.loads(ids) for ids in member_gene_ids]
    per_gene = gene_ninemers([g for ids in members for g in ids], release)
    return [frozenset().union(*(per_gene[g] for g in ids)) for ids in members]


def unsafe_overlap(
    ninemers: Sequence[frozenset], index: pd.Index
) -> Callable[[pd.Series], pd.Series]:
    """Return ``mask -> fraction of each group's 9-mers present in any masked group``.

    The result is aligned to ``index`` and NaN where a group has no 9-mer. A group
    does not count toward its own overlap unless it is itself in the mask.
    """
    sizes = np.array([len(kmers) for kmers in ninemers], dtype=float)
    cache: dict[bytes, pd.Series] = {}

    def overlap(unsafe: pd.Series) -> pd.Series:
        key = unsafe.to_numpy(bool).tobytes()
        if key not in cache:
            blacklist = frozenset().union(*(k for k, flag in zip(ninemers, unsafe) if flag))
            shared = np.array([len(k & blacklist) for k in ninemers], dtype=float)
            with np.errstate(invalid="ignore", divide="ignore"):
                fraction = np.where(sizes > 0, shared / sizes, np.nan)
            cache[key] = pd.Series(fraction, index=index)
        return cache[key].copy()

    return overlap
