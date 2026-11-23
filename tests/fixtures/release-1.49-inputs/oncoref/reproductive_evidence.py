"""Positive Reproductive evidence with source, lineage and identity constraints.

The HPA v23 tissue-cluster vocabulary is pinned. A generic resident cell supplies
only indirect sampling-site evidence. Its
qualification additionally requires complete outside-cluster coverage, no outside
cluster reaching 5 nTPM, and no curated normal-induction counterevidence.
These mappings attest to positive association; all normal and concentration
gates remain independent. See https://v23.proteinatlas.org/about/download.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .expression_diffuseness import HPA_REPRODUCTIVE_BRANCHES

HPA_REPRODUCTIVE_CELLS = {
    "testis": {
        "leydig cells",
        "spermatocytes",
        "late spermatids",
        "early spermatids",
        "spermatogonia",
        "peritubular cells",
        "sertoli cells",
    },
    "placenta": {"cytotrophoblasts", "extravillous trophoblasts", "syncytiotrophoblasts"},
    "endometrium": {"glandular and luminal cells", "endometrial stromal cells", "ciliated cells"},
    "prostate": {"basal prostatic cells", "prostatic glandular cells"},
    "breast": {"breast glandular cells", "breast myoepithelial cells"},
    "fallopian tube": {"secretory cells", "ciliated cells"},
    "ovary": {"ovarian stromal cells", "granulosa cells", "oocytes"},
}


def bulk_reference(facts: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Use the current full Reproductive panel, never a frozen narrower maximum.

    The reference is a measured same-source protein sum within a tissue, followed
    by a maximum across mapped Reproductive tissues. Missing loci/values stay
    missing. Partial panel coverage can attest to presence but not concentration.
    """
    if (
        scores.entity_id.duplicated().any()
        or set(scores.entity_id) != set(facts.protein_equivalence_group_id)
        or not (
            scores.source.eq("hpa_rna_consensus")
            & scores.version.eq("v23")
            & scores.unit.eq("nTPM")
            & scores.taxonomy_id.eq(9606)
            & scores.summary.eq("consensus")
        ).all()
    ):
        raise ValueError("Require one human HPA v23 consensus reference per protein")
    values = pd.to_numeric(scores.reference_value, errors="raise")
    if values.lt(0).any() or np.isinf(values).any():
        raise ValueError("Require finite nonnegative RNA or missing")
    result = facts.copy()
    result["reproductive_ntpm"] = result.protein_equivalence_group_id.map(
        scores.set_index("entity_id").reference_value
    )
    return result


def cell_observations(contexts: pd.DataFrame) -> pd.DataFrame:
    """Assess native >=5 nTPM, >=20-read reproductive-lineage observations.

    All same-protein loci must be measured in the same cluster before summing;
    the caller supplies complete_members from the identity-aware aggregation.
    Generic vascular, immune, fibroblast and mixed populations are marked as
    indirect candidates, requiring the separate outside-panel and induction checks.
    An unlisted source/tissue/cell combination remains nonqualifying.
    """
    d = contexts.copy()
    if d.value.lt(0).any() or np.isinf(d.value).any():
        raise ValueError("Require finite nonnegative RNA or missing")
    site = d.tissue.isin({t for ts in HPA_REPRODUCTIVE_BRANCHES.values() for t in ts})
    source = (
        d.source.eq("hpa_single_cell_tissue")
        & d.version.eq("v23")
        & d.taxonomy_id.eq(9606)
        & d.unit.eq("nTPM")
        & d.mapping_status.eq("mapped")
        & site
    )
    lineage = pd.Series(
        [
            cell in HPA_REPRODUCTIVE_CELLS.get(tissue, set())
            for tissue, cell in zip(d.tissue, d.cell_type, strict=True)
        ],
        index=d.index,
    )
    complete = d.complete_members.eq(True) & d.member_count.ge(1)
    supported = d["Read count"].ge(20) & d.value.notna()
    d["reproductive_cell_usable"] = source & lineage & complete & supported
    d["reproductive_cell_positive"] = d.reproductive_cell_usable & d.value.ge(5)
    d["reproductive_sampling_site"] = site
    d["reproductive_cell_indirect_candidate"] = (
        source & ~lineage & complete & supported & d.value.ge(5)
    )
    d["reproductive_cell_reason"] = np.select(
        [~source, ~lineage, ~complete, ~supported, d.value.lt(5)],
        [
            "source_or_mapping",
            "generic_or_unreviewed_cell",
            "incomplete_members",
            "insufficient_read_support",
            "below_5_ntpm",
        ],
        default="qualifying",
    )
    return d


def indirect_association(facts: pd.DataFrame) -> pd.Series:
    """Admit localized indirect RNA only with observed outside-cluster coverage.

    This establishes sampling-site association, not the producing cell. Curated
    induction in nonreproductive cells contradicts an indirect-only localization.
    It does not itself establish Unsafe or invalidate direct lineage evidence.
    """
    return (
        facts.reproductive_cell_indirect_ntpm.ge(5)
        & facts.outside_cluster_complete.eq(True)
        & facts.outside_cluster_max_ntpm.lt(5)
        & facts.normal_induction_evidence.eq(False)
    )
