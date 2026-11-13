# Licensed under the Apache License, Version 2.0.
"""Optional tightening of the Target candidate catalog, one named knob at a time.

The catalog in :mod:`oncoref.target_selection` is the *frame*: IHC-dependent bulk
RNA ceilings, single-cell rules, coverage, identity, sequence and expression bounds
have already run. Nothing here changes a catalog category. ``TargetFilter`` follows
the catalog's own shape and can only remove proteins:

1. **off_target**  lower the ceiling on expression outside tumors and reproductive tissue
2. **warnings**    act on flags the catalog records but does not exclude on
3. **curation**    act on cited literature: normal function by tissue tier, read-through loci
4. **sets**        raise the evidence each antigen set needs, then pick a set
5. **papers**      require published support
6. **exclusions**  drop named gene families (a judgment, not a measurement)

A protein is kept only when every layer passes. ``filter_failed`` names the first
layer that removed it. Two conventions hold for every knob: ``max_*`` means
*below* the value and ``min_*`` means *at least* the value; a missing measurement
never passes. Filter modes (``FILTER_MODES``) are starting points, not policy.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache

import pandas as pd

FILTER_POLICY_ID = "oncoref.target-filters.2026-11-12.4"

#: Which set a protein must belong to, after the evidence knobs have tightened both.
MEMBERSHIPS = {
    "target_candidates": lambda ca, ra: ca | ra,
    "cancer_antigen": lambda ca, ra: ca,
    "reproductive_antigen": lambda ca, ra: ra,
    "cta": lambda ca, ra: ca & ra,
    "cancer_antigen_only": lambda ca, ra: ca & ~ra,
    "reproductive_antigen_only": lambda ca, ra: ra & ~ca,
}

#: Paper families, from the source-level roles recorded by ``cta_literature``.
PAPER_FAMILIES = {
    "cta": "catalog_nomination",  # cancer/testis, cancer-germline, cancer-placenta lists
    "oncofetal": "cancer_embryo_program",  # embryonic gene programs re-expressed in cancer
    "tumor": "tumor_expression",  # single-gene tumor expression reports
    "reproductive": "reproductive_developmental",  # normal reproductive enrichment only
}
#: Literature verdicts that ``literature="no_conflict"`` removes.
LITERATURE_CONFLICTS = ("conflicting_or_nonspecific", "specialized_tissue_expression")
LAYERS = ("off_target", "warnings", "curation", "sets", "papers", "exclusions")
#: ``normal_function`` values and the curated tiers each one removes.
NORMAL_FUNCTION = {
    "any": (),
    "no_vital": ("vital",),
    "no_protected": ("vital", "protected"),
    "no_manageable": ("vital", "protected", "manageable"),
}
_REPRODUCTIVE_RNA = ("reproductive_ntpm", "reproductive_cell_ntpm", "developmental_mean_tpm")


@dataclass(frozen=True)
class TargetFilter:
    """Every knob, with the loosest value as its default.

    ``max_*`` knobs keep values strictly below the knob; ``min_*`` knobs keep values
    at or above it. TPM and HPA nTPM are separate scales and are never converted.
    """

    # 1. off target: expression outside tumors and reproductive tissue -------------
    #: Tissue level, native TPM: max(normal reference N, Disposable) must be below this.
    max_off_target_tpm: float | None = None
    #: HPA nTPM, which also takes in supported single-cell clusters, must be below this.
    max_off_target_ntpm: float | None = None
    #: Disposable tissue alone, on both scales (the catalog already caps it at 32).
    max_disposable: float | None = None
    # 2. warnings: recorded by the catalog, not acted on by it ----------------------
    rare_cell: str = "any"  # "any" | "no_dismissed"
    literature: str = "any"  # "any" | "no_conflict"
    mapping: str = "any"  # "any" | "no_risk"
    # 3. curation: cited literature, see oncoref.target_curation ---------------------
    #: Drop proteins with a curated normal function in a tissue at or above this tier.
    normal_function: str = "any"  # "any" | "no_vital" | "no_protected" | "no_manageable"
    readthrough: str = "any"  # "any" | "exclude": read-through loci of other genes
    # 4. sets: tighten each set's evidence, then choose a set ----------------------
    min_patients: int = 0  # Cancer Antigen: independent tumors above the normal reference
    min_tumor_fold: float = 0.0  # Cancer Antigen: highest tumor / max(normal reference, 1)
    min_reproductive: float = 0.0  # Reproductive Antigen: best reproductive RNA route
    membership: str = "target_candidates"
    # 5. papers: distinct DOIs; a DOI counts once however many tables it has -------
    min_papers: int = 0
    families: tuple = ("cta", "oncofetal", "tumor")
    min_families: int = 0
    exact_members_only: bool = False  # drop nominations naming only some member genes
    # 6. exclusions: gene families dropped by name ---------------------------------
    exclude_prefixes: tuple = ()
    except_labels: tuple = ()

    def __post_init__(self):
        if self.membership not in MEMBERSHIPS:
            raise ValueError(f"membership must be one of {sorted(MEMBERSHIPS)}")
        if unknown := set(self.families) - set(PAPER_FAMILIES):
            raise ValueError(f"unknown paper families {sorted(unknown)}")
        if isinstance(self.exclude_prefixes, str) or isinstance(self.except_labels, str):
            raise ValueError("exclude_prefixes and except_labels must be tuples, not strings")
        if self.rare_cell not in ("any", "no_dismissed"):
            raise ValueError("rare_cell must be 'any' or 'no_dismissed'")
        if self.literature not in ("any", "no_conflict"):
            raise ValueError("literature must be 'any' or 'no_conflict'")
        if self.mapping not in ("any", "no_risk"):
            raise ValueError("mapping must be 'any' or 'no_risk'")
        if self.normal_function not in NORMAL_FUNCTION:
            raise ValueError(f"normal_function must be one of {list(NORMAL_FUNCTION)}")
        if self.readthrough not in ("any", "exclude"):
            raise ValueError("readthrough must be 'any' or 'exclude'")
        if self.min_families > len(self.families):
            raise ValueError("min_families exceeds the number of counted families")

    def changed(self) -> dict:
        """Knobs that differ from the loosest setting; empty means the plain catalog."""
        base = asdict(TargetFilter())
        return {k: v for k, v in asdict(self).items() if v != base[k]}


#: Proteins every filter mode must keep; a test enforces it. PRAME and NY-ESO-1 are
#: targets of T-cell therapies given to patients, so a mode that drops one is miscalibrated.
POSITIVE_CONTROLS = ("PRAME", "NY-ESO-1", "MAGEA4")

#: Named filter modes, loosest to strictest. Each keeps a subset of the one before it.
FILTER_MODES = {
    "all": TargetFilter(),
    "loose": TargetFilter(readthrough="exclude", min_papers=1),
    "standard": TargetFilter(
        max_off_target_tpm=5.0,
        normal_function="no_vital",
        readthrough="exclude",
        min_patients=2,
        min_papers=2,
    ),
    "strict": TargetFilter(
        max_off_target_tpm=4.0,
        normal_function="no_protected",
        readthrough="exclude",
        min_patients=5,
        min_tumor_fold=10.0,
        membership="cta",
        min_papers=3,
        exact_members_only=True,
        exclude_prefixes=("MAGE",),
        except_labels=("MAGEA4",),
    ),
}


@lru_cache(maxsize=1)
def paper_table() -> pd.DataFrame:
    """One row per protein, family and paper that counts for that protein.

    A CTA-list paper counts only for proteins it nominates: appearing in such a
    paper's background table is not a nomination. Oncofetal and tumor papers have
    no nomination tables, so membership of the program or the report itself counts.
    """
    from .cta_literature import attestations

    links = attestations()
    rows = []
    columns = ["entity_id", "source_id", "source_roles", "published_nomination"]
    for record in links[columns].itertuples(index=False):
        for family, role in PAPER_FAMILIES.items():
            if role not in record.source_roles:
                continue
            if family == "cta" and not record.published_nomination:
                continue
            rows.append((record.entity_id, family, record.source_id))
    return pd.DataFrame(rows, columns=["entity_id", "family", "source_id"]).drop_duplicates()


def paper_counts(families=tuple(PAPER_FAMILIES)) -> pd.DataFrame:
    """One row per entity: distinct papers in each family and across ``families``."""
    table = paper_table()
    counts = table.groupby(["entity_id", "family"]).source_id.nunique().unstack(fill_value=0)
    counts = counts.reindex(columns=list(PAPER_FAMILIES), fill_value=0)
    counts = counts.add_prefix("papers_")
    chosen = table[table.family.isin(families)]
    counts["papers_counted"] = chosen.groupby("entity_id").source_id.nunique()
    counts["families_counted"] = chosen.groupby("entity_id").family.nunique()
    return counts.fillna(0).astype(int)


def _frame(table: pd.DataFrame | None) -> pd.DataFrame:
    if table is not None:
        return table.copy()
    from . import target_selection

    return target_selection.results()


def _reference_context(value) -> str:
    """Tissue behind the catalog's normal reference, from its stored JSON receipt."""
    try:
        records = json.loads(value) if isinstance(value, str) else value
        return str(records[0].get("tissue", "")) if records else ""
    except (TypeError, ValueError, AttributeError, IndexError, KeyError):
        return ""


def measures(table: pd.DataFrame | None = None) -> pd.DataFrame:
    """The catalog plus every quantity a knob reads, before any knob is applied.

    ``off_target_tpm``   max(normal reference N, Disposable), native TPM. N is the
                         catalog's own reference: bulk tissues, GTEx donor medians,
                         brain subregions and immune-cell summaries.
    ``off_target_tpm_tissue`` where that TPM value was measured ("" when it is 0).
    ``off_target_ntpm``  max(healthy maximum, Disposable), HPA nTPM. The healthy
                         maximum also takes in supported single-cell clusters.
    ``tumor_fold``       highest tumor / max(N, 1), the catalog's one-tumor ratio.
    ``reproductive_best`` highest of adult reproductive nTPM, reproductive cell
                         nTPM and developmental mean TPM, each in its own unit.
    """
    d = _frame(table)
    derived = [f"papers_{family}" for family in PAPER_FAMILIES]
    d = d.drop(columns=[*derived, "papers_counted", "families_counted"], errors="ignore")
    for column in ("is_cancer_antigen", "is_reproductive_antigen", "is_cta", "is_target_candidate"):
        d[column] = d[column].astype(bool)
    d["off_target_tpm"] = d[["normal_reference_tpm", "disposable_max_tpm"]].max(
        axis=1, skipna=False
    )
    d["off_target_ntpm"] = d[["healthy_maximum_ntpm", "disposable_max_ntpm"]].max(
        axis=1, skipna=False
    )
    tissue = d.normal_reference_contexts.map(_reference_context)
    tissue = tissue.where(d.normal_reference_tpm.ge(d.disposable_max_tpm), "Disposable")
    d["off_target_tpm_tissue"] = tissue.where(d.off_target_tpm.gt(0), "")
    d["tumor_fold"] = d.tumor_max_tpm / d.normal_reference_tpm.clip(lower=1.0)
    routes = d[list(_REPRODUCTIVE_RNA)].apply(pd.to_numeric, errors="coerce")
    d["reproductive_best"] = routes.max(axis=1).fillna(0.0)
    d["mapping_risk"] = d.paralog_ambiguous.eq(True) | d.host_overlapped.eq(True)
    from . import target_curation

    curated = target_curation.table(d.protein_group_label)
    d = d.drop(columns=list(curated.columns), errors="ignore").join(
        curated, on="protein_group_label"
    )
    counts = paper_counts()
    d = d.join(counts.drop(columns=["papers_counted", "families_counted"]), on="entity_id")
    for column in counts.columns.drop(["papers_counted", "families_counted"]):
        d[column] = d[column].fillna(0).astype(int)
    return d


def apply(f: TargetFilter, table: pd.DataFrame | None = None) -> pd.DataFrame:
    """Evaluate every layer for every protein; nothing is dropped.

    Returns :func:`measures` plus ``pass_<layer>`` for each layer, ``keep`` (all
    layers pass), ``filter_failed`` (first failing layer, or "") and
    ``filtered_category`` (the category under the tightened set evidence).
    """
    d = measures(table)

    low = pd.Series(True, index=d.index)
    if f.max_off_target_tpm is not None:
        low &= d.off_target_tpm.lt(f.max_off_target_tpm)
    if f.max_off_target_ntpm is not None:
        low &= d.off_target_ntpm.lt(f.max_off_target_ntpm)
    if f.max_disposable is not None:
        low &= d.disposable_max_tpm.lt(f.max_disposable) & d.disposable_max_ntpm.lt(
            f.max_disposable
        )
    d["pass_off_target"] = low

    clean = pd.Series(True, index=d.index)
    if f.rare_cell == "no_dismissed":
        clean &= d.rare_cell_status.ne("dismissed")
    if f.literature == "no_conflict":
        clean &= ~d.literature_verdict.isin(LITERATURE_CONFLICTS)
    if f.mapping == "no_risk":
        clean &= ~d.mapping_risk
    d["pass_warnings"] = clean

    curated = ~d.curated_tier.isin(NORMAL_FUNCTION[f.normal_function])
    if f.readthrough == "exclude":
        curated &= ~d.is_readthrough
    d["pass_curation"] = curated

    # The evidence knobs tighten each set separately; membership is decided afterwards.
    cancer = d.is_cancer_antigen
    if f.min_patients:
        cancer = cancer & d.independent_patients_above_normal.ge(f.min_patients)
    if f.min_tumor_fold:
        cancer = cancer & d.tumor_fold.ge(f.min_tumor_fold)
    reproductive = d.is_reproductive_antigen
    if f.min_reproductive:
        reproductive = reproductive & d.reproductive_best.ge(f.min_reproductive)
    d["filtered_cancer_antigen"] = cancer
    d["filtered_reproductive_antigen"] = reproductive
    d["filtered_category"] = ""
    d.loc[cancer & reproductive, "filtered_category"] = "CTA"
    d.loc[cancer & ~reproductive, "filtered_category"] = "Cancer Antigen only"
    d.loc[reproductive & ~cancer, "filtered_category"] = "Reproductive Antigen only"
    d["pass_sets"] = MEMBERSHIPS[f.membership](cancer, reproductive).astype(bool)

    counted = paper_counts(f.families)[["papers_counted", "families_counted"]]
    d = d.join(counted, on="entity_id")
    for column in counted:
        d[column] = d[column].fillna(0).astype(int)
    papers = d.papers_counted.ge(f.min_papers) & d.families_counted.ge(f.min_families)
    if f.exact_members_only:
        papers &= d.nomination_scope.ne("some_members")
    d["pass_papers"] = papers

    named = pd.Series(False, index=d.index)
    if f.exclude_prefixes:
        label = d.protein_group_label.astype(str)
        excepted = label.isin(f.except_labels) | label.str.split("/").str[0].isin(f.except_labels)
        named = label.str.startswith(tuple(f.exclude_prefixes)) & ~excepted
    d["pass_exclusions"] = ~named

    passes = d[[f"pass_{layer}" for layer in LAYERS]].astype(bool)
    d["keep"] = passes.all(axis=1)
    first = passes.to_numpy().argmin(axis=1)
    d["filter_failed"] = [("" if kept else LAYERS[i]) for kept, i in zip(d.keep, first)]
    return d


def select(f: TargetFilter | str = "all", table: pd.DataFrame | None = None) -> pd.DataFrame:
    """Proteins that survive every layer. ``f`` is a ``TargetFilter`` or a filter mode name."""
    f = FILTER_MODES[f] if isinstance(f, str) else f
    d = apply(f, table)
    return d[d.keep].reset_index(drop=True)


def funnel(f: TargetFilter | str, table: pd.DataFrame | None = None) -> pd.DataFrame:
    """Proteins remaining after each layer, in order, starting from the candidates."""
    f = FILTER_MODES[f] if isinstance(f, str) else f
    d = apply(f, table)
    alive = d.is_target_candidate.copy()
    rows = [("Target candidates", int(alive.sum()), None)]
    for layer in LAYERS:
        before = int(alive.sum())
        alive &= d[f"pass_{layer}"]
        rows.append((layer, int(alive.sum()), before - int(alive.sum())))
    out = pd.DataFrame(rows, columns=["layer", "remaining", "removed"])
    out.attrs["remaining_by_category"] = d[alive].filtered_category.value_counts().to_dict()
    return out


def describe(f: TargetFilter | str) -> str:
    """Plain-language statement of exactly what ``f`` keeps, one line per active knob."""
    f = FILTER_MODES[f] if isinstance(f, str) else f
    lines = ["Start from the catalog's Target candidates; catalog rules are unchanged."]
    if f.max_off_target_tpm is not None:
        lines.append(
            f"Off target: the highest tissue value outside tumors and reproductive tissue, "
            f"Disposable tissue included, is below {f.max_off_target_tpm:g} TPM."
        )
    if f.max_off_target_ntpm is not None:
        lines.append(
            f"Off target: the highest HPA value outside tumors and reproductive tissue, "
            f"single-cell clusters and Disposable tissue included, is below "
            f"{f.max_off_target_ntpm:g} nTPM."
        )
    if f.max_disposable is not None:
        lines.append(
            f"Off target: Disposable tissue alone is below {f.max_disposable:g} on both scales."
        )
    if f.rare_cell == "no_dismissed":
        lines.append("Warnings: no dismissed rare-cell signal.")
    if f.literature == "no_conflict":
        lines.append("Warnings: literature review is not conflicting or tissue-specialized.")
    if f.mapping == "no_risk":
        lines.append("Warnings: no paralog or host-gene mapping risk.")
    if f.normal_function != "any":
        tiers = " or ".join(t.capitalize() for t in NORMAL_FUNCTION[f.normal_function])
        lines.append(f"Curation: no cited normal function or expression in a {tiers} tissue.")
    if f.readthrough == "exclude":
        lines.append("Curation: not a read-through of other protein-coding genes.")
    if f.min_patients:
        lines.append(
            f"Sets: a Cancer Antigen needs at least {f.min_patients} independent tumors above "
            f"its normal reference."
        )
    if f.min_tumor_fold:
        lines.append(
            f"Sets: a Cancer Antigen's highest tumor is at least {f.min_tumor_fold:g} x "
            f"max(normal reference, 1)."
        )
    if f.min_reproductive:
        lines.append(
            f"Sets: a Reproductive Antigen has a reproductive RNA route of at least "
            f"{f.min_reproductive:g} in its own unit."
        )
    lines.append(
        f"Sets: keep proteins that are then in {f.membership.replace('_', ' ')}."
        if f.membership != "target_candidates"
        else "Sets: keep proteins that are then in either antigen set."
    )
    if f.min_papers or f.min_families:
        text = f"Papers: at least {f.min_papers} distinct paper(s) in {'/'.join(f.families)}"
        if f.min_families:
            text += f", spanning at least {f.min_families} families"
        lines.append(text + ".")
    if f.exact_members_only:
        lines.append("Papers: no nomination that names only some of the protein's member genes.")
    if f.exclude_prefixes:
        text = f"Exclusions: drop labels starting with {'/'.join(f.exclude_prefixes)}"
        if f.except_labels:
            text += f", except {', '.join(f.except_labels)}"
        lines.append(text + ".")
    return "\n".join(lines)
