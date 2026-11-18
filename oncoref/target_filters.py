# Licensed under the Apache License, Version 2.0.
"""Optional tightening of the Target candidate catalog, one named knob at a time.

The catalog in :mod:`oncoref.target_selection` is the *frame*. A Cancer Antigen
already has a tumor sample at 10 x every non-reproductive normal value, and the
normal screens, published normal expression, coverage, identity and sequence have
run. Catalog knobs read that catalog; optional host constraints read the separate
versioned exposure evidence. Neither rewrites biological catalog membership:

1. **off_target**  a lower ceiling on expression outside tumors and reproductive tissue
2. **curation**    the tier of a recorded normal home (literature finding, family, reference)
3. **host**        explicit protection, absence and accepted-loss choices
4. **sets**        enrichment in each compartment, tumor replication, then which set to keep
5. **papers**      published support
6. **exclusions**  gene families dropped by name (a judgment, not a measurement)

A protein is kept only when every layer passes. ``filter_failed`` names the first
layer that removed it. ``max_*`` means *below* the value and ``min_*`` means *at
least* the value; a missing measurement never passes. Filter modes
(``FILTER_MODES``) are starting points, not policy.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from . import host_context, target_selection

FILTER_POLICY_ID = "oncoref.target-filters.2026-11-17.5"

#: Compatibility names describe assumed anatomy, never accepted organ loss.
HOST_PROFILES = {
    "all": (),
    **{
        name: host_context.from_legacy(name).protected
        for name in ("male-patient", "female-patient")
    },
}


def normalize_compartments(values) -> tuple[str, ...]:
    return host_context.normalize(values)


#: Which set a protein must belong to after tightening compartment evidence.
MEMBERSHIPS = {
    "target_candidates": lambda ca, ra: ca | ra,
    "cancer_antigen": lambda ca, ra: ca,
    "reproductive_antigen": lambda ca, ra: ra,
    "cta": lambda ca, ra: ca & ra,
    "cancer_antigen_only": lambda ca, ra: ca & ~ra,
    "reproductive_antigen_only": lambda ca, ra: ra & ~ca,
}

#: Legacy source-role labels; eligibility now consumes per-gene PAPER_USES.
PAPER_FAMILIES = {
    "cta": "catalog_nomination",  # cancer/testis, cancer-germline, cancer-placenta lists
    "oncofetal": "cancer_embryo_program",  # embryonic gene programs re-expressed in cancer
    "tumor": "tumor_expression",  # single-gene tumor expression reports
    "reproductive": "reproductive_developmental",  # normal reproductive enrichment only
}
#: Families whose papers ``min_papers`` counts: target nominations and tumor reports.
COUNTED_FAMILIES = ("cta", "oncofetal", "tumor")
LAYERS = ("off_target", "curation", "host", "sets", "papers", "exclusions")
#: ``normal_function`` values and the curated tiers each one removes.
NORMAL_FUNCTION = {
    "any": (),
    "no_vital": ("vital",),
    "no_protected": ("vital", "protected"),
    "no_manageable": ("vital", "protected", "manageable"),
}


@dataclass(frozen=True)
class TargetFilter:
    """Every knob, with the loosest value as its default."""

    #: 1. Highest tissue value outside tumors and reproductive tissue, Disposable
    #: included, native TPM, must be below this.
    max_off_target_tpm: float | None = None
    #: 2. Drop proteins whose recorded normal home is at or above this tier.
    normal_function: str = "any"  # "any" | "no_vital" | "no_protected" | "no_manageable"
    #: 3. A Cancer Antigen's highest tumor / max(non-reproductive normal, 1). The
    #: catalog already requires 10; a higher value moves weaker ones out of the set.
    min_tumor_fold: float = 10.0
    min_reproductive_fold: float = 10.0
    #: Known independent patients at min_tumor_fold; anonymous samples do not count.
    min_tumor_patients: int = 0
    membership: str = "target_candidates"
    #: Host compartment profile and explicit physiological/surgical overrides.
    host: host_context.HostContext | None = None
    absent_tissues: tuple[str, ...] = ()
    host_profile: str = "all"  # "all" | "male-patient" | "female-patient"
    lactating: bool = False
    pregnant: bool = False
    allow_tissues: tuple[str, ...] = ()
    protect_tissues: tuple[str, ...] = ()
    #: 4. Distinct papers that nominate the protein or report it in tumors.
    min_papers: int = 0
    #: 5. Labels dropped by prefix, and labels kept despite it.
    exclude_prefixes: tuple = ()
    except_labels: tuple = ()

    def __post_init__(self):
        if self.membership not in MEMBERSHIPS:
            raise ValueError(f"membership must be one of {sorted(MEMBERSHIPS)}")
        if self.host_profile == "core":
            host_context.from_legacy("core")
        if self.host_profile not in HOST_PROFILES:
            raise ValueError(f"host_profile must be one of {list(HOST_PROFILES)}")
        for flag in ("lactating", "pregnant"):
            if not isinstance(getattr(self, flag), bool):
                raise ValueError(f"{flag} must be a boolean")
        for field in (
            "exclude_prefixes",
            "except_labels",
            "allow_tissues",
            "protect_tissues",
            "absent_tissues",
        ):
            if not isinstance(getattr(self, field), tuple):
                raise ValueError(
                    f"{field} must be a tuple, not {type(getattr(self, field)).__name__}"
                )
        context(self)
        if self.max_off_target_tpm is not None and (
            isinstance(self.max_off_target_tpm, bool)
            or not isinstance(self.max_off_target_tpm, (int, float))
            or not math.isfinite(self.max_off_target_tpm)
            or self.max_off_target_tpm <= 0
        ):
            raise ValueError("max_off_target_tpm must be finite and positive")
        if self.normal_function not in NORMAL_FUNCTION:
            raise ValueError(f"normal_function must be one of {list(NORMAL_FUNCTION)}")
        for name in ("min_tumor_fold", "min_reproductive_fold"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 10.0:
                raise ValueError(f"{name} must be finite and at least the catalog's 10")
        if (
            isinstance(self.min_tumor_patients, bool)
            or not isinstance(self.min_tumor_patients, int)
            or self.min_tumor_patients < 0
        ):
            raise ValueError("min_tumor_patients must be a nonnegative integer")

    def changed(self) -> dict:
        """Knobs that differ from the loosest setting; empty means the plain catalog."""
        base = asdict(TargetFilter())
        return {k: v for k, v in asdict(self).items() if v != base[k]}


def context(f):
    if isinstance(f, str):
        f = FILTER_MODES[f]
    if f.host is not None:
        if not isinstance(f.host, host_context.HostContext):
            raise ValueError("host must be a HostContext")
        if any(
            (
                f.host_profile != "all",
                f.pregnant,
                f.lactating,
                f.allow_tissues,
                f.protect_tissues,
                f.absent_tissues,
            )
        ):
            raise ValueError("Use host or legacy profile/tissue arguments, not both")
        return f.host
    return host_context.from_legacy(
        f.host_profile,
        pregnant=f.pregnant,
        lactating=f.lactating,
        allow_tissues=f.allow_tissues,
        protect_tissues=f.protect_tissues,
        absent_tissues=f.absent_tissues,
    )


def protected_compartments(f: TargetFilter | str) -> tuple[str, ...]:
    host = context(f)
    return host.protected if host else ()


def allowed_compartments(f: TargetFilter | str) -> tuple[str, ...]:
    """Compatibility view of non-protected compartments; inspect context for states."""
    host = context(f)
    return tuple(
        c for c in target_selection.REPRODUCTIVE_COMPARTMENTS if not host or c not in host.protected
    )


#: Baseline catalog regression controls. Host constraints may exclude any of them.
POSITIVE_CONTROLS = ("PRAME", "NY-ESO-1", "MAGEA4")

#: Named filter modes, loosest to strictest. Each keeps a subset of the one before it.
FILTER_MODES = {
    "all": TargetFilter(),
    "loose": TargetFilter(min_papers=1),
    "standard": TargetFilter(
        max_off_target_tpm=5.0,
        normal_function="no_vital",
        min_papers=2,
    ),
    "strict": TargetFilter(
        max_off_target_tpm=4.0,
        normal_function="no_protected",
        membership="cta",
        min_papers=3,
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
    Exact loci map each record to the current protein group. Rejected assay,
    related-species, fusion-only and other ineligible coverage cannot earn credit.
    """
    from .cta_literature import PAPER_USES, current_attestations

    links = current_attestations()
    links = links[links.eligible_for_coverage]
    rows = []
    columns = ["entity_id", "source_id", "allowed_uses"]
    for record in links[columns].itertuples(index=False):
        for family, use in PAPER_USES.items():
            if use not in record.allowed_uses:
                continue
            rows.append((record.entity_id, family, record.source_id))
    return pd.DataFrame(rows, columns=["entity_id", "family", "source_id"]).drop_duplicates()


def paper_counts() -> pd.DataFrame:
    """One row per entity: distinct papers in each family, and ``papers_counted`` across
    the families :data:`COUNTED_FAMILIES` (a paper counts once however many tables it has)."""
    table = paper_table()
    counts = table.groupby(["entity_id", "family"]).source_id.nunique().unstack(fill_value=0)
    counts = counts.reindex(columns=list(PAPER_FAMILIES), fill_value=0)
    counts = counts.add_prefix("papers_")
    chosen = table[table.family.isin(COUNTED_FAMILIES)]
    counts["papers_counted"] = chosen.groupby("entity_id").source_id.nunique()
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
                         catalog's own reference: HPA samples, bulk summaries, GTEx medians,
                         brain subregions and immune-cell summaries.
    ``off_target_tpm_tissue`` where that TPM value was measured ("" when it is 0).
    ``tumor_fold``       the catalog's highest tumor / max(N, Disposable, 1).
    ``reproductive_fold`` highest qualified reproductive observation / max(N, D, 1).
    ``reproductive_best`` compatibility alias for reproductive_max_sample_tpm;
                          native TPM only. It no longer mixes assay scales.
    """
    d = _frame(table)
    derived = [f"papers_{family}" for family in PAPER_FAMILIES]
    d = d.drop(columns=[*derived, "papers_counted"], errors="ignore")
    for column in ("is_cancer_antigen", "is_reproductive_antigen", "is_cta", "is_target_candidate"):
        d[column] = d[column].astype(bool)
    d["off_target_tpm"] = d[["normal_reference_tpm", "disposable_max_tpm"]].max(
        axis=1, skipna=False
    )
    tissue = d.normal_reference_contexts.map(_reference_context)
    tissue = tissue.where(d.normal_reference_tpm.ge(d.disposable_max_tpm), "Disposable")
    d["off_target_tpm_tissue"] = tissue.where(d.off_target_tpm.gt(0), "")
    d["reproductive_best"] = d.reproductive_max_sample_tpm
    from . import target_curation

    d = d.drop(columns=[c for c in d if c.startswith("curated_")])
    d = d.join(target_curation.table(d))
    counts = paper_counts()
    d = d.join(counts, on="entity_id")
    for column in counts.columns:
        d[column] = d[column].fillna(0).astype(int)
    return d


def apply(
    f: TargetFilter, table: pd.DataFrame | None = None, *, host_evidence=None
) -> pd.DataFrame:
    """Evaluate every layer for every protein; nothing is dropped.

    Returns :func:`measures` plus ``pass_<layer>`` for each layer, ``keep`` (all
    layers pass), ``filter_failed`` (first failing layer, or "") and
    ``filtered_category`` (the category after tightening compartment evidence).
    """
    d = measures(table)
    host = context(f)
    assessment = host_context.assess(d, host, evidence=host_evidence)
    for name in assessment:
        d[name] = assessment[name]
    h0 = d[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1, skipna=False).clip(lower=1.0)
    d["effective_off_target_tpm"] = np.maximum(d.off_target_tpm, d.host_protected_max_tpm)
    host_denom = np.maximum(h0, d.host_protected_max_tpm)
    d["effective_tumor_fold"] = d.tumor_fold if host is None else d.tumor_max_tpm / host_denom
    d["effective_reproductive_fold"] = d.reproductive_fold
    d["pass_off_target"] = (
        d.effective_off_target_tpm.lt(f.max_off_target_tpm)
        if f.max_off_target_tpm is not None
        else True
    )
    d["pass_curation"] = ~d.curated_tier.isin(NORMAL_FUNCTION[f.normal_function])
    # Tightened biological membership and host eligibility remain separate axes.
    cancer = d.is_cancer_antigen & d.tumor_fold.ge(f.min_tumor_fold)
    host_tumor_pass = d.effective_tumor_fold.ge(f.min_tumor_fold)
    if f.min_tumor_patients:

        def patient_fold(value):
            if not isinstance(value, str):
                return float("nan")
            folds = json.loads(value)
            if (
                not isinstance(folds, list)
                or any(
                    isinstance(v, bool)
                    or not isinstance(v, (int, float))
                    or not math.isfinite(v)
                    or v < 10
                    for v in folds
                )
                or folds != sorted(folds, reverse=True)
            ):
                raise ValueError("Patient folds must be descending finite values >=10")
            return (
                folds[f.min_tumor_patients - 1]
                if len(folds) >= f.min_tumor_patients
                else float("nan")
            )

        folds = d.get("independent_patient_tumor_folds", pd.Series(index=d.index, dtype=object))
        kth_patient_fold = folds.map(patient_fold)
        cancer &= kth_patient_fold.ge(f.min_tumor_fold)
        host_tumor_pass &= (kth_patient_fold * (h0 / host_denom)).ge(f.min_tumor_fold)
    reproductive = d.is_reproductive_antigen & d.reproductive_fold.ge(f.min_reproductive_fold)
    if host is not None:
        # A stricter cancer-evidence requirement may leave a protein in RA.
        # That alone is not host exposure and must not enter the background.
        failed_margin = cancer & ~host_tumor_pass
        d.loc[failed_margin, "host_reasons"] += (
            " | Insufficient independent tumor expression over the protected host reference"
        )
        d.loc[failed_margin & d.host_status.eq("eligible"), "host_status"] = "excluded"
        d["pass_host"] &= ~failed_margin
    d.attrs["host_context"] = host.as_dict() if host else None
    d.attrs["selection_policy_id"] = FILTER_POLICY_ID
    d.attrs["filters"] = asdict(f)
    d["filtered_cancer_antigen"] = cancer
    d["filtered_reproductive_antigen"] = reproductive
    d["filtered_category"] = ""
    d.loc[cancer & reproductive, "filtered_category"] = "CTA"
    d.loc[cancer & ~reproductive, "filtered_category"] = "Cancer Antigen only"
    d.loc[reproductive & ~cancer, "filtered_category"] = "Reproductive Antigen only"
    d["pass_sets"] = MEMBERSHIPS[f.membership](cancer, reproductive).astype(bool)
    d["pass_papers"] = d.papers_counted.ge(f.min_papers)
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
    lines = [
        "Start from the catalog's Target candidates: CA and RA independently require "
        "10 x the non-reproductive normal reference in their respective compartments."
    ]
    if f.max_off_target_tpm is not None:
        lines.append(
            f"Off target: the highest tissue value outside tumors and reproductive tissue, "
            f"Disposable and any protected host compartments included, is below {f.max_off_target_tpm:g} TPM."
        )
    if f.normal_function != "any":
        tiers = " or ".join(t.capitalize() for t in NORMAL_FUNCTION[f.normal_function])
        lines.append(f"Curation: no recorded normal home in a {tiers} tissue.")
    protected = protected_compartments(f)
    if context(f) is not None:
        host = context(f)
        prot_text = ", ".join(protected) if protected else "none"
        absent_text = ", ".join(host.absent) or "none"
        loss_text = ", ".join(host.accept_loss) or "none"
        lines.append(
            f"Host: protect [{prot_text}]; absent [{absent_text}]; accept loss [{loss_text}]. "
            "Other compartments are not protected. "
            f"Protected-tissue RNA limits: < {host.max_tpm:g} TPM and < {host.max_ntpm:g} nTPM. "
            "Missing protected-tissue coverage and unresolved protein findings fail; membership stays biological."
        )
    if f.min_tumor_fold > 10.0:
        lines.append(
            f"Sets: a Cancer Antigen's highest tumor is at least {f.min_tumor_fold:g} x "
            f"max(non-reproductive normal, Disposable, 1 TPM)."
        )
    if f.min_reproductive_fold > 10.0:
        lines.append(
            f"Sets: reproductive evidence reaches at least {f.min_reproductive_fold:g} x max(non-reproductive normal, Disposable, 1 TPM)."
        )
    if f.min_tumor_patients:
        lines.append(
            f"Sets: cancer evidence includes at least {f.min_tumor_patients} known patients at {f.min_tumor_fold:g} x; anonymous samples do not count."
        )
    lines.append(
        f"Sets: keep proteins that are then in {f.membership.replace('_', ' ')}."
        if f.membership != "target_candidates"
        else "Sets: keep proteins that are then in either antigen set."
    )
    if f.min_papers:
        lines.append(
            f"Papers: at least {f.min_papers} distinct paper(s) nominating it or reporting it "
            "in tumors."
        )
    if f.exclude_prefixes:
        text = f"Exclusions: drop labels starting with {'/'.join(f.exclude_prefixes)}"
        if f.except_labels:
            text += f", except {', '.join(f.except_labels)}"
        lines.append(text + ".")
    return "\n".join(lines)


def profile_summary(table: pd.DataFrame | None = None) -> pd.DataFrame:
    """Counts of CA, RA, CTA and their union across canonical host/toxicity profiles."""
    measured = measures(table)
    scenarios = [
        ("all", TargetFilter(host_profile="all")),
        ("female-patient", TargetFilter(host_profile="female-patient")),
        (
            "female-patient (protect ovary)",
            TargetFilter(host_profile="female-patient", protect_tissues=("ovary",)),
        ),
        (
            "female-patient (lactating)",
            TargetFilter(host_profile="female-patient", lactating=True),
        ),
        (
            "female-patient (pregnant)",
            TargetFilter(host_profile="female-patient", pregnant=True),
        ),
        ("male-patient", TargetFilter(host_profile="male-patient")),
        (
            "male-patient (protect testis)",
            TargetFilter(
                host_profile="male-patient",
                protect_tissues=("testis",),
            ),
        ),
        (
            "male-patient (protect adult tissues)",
            TargetFilter(
                host_profile="male-patient",
                protect_tissues=(
                    "testis",
                    "epididymis",
                    "prostate",
                    "seminal vesicle",
                    "resting_breast",
                ),
            ),
        ),
    ]
    rows = []
    for label, filt in scenarios:
        evaluated = apply(filt, measured)
        kept = evaluated[evaluated.keep]
        ca = int(kept.filtered_cancer_antigen.sum())
        ra = int(kept.filtered_reproductive_antigen.sum())
        cta = int((kept.filtered_cancer_antigen & kept.filtered_reproductive_antigen).sum())
        rows.append(
            {
                "profile": label,
                "protected_compartments": ", ".join(protected_compartments(filt)) or "none",
                "allowed_compartments": ", ".join(allowed_compartments(filt)),
                "host_states": json.dumps(context(filt).states() if context(filt) else {}),
                "unresolved": int(
                    (evaluated.is_target_candidate & evaluated.host_status.eq("unresolved")).sum()
                ),
                "CTA_unresolved": int(
                    (evaluated.is_cta & evaluated.host_status.eq("unresolved")).sum()
                ),
                "CTA_excluded": int(
                    (evaluated.is_cta & evaluated.host_status.eq("excluded")).sum()
                ),
                "CTA": cta,
                "CA_only": ca - cta,
                "RA_only": ra - cta,
                "CA": ca,
                "RA": ra,
                "Union": len(kept),
            }
        )
    return pd.DataFrame(rows)
