"""Explicit tissue-loss choices, separate from biological antigen membership.

These are research screening assumptions, not clinical safety determinations.
An absent compartment contributes no local target tissue; accepted loss still
has possible acute inflammatory consequences. Neither changes global Unsafe.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass

import pandas as pd

from . import target_selection as catalog

POLICY_ID = "oncoref.host-selection.2026-11-25.1"
COMPARTMENTS = catalog.REPRODUCTIVE_COMPARTMENTS
SEX_RESTRICTION_CODES: tuple[str, ...] = (
    "CTA-M",
    "CTA-F",
    "CTA-MF",
    "CTA-P",
    "RA-M-only",
    "RA-F-only",
    "RA-MF-only",
    "RA-P-only",
)


def compartments(species="human") -> tuple[str, ...]:
    """Return canonical reproductive and mammary HostContext compartments for ``species``."""
    from .species import resolve_species
    from .species_target_selection import anatomy_for_species

    return anatomy_for_species(resolve_species(species)).reproductive_compartments


_compartments_for_species = compartments


def normalize(values, *, species="human") -> tuple[str, ...]:
    from .species import resolve_species
    from .species_target_selection import anatomy_for_species

    return anatomy_for_species(resolve_species(species)).normalize_compartments(values)


@dataclass(frozen=True)
class HostContext:
    """Protect only explicitly requested compartments; default protection is empty.

    The RNA limits are independent research cutoffs in their native units. They
    do not convert nTPM to TPM or establish tolerability. A generic compartment
    label always retains whole-organ scope, including endocrine/support cells.

    Pregnancy & developmental completeness (#790)
    ---------------------------------------------
    Protecting ``"embryo"`` (including via ``for_host(..., pregnant=True)``) is
    intentionally fail-closed (0 eligible candidates across human, dog, and mouse)
    by default because ``embryo`` carries ``scope_complete=False`` on every entity:
    integrated blastocyst/gastrulation/early-embryo libraries do not cover 1st–3rd
    trimester fetal organogenesis. Callers have two explicit screening options:

    1. Screen the adult maternal-fetal & lactational interface directly via
       ``HostContext(protect=("placenta", "lactating_breast"))`` (or
       ``"lactating_mammary"`` in dog/mouse), where ``scope_complete=True``.
    2. Opt in via ``allow_incomplete_scope=True`` to evaluate measured early-embryo,
       placenta, and lactating-mammary expression (``< max_tpm`` / ``< max_ntpm``)
       without failing on ``sampled scope incomplete`` (unmeasured ``coverage !=
       "measured"`` and unresolved/positive protein findings still fail closed).
    """

    absent: tuple[str, ...] = ()
    accept_loss: tuple[str, ...] = ()
    max_tpm: float = 8.0
    max_ntpm: float = 8.0
    assumptions: tuple[str, ...] = ()
    protect: tuple[str, ...] = ()
    species: str = "human"
    negative_ihc_multiplier: float | None = None
    allow_incomplete_scope: bool = False

    def __post_init__(self):
        from .species import resolve_species

        resolved_species = resolve_species(self.species)
        object.__setattr__(self, "species", resolved_species.name)
        if not isinstance(self.allow_incomplete_scope, bool):
            raise ValueError("allow_incomplete_scope must be a boolean")
        if self.negative_ihc_multiplier is None:
            from .species_policies import species_policy

            ihc_mult = species_policy(resolved_species).thresholds.host_negative_ihc_multiplier
        else:
            ihc_mult = self.negative_ihc_multiplier
        if (
            isinstance(ihc_mult, bool)
            or not isinstance(ihc_mult, (int, float))
            or not math.isfinite(ihc_mult)
            or ihc_mult < 1.0
        ):
            raise ValueError("negative_ihc_multiplier must be finite and >= 1.0")
        object.__setattr__(self, "negative_ihc_multiplier", float(ihc_mult))
        for field in ("absent", "accept_loss", "protect"):
            object.__setattr__(
                self, field, normalize(getattr(self, field), species=resolved_species)
            )
        if set(self.absent) & set(self.accept_loss):
            raise ValueError("A compartment cannot be both absent and accepted for loss")
        if set(self.protect) & (set(self.absent) | set(self.accept_loss)):
            raise ValueError("A protected compartment cannot be absent or accepted for loss")
        for field in ("max_tpm", "max_ntpm"):
            value = getattr(self, field)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ValueError(f"{field} must be finite and positive")
        if not isinstance(self.assumptions, tuple) or any(
            not isinstance(v, str) for v in self.assumptions
        ):
            raise ValueError("assumptions must be a tuple of strings")

    @classmethod
    def for_host(
        cls,
        profile="all",
        *,
        pregnant: bool = False,
        lactating: bool = False,
        allow_tissues=(),
        protect_tissues=(),
        absent_tissues=(),
        species="human",
        allow_incomplete_scope: bool = False,
        max_tpm: float = 8.0,
        max_ntpm: float = 8.0,
        negative_ihc_multiplier: float | None = None,
    ) -> HostContext | None:
        """Build a species-aware ``HostContext`` from host profile and state choices (#790)."""
        return for_host(
            profile,
            pregnant=pregnant,
            lactating=lactating,
            allow_tissues=allow_tissues,
            protect_tissues=protect_tissues,
            absent_tissues=absent_tissues,
            species=species,
            allow_incomplete_scope=allow_incomplete_scope,
            max_tpm=max_tpm,
            max_ntpm=max_ntpm,
            negative_ihc_multiplier=negative_ihc_multiplier,
        )

    @property
    def protected(self) -> tuple[str, ...]:
        return self.protect

    def as_dict(self) -> dict:
        return {"policy_id": POLICY_ID, **asdict(self), "protected": self.protected}

    def states(self) -> dict[str, str]:
        return {
            c: "Absent"
            if c in self.absent
            else "Accept loss"
            if c in self.accept_loss
            else "Protect"
            if c in self.protect
            else "Not protected"
            for c in compartments(self.species)
        }


def for_host(
    profile="all",
    *,
    pregnant: bool = False,
    lactating: bool = False,
    allow_tissues=(),
    protect_tissues=(),
    absent_tissues=(),
    species="human",
    allow_incomplete_scope: bool = False,
    max_tpm: float = 8.0,
    max_ntpm: float = 8.0,
    negative_ihc_multiplier: float | None = None,
) -> HostContext | None:
    """Build a species-aware ``HostContext`` from profile/physiological arguments (#790).

    Retains every anatomy and physiological-state assumption on ``HostContext.assumptions``.
    When ``pregnant=True``, protects ``("embryo", "placenta", <lactating_mammary>)``;
    because ``embryo`` carries ``scope_complete=False`` (no 1st–3rd trimester fetal
    organogenesis panel), ``assess()`` fails closed (0 eligible) unless
    ``allow_incomplete_scope=True`` is explicitly requested or ``embryo`` is omitted.
    """
    from .species import resolve_species

    if not isinstance(allow_incomplete_scope, bool):
        raise ValueError("allow_incomplete_scope must be a boolean")
    resolved_species = resolve_species(species)
    allow, protect, absent_override = (
        normalize(values, species=resolved_species)
        for values in (allow_tissues, protect_tissues, absent_tissues)
    )
    if (set(allow) & set(protect)) or (set(absent_override) & (set(allow) | set(protect))):
        raise ValueError("Compartment choices conflict")
    if profile == "core":
        raise ValueError(
            "core describes biological origin, not host protection; use explicit HostContext protection choices"
        )
    if profile not in {"all", "male-patient", "female-patient"}:
        raise ValueError("Unknown host profile")
    if not any((profile != "all", pregnant, lactating, allow, protect, absent_override)):
        return None
    from .species_target_selection import anatomy_for_species

    anat = anatomy_for_species(resolved_species)
    female_tract = anat.female_tract_compartments
    male_accessory = anat.male_accessory_compartments
    lactating_comp = anat.lactating_mammary_compartment
    absent, accepted, protected, notes = set(), set(), set(protect), []
    if profile in {"male-patient", "female-patient"}:
        opposite = (
            (
                *female_tract,
                "ovary",
                *((anat.resting_mammary_compartment,) if anat.taxon_id == 10090 else ()),
            )
            if profile == "male-patient"
            else (*male_accessory, "testis")
        )
        absent.update(opposite)
        notes.append(
            f"{profile}: assumed anatomy; verify actual organs rather than inferring them from sex or chromosomes"
        )
        if not pregnant:
            absent.update(("embryo", "placenta"))
            notes.append(
                "Assumes no current pregnancy; does not assess future pregnancy during persistent treatment"
            )
        if not lactating and not pregnant:
            absent.add(lactating_comp)
            notes.append(
                "Assumes inactive milk-secretory state, not an absent breast; state inference has limited sampling"
            )
    else:
        notes.append("Custom context: only explicitly requested compartments protected")
    notes.append("Unspecified reproductive compartments are not protected by this screen")
    if pregnant:
        absent.difference_update(("embryo", "placenta", lactating_comp))
        protected.update(("embryo", "placenta", lactating_comp))
    if lactating:
        absent.discard(lactating_comp)
        protected.add(lactating_comp)
        accepted.discard(lactating_comp)
    absent.difference_update((*allow, *protect))
    accepted.difference_update((*protect, *absent_override))
    absent.update(absent_override)
    accepted.update(allow)
    protected.difference_update((*allow, *absent_override))
    if allow_incomplete_scope:
        notes.append(
            "allow_incomplete_scope=True: waives sampled-scope completeness holds (e.g. early-embryo "
            "atlases lack 1st-3rd trimester fetal organogenesis coverage)"
        )
    return HostContext(
        absent=tuple(absent),
        accept_loss=tuple(accepted),
        protect=tuple(protected),
        assumptions=tuple(notes),
        species=resolved_species.name,
        max_tpm=max_tpm,
        max_ntpm=max_ntpm,
        negative_ihc_multiplier=negative_ihc_multiplier,
        allow_incomplete_scope=allow_incomplete_scope,
    )


#: Non-deprecated compatibility alias for :func:`for_host`.
from_legacy = for_host


def _evidence(species="human") -> dict:
    from .species_target_selection import _packaged_species_host_exposure_ref

    return _packaged_species_host_exposure_ref(species)


def exposure_evidence(species="human") -> dict:
    """Detached per-compartment observations and coverage, never clinical clearance."""
    from .species_target_selection import packaged_species_host_exposure

    return packaged_species_host_exposure(species)


def assess(table, host: HostContext | None, *, evidence=None) -> pd.DataFrame:
    """Evaluate protection independently of target-set membership and paper counts.

    Protected compartments require ``coverage == "measured"``, finite RNA below
    ``max_tpm`` / ``max_ntpm`` (scaled by ``negative_ihc_multiplier`` when routine
    IHC is ``not_detected``), no positive/unresolved protein finding, and
    ``scope_complete is True`` (unless ``host.allow_incomplete_scope=True`` is
    explicitly set, e.g. to screen measured early-embryo stages despite missing
    1st–3rd trimester fetal organogenesis coverage).
    """
    if host is not None and not isinstance(host, HostContext):
        raise ValueError("host must be a HostContext")
    if host is not None and host.protected and evidence is None:
        try:
            evidence = _evidence(host.species)
        except TypeError:
            evidence = _evidence()
    else:
        evidence = evidence or {}
    is_female_host = False
    is_male_host = False
    if host is not None:
        species_comps = set(compartments(host.species))
        fem_tract = tuple(
            c
            for c in ("uterus", "uterine_horn", "vagina", "cervix", "endometrium", "fallopian tube")
            if c in species_comps
        )
        male_tract = tuple(
            c
            for c in ("prostate", "epididymis", "seminal vesicle", "coagulating_gland")
            if c in species_comps
        )
        male_all_comps = tuple(c for c in ("testis", *male_tract) if c in species_comps)
        fem_all_comps = tuple(c for c in ("ovary", *fem_tract) if c in species_comps)
        lact_comps = tuple(
            c for c in ("lactating_breast", "lactating_mammary") if c in species_comps
        )
        declared_female = any(a.startswith("female-patient:") for a in host.assumptions)
        declared_male = any(a.startswith("male-patient:") for a in host.assumptions)
        active_host_comps = set(host.protect) | set(host.accept_loss)
        is_female_host = declared_female or (
            not declared_male
            and bool(male_all_comps)
            and all(c in host.absent for c in male_all_comps)
            and not any(c in active_host_comps for c in male_all_comps)
            and (
                any(c in active_host_comps for c in (*fem_all_comps, *lact_comps))
                or any(c not in host.absent for c in fem_all_comps)
            )
        )
        is_male_host = not is_female_host and (
            declared_male
            or (
                bool(fem_all_comps)
                and all(c in host.absent for c in (*fem_all_comps, *lact_comps))
                and not any(c in active_host_comps for c in (*fem_all_comps, *lact_comps))
                and (
                    any(c in active_host_comps for c in male_all_comps)
                    or any(c not in host.absent for c in male_all_comps)
                )
            )
        )
    rows = []
    for row in table.itertuples():
        rejected, unknown, observed = [], [], []
        for tissue in () if host is None else host.protected:
            item = evidence.get(row.entity_id, {}).get(tissue)
            if not isinstance(item, dict) or item.get("coverage") != "measured":
                unknown.append(tissue)
                continue
            value = item.get("max_tpm")
            if tissue in {"resting_breast", "resting_mammary"} and host is not None:
                if is_female_host and item.get("female_max_tpm") is not None:
                    value = item.get("female_max_tpm")
                elif is_male_host and item.get("male_max_tpm") is not None:
                    value = item.get("male_max_tpm")
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                unknown.append(tissue)
                continue
            observed.append((value, tissue))
            ihc_mult = (
                float(host.negative_ihc_multiplier)
                if item.get("protein_status") == "not_detected"
                else 1.0
            )
            eff_max_tpm = host.max_tpm * ihc_mult
            eff_max_ntpm = host.max_ntpm * ihc_mult
            if value >= eff_max_tpm:
                rejected.append(f"{tissue}: RNA {value:g} TPM >= {eff_max_tpm:g}")
            ntpm = item.get("max_ntpm")
            if tissue in {"resting_breast", "resting_mammary"} and host is not None:
                if is_female_host and item.get("female_max_ntpm") is not None:
                    ntpm = item.get("female_max_ntpm")
                elif is_male_host and item.get("male_max_ntpm") is not None:
                    ntpm = item.get("male_max_ntpm")
            if ntpm is not None:
                if (
                    isinstance(ntpm, bool)
                    or not isinstance(ntpm, (int, float))
                    or not math.isfinite(ntpm)
                    or ntpm < 0
                ):
                    unknown.append(tissue)
                elif ntpm >= eff_max_ntpm:
                    rejected.append(f"{tissue}: RNA {ntpm:g} nTPM >= {eff_max_ntpm:g}")
            if item.get("protein_status") in {"positive", "unresolved"}:
                female_only_mammary_ihc = (
                    tissue in {"resting_breast", "resting_mammary"}
                    and is_male_host
                    and item.get("protein_status") == "unresolved"
                    and value < eff_max_tpm
                    and (ntpm is None or ntpm < eff_max_ntpm)
                )
                if not female_only_mammary_ihc:
                    unknown.append(tissue + ": protein finding requires host-specific review")
            if item.get("scope_complete") is not True and not host.allow_incomplete_scope:
                unknown.append(tissue + ": sampled scope incomplete")
        maximum, limiting = max(observed, default=(0.0, ""))
        status = (
            "not_applied"
            if host is None
            else "excluded"
            if rejected
            else "unresolved"
            if unknown
            else "eligible"
        )
        rows.append(
            {
                "host_status": status,
                "host_protected_max_tpm": maximum if not unknown else float("nan"),
                "host_limiting_compartment": limiting,
                "host_reasons": " | ".join(rejected + ["Unresolved: " + t for t in unknown]),
                "host_unknown_compartments": json.dumps(
                    sorted({reason.split(": ", 1)[0] for reason in unknown})
                ),
                "pass_host": status in {"eligible", "not_applied"},
            }
        )
    return pd.DataFrame(
        rows,
        index=table.index,
        columns=[
            "host_status",
            "host_protected_max_tpm",
            "host_limiting_compartment",
            "host_reasons",
            "host_unknown_compartments",
            "pass_host",
        ],
    )


__all__ = [
    "COMPARTMENTS",
    "POLICY_ID",
    "SEX_RESTRICTION_CODES",
    "HostContext",
    "assess",
    "compartments",
    "exposure_evidence",
    "for_host",
    "from_legacy",
    "normalize",
]
