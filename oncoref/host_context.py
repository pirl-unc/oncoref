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

POLICY_ID = "oncoref.host-selection.2026-11-18.2"
COMPARTMENTS = catalog.REPRODUCTIVE_COMPARTMENTS


def _compartments_for_species(species="human") -> tuple[str, ...]:
    from .species import resolve_species
    from .species_target_selection import anatomy_for_species

    return anatomy_for_species(resolve_species(species)).reproductive_compartments


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
    """

    absent: tuple[str, ...] = ()
    accept_loss: tuple[str, ...] = ()
    max_tpm: float = 4.0
    max_ntpm: float = 4.0
    assumptions: tuple[str, ...] = ()
    protect: tuple[str, ...] = ()
    species: str = "human"

    def __post_init__(self):
        from .species import resolve_species

        resolved_species = resolve_species(self.species)
        object.__setattr__(self, "species", resolved_species.name)
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
            for c in _compartments_for_species(self.species)
        }


def from_legacy(
    profile="all",
    *,
    pregnant=False,
    lactating=False,
    allow_tissues=(),
    protect_tissues=(),
    absent_tissues=(),
    species="human",
):
    """Translate compatibility arguments; retain every anatomy/state assumption."""
    from .species import resolve_species

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
            (*female_tract, "ovary") if profile == "male-patient" else (*male_accessory, "testis")
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
    return HostContext(
        absent=tuple(absent),
        accept_loss=tuple(accepted),
        protect=tuple(protected),
        assumptions=tuple(notes),
        species=resolved_species.name,
    )


def _evidence(species="human") -> dict:
    from .species_target_selection import _packaged_species_host_exposure_ref

    return _packaged_species_host_exposure_ref(species)


def exposure_evidence(species="human") -> dict:
    """Detached per-compartment observations and coverage, never clinical clearance."""
    from .species_target_selection import packaged_species_host_exposure

    return packaged_species_host_exposure(species)


def assess(table, host: HostContext | None, *, evidence=None) -> pd.DataFrame:
    """Evaluate protection independently of target-set membership and paper counts."""
    if host is not None and not isinstance(host, HostContext):
        raise ValueError("host must be a HostContext")
    if host is not None and host.protected and evidence is None:
        try:
            evidence = _evidence(host.species)
        except TypeError:
            evidence = _evidence()
    else:
        evidence = evidence or {}
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
                if (
                    "testis" in host.absent
                    and "ovary" not in host.absent
                    and item.get("female_max_tpm") is not None
                ):
                    value = item.get("female_max_tpm")
                elif (
                    "ovary" in host.absent
                    and "testis" not in host.absent
                    and item.get("male_max_tpm") is not None
                ):
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
            if value >= host.max_tpm:
                rejected.append(f"{tissue}: RNA {value:g} TPM >= {host.max_tpm:g}")
            ntpm = item.get("max_ntpm")
            if ntpm is not None:
                if (
                    isinstance(ntpm, bool)
                    or not isinstance(ntpm, (int, float))
                    or not math.isfinite(ntpm)
                    or ntpm < 0
                ):
                    unknown.append(tissue)
                elif ntpm >= host.max_ntpm:
                    rejected.append(f"{tissue}: RNA {ntpm:g} nTPM >= {host.max_ntpm:g}")
            if item.get("protein_status") in {"positive", "unresolved"}:
                unknown.append(tissue + ": protein finding requires host-specific review")
            if item.get("scope_complete") is not True:
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
