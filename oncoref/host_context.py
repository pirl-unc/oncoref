"""Explicit tissue-loss choices, separate from biological antigen membership.

These are research screening assumptions, not clinical safety determinations.
An absent compartment contributes no local target tissue; accepted loss still
has possible acute inflammatory consequences. Neither changes global Unsafe.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from functools import cache
from importlib.resources import files

import pandas as pd

from . import target_selection as catalog

POLICY_ID = "oncoref.host-selection.2026-11-17.2"
COMPARTMENTS = catalog.REPRODUCTIVE_COMPARTMENTS


def normalize(values) -> tuple[str, ...]:
    if not isinstance(values, tuple):
        raise ValueError("Tissue choices must be tuples of names")
    found = set()
    aliases = {c.replace("_", " "): (c,) for c in COMPARTMENTS}
    aliases.update(breast=catalog.MAMMARY_COMPARTMENTS, mammary=catalog.MAMMARY_COMPARTMENTS)
    for name in values:
        if not isinstance(name, str):
            raise ValueError("Tissue names must be strings")
        key = " ".join(name.lower().replace("_", " ").split())
        if key not in aliases:
            raise ValueError(f"Unknown compartment {name!r}")
        found.update(aliases[key])
    return tuple(c for c in COMPARTMENTS if c in found)


@dataclass(frozen=True)
class HostContext:
    """Protect every compartment except explicitly absent or accepted-loss ones.

    The RNA limits are independent research cutoffs in their native units. They
    do not convert nTPM to TPM or establish tolerability. A generic compartment
    label always retains whole-organ scope, including endocrine/support cells.
    """

    absent: tuple[str, ...] = ()
    accept_loss: tuple[str, ...] = ()
    max_tpm: float = 4.0
    max_ntpm: float = 4.0
    assumptions: tuple[str, ...] = ()

    def __post_init__(self):
        for field in ("absent", "accept_loss"):
            object.__setattr__(self, field, normalize(getattr(self, field)))
        if set(self.absent) & set(self.accept_loss):
            raise ValueError("A compartment cannot be both absent and accepted for loss")
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
        return tuple(c for c in COMPARTMENTS if c not in {*self.absent, *self.accept_loss})

    def as_dict(self) -> dict:
        return {"policy_id": POLICY_ID, **asdict(self), "protected": self.protected}

    def states(self) -> dict[str, str]:
        return {
            c: "Absent"
            if c in self.absent
            else "Accept loss"
            if c in self.accept_loss
            else "Protect"
            for c in COMPARTMENTS
        }


def from_legacy(
    profile="all",
    *,
    pregnant=False,
    lactating=False,
    allow_tissues=(),
    protect_tissues=(),
    absent_tissues=(),
):
    """Translate compatibility arguments; retain every anatomy/state assumption."""
    allow, protect, absent_override = map(
        normalize, (allow_tissues, protect_tissues, absent_tissues)
    )
    if (set(allow) & set(protect)) or (set(absent_override) & (set(allow) | set(protect))):
        raise ValueError("Compartment choices conflict")
    if profile == "core":
        raise ValueError(
            "core describes biological origin, not host protection; use an explicit HostContext and accept_loss choices"
        )
    if profile not in {"all", "male-patient", "female-patient"}:
        raise ValueError("Unknown host profile")
    if not any((profile != "all", pregnant, lactating, allow, protect, absent_override)):
        return None
    absent, accepted, notes = set(), set(), []
    if profile in {"male-patient", "female-patient"}:
        opposite = (
            (*catalog.FEMALE_TRACT_COMPARTMENTS, "ovary")
            if profile == "male-patient"
            else (*catalog.MALE_ACCESSORY_COMPARTMENTS, "testis")
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
            absent.add("lactating_breast")
            notes.append(
                "Assumes inactive milk-secretory state, not an absent breast; state inference has limited sampling"
            )
    else:
        # A custom policy must name all exceptions. Unspecified organs stay protected.
        notes.append("Custom context: all unspecified compartments protected")
    if pregnant:
        absent.difference_update(("embryo", "placenta", "lactating_breast"))
        accepted.difference_update(("embryo", "placenta", "lactating_breast"))
    if lactating:
        absent.discard("lactating_breast")
        accepted.discard("lactating_breast")
    absent.difference_update((*allow, *protect))
    accepted.difference_update((*protect, *absent_override))
    absent.update(absent_override)
    accepted.update(allow)
    return HostContext(absent=tuple(absent), accept_loss=tuple(accepted), assumptions=tuple(notes))


@cache
def _evidence():
    import hashlib

    root = files("oncoref") / "data"
    path = root / "host-exposure.csv.gz"
    meta = json.loads((root / "host-exposure.json").read_text())
    if meta["policy_id"] != POLICY_ID or meta["catalog_policy_id"] != catalog.POLICY_ID:
        raise ValueError("Host exposure belongs to another policy")
    if meta.get("catalog_data_sha256") != catalog.metadata()["data_sha256"]:
        raise ValueError("Host exposure belongs to another catalog snapshot")
    if hashlib.sha256(path.read_bytes()).hexdigest() != meta["sha256"]:
        raise ValueError("Host exposure does not match its receipt")
    return (
        pd.read_csv(path, keep_default_na=False)
        .set_index("entity_id")
        .exposures.map(json.loads)
        .to_dict()
    )


def exposure_evidence() -> dict:
    """Detached per-compartment observations and coverage, never clinical clearance."""
    from copy import deepcopy

    return deepcopy(_evidence())


def assess(table, host: HostContext | None, *, evidence=None) -> pd.DataFrame:
    """Evaluate protection independently of target-set membership and paper counts."""
    if host is not None and not isinstance(host, HostContext):
        raise ValueError("host must be a HostContext")
    evidence = _evidence() if host is not None and evidence is None else (evidence or {})
    rows = []
    for row in table.itertuples():
        rejected, unknown, observed = [], [], []
        for tissue in () if host is None else host.protected:
            item = evidence.get(row.entity_id, {}).get(tissue)
            if not isinstance(item, dict) or item.get("coverage") != "measured":
                unknown.append(tissue)
                continue
            value = item.get("max_tpm")
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
