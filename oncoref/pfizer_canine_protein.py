# Licensed under the Apache License, Version 2.0.
"""Pfizer MSV000099744 canine bulk protein abundance.

Use `source()` to inspect provenance and acquisition status without a download.
"""

from __future__ import annotations

from . import expression_registry as _registry
from . import species_expression as _expression
from .expression_registry import ExpressionSource
from .species_assets import DEFAULT_MAX_BYTES
from .species_expression import ExpressionDataset

__all__ = ["read", "source"]


def source() -> ExpressionSource:
    """Return the bundled source record, including its measurement contract."""
    record = _registry.expression_source("pfizer-dog-protein", species=9615)
    if record is None:
        raise KeyError("Missing bundled canine source: pfizer-dog-protein")
    return record


def read(
    *,
    offline: bool = False,
    cache_root=None,
    local_overrides=None,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> ExpressionDataset:
    """Native protein-group iBAQ across tissue pools, not independent donors.

    Uses the pinned ORFdata dataset and shared managed cache. Requires
    oncoref[quantitative] on Python 3.11+. An override names
    ``9615:pfizer-atlas-dog-protein`` and the released dataset.sqlite file.
    """
    return _expression.read_species_expression(
        "pfizer-dog-protein",
        species=9615,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )
