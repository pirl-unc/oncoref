# Licensed under the Apache License, Version 2.0.
"""Bgee SRP181218 canine gene TPM.

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
    record = _registry.expression_source("bgee-srp181218", species=9615)
    if record is None:
        raise KeyError("Missing bundled canine source: bgee-srp181218")
    return record


def read(
    *,
    offline: bool = False,
    cache_root=None,
    local_overrides=None,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> ExpressionDataset:
    """Native gene TPM with evidence-backed sample tissue corrections.

    Assets use the shared verified cache or taxon-qualified local overrides.
    """
    return _expression.read_species_expression(
        "bgee-srp181218",
        species=9615,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )
