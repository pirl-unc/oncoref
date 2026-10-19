# Licensed under the Apache License, Version 2.0.
"""GSE304066 canine osteosarcoma metadata; measurement provenance is unresolved.

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
    record = _registry.expression_source("gse304066-osteosarcoma", species=9615)
    if record is None:
        raise KeyError("Missing bundled canine source: gse304066-osteosarcoma")
    return record


def read(
    *,
    offline: bool = False,
    cache_root=None,
    local_overrides=None,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> ExpressionDataset:
    """Raises with the recorded provenance conflict before downloading.

    Assets use the shared verified cache or taxon-qualified local overrides.
    """
    return _expression.read_species_expression(
        "gse304066-osteosarcoma",
        species=9615,
        offline=offline,
        cache_root=cache_root,
        local_overrides=local_overrides,
        max_bytes=max_bytes,
    )
