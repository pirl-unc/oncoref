"""Versioned direct GTEx tissue medians in native TPM.

The default is the newest registered release. The source's original tissue
columns remain intact, including microdissected compartments in newer releases.
Medians are population summaries, not bounds on individual-donor expression.
"""

from __future__ import annotations

import pandas as pd

from . import reference_data

SOURCE_ID = "gtex_tissue_rna"


def tissue_rna(version: str | None = None) -> pd.DataFrame:
    """Return the native gene × tissue matrix, with concrete source provenance.

    ``version=None`` and ``"latest"`` resolve to the newest registered release.
    Explicit versions such as ``"v8"`` remain available. ``Name`` retains the
    versioned Ensembl identifier and ``Description`` the original source symbol.
    """
    concrete = reference_data.resolve_version(SOURCE_ID, version)
    path = reference_data.ensure(SOURCE_ID, concrete)
    result = pd.read_csv(path, sep="\t", skiprows=2)
    result.attrs.update(reference_data.provenance(SOURCE_ID, concrete))
    result.attrs["unit"] = "TPM"
    result.attrs["statistic"] = "donor median"
    return result
