# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Compatibility helpers; hazard assignments live in :mod:`tissue_hazards`.

``differentiated_airway`` identifies an anatomical family, not a common tier.
Club cells are Vital; mature ciliated and goblet cells are Protected.
"""

from __future__ import annotations

import re

import pandas as pd

from .tissue_hazards import assess

AIRWAY_TISSUES = frozenset({"lung", "bronchus"})
DIFFERENTIATED_AIRWAY = re.compile(
    r"ciliated|club cell|goblet|ionocyte|secretory|serous|mucous", re.I
)


def differentiated_airway(tissue: pd.Series, cell_type: pd.Series) -> pd.Series:
    """True for differentiated (non-basal) airway epithelial cell types."""
    return tissue.fillna("").str.lower().isin(AIRWAY_TISSUES) & cell_type.fillna("").str.contains(
        DIFFERENTIATED_AIRWAY
    )


def cell_tiers(tissue: pd.Series, cell_type: pd.Series) -> pd.Series:
    """Resolve once per unique anatomy/cell label, retaining the caller's index."""
    pairs = list(zip(tissue.fillna(""), cell_type.fillna("")))
    mapping = {p: assess(*p, resolution="cell")["tier"] for p in set(pairs)}
    return pd.Series([mapping[p] for p in pairs], index=tissue.index)
