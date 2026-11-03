# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""One cell-type tier table shared by the RNA cascade and the rare-cell screen.

Tier follows how much harm losing a cell type does. Neurons, glia, cardiomyocytes,
alveolar epithelium and the endothelium of vital organs are Vital: post-mitotic or
needed for gas exchange and perfusion. Differentiated airway epithelium (ciliated,
club, goblet, secretory cells and ionocytes) is a renewing barrier replaced from
basal stem cells, so it is Manageable, like fibroblasts and bowel epithelium.
Airway basal cells are the stem pool that makes that renewal possible and stay Vital.
This is a policy judgement, not a measurement of tolerability.
"""

from __future__ import annotations

import re

import pandas as pd

AIRWAY_TISSUES = frozenset({"lung", "bronchus"})
DIFFERENTIATED_AIRWAY = re.compile(
    r"ciliated|club cell|goblet|ionocyte|secretory|serous|mucous", re.I
)


def differentiated_airway(tissue: pd.Series, cell_type: pd.Series) -> pd.Series:
    """True for differentiated (non-basal) airway epithelial cell types."""
    return tissue.fillna("").str.lower().isin(AIRWAY_TISSUES) & cell_type.fillna("").str.contains(
        DIFFERENTIATED_AIRWAY
    )
