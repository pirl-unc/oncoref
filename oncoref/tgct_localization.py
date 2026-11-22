# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""TCGA-TGCT somatic testicular admixture vs malignant germ-cell localization audit (#713).

Exposes the packaged per-library somatic niche / spermatid stoichiometry table
(:func:`tgct_reproductive_localization_audit`) across all 154 ``TCGA-TGCT``
libraries and the per-protein cross-cohort witness sensitivity table
(:func:`tgct_aafl_protein_support`) for the 34 catalog proteins whose peak tumor
measurement occurs in ``TCGA-2G-AAFL-01``.
"""

from __future__ import annotations

import pandas as pd

from .load_dataset import get_data

TGCT_LOCALIZATION_AUDIT_ID = "oncoref.tgct-localization-audit.2026-11-22.1"
AAFL_ADMIXED_SOLE_WITNESS_TIER = "aafl_sole_witness_admixed_parenchyma_candidate"
R_TGCT_MULTI_LIBRARY_ONLY_TIER = "r_tgct_multi_library_only"


def tgct_reproductive_localization_audit(*, r_only: bool = False) -> pd.DataFrame:
    """Per-library somatic testicular niche (Leydig + Sertoli) and spermatid stoichiometry
    across all 154 ``TCGA-TGCT`` libraries.

    Parameters
    ----------
    r_only : bool, default False
        When ``True``, restrict to the 7 ``library_qc`` Note ``R`` libraries
        (``testis_tissue_content >= 0.001``).
    """
    df = get_data("tgct-reproductive-localization-audit").copy()
    if r_only:
        df = df[df["library_qc_note_R"].astype(bool)].reset_index(drop=True)
    return df


def tgct_aafl_protein_support(*, support_tier: str | None = None) -> pd.DataFrame:
    """Cross-cohort witness audit for all 34 catalog proteins peaking in ``TCGA-2G-AAFL-01``.

    Distinguishes CTAs independently confirmed outside ``R``-flagged ``TGCT`` libraries
    (``confirmed_outside_r_tgct``) from spermatid candidates whose qualifying tumor
    evidence comes solely from ``TCGA-2G-AAFL-01`` (``aafl_sole_witness_admixed_parenchyma_candidate``)
    or solely from ``R``-flagged ``TGCT`` libraries (``r_tgct_multi_library_only``).
    """
    df = get_data("tgct-aafl-protein-support").copy()
    if support_tier is not None:
        df = df[df["support_tier"].astype(str).eq(str(support_tier))].reset_index(drop=True)
    return df


def tgct_admixture_sole_witness_entity_ids() -> frozenset[str]:
    """Entity IDs (11 CTAs) whose Cancer Antigen qualification depends solely on ``TCGA-2G-AAFL-01``."""
    df = tgct_aafl_protein_support(support_tier=AAFL_ADMIXED_SOLE_WITNESS_TIER)
    return frozenset(df["entity_id"].astype(str))


def tgct_admixture_restricted_entity_ids(*, include_multi_r_tgct: bool = True) -> frozenset[str]:
    """Entity IDs whose Cancer Antigen qualification lacks any witness outside ``R``-flagged TGCT.

    Includes the 11 ``TCGA-2G-AAFL-01`` sole-witness candidates plus, when
    ``include_multi_r_tgct=True`` (default), ``PRM1`` and ``PRM2`` (13 total).
    """
    tiers = {AAFL_ADMIXED_SOLE_WITNESS_TIER}
    if include_multi_r_tgct:
        tiers.add(R_TGCT_MULTI_LIBRARY_ONLY_TIER)
    df = tgct_aafl_protein_support()
    return frozenset(df.loc[df["support_tier"].isin(tiers), "entity_id"].astype(str))
