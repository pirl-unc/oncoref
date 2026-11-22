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

"""Centralized normal-library (#712) and rare-cell attribution (#711) adjudications."""

from __future__ import annotations

import json
from copy import deepcopy
from importlib.resources import files

import pandas as pd

from .load_dataset import get_data

ADJUDICATION_DISPOSITIONS_712: tuple[str, ...] = (
    "restored_to_cta_under_active_p95",
    "superseded_by_independent_gtex_p95",
    "retained_multi_tissue_cross_atlas_concordant",
    "retained_paralog_multimapping_unverified",
    "retained_sole_adult_tissue_panel_biopsy",
)

ATTRIBUTION_STATUSES_711: tuple[str, ...] = (
    "confounder_rejected",
    "single_study_unconfirmed",
    "unresolved",
)

RESOLUTION_REQUIREMENTS_711: tuple[str, ...] = (
    "resolved_unsafe_in_forbidden_background",
    "independent_second_donor_study_or_targeted_ish_ihc",
    "unique_locus_multimapper_excluded_umi_realignment",
    "spliced_exon_junction_or_polya_bam_attribution",
)


def normal_library_peak_adjudications(*, disposition: str | None = None) -> pd.DataFrame:
    """Return per-protein and per-library adjudications for the 8 #712 HPA-limited candidates."""
    df = get_data("normal-library-peak-adjudications")
    if disposition is not None:
        if disposition not in ADJUDICATION_DISPOSITIONS_712:
            raise ValueError(
                f"Unknown disposition {disposition!r}; expected one of {ADJUDICATION_DISPOSITIONS_712}"
            )
        df = df[df["adjudication_disposition"].eq(disposition)].reset_index(drop=True)
    return df


def hpa_normal_statistic_sensitivity() -> pd.DataFrame:
    """Return genome-wide catalog sensitivity to HPA v23 individual-biopsy statistics (#712)."""
    return get_data("hpa-normal-statistic-sensitivity")


def rare_cell_attribution_reviews() -> dict:
    """Return the validated Ensembl 115 rare-cell host/paralog attribution review ledger (#711)."""
    payload = json.loads(
        files("oncoref").joinpath("data/rare-cell-attribution-reviews.json").read_text()
    )
    if payload.get("schema") != "oncoref.rare-cell-attribution-reviews":
        raise ValueError("Unexpected rare-cell attribution review schema")
    return deepcopy(payload)


def rare_cell_candidate_adjudications(
    *,
    attribution_status: str | None = None,
    resolution_requirement: str | None = None,
) -> pd.DataFrame:
    """Return multi-modal rare-cell, bulk, and background-propagation adjudications (#711)."""
    df = get_data("rare-cell-candidate-adjudications")
    if attribution_status is not None:
        if attribution_status not in ATTRIBUTION_STATUSES_711:
            raise ValueError(
                f"Unknown attribution_status {attribution_status!r}; expected one of {ATTRIBUTION_STATUSES_711}"
            )
        df = df[df["attribution_status"].eq(attribution_status)].reset_index(drop=True)
    if resolution_requirement is not None:
        if resolution_requirement not in RESOLUTION_REQUIREMENTS_711:
            raise ValueError(
                f"Unknown resolution_requirement {resolution_requirement!r}; "
                f"expected one of {RESOLUTION_REQUIREMENTS_711}"
            )
        df = df[df["resolution_requirement"].eq(resolution_requirement)].reset_index(drop=True)
    return df


def unresolved_tenfold_rare_cell_entity_ids() -> frozenset[str]:
    """Return ``entity_id``s of #711 ``Other:Review`` candidates still held on normal-cell attribution."""
    df = rare_cell_candidate_adjudications()
    held = df[df["active_classification"].eq("Other:Review")]
    return frozenset(held["entity_id"].astype(str))
