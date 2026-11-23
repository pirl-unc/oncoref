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

"""Quantitative adult human normal-tissue panels and adjudications for unmeasured home tissues (#735).

Integrates unit-separated adult human bulk RNA-seq (linear ``TPM``) and donor-resolved
single-cell RNA-seq (``UMI``) panels across the five ``UNMEASURED_HOME_TISSUES``
(``olfactory epithelium``, ``hair follicle``, ``taste bud``, ``dental pulp``, ``bone``)
plus ``nasal mucosa``, without pooling ``UMI`` counts or targeted qPCR copy numbers
into bulk ``TPM`` baselines.
"""

from __future__ import annotations

import pandas as pd

from .load_dataset import get_data

HOME_TISSUES_735: tuple[str, ...] = (
    "olfactory epithelium",
    "nasal mucosa",
    "hair follicle",
    "taste bud",
    "dental pulp",
    "bone",
)

PANEL_MODALITIES_735: tuple[str, ...] = (
    "bulk_rnaseq_tpm",
    "scrnaseq_umi",
    "targeted_qpcr_copies",
)

HOME_PANEL_DISPOSITIONS_735: tuple[str, ...] = (
    "home_expressed_above_tier_ceiling",
    "home_expressed_insufficient_tumor_fold",
    "home_cleared_panel_passes_tenfold",
    "home_undetected_in_sampled_panel",
    "home_disposable_follicle_expressed_tenfold",
    "already_excluded_upstream",
)


def unmeasured_home_tissue_panels(
    *,
    home_tissue: str | None = None,
    modality: str | None = None,
) -> pd.DataFrame:
    """Return checksum-pinned adult human normal panels covering unmeasured home tissues (#735)."""
    df = get_data("unmeasured-home-tissue-panels")
    if home_tissue is not None:
        if home_tissue not in HOME_TISSUES_735:
            raise ValueError(
                f"Unknown home_tissue {home_tissue!r}; expected one of {HOME_TISSUES_735}"
            )
        df = df[df["home_tissue"].eq(home_tissue)].reset_index(drop=True)
    if modality is not None:
        if modality not in PANEL_MODALITIES_735:
            raise ValueError(
                f"Unknown modality {modality!r}; expected one of {PANEL_MODALITIES_735}"
            )
        df = df[df["modality"].eq(modality)].reset_index(drop=True)
    return df


def unmeasured_home_protein_adjudications(
    *,
    home_tissue: str | None = None,
    disposition: str | None = None,
    only_ca_cta_gated: bool = False,
) -> pd.DataFrame:
    """Return quantitative home-panel adjudications across all 568 ``home_tissue`` proteins (#735)."""
    df = get_data("unmeasured-home-protein-adjudications")
    if home_tissue is not None:
        if home_tissue not in HOME_TISSUES_735:
            raise ValueError(
                f"Unknown home_tissue {home_tissue!r}; expected one of {HOME_TISSUES_735}"
            )
        df = df[df["home_tissue"].eq(home_tissue)].reset_index(drop=True)
    if disposition is not None:
        if disposition not in HOME_PANEL_DISPOSITIONS_735:
            raise ValueError(
                f"Unknown disposition {disposition!r}; expected one of {HOME_PANEL_DISPOSITIONS_735}"
            )
        df = df[df["home_panel_disposition"].eq(disposition)].reset_index(drop=True)
    if only_ca_cta_gated:
        df = df[df["is_ca_cta_home_gated"].astype(bool)].reset_index(drop=True)
    return df


def unmeasured_home_gene_expression(*, gene: str | None = None) -> pd.DataFrame:
    """Return genome-wide Ensembl 115 quantitative expression across the #735 home-tissue panels."""
    df = get_data("unmeasured-home-gene-expression")
    if gene is not None:
        q = str(gene).strip()
        df = df[df["gene_id"].eq(q) | df["symbol"].eq(q)].reset_index(drop=True)
    return df


def unmeasured_home_tissue_summary() -> pd.DataFrame:
    """Summarize #735 home-tissue panel coverage and gated-protein adjudications by tissue."""
    panels = unmeasured_home_tissue_panels()
    adj = unmeasured_home_protein_adjudications()
    rows = []
    for tissue in HOME_TISSUES_735:
        tp = panels[panels["home_tissue"].eq(tissue)]
        ta = adj[adj["home_tissue"].eq(tissue)]
        tg = ta[ta["is_ca_cta_home_gated"].astype(bool)]
        rows.append(
            {
                "home_tissue": tissue,
                "tissue_tier": tp["tissue_tier"].iloc[0] if not tp.empty else "",
                "n_panels": len(tp),
                "has_bulk_tpm": bool(tp["modality"].eq("bulk_rnaseq_tpm").any()),
                "has_scrnaseq_umi": bool(tp["modality"].eq("scrnaseq_umi").any()),
                "total_home_proteins": len(ta),
                "ca_cta_gated_proteins": len(tg),
                "gated_above_tier_ceiling": int(
                    tg["home_panel_disposition"].eq("home_expressed_above_tier_ceiling").sum()
                ),
                "gated_insufficient_tumor_fold": int(
                    tg["home_panel_disposition"].eq("home_expressed_insufficient_tumor_fold").sum()
                ),
                "gated_disposable_follicle_expressed_tenfold": int(
                    tg["home_panel_disposition"]
                    .eq("home_disposable_follicle_expressed_tenfold")
                    .sum()
                ),
                "gated_cleared_panel_passes_tenfold": int(
                    tg["home_panel_disposition"].eq("home_cleared_panel_passes_tenfold").sum()
                ),
                "gated_undetected_in_sampled_panel": int(
                    tg["home_panel_disposition"].eq("home_undetected_in_sampled_panel").sum()
                ),
            }
        )
    return pd.DataFrame(rows)
