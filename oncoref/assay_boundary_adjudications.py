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

"""Assay-boundary adjudications (#768) & 1-mismatch Hamming near-self / tryptic-MS audit (#769).

Implements four orthogonal assay-boundary rules (#768) and packaged 1-mismatch Hamming
cross-reactivity plus isobaric (`I == L`) tryptic-peptide MS paralog auditing (#769):

1. **Non-admixed tumor witness gating (`#713` / `#768`)**:
   Requires at least one qualifying ``>= 10x`` tumor witness outside ``library_qc`` Note ``R``
   (``testis_tissue_content >= 0.001``) ``TCGA-TGCT`` parenchymal-admixture libraries. Moves
   ``11`` post-meiotic spermatid markers lacking non-admixed somatic witnesses (``CETN1``,
   ``PGK2``, ``PRM2``, + ``8`` ``AAFL`` sole witnesses) from ``CTA -> Reproductive Antigen
   only``, while ``PRM1`` stays ``CTA`` via its female ``BRCA`` witness (``58.97 >= 57.11 TPM``
   under ``#764`` bleed QC) and ``SPATA19`` stays ``CTA`` via its ``#737`` ``cSCC`` somatic
   witness (``GSM3567709 = 13.28 >= 10.0 TPM``).

2. **RNA-zero single-antibody HPA IHC reconciliation (`#768`)**:
   Reconciles curated cross-reactive HPA IHC exclusions (``cta-ihc-unreliable.csv``) only when
   ``required_coverage_complete == True``, somatic normal RNA ``< 0.5 TPM`` and ``< 1.0 nTPM``,
   and no other ``Unsafe`` gate fires — advancing ``MAGEB6`` (``Unsafe -> CTA``) while keeping
   somatically expressed ``MAGEB4``, ``XAGE2``, and ``NLRP9`` fail-closed.

3. **Adult cornea multi-UMI cell support floor (`#768`)**:
   Requires ``>= 10`` cells with ``>= 2`` UMIs across ``GSE155683`` adult cornea (``21,343``
   cells) before ``Adult cornea`` alone marks a candidate ``Unsafe`` — advancing ``DKKL1``
   (``8`` multi-UMI cells total, ``Unsafe -> CTA``) while retaining ``INSL4`` (``64``),
   ``GPHB5`` (``24``), ``QRFP`` (``40``), and ``C20orf202`` (``42``) in ``Unsafe``.

4. **Protamine-diluted germline transcription-factor qualification (`#768`)**:
   Qualifies ``reproductive_evidence`` when somatic normal ``<= 0.10 TPM``, germline bulk
   ``>= 5.0 TPM`` (``>= 50x`` over somatic), ``>= 2`` independent germline donors reach
   ``>= 3.0 TPM``, AND germline single-cell ``nTPM >= 1.0`` with ``0.0`` somatic single-cell
   ``nTPM`` — advancing ``FOXR2``, ``SPATA31A3``, and ``PPP4R3C`` from ``Cancer Antigen only
   -> CTA`` while rejecting single-donor ``DEFB125``.

5. **1-mismatch Hamming near-self & tryptic-MS paralog audit (`#769`)**:
   Exposes ``exact_unique_9mers``, ``safe_2mm_9mers``, ``unsafe_near_self_1mm_pct``, and
   ``unsafe_shared_tryptic_peptides`` across ``target_selection.results()``, ``antigen_sets``,
   ``sequence_bundle.target_selection_dispositions()``, and ``evaluate_human_catalog()``, plus
   ``hamming_1mm_policy()`` for ``PeptideIndex`` and ``tryptic_ms_paralog_audit()`` for
   ``CTAGE1`` vs ``MIA2``/``CTAGE4``/``CTAGE9``, ``TGM4`` vs ``TGM2``, and ``SPATA31A3`` vs
   ``SPATA31C1``/``SPATA31C2``.
"""

from __future__ import annotations

import hashlib
import json
import re
from copy import deepcopy
from importlib.resources import files

import numpy as np
import pandas as pd

from .load_dataset import get_data
from .peptide_subtraction import AMINO_ACIDS, MatchPolicy

ASSAY_BOUNDARY_POLICY_ID = "oncoref.assay-boundary-and-near-self.2026-11-23.1"

ADULT_CORNEA_MIN_MULTI_UMI_CELLS = 10
IHC_RECONCILE_MAX_SOMATIC_TPM = 0.5
IHC_RECONCILE_MAX_SOMATIC_NTPM = 1.0
GERMLINE_TF_MAX_SOMATIC_TPM = 0.10
GERMLINE_TF_MIN_BULK_TPM = 5.0
GERMLINE_TF_MIN_FOLD = 50.0
GERMLINE_TF_MIN_DONORS_GE3 = 2
GERMLINE_TF_MIN_SC_NTPM = 1.0
SC_CIS_SPILLOVER_MAX_DIST_BP = 100_000
SC_CIS_SPILLOVER_MIN_NEIGHBOR_NTPM = 100_000.0
SC_CIS_SPILLOVER_MAX_RATIO = 1e-4
SC_CIS_SPILLOVER_MAX_SOMATIC_BULK_TPM = 1.0

NEAR_SELF_1MM_COLUMNS: tuple[str, ...] = (
    "exact_unique_9mers",
    "safe_2mm_9mers",
    "unsafe_near_self_1mm_pct",
    "unsafe_shared_tryptic_peptides",
)

ASSAY_BOUNDARY_DISPOSITIONS_768: tuple[str, ...] = (
    "reclassified_ra_only_admixed_tgct_witness",
    "retained_cta_independent_somatic_witness",
    "ihc_rna_zero_cleared",
    "retained_unsafe_ihc_or_rna_discordant_guard_failed",
    "adult_cornea_sparse_umi_cleared",
    "retained_unsafe_adult_cornea_multi_umi",
    "germline_tf_protamine_diluted_qualified",
    "retained_ca_only_germline_guard_failed",
    "sc_soup_and_cis_spillover_cleared",
    "retained_unsafe_sc_cluster",
)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _verify_csv_receipt(stem: str) -> dict:
    data_root = files("oncoref").joinpath("data")
    csv_bytes = data_root.joinpath(f"{stem}.csv").read_bytes()
    receipt = json.loads(data_root.joinpath(f"{stem}.json").read_text(encoding="utf-8"))
    if _sha256_bytes(csv_bytes) != receipt["sha256"]:
        raise ValueError(f"{stem}.csv does not match its SHA-256 receipt")
    return receipt


def hamming_1mm_policy() -> MatchPolicy:
    """Return a pure 1-substitution Hamming ``MatchPolicy`` (indels forbidden at ``cost = 2.0``)."""
    return MatchPolicy(max_distance=1.0, insertion_cost=2.0, deletion_cost=2.0)


def tryptic_peptides(
    sequence: str,
    *,
    min_len: int = 6,
    max_len: int = 30,
    max_missed_cleavages: int = 1,
) -> set[str]:
    """In-silico trypsin digest (cleave after ``K``/``R`` not before ``P``, ``0..max_missed_cleavages``)."""
    seq = str(sequence or "").rstrip("*")
    if not seq:
        return set()
    cuts = [0, *[m.end() for m in re.finditer(r"[KR](?!P)", seq)], len(seq)]
    frags = [seq[cuts[i] : cuts[i + 1]] for i in range(len(cuts) - 1) if cuts[i] < cuts[i + 1]]
    out: set[str] = set()
    for i in range(len(frags)):
        for span in range(1, max_missed_cleavages + 2):
            if i + span > len(frags):
                break
            pep = "".join(frags[i : i + span])
            if min_len <= len(pep) <= max_len and set(pep) <= AMINO_ACIDS:
                out.add(pep)
    return out


def evaluate_non_admixed_tumor_witness(
    *,
    entity_id: str,
    max_non_tgct_tpm: float,
    single_tumor_threshold_tpm: float,
    admixture_restricted_ids: set[str] | frozenset[str],
) -> dict:
    """Evaluate whether a CTA relies solely on ``R``-flagged ``TCGA-TGCT`` admixture (#713/#768)."""
    if str(entity_id) not in admixture_restricted_ids:
        return {
            "admixture_sole_witness_excluded": False,
            "disposition": "retained_cta_independent_somatic_witness",
        }
    if float(max_non_tgct_tpm or 0.0) >= float(single_tumor_threshold_tpm) - 1e-9:
        return {
            "admixture_sole_witness_excluded": False,
            "disposition": "retained_cta_independent_somatic_witness",
        }
    return {
        "admixture_sole_witness_excluded": True,
        "disposition": "reclassified_ra_only_admixed_tgct_witness",
    }


def evaluate_ihc_rna_reconciliation(
    *,
    unreliable_ihc_curated: bool,
    sole_ihc_unsafe_gate: bool,
    required_coverage_complete: bool,
    somatic_rna_max_tpm: float,
    somatic_rna_max_ntpm: float,
) -> dict:
    """Evaluate RNA-zero single-antibody HPA IHC reconciliation with complete-RNA guard (#768)."""
    if (
        bool(unreliable_ihc_curated)
        and bool(sole_ihc_unsafe_gate)
        and bool(required_coverage_complete)
        and float(somatic_rna_max_tpm) < IHC_RECONCILE_MAX_SOMATIC_TPM
        and float(somatic_rna_max_ntpm) < IHC_RECONCILE_MAX_SOMATIC_NTPM
    ):
        return {
            "cleared": True,
            "disposition": "ihc_rna_zero_cleared",
        }
    return {
        "cleared": False,
        "disposition": "retained_unsafe_ihc_or_rna_discordant_guard_failed",
    }


def evaluate_adult_cornea_multi_umi(
    *,
    sole_cornea_unsafe_gate: bool,
    rare_cell_established: bool,
    n_cells_ge2_umi: int,
    min_multi_umi_cells: int = ADULT_CORNEA_MIN_MULTI_UMI_CELLS,
) -> dict:
    """Require ``>= min_multi_umi_cells`` multi-UMI cells in ``GSE155683`` adult cornea (#768)."""
    if (
        bool(sole_cornea_unsafe_gate)
        and not bool(rare_cell_established)
        and int(n_cells_ge2_umi) < int(min_multi_umi_cells)
    ):
        return {
            "cleared": True,
            "disposition": "adult_cornea_sparse_umi_cleared",
        }
    return {
        "cleared": False,
        "disposition": "retained_unsafe_adult_cornea_multi_umi",
    }


def evaluate_germline_tf_qualification(
    *,
    baseline_classification: str,
    normal_reference_tpm: float,
    germline_max_tpm: float,
    germline_donors_ge3_tpm: int,
    germline_sc_ntpm: float,
    somatic_sc_ntpm: float,
) -> dict:
    """Evaluate protamine-diluted germline transcription-factor qualification (#768)."""
    n_val = float(normal_reference_tpm)
    g_val = float(germline_max_tpm)
    fold_ok = g_val >= GERMLINE_TF_MIN_FOLD * max(n_val, 1e-6)
    if (
        str(baseline_classification) == "Cancer Antigen only"
        and n_val <= GERMLINE_TF_MAX_SOMATIC_TPM
        and g_val >= GERMLINE_TF_MIN_BULK_TPM
        and fold_ok
        and int(germline_donors_ge3_tpm) >= GERMLINE_TF_MIN_DONORS_GE3
        and float(germline_sc_ntpm) >= GERMLINE_TF_MIN_SC_NTPM
        and float(somatic_sc_ntpm) <= 0.0
    ):
        return {
            "qualifies": True,
            "disposition": "germline_tf_protamine_diluted_qualified",
        }
    return {
        "qualifies": False,
        "disposition": "retained_ca_only_germline_guard_failed",
    }


def evaluate_sc_cluster_soup_and_spillover(
    *,
    sole_sc_cluster_unsafe_gate: bool,
    rare_cell_established: bool,
    unexplained_somatic_clusters_above_ceiling: int,
    max_cis_spillover_ratio: float,
    flagged_organ_hpa_bulk_tpm: float,
    flagged_organ_ihc_positive: bool,
    raw_somatic_clusters_above_ceiling: int = 1,
    soup_cleared_clusters: int = 0,
) -> dict:
    """Evaluate droplet ambient-soup & <=100 kb cis-chromatin spillover reconciliation (#779)."""
    if (
        bool(sole_sc_cluster_unsafe_gate)
        and not bool(rare_cell_established)
        and int(raw_somatic_clusters_above_ceiling) >= 1
        and int(unexplained_somatic_clusters_above_ceiling) == 0
        and (int(soup_cleared_clusters) >= 1 or float(max_cis_spillover_ratio) > 0.0)
        and float(max_cis_spillover_ratio) <= SC_CIS_SPILLOVER_MAX_RATIO
        and float(flagged_organ_hpa_bulk_tpm) < SC_CIS_SPILLOVER_MAX_SOMATIC_BULK_TPM
        and not bool(flagged_organ_ihc_positive)
    ):
        return {
            "cleared": True,
            "disposition": "sc_soup_and_cis_spillover_cleared",
        }
    return {
        "cleared": False,
        "disposition": "retained_unsafe_sc_cluster",
    }


def assay_boundary_adjudications() -> pd.DataFrame:
    """Return SHA-256-verified ``#768`` assay-boundary adjudications."""
    _verify_csv_receipt("assay-boundary-adjudications")
    return get_data("assay-boundary-adjudications").copy()


def candidate_near_self_1mm() -> pd.DataFrame:
    """Return SHA-256-verified ``#769`` 1-mismatch Hamming near-self & tryptic-MS metrics."""
    _verify_csv_receipt("candidate-near-self-1mm")
    return get_data("candidate-near-self-1mm").copy()


def tryptic_ms_paralog_audit() -> pd.DataFrame:
    """Return SHA-256-verified ``I == L`` tryptic-peptide MS paralog sharing audit (#769)."""
    _verify_csv_receipt("tryptic-ms-paralog-audit")
    return get_data("tryptic-ms-paralog-audit").copy()


def assay_boundary_comparison_observations() -> dict:
    """Return verified native-TPM comparison observations for cleared ``#768`` candidates."""
    payload = json.loads(
        files("oncoref")
        .joinpath("data/assay-boundary-comparison-observations.json")
        .read_text(encoding="utf-8")
    )
    if payload.get("schema") != "oncoref.assay-boundary-comparison-observations.v1":
        raise ValueError("Unexpected assay-boundary comparison observations schema")
    return deepcopy(payload)


def annotate_catalog(
    frame: pd.DataFrame, *, enabled: bool = True, integrated: bool = True
) -> pd.DataFrame:
    """Attach ``tsarina`` ``#769`` 1-mismatch Hamming & tryptic-MS columns to a catalog DataFrame."""
    out = frame.copy()
    if not enabled:
        out["exact_unique_9mers"] = 0
        out["safe_2mm_9mers"] = 0
        out["unsafe_near_self_1mm_pct"] = np.nan
        out["unsafe_shared_tryptic_peptides"] = 0
        return out

    ns_df = candidate_near_self_1mm()
    use_int = bool(integrated) and "exact_unique_9mers_integrated" in ns_df.columns
    by_eid: dict[str, tuple[int, int, float, int]] = {}
    for r in ns_df.itertuples(index=False):
        if use_int:
            ex_v = int(r.exact_unique_9mers_integrated)
            s2_v = int(r.safe_2mm_9mers_integrated)
            pct_v = (
                float(r.unsafe_near_self_1mm_pct_integrated)
                if pd.notna(r.unsafe_near_self_1mm_pct_integrated)
                else np.nan
            )
            tr_v = int(r.unsafe_shared_tryptic_peptides_integrated)
        else:
            ex_v = int(r.exact_unique_9mers)
            s2_v = int(r.safe_2mm_9mers)
            pct_v = (
                float(r.unsafe_near_self_1mm_pct)
                if pd.notna(r.unsafe_near_self_1mm_pct)
                else np.nan
            )
            tr_v = int(r.unsafe_shared_tryptic_peptides)
        by_eid[str(r.entity_id)] = (ex_v, s2_v, pct_v, tr_v)

    is_unsafe_col = (
        out["category"].astype(str).eq("Unsafe")
        if "category" in out.columns
        else pd.Series(False, index=out.index)
    )
    ex_col, s2_col, pct_col, tr_col = [], [], [], []
    for eid, is_uns in zip(out["entity_id"].astype(str), is_unsafe_col):
        if is_uns:
            ex, s2, pct, tr = (0, 0, np.nan, 0)
        else:
            ex, s2, pct, tr = by_eid.get(eid, (0, 0, np.nan, 0))
        ex_col.append(ex)
        s2_col.append(s2)
        pct_col.append(pct)
        tr_col.append(tr)

    out["exact_unique_9mers"] = ex_col
    out["safe_2mm_9mers"] = s2_col
    out["unsafe_near_self_1mm_pct"] = pct_col
    out["unsafe_shared_tryptic_peptides"] = tr_col
    return out


__all__ = [
    "ADULT_CORNEA_MIN_MULTI_UMI_CELLS",
    "ASSAY_BOUNDARY_DISPOSITIONS_768",
    "ASSAY_BOUNDARY_POLICY_ID",
    "GERMLINE_TF_MAX_SOMATIC_TPM",
    "GERMLINE_TF_MIN_BULK_TPM",
    "GERMLINE_TF_MIN_DONORS_GE3",
    "GERMLINE_TF_MIN_FOLD",
    "GERMLINE_TF_MIN_SC_NTPM",
    "IHC_RECONCILE_MAX_SOMATIC_NTPM",
    "IHC_RECONCILE_MAX_SOMATIC_TPM",
    "NEAR_SELF_1MM_COLUMNS",
    "SC_CIS_SPILLOVER_MAX_DIST_BP",
    "SC_CIS_SPILLOVER_MAX_RATIO",
    "SC_CIS_SPILLOVER_MAX_SOMATIC_BULK_TPM",
    "SC_CIS_SPILLOVER_MIN_NEIGHBOR_NTPM",
    "annotate_catalog",
    "assay_boundary_adjudications",
    "assay_boundary_comparison_observations",
    "candidate_near_self_1mm",
    "evaluate_adult_cornea_multi_umi",
    "evaluate_germline_tf_qualification",
    "evaluate_ihc_rna_reconciliation",
    "evaluate_non_admixed_tumor_witness",
    "evaluate_sc_cluster_soup_and_spillover",
    "hamming_1mm_policy",
    "tryptic_ms_paralog_audit",
    "tryptic_peptides",
]
