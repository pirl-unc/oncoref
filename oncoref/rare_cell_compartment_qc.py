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

"""Nuclear/intronic snRNA-seq compartment & transcriptional spillover QC (#765).

Separates intron-inclusive single-nucleus RNA-seq (``assay_compartment == "nucleus"``)
and proximal promoter / termination-readthrough spillover from spliced cytoplasmic
``CDS`` mRNA in pMHC safety gates, enforcing downstream cross-repo guardrails:

1. **Strict bulk-dilution & literature/IHC guardrails (``MAGEA12`` regression anchor):**
   A rare-cell signal is never cleared on ``finding_bulk_tpm < 1.0`` alone (``MAGEA12``
   sits at ``0.30 TPM`` whole-brain bulk with ``literature_effect == "confirm_detection"``
   and remains unconditionally ``Unsafe``). Clearance requires ``assay_compartment ==
   "nucleus"`` with zero cytoplasmic ``scRNA-seq`` multi-UMI replication,
   ``finding_bulk_tpm < 0.05``, negative somatic organ IHC, ``literature_effect`` not in
   ``{"exclude", "confirm_detection"}``, and verified Ensembl 115 spillover geometry
   (or a verified closed repeat array under ``#767``).
2. **Intronless / single-CDS-exon support (``DUX4``, ``PLAC1``):**
   Explicitly records ``is_intronless_cds`` and accepts ``cds_exonic_reads`` alongside
   ``spliced_exon_junction`` for coding-sequence attribution so single-exon ORFs are
   neither penalized nor bypassed.
3. **Explicit downstream ``nuclear_only`` disposition:**
   Cleared nuclear-spillover loci receive ``effective_rare_cell_status = "nuclear_only"``
   and ``nuclear_rna_caution = "snrnaseq_intronic_spillover_cleared"``, allowing loci
   that pass all ``10x`` normal/reproductive/tumor gates (such as ``PLAC1`` and ``TENT5D``)
   to classify into ``CTA`` / ``Reproductive Antigen only`` with an auditable machine-
   readable caution column rather than being trapped in ``Unsafe`` or ``Other:Review``.
"""

from __future__ import annotations

import hashlib
import json
import math
from importlib.resources import files

import pandas as pd

from .load_dataset import get_data

NUCLEAR_SPILLOVER_POLICY_ID = "oncoref.nuclear-spillover-qc.2026-11-23.1"

ASSAY_COMPARTMENTS: tuple[str, ...] = ("nucleus", "cell", "mixed")

NUCLEAR_CLEARANCE_BULK_CEILING_TPM = 0.05
BIDIRECTIONAL_PROMOTER_MAX_BP = 1000
BIDIRECTIONAL_PARTNER_MIN_TPM = 10.0
SAME_STRAND_UPSTREAM_MAX_BP = 10000
SAME_STRAND_PARTNER_MIN_TPM = 5.0
NONCODING_5PRIME_SPAN_MIN_BP = 50000
CDS_TO_SPAN_MAX_RATIO = 0.05

NUCLEAR_DISPOSITIONS_765: tuple[str, ...] = (
    "nuclear_spillover_cleared",
    "promoter_spillover_immune_explained",
    "retained_established_cytoplasmic_or_bulk",
    "retained_literature_or_ihc_supported",
    "retained_unresolved_paralog_or_exonic_host",
    "retained_nuclear_no_spillover_geometry",
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


def assess_transcriptional_spillover_geometry(
    *,
    bidirectional_promoter_bp: int | float | None = None,
    bidirectional_partner_tpm: float = 0.0,
    same_strand_upstream_bp: int | float | None = None,
    same_strand_partner_tpm: float = 0.0,
    noncoding_5prime_span_bp: int | float = 0,
    cds_span_ratio: float = 1.0,
    n_cds_exons: int = 1,
) -> dict:
    """Evaluate Ensembl 115 promoter, termination-readthrough, and UTR-intron spillover geometry."""
    mechanisms: list[str] = []
    bidir_tpm = float(bidirectional_partner_tpm) if pd.notna(bidirectional_partner_tpm) else 0.0
    if (
        pd.notna(bidirectional_promoter_bp)
        and math.isfinite(float(bidirectional_promoter_bp))
        and 0 <= float(bidirectional_promoter_bp) <= BIDIRECTIONAL_PROMOTER_MAX_BP
        and math.isfinite(bidir_tpm)
        and bidir_tpm >= BIDIRECTIONAL_PARTNER_MIN_TPM
    ):
        mechanisms.append("bidirectional_promoter_prompt")
    ss_tpm = float(same_strand_partner_tpm) if pd.notna(same_strand_partner_tpm) else 0.0
    if (
        pd.notna(same_strand_upstream_bp)
        and math.isfinite(float(same_strand_upstream_bp))
        and -5000 <= float(same_strand_upstream_bp) <= SAME_STRAND_UPSTREAM_MAX_BP
        and math.isfinite(ss_tpm)
        and ss_tpm >= SAME_STRAND_PARTNER_MIN_TPM
    ):
        mechanisms.append("same_strand_upstream_readthrough")
    nc_span = float(noncoding_5prime_span_bp) if pd.notna(noncoding_5prime_span_bp) else 0.0
    c_ratio = float(cds_span_ratio) if pd.notna(cds_span_ratio) else 1.0
    if (
        math.isfinite(nc_span)
        and math.isfinite(c_ratio)
        and nc_span >= NONCODING_5PRIME_SPAN_MIN_BP
        and c_ratio <= CDS_TO_SPAN_MAX_RATIO
    ):
        mechanisms.append("long_noncoding_5prime_intron_span")

    intronless_cds = int(n_cds_exons) == 1
    accepted_cds_assays = (
        ("cds_exonic_reads",) if intronless_cds else ("spliced_exon_junction", "cds_exonic_reads")
    )
    return {
        "transcriptional_spillover_risk": bool(mechanisms),
        "spillover_mechanisms": tuple(mechanisms),
        "is_intronless_cds": intronless_cds,
        "accepted_cds_confirmation_assays": accepted_cds_assays,
    }


def evaluate_nuclear_rare_cell_disposition(
    *,
    baseline_rare_cell_status: str,
    nucleus_only_signal: bool,
    finding_bulk_tpm: float,
    somatic_ihc_positive: bool,
    literature_effect: str,
    other_unsafe_gates: bool,
    transcriptional_spillover_risk: bool,
    paralog_ambiguous: bool,
    host_overlapped: bool,
    attribution_review_status: str = "",
    closed_repeat_family_cleared: bool = False,
) -> dict:
    """Determine whether an snRNA-seq rare-cell finding clears to ``nuclear_only`` (#765)."""
    lit = str(literature_effect or "not_reviewed").strip()
    if str(attribution_review_status or "").strip() == "confounder_rejected":
        return {
            "cleared": False,
            "disposition": "retained_established_cytoplasmic_or_bulk",
            "effective_rare_cell_status": "established",
            "effective_rare_cell_action": "unsafe",
            "nuclear_rna_caution": "",
        }
    if somatic_ihc_positive or lit in {"exclude", "confirm_detection"} or other_unsafe_gates:
        return {
            "cleared": False,
            "disposition": "retained_literature_or_ihc_supported",
            "effective_rare_cell_status": baseline_rare_cell_status,
            "effective_rare_cell_action": (
                "unsafe" if baseline_rare_cell_status == "established" else "review"
            ),
            "nuclear_rna_caution": "",
        }
    bulk = float("nan") if finding_bulk_tpm is None else float(finding_bulk_tpm)
    if (
        not nucleus_only_signal
        or not math.isfinite(bulk)
        or bulk < 0.0
        or bulk >= NUCLEAR_CLEARANCE_BULK_CEILING_TPM
    ):
        return {
            "cleared": False,
            "disposition": "retained_established_cytoplasmic_or_bulk",
            "effective_rare_cell_status": baseline_rare_cell_status,
            "effective_rare_cell_action": (
                "unsafe" if baseline_rare_cell_status == "established" else "review"
            ),
            "nuclear_rna_caution": "",
        }
    effective_paralog_ambig = bool(paralog_ambiguous) and not bool(closed_repeat_family_cleared)
    if effective_paralog_ambig or bool(host_overlapped):
        return {
            "cleared": False,
            "disposition": "retained_unresolved_paralog_or_exonic_host",
            "effective_rare_cell_status": "unresolved",
            "effective_rare_cell_action": "review",
            "nuclear_rna_caution": "snrnaseq_nucleus_only_pending_locus_attribution",
        }
    if not (bool(transcriptional_spillover_risk) or bool(closed_repeat_family_cleared)):
        return {
            "cleared": False,
            "disposition": "retained_nuclear_no_spillover_geometry",
            "effective_rare_cell_status": baseline_rare_cell_status,
            "effective_rare_cell_action": (
                "unsafe" if baseline_rare_cell_status == "established" else "review"
            ),
            "nuclear_rna_caution": "",
        }
    return {
        "cleared": True,
        "disposition": "nuclear_spillover_cleared",
        "effective_rare_cell_status": "nuclear_only",
        "effective_rare_cell_action": "caution",
        "nuclear_rna_caution": "snrnaseq_intronic_spillover_cleared",
    }


def census_collection_compartments() -> pd.DataFrame:
    """Return SHA-256-verified ``assay_compartment`` annotations across all 93 Census collections (#765)."""
    _verify_csv_receipt("census-collection-compartments")
    return get_data("census-collection-compartments").copy()


def nuclear_spillover_adjudications(*, only_cleared: bool = False) -> pd.DataFrame:
    """Return SHA-256-verified Ensembl 115 spillover & compartment adjudications (#765)."""
    _verify_csv_receipt("nuclear-spillover-adjudications")
    df = get_data("nuclear-spillover-adjudications").copy()
    if only_cleared:
        df = df[df["nuclear_spillover_cleared"].astype(bool)].reset_index(drop=True)
    return df


def cleared_candidate_comparison_observations() -> dict:
    """Return verified native-TPM normal, reproductive, and tumor observations for cleared Unsafe candidates."""
    payload = json.loads(
        files("oncoref")
        .joinpath("data/cleared-candidate-comparison-observations.json")
        .read_text(encoding="utf-8")
    )
    if payload.get("schema") != "oncoref.cleared-candidate-comparison-observations.v1":
        raise ValueError("Unexpected cleared-candidate comparison observations schema")
    return payload


def nuclear_spillover_receipt() -> dict:
    """Return the SHA-256 provenance receipt for #765 nuclear spillover QC."""
    return _verify_csv_receipt("nuclear-spillover-adjudications")


__all__ = [
    "ASSAY_COMPARTMENTS",
    "BIDIRECTIONAL_PARTNER_MIN_TPM",
    "BIDIRECTIONAL_PROMOTER_MAX_BP",
    "CDS_TO_SPAN_MAX_RATIO",
    "NONCODING_5PRIME_SPAN_MIN_BP",
    "NUCLEAR_CLEARANCE_BULK_CEILING_TPM",
    "NUCLEAR_DISPOSITIONS_765",
    "NUCLEAR_SPILLOVER_POLICY_ID",
    "SAME_STRAND_PARTNER_MIN_TPM",
    "SAME_STRAND_UPSTREAM_MAX_BP",
    "assess_transcriptional_spillover_geometry",
    "census_collection_compartments",
    "cleared_candidate_comparison_observations",
    "evaluate_nuclear_rare_cell_disposition",
    "nuclear_spillover_adjudications",
    "nuclear_spillover_receipt",
]
