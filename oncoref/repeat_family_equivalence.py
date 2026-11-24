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

"""Closed homologous repeat-family sequence equivalence & symmetric quantification (#767).

Resolves symmetric multimapping across closed germline/cleavage repeat arrays (such as
``DUX4`` and the ``4q35``/``10q26`` ``D4Z4`` macrosatellite ``DUX4L*`` copies, ``PRAMEF*``,
``CT47*``, ``TSPY*``, ``SPANX*``, and ``CSAG*``) without opening loopholes for heterologous
somatic paralogs:

1. **Sequence-partner partitioning**:
   Partitions primary-assembly sequence partners of a candidate protein group ``g`` into:
   - ``P_self(g)``: closed self-family homologs / ``unprocessed_pseudogene`` repeat copies
     sharing ``>= 85%`` ``CDS`` 31-mers or ``>= 90%`` coding amino-acid sequence with ``g``
     and exhibiting zero ``Unsafe`` somatic expression.
   - ``P_ext(g)``: external partners — any ``Unsafe`` coding paralog (``unsafe_sequence_overlap_pct
     >= 5.0%``, e.g., ``PSG8 <-> PSG1``, ``POTEC <-> POTED/E``, ``GAGE2E <-> GAGE12F/G``),
     somatic ``processed_pseudogene`` retrocopy, or somatic ``lncRNA``/chimeric host.

2. **Symmetric family-summed quantification ``(N_fam, R_fam, T_fam)``**:
   Clears ``paralog_ambiguous`` only when ``P_self(g)`` is non-empty and ``P_ext(g)`` is
   empty, then evaluates ``(N_fam, R_fam, T_fam)`` summed symmetrically across
   ``{g} U P_self(g)`` on locus/macrosatellite-resolved genome alignments (GTEx v11,
   HPA genome-aligned ``rna_tissue_hpa``, preimplantation embryo atlases, and linear-TPM
   tumor cohorts) so any leaky somatic repeat copy inflates ``N_fam`` equally.
   For ``DUX4``, clearing ``paralog_ambiguous`` interacts cleanly with ``#765`` (nuclear
   brain ``snRNA-seq`` gating): neither ``#765`` alone nor ``#767`` alone moves ``DUX4``
   out of ``Other:Review``, whereas joint ``#765 + #767`` clears both holds.

3. **Downstream coding closed-family contract for ``tsarina``**:
   Exposes ``closed_homolog_family_id``, ``co_silenced_coding_sibling_ids``, and
   ``closed_self_paralog_ids`` across ``target_selection.results()``, ``antigen_sets``,
   ``sequence_bundle.target_selection_dispositions()``, and ``evaluate_human_catalog()``
   so downstream peptide subtraction (``cta_exclusive_peptides()``) does not wipe out
   valid family epitopes against germline-restricted ``RA only`` / ``CA only`` coding
   siblings (e.g., ``CT47C1``, ``TSPY9``, ``SPANXN2/5``, ``CSAG2/3``) while strictly
   preserving subtraction against ``Unsafe`` or ``Other:Review`` somatic paralogs.
"""

from __future__ import annotations

import hashlib
import json
import math
from copy import deepcopy
from importlib.resources import files

import pandas as pd

from .load_dataset import get_data

REPEAT_FAMILY_POLICY_ID = "oncoref.repeat-family-equivalence.2026-11-23.1"

MIN_SELF_CDS_31MER_RECALL = 0.85
MIN_SELF_AA_IDENTITY = 0.90
MAX_CLOSED_UNSAFE_OVERLAP_PCT = 5.0
SOMATIC_FAMILY_NORMAL_CEILING_TPM = 5.0

REPEAT_FAMILY_DISPOSITIONS_767: tuple[str, ...] = (
    "closed_repeat_family_cleared",
    "retained_external_unsafe_or_somatic_paralog",
    "retained_processed_pseudogene_or_host_partner",
    "retained_family_somatic_ceiling_exceeded",
    "retained_insufficient_self_homology",
)

REPEAT_FAMILY_COLUMNS: tuple[str, ...] = (
    "closed_homolog_family_id",
    "co_silenced_coding_sibling_ids",
    "closed_self_paralog_ids",
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


def partition_sequence_partners(
    *,
    partners: list[dict],
    unsafe_sequence_overlap_pct: float = 0.0,
    family_somatic_max_tpm: float = 0.0,
) -> dict:
    """Partition sequence partners into ``P_self`` and ``P_ext`` under ``#767``.

    Each partner dict carries ``gene_id``, ``symbol``, ``biotype``, ``cds_31mer_recall``
    (or ``aa_identity``), and ``is_unsafe`` (True when the partner is ``Unsafe`` or
    somatically expressed above ``5.0 TPM``).
    """
    p_self: list[str] = []
    p_ext: list[str] = []
    has_processed_or_host = False

    for p in partners:
        gid = str(p.get("gene_id") or p.get("symbol") or "")
        biotype = str(p.get("biotype") or "")
        raw_kmer = p.get("cds_31mer_recall", 0.0)
        raw_aa = p.get("aa_identity", 0.0)
        kmer_rec = float(raw_kmer) if pd.notna(raw_kmer) else 0.0
        aa_id = float(raw_aa) if pd.notna(raw_aa) else 0.0
        if not math.isfinite(kmer_rec):
            kmer_rec = 0.0
        if not math.isfinite(aa_id):
            aa_id = 0.0
        is_unsafe = bool(p.get("is_unsafe", False))

        if biotype not in {"protein_coding", "unprocessed_pseudogene"}:
            p_ext.append(gid)
            has_processed_or_host = True
            continue
        if is_unsafe:
            p_ext.append(gid)
            continue
        if kmer_rec >= MIN_SELF_CDS_31MER_RECALL or aa_id >= MIN_SELF_AA_IDENTITY:
            p_self.append(gid)
        else:
            p_ext.append(gid)

    unsafe_ov = (
        float("nan") if unsafe_sequence_overlap_pct is None else float(unsafe_sequence_overlap_pct)
    )
    if (
        not math.isfinite(unsafe_ov)
        or unsafe_ov < 0.0
        or unsafe_ov >= MAX_CLOSED_UNSAFE_OVERLAP_PCT
    ):
        return {
            "closed_repeat_family_cleared": False,
            "disposition": "retained_external_unsafe_or_somatic_paralog",
            "p_self": tuple(sorted(set(p_self))),
            "p_ext": tuple(sorted(set(p_ext))),
        }
    if p_ext:
        disp = (
            "retained_processed_pseudogene_or_host_partner"
            if has_processed_or_host
            else "retained_external_unsafe_or_somatic_paralog"
        )
        return {
            "closed_repeat_family_cleared": False,
            "disposition": disp,
            "p_self": tuple(sorted(set(p_self))),
            "p_ext": tuple(sorted(set(p_ext))),
        }
    if not p_self:
        return {
            "closed_repeat_family_cleared": False,
            "disposition": "retained_insufficient_self_homology",
            "p_self": (),
            "p_ext": (),
        }
    som_val = float("nan") if family_somatic_max_tpm is None else float(family_somatic_max_tpm)
    if not math.isfinite(som_val) or som_val < 0.0 or som_val >= SOMATIC_FAMILY_NORMAL_CEILING_TPM:
        return {
            "closed_repeat_family_cleared": False,
            "disposition": "retained_family_somatic_ceiling_exceeded",
            "p_self": tuple(sorted(set(p_self))),
            "p_ext": (),
        }
    return {
        "closed_repeat_family_cleared": True,
        "disposition": "closed_repeat_family_cleared",
        "p_self": tuple(sorted(set(p_self))),
        "p_ext": (),
    }


def repeat_family_adjudications(*, only_cleared: bool = False) -> pd.DataFrame:
    """Return SHA-256-verified closed repeat-family adjudications across ``paralog_ambiguous`` candidates."""
    _verify_csv_receipt("repeat-family-adjudications")
    df = get_data("repeat-family-adjudications").copy()
    if only_cleared:
        df = df[df["closed_repeat_family_cleared"].astype(bool)].reset_index(drop=True)
    return df


def closed_homolog_families() -> pd.DataFrame:
    """Return SHA-256-verified closed coding & repeat family memberships for downstream subtraction (#767)."""
    _verify_csv_receipt("closed-homolog-families")
    return get_data("closed-homolog-families").copy()


def repeat_family_comparison_observations() -> dict:
    """Return symmetric family-summed ``(N_fam, R_fam, T_fam)`` observations for cleared repeat families."""
    payload = json.loads(
        files("oncoref")
        .joinpath("data/repeat-family-comparison-observations.json")
        .read_text(encoding="utf-8")
    )
    if payload.get("schema") != "oncoref.repeat-family-comparison-observations.v1":
        raise ValueError("Unexpected repeat-family comparison observations schema")
    return deepcopy(payload)


def repeat_family_receipt() -> dict:
    """Return the SHA-256 provenance receipt for ``#767`` repeat-family adjudications."""
    return _verify_csv_receipt("repeat-family-adjudications")


def annotate_catalog(frame: pd.DataFrame, *, enabled: bool = True) -> pd.DataFrame:
    """Attach ``tsarina`` closed-homolog-family columns to a catalog DataFrame."""
    out = frame.copy()
    if not enabled:
        out["closed_homolog_family_id"] = ""
        out["co_silenced_coding_sibling_ids"] = "[]"
        out["closed_self_paralog_ids"] = "[]"
        return out

    fam_df = closed_homolog_families()
    by_eid: dict[str, tuple[str, str, str]] = {}
    for r in fam_df.itertuples(index=False):
        by_eid[str(r.entity_id)] = (
            str(r.closed_homolog_family_id) if pd.notna(r.closed_homolog_family_id) else "",
            str(r.co_silenced_coding_sibling_ids)
            if pd.notna(r.co_silenced_coding_sibling_ids)
            else "[]",
            str(r.closed_self_paralog_ids) if pd.notna(r.closed_self_paralog_ids) else "[]",
        )

    fids, sibs, pars = [], [], []
    for eid in out["entity_id"].astype(str):
        fid, sib_json, par_json = by_eid.get(eid, ("", "[]", "[]"))
        fids.append(fid)
        sibs.append(sib_json)
        pars.append(par_json)

    out["closed_homolog_family_id"] = fids
    out["co_silenced_coding_sibling_ids"] = sibs
    out["closed_self_paralog_ids"] = pars
    return out


__all__ = [
    "MAX_CLOSED_UNSAFE_OVERLAP_PCT",
    "MIN_SELF_AA_IDENTITY",
    "MIN_SELF_CDS_31MER_RECALL",
    "REPEAT_FAMILY_COLUMNS",
    "REPEAT_FAMILY_DISPOSITIONS_767",
    "REPEAT_FAMILY_POLICY_ID",
    "SOMATIC_FAMILY_NORMAL_CEILING_TPM",
    "annotate_catalog",
    "closed_homolog_families",
    "partition_sequence_partners",
    "repeat_family_adjudications",
    "repeat_family_comparison_observations",
    "repeat_family_receipt",
]
