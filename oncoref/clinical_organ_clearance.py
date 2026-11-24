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

"""Clinical pMHC-directed organ-safety evidence tier & Xq28 CTAG2 attribution (#766).

Validates and applies epitope-, HLA-, and Ensembl-115-isoform-scoped human clinical
organ-safety clearances from high-potency pMHC-directed therapies (``TCR-T`` /
``ImmTAC``) anchored to ``oncoref.tcr_trials`` and ``oncoref.reference_data``:

1. **Qualifying clinical modality & potency:** Requires adoptive ``TCR-T`` or
   ``ImmTAC`` targeting an exact ``(peptide, HLA)`` complex with published objective
   clinical tumor responses and ``>= 20`` reviewed HLA-matched infused patients.
   Peptide/mRNA vaccines cannot supply organ-safety clearance.
2. **Ensembl 115 isoform & epitope scope (``CTAG2`` ``LAGE-1a`` vs. ``LAGE-1b``):**
   Verifies exact coding-transcript sequence containment. For ``CTAG2``
   (``ENSG00000126890``), ``SLLMWITQC`` is present in ``ENST00000369585``
   (``LAGE-1a``, 180 aa) at positions 157–165 and **absent** from
   ``ENST00000247306`` (``LAGE-1b``, 210 aa); clearance is explicitly scoped to
   ``clinical_clearance_scope = "epitope"`` and ``clinical_cleared_transcripts =
   ["ENST00000369585"]``.
3. **Vascular perfusion & blood–brain barrier guardrails:** Systemic intravenous
   infusion clears only vascularly perfused peripheral organs (e.g., ``lung``,
   ``urinary bladder``) with zero attributed severe specificity toxicity, and fails
   closed on barrier-protected compartments (``brain``, ``spinal cord``, ``eye``,
   ``retina``, ``testis``) unless direct compartment T-cell exposure is proven.
4. **Negative-clearance clinical toxicity controls:** Receptor/target records with
   suspected, confirmed, or unresolved specificity toxicity (``magea3_a3a`` cardiac,
   ``magea3_a2_neuro`` brain, ``tbi1301`` unresolved lung injury) fail closed and can
   never supply organ clearance.
"""

from __future__ import annotations

import json
import re
from copy import deepcopy
from importlib.resources import files

import pandas as pd

from . import reference_data, tcr_trials

CLINICAL_ORGAN_CLEARANCE_POLICY_ID = "oncoref.clinical-organ-clearance.2026-11-23.1"
CLEARANCE_FILE = "clinical-organ-clearances.json"

ALLOWED_MODALITIES: frozenset[str] = frozenset({"TCR-T", "ImmTAC"})
MIN_HLA_MATCHED_PATIENTS = 20

VASCULARLY_PERFUSED_ORGANS: frozenset[str] = frozenset(
    {
        "lung",
        "urinary bladder",
        "kidney",
        "liver",
        "colon",
        "small intestine",
        "stomach",
        "esophagus",
        "skin",
        "bone marrow",
        "heart",
    }
)

BARRIER_PROTECTED_ORGANS: frozenset[str] = frozenset(
    {
        "brain",
        "spinal cord",
        "eye",
        "retina",
        "testis",
    }
)

CLINICAL_CLEARANCE_COLUMNS: tuple[str, ...] = (
    "clinical_clearance_scope",
    "clinical_cleared_organs",
    "clinical_cleared_epitopes",
    "clinical_hla_alleles",
    "clinical_cleared_transcripts",
    "clinical_unverified_isoforms",
    "clinical_trial_source_ids",
)


def validate(record: dict) -> None:
    """Validate clinical organ-clearance ledger against ``tcr_trials`` & perfusion rules."""
    if (record.get("schema"), record.get("schema_version")) != (
        "oncoref.clinical-organ-clearances",
        1,
    ):
        raise ValueError("Unsupported clinical-organ-clearances schema")

    tcr_doc = tcr_trials.evidence()
    products_by_id = {p["product_id"]: p for p in tcr_doc["products"]}
    cohorts_by_id = {c["cohort_id"]: c for c in tcr_doc["cohorts"]}
    events_by_product: dict[str, list[dict]] = {}
    for ev in tcr_doc["events"]:
        events_by_product.setdefault(ev["product_id"], []).append(ev)
    sources_by_key = tcr_doc["sources"]

    seen_ids: set[str] = set()
    for entry in record.get("entries", []):
        cid = entry.get("clearance_id", "")
        if not cid or cid in seen_ids:
            raise ValueError(f"Duplicate or missing clearance_id: {cid!r}")
        seen_ids.add(cid)

        status = entry.get("status")
        if status not in {"cleared", "blocked_clinical_toxicity", "blocked_unresolved_event"}:
            raise ValueError(f"Invalid clearance status {status!r} on {cid!r}")

        pid = entry.get("product_id")
        if pid not in products_by_id:
            raise ValueError(f"Unknown tcr_trials product_id {pid!r} on {cid!r}")
        prod = products_by_id[pid]

        modality = entry.get("modality")
        if modality not in ALLOWED_MODALITIES or modality != prod["modality"]:
            raise ValueError(
                f"Clinical organ clearance requires {sorted(ALLOWED_MODALITIES)}, got {modality!r}"
            )

        for skey in entry.get("source_keys", []):
            if skey not in sources_by_key:
                raise ValueError(f"Unknown tcr_trials source_key {skey!r}")
            src_info = sources_by_key[skey]
            sid = src_info["source_id"]
            ver = src_info["version"]
            if sid not in reference_data.REFERENCE_SOURCES:
                raise ValueError(f"Unpinned reference source_id {sid!r}")
            if reference_data.REFERENCE_SOURCES[sid]["pins"][ver]["sha256"] != src_info["sha256"]:
                raise ValueError(f"SHA-256 mismatch on {sid}:{ver}")

        if status != "cleared":
            if entry.get("cleared_organs"):
                raise ValueError(f"Blocked entry {cid!r} cannot list cleared_organs")
            continue

        if entry.get("clearance_scope") != "epitope":
            raise ValueError(
                f"Cleared entry {cid!r} must declare clearance_scope='epitope' (never blanket gene)"
            )
        if not bool(entry.get("objective_tumor_responses", False)):
            raise ValueError(f"Cleared entry {cid!r} requires objective_tumor_responses=True")

        cohort_ids = entry.get("cohort_ids", [])
        if not cohort_ids:
            raise ValueError(f"Cleared entry {cid!r} requires reviewed cohort_ids")
        total_patients = 0
        for coh_id in cohort_ids:
            if coh_id not in cohorts_by_id:
                raise ValueError(f"Unknown cohort_id {coh_id!r}")
            coh = cohorts_by_id[coh_id]
            if coh["product_id"] != pid or coh["review_status"] != "reviewed":
                raise ValueError(f"Cohort {coh_id!r} not reviewed for product {pid!r}")
            if coh["grade45_specificity"] != "none_identified":
                raise ValueError(
                    f"Cohort {coh_id!r} has specificity finding {coh['grade45_specificity']!r}"
                )
            total_patients += int(coh["patients_treated"] or 0)
        if total_patients < MIN_HLA_MATCHED_PATIENTS:
            raise ValueError(
                f"Cleared entry {cid!r} has {total_patients} patients < {MIN_HLA_MATCHED_PATIENTS}"
            )
        if int(entry.get("hla_matched_patients_treated", 0)) != total_patients:
            raise ValueError(f"Denominator mismatch on {cid!r}")

        # Verify target gene & peptide match in tcr-trial-evidence
        gid = entry.get("gene_id")
        epitopes = entry.get("cleared_epitopes", [])
        if not epitopes:
            raise ValueError(f"Cleared entry {cid!r} requires cleared_epitopes")
        prod_targets = [t for t in prod["targets"] if t["gene_id"] == gid]
        if not prod_targets:
            raise ValueError(f"Gene {gid!r} is not a verified target of {pid!r}")
        if set(epitopes) != {t["peptide"] for t in prod_targets}:
            raise ValueError(f"Epitope mismatch between {cid!r} and tcr-trial-evidence")

        hlas = entry.get("hla_alleles", [])
        if not hlas or not (set(hlas) <= set(prod.get("hla", []))):
            raise ValueError(
                f"Cleared entry {cid!r} requires reviewed hla_alleles subset of product {pid!r} HLA"
            )

        txs = entry.get("cleared_transcripts", [])
        unc_txs = entry.get("uncleared_coding_transcripts", [])
        if not txs or any(not re.fullmatch(r"ENST\d+", str(t)) for t in [*txs, *unc_txs]):
            raise ValueError(
                f"Cleared entry {cid!r} requires canonical unversioned ENST transcript IDs"
            )
        overlap_isoforms = set(txs) & set(unc_txs)
        if overlap_isoforms:
            raise ValueError(f"Transcript cannot be both cleared and uncleared: {overlap_isoforms}")

        # Perfusion / barrier & product event checks
        barrier_ok = bool(entry.get("barrier_compartment_exposure_verified", False))
        prod_events = events_by_product.get(pid, [])
        # Any product with unresolved/suspected/confirmed specificity event anywhere blocks product credit
        for pev in prod_events:
            if pev["specificity_attribution"] in {
                "suspected",
                "confirmed",
                "unresolved",
                "unassessed",
            }:
                raise ValueError(
                    f"Product {pid!r} has {pev['specificity_attribution']!r} event {pev['event_id']!r}"
                )
        for organ in entry.get("cleared_organs", []):
            if organ in BARRIER_PROTECTED_ORGANS and not barrier_ok:
                raise ValueError(
                    f"Systemic IV exposure on {cid!r} cannot clear barrier-protected organ {organ!r}"
                )
            if organ not in VASCULARLY_PERFUSED_ORGANS and not barrier_ok:
                raise ValueError(f"Organ {organ!r} is not in VASCULARLY_PERFUSED_ORGANS")
            # Even a non-established fatal/confounded event in the same organ precludes claiming clearance there
            if any(pev.get("affected_tissue") == organ for pev in prod_events):
                raise ValueError(
                    f"Entry {cid!r} cannot claim organ clearance in {organ!r} with recorded event in that organ"
                )


def ledger() -> dict:
    """Return the validated clinical organ-clearance ledger (`#766`), offline."""
    payload = json.loads(
        files("oncoref").joinpath("data", CLEARANCE_FILE).read_text(encoding="utf-8")
    )
    validate(payload)
    return deepcopy(payload)


def clearances_table(*, only_cleared: bool = False) -> pd.DataFrame:
    """Return one row per reviewed clinical organ-clearance or negative-control entry."""
    doc = ledger()
    rows = []
    for e in doc["entries"]:
        rows.append(
            {
                "clearance_id": e["clearance_id"],
                "entity_id": e["entity_id"],
                "protein_group_label": e["protein_group_label"],
                "gene_id": e["gene_id"],
                "product_id": e["product_id"],
                "modality": e["modality"],
                "status": e["status"],
                "clearance_scope": e.get("clearance_scope", ""),
                "hla_matched_patients_treated": int(e.get("hla_matched_patients_treated", 0)),
                "clinical_cleared_organs": json.dumps(e.get("cleared_organs", [])),
                "clinical_cleared_epitopes": json.dumps(e.get("cleared_epitopes", [])),
                "clinical_hla_alleles": json.dumps(e.get("hla_alleles", [])),
                "clinical_cleared_transcripts": json.dumps(e.get("cleared_transcripts", [])),
                "uncleared_coding_transcripts": json.dumps(
                    e.get("uncleared_coding_transcripts", [])
                ),
                "clinical_trial_source_ids": json.dumps(e.get("trial_source_ids", [])),
                "rationale": e.get("rationale", ""),
            }
        )
    df = pd.DataFrame(rows)
    if only_cleared:
        df = df[df["status"].eq("cleared")].reset_index(drop=True)
    return df


def ctag2_xq28_attribution() -> dict:
    """Return the Ensembl 115 Xq28 locus & self-family paralog attribution for CTAG2 (#711/#766)."""
    doc = ledger()
    return deepcopy(doc["locus_attributions"]["CTAG2"])


def annotate_catalog(frame: pd.DataFrame, *, enabled: bool = True) -> pd.DataFrame:
    """Attach structured ``tsarina`` ``clinical_cleared_*`` columns to a catalog DataFrame."""
    out = frame.copy()
    if not enabled:
        out["clinical_clearance_scope"] = ""
        for col in CLINICAL_CLEARANCE_COLUMNS[1:]:
            out[col] = "[]"
        return out

    doc = ledger()
    by_eid: dict[str, dict] = {}
    for e in doc["entries"]:
        if e["status"] != "cleared":
            continue
        eid = e["entity_id"]
        cur = by_eid.setdefault(
            eid,
            {
                "scope": e["clearance_scope"],
                "organs": set(),
                "epitopes": set(),
                "hlas": set(),
                "transcripts": set(),
                "unverified": set(),
                "sources": set(),
            },
        )
        cur["organs"].update(e.get("cleared_organs", []))
        cur["epitopes"].update(e.get("cleared_epitopes", []))
        cur["hlas"].update(e.get("hla_alleles", []))
        cur["transcripts"].update(e.get("cleared_transcripts", []))
        cur["unverified"].update(e.get("uncleared_coding_transcripts", []))
        cur["sources"].update(e.get("trial_source_ids", []))

    scopes, organs, epitopes, hlas, txs, unv_txs, srcs = [], [], [], [], [], [], []
    for eid in out["entity_id"].astype(str):
        rec = by_eid.get(eid)
        if rec is None:
            scopes.append("")
            organs.append("[]")
            epitopes.append("[]")
            hlas.append("[]")
            txs.append("[]")
            unv_txs.append("[]")
            srcs.append("[]")
        else:
            scopes.append(rec["scope"])
            organs.append(json.dumps(sorted(rec["organs"])))
            epitopes.append(json.dumps(sorted(rec["epitopes"])))
            hlas.append(json.dumps(sorted(rec["hlas"])))
            txs.append(json.dumps(sorted(rec["transcripts"])))
            unv_txs.append(json.dumps(sorted(rec["unverified"] - rec["transcripts"])))
            srcs.append(json.dumps(sorted(rec["sources"])))

    out["clinical_clearance_scope"] = scopes
    out["clinical_cleared_organs"] = organs
    out["clinical_cleared_epitopes"] = epitopes
    out["clinical_hla_alleles"] = hlas
    out["clinical_cleared_transcripts"] = txs
    out["clinical_unverified_isoforms"] = unv_txs
    out["clinical_trial_source_ids"] = srcs
    return out


__all__ = [
    "ALLOWED_MODALITIES",
    "BARRIER_PROTECTED_ORGANS",
    "CLEARANCE_FILE",
    "CLINICAL_CLEARANCE_COLUMNS",
    "CLINICAL_ORGAN_CLEARANCE_POLICY_ID",
    "MIN_HLA_MATCHED_PATIENTS",
    "VASCULARLY_PERFUSED_ORGANS",
    "annotate_catalog",
    "clearances_table",
    "ctag2_xq28_attribution",
    "ledger",
    "validate",
]
