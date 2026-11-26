# Licensed under the Apache License, Version 2.0.
"""Reproductive-age female GTEx stratification & adult human oocyte/follicle RNA panel.

Resolves two systematic blind spots in bulk postmenopausal reproductive tissue atlases:

1. **GTEx v11 donor age stratification (ages 20–49 vs. 50–79):** Nearly 60% of GTEx
   ovary, uterus, and vagina biopsies derive from postmenopausal donors (ages 50–79)
   with follicularly depleted ovarian cortex and atrophic endometrium. Computing female
   ``HostContext`` organ exposures across reproductive-age ``XX`` donors (``20–49``)
   captures active premenopausal physiology without diluting cyclical/follicular signals.
2. **Adult human folliculogenesis RNA-seq (GSE107746, Zhang et al. 2018, PMID:30420494):**
   Even in premenopausal ovarian cortex, oocytes represent ``<0.1%`` of stromal tissue
   mass. Stage-replicated linear ``TPM`` across ``128`` QC-passing (``<= 2.5%`` ectopic ORs)
   adult human primordial, primary, secondary, and antral oocyte libraries plus
   primordial-to-preovulatory granulosa libraries (`9` cell-type × stage strata)
   measures genuine adult female germ-cell expression on unambiguous,
   non-paralog-confounded Ensembl 115 coding loci.
"""

from __future__ import annotations

import hashlib
import json
from functools import cache
from importlib.resources import files

import pandas as pd

from . import gtex
from .load_dataset import get_data

ADULT_FOLLICLE_PANEL_ID = "oncoref.human-adult-oocyte-follicle-rna.2026-11-26.1"


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


@cache
def _verified_follicle_payload() -> tuple[pd.DataFrame, dict]:
    data_root = files("oncoref").joinpath("data")
    csv_res = data_root.joinpath("human-adult-oocyte-follicle-rna.csv.gz")
    json_res = data_root.joinpath("human-adult-oocyte-follicle-rna.json")
    receipt = json.loads(json_res.read_text(encoding="utf-8"))
    if _sha256_bytes(csv_res.read_bytes()) != receipt["sha256"]:
        raise ValueError(
            "human-adult-oocyte-follicle-rna.csv.gz does not match its SHA-256 receipt"
        )
    df = get_data("human-adult-oocyte-follicle-rna").copy()
    return df, receipt


def adult_oocyte_follicle_rna() -> pd.DataFrame:
    """Return SHA-256-verified stage-replicated adult human oocyte & granulosa TPMs (GSE107746)."""
    df, _ = _verified_follicle_payload()
    return df.copy()


def adult_oocyte_follicle_receipt() -> dict:
    """Return SHA-256 provenance and library-QC metadata for the GSE107746 adult follicle panel."""
    _, receipt = _verified_follicle_payload()
    return json.loads(json.dumps(receipt))


def reproductive_age_female_donor_ids() -> frozenset[str]:
    """Return chromosomal ``XX`` GTEx v11 donor IDs in reproductive age brackets (``20–49``)."""
    age_df = gtex.donor_age(reproductive_age_only=True, sex="XX")
    return frozenset(age_df["donor_id"].astype(str))


__all__ = [
    "ADULT_FOLLICLE_PANEL_ID",
    "adult_oocyte_follicle_receipt",
    "adult_oocyte_follicle_rna",
    "reproductive_age_female_donor_ids",
]
