"""Individual reviews of the nonreproductive tumor-expression screen."""

from __future__ import annotations

import gzip
import json
from importlib.resources import files

import pandas as pd


def audit() -> dict:
    """Return the dated selection definition, counts and all 47 source-linked reviews.

    The set is selected by zero recorded adult reproductive bulk RNA, not proven
    absence in every normal tissue. It contains no confirmed Cancer-Only member.
    ``classification`` includes accepted pending normal-evidence curation;
    ``released_classification`` identifies its explicitly named baseline snapshot.
    The audit does not add a fourth classification category.
    """
    path = files("oncoref").joinpath("data/cancer-only-audit.json.gz")
    data = json.loads(gzip.decompress(path.read_bytes()))
    if data.get("schema") != "oncoref.cancer-only-audit" or data.get("version") != 1:
        raise ValueError("Unsupported cancer-only audit")
    return data


def records(gene: str | None = None) -> pd.DataFrame:
    """One row per reviewed protein; exact current symbol or Ensembl gene ID."""
    result = pd.DataFrame(audit()["records"])
    if gene is not None:
        if not isinstance(gene, str) or not gene:
            raise ValueError("gene must be a nonempty exact symbol or Ensembl gene ID")
        result = result[result.protein.eq(gene) | result.gene_id.eq(gene)]
    return result.reset_index(drop=True)
