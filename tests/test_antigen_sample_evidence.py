"""Individual-sample grouping must preserve locus coverage and native units."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from antigen_sample_evidence import sum_transcript_chunks, tissue_samples

from oncoref.target_selection import developmental_donor_requirement


def test_transcript_sum_across_chunks_preserves_missing_members():
    chunks = [
        pd.DataFrame({"ensgid": ["g1", "g2"], "TPM.testis.1": [2, 4], "TPM.liver.2": [np.nan, 1]}),
        pd.DataFrame({"ensgid": ["g1"], "TPM.testis.1": [3], "TPM.liver.2": [5]}),
    ]
    gene = sum_transcript_chunks(chunks)
    assert gene.loc["g1", "TPM.testis.1"] == 5
    assert pd.isna(gene.loc["g1", "TPM.liver.2"])
    facts = pd.DataFrame(
        {"entity_id": ["complete", "missing"], "member_gene_ids": ['["g1", "g2"]', '["g1", "g3"]']}
    )
    ids = pd.DataFrame({"source_gene_id": ["g1", "g2"], "canonical_gene_id": ["g1", "g2"]})
    d = tissue_samples(facts, ids, genes=gene).set_index(["entity_id", "sample_id"])
    assert d.loc[("complete", "testis.1"), "value"] == 9
    assert pd.isna(d.loc[("complete", "liver.2"), "value"])
    assert not d.loc["missing", "identity_resolved"].any()
    assert not d.loc["missing", "qc_pass"].any()
    assert d.participant_id.eq("").all()  # Sample identifiers are not donor IDs.


@pytest.mark.parametrize("value", [-1, float("inf")])
def test_invalid_native_tpm_is_rejected(value):
    with pytest.raises(ValueError, match="finite and nonnegative"):
        sum_transcript_chunks([pd.DataFrame({"ensgid": ["g1"], "TPM.testis.1": [value]})])


@pytest.mark.parametrize(
    "cells,positive,required",
    [
        (8, 4, 2),
        (8, 3, 0),
        (5, 3, 2),
        (4, 4, 0),
        (8, 8, 2),
        (20, 5, 1),
        (60, 5, 0),
        (60, 6, 1),
    ],
)
def test_cell_support_sets_the_required_donor_replication(cells, positive, required):
    assert developmental_donor_requirement(cells, positive) == required


def test_mammary_state_is_inferred_and_incomplete_marker_panel_is_unknown(monkeypatch):
    from antigen_sample_evidence import LACTATION_MARKERS, mammary_sample_states

    from oncoref import gene_matching

    ids = pd.DataFrame(
        {"canonical_gene_id": ["m1", "m2", "m3", "m4"], "canonical_symbol": LACTATION_MARKERS}
    )
    monkeypatch.setattr(gene_matching, "gene_identity_catalog", lambda: ids)
    accepted = ids[["canonical_gene_id"]].assign(source_gene_id=ids.canonical_gene_id)
    values = pd.DataFrame(
        {
            "TPM.breast.1": [100, 100, 100, 0],
            "TPM.breast.2": [0, 0, 0, 0],
            "TPM.breast.3": [0, 0, 0, np.nan],
        },
        index=ids.canonical_gene_id,
    )
    states = mammary_sample_states(values, accepted)
    assert [states[f"breast.{i}"]["state"] for i in (1, 2, 3)] == [
        "secretory_signature",
        "no_high_secretory_signature",
        "unknown",
    ]
    assert not any(v["provider_confirmed"] for v in states.values())
