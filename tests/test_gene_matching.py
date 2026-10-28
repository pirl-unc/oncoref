import json

import numpy as np
import pandas as pd
import pytest

from oncoref import gene_ids
from oncoref.gene_matching import (
    attach_gene_identities,
    collapse_reference_expression,
    gene_identity_catalog,
    harmonize_reference_expression,
    match_reference_genes,
)


def match(frame):
    return match_reference_genes(
        frame, id_column="id", symbol_column="name", source="fixture", source_version="1"
    )


def test_stable_ids_names_and_ambiguity_are_not_interchangeable():
    data = pd.DataFrame(
        {
            "id": [
                "ENSG00000181222.5",
                "ENSG00000287363",
                "ENSG00000285437",
                "",
                "ENSG00000148362",
                "ENSG99999999999",
                "ENSG00000181222",
            ],
            "name": ["POLR2A", "GUCA1ANB", "POLR2J3", "POLR2J3", "PAXX", "TP53", "TP53"],
        },
        index=[9, 3, 3, 2, 8, 4, 1],
    )
    result = match(data)
    assert result.index.equals(data.index)
    assert result.identity_status.tolist() == ["matched"] * 3 + [
        "ambiguous_name",
        "matched",
        "unmapped_id",
        "id_name_conflict",
    ]
    assert result.canonical_gene_id.iloc[4] == "ENSG00000310560"
    assert json.loads(result.historical_ensembl_gene_ids.iloc[4]) == ["ENSG00000148362"]
    assert len(json.loads(result.candidate_gene_ids.iloc[3])) == 2
    assert not result.canonical_gene_id.iloc[5]
    assert gene_ids.canonical_gene_id("POLR2J3") is None


def test_unknown_source_label_does_not_erase_verified_id():
    result = match(pd.DataFrame({"id": ["ENSG00000181222"], "name": ["old-unregistered-label"]}))
    assert result.identity_status.iloc[0] == "matched"
    assert result.source_name_unrecognized.iloc[0]
    assert result.source_gene_name.iloc[0] == "old-unregistered-label"


def test_missing_name_sentinels_cannot_become_scn11a_synonyms():
    result = match(
        pd.DataFrame(
            {"id": ["ENSG00000181222", "", "ENSG00000168356"], "name": ["NaN", "NaN", "NaN"]}
        )
    )
    assert result.identity_status.tolist() == ["matched", "unmapped", "matched"]
    assert result.canonical_gene_id.tolist() == ["ENSG00000181222", "", "ENSG00000168356"]
    assert result.source_gene_name.tolist() == ["NaN"] * 3


def test_symbol_only_and_empty_inputs():
    result = match(pd.DataFrame({"id": [None, ""], "name": ["GNB2L1", "TP53"]}))
    assert result.canonical_gene_id.tolist() == ["ENSG00000204628", "ENSG00000141510"]
    assert match(pd.DataFrame({"id": [], "name": []})).empty


def test_identity_bundle_survives_stage_rows():
    data = pd.DataFrame(
        {"gene_id": ["ENSG00000310560"] * 2, "stage": ["bulk", "tumor"]}, index=[8, 4]
    )
    result = attach_gene_identities(data)
    assert result.index.equals(data.index)
    assert result.historical_ensembl_gene_ids.nunique() == 1
    assert "PAXX" in json.loads(result.gene_name_synonyms.iloc[0])
    with pytest.raises(ValueError, match="outside"):
        attach_gene_identities(pd.DataFrame({"gene_id": ["ENSG99999999999"]}))


def test_alias_duplicates_never_become_extra_expression():
    data = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG00000148362", "ENSG00000310560"],
            "Symbol": ["PAXX", "PAXX"],
            "sample": [3.0, 3.0],
        }
    )
    kw = {"sample_cols": ["sample"], "source": "fixture", "source_version": "1"}
    with pytest.raises(ValueError, match="Conflicting/repeated"):
        harmonize_reference_expression(data, **kw)
    result, audit = harmonize_reference_expression(data, duplicate_policy="equal_alias_rows", **kw)
    assert result["sample"].tolist() == [3.0]
    assert audit.identity_status.tolist() == ["matched", "duplicate_alias_not_added"]
    data.loc[1, "sample"] = 4
    with pytest.raises(ValueError, match="Conflicting/repeated"):
        harmonize_reference_expression(data, duplicate_policy="equal_alias_rows", **kw)
    quarantined, audit = harmonize_reference_expression(data, duplicate_policy="quarantine", **kw)
    assert quarantined.empty
    assert audit.identity_status.eq("conflicting_duplicate_measurements").all()


def test_same_sample_proteoform_sum_precedes_maximum_and_tracks_missingness():
    frame = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG00000268651", "ENSG00000184033"],
            "Symbol": ["CTAG1A", "CTAG1B"],
            "a": [1.2, 1.3],
            "b": [10, 0],
            "partial": [np.nan, 0.5],
            "absent": [np.nan, np.nan],
        }
    )
    kw = {
        "sample_cols": ["a", "b", "partial", "absent"],
        "source": "fixture",
        "source_version": "1",
        "unit": "TPM",
    }
    result, coverage, audit = collapse_reference_expression(frame, **kw)
    assert len(result) == 1 and result.proteoform_key.iloc[0] == "NY-ESO-1"
    assert result.a.iloc[0] == 2.5
    assert result.b.iloc[0] == 10  # never sum different samples' maxima (11.3)
    assert result.partial.iloc[0] == 0.5
    assert np.isnan(result.absent.iloc[0])
    assert coverage.complete.tolist() == [True, True, False, False]
    assert coverage.measured_members.tolist() == [2, 2, 1, 0]
    assert len(json.loads(result.member_identities.iloc[0])) == 2
    assert audit.identity_status.eq("matched").all()
    swapped = frame.copy()
    swapped[kw["sample_cols"]] = frame[kw["sample_cols"]].iloc[::-1].to_numpy()
    other, _, _ = collapse_reference_expression(swapped, **kw)
    pd.testing.assert_frame_equal(result, other)


def test_missing_locus_is_not_an_observed_zero():
    frame = pd.DataFrame({"Ensembl_Gene_ID": ["ENSG00000204379"], "Symbol": ["XAGE1A"], "s": [0.0]})
    result, coverage, _ = collapse_reference_expression(
        frame, sample_cols=["s"], source="fixture", source_version="1", unit="TPM"
    )
    assert result.s.iloc[0] == 0
    assert not coverage.complete.iloc[0]
    assert json.loads(result.missing_member_gene_ids.iloc[0]) == ["ENSG00000204382"]


def test_catalog_is_defensive_and_sources_and_units_are_explicit():
    first = gene_identity_catalog()
    first.loc[0, "canonical_symbol"] = "changed"
    assert gene_identity_catalog().canonical_symbol.iloc[0] != "changed"
    with pytest.raises(ValueError, match="taxon"):
        match_reference_genes(
            pd.DataFrame(), id_column=None, source="fixture", source_version="1", taxonomy_id=10090
        )
    with pytest.raises(ValueError, match="unit"):
        collapse_reference_expression(
            pd.DataFrame(), sample_cols=["s"], source="fixture", source_version="1", unit="logTPM"
        )
