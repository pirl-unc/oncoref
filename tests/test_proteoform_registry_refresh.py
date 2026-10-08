# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

"""Registry refreshes must reach expression summaries, including cached shards (#567)."""

import importlib

import pandas as pd
import pytest

from oncoref import expression
from oncoref.expression_builders import (
    cohort_percentile_vectors,
    sum_proteoform_tpm,
    within_sample_top_fractions,
)
from oncoref.proteoforms import (
    _proteoform_registry_metadata,
    collapse_to_proteoforms,
    proteoform_group_map,
    proteoform_key,
    proteoform_symbol,
)

_RBMY_IDS = ["ENSG00000169800", "ENSG00000226941"]


@pytest.fixture
def rbmy_expression():
    # Allocated per-gene TPM contributions. RBMY sums to 12 in both samples;
    # HBA sums to 16 only in genome scope, changing the within-sample rank axis.
    return pd.DataFrame(
        {
            "Ensembl_Gene_ID": [
                *_RBMY_IDS,
                "ENSG00000185686",
                "ENSG00000206172",
                "ENSG00000188536",
            ],
            "Symbol": ["RBMY1F", "RBMY1J", "PRAME", "HBA1", "HBA2"],
            "s1": [6.0, 6.0, 11.0, 8.0, 8.0],
            "s2": [0.0, 12.0, 11.0, 8.0, 8.0],
        }
    )


def test_rbmy_identity_and_per_sample_collapse(rbmy_expression):
    assert {proteoform_key(gene_id) for gene_id in _RBMY_IDS} == {"RBMY1F/J"}
    collapsed = collapse_to_proteoforms(rbmy_expression).set_index("proteoform_key")
    assert not set(_RBMY_IDS).intersection(collapsed.index)
    assert collapsed.loc["RBMY1F/J", ["s1", "s2"]].tolist() == [12.0, 12.0]
    assert collapsed.loc["RBMY1F/J", "proteoform_members"] == "RBMY1F/RBMY1J"


@pytest.mark.parametrize("scope,expected_fraction", [("cta", 1.0), ("genome", 0.0)])
def test_rbmy_sums_before_ranking_on_each_scope_axis(rbmy_expression, scope, expected_fraction):
    per_gene = within_sample_top_fractions(rbmy_expression).set_index("Symbol")
    assert per_gene.loc["RBMY1F", "frac_samples_top5pct"] == 0.0
    assert per_gene.loc["RBMY1J", "frac_samples_top5pct"] == 0.5
    collapsed = collapse_to_proteoforms(rbmy_expression, scope=scope)
    ranked = within_sample_top_fractions(collapsed).set_index("Symbol")
    assert ranked.loc["RBMY1F/J", "frac_samples_top5pct"] == expected_fraction


@pytest.mark.parametrize("scope", ["cta", "genome"])
@pytest.mark.parametrize(
    "module_name", ["generate_cohort_percentiles", "generate_within_sample_top5"]
)
def test_generators_persist_registry_identity(tmp_path, rbmy_expression, scope, module_name):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    rbmy_expression.to_parquet(input_dir / "PRAD.parquet", index=False)
    out_dir = tmp_path / "output"
    importlib.import_module(module_name).build(
        input_dir, drop_genes=set(), out_dir=out_dir, proteoform=True, scope=scope
    )
    out = pd.read_parquet(out_dir / "PRAD.parquet")
    assert out.attrs == _proteoform_registry_metadata(scope)
    assert "RBMY1F/J" in set(out["proteoform_key"])


@pytest.fixture
def old_registry_shards(tmp_path, monkeypatch, rbmy_expression):
    # Simulate the omitted RBMY pair in the previous registry. Old per-gene
    # percentiles/fractions must never be combined to repair these artifacts.
    old_map = {k: v for k, v in proteoform_group_map().items() if k != "RBMY1F/RBMY1J"}
    old = sum_proteoform_tpm(
        rbmy_expression, old_map, group_symbols={k: proteoform_symbol(k) for k in old_map}
    )
    for dataset_name, build in [
        ("percentiles", cohort_percentile_vectors),
        ("within_sample", within_sample_top_fractions),
    ]:
        directory = tmp_path / expression.SHARD_DATASETS[dataset_name].subdir(proteoform=True)
        directory.mkdir()
        build(old).to_parquet(directory / "PRAD.parquet", index=False)
    monkeypatch.setenv("CANCERDATA_BUNDLED_DATA", str(tmp_path))
    monkeypatch.setattr(
        expression, "per_sample_expression", lambda code, **kw: rbmy_expression.copy()
    )
    return tmp_path


@pytest.mark.parametrize("metadata", [None, "old_digest", "wrong_scope"])
@pytest.mark.parametrize("dataset_name", ["percentiles", "within_sample"])
def test_stale_or_unverified_shards_recompute(old_registry_shards, metadata, dataset_name):
    path = old_registry_shards / expression.SHARD_DATASETS[dataset_name].subdir(proteoform=True)
    if metadata:
        stale = pd.read_parquet(path / "PRAD.parquet")
        stale.attrs.update(_proteoform_registry_metadata("cta"))
        if metadata == "old_digest":
            stale.attrs["proteoform_registry_sha256"] = "old"
        else:
            stale.attrs["proteoform_scope"] = "genome"
        stale.to_parquet(path / "PRAD.parquet", index=False)
    if dataset_name == "percentiles":
        out = expression.proteoform_cohort_percentiles("PRAD", auto_fetch=False)
        assert out.set_index("Symbol").loc["RBMY1F/J", "p95"] == pytest.approx(12.0, rel=0.002)
    else:
        out = expression.proteoform_within_sample_top_fraction("PRAD", auto_fetch=False)
        assert out.set_index("Symbol").loc["RBMY1F/J", "frac_samples_top5pct"] == 1.0
    assert out.attrs["source"] == "recomputed"
    assert not set(_RBMY_IDS).intersection(out["proteoform_key"])
    assert (
        out.attrs["proteoform_registry_sha256"]
        == _proteoform_registry_metadata("cta")["proteoform_registry_sha256"]
    )


@pytest.mark.parametrize("dataset_name", ["percentiles", "within_sample"])
def test_stale_shard_without_source_fails_clearly(old_registry_shards, monkeypatch, dataset_name):
    def no_matrix(*args, **kwargs):
        raise FileNotFoundError("no source matrix")

    monkeypatch.setattr(expression, "per_sample_expression", no_matrix)
    read = (
        expression.proteoform_cohort_percentiles
        if dataset_name == "percentiles"
        else expression.proteoform_within_sample_top_fraction
    )
    with pytest.raises(ValueError, match="stale or unverified proteoform registry"):
        read("PRAD", auto_fetch=False)


def test_recomputed_percentiles_do_not_claim_stale_artifact_qc(old_registry_shards, monkeypatch):
    def no_stale_qc(*args, **kwargs):
        pytest.fail("recomputed values must not use the stale artifact's QC metadata")

    monkeypatch.setattr(expression, "_require_expression_artifact_sample_qc", no_stale_qc)
    out = expression.proteoform_cohort_percentiles("PRAD", auto_fetch=False)
    assert out.attrs["artifact_sample_qc_verified"] is False


@pytest.mark.parametrize("dataset_name", ["percentiles", "within_sample"])
def test_local_availability_excludes_stale_shards_without_sources(
    old_registry_shards, monkeypatch, dataset_name
):
    from oncoref import source_matrices

    monkeypatch.setattr(source_matrices, "available_cohorts", lambda: ["PRAD"])
    monkeypatch.setattr(source_matrices, "is_cached", lambda code: False)
    available = (
        expression.locally_available_percentile_cohorts
        if dataset_name == "percentiles"
        else expression.locally_available_within_sample_cohorts
    )
    assert available(proteoform=True) == []
    monkeypatch.setattr(source_matrices, "is_cached", lambda code: True)
    assert available(proteoform=True) == ["PRAD"]
    assert available(proteoform=True, include_recomputable=False) == []


def test_registry_fingerprint_is_order_independent_and_tracks_membership(monkeypatch):
    from oncoref import proteoforms

    groups = proteoform_group_map()
    original = _proteoform_registry_metadata("cta")
    reordered = {label: tuple(reversed(ids)) for label, ids in reversed(list(groups.items()))}
    monkeypatch.setattr(proteoforms, "proteoform_group_map", lambda **kw: reordered)
    assert _proteoform_registry_metadata("cta") == original
    del reordered["RBMY1F/RBMY1J"]
    assert _proteoform_registry_metadata("cta") != original
