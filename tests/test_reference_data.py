# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import pytest

from cancerdata import reference_data


def test_cache_dir_honors_env(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    assert reference_data.cache_dir() == tmp_path / "sources"
    assert reference_data.cache_dir().exists()


def test_resolve_version_default_and_errors():
    assert reference_data.resolve_version("hpa_rna_consensus") == "v23"
    assert reference_data.resolve_version("hpa_rna_consensus", "latest") == "latest"
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.resolve_version("hpa_normal_tissue", "v99")
    with pytest.raises(reference_data.ReferenceDataError):
        reference_data.resolve_version("not_a_source")


def test_local_path_and_is_cached(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    p = reference_data.local_path("hpa_rna_consensus")
    assert p == tmp_path / "sources" / "hpa_rna_consensus" / "v23" / "rna_tissue_consensus.tsv"
    assert not reference_data.is_cached("hpa_rna_consensus")
    p.parent.mkdir(parents=True)
    p.write_text("Gene\tTissue\tnTPM\n")
    assert reference_data.is_cached("hpa_rna_consensus")


def test_status_shape(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    rows = reference_data.status()
    names = {r["name"] for r in rows}
    assert names == {"hpa_rna_consensus", "hpa_normal_tissue", "hpa_single_cell"}
    assert all(r["cached"] is False for r in rows)  # nothing downloaded


def test_single_cell_registered():
    # The single-cell nTPM source is part of the registry.
    spec = reference_data.REFERENCE_SOURCES["hpa_single_cell"]
    assert spec["filename"] == "rna_single_cell_type.tsv"
    assert "v23" in spec["urls"]
