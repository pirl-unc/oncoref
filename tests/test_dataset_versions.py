"""Concrete release selection across the catalog, with no version substitution."""

import json
from dataclasses import replace

import pytest

from oncoref import catalog, cli, data_bundle, gtex, reference_data, source_matrices, species_assets


def test_latest_is_concrete_and_old_reference_paths_are_distinct(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    assert catalog.versions("gtex_tissue_rna") == ("v11", "v10", "v8")
    assert catalog.resolve_version("gtex_tissue_rna", "latest") == "v11"
    paths = {
        reference_data.local_path("gtex_tissue_rna", v) for v in catalog.versions("gtex_tissue_rna")
    }
    assert len(paths) == 3
    assert all("latest" not in p.parts for p in paths)
    assert catalog.resolve_version("hpa_rna_consensus") == "v25.1"
    assert catalog.resolve_version("hpa_rna_consensus", "v23") == "v23"
    assert (
        catalog.resolve_version("hpa_normal_tissue") == "v23"
    )  # newest available legacy IHC table
    assert all("latest" not in spec["urls"] for spec in reference_data.REFERENCE_SOURCES.values())
    with pytest.raises(reference_data.ReferenceDataError, match="available"):
        catalog.ensure("gtex_tissue_rna", version="v999")


def test_gtex_reader_records_annotation_and_version_without_rewriting_identity(
    monkeypatch, tmp_path
):
    fixture = tmp_path / "gtex.gct"
    fixture.write_text("#1.2\n1\t1\nName\tDescription\tPituitary\nENSG00000287363.4\tGUCA1ANB\t0\n")
    calls = []

    def ensure(name, version):
        calls.append((name, version))
        return fixture

    monkeypatch.setattr(reference_data, "ensure", ensure)
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path / "cache"))
    current = gtex.tissue_rna()
    old = gtex.tissue_rna("v8")
    assert calls == [("gtex_tissue_rna", "v11"), ("gtex_tissue_rna", "v8")]
    assert current.attrs["version"] == "v11" and old.attrs["version"] == "v8"
    assert current.attrs["reference"] == "GENCODE v47"
    assert current.attrs["unit"] == "TPM"
    assert current.iloc[0].Name == "ENSG00000287363.4"
    assert current.iloc[0].Description == "GUCA1ANB" and current.iloc[0].Pituitary == 0


def test_catalog_forwards_exact_release_to_each_backend(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(reference_data, "ensure", lambda n, v: calls.append((n, v)) or tmp_path)
    monkeypatch.setattr(
        source_matrices, "ensure", lambda n, **kw: calls.append((n, kw["version"])) or tmp_path
    )
    monkeypatch.setattr(
        data_bundle,
        "ensure_local",
        lambda **kw: calls.append(("bundle", kw["version"])) or tmp_path,
    )
    monkeypatch.setattr(data_bundle, "find", lambda n, **kw: tmp_path / kw["version"] / n)
    catalog.ensure("gtex_tissue_rna", "v10")
    catalog.ensure("per-sample:LUAD", "5.22.13")
    result = catalog.ensure("pan-cancer-expression.csv", "5.23.25")
    assert calls == [("gtex_tissue_rna", "v10"), ("LUAD", "5.22.13"), ("bundle", "5.23.25")]
    assert result == tmp_path / "5.23.25" / "pan-cancer-expression.csv"
    with pytest.raises(ValueError, match="one dataset"):
        catalog.fetch("all", version="v8")


def test_historical_matrix_never_reuses_current_cache(monkeypatch, tmp_path):
    monkeypatch.setenv(source_matrices.CACHE_DIR_ENV_VAR, str(tmp_path))
    current = source_matrices.local_path("LUAD")
    current.write_bytes(b"current")
    old = source_matrices.local_path("LUAD", "5.22.13")
    assert old != current and not old.exists()
    assert catalog.path("per-sample:LUAD", "5.22.13") is None
    assert "/source-v5.22.13/" in source_matrices.release_url("LUAD", "5.22.13")
    with pytest.raises(source_matrices.SourceMatrixError):
        source_matrices.local_path("LUAD", "../../other")


def test_historical_bundle_validates_requested_release(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_BUNDLED_DATA", str(tmp_path / "current"))
    old = "5.23.25"
    seen = {}

    def manifest(source, **kwargs):
        seen.update(kwargs)
        assert "/v5.23.25/" in source["url"]
        return {"tarball": {"sha256": "a" * 64}}

    def extract(url, root, **kwargs):
        seen["written_version"] = kwargs["data_version"]
        assert root != data_bundle.cache_dir()

    monkeypatch.setattr(data_bundle, "_fetch_release_manifest", manifest)
    monkeypatch.setattr(data_bundle, "_download_and_extract", extract)
    data_bundle.fetch(version=old, verbose=False)
    assert seen == {
        "expected_data_version": old,
        "expected_source_matrix_version": None,
        "written_version": old,
    }


def test_species_versions_select_exact_asset_and_content_path(monkeypatch, tmp_path):
    first = species_assets.species_assets("dog")[0]
    second = replace(first, source_version="next-release", sha256="a" * 64)
    monkeypatch.setattr(
        species_assets, "_manifest", lambda species: {"assets": [first.as_dict(), second.as_dict()]}
    )
    monkeypatch.setenv("ONCOREF_ASSET_CACHE", str(tmp_path))
    assert catalog.versions(first.id, species="dog") == ("next-release", first.source_version)
    assert catalog.resolve_version(first.id, species="dog") == "next-release"
    assert (
        species_assets.species_asset(first.id, species="dog", version=first.source_version) == first
    )
    assert species_assets.asset_path(first) != species_assets.asset_path(second)
    with pytest.raises(ValueError, match="available"):
        catalog.resolve_version(first.id, "unknown", species="dog")


def test_packaged_tables_expose_version_and_content_identity():
    record = catalog.provenance("cta-selection-human")
    assert record["version"].startswith("oncoref-")
    assert len(record["sha256"]) == 64 and record["verification_state"] == "content_recorded"
    assert catalog.path("cta-selection-human", version="oncoref-0.0.1") is None


def test_cli_versions_and_provenance(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    assert cli.main(["data", "versions", "gtex_tissue_rna"]) == 0
    assert json.loads(capsys.readouterr().out) == ["v11", "v10", "v8"]
    assert cli.main(["data", "provenance", "gtex_tissue_rna", "--version", "v10"]) == 0
    assert json.loads(capsys.readouterr().out)["version"] == "v10"


@pytest.mark.parametrize(
    "name,species,relative,content",
    [
        ("cta-selection-human", "human", "cta-selection-human.csv.gz", b"historical table"),
        (
            "cta-evidence",
            "dog",
            "species/dog/cta-evidence.json",
            b'{"source_version":"2026-10-01"}',
        ),
    ],
)
def test_historical_curated_data_verify_wheel_without_executing_code(
    monkeypatch, tmp_path, name, species, relative, content
):
    import hashlib
    import io
    import zipfile

    from oncoref import packaged_data

    monkeypatch.setenv("CANCERDATA_DATA_DIR", str(tmp_path))
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr("oncoref/data/" + relative, content)
        archive.writestr("oncoref/__init__.py", "raise RuntimeError('must never execute')")
    wheel_bytes = payload.getvalue()
    metadata = {
        "urls": [
            {
                "filename": "oncoref-1.20.0-py3-none-any.whl",
                "url": "https://files.pythonhosted.org/test.whl",
                "digests": {"sha256": hashlib.sha256(wheel_bytes).hexdigest()},
            }
        ]
    }
    monkeypatch.setattr(
        packaged_data.urllib.request,
        "urlopen",
        lambda *a, **kw: io.BytesIO(json.dumps(metadata).encode()),
    )
    downloads = []

    def fetch(url, destination, **kwargs):
        downloads.append(url)
        destination.write_bytes(wheel_bytes)

    monkeypatch.setattr(packaged_data._downloads, "fetch_file", fetch)
    assert catalog.resolve_version(name, species=species).startswith("oncoref-")
    dest = catalog.ensure(name, species=species, version="oncoref-1.20.0")
    assert dest.read_bytes() == content
    assert "oncoref-1.20.0" in catalog.versions(name, species=species)
    assert (
        catalog.provenance(name, "oncoref-1.20.0", species=species, verify_content=True)[
            "verification_state"
        ]
        == "verified"
    )
    if species == "dog":
        assert (
            catalog.provenance(name, "oncoref-1.20.0", species=species)["source_version"]
            == "2026-10-01"
        )
    catalog.ensure(name, species=species, version="oncoref-1.20.0")
    assert len(downloads) == 1
    dest.write_bytes(b"corrupted")
    assert (
        catalog.provenance(name, "oncoref-1.20.0", species=species, verify_content=True)[
            "verification_state"
        ]
        == "checksum_mismatch"
    )
    catalog.ensure(name, species=species, version="oncoref-1.20.0")
    assert len(downloads) == 2 and dest.read_bytes() == content
    with pytest.raises(FileNotFoundError, match="absent"):
        packaged_data.ensure("not-in-wheel.csv", "oncoref-1.20.0")
    assert not (dest.parent / "__init__.py").exists()


def test_bundle_override_named_like_older_release_does_not_collapse_versions(monkeypatch, tmp_path):
    monkeypatch.setenv("CANCERDATA_BUNDLED_DATA", str(tmp_path / "v5.23.25"))
    assert data_bundle.cache_dir() != data_bundle.cache_dir("5.23.25")
