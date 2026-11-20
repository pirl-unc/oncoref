"""Build ordering and receipt integrity are independent of available raw datasets."""

import hashlib
import json

import pytest

from oncoref import curation, reference_data


def test_dependency_order_and_invalid_graphs():
    steps = {
        "report": {"needs": ["catalog"]},
        "catalog": {"needs": ["rare"]},
        "rare": {"needs": []},
    }
    assert curation.ordered_steps(steps) == ["rare", "catalog", "report"]
    with pytest.raises(ValueError, match="Missing or cyclic"):
        curation.ordered_steps({"a": {"needs": ["b"]}})
    with pytest.raises(ValueError, match="Missing or cyclic"):
        curation.ordered_steps({"a": {"needs": ["b"]}, "b": {"needs": ["a"]}})


def test_receipt_detects_changed_implementation_and_missing_output(tmp_path, monkeypatch):
    monkeypatch.setattr(curation, "manifest", lambda: {"steps": {"catalog": {"needs": []}}})
    receipt = {
        "policy_id": curation.target_selection.POLICY_ID,
        "steps": ["catalog"],
        "implementation": {"builder.py": hashlib.sha256(b"code").hexdigest()},
        "outputs": {"catalog.csv": hashlib.sha256(b"data").hexdigest()},
    }
    path = tmp_path / curation.RECEIPT
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(receipt))
    (tmp_path / "builder.py").write_bytes(b"code")
    (tmp_path / "catalog.csv").write_bytes(b"data")
    assert curation.check(tmp_path)["verified_outputs"] == 1
    (tmp_path / "builder.py").write_bytes(b"new code")
    with pytest.raises(ValueError, match=r"builder\.py"):
        curation.check(tmp_path)
    (tmp_path / "catalog.csv").unlink()
    with pytest.raises(ValueError, match=r"catalog\.csv"):
        curation.check(tmp_path)


def test_content_addressed_import_validates_and_repairs_bytes(tmp_path, monkeypatch):
    monkeypatch.setenv(reference_data.CACHE_DIR_ENV_KEY, str(tmp_path / "cache"))
    data = b"pinned reference"
    digest = hashlib.sha256(data).hexdigest()
    spec = {
        "filename": "input.tsv",
        "urls": {"v1": "https://example.invalid/input"},
        "pins": {"v1": {"sha256": digest, "bytes": len(data)}},
    }
    monkeypatch.setattr(reference_data, "REFERENCE_SOURCES", {"test": spec})
    source = tmp_path / "source.tsv"
    source.write_bytes(data)
    target = reference_data.pinned_path("test", "v1", import_path=source)
    assert target.parent.name == digest and target.read_bytes() == data
    target.write_bytes(b"damaged")
    assert reference_data.pinned_path("test", "v1", import_path=source).read_bytes() == data
    target.unlink()
    source.write_bytes(b"changed reference")
    with pytest.raises(reference_data.ReferenceDataError, match="does not match"):
        reference_data.pinned_path("test", "v1", import_path=source)


def test_required_curation_inputs_are_all_registered_and_portable():
    from pathlib import Path

    from oncoref import catalog

    refs = curation.required_reference_inputs()
    ref_sources = {source for source, _ in refs}
    assert {
        "hpa_individual_tissue_tpm",
        "hpa_rna_consensus",
        "hpa_single_cell_tissue",
        "hpa_normal_tissue",
        "rna_tissue_hpa",
        "rna_tissue_gtex",
        "rna_brain_gtex",
        "rna_brain_hpa",
        "gtex_tissue_rna",
        "gtex_individual_tissue_tpm",
        "hpa_immune_cell_samples",
        "hpa_immune_monaco",
        "petropoulos_embryonic_tpm",
        "tyser_embryonic_tpm",
    } <= ref_sources
    for source, version in refs:
        assert source in reference_data.REFERENCE_SOURCES
        assert version in reference_data.REFERENCE_SOURCES[source]["urls"]

    root = Path(__file__).resolve().parents[1]
    tumors = curation.required_tumor_matrices(root)
    assert len(tumors) == 132
    for code, version, sha256 in tumors:
        assert "/" not in code and not version.startswith("v") and len(sha256) == 64
        assert catalog.resolve_version(f"per-sample:{code}", version) == version
