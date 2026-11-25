"""Verified policy-bound sequence consumers retain occurrences and explicit gaps."""

import gzip
import json
import zipfile
from pathlib import Path

import pytest

from oncoref import data_file_renames, target_selection
from oncoref import sequence_bundle as bundle
from tests import released_inputs


def test_packaged_asset_pin_matches_current_policy_and_complete_inventory():
    root = Path(__file__).resolve().parents[1]
    audit = root / "docs/audits/sequence-bundle-20261121"
    manifest = json.loads((audit / "manifest.json").read_text())
    pin = bundle.descriptor()
    assert pin["bundle_id"] == manifest["bundle_id"] == bundle.content_id(manifest)
    assert manifest["selection_policy_id"] == target_selection.BASELINE_2026_11_20_POLICY_ID
    assert manifest["role_policy_id"] == bundle.ROLE_POLICY_ID
    assert manifest["record_count"] == sum(manifest["record_roles"].values())
    assert manifest["record_count"] == 245535
    assert manifest["occurrence_count"] >= manifest["record_count"]
    for name, expected in manifest["package_dependencies"].items():
        try:
            raw = released_inputs.read(f"oncoref/data/{name}", expected)
        except AssertionError:
            raw = data_file_renames.as_released(name)
        assert bundle.digest(raw) == expected

    historical_audit = root / "docs/audits/sequence-bundle-20261104"
    historical_manifest = json.loads((historical_audit / "manifest.json").read_text())
    assert historical_manifest["bundle_id"] == bundle.content_id(historical_manifest)
    transition = json.loads((historical_audit / "historical-transition.json").read_text())
    receipt = json.loads(
        (root / "docs/audits/cta-definition-20261104/release-receipt.json").read_text()
    )
    legacy_pin = json.loads(
        released_inputs.read(
            "oncoref/data/sequence-bundle.json",
            receipt["inputs"]["oncoref/data/sequence-bundle.json"],
        )
    )
    assert transition["current_bundle_id"] == legacy_pin["bundle_id"]
    assert transition["historical_policy_id"]


def test_target_selection_dispositions_cover_all_candidates_and_exclude_unsafe():
    dispositions = bundle.target_selection_dispositions()
    assert len(dispositions) == len(target_selection.results())

    unsafe = dispositions[dispositions.category.eq("Unsafe")]
    assert (unsafe.shared_sequence == "exclude").all()
    assert (~unsafe.direct_target_candidate).all()
    assert (~unsafe.cancer_antigen_candidate).all()

    ctas = dispositions[dispositions.category.eq("CTA")]
    assert len(ctas) == 323
    assert (ctas.shared_sequence == "allow").all()
    assert ctas.direct_target_candidate.all()
    assert ctas.cancer_antigen_candidate.all()
    assert ctas.reproductive_antigen_candidate.all()
    assert ctas.representative_sequence_sha256s.map(bool).all()

    ca_only = dispositions[dispositions.category.eq("Cancer Antigen only")]
    assert len(ca_only) == 97
    assert (ca_only.shared_sequence == "allow").all()
    assert (~ca_only.direct_target_candidate).all()
    assert ca_only.cancer_antigen_candidate.all()
    assert ca_only.representative_sequence_sha256s.map(bool).all()

    ra_only = dispositions[dispositions.category.eq("Reproductive Antigen only")]
    assert len(ra_only) == 142
    assert (ra_only.shared_sequence == "allow").all()
    assert (~ra_only.direct_target_candidate).all()
    assert ra_only.reproductive_antigen_candidate.all()

    for subtype in ("Other:Review", "Other:Sequence", "Other:Unscreened tissue"):
        sub = dispositions[dispositions.classification.eq(subtype)]
        assert not sub.empty
        assert (sub.shared_sequence == "review").all()

    by_label = dispositions.set_index("protein_group_label")
    for rescued in (
        "MAGEC2",
        "CSAG1",
        "PLAC1",
        "DKKL1",
        "CTAG2",
        "DUX4",
        "MAGEB6",
        "FOXR2",
        "CSN2",
    ):
        assert by_label.loc[rescued, "shared_sequence"] == "allow"
        assert bool(by_label.loc[rescued, "direct_target_candidate"])
    for newly_unsafe in ("OBP2A", "SLCO6A1", "SPATA45", "COX7B2", "SP7", "GFY", "DAZL"):
        assert by_label.loc[newly_unsafe, "shared_sequence"] == "exclude"
    for unscreened in ("DSPP", "CALHM3"):
        assert by_label.loc[unscreened, "shared_sequence"] == "review"


def test_cached_current_bundle_enforces_target_selection_roles():
    try:
        path = bundle.ensure(download=False)
    except FileNotFoundError:
        pytest.skip("Current sequence bundle is not staged in local cache")
    loaded = bundle.load(path)
    ts = target_selection.results()
    unsafe_genes = {
        g
        for members in ts.loc[ts.category.eq("Unsafe"), "member_gene_ids"]
        for g in json.loads(members)
    }
    cta_entities = set(ts.loc[ts.category.eq("CTA"), "entity_id"])
    ca_entities = set(ts.loc[ts.category.eq("Cancer Antigen only"), "entity_id"])
    covered_cta = set()
    covered_ca = set()
    for record in loaded.records():
        if record["shared_sequence"] == "allow":
            for occ in record["occurrences"]:
                assert occ["gene_id"] not in unsafe_genes
                if occ["direct_target_candidate"]:
                    covered_cta.update(occ["decision_ids"])
                if occ.get("cancer_antigen_candidate"):
                    covered_ca.update(occ["decision_ids"])
    assert cta_entities <= covered_cta
    assert ca_entities <= covered_ca
    seq_tbl = loaded.to_sequences_table()
    label_to_gid = {
        row.protein_group_label: json.loads(row.member_gene_ids)[0]
        for row in ts.itertuples(index=False)
    }
    roles_by_gid = seq_tbl.groupby("gene_id")["shared_sequence"].unique().to_dict()
    for unsafe_lbl in ("SP7", "GFY", "DAZL"):
        assert set(roles_by_gid[label_to_gid[unsafe_lbl]]) == {"exclude"}
    for allow_lbl in ("PLAC1", "ZAN", "DKKL1", "MAGEB6", "CTAG2", "DUX4"):
        assert "allow" in set(roles_by_gid[label_to_gid[allow_lbl]])


def make_bundle(folder, *, role="allow", second=False):
    folder.mkdir(parents=True, exist_ok=True)
    sequence = "MARTXACD*"
    segments = [
        {
            "start": start,
            "end": end,
            "sequence": sequence[start:end],
            "sha256": bundle.digest(sequence[start:end].encode()),
        }
        for start, end in [(0, 4), (5, 8)]
    ]
    records = []
    for i, r in enumerate([role, "exclude"] if second else [role]):
        records.append(
            {
                "record_id": str(i),
                "reference_id": "reference",
                "sequence": sequence,
                "sequence_sha256": bundle.digest(sequence.encode()),
                "valid_segments": segments,
                "shared_sequence": r,
                "occurrences": [{"occurrence_id": str(i), "gene_id": str(i), "shared_sequence": r}],
                "source_record": {"original_accession": str(i)},
            }
        )
    files = {
        "records.jsonl.gz": gzip.compress(
            ("\n".join(json.dumps(r) for r in records) + "\n").encode(), mtime=0
        ),
        "coverage.json": json.dumps(
            {"reference_id": "reference", "missing_translations": ["not_measured"]}
        ).encode(),
        "decisions.csv.gz": gzip.compress(b"gene,role\n", mtime=0),
        "orfdata-proteome.json.gz": gzip.compress(b"{}", mtime=0),
    }
    manifest = {
        "schema": bundle.SCHEMA,
        "schema_version": 1,
        "reference": {"reference_id": "reference"},
        "record_count": len(records),
        "occurrence_count": len(records),
        "package_dependencies": {},
        "payloads": {n: {"bytes": len(v), "sha256": bundle.digest(v)} for n, v in files.items()},
    }
    manifest["bundle_id"] = bundle.content_id(manifest)
    for name, value in files.items():
        (folder / name).write_bytes(value)
    (folder / "manifest.json").write_text(json.dumps(manifest))
    return manifest


def test_originals_gaps_segments_and_duplicate_occurrences_survive(tmp_path):
    m = make_bundle(tmp_path / "bundle", second=True)
    b = bundle.load(tmp_path / "bundle", expected_id=m["bundle_id"], require_current=False)
    assert b.identity == m["bundle_id"]
    assert b.coverage["missing_translations"] == ["not_measured"]
    groups = list(b.groups())
    assert (
        len(groups) == 1 and groups[0]["record_count"] == 2 and groups[0]["occurrence_count"] == 2
    )
    assert groups[0]["shared_sequence"] == "exclude"
    segments = list(b.segments(role="exclude", min_length=4))
    assert len(segments) == 1 and segments[0]["sequence"] == "MART"
    assert next(b.records())["sequence"] == "MARTXACD*"
    first = b.manifest
    first["record_count"] = 99
    assert b.manifest["record_count"] == 2


def test_missing_trusted_pin_mixed_and_tampered_payloads_fail(tmp_path):
    m = make_bundle(tmp_path / "a")
    with pytest.raises(ValueError, match="trusted"):
        bundle.load(tmp_path / "a", require_current=False)
    n = make_bundle(tmp_path / "b", role="review")
    with pytest.raises(ValueError, match="identity"):
        bundle.load(tmp_path / "b", expected_id=m["bundle_id"], require_current=False)
    (tmp_path / "b/records.jsonl.gz").write_bytes((tmp_path / "a/records.jsonl.gz").read_bytes())
    with pytest.raises(ValueError, match="checksum"):
        bundle.load(tmp_path / "b", expected_id=n["bundle_id"], require_current=False)


def test_cached_default_requires_current_package_dependencies(tmp_path, monkeypatch):
    m = make_bundle(tmp_path / "a")
    monkeypatch.setattr(bundle, "descriptor", lambda: {"bundle_id": m["bundle_id"]})
    assert bundle.load(tmp_path / "a").identity == m["bundle_id"]
    with pytest.raises(ValueError, match="current package"):
        bundle.load(tmp_path / "a", expected_id="other")


@pytest.mark.parametrize("mutation", ["missing_segment", "unknown_role"])
def test_even_pinned_records_must_preserve_segments_and_explicit_roles(tmp_path, mutation):
    folder = tmp_path / "bundle"
    manifest = make_bundle(folder)
    path = folder / "records.jsonl.gz"
    record = json.loads(gzip.decompress(path.read_bytes()))
    if mutation == "missing_segment":
        record["valid_segments"].pop()
    else:
        record["occurrences"][0]["shared_sequence"] = "unknown"
    payload = gzip.compress((json.dumps(record) + "\n").encode(), mtime=0)
    path.write_bytes(payload)
    manifest["payloads"][path.name] = {"bytes": len(payload), "sha256": bundle.digest(payload)}
    manifest["bundle_id"] = bundle.content_id(manifest)
    (folder / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match=r"AA20 run|explicit sequence role"):
        bundle.load(folder, expected_id=manifest["bundle_id"], require_current=False)


def test_opt_in_bounded_acquisition_and_corrupt_cache_recovery(tmp_path, monkeypatch):
    source = tmp_path / "source"
    m = make_bundle(source)
    archive = tmp_path / "source.zip"
    with zipfile.ZipFile(archive, "w") as z:
        for p in source.iterdir():
            z.write(p, p.name)
    data = archive.read_bytes()
    pin = {
        "bundle_id": m["bundle_id"],
        "url": "https://example.org/pinned.zip",
        "archive_bytes": len(data),
        "archive_sha256": bundle.digest(data),
        "unpacked_bytes": sum(p.stat().st_size for p in source.iterdir()),
    }
    monkeypatch.setattr(bundle, "descriptor", lambda: dict(pin))
    cache = tmp_path / "cache"
    with pytest.raises(FileNotFoundError):
        bundle.ensure(cache_dir=cache)
    assert not cache.exists()

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def raise_for_status(self):
            pass

        def iter_content(self, size):
            yield data

    import requests

    monkeypatch.setattr(requests, "get", lambda *a, **k: Response())
    destination = bundle.ensure(cache_dir=cache, download=True)
    assert bundle.load(destination).identity == m["bundle_id"]
    (destination / "coverage.json").write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        bundle.ensure(cache_dir=cache)
    assert bundle.ensure(cache_dir=cache, download=True) == destination
    assert bundle.load(destination).coverage["missing_translations"]
    pin["archive_bytes"] -= 1
    with pytest.raises(ValueError, match="bound"):
        bundle.ensure(cache_dir=tmp_path / "small", download=True)
