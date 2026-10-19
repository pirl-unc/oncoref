"""Species listings stay readable, lossless and free of downloads."""

import json

import pytest

from oncoref import catalog, cli, data_bundle, expression_registry, reference_data, source_matrices
from oncoref import species_assets as assets


@pytest.fixture(autouse=True)
def no_downloads(monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("A listing tried to download data")

    for module, method in (
        (data_bundle, "fetch"),
        (reference_data, "download"),
        (source_matrices, "fetch"),
        (assets, "ensure_asset"),
    ):
        monkeypatch.setattr(module, method, unexpected)
    monkeypatch.setenv("COLUMNS", "100")


@pytest.mark.parametrize("command", [["data", "list"], ["expression-sources"]])
@pytest.mark.parametrize("json_output", [False, True])
def test_omitted_species_is_exactly_explicit_human(command, json_output, capsys):
    options = ["--json"] if json_output else []
    assert cli.main([*command, *options]) == 0
    default = capsys.readouterr()
    assert not default.err
    assert cli.main([*command, "--species", "human", *options]) == 0
    assert capsys.readouterr() == default


@pytest.mark.parametrize("command", [["data", "list"], ["expression-sources"]])
@pytest.mark.parametrize("alias", ["dog", "9615", "Canis lupus familiaris"])
def test_canine_aliases_are_equivalent(command, alias, capsys):
    assert cli.main([*command, "--species", "canine"]) == 0
    expected = capsys.readouterr()
    assert "Canine" in expected.out
    assert "NCBI taxon 9615" in expected.out
    assert cli.main([*command, "--species", alias]) == 0
    assert capsys.readouterr() == expected


@pytest.mark.parametrize("width", [40, 80, 100])
@pytest.mark.parametrize("species", ["human", "canine"])
@pytest.mark.parametrize("command", [["data", "list"], ["expression-sources"]])
def test_listing_fits_terminal_and_keeps_every_id(command, species, width, capsys, monkeypatch):
    monkeypatch.setenv("COLUMNS", str(width))
    assert cli.main([*command, "--species", species]) == 0
    result = capsys.readouterr()
    assert not result.err
    assert all(len(line) <= width for line in result.out.splitlines())
    if command[0] == "data":
        ids = [r["name"] for r in catalog.inventory(species=species)]
    else:
        sources = expression_registry.expression_sources(species=species)
        ids = [s.id for s in sources] + [code for s in sources for code in s.cancer_codes]
    compact = "".join(result.out.split())
    for identifier in ids:
        assert identifier in compact
    assert "..." not in result.out


def test_canine_source_output_explains_pending_reads(capsys):
    assert cli.main(["expression-sources", "--species", "canine"]) == 0
    output = capsys.readouterr().out
    assert "Normal / bgee" in output and "Tumor / pending-native-source" in output
    assert "metadata only" in output and "blocked provenance" in output
    assert "Read unavailable:" in output
    assert "unknown" in output
    assert "Cancer codes" not in output


@pytest.mark.parametrize("species", ["human", "canine"])
def test_inventory_json_preserves_catalog_fields(species, capsys):
    assert cli.main(["data", "list", "--species", species, "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == catalog.inventory(species=species)


def test_source_json_preserves_complete_codes_and_contracts(capsys):
    assert cli.main(["expression-sources", "--json"]) == 0
    rows = {r["id"]: r for r in json.loads(capsys.readouterr().out)}
    for source in expression_registry.expression_sources():
        assert rows[source.id]["cancer_codes"] == list(source.cancer_codes)
        assert rows[source.id]["contract"] == source.expression_contract.as_dict()


@pytest.mark.parametrize("command", [["data", "list"], ["expression-sources"]])
def test_empty_species_and_invalid_species_have_distinct_results(command, capsys):
    assert cli.main([*command, "--species", "mouse"]) == 0
    output = capsys.readouterr()
    assert not output.err
    assert "Mouse" in output.out and "No " in output.out
    assert cli.main([*command, "--species", "mouse", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == []
    assert cli.main([*command, "--species", "wolf"]) == 1
    error = capsys.readouterr()
    assert not error.out
    assert "Unknown species" in error.err and "Traceback" not in error.err


def test_filters_and_empty_selection(capsys):
    assert cli.main(["data", "list", "--species", "canine", "--group", "references", "--json"]) == 0
    assert {r["group"] for r in json.loads(capsys.readouterr().out)} == {"references"}
    assert (
        cli.main(["expression-sources", "--species", "canine", "--context", "tumor", "--json"]) == 0
    )
    assert [r["id"] for r in json.loads(capsys.readouterr().out)] == ["gse304066-osteosarcoma"]
    assert cli.main(["expression-sources", "--species", "human", "--context", "normal"]) == 0
    assert "No matching" in capsys.readouterr().out
    assert cli.main(["expression-sources", "--type", "gdc", "--code", "ACC", "--json"]) == 0
    assert [r["id"] for r in json.loads(capsys.readouterr().out)] == ["tcga-acc"]


def test_asset_listing_uses_explicit_overrides_without_claiming_integrity(
    tmp_path, capsys, monkeypatch
):
    # A listing reports presence only. Arbitrary bytes must not be called verified.
    asset = assets.species_asset("ensembl98-gtf", species="canine")
    local = tmp_path / "local.gtf.gz"
    local.write_bytes(b"not verified")
    overrides = tmp_path / "overrides.json"
    overrides.write_text(json.dumps({f"9615:{asset.id}": str(local)}))
    monkeypatch.setenv("ONCOREF_ASSET_OVERRIDES", str(tmp_path / "absent-default.json"))
    assert (
        cli.main(
            [
                "data",
                "list",
                asset.id,
                "--species",
                "canine",
                "--asset-overrides",
                str(overrides),
                "--json",
            ]
        )
        == 0
    )
    (row,) = json.loads(capsys.readouterr().out)
    assert row["available"] is True
    assert row["integrity"] == "not_checked"


def test_source_matrix_local_count_is_explicit(monkeypatch, capsys):
    monkeypatch.setattr(source_matrices, "available_cohorts", lambda: ["ACC", "LUAD"])
    monkeypatch.setattr(source_matrices, "is_cached", lambda code: code == "ACC")
    assert cli.main(["data", "list", "source"]) == 0
    output = capsys.readouterr().out
    assert "1/2" in output
    assert "1 of 2 cohort matrices present locally" in output
