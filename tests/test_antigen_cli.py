"""End-user commands must use the current catalog and fail closed on background gaps."""

import json

import pandas as pd

from oncoref import antigen_sets, cli, target_selection


def test_current_list_matches_owner_and_backgrounds(capsys):
    assert cli.main(["antigens", "list", "--set", "CTA", "--ids", "--format", "json"]) == 0
    assert {r["gene_id"] for r in json.loads(capsys.readouterr().out)} == antigen_sets.gene_ids()
    assert cli.main(["antigens", "list", "--set", "Other", "--format", "json"]) == 0
    assert len(json.loads(capsys.readouterr().out)) == len(antigen_sets.proteins("Other"))
    assert cli.main(["antigens", "list", "--set", "Unsafe", "--mode", "strict"]) == 1
    assert "Candidate filters" in capsys.readouterr().err


def test_combined_cli_fold_and_patient_constraints(capsys):
    assert (
        cli.main(
            [
                "antigens",
                "list",
                "--set",
                "CA",
                "--min-tumor-fold",
                "30",
                "--min-patients",
                "2",
                "--format",
                "json",
            ]
        )
        == 0
    )
    result = json.loads(capsys.readouterr().out)
    rows = target_selection.results().set_index("entity_id")
    assert result
    assert all(
        json.loads(rows.loc[r["entity_id"], "independent_patient_tumor_folds"])[1] >= 30
        for r in result
    )


def test_explain_held_gene_and_member_symbol(capsys):
    assert cli.main(["antigens", "explain", "CAGE1", "--format", "json"]) == 0
    row = json.loads(capsys.readouterr().out)[0]
    assert row["classification"] == "Other:Review"
    assert row["rare_cell_status"] == "unresolved"
    assert row["strict"].startswith("fail:")
    assert cli.main(["antigens", "explain", "CTAG2", "--format", "json"]) == 0
    ctag2_row = json.loads(capsys.readouterr().out)[0]
    assert ctag2_row["classification"] == "CTA"
    assert cli.main(["antigens", "explain", "SSX2", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["protein_group_label"] == "SSX2/B"
    assert cli.main(["antigens", "explain", "not-a-gene"]) == 1


def test_explain_distinguishes_all_papers_from_tumor_support(capsys):
    assert cli.main(["antigens", "explain", "DPRX", "--format", "json"]) == 0
    row = json.loads(capsys.readouterr().out)[0]
    assert row["paper_source_count"] == 2
    assert row["papers_reproductive"] == 2
    assert row["papers_counted"] == 0


def test_explain_papers_retains_original_gene_name_and_tissue(capsys):
    assert cli.main(["antigens", "explain", "FOXL3", "--papers", "--format", "json"]) == 0
    row = json.loads(capsys.readouterr().out)[0]
    claim = row["paper_claims"][0]
    assert claim["reported_name"] == "WI2-2373I1.1"
    assert claim["provenance"]["source_tissue"] == "Testis"
    assert claim["allowed_uses"] == ["reproductive_support"]


def test_subtraction_exact_equivalence_fasta_and_missing_coverage(tmp_path, capsys, monkeypatch):
    background = pd.DataFrame({"gene_id": ["background"], "sequence": ["LLLLLLLL"]})
    monkeypatch.setattr(antigen_sets, "sequences", lambda _: background)
    query = tmp_path / "query.fa"
    query.write_text(">first\nIIIIIIII\n>second\nAAAAAAAA\n")
    args = [
        "antigens",
        "subtract",
        str(query),
        "--min-length",
        "8",
        "--max-length",
        "8",
        "--format",
        "json",
    ]
    assert cli.main(args) == 0
    assert len(json.loads(capsys.readouterr().out)) == 2
    assert cli.main([*args, "--equiv", "ILV"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert [r["query_id"] for r in rows] == ["second"]
    background.attrs["missing_gene_ids"] = ["missing"]
    assert cli.main(args) == 1
    assert "missing gene translations" in capsys.readouterr().err


def test_funnel_and_legacy_notice(capsys):
    assert cli.main(["antigens", "funnel", "--mode", "strict", "--format", "json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    sequential = [r["remaining"] for r in rows if r["kind"] == "sequential filter"]
    assert sequential == sorted(sequential, reverse=True)
    catalog_steps = [r["remaining"] for r in rows if r["kind"] == "sequential catalog"]
    assert catalog_steps[-1] == len(antigen_sets.proteins("CA|RA"))
    assert catalog_steps == sorted(catalog_steps, reverse=True)
    assert cli.main(["cta", "--count"]) == 0
    assert "Legacy CTA snapshot" in capsys.readouterr().err


def test_cli_host_profiles_and_tissue_toggles(capsys):
    args = [
        "antigens",
        "list",
        "--set",
        "CTA",
        "--profile",
        "female-patient",
        "--lactating",
        "--format",
        "json",
    ]
    assert cli.main(args) == 0
    capture = capsys.readouterr()
    rows = json.loads(capture.out)
    metadata = json.loads(capture.err)
    assert {r["entity_id"] for r in rows} == set(
        antigen_sets.proteins("CTA", profile="female-patient", lactating=True).entity_id
    )
    assert metadata["host_context"]["protected"] == ["lactating_breast"]
    assert (
        cli.main(
            [
                "antigens",
                "explain",
                "PRAME",
                "--profile",
                "female-patient",
                "--protect-tissue",
                "endometrium",
                "--format",
                "json",
            ]
        )
        == 0
    )
    row = json.loads(capsys.readouterr().out)[0]
    assert row["classification"] == "CTA" and not row["keep"]
    assert row["host_status"] in {"excluded", "unresolved"}
    assert row["host_reasons"]
    args = [
        "antigens",
        "list",
        "--profile",
        "male-patient",
        "--absent-tissue",
        "prostate",
        "--accept-loss",
        "testis",
        "--format",
        "json",
    ]
    assert cli.main(args) == 0
    meta = json.loads(capsys.readouterr().err)["host_context"]
    assert "prostate" in meta["absent"] and "testis" in meta["accept_loss"]


def test_cli_subtraction_uses_the_host_selection_background(capsys, monkeypatch):
    captured = {}
    records = pd.DataFrame({"gene_id": ["host_excluded_cta"], "sequence": ["LLLLLLLL"]})
    records.attrs["host_context"] = {"protected": ["testis"]}

    class Selection:
        def background_sequences(self, background):
            captured["background"] = background
            return records

    def select(group, **kwargs):
        captured.update(group=group, **kwargs)
        return Selection()

    monkeypatch.setattr(antigen_sets, "select", select)
    assert (
        cli.main(
            ["antigens", "subtract", "LLLLLLLL", "--profile", "male-patient", "--format", "json"]
        )
        == 0
    )
    output = capsys.readouterr()
    assert json.loads(output.out) == []
    assert captured["host_profile"] == "male-patient"
    assert captured["background"] == "Unsafe|Other"
    assert json.loads(output.err)["selection"]["host_context"]["protected"] == ["testis"]
