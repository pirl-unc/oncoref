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
    assert cli.main(["antigens", "explain", "CTAG2", "--format", "json"]) == 0
    row = json.loads(capsys.readouterr().out)[0]
    assert row["classification"] == "Other:Review"
    assert row["rare_cell_status"] == "unresolved"
    assert row["strict"].startswith("fail:")
    assert cli.main(["antigens", "explain", "SSX2", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["protein_group_label"] == "SSX2/B"
    assert cli.main(["antigens", "explain", "not-a-gene"]) == 1


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
    assert (
        cli.main(
            [
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
        )
        == 0
    )
    lactating_rows = json.loads(capsys.readouterr().out)
    assert len(lactating_rows) == 283
    assert "LALBA" not in {r["protein_group_label"] for r in lactating_rows}

    assert (
        cli.main(
            [
                "antigens",
                "list",
                "--set",
                "CTA",
                "--profile",
                "male-patient",
                "--allow-tissue",
                "prostate",
                "--allow-tissue",
                "seminal vesicle",
                "--format",
                "json",
            ]
        )
        == 0
    )
    male_post_rp = json.loads(capsys.readouterr().out)
    assert len(male_post_rp) == 264
    assert {"LALBA", "TGM4", "SEMG2"} <= {r["protein_group_label"] for r in male_post_rp}

    assert cli.main(["antigens", "explain", "LALBA", "--format", "json"]) == 0
    explained = json.loads(capsys.readouterr().out)[0]
    assert explained["cta_tier"] == "extended_only"
    assert explained["reproductive_qualifying_compartments"] == "lactating_breast"
    assert explained["profile:male-patient"] == "pass"
    assert explained["profile:female-patient"] == "pass"
    assert explained["profile:female-lactating"].startswith("fail:")
    assert explained["profile:core"].startswith("fail:")
