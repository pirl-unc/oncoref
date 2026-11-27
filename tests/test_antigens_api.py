# Licensed under the Apache License, Version 2.0.
"""Tests for the unified oncoref.antigens API, CLI upgrades, pre-baked catalog, and explorer."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import oncoref
from oncoref import antigens, cli, target_selection

ROOT = Path(__file__).resolve().parents[1]


def test_antigens_top_level_export_and_summary():
    assert hasattr(oncoref, "antigens")
    summ = antigens.summary()
    assert list(summ["species"]) == ["human", "dog", "mouse"]
    by_sp = summ.set_index("species")

    assert by_sp.loc["human", "CTA"] == 322
    assert by_sp.loc["human", "CTA_ectopic"] == 311
    assert (
        by_sp.loc["human", "CTA_M"]
        + by_sp.loc["human", "CTA_P"]
        + by_sp.loc["human", "CTA_F"]
        + by_sp.loc["human", "CTA_MF"]
        == 322
    )
    assert by_sp.loc["dog", "CTA"] == 118
    assert by_sp.loc["dog", "tumor_cohorts"] == 19
    assert by_sp.loc["dog", "cell_lines"] == 12
    assert by_sp.loc["mouse", "CTA"] == 8
    assert by_sp.loc["mouse", "tumor_cohorts"] == 14
    assert by_sp.loc["mouse", "cell_lines"] == 9


def test_antigens_sex_protection_differs_from_empty_legacy_profile():
    # sex="female" protects all non-pregnant, non-lactating adult female reproductive
    # and resting mammary organs AND requires >=10x effective tumor fold over protected
    # female exposure -> 258 of the 270 (CTA-M + CTA-P) CTAs pass both organ ceiling and >=10x margin
    fem = antigens.table("CTA", sex="female")
    assert len(fem) == 258
    assert set(fem["sex_restriction_code"]) == {"CTA-M", "CTA-P"}
    assert (fem["effective_tumor_fold"] >= 10.0).all()
    assert "MAGEA4" in set(fem["protein_group_label"])
    assert "DUX4" in set(fem["protein_group_label"])
    assert "LALBA" not in set(fem["protein_group_label"])
    assert "CSN2" not in set(fem["protein_group_label"])

    # sex="male" protects all adult male reproductive & resting mammary organs + >=10x margin -> 58 CTAs
    # (57 CTA-F / CTA-P plus FOXR2 at 4.96 TPM testis < 8 TPM Protected ceiling with 28.8x tumor fold)
    male = antigens.table("CTA", sex="male")
    assert len(male) == 58
    assert set(male["sex_restriction_code"]) == {"CTA-F", "CTA-P", "CTA-M"}
    assert set(male.loc[male["sex_restriction_code"].eq("CTA-M"), "protein_group_label"]) == {
        "FOXR2"
    }
    assert (male["effective_tumor_fold"] >= 10.0).all()
    assert "LALBA" in set(male["protein_group_label"])
    assert "CSN2" in set(male["protein_group_label"])
    assert "MAGEA4" not in set(male["protein_group_label"])

    # Surgical absence / accepted loss unlocks homologous organ targets while preserving XY resting_breast
    male_post_prostatectomy = antigens.table(
        "CTA", sex="male", absent=("prostate", "seminal vesicle", "epididymis", "testis")
    )
    male_post_labels = set(male_post_prostatectomy["protein_group_label"])
    assert {"SEMG2", "CSN2", "LALBA"} <= male_post_labels

    # Conflicting compartment choices & male pregnancy/lactation fail closed
    with pytest.raises(ValueError, match="conflict"):
        antigens.table("CTA", protect=("ovary",), accept_loss=("ovary",))
    with pytest.raises(ValueError, match="incompatible"):
        antigens.table("CTA", sex="male", pregnant=True)
    # pregnant=True respects explicit surgical absent / accept_loss without crashing
    preg_no_lact = antigens._build_host_context(
        species="human", sex="female", pregnant=True, absent=("lactating_breast",)
    )
    assert preg_no_lact is not None and "lactating_breast" in preg_no_lact.absent


def test_antigens_ectopic_and_replicated_filters_across_species():
    ect_human = antigens.table("CTA", ectopic_only=True)
    assert len(ect_human) == 311
    assert "SEMG2" not in set(ect_human["protein_group_label"])
    assert "CSN2" in set(ect_human["protein_group_label"])

    rep_human = antigens.table("CTA", replicated_only=True)
    assert (rep_human["tumor_support"] == "replicated").all()
    assert 0 < len(rep_human) < 322

    dog_rep = antigens.table("CTA", species="dog", ectopic_only=True, replicated_only=True)
    assert not dog_rep.empty
    assert (dog_rep["tumor_support"] == "replicated").all()
    assert set(dog_rep["tumor_activation_mode"]) <= {"ectopic", "ectopic_and_lineage_retained"}

    mouse_ect = antigens.table("CTA", species="mouse", ectopic_only=True)
    assert len(mouse_ect) == 6
    assert "Pbsn" not in set(mouse_ect["protein_group_label"])
    assert "Hmga2" in set(mouse_ect["protein_group_label"])

    gids = antigens.gene_ids("CTA-P")
    assert len(gids) >= 48

    with pytest.raises(ValueError, match="Candidate filters"):
        antigens.table("Unsafe", sex="female")


def test_antigens_explain_and_funnel():
    prame = antigens.explain("PRAME", sex="female")
    # Strict JSON compliance (no NaN values anywhere in output dict)
    assert json.dumps(prame, allow_nan=False)
    assert prame["category"] == "CTA"
    assert prame["sex_subset"] == "CTA-MF"
    assert prame["host"]["pass_host"] is False
    assert prame["clinical_tcr_benchmarks"]

    # Cross-recognized normal off-target protein (TTN) surfaces its TCR benchmark
    ttn = antigens.explain("TTN")
    assert json.dumps(ttn, allow_nan=False)
    assert ttn["clinical_tcr_benchmarks"]

    # DKKL1 passes organ ceiling (<8 TPM) on sex="female" but fails >=10x effective_tumor_fold,
    # matching its exclusion from antigens.table("CTA", sex="female")
    dkkl1 = antigens.explain("DKKL1", sex="female")
    assert dkkl1["host"]["organ_ceiling_pass"] is True
    assert dkkl1["host"]["keep"] is False
    assert dkkl1["host"]["pass_host"] is False

    magea4 = antigens.explain("MAGEA4", sex="female")
    assert magea4["sex_subset"] == "CTA-M"
    assert magea4["host"]["pass_host"] is True
    assert magea4["host"]["keep"] is True

    hmga2 = antigens.explain("Hmga2", species="mouse")
    assert hmga2["category"] == "CTA"
    assert hmga2["cell_line_peak_tpm"] is not None and hmga2["cell_line_peak_tpm"] > 10.0

    fun = antigens.funnel("human", sex="female", ectopic_only=True, replicated_only=True)
    rem = list(fun["remaining"])
    assert rem == sorted(rem, reverse=True)
    assert rem[0] == 19893


def test_antigens_subtract_symbol_long_label_and_invalid_queries(monkeypatch):
    bg = pd.DataFrame({"gene_id": ["bg1"], "sequence": ["LLLLLLLLL"]})
    monkeypatch.setattr(oncoref.antigen_sets, "sequences", lambda *a, **kw: bg)
    monkeypatch.setattr(
        oncoref.antigen_sets,
        "_sequences",
        lambda rows, *a, **kw: pd.DataFrame(
            {
                "transcript_id": ["ENST_A", "ENST_B"],
                "sequence": ["MAGEAFIRSTSEQ", "MAGEASECNDSEQ"],
            }
        ),
    )
    # >=30-char catalog label resolves cleanly and preserves both isoforms
    sub_ct47 = antigens.subtract("CT47A1/10/11/12/2/3/4/5/6/7/8/9", lengths=[9])
    assert set(sub_ct47["query_id"]) == {"ENST_A", "ENST_B"}

    # Unresolved symbols/aliases or non-AA identifiers fail closed with ValueError
    for bad_q in ("CTAG1A/B", "MAGE-A4", "PRAM", "METASTASIS_1"):
        with pytest.raises(ValueError, match="neither a known gene"):
            antigens.subtract(bad_q, lengths=[9])


def test_antigens_cli_summary_sex_ectopic_and_replicated(capsys):
    assert cli.main(["antigens", "summary", "--format", "json"]) == 0
    out_summ = json.loads(capsys.readouterr().out)
    assert len(out_summ) == 3

    assert (
        cli.main(
            [
                "antigens",
                "list",
                "--sex",
                "female",
                "--ectopic-only",
                "--replicated",
                "--format",
                "json",
            ]
        )
        == 0
    )
    cap = capsys.readouterr()
    rows = json.loads(cap.out)
    meta = json.loads(cap.err)
    assert rows
    assert all(r["sex_restriction_code"] in {"CTA-M", "CTA-P"} for r in rows)
    assert all(r["tumor_support"] == "replicated" for r in rows)
    assert "ovary" in meta["host_context"]["protected"]

    # Contradictory --sex and --profile fails closed
    assert cli.main(["antigens", "list", "--sex", "female", "--profile", "male-patient"]) == 1
    assert "not both" in capsys.readouterr().err

    # CLI funnel supports --ectopic-only and --replicated
    assert (
        cli.main(
            [
                "antigens",
                "funnel",
                "--sex",
                "female",
                "--ectopic-only",
                "--replicated",
                "--format",
                "json",
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)

    # Table output has clean stderr (no noisy JSON metadata dump)
    assert cli.main(["antigens", "list", "--set", "CTA-F", "--format", "table"]) == 0
    cap_tbl = capsys.readouterr()
    assert "LALBA" in cap_tbl.out
    assert cap_tbl.err == ""


def test_packaged_human_catalog_has_prebaked_columns_without_sidecars():
    with (ROOT / "oncoref/data/target-selection.csv.gz").open("rb") as fh:
        raw = pd.read_csv(fh, compression="gzip", nrows=5)
    for col in (
        "sex_restriction_code",
        "tumor_activation_mode",
        "tumor_support",
        "clinical_clearance_scope",
        "closed_homolog_family_id",
        "exact_unique_9mers",
        "safe_2mm_9mers",
        "unsafe_near_self_1mm_pct",
    ):
        assert col in raw.columns, f"Expected {col} pre-baked in target-selection.csv.gz"

    ts = target_selection.results()
    assert set(ts.loc[ts["is_cancer_antigen"], "tumor_support"]) == {"replicated", "single_sample"}
    assert (ts.loc[~ts["is_cancer_antigen"], "tumor_support"] == "").all()


def test_interactive_antigen_explorer_html_is_live_and_current():
    explorer = ROOT / "docs/antigen-explorer.html"
    assert explorer.is_file()
    text = explorer.read_text(encoding="utf-8")
    assert 'id="explorerData"' in text and 'id=" explorerData"' not in text
    assert target_selection.POLICY_ID in text
    assert '"human":' in text and '"dog":' in text and '"mouse":' in text
    for gene in ("MAGEA4", "PRAME", "LALBA", "CSN2", "DUX4", "NT5C1B", "Hmga2", "Trap1a"):
        assert gene in text
