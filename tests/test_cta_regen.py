# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""CTA regeneration from HPA — parity with the shipped table (#35, Phase C)."""

import io
import os
from pathlib import Path

import pandas as pd
import pytest

from oncoref import cta_regen, hpa, reference_data

_CSV = Path(__file__).resolve().parents[1] / "oncoref" / "data" / "cancer-testis-antigens.csv"

_HPA_READY = reference_data.is_cached("hpa_rna_consensus") and reference_data.is_cached(
    "hpa_normal_tissue"
)


def _require_hpa_v23():
    if _HPA_READY:
        return
    if os.environ.get("CI"):
        pytest.fail("CI must stage the pinned HPA v23 RNA and normal-tissue tables")
    pytest.skip("HPA v23 not cached")


def test_fraction_deflation_semantics():
    # All-below-1 nTPM -> deflated total 0 -> 1.0 (restricted by the +1 pseudocount).
    f = cta_regen._fraction
    repro = frozenset({"testis"})
    assert f({"testis": 0.5, "liver": 0.4}, repro, deflate=True) == 1.0
    assert f({"testis": 0.5, "liver": 0.4}, repro, deflate=False) < 1.0
    # A real somatic signal lowers the deflated reproductive fraction.
    val = f({"testis": 10.0, "liver": 10.0}, repro, deflate=True)
    assert 0.0 < val < 1.0


@pytest.mark.parametrize("tissue_scope", ["core", "extended"])
def test_regeneration_reproduces_shipped_table(tissue_scope):
    _require_hpa_v23()

    # Regenerating the HPA-derived columns from HPA v23 must reproduce the shipped,
    # oncoref-owned table exactly — proof the regenerator is the source of truth.
    path = _CSV if tissue_scope == "core" else _CSV.with_name("cancer-testis-antigens-extended.csv")
    regen = cta_regen.regenerate_cta_columns(pd.read_csv(path), tissue_scope=tissue_scope)
    buf = io.StringIO()
    regen.to_csv(buf, index=False)
    from_regen = pd.read_csv(io.StringIO(buf.getvalue()), dtype=str, keep_default_na=False)
    shipped = pd.read_csv(path, dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(from_regen, shipped)


def test_core_and_extended_snapshots_share_candidates_and_protein_evidence():
    core = pd.read_csv(_CSV)
    extended = pd.read_csv(_CSV.with_name("cancer-testis-antigens-extended.csv"))
    cols = cta_regen.PRESERVED_COLUMNS + cta_regen._PROTEIN_COLUMNS
    pd.testing.assert_frame_equal(core[cols], extended[cols])
    observed = core.rna_deflated_reproductive_frac.notna()
    assert observed.equals(extended.rna_deflated_reproductive_frac.notna())
    assert (
        extended.loc[observed, "rna_deflated_reproductive_frac"]
        >= core.loc[observed, "rna_deflated_reproductive_frac"]
    ).all()
    with pytest.raises(ValueError, match="Unknown CTA tissue scope"):
        cta_regen.regenerate_cta_columns(core, tissue_scope="typo")


def test_regeneration_preserves_identity_columns():
    _require_hpa_v23()

    old = pd.read_csv(_CSV)
    regen = cta_regen.regenerate_cta_columns(old)
    for col in ("Symbol", "Ensembl_Gene_ID", "Aliases", "source_databases", "biotype"):
        pd.testing.assert_series_equal(regen[col], old[col], check_dtype=False)


def test_hpa_v23_safety_tissue_mappings_use_native_ihc_labels():
    _require_hpa_v23()
    native_labels = set(hpa.hpa_normal_tissue("v23")["Tissue"].dropna().astype(str))
    for group in ("brain", "heart", "liver", "lung", "pancreas"):
        resolution = hpa.resolve_safety_tissue_group(group, require_complete=False)
        assert resolution.source_version == "v23"
        assert set(resolution.source_tissues) <= native_labels


def test_rna_only_below_floor_confidence_is_capped():
    # tsarina#114 cap: no gene with no protein and rna_max_ntpm < 2.0 may keep a
    # HIGH restriction_confidence (near-noise RNA must not earn HIGH).
    df = pd.read_csv(_CSV, dtype=str, keep_default_na=False)
    rna_max = pd.to_numeric(df["rna_max_ntpm"], errors="coerce")
    offenders = df[
        (df["restriction_confidence"] == "HIGH")
        & (df["protein_restriction"] == "NO_DATA")
        & (rna_max < 2.0)
    ]
    assert offenders.empty, sorted(offenders["Symbol"])


def test_reproductive_tier_audit_separates_accessory_expression_and_missing_data(monkeypatch):
    monkeypatch.setattr(
        hpa,
        "hpa_rna_consensus",
        lambda: pd.DataFrame(
            {
                "Gene": ["ENSG1", "ENSG1", "ENSG1"],
                "Tissue": ["placenta", "epididymis", "smooth muscle"],
                "nTPM": [101.0, 101.0, 3.0],
            }
        ),
    )
    table = pd.DataFrame(
        [
            {
                "Symbol": symbol,
                "Ensembl_Gene_ID": gid,
                "biotype": "protein_coding",
                "protein_reliability": "Enhanced",
                "protein_reproductive": True,
                "rna_max_somatic_tissue": "smooth muscle",
                "rna_max_somatic_ntpm": 3.0,
            }
            for symbol, gid in [("observed", "ENSG1"), ("missing", "ENSG2")]
        ]
    )
    audit = cta_regen.reproductive_tier_audit(table).set_index("Symbol")
    observed = audit.loc["observed"]
    assert observed.core_deflated_frac == round(100 / 202, 4)
    assert observed.extended_deflated_frac == round(200 / 202, 4)
    assert not observed.core_hpa_pass and observed.extended_hpa_pass
    assert observed.max_outside_core_tissue == "epididymis"
    assert observed.max_somatic_tissue == "smooth muscle"
    missing = audit.loc["missing"]
    assert not missing.rna_observed
    assert pd.isna(missing.core_deflated_frac)
    assert not missing.core_hpa_pass and not missing.extended_hpa_pass
    assert missing.max_somatic_tissue == ""
    assert pd.isna(missing.max_somatic_ntpm)


@pytest.mark.parametrize("fraction", [float("nan"), "nan", None])
def test_missing_fraction_never_passes_the_hpa_gate(fraction):
    row = pd.Series(
        {
            "biotype": "protein_coding",
            "protein_reliability": "Enhanced",
            "protein_reproductive": True,
            "rna_deflated_reproductive_frac": fraction,
        }
    )
    assert cta_regen.restriction_failure_reasons(row) == "missing_rna"
    assert not cta_regen._passes_filters_rule(row, 0.97)


def test_extended_gate_failure_names_its_scope():
    row = pd.Series(
        {
            "biotype": "protein_coding",
            "protein_reliability": "Enhanced",
            "protein_reproductive": True,
            "rna_deflated_reproductive_frac": 0.5,
            "rna_reproductive_tissue_scope": "extended",
        }
    )
    assert cta_regen.restriction_failure_reasons(row) == "extended_rna_fraction_below_threshold"


@pytest.mark.parametrize("source_scope", ["core", "extended", "legacy"])
@pytest.mark.parametrize("target_scope", ["core", "extended"])
@pytest.mark.parametrize("protein_observed", [False, True])
def test_regeneration_clears_missing_rna_instead_of_reusing_old_evidence(
    monkeypatch, source_scope, target_scope, protein_observed
):
    path = (
        _CSV.with_name("cancer-testis-antigens-extended.csv")
        if source_scope == "extended"
        else _CSV
    )
    table = pd.read_csv(path)
    table = table[table.Symbol.isin(["RNASE12", "PAGE4"])].copy()
    if source_scope == "legacy":
        table = table.drop(columns="rna_reproductive_tissue_scope")
    original = table.copy(deep=True)
    page4_id = table.set_index("Symbol").loc["PAGE4", "Ensembl_Gene_ID"].split(".")[0]
    rnase12_id = table.set_index("Symbol").loc["RNASE12", "Ensembl_Gene_ID"].split(".")[0]
    # RNASE12 is entirely absent; PAGE4 has measured zero expression. These
    # must stay distinguishable, even when neither belongs in the final panel.
    monkeypatch.setattr(
        hpa,
        "hpa_rna_consensus",
        lambda: pd.DataFrame(
            {"Gene": [page4_id] * 3, "Tissue": ["testis", "epididymis", "liver"], "nTPM": [0.0] * 3}
        ),
    )
    monkeypatch.setattr(
        hpa,
        "hpa_normal_tissue",
        lambda: pd.DataFrame(
            [[rnase12_id, "epididymis", "High", "Enhanced"]] if protein_observed else [],
            columns=["Gene", "Tissue", "Level", "Reliability"],
        ),
    )
    with pytest.warns(UserWarning, match="absent from HPA.*RNASE12"):
        result = cta_regen.regenerate_cta_columns(table, tissue_scope=target_scope)
    pd.testing.assert_frame_equal(table, original)
    rows = result.set_index("Symbol")
    missing = rows.loc["RNASE12"]
    assert missing.rna_reproductive_tissue_scope == target_scope
    measurements = [
        c for c in result if c.startswith("rna_") and (c.endswith("_frac") or c.endswith("_ntpm"))
    ] + ["rna_somatic_detected_count"]
    assert missing[measurements].isna().all()
    assert not missing.passes_filters and not missing.rna_reproductive
    assert not missing[[f"rna_{pct}_pct_filter" for pct, _ in cta_regen._PCT_FILTERS]].any()
    assert missing.rna_restriction == missing.rna_restriction_level == "NO_DATA"
    assert cta_regen.restriction_failure_reasons(missing) == "missing_rna"
    assert missing.protein_reliability == ("Enhanced" if protein_observed else "no data")
    assert missing.rna_max_somatic_tissue == ""
    assert pd.isna(missing.rna_thymus)
    observed = rows.loc["PAGE4"]
    assert observed.rna_max_ntpm == observed.rna_testis_ntpm == 0
    assert observed.rna_somatic_detected_count == 0
    assert observed.rna_reproductive_frac == 0
    assert observed.rna_deflated_reproductive_frac == 1
    assert observed.never_expressed

    # Persisted results remain auditable and a second regeneration is stable.
    from oncoref.cta_sources import intake_membership

    reloaded = pd.read_csv(io.StringIO(result.to_csv(index=False)))
    intake = intake_membership(table=reloaded, tissue_scope=target_scope)
    rnase12 = intake[intake.Symbol == "RNASE12"]
    assert not rnase12.empty
    assert not rnase12.hpa_restriction.any()
    assert not rnase12.default_panel.any()
    with pytest.warns(UserWarning, match="absent from HPA.*RNASE12"):
        repeated = cta_regen.regenerate_cta_columns(reloaded, tissue_scope=target_scope)
    assert repeated.to_csv(index=False) == result.to_csv(index=False)


def test_fraction_is_stable_across_tissue_order_and_python_sum_versions():
    # These values expose the difference between ordinary iterative summation
    # and compensated summation (Python 3.12 changed built-in sum()).
    from math import fsum

    values = {"testis": 1.0, **{f"t{i}": 0.1 for i in range(50)}}
    expected = 1.0 / fsum(values.values())
    assert cta_regen._fraction(values, frozenset({"testis"}), False) == expected
    assert (
        cta_regen._fraction(dict(reversed(list(values.items()))), frozenset({"testis"}), False)
        == expected
    )
