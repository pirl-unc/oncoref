"""The category audit preserves calls and ranks supported source-native values."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from oncoref import target_selection

ROOT = Path(__file__).parents[1]
OUT = ROOT / "docs/audits/target-selection-20261109/category-profiles"


@pytest.fixture(scope="module")
def profiles():
    return pd.read_csv(OUT / "all-protein-profiles.csv.gz", low_memory=False)


def test_audit_does_not_reclassify_proteins(profiles):
    current = target_selection.results().set_index("entity_id")
    observed = profiles.set_index("entity_id")
    assert observed.index.is_unique and set(observed.index) == set(current.index)
    assert observed.classification.equals(current.classification.reindex(observed.index))
    published = pd.read_csv(OUT / "published-targets.csv.gz", low_memory=False)
    assert set(published.entity_id) == set(current.index[current.published_nomination])
    assert len(published) == 2233


def test_exclusive_sets_preserve_every_gene_and_reference_record(profiles):
    for category, expected_proteins, expected_genes in [
        ("Cancer Antigen only", 260, 261),
        ("Reproductive Antigen only", 22, 22),
    ]:
        slug = category.lower().replace(" ", "-")
        proteins = pd.read_csv(OUT / f"{slug}.csv")
        genes = pd.read_csv(OUT / f"{slug}-genes.csv.gz")
        assert len(proteins) == expected_proteins and len(genes) == expected_genes
        assert genes.gene_id.is_unique
        expected = {g for ids in proteins.member_gene_ids for g in json.loads(ids)}
        assert set(genes.gene_id) == expected
        assert set(proteins.entity_id) == set(
            profiles.loc[profiles.category.eq(category), "entity_id"]
        )


def test_peaks_are_finite_and_keep_disposable_tissues(profiles):
    peaks = pd.read_csv(OUT / "normal-peaks.csv.gz", low_memory=False)
    assert np.isfinite(peaks.value).all() and peaks.value.ge(0).all()
    assert set(peaks.scope) == {"normal", "disposable"}
    assert not peaks.duplicated(["entity_id", "unit"]).any()
    cornea = peaks[peaks.unit.eq("CP10K")]
    assert cornea.source_label.eq("Adult cornea").all()
    assert cornea.source_url.eq(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683"
    ).all()
    by_name = profiles.set_index("protein")
    assert by_name.loc["LCE3C", "normal_peak_tpm"] == 21.4
    assert by_name.loc["LCE3C", "normal_peak_ntpm"] == 23.4
    assert by_name.loc["LCE3C", "normal_scope_tpm"] == "disposable"
    assert by_name.loc["PRAME", "normal_peak_ntpm"] == 20.3
    assert by_name.loc["PRAME", "normal_ceiling_ntpm"] == 32
    assert by_name.loc["PRAME", "normal_ihc_ntpm"] == "Enhanced"
    assert by_name.loc["CTAG2", "normal_peak_ntpm"] == 24.9


def test_each_top_twenty_is_ranked_within_its_own_category_and_unit(profiles):
    for category in target_selection.CATEGORIES:
        slug = category.lower().replace(" ", "-")
        for key in ["tpm", "ntpm", "cp10k"]:
            ranking = pd.read_csv(OUT / f"top20-{slug}-{key}.csv")
            field = "normal_peak_" + key
            expected = profiles[profiles.category.eq(category) & profiles[field].notna()]
            expected = expected.sort_values([field, "protein"], ascending=[False, True]).head(20)
            assert ranking.entity_id.tolist() == expected.entity_id.tolist()
            assert ranking[field].is_monotonic_decreasing
            assert len(ranking) <= 20


def test_named_examples_and_reason_for_reproductive_only(profiles):
    names = profiles.set_index("protein")
    for name in ["MAGEA4", "PRAME", "NY-ESO-1", "XAGE1A/B", "CT83"]:
        assert names.loc[name, "category"] == "CTA"
    assert names.loc["MAGEC2", "classification"] == "Other:Review"
    assert names.loc["MAGEA3", "category"] == names.loc["MAGEA12", "category"] == "Unsafe"
    reproductive = profiles[profiles.category.eq("Reproductive Antigen only")]
    assert reproductive.cancer_evidence.value_counts().to_dict() == {
        "unassessed": 15,
        "does_not_qualify": 7,
    }


def test_provenance_covers_every_committed_output():
    manifest = json.loads((OUT / "manifest.json").read_text())
    assert manifest["classifications_changed"] is False
    for section, base in [("implementation_sha256", ROOT), ("outputs_sha256", OUT)]:
        for path, expected in manifest[section].items():
            assert hashlib.sha256((base / path).read_bytes()).hexdigest() == expected, path
