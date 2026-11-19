"""Public source access, independent literature controls and normal percentile contracts."""

import json
from pathlib import Path

import pandas as pd
import pytest

from oncoref import antigen_data, catalog, reference_data, target_selection

ROOT = Path(__file__).resolve().parents[1]


def test_cited_lists_keep_identity_and_scope():
    lists = antigen_data.gene_lists()
    assert lists.list_id.is_unique
    negative = antigen_data.gene_list("benchmark/not-CTA").set_index("gene_name")
    assert set(negative.index) == {"OR4M2B", "OR52D1", "OR2T4", "OR6F1", "OR7E24"}
    assert negative.at["OR4M2B", "control_basis"] == "withhold"
    assert negative.gene_id.is_unique
    actual = target_selection.results().set_index("protein_group_label")
    assert not actual.loc[negative.index, "is_cta"].any()
    positive = antigen_data.gene_list("benchmark/CTA")
    assert actual.loc[positive.gene_name, "is_cta"].all()
    claim = antigen_data.gene_list("paper/10.1186/s40425-019-0678-x").iloc[0]
    assert claim.gene_id == "ENSG00000204019" and "tumor_support" in claim.allowed_uses
    table = catalog.read_table("marcinkowski2019-ct83-expression")
    ih = table[table.assay.eq("RNAscope") & table.statistic.eq("any positive tumor cells")]
    assert list(zip(ih.positive_specimens, ih.assessed_specimens)) == [(9, 13), (4, 9)]
    assert catalog.ensure("target-selection").name == "target-selection.csv.gz"
    assert catalog.ensure("target-selection.json").name == "target-selection.json"
    assert antigen_data.download("cta-literature-benchmarks.csv").is_file()
    with pytest.raises(ValueError, match="Choose a version"):
        antigen_data.download("gtex_individual_tissue_tpm")


def test_all_catalog_inputs_resolve_to_the_recorded_source_digest():
    sources = antigen_data.sources()
    assert not sources.duplicated(["dataset", "version"]).any()
    assert sources.sha256.str.fullmatch("[a-f0-9]{64}").all()
    assert sources.url.str.startswith("https://").all()
    for row in sources.itertuples():
        version = row.version or None
        assert catalog.resolve_version(row.dataset, version)
        if row.dataset in reference_data.REFERENCE_SOURCES:
            pin = reference_data.REFERENCE_SOURCES[row.dataset].get("pins", {}).get(version)
            if pin:
                assert pin["sha256"] == row.sha256
    indexed = set(sources.sha256)
    inputs = json.loads(
        (ROOT / "docs/audits/target-selection-20261110/rebuild/inputs.json").read_text()
    )
    assert {v for k, v in inputs.items() if not k.endswith(".py")} <= indexed


def test_comparison_and_atlas_views_preserve_native_rows(monkeypatch):
    sources = antigen_data.sources().set_index("dataset")

    def local(name, version=None):
        return ROOT / sources.loc[name, "source_file"]

    monkeypatch.setattr(catalog, "ensure", local)
    entity = target_selection.results("CT83").entity_id.item()
    observed = antigen_data.expression("tumor", entity_ids=entity)
    assert set(observed.entity_id) == {entity} and set(observed.unit) == {"TPM"}
    assert observed.loc[observed.single_sample_pass, "participant_id"].nunique() == 22
    atlases = antigen_data.cell_atlases()
    assert atlases.dataset_version_id.is_unique and len(atlases) == 452


def test_biopsy_p95_and_anonymous_witnesses_are_distinct():
    from scripts.bulk_normal_percentiles import tissue_percentiles

    raw = pd.DataFrame(
        {
            "entity_id": ["a"] * 4,
            "value": [2.0, 8.0, 100.0, 23.0],
            "tissue": ["lung"] * 3 + ["testis"],
            "scope": ["normal"] * 3 + ["reproductive"],
            "sample_id": ["one", "two", "three", "four"],
        }
    )
    frame = raw.assign(
        source_id="hpa_tissue_samples:v23",
        identity_resolved=[True, True, False, True],
        qc_pass=True,
        participant_id="",
        library_qc_review="",
        mammary_state_inference="",
    )
    result = tissue_percentiles(frame).set_index("tissue")
    assert result.loc["lung", "value"] == pytest.approx(7.7)
    assert result.loc["lung", "n_observations"] == 2 and result.loc["lung", "sample_id"] == ""
    assert result.loc["testis", "sample_id"] == "four" and result.loc["testis", "value"] == 23
