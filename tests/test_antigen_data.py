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
    assert set(negative.index) == {"OR52D1", "OR2T4", "OR6F1", "OR7E24"}
    withheld = antigen_data.gene_list("benchmark/withhold")
    assert withheld.gene_name.tolist() == ["OR4M2B"]
    assert negative.control_basis.eq("not-CTA").all()
    assert negative.gene_id.is_unique
    actual = target_selection.results().set_index("protein_group_label")
    assert not actual.loc[[*negative.index, *withheld.gene_name], "is_cta"].any()
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
        assert catalog.dataset(row.dataset).description
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


def test_dice_uses_native_ranks_and_bounds_unpaired_loci(monkeypatch, tmp_path):
    from scripts import bulk_normal_percentiles as bulk

    path = tmp_path / "dice.csv"
    pd.DataFrame(
        {"Feature_name": ["A.1", "B.1", "C.1", "D.1"], "1": [2, 6, 9, 1], "2": [8, 7, 10, 2]}
    ).to_csv(path, index=False)
    facts = pd.DataFrame(
        {
            "entity_id": ["single", "pair", "missing"],
            "member_count": [1, 2, 2],
            "member_gene_ids": ['["A"]', '["B","C"]', '["D","E"]'],
        }
    )
    monkeypatch.setattr(
        bulk.reference_data, "REFERENCE_SOURCES", {"dice_x": {"version": "v1", "cell_type": "T"}}
    )
    monkeypatch.setattr(bulk, "_input", lambda name, receipts: path)
    d = bulk.dice_percentiles(facts, {}).set_index("entity_id")
    assert d.loc["single", "value"] == pytest.approx(7.7) and d.loc["single", "n_observations"] == 2
    assert d.loc["pair", "value"] == 17 and pd.isna(d.loc["pair", "n_observations"])
    assert not d.loc["missing", "qc_pass"] and pd.isna(d.loc["missing", "value"])


def test_immune_percentile_sums_loci_before_each_donor_votes(monkeypatch, tmp_path):
    from scripts import bulk_normal_percentiles as bulk

    path = tmp_path / "immune.tsv"
    rows = [
        (locus, donor, sample, "NK", value)
        for donor, sample, values in [
            ("one", "a", (1, 2)),
            ("one", "b", (3, 4)),
            ("two", "c", (5, 8)),
            ("three", "d", (100, None)),
        ]
        for locus, value in zip(["A", "B"], values)
    ]
    pd.DataFrame(rows, columns=["ENSG ID", "Donor", "Sample ID", "Immune cell", "TPM"]).to_csv(
        path, sep="\t", index=False
    )
    facts = pd.DataFrame(
        {"entity_id": ["pair"], "member_count": [2], "member_gene_ids": ['["A","B"]']}
    )
    ids = pd.DataFrame({"source_gene_id": ["A", "B"], "canonical_gene_id": ["A", "B"]})
    monkeypatch.setattr(bulk, "_input", lambda name, receipts: path)
    d = bulk.hpa_immune_percentiles(facts, ids, {}).iloc[0]
    assert d.value == pytest.approx(12.6) and d.n_observations == 2


def test_baseline_panels_expose_only_their_own_scope(monkeypatch, tmp_path):
    raw = pd.DataFrame(
        {
            "entity_id": ["a"] * 4 + ["b"],
            "scope": ["normal", "disposable", "reproductive", "thymic", "normal"],
            "value": [1, 2, 1000, 30, 5],
        }
    )
    path = tmp_path / "healthy.csv"
    raw.to_csv(path, index=False)
    monkeypatch.setattr(catalog, "ensure", lambda *args, **kwargs: path)
    for panel, expected in [
        ("normal", [1]),
        ("disposable", [2]),
        ("thymic", [30]),
        ("healthy", [1, 2, 1000, 30]),
    ]:
        assert antigen_data.expression(panel, entity_ids="a").value.tolist() == expected


def test_each_locus_in_a_publication_claim_has_its_own_symbol(monkeypatch):
    claim = {"protein": "ABC group", "gene_ids": ["id1", "id2"], "reported_name": "ABC"}
    registry = {"sources": [{"source_id": "test", "url": "source", "gene_claims": [claim]}]}
    monkeypatch.setattr(antigen_data.cta_literature, "publication_registry", lambda: registry)
    monkeypatch.setattr(
        antigen_data,
        "gene_identity_catalog",
        lambda: pd.DataFrame({"canonical_gene_id": ["id1", "id2"], "canonical_symbol": ["A", "B"]}),
    )
    data = antigen_data.gene_list("paper/test")
    assert data.gene_name.tolist() == ["A", "B"]
    assert data.reported_protein.eq("ABC group").all()
