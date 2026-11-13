"""Target-filter knobs: defaults are inert, modes nest, unmeasured never passes."""

import itertools

import pandas as pd
import pytest

from oncoref import target_filters as tf
from oncoref import target_selection


@pytest.fixture(scope="module")
def catalog():
    return tf.measures(target_selection.results())


def kept(f, catalog):
    return set(tf.select(f, catalog).entity_id)


def test_default_filter_is_the_plain_catalog(catalog):
    d = tf.apply(tf.TargetFilter(), catalog)
    assert d.keep.sum() == catalog.is_target_candidate.sum() == 902
    assert (d.keep == d.is_target_candidate).all()
    assert (d[d.keep].filtered_category == d[d.keep].category).all()
    assert tf.TargetFilter().changed() == {}


def test_filter_modes_are_nested(catalog):
    assert list(tf.FILTER_MODES) == ["all", "loose", "standard", "strict"]
    ladder = [kept(name, catalog) for name in tf.FILTER_MODES]
    for wider, narrower in itertools.pairwise(ladder):
        assert narrower < wider


def test_no_knob_can_add_a_protein(catalog):
    candidates = set(catalog[catalog.is_target_candidate].entity_id)
    for name in tf.FILTER_MODES:
        assert kept(name, catalog) <= candidates


def test_every_mode_keeps_the_positive_controls(catalog):
    for name in tf.FILTER_MODES:
        labels = set(tf.select(name, catalog).protein_group_label)
        assert set(tf.POSITIVE_CONTROLS) <= labels, name


def test_layers_report_first_failure_in_order(catalog):
    d = tf.apply(tf.TargetFilter(min_papers=2, max_off_target_tpm=0.5), catalog)
    both = d[d.is_target_candidate & ~d.pass_papers & ~d.pass_off_target]
    assert len(both) and (both.filter_failed == "off_target").all()
    only_papers = d[d.is_target_candidate & ~d.pass_papers & d.pass_off_target]
    assert len(only_papers) and (only_papers.filter_failed == "papers").all()


def test_off_target_knobs_are_per_scale_and_include_disposable(catalog):
    tissue = tf.select(tf.TargetFilter(max_off_target_tpm=1.0), catalog)
    assert tissue.normal_reference_tpm.lt(1.0).all() and tissue.disposable_max_tpm.lt(1.0).all()
    assert tissue.healthy_maximum_ntpm.ge(1.0).any()  # the TPM knob does not read nTPM
    both = tf.select(tf.TargetFilter(max_off_target_tpm=1.0, max_off_target_ntpm=1.0), catalog)
    assert both.healthy_maximum_ntpm.lt(1.0).all() and both.disposable_max_ntpm.lt(1.0).all()
    assert set(both.entity_id) < set(tissue.entity_id)
    # PRAME is low at tissue level and high only on a single-cell cluster.
    prame = catalog[catalog.protein_group_label == "PRAME"].iloc[0]
    assert prame.off_target_tpm < 5 < prame.off_target_ntpm
    # The catalog's reference N is never below the legacy bulk maximum.
    candidates = catalog[catalog.is_target_candidate]
    assert candidates.normal_reference_tpm.ge(candidates.healthy_maximum_tpm - 1e-9).all()


def test_evidence_knobs_move_proteins_between_sets_not_out_of_the_union(catalog):
    d = tf.apply(tf.TargetFilter(min_patients=5), catalog)
    weak = d.is_cta & d.independent_patients_above_normal.lt(5)
    assert weak.any()
    assert (d[weak].filtered_category == "Reproductive Antigen only").all()
    assert d[weak].keep.all()
    # Reproductive-Antigen-only proteins are untouched by a tumor knob.
    ra_only = d.category.eq("Reproductive Antigen only")
    assert d[ra_only].keep.all()
    # Requiring both sets is the membership knob's job.
    both = tf.apply(tf.TargetFilter(min_patients=5, membership="cta"), catalog)
    assert not both[weak].keep.any()


def test_cta_list_papers_count_only_where_they_nominate(catalog):
    papers = tf.paper_table()
    from oncoref.cta_literature import attestations

    links = attestations()
    nominated = set(
        map(tuple, links[links.published_nomination][["entity_id", "source_id"]].to_numpy())
    )
    cta = papers[papers.family == "cta"]
    assert set(map(tuple, cta[["entity_id", "source_id"]].to_numpy())) <= nominated
    candidates = catalog[catalog.is_target_candidate]
    assert candidates.papers_cta.gt(0).sum() == candidates.published_nomination.astype(bool).sum()


def test_exact_members_only_removes_partial_nominations_and_nothing_else(catalog):
    d = tf.apply(tf.TargetFilter(exact_members_only=True), catalog)
    removed = d[d.is_target_candidate & ~d.keep]
    assert len(removed) == 6
    assert (removed.nomination_scope == "some_members").all()


def test_missing_values_never_pass_a_threshold(catalog):
    row = catalog[catalog.protein_group_label == "MAGEA4"].copy()
    assert tf.apply(tf.TargetFilter(max_off_target_tpm=1.0), row).keep.all()
    row["normal_reference_tpm"] = float("nan")
    d = tf.apply(tf.TargetFilter(max_off_target_tpm=1.0), row)
    assert d.filter_failed.iloc[0] == "off_target"


def test_bad_knobs_are_rejected():
    for bad in (
        {"membership": "everything"},
        {"families": ("nope",)},
        {"rare_cell": "some"},
        {"mapping": "clean"},
        {"normal_function": "no_disposable"},
        {"readthrough": "drop"},
        {"min_families": 5},
        {"exclude_prefixes": "MAGE"},
    ):
        with pytest.raises(ValueError):
            tf.TargetFilter(**bad)


def test_describe_states_every_active_knob():
    text = tf.describe("strict")
    for phrase in (
        "below 4 TPM",
        "at least 5 independent tumors",
        "10 x",
        "in cta",
        "at least 3 distinct paper",
        "MAGE",
        "except MAGEA4",
    ):
        assert phrase in text
    assert len(tf.describe("all").splitlines()) == 2


def test_mage_rule_drops_the_family_but_keeps_magea4(catalog):
    f = tf.TargetFilter(exclude_prefixes=("MAGE",), except_labels=("MAGEA4",))
    d = tf.apply(f, catalog)
    removed = d[d.is_target_candidate & ~d.keep]
    assert len(removed) == 17
    assert removed.protein_group_label.str.startswith("MAGE").all()
    assert (removed.filter_failed == "exclusions").all()
    assert d[d.protein_group_label == "MAGEA4"].keep.all()
    assert isinstance(d, pd.DataFrame)


def test_readthrough_knob_removes_the_listed_loci_and_nothing_else(catalog):
    from oncoref import target_curation

    d = tf.apply(tf.TargetFilter(readthrough="exclude"), catalog)
    removed = set(d[d.is_target_candidate & ~d.keep].protein_group_label)
    assert removed == set(target_curation.readthrough().protein_group_label)
    assert {"ATF7-NPFF", "ENSG00000273155", "ENSG00000267618"} <= removed
    assert (d[d.is_target_candidate & ~d.keep].filter_failed == "curation").all()


def test_normal_function_knob_follows_the_tier_order(catalog):
    sizes = [
        int(tf.apply(tf.TargetFilter(normal_function=value), catalog).keep.sum())
        for value in tf.NORMAL_FUNCTION
    ]
    assert sizes[0] == 902 and sizes == sorted(sizes, reverse=True) and sizes[-1] < sizes[1]
    strict = tf.apply(tf.TargetFilter(normal_function="no_protected"), catalog)
    kept = strict[strict.keep]
    assert not kept.curated_tier.isin(["vital", "protected"]).any()
    # Hair follicle and cutaneous melanocytes are Disposable in the registry, so
    # hair keratins and PRAME are recorded but never removed by this knob.
    loosest = tf.apply(tf.TargetFilter(normal_function="no_manageable"), catalog)
    disposable = loosest[loosest.is_target_candidate & loosest.curated_tier.eq("disposable")]
    assert len(disposable) and disposable.keep.all()
    assert loosest[loosest.protein_group_label == "PRAME"].curated_tier.item() == "disposable"


def test_every_acting_curation_record_keeps_a_citation_and_a_tier():
    from oncoref import target_curation as tc

    found = tc.findings()
    acting = found[found.acts.astype(str).str.lower() == "true"]
    assert len(acting) and acting.best_source.str.len().gt(0).all()
    assert acting.tier.isin(tc.TIERS).all()
    assert acting.verdict.eq("confirmed").all()
    families = tc.families()
    assert families[families.acts].sources.str.len().gt(0).all()
    compartments = tc.compartments()
    proposed = compartments[compartments.tier_source == "proposed"]
    assert len(proposed) == 7 and proposed.tier_references.str.contains("PMID").all()
    resolved = compartments[compartments.tier_source != "proposed"]
    assert resolved.tier_source.str.startswith("tissue_hazards:").all()
