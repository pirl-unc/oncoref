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
    assert d.keep.sum() == catalog.is_target_candidate.sum() == 534
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


def test_positive_controls_follow_the_same_thresholds(catalog):
    expected = {
        "all": set(tf.POSITIVE_CONTROLS),
        "loose": set(tf.POSITIVE_CONTROLS),
        "standard": {"MAGEA4", "NY-ESO-1"},
        "strict": {"MAGEA4", "NY-ESO-1"},
    }
    for name in tf.FILTER_MODES:
        labels = set(tf.select(name, catalog).protein_group_label)
        assert set(tf.POSITIVE_CONTROLS) & labels == expected[name], name


def test_layers_report_first_failure_in_order(catalog):
    d = tf.apply(tf.TargetFilter(min_papers=2, max_off_target_tpm=0.5), catalog)
    both = d[d.is_target_candidate & ~d.pass_papers & ~d.pass_off_target]
    assert len(both) and (both.filter_failed == "off_target").all()
    only_papers = d[d.is_target_candidate & ~d.pass_papers & d.pass_off_target]
    assert len(only_papers) and (only_papers.filter_failed == "papers").all()


def test_off_target_knob_reads_tissue_tpm_and_includes_disposable(catalog):
    tissue = tf.select(tf.TargetFilter(max_off_target_tpm=1.0), catalog)
    assert tissue.normal_reference_tpm.lt(1.0).all() and tissue.disposable_max_tpm.lt(1.0).all()
    assert tissue.healthy_maximum_ntpm.ge(1.0).any()  # the knob does not read nTPM
    # Reconstruct the filter from the two current source-native references.
    candidates = catalog[catalog.is_target_candidate]
    expected = candidates.normal_reference_tpm.lt(1.0) & candidates.disposable_max_tpm.lt(1.0)
    assert set(tissue.entity_id) == set(candidates.loc[expected, "entity_id"])


def test_every_cancer_antigen_has_a_tumor_at_tenfold(catalog):
    antigens = catalog[catalog.is_cancer_antigen]
    assert antigens.tumor_fold.ge(target_selection.TUMOR_FOLD).all()
    floor = antigens[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1).clip(lower=1)
    assert (antigens.tumor_max_tpm / floor - antigens.tumor_fold).abs().lt(1e-6).all()
    assert (
        not catalog[catalog.is_target_candidate]
        .tumor_top_cancer_code.isin(target_selection.GERM_CELL_TUMORS)
        .any()
    )


def test_fold_knob_moves_proteins_between_sets_not_out_of_the_union(catalog):
    d = tf.apply(tf.TargetFilter(min_tumor_fold=50.0), catalog)
    weak = d.is_cta & d.tumor_fold.lt(50)
    assert weak.any()
    assert (d[weak].filtered_category == "Reproductive Antigen only").all()
    assert d[weak].keep.all()
    # Reproductive-Antigen-only proteins are untouched by a tumor knob.
    ra_only = d.category.eq("Reproductive Antigen only")
    assert d[ra_only].keep.all()
    # Requiring both sets is the membership knob's job.
    both = tf.apply(tf.TargetFilter(min_tumor_fold=50.0, membership="cta"), catalog)
    assert not both[weak].keep.any()


def test_reproductive_fold_tightens_its_own_set(catalog):
    d = tf.apply(tf.TargetFilter(min_reproductive_fold=50), catalog)
    weak = d.is_cta & d.reproductive_fold.lt(50)
    assert weak.any()
    assert d[weak].filtered_category.eq("Cancer Antigen only").all()
    assert d[weak].keep.all()
    assert not tf.apply(tf.TargetFilter(min_reproductive_fold=50, membership="cta"), catalog)[
        weak
    ].keep.any()


def test_patient_replication_uses_tenfold_known_patients_not_above_normal(catalog):
    row = catalog[catalog.protein_group_label.eq("MAGEA4")].copy()
    f = tf.TargetFilter(membership="cancer_antigen", min_tumor_patients=2)
    row["independent_patients_at_tenfold"] = 1
    row["independent_patient_tumor_folds"] = "[100]"
    row["independent_patients_above_normal"] = 100
    row["anonymous_samples_at_tenfold"] = 100
    assert not tf.apply(f, row).keep.any()
    row["independent_patients_at_tenfold"] = 2
    row["independent_patient_tumor_folds"] = "[100, 10]"
    assert tf.apply(f, row).keep.all()
    assert not tf.apply(
        tf.TargetFilter(membership="cancer_antigen", min_tumor_patients=2, min_tumor_fold=30), row
    ).keep.any()
    row["independent_patient_tumor_folds"] = "[100, 30]"
    assert tf.apply(
        tf.TargetFilter(membership="cancer_antigen", min_tumor_patients=2, min_tumor_fold=30), row
    ).keep.all()
    row["independent_patients_at_tenfold"] = float("nan")
    row["independent_patient_tumor_folds"] = float("nan")
    assert not tf.apply(f, row).keep.any()
    assert "known patients" in tf.describe(f)


def test_cta_list_papers_count_only_where_they_nominate(catalog):
    papers = tf.paper_table()
    from oncoref.cta_literature import current_attestations

    links = current_attestations()
    nominated = set(
        map(tuple, links[links.published_nomination][["entity_id", "source_id"]].to_numpy())
    )
    cta = papers[papers.family == "cta"]
    assert set(map(tuple, cta[["entity_id", "source_id"]].to_numpy())) <= nominated
    candidates = catalog[catalog.is_target_candidate]
    assert candidates.papers_cta.gt(0).sum() == candidates.published_nomination.astype(bool).sum()


def test_missing_values_never_pass_a_threshold(catalog):
    row = catalog[catalog.protein_group_label == "MAGEA4"].copy()
    row["normal_reference_tpm"] = 1.0
    assert tf.apply(tf.TargetFilter(max_off_target_tpm=4.0), row).keep.all()
    row["normal_reference_tpm"] = float("nan")
    d = tf.apply(tf.TargetFilter(max_off_target_tpm=4.0), row)
    assert d.filter_failed.iloc[0] == "off_target"


def test_bad_knobs_are_rejected():
    for bad in (
        {"membership": "everything"},
        {"normal_function": "no_disposable"},
        {"readthrough": "exclude"},
        {"min_patients": 5},
        {"min_tumor_fold": 5.0},
        {"min_reproductive_fold": 5.0},
        {"min_reproductive_fold": float("nan")},
        {"min_tumor_patients": -1},
        {"min_tumor_patients": 1.5},
        {"exclude_prefixes": "MAGE"},
    ):
        with pytest.raises((ValueError, TypeError)):
            tf.TargetFilter(**bad)


def test_describe_states_every_active_knob():
    text = tf.describe("strict")
    for phrase in (
        "below 4 TPM",
        "Vital or Protected",
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
    assert len(removed) == 18
    assert "MAGEC2" in set(removed.protein_group_label)
    assert removed.protein_group_label.str.startswith("MAGE").all()
    assert (removed.filter_failed == "exclusions").all()
    assert d[d.protein_group_label == "MAGEA4"].keep.all()
    assert isinstance(d, pd.DataFrame)


def test_readthrough_loci_are_resolved_by_the_catalog_not_by_a_filter(catalog):
    from oncoref import target_curation

    loci = set(target_curation.readthrough().entity_id)
    rows = catalog[catalog.entity_id.isin(loci)]
    assert len(rows) == len(loci) and rows.readthrough_locus.all()
    assert not rows.is_target_candidate.any() and not rows.sequence_pass.any()
    assert not catalog[~catalog.entity_id.isin(loci)].readthrough_locus.any()
    named = catalog[
        catalog.protein_group_label.isin(["ATF7-NPFF", "ENSG00000273155", "ENSG00000267618"])
    ]
    assert len(named) == 3 and named.classification.isin(["Other:Sequence", "Other:Review"]).all()
    assert not named.sequence_pass.any()  # Review precedence does not erase sequence failure.
    assert rows.readthrough_parents.str.len().gt(0).sum() >= len(rows) - 1


def test_normal_function_knob_follows_the_tier_order(catalog):
    sizes = [
        int(tf.apply(tf.TargetFilter(normal_function=value), catalog).keep.sum())
        for value in tf.NORMAL_FUNCTION
    ]
    assert sizes[0] == 534 and sizes == sorted(sizes, reverse=True) and sizes[-1] < sizes[1]
    strict = tf.apply(tf.TargetFilter(normal_function="no_protected"), catalog)
    kept = strict[strict.keep]
    assert not kept.curated_tier.isin(["vital", "protected"]).any()
    # Hair follicle and cutaneous melanocytes are Disposable in the registry, so
    # hair keratins and PRAME are recorded but never removed by this knob.
    loosest = tf.apply(tf.TargetFilter(normal_function="no_manageable"), catalog)
    disposable = loosest[loosest.is_target_candidate & loosest.curated_tier.eq("disposable")]
    assert len(disposable) and disposable.keep.all()
    assert loosest[loosest.protein_group_label == "PRAME"].curated_tier.item() == "disposable"


def test_filters_act_only_on_records_the_catalog_carries(catalog):
    from oncoref import target_curation as tc
    from oncoref import tissue_tiers

    measured = tf.measures(catalog)
    candidates = measured[measured.is_target_candidate]
    # An excluding entry is Unsafe in the catalog, so no filter ever sees it.
    entries = tc.entries()
    excluded = entries[entries.effect.eq("exclude")]
    assert catalog[catalog.entity_id.isin(excluded.entity_id)].category.eq("Unsafe").all()
    # What layer 3 reads is the catalog's own recorded tier or an acting family rule.
    families = tc.families().set_index("id")
    acting = set(families.index[families.filter_acts])
    for row in candidates[candidates.curated_tier.ne("")].itertuples():
        own = row.literature_recorded_tier if isinstance(row.literature_recorded_tier, str) else ""
        family = row.literature_family_tier if row.literature_family in acting else ""
        measured = (
            row.reference_recorded_tier if isinstance(row.reference_recorded_tier, str) else ""
        )
        assert row.curated_tier == tc.most_protective({own, family, measured} - {""})
        assert row.curated_sources
    assert candidates[candidates.curated_tier.eq("")].literature_recorded_tier.isna().all()
    assert candidates[candidates.curated_tier.eq("")].reference_recorded_tier.isna().all()
    # The nose and skin reference records a tissue only on agreement of people or datasets.
    reference = tc.nose_skin()
    rules = tc.nose_skin_provenance()["rules"]
    recorded = reference[reference.recorded_tissue.ne("")]
    assert set(candidates[candidates.reference_recorded_tier.notna()].entity_id) <= set(
        recorded.entity_id
    )
    for row in recorded.itertuples():
        assert (
            row.olfactory_samples_ge1 >= rules["olfactory_min_samples"]
            or row.nose_datasets_detected >= rules["cell_min_datasets"]
            or row.skin_datasets_detected >= rules["cell_min_datasets"]
        )
    # Every family and recorded tissue resolves through the tissue-tier registry or its updates.
    for rule in families.itertuples():
        assert tissue_tiers.assess(rule.tissue)["tier"] == rule.tier
    assert families[families.filter_acts].sources.str.len().gt(0).all()
    added = tc.compartments()
    assert len(added) == 7 and added.references.map(len).ge(2).all()


def test_host_profile_invariants_on_published_catalog(catalog):
    summary = tf.profile_summary(catalog).set_index("profile")
    assert summary.loc["all", "CTA"] == catalog.is_cta.sum()
    for row in summary.itertuples():
        assert catalog.is_cta.sum() >= row.CTA
        assert row.CTA + row.CTA_unresolved + row.CTA_excluded == catalog.is_cta.sum()
    assert summary.loc["female-patient (pregnant)", "Union"] == 0
    male = tf.apply(tf.TargetFilter(host_profile="male-patient"), catalog)
    female = tf.apply(tf.TargetFilter(host_profile="female-patient"), catalog)
    baseline = tf.apply(tf.TargetFilter(), catalog)
    assert male.keep.equals(baseline.keep) and female.keep.equals(baseline.keep)
    assert tf.protected_compartments(tf.TargetFilter(host_profile="male-patient")) == ()
    protected = tf.apply(
        tf.TargetFilter(host_profile="male-patient", protect_tissues=("testis",)), catalog
    )
    assert not protected[protected.protein_group_label.eq("MAGEA4")].keep.any()
    assert set(protected.loc[protected.keep, "entity_id"]) < set(male.loc[male.keep, "entity_id"])
    # No host policy changes biological category, even when the host rejects a CTA.
    assert male.filtered_category[male.is_cta].eq("CTA").all()
    for bad in (
        {"host_profile": "unknown"},
        {"lactating": "yes"},
        {"allow_tissues": "prostate"},
        {"protect_tissues": ("not_a_compartment",)},
        {"allow_tissues": ("prostate",), "protect_tissues": ("prostate",)},
        {"allow_tissues": ("breast",), "protect_tissues": ("resting_breast",)},
    ):
        with pytest.raises(ValueError):
            tf.TargetFilter(**bad)
