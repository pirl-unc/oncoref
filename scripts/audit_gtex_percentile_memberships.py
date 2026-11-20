"""Compare GTEx donor percentiles against the November 18 protein-level audit."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from audit_gtex_cross_tissue_signals import source_panels
from gtex_sample_evidence import SOURCE, donor_percentile, sample_matrix

from oncoref import reference_data
from oncoref import target_selection as ts

AUDITS = Path(__file__).resolve().parents[1] / "docs/audits"
OUT = AUDITS / "gtex-cross-tissue-signals"


def qc_scenarios(facts):
    """Recompute P99 for source-panel proteins in matched libraries, retaining all other evidence."""
    rebuild = AUDITS / "target-selection-20261110/rebuild"
    normal = pd.read_csv(
        rebuild / "normal-comparison.csv.gz",
        usecols=[
            "entity_id",
            "source_context",
            "sample_id",
            "tissue",
            "source_id",
            "scope",
            "value",
            "qc_pass",
            "identity_resolved",
        ],
    )
    normal = normal.loc[
        normal.qc_pass & normal.identity_resolved & normal.scope.isin(["normal", "disposable"])
    ].copy()
    hpa_qc = pd.read_csv(rebuild / "normal-library-qc.csv")
    flagged = set(pd.read_csv(OUT / "hpa-signals.csv").sample_id)
    flagged.update(
        hpa_qc.loc[hpa_qc[["carryover_P", "carryover_E"]].eq(True).any(axis=1), "sample_id"]
    )
    biopsies = normal.source_id.eq("hpa_tissue_samples:v23") & normal.scope.eq("normal")
    suspect = biopsies & normal.sample_id.isin(flagged)
    # Only perturb a biopsy where the same protein still has another biopsy of this tissue.
    sites = ["entity_id", "tissue"]
    retained = pd.MultiIndex.from_frame(normal.loc[biopsies & ~suspect, sites]).unique()
    omit = suspect & pd.MultiIndex.from_frame(normal[sites]).isin(retained)
    hpa = normal.loc[~omit].groupby("entity_id").value.max().reindex(facts.index)
    native = normal.source_id.eq("gtex_donor_p99:v11")
    sites = pd.read_csv(
        rebuild / "gtex-library-qc.csv.gz", usecols=["source_tissue"]
    ).source_tissue.unique()
    _, panels = source_panels(sorted(sites))
    names = pd.read_csv(
        reference_data.pinned_path("gtex_tissue_rna", "v11"), sep="\t", skiprows=2, usecols=["Name"]
    )
    chosen = pd.read_csv(
        AUDITS / "cta-selection-20261028/alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz"
    ).dropna(subset="row_position")
    versioned = pd.Series(
        names.Name.iloc[chosen.row_position.astype(int)].to_numpy(), index=chosen.gene_id
    )
    groups = facts.member_gene_ids.map(lambda text: frozenset(json.loads(text)))
    matches = pd.read_csv(OUT / "batch-matched.csv")
    for site, hits in matches.groupby("tissue"):
        measured = native & normal.source_context.eq(site)
        if not measured.any():
            continue
        removals = {}
        for hit in hits.itertuples():
            panel = set(versioned.index[versioned.isin(panels[hit.source_panel])])
            for entity, members in groups.items():
                if members and members <= panel:
                    removals.setdefault(entity, set()).add(hit.sample_id)
        if not removals:
            continue
        group = {gene: entity for entity in removals for gene in groups.at[entity]}
        exact = {versioned.at[g].encode(): g for g in group}
        values, donors = sample_matrix(
            reference_data.pinned_path(SOURCE, "v11-" + site.lower()), exact
        )
        grouped = values.groupby(values.index.map(group))
        sizes = groups.loc[list(removals)].map(len)
        sums = grouped.sum(min_count=1).where(grouped.count().eq(sizes, axis=0))
        rows = normal.loc[measured & normal.entity_id.isin(removals)]
        original = donor_percentile(sums, donors, 0.99).reindex(rows.entity_id)
        np.testing.assert_allclose(original, rows.value)
        for entity, samples in removals.items():
            if entity in sums.index:
                sums.loc[entity, list(samples)] = np.nan
        updated = donor_percentile(sums, donors, 0.99).reindex(rows.entity_id).fillna(original)
        normal.loc[rows.index, "value"] = updated.to_numpy()
    return normal.groupby("entity_id").value.max().reindex(facts.index), hpa


def main():
    if ts.POLICY_ID != "oncoref.target-selection.2026-11-18.4":
        raise SystemExit(
            "This historical P99 audit requires the 2026-11-18.4 catalog (commit 8e85f037)."
        )
    ceilings = pd.read_csv(OUT / "protein-sensitivity.csv.gz").set_index("entity_id")
    fields = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "classification",
        "unsafe",
        *ts.SHARED_CHECKS,
        "home_tissue_unmeasured",
        "expression_bounds_status",
        "disposable_cap_pass",
        "target_peak_tpm",
        "tumor_max_tpm",
        "reproductive_supported_tpm",
        "cancer_evidence",
        "reproductive_evidence",
        "normal_reference_tpm",
        "disposable_max_tpm",
    ]
    catalog = pd.read_csv(
        AUDITS / "target-selection-20261110/rebuild/protein-decisions.csv.gz", usecols=fields
    ).set_index("entity_id")
    facts = catalog.loc[ceilings.index]
    if set(catalog.index[catalog.normal_reference_tpm.notna()]) != set(facts.index):
        raise ValueError("Assessed protein universe has changed")
    np.testing.assert_allclose(
        facts[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1), ceilings.ceiling_p99
    )
    summary = []
    changes = facts[["protein_group_label"]].copy()
    matched, hpa = qc_scenarios(facts)
    for key, ceiling in [
        ("p99", ceilings.ceiling_p99),
        ("matched_p99", matched),
        ("hpa_flag_sensitivity_p99", hpa),
        ("p95", ceilings.ceiling_p95),
    ]:
        threshold = ts.ENRICHMENT_FOLD * ceiling.clip(lower=ts.DETECTION_FLOOR)
        d = facts.copy()
        d["expression_bounds_pass"] = (
            d.expression_bounds_status.ne("unassessed")
            & d.disposable_cap_pass
            & d.target_peak_tpm.ge(ts.DETECTION_FLOOR)
            & d.target_peak_tpm.gt(ceiling)
        )
        for field, peak in [
            ("cancer_evidence", d.tumor_max_tpm),
            ("reproductive_evidence", d.reproductive_supported_tpm),
        ]:
            measured = d[field].ne("unassessed")
            d.loc[measured, field] = np.where(
                peak[measured].ge(threshold[measured]), "qualifies", "does_not_qualify"
            )
        results = ts.classify(d)
        if key == "p99" and not results.classification.equals(facts.classification):
            raise ValueError("P99 did not reproduce catalog decisions")
        summary.append(
            {
                "normal_policy": key,
                "assessed": len(d),
                "changed_thresholds": (
                    ~np.isclose(
                        threshold,
                        ts.ENRICHMENT_FOLD * ceilings.ceiling_p99.clip(lower=ts.DETECTION_FLOOR),
                    )
                ).sum(),
                "changed_classifications": results.classification.ne(facts.classification).sum(),
                "tumor_at_tenfold": d.tumor_max_tpm.ge(threshold).sum(),
                "CA": results.is_cancer_antigen.sum(),
                "RA": results.is_reproductive_antigen.sum(),
                "CTA": results.is_cta.sum(),
                "union": results.is_target_candidate.sum(),
            }
        )
        changes[key] = results.classification
        changes[key + "_threshold_tpm"] = threshold
    changes[
        changes.p99.ne(changes.p95)
        | changes.p99.ne(changes.matched_p99)
        | changes.p99.ne(changes.hpa_flag_sensitivity_p99)
    ].sort_values("protein_group_label").to_csv(OUT / "percentile-transitions.csv")
    pd.DataFrame(summary).to_csv(OUT / "percentile-memberships.csv", index=False)


if __name__ == "__main__":
    main()
