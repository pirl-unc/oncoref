#!/usr/bin/env python3
"""Sum distinct oncoref burden categories represented by each primary CTA's p90 hits."""

import argparse
import json
from pathlib import Path

import pandas as pd
from cta_mortality_coverage_report import report_burden_category
from cta_report_common import seal_stage, verify_analysis, verify_stage

import oncoref as od

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/cta_proteoform_report_20260917"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    out = parser.parse_args().out.resolve()
    verify_analysis(out)
    verify_stage(out, "primary")
    dest = out / "primary_panel"
    primary = pd.read_csv(dest / "primary_proteoforms.csv")
    reference = od.cancer_burden_df().set_index("burden_category")
    reference.to_csv(dest / "burden_reference.csv")
    groups = pd.read_csv(out / "cohort_overlap_groups.csv").set_index("cancer_code")
    cohorts = pd.read_csv(out / "cohort_audit.csv").set_index("cancer_code")
    mapping = {code: report_burden_category(code) for code in cohorts.index}
    metrics = pd.read_csv(out / "all_proteoform_cohort_metrics.csv.gz")
    summary, detail = [], []
    for r in primary.itertuples():
        good = metrics[
            (metrics.proteoform_key == r.proteoform_key)
            & (metrics.transcriptome_percentile == 90)
            & metrics.complete_measurement
            & metrics.n_patients.ge(20)
            & (100 * metrics.n_expressing > 10 * metrics.n_patients)
        ]
        codes = set(good.cancer_code)
        assert codes == set(str(r.p90_cohort_views_gt10).split(";")) - {"nan", ""}
        assert (
            groups.loc[list(codes)].cancer_type_group.nunique() == r.p90_n_cancer_type_groups_gt10
        )
        categories = sorted({mapping[c] for c in codes})
        included = []
        for cat in categories:
            assert cat in reference.index
            ref = reference.loc[cat]
            use = ref.derivation_basis != "residual"
            if use:
                included.append(cat)
            hit_codes = sorted(c for c in codes if mapping[c] == cat)
            detail.append(
                {
                    "proteoform_key": r.proteoform_key,
                    "Symbol": r.Symbol,
                    "burden_category": cat,
                    "included_in_sum": use,
                    "qualifying_cohort_codes": ";".join(hit_codes),
                    "qualifying_cancer_type_groups": ";".join(
                        sorted(set(groups.loc[hit_codes].cancer_type_group))
                    ),
                    "world_incidence_pct": ref.world_incidence_pct,
                    "world_mortality_pct": ref.world_mortality_pct,
                    "source": ref.source,
                    "mapping_scope_note": "ADCC salivary/head-neck scope under audit (#543)"
                    if "ADCC" in hit_codes
                    else "CHOL mixed biliary scope under audit (#543)"
                    if "CHOL" in hit_codes
                    else "",
                    "interpretation": "Full category represented by at least one cohort hit; not patient expression coverage",
                }
            )
        summary.append(
            {
                "proteoform_key": r.proteoform_key,
                "Symbol": r.Symbol,
                "world_incidence_pct_represented": reference.loc[
                    included
                ].world_incidence_pct.sum(),
                "world_mortality_pct_represented": reference.loc[
                    included
                ].world_mortality_pct.sum(),
                "n_burden_categories_represented": len(included),
                "burden_categories_represented": ";".join(included),
            }
        )
    pd.DataFrame(summary).to_csv(dest / "primary_burden_coverage.csv", index=False)
    pd.DataFrame(detail).to_csv(dest / "primary_burden_coverage_details.csv", index=False)
    audit = cohorts.reset_index().merge(
        groups.reset_index(), on="cancer_code", suffixes=("", "_group")
    )
    audit["eligible_n20"] = audit.n_patients.ge(20)
    audit["burden_category"] = audit.cancer_code.map(mapping)
    audit[
        [
            "cancer_code",
            "cancer_name",
            "n_patients",
            "eligible_n20",
            "overlap_group",
            "cancer_type_group",
            "burden_category",
        ]
    ].to_csv(dest / "coverage_denominator_audit.csv", index=False)
    (dest / "burden_coverage_method.json").write_text(
        json.dumps(
            {
                "source": "oncoref.cancer_burden_df(), internal GLOBOCAN2022 reference",
                "mapping": "oncoref.burden_category() applied to actual qualifying cohort codes",
                "definition": "For each protein separately, sum each distinct non-residual burden category once when any eligible fully measured cohort has >10% p90-positive patients.",
                "not_a_running_sum_across_proteins": True,
                "not_patient_level_global_coverage": True,
                "no_renormalization": True,
                "reference_totals_pct": reference[["world_incidence_pct", "world_mortality_pct"]]
                .sum()
                .to_dict(),
                "n_cohort_views": len(cohorts),
                "n_eligible_cohort_views": int(audit.eligible_n20.sum()),
                "n_eligible_cancer_type_groups": int(
                    audit.loc[audit.eligible_n20].cancer_type_group.nunique()
                ),
                "limitations": "Any subtype hit represents the full parent burden category. Pediatric and rare subtypes do not have separate global weights. Rounded internal shares sum approximately to 100%; category-scope audit #543 remains relevant.",
            },
            indent=2,
        )
        + "\n"
    )
    seal_stage(
        out,
        "burden",
        [Path(__file__), out / "analysis_receipt.json", dest / "primary_proteoforms.csv"],
        [
            dest / n
            for n in [
                "primary_burden_coverage.csv",
                "primary_burden_coverage_details.csv",
                "coverage_denominator_audit.csv",
                "burden_coverage_method.json",
            ]
        ],
    )
    print(
        pd.DataFrame(summary)[
            ["Symbol", "world_incidence_pct_represented", "world_mortality_pct_represented"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
