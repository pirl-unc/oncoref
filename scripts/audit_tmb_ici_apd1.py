#!/usr/bin/env python3
"""Inventory every TMB/ICI/aPD1 number without mistaking consistency for validation.

This offline audit detects arithmetic, endpoint, regimen and provenance risks.
It does not establish that a paper supports a value: source reviews and unresolved
items remain explicit in the output. Run from any directory with --output PATH.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLES = {
    "cancer-tmb": "median_tmb_mut_mb",
    "cancer-ici-response": "orr_pct",
    "cancer-apd1-response": "apd1_orr_pct",
    "cancer-ici-response-estimates": "value",
}
FIELDS = (
    "dataset",
    "row_key",
    "cancer_code",
    "regimen",
    "metric",
    "value",
    "reference",
    "source_review_status",
    "source_locator",
    "findings",
)
CONTEXT_BASES = {
    "reported_context",
    "derived_blend",
    "derived_cross_cohort",
    "inferred_from_outcomes",
}


def read_table(data_dir, name):
    with (Path(data_dir) / f"{name}.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def regimen_findings(row):
    """Check explicit drug identities; never infer a regimen from response size."""
    regimen = row.get("regimen", row.get("drug_target", ""))
    drug = row.get("drug", "").lower()
    findings = []
    if regimen in {"PD-1", "PD-L1"}:
        if any(
            x in drug
            for x in (
                "ipilimumab",
                "tremelimumab",
                "bevacizumab",
                "tiragolumab",
                "chemotherapy",
                "everolimus",
                "sunitinib",
                "sorafenib",
            )
        ):
            findings.append("drug_regimen_mismatch")
        if "+/-" in drug or "pooled anti-pd-1/pd-l1" in drug:
            findings.append("mixed_regimen")
    if regimen == "PD-1+CTLA-4" and not ("nivolumab" in drug and "ipilimumab" in drug):
        findings.append("drug_regimen_mismatch")
    return findings


def audit(data_dir=ROOT / "oncoref" / "data"):
    tables = {name: read_table(data_dir, name) for name in TABLES}
    estimates = tables["cancer-ici-response-estimates"]
    primary = {
        (r["cancer_code"], r["regimen"]): r
        for r in estimates
        if r["role"] == "primary" and r["metric"] == "ORR"
    }
    tmb_audit_path = Path(data_dir) / "cancer-tmb-source-audit.csv"
    tmb_reviews = (
        {r["cancer_code"]: r for r in read_table(data_dir, "cancer-tmb-source-audit")}
        if tmb_audit_path.exists()
        else {}
    )
    records = []
    for dataset, rows in tables.items():
        for row in rows:
            value = number(row[TABLES[dataset]])
            regimen = row.get("regimen", row.get("drug_target", ""))
            key = row.get("estimate_id", row["cancer_code"] + (f"/{regimen}" if regimen else ""))
            findings = []
            source = row
            metric = row.get("metric", "TMB" if dataset == "cancer-tmb" else "ORR")
            ref = row.get("ref", row.get("pmid_doi", ""))
            if row[TABLES[dataset]] and value is None:
                findings.append("invalid_numeric_value")
            if value is not None:
                if value < 0 and metric != "TUMOR_SHRINKAGE":
                    findings.append("negative_value")
                if metric in {"ORR", "CRR", "PR", "PRR", "DCR", "CBR"} and value > 100:
                    findings.append("percentage_out_of_bounds")
                if not ref:
                    findings.append("no_direct_reference")
            if dataset in {"cancer-ici-response", "cancer-apd1-response"} and value is not None:
                source = primary.get((row["cancer_code"], regimen), {})
                if not source:
                    findings.append("missing_primary_endpoint")
                else:
                    if ref != source["ref"]:
                        findings.append("anchor_reference_mismatch")
                    if (
                        number(source["value"]) is not None
                        and abs(value - number(source["value"])) > 0.05
                    ):
                        findings.append("anchor_value_differs_from_endpoint")
                    if source["value_basis"] in CONTEXT_BASES:
                        findings.append("nonreported_population_anchor")
                if re.search(
                    r"case report|single.patient|exceptional respon",
                    row.get("setting", "") + row.get("notes", ""),
                    re.I,
                ):
                    findings.append("selected_case_as_population_anchor")
            if dataset == "cancer-tmb":
                review = tmb_reviews.get(row["cancer_code"], {})
                status = review.get(
                    "source_review_status",
                    "needs_source_review" if value is not None else "audited_gap",
                )
                locator = review.get("source_locator", "")
                if value is not None:
                    if status != "source_checked":
                        findings.append("tmb_source_not_revalidated")
                    notes = row["notes"].lower()
                    if re.search(
                        r"\bmean\b|inferred|order.of.magnitude|approximate|no published per.mb median",
                        notes,
                    ):
                        findings.append("tmb_statistic_or_derivation_requires_review")
                    if not row["n_samples"]:
                        findings.append("sample_size_not_curated")
                    if not review.get("tmb_assay"):
                        findings.append("assay_not_structured")
            else:
                status = (
                    "prior_locator_" + source.get("source_locator_status", "missing")
                    if value is not None
                    else "audited_gap"
                )
                locator = source.get("source_locator", "")
                findings += regimen_findings(row)
                if value is not None and source.get("source_locator_status") != "verified":
                    findings.append("numeric_source_location_not_verified")
                n, k = number(source.get("metric_n")), number(source.get("responders"))
                lo, hi = number(source.get("ci_low")), number(source.get("ci_high"))
                if n is not None and (n <= 0 or n != int(n)):
                    findings.append("invalid_denominator")
                if n is not None and 0 < n < 10:
                    findings.append("small_denominator")
                if k is not None and (k < 0 or k != int(k) or (n is not None and k > n)):
                    findings.append("invalid_numerator")
                if lo is not None and hi is not None and lo > hi:
                    findings.append("reversed_ci")
                source_value = number(source.get("value"))
                if source_value is not None and (
                    (lo is not None and source_value < lo - 0.1)
                    or (hi is not None and source_value > hi + 0.1)
                ):
                    findings.append("value_outside_ci")
                if (
                    source.get("unit") == "percent"
                    and n
                    and k is not None
                    and source_value is not None
                    and abs(source_value - 100 * k / n) > 0.6
                ):
                    findings.append("count_rate_mismatch")
                if metric in {"ORR", "CRR", "PR", "PRR", "DCR"} and value is not None:
                    if n is None:
                        findings.append("response_denominator_missing")
                    if source.get("value_basis") in CONTEXT_BASES:
                        findings.append("excluded_from_pooling")
            records.append(
                dict(
                    zip(
                        FIELDS,
                        [
                            dataset,
                            key,
                            row["cancer_code"],
                            regimen,
                            metric,
                            row[TABLES[dataset]],
                            ref,
                            status,
                            locator,
                            ";".join(sorted(set(findings))),
                        ],
                    )
                )
            )
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "oncoref" / "data")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = audit(args.data_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Inventoried {len(rows)} rows; consistency checks do not establish source validity.")
    for flag, count in sorted(
        Counter(flag for row in rows for flag in row["findings"].split(";") if flag).items()
    ):
        print(f"{flag}: {count}")


if __name__ == "__main__":
    main()
