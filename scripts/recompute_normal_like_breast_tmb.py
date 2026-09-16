#!/usr/bin/env python3
"""Summarize normal-like TCGA breast TMB from a pinned cBioPortal data release.

Download the three SOURCES separately, then pass their paths. Uses the portal's
nonsynonymous TMB rates without inventing a callable denominator. No network I/O.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, median

REVISION = "0cc9138746c08b304f8dac92c31983e0ef44af1d"
STUDY = "brca_tcga_pan_can_atlas_2018"
SOURCES = {
    "sample": {
        "url": f"https://media.githubusercontent.com/media/cBioPortal/datahub/{REVISION}/public/{STUDY}/data_clinical_sample.txt",
        "sha256": "e96e6e4895a6ddcdb31d508a9e8f5676ca78ce396d2f82cd111a286a17e75ca6",
    },
    "patient": {
        "url": f"https://media.githubusercontent.com/media/cBioPortal/datahub/{REVISION}/public/{STUDY}/data_clinical_patient.txt",
        "sha256": "13c939025794ce21905ccb3409865fa08506b2283b432f8e3a3022f60a4c7682",
    },
    "sequenced": {
        "url": f"https://raw.githubusercontent.com/cBioPortal/datahub/{REVISION}/public/{STUDY}/case_lists/cases_sequenced.txt",
        "sha256": "ae464df8154970f13b1cfe34e5538f9286abbe8bd13e1105474638ff2b45c19a",
    },
}
FIELDS = ("sample_id", "patient_id", "subtype", "nonsynonymous_mut_mb")


def summarize(rows):
    """Unweighted sample statistics; reject duplicates, negative and missing rates."""
    if len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate samples")
    if len({r["patient_id"] for r in rows}) != len(rows):
        raise ValueError("Repeated patients require an explicit sampling policy")
    values = [float(r["nonsynonymous_mut_mb"]) for r in rows]
    if not values or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("Missing or invalid TMB")
    return {
        "n": len(values),
        "median": median(values),
        "mean": mean(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def extract(sample, patient, sequenced):
    paths = {"sample": sample, "patient": patient, "sequenced": sequenced}
    for key, path in paths.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != SOURCES[key]["sha256"]:
            raise ValueError(f"{key} hash differs from the reviewed source")

    def read_tsv(path):
        with Path(path).open() as handle:
            return list(
                csv.DictReader((s for s in handle if not s.startswith("#")), delimiter="\t")
            )

    patients = read_tsv(patient)
    subtypes = {r["PATIENT_ID"]: r["SUBTYPE"] for r in patients}
    if len(subtypes) != len(patients):
        raise ValueError("Duplicate patient metadata")
    metadata = dict(s.split(":", 1) for s in Path(sequenced).read_text().splitlines() if ":" in s)
    eligible = set(metadata["case_list_ids"].split())
    rows = [
        dict(zip(FIELDS, (r["SAMPLE_ID"], r["PATIENT_ID"], "BRCA_Normal", r["TMB_NONSYNONYMOUS"])))
        for r in read_tsv(sample)
        if subtypes.get(r["PATIENT_ID"]) == "BRCA_Normal" and r["SAMPLE_ID"] in eligible
    ]
    rows.sort(key=lambda r: r["sample_id"])
    summarize(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in SOURCES:
        parser.add_argument(f"--{key}", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = extract(args.sample, args.patient, args.sequenced)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "tmb-normal-like-sample-rates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "study": STUDY,
        "revision": REVISION,
        "sources": SOURCES,
        "reference": "PMID:29625048",
        "method": "Join sample PATIENT_ID to patient SUBTYPE=BRCA_Normal; require mutation-profile eligibility in cases_sequenced; summarize portal TMB_NONSYNONYMOUS. All 36 subtype samples pass eligibility and have finite rates; one sample per patient.",
        "limitations": "Derived sample statistic, not a published subgroup summary. Uses portal-provided nonsynonymous TMB, without independent reconstruction of callable bases. Normal-like is an expression classification, susceptible to nonneoplastic admixture; it is not TNBC. This 36-case release differs from the eight normal-like cases excluded in TCGA 2012.",
        "results": summarize(rows),
    }
    (args.output_dir / "tmb-normal-like-recomputed.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
