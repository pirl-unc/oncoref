#!/usr/bin/env python3
"""Build a receptor-defined primary TCGA TNBC subset and a complete membership audit.

Inputs are immutable public source files, checked before selection. Expression
columns are copied from the released BRCA matrix without gene filtering or
renormalization. PAM50 is an annotation for overlap analysis, never a selector.
Run with --download-sources to populate --source-dir from the pinned URLs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REVISION = "0cc9138746c08b304f8dac92c31983e0ef44af1d"
SOURCE_COHORT = "TREEHOUSE_POLYA_25_01_TCGA_BRCA_TNBC"
PARENT = {
    "url": "https://github.com/pirl-unc/oncoref/releases/download/source-v5.22.10/BRCA_per_sample_tpm.parquet",
    "sha256": "e07f737c4f8d068f5aaa4096129ab3be5253b112cdced3969f2514c8afc48fec",
}
_MEDIA = f"https://media.githubusercontent.com/media/cBioPortal/datahub/{REVISION}/public"
_ATLAS = "brca_tcga_pan_can_atlas_2018"
SOURCES = {
    "receptors": {
        "filename": "brca_tcga_patient.tsv",
        "url": f"{_MEDIA}/brca_tcga/data_clinical_patient.txt",
        "sha256": "973b6e03436330e0563a15989d04fb3566a2c7cdc66db9aaa3badbb32c58badc",
    },
    "pam50": {
        "filename": f"{_ATLAS}_patient.tsv",
        "url": f"{_MEDIA}/{_ATLAS}/data_clinical_patient.txt",
        "sha256": "13c939025794ce21905ccb3409865fa08506b2283b432f8e3a3022f60a4c7682",
    },
    "tmb": {
        "filename": f"{_ATLAS}_sample.tsv",
        "url": f"{_MEDIA}/{_ATLAS}/data_clinical_sample.txt",
        "sha256": "e96e6e4895a6ddcdb31d508a9e8f5676ca78ce396d2f82cd111a286a17e75ca6",
    },
    "sequenced": {
        "filename": f"{_ATLAS}_cases_sequenced.txt",
        "url": f"https://raw.githubusercontent.com/cBioPortal/datahub/{REVISION}/public/{_ATLAS}/case_lists/cases_sequenced.txt",
        "sha256": "ae464df8154970f13b1cfe34e5538f9286abbe8bd13e1105474638ff2b45c19a",
    },
}
RECEPTOR_COLUMNS = ("ER_STATUS_BY_IHC", "PR_STATUS_BY_IHC", "IHC_HER2", "HER2_FISH_STATUS")
DEFINITE = {"Positive", "Negative"}


def receptor_call(er: str, pr: str, ihc: str, fish: str) -> tuple[str, str, str]:
    """Return TNBC state, resolved HER2 state, and reason from reported categories.

    A definite FISH call resolves missing/equivocal/indeterminate IHC. Opposing
    definite IHC/FISH calls remain indeterminate. With no definite FISH, use a
    definite IHC call. Any definite positive receptor excludes TNBC; otherwise
    all three receptors must be explicitly negative. No expression, PAM50,
    amplification proxy, or historic percentage bin supplies a missing call.
    """
    conflict = ihc in DEFINITE and fish in DEFINITE and ihc != fish
    her2 = "Indeterminate" if conflict else fish if fish in DEFINITE else ihc
    if her2 not in DEFINITE:
        her2 = "Indeterminate"
    if "Positive" in (er, pr, her2):
        return "non_TNBC", her2, "positive_receptor"
    if er == pr == her2 == "Negative":
        return "TNBC", her2, "three_negative_receptors"
    return "indeterminate", her2, "conflicting_HER2" if conflict else "unresolved_receptor"


def membership(sample_ids, receptors: pd.DataFrame, pam50: pd.DataFrame) -> pd.DataFrame:
    """Audit every expression column; only primary, definitively TNBC cases enter."""
    for name, frame in (("receptor", receptors), ("PAM50", pam50)):
        if frame["PATIENT_ID"].isna().any() or frame["PATIENT_ID"].duplicated().any():
            raise ValueError(f"Missing or duplicate {name} patient metadata")
    receptors = receptors.set_index("PATIENT_ID")
    pam50 = pam50.set_index("PATIENT_ID")["SUBTYPE"]
    sample_ids = list(sample_ids)
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Duplicate expression sample IDs")
    records = []
    for sample in sample_ids:
        if not re.fullmatch(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}", sample):
            raise ValueError(f"Unexpected TCGA sample barcode: {sample}")
        patient = sample.rsplit("-", 1)[0]
        primary = sample.endswith("-01")
        row = receptors.loc[patient] if patient in receptors.index else {}
        raw = {key: row.get(key, "") for key in RECEPTOR_COLUMNS}
        state, her2, reason = receptor_call(*(str(raw[c]) for c in RECEPTOR_COLUMNS))
        if patient not in receptors.index:
            reason = "missing_receptor_metadata"
        records.append(
            {
                "sample_id": sample,
                "patient_id": patient,
                "sample_type": "primary" if primary else "non_primary",
                **raw,
                "HER2_resolved": her2,
                "HER2_conflict": raw["IHC_HER2"] in DEFINITE
                and raw["HER2_FISH_STATUS"] in DEFINITE
                and raw["IHC_HER2"] != raw["HER2_FISH_STATUS"],
                "receptor_state": state,
                "receptor_reason": reason,
                "pam50": pam50.get(patient, "unavailable"),
                "included": primary and state == "TNBC",
                "selection_reason": reason if primary else "non_primary_sample",
            }
        )
    result = pd.DataFrame(records)
    result["pam50"] = result["pam50"].fillna("unavailable")
    if result.loc[result["included"], "patient_id"].duplicated().any():
        raise ValueError("Multiple selected expression samples per patient")
    return result


def tmb_rows(audit: pd.DataFrame, samples: pd.DataFrame, eligible: set[str]) -> pd.DataFrame:
    """TMB for the same primary RNA cohort, with explicit per-sample availability."""
    if samples["SAMPLE_ID"].duplicated().any():
        raise ValueError("Duplicate TMB sample metadata")
    selected = audit.loc[audit["included"], ["sample_id", "patient_id"]]
    joined = selected.merge(
        samples[["SAMPLE_ID", "PATIENT_ID", "TMB_NONSYNONYMOUS"]],
        left_on="sample_id",
        right_on="SAMPLE_ID",
        how="left",
        validate="one_to_one",
    )
    matched = joined["SAMPLE_ID"].notna()
    if (joined.loc[matched, "patient_id"] != joined.loc[matched, "PATIENT_ID"]).any():
        raise ValueError("TMB sample/patient identity mismatch")
    values = pd.to_numeric(joined["TMB_NONSYNONYMOUS"], errors="coerce")
    mutation_eligible = joined["sample_id"].isin(eligible)
    valid = values.notna() & np.isfinite(values) & values.ge(0)
    joined["tmb_included"] = matched & mutation_eligible & valid
    joined["tmb_reason"] = np.select(
        [~matched, ~mutation_eligible, ~valid],
        ["missing_sample_metadata", "not_mutation_profile_eligible", "missing_or_invalid_TMB"],
        default="included",
    )
    joined["nonsynonymous_mut_mb"] = values.where(joined["tmb_included"])
    return joined[["sample_id", "patient_id", "nonsynonymous_mut_mb", "tmb_included", "tmb_reason"]]


def verified_bytes(path: Path, expected: str) -> bytes:
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"SHA-256 differs from the reviewed source: {path}")
    return data


def build(parent: Path, source_dir: Path, output_dir: Path) -> dict:
    verified_bytes(parent, PARENT["sha256"])
    for source in SOURCES.values():
        verified_bytes(source_dir / source["filename"], source["sha256"])
    frames = {
        key: pd.read_csv(source_dir / SOURCES[key]["filename"], sep="\t", comment="#")
        for key in ("receptors", "pam50", "tmb")
    }
    matrix = pd.read_parquet(parent)
    audit = membership(
        matrix.columns.drop(["Ensembl_Gene_ID", "Symbol"]), frames["receptors"], frames["pam50"]
    )
    ids = audit.loc[audit["included"], "sample_id"].tolist()
    if not ids:
        raise ValueError("No definitive primary TNBC samples")
    subset = matrix[["Ensembl_Gene_ID", "Symbol", *ids]]
    derived = output_dir / "derived"
    derived.mkdir(parents=True, exist_ok=True)
    matrix_path = derived / "tcga_brca_tnbc_per_sample_tpm.parquet"
    subset.to_parquet(matrix_path, index=False)
    audit.to_csv(output_dir / "brca-receptor-membership.csv", index=False)
    primary = audit.loc[audit["sample_type"].eq("primary")]
    overlap = pd.crosstab(primary["pam50"], primary["receptor_state"]).reindex(
        columns=["TNBC", "non_TNBC", "indeterminate"], fill_value=0
    )
    overlap.to_csv(output_dir / "brca-tnbc-pam50-overlap.csv")
    case_text = (source_dir / SOURCES["sequenced"]["filename"]).read_text()
    cases = dict(line.split(":", 1) for line in case_text.splitlines() if ":" in line)
    rates = tmb_rows(audit, frames["tmb"], set(cases["case_list_ids"].split()))
    rates.to_csv(output_dir / "brca-tnbc-tmb.csv", index=False)
    measured = rates.loc[rates["tmb_included"], "nonsynonymous_mut_mb"]
    if measured.empty:
        raise ValueError("No mutation-profile-eligible TNBC TMB measurements")
    report = {
        "source_cohort": SOURCE_COHORT,
        "revision": REVISION,
        "parent": PARENT,
        "sources": SOURCES,
        "matrix_sha256": hashlib.sha256(matrix_path.read_bytes()).hexdigest(),
        "n_genes": len(subset),
        "n_samples": len(ids),
        "n_parent_samples": len(audit),
        "n_primary_samples": len(primary),
        "primary_receptor_states": primary["receptor_state"].value_counts().to_dict(),
        "tnbc_pam50": audit.loc[audit["included"], "pam50"].value_counts().to_dict(),
        "tmb": {
            "n": len(measured),
            "median": measured.median(),
            "mean": measured.mean(),
            "minimum": measured.min(),
            "maximum": measured.max(),
        },
        "limitations": "Historical categorical receptor calls, not modern threshold re-adjudication. "
        "Primary TCGA tissue is not the metastatic trial population. TMB is a portal-rate "
        "reanalysis of the same receptor-defined RNA samples with mutation profiles, not a "
        "published subgroup statistic or a reconstruction of callable bases. PAM50 is independent.",
    }
    (output_dir / "brca-tnbc-build.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--brca-matrix", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--download-sources", action="store_true")
    args = parser.parse_args()
    if args.download_sources:
        args.source_dir.mkdir(parents=True, exist_ok=True)
        for source in SOURCES.values():
            path = args.source_dir / source["filename"]
            if not path.exists():
                data = urllib.request.urlopen(source["url"], timeout=60).read()
                if hashlib.sha256(data).hexdigest() != source["sha256"]:
                    raise ValueError(f"Download checksum mismatch: {source['url']}")
                path.write_bytes(data)
    print(json.dumps(build(args.brca_matrix, args.source_dir, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
