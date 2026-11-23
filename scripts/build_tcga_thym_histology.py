#!/usr/bin/env python3
"""Build WHO-histology-defined TCGA thymoma (THYM) and thymic carcinoma (THYMCA) cohorts.

Inputs are immutable public source files, verified by SHA-256 before selection.
Expression columns are copied from the released THYM matrix without gene filtering
or renormalization. Only primary tumor barcodes (-01) with concordant WHO histology
and ICD-O-3 morphology enter THYM (WHO Types A/AB/B1/B2/B3) or THYMCA (WHO Type C,
ICD-O-3 8586/3). Run with --download-sources to populate --source-dir from pinned URLs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

REVISION = "0cc9138746c08b304f8dac92c31983e0ef44af1d"
SOURCE_COHORT = "TREEHOUSE_POLYA_25_01_TCGA_THYM_HISTOLOGY"
PARENT = {
    "url": "https://github.com/pirl-unc/oncoref/releases/download/source-v5.22.10/THYM_per_sample_tpm.parquet",
    "sha256": "14d8eebcd8e1207b8fd6f54f3728f19d1f6c8519bd027cd23c71d51b543c5707",
    "md5": "063f6729be3ca25f45dc4d09fa5d8603",
}
_MEDIA = f"https://media.githubusercontent.com/media/cBioPortal/datahub/{REVISION}/public"
_ATLAS = "thym_tcga_pan_can_atlas_2018"
SOURCES = {
    "pan_sample": {
        "filename": f"{_ATLAS}_sample.tsv",
        "url": f"{_MEDIA}/{_ATLAS}/data_clinical_sample.txt",
        "sha256": "924a652da15705672d4dba11f7fb9a040d39fd42d63cfcce6cdc4008d6c20a71",
    },
    "pan_patient": {
        "filename": f"{_ATLAS}_patient.tsv",
        "url": f"{_MEDIA}/{_ATLAS}/data_clinical_patient.txt",
        "sha256": "6d2d8a0bd6ca5c8795475677eda5e4a42db1cbe29ad1d1749dc2abd9ff3438e4",
    },
    "tcga_patient": {
        "filename": "thym_tcga_patient.tsv",
        "url": f"{_MEDIA}/thym_tcga/data_clinical_patient.txt",
        "sha256": "474a4cceb0a3f1d4fe76290a435d737a54cfa71511a7b32057d94ef3c239dd66",
    },
}

THYMCA_WHO_LABELS = frozenset({"Thymoma, Type C", "Thymoma; Type C"})
THYMOMA_WHO_LABELS = frozenset(
    {
        "Thymoma, Type A",
        "Thymoma; Type A",
        "Thymoma, Type Ab",
        "Thymoma; Type AB",
        "Thymoma, Type A|Thymoma, Type Ab",
        "Thymoma; Type A|Thymoma; Type AB",
        "Thymoma, Type B1",
        "Thymoma; Type B1",
        "Thymoma, Type B1|Thymoma, Type B2",
        "Thymoma; Type B1|Thymoma; Type B2",
        "Thymoma, Type B2",
        "Thymoma; Type B2",
        "Thymoma, Type B2|Thymoma, Type B3",
        "Thymoma; Type B2|Thymoma; Type B3",
        "Thymoma, Type B3",
        "Thymoma; Type B3",
    }
)
THYMCA_ICDO_CODES = frozenset({"8586/3"})
THYMOMA_ICDO_CODES = frozenset(
    {
        "8581/1",
        "8581/3",
        "8582/1",
        "8582/3",
        "8583/1",
        "8583/3",
        "8584/1",
        "8584/3",
        "8585/3",
        "8521/1",
    }
)


def _clean_field(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if text in {"", "nan", "[Not Available]", "[Unknown]", "NA"}:
        return ""
    return text


def histology_call(
    pan_tumor_type: object,
    tcga_hist_diag: object,
    pan_icdo: object,
    tcga_icdo: object,
) -> tuple[str, str, str]:
    """Return ('THYM' | 'THYMCA' | 'indeterminate', normalized_who_subtype, reason).

    WHO Type C / ICD-O-3 8586/3 defines thymic carcinoma (THYMCA). WHO Types
    A, AB, B1, B2, and B3 (ICD-O-3 8581-8585 and 8521/1) define thymoma (THYM).
    Any missing metadata across all four fields, unrecognized histology label,
    or disagreement between carcinoma and thymoma fields yields 'indeterminate'.
    """
    pt = _clean_field(pan_tumor_type)
    hd = _clean_field(tcga_hist_diag)
    pi = _clean_field(pan_icdo)
    ti = _clean_field(tcga_icdo)

    if not any((pt, hd, pi, ti)):
        return "indeterminate", "unavailable", "missing_histology_metadata"

    calls: list[str] = []
    for label in (pt, hd):
        if not label:
            continue
        if label in THYMCA_WHO_LABELS:
            calls.append("THYMCA")
        elif label in THYMOMA_WHO_LABELS:
            calls.append("THYM")
        else:
            return "indeterminate", label, "unrecognized_who_histology"

    for icdo in (pi, ti):
        if not icdo:
            continue
        if icdo in THYMCA_ICDO_CODES:
            calls.append("THYMCA")
        elif icdo in THYMOMA_ICDO_CODES:
            calls.append("THYM")
        else:
            return "indeterminate", icdo, "unrecognized_icdo_morphology"

    unique_calls = set(calls)
    if len(unique_calls) != 1:
        return "indeterminate", pt or hd or pi or ti, "conflicting_histology_fields"

    resolved = calls[0]
    subtype = (
        (hd or pt).replace("Thymoma; ", "").replace("Thymoma, ", "").replace("Type Ab", "Type AB")
    )
    if resolved == "THYMCA":
        return "THYMCA", "WHO Type C (Thymic Carcinoma, 8586/3)", "concordant_thymic_carcinoma"
    return "THYM", f"WHO {subtype}", "concordant_thymoma"


def membership(
    sample_ids,
    pan_sample: pd.DataFrame,
    pan_patient: pd.DataFrame,
    tcga_patient: pd.DataFrame,
) -> pd.DataFrame:
    """Audit every parent expression column and route primary THYM vs THYMCA."""
    if pan_sample["SAMPLE_ID"].isna().any() or pan_sample["SAMPLE_ID"].duplicated().any():
        raise ValueError("Missing or duplicate PanCancer sample metadata")
    for name, frame in (("PanCancer patient", pan_patient), ("TCGA patient", tcga_patient)):
        if frame["PATIENT_ID"].isna().any() or frame["PATIENT_ID"].duplicated().any():
            raise ValueError(f"Missing or duplicate {name} metadata")

    pan_s = pan_sample.set_index("SAMPLE_ID")
    pan_p = pan_patient.set_index("PATIENT_ID")
    tcga_p = tcga_patient.set_index("PATIENT_ID")

    sample_ids = list(sample_ids)
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Duplicate expression sample IDs")

    records = []
    for sample in sample_ids:
        if not re.fullmatch(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}", sample):
            raise ValueError(f"Unexpected TCGA sample barcode: {sample}")
        patient = sample.rsplit("-", 1)[0]
        primary = sample.endswith("-01")
        pan_tumor_type = pan_s.loc[sample, "TUMOR_TYPE"] if sample in pan_s.index else ""
        tcga_hist_diag = (
            tcga_p.loc[patient, "HISTOLOGICAL_DIAGNOSIS"] if patient in tcga_p.index else ""
        )
        pan_icdo = pan_p.loc[patient, "ICD_O_3_HISTOLOGY"] if patient in pan_p.index else ""
        tcga_icdo = tcga_p.loc[patient, "ICD_O_3_HISTOLOGY"] if patient in tcga_p.index else ""
        code, who_subtype, reason = histology_call(
            pan_tumor_type, tcga_hist_diag, pan_icdo, tcga_icdo
        )
        included = primary and code in {"THYM", "THYMCA"}
        records.append(
            {
                "sample_id": sample,
                "patient_id": patient,
                "sample_type": "primary" if primary else "non_primary",
                "pan_tumor_type": _clean_field(pan_tumor_type),
                "tcga_histological_diagnosis": _clean_field(tcga_hist_diag),
                "pan_icd_o_3_histology": _clean_field(pan_icdo),
                "tcga_icd_o_3_histology": _clean_field(tcga_icdo),
                "histology_state": code,
                "who_subtype": who_subtype,
                "histology_reason": reason,
                "routed_cancer_code": code if included else "",
                "included": included,
                "selection_reason": reason if primary else "non_primary_sample",
            }
        )
    result = pd.DataFrame(records)
    if result.loc[result["included"], "patient_id"].duplicated().any():
        raise ValueError("Multiple selected expression samples per patient")
    return result


def sample_manifest_rows(audit: pd.DataFrame, parent_path: Path) -> pd.DataFrame:
    """Emit the 20-column public sample-manifest rows for THYM and THYMCA."""
    parent_bytes = parent_path.stat().st_size
    parent_md5 = hashlib.md5(parent_path.read_bytes()).hexdigest()
    pipeline = (
        "treehouse_polya_25_01_tcga_thym_histology_tpm_to_tpm_oncoref_canonical_clean_tpm_16_9_75"
    )
    rows = []
    for row in audit.itertuples(index=False):
        code = str(row.routed_cancer_code) if row.included else "THYM_EPITHELIAL"
        rows.append(
            {
                "cancer_code": code,
                "source_cohort": SOURCE_COHORT,
                "source_project": "Treehouse (TCGA-THYM) x cBioPortal WHO histology",
                "case_id": row.patient_id,
                "sample_id": row.sample_id,
                "source_file_id": row.sample_id,
                "source_file_name": "THYM_per_sample_tpm.parquet",
                "source_project_id": "TCGA-THYM",
                "sample_type": "Primary Tumor" if row.sample_type == "primary" else "Non-primary",
                "primary_diagnosis": row.who_subtype,
                "md5sum": parent_md5,
                "file_size": parent_bytes,
                "workflow_type": "Treehouse 25.01 PolyA RSEM TPM split by TCGA WHO histology",
                "raw_unit": "TPM",
                "processing_pipeline": pipeline,
                "source_url": PARENT["url"],
                "lineage_evidence_source": (
                    f"cBioPortal {REVISION}: {row.tcga_histological_diagnosis or row.pan_tumor_type}; "
                    f"ICD-O-3={row.tcga_icd_o_3_histology or row.pan_icd_o_3_histology}"
                ),
                "included": bool(row.included),
                "exclusion_reason": "" if row.included else str(row.selection_reason),
                "lineage_label": code if row.included else "",
            }
        )
    return pd.DataFrame(rows)


def verified_bytes(path: Path, expected: str) -> bytes:
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"SHA-256 differs from the reviewed source: {path}")
    return data


def build(parent: Path, source_dir: Path, output_dir: Path, audit_dir: Path | None = None) -> dict:
    verified_bytes(parent, PARENT["sha256"])
    for source in SOURCES.values():
        verified_bytes(source_dir / source["filename"], source["sha256"])
    frames = {
        key: pd.read_csv(source_dir / SOURCES[key]["filename"], sep="\t", comment="#")
        for key in ("pan_sample", "pan_patient", "tcga_patient")
    }
    matrix = pd.read_parquet(parent)
    audit = membership(
        matrix.columns.drop(["Ensembl_Gene_ID", "Symbol"]),
        frames["pan_sample"],
        frames["pan_patient"],
        frames["tcga_patient"],
    )
    thym_ids = audit.loc[audit["routed_cancer_code"].eq("THYM"), "sample_id"].tolist()
    thymca_ids = audit.loc[audit["routed_cancer_code"].eq("THYMCA"), "sample_id"].tolist()
    if not thym_ids or not thymca_ids:
        raise ValueError("Missing primary THYM or THYMCA samples after histology routing")

    derived = output_dir / "derived"
    derived.mkdir(parents=True, exist_ok=True)
    thym_path = derived / "THYM_per_sample_tpm.parquet"
    thymca_path = derived / "THYMCA_per_sample_tpm.parquet"
    matrix[["Ensembl_Gene_ID", "Symbol", *thym_ids]].to_parquet(thym_path, index=False)
    matrix[["Ensembl_Gene_ID", "Symbol", *thymca_ids]].to_parquet(thymca_path, index=False)

    manifest_df = sample_manifest_rows(audit, parent)
    manifest_path = derived / "thym_histology_sample_manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)

    target_audit_dir = audit_dir if audit_dir is not None else output_dir
    target_audit_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(target_audit_dir / "thym-histology-membership.csv", index=False)

    report = {
        "source_cohort": SOURCE_COHORT,
        "revision": REVISION,
        "parent": PARENT,
        "sources": SOURCES,
        "matrix_sha256": {
            "THYM": hashlib.sha256(thym_path.read_bytes()).hexdigest(),
            "THYMCA": hashlib.sha256(thymca_path.read_bytes()).hexdigest(),
        },
        "n_genes": len(matrix),
        "n_parent_samples": len(audit),
        "n_samples_by_code": {
            "THYM": len(thym_ids),
            "THYMCA": len(thymca_ids),
        },
        "n_indeterminate": int(audit["histology_state"].eq("indeterminate").sum()),
        "who_subtype_counts": audit.loc[audit["included"], "who_subtype"].value_counts().to_dict(),
        "limitations": (
            "TCGA-THYM splits cleanly into 108 primary WHO Type A/AB/B1/B2/B3 thymomas (THYM) "
            "and 11 primary WHO Type C / ICD-O-3 8586/3 thymic carcinomas (THYMCA) with zero "
            "cross-table conflicts. THYMCA n=11 has thin upper-tail quantile precision; "
            "THYM_EPITHELIAL unions both disjoint cohorts (n=119)."
        ),
    }
    (target_audit_dir / "thym-histology-build.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--thym-matrix", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, default=None)
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
    print(
        json.dumps(
            build(args.thym_matrix, args.source_dir, args.output_dir, args.audit_dir), indent=2
        )
    )


if __name__ == "__main__":
    main()
