# Licensed under the Apache License, Version 2.0.
"""Per-sample expression QC, housekeeping coverage, and artifact build manifests (#756)."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from . import data_bundle
from .expression_engine import id_columns, sample_columns
from .gene_ids import unversioned
from .normalization import clean_tpm
from .version import DATA_VERSION, SOURCE_MATRIX_VERSION


def _expr():
    from . import expression

    return expression


_PER_SAMPLE_NORMALIZE = ("tpm_raw", "tpm_clean", "tpm_clean_log1p", "tpm_clean_hk")
_SAMPLE_QC_MODES = ("all", "pass", "pass_or_warn")
_ARTIFACT_SAMPLE_QC_MODES = ("artifact", *_SAMPLE_QC_MODES)
SAMPLE_EXPRESSION_QC_POLICY_VERSION = "sample_expression_qc_v5"
EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION = "expression_artifact_build_metadata_v3"
SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH = "source-matrix-sample-qc.csv"
EXPRESSION_ARTIFACT_BUILD_METADATA_PATH = "expression-artifact-build-metadata.csv"
EXPRESSION_ARTIFACT_BUILD_METADATA_JSON_PATH = "expression-artifact-build-metadata.json"

_SOURCE_MATRIX_SAMPLE_QC_MANIFEST_COLUMNS = [
    "cancer_code",
    "source_cohort",
    "source_type",
    "unit",
    "source_scale_class",
    "linear_tpm_comparable",
    "sample_id",
    "n_measured_genes",
    "n_detected_genes",
    "n_detected_raw",
    "n_detected_clean",
    "n_detected_clean_biological",
    "detected_gene_fraction",
    "zero_fraction_raw",
    "parse_missing_fraction",
    "total_tpm",
    "top_gene_id",
    "top_gene_symbol",
    "top_gene_tpm",
    "top_gene_fraction",
    "top10_tpm",
    "top10_fraction",
    "top1_fraction_clean",
    "top10_fraction_clean",
    "housekeeping_genes_present",
    "housekeeping_genes_detected",
    "housekeeping_detection_floor_tpm",
    "housekeeping_genes_above_floor",
    "housekeeping_fraction_above_floor",
    "housekeeping_zero_fraction",
    "tpm_proxy",
    "sample_qc_status",
    "sample_qc_reasons",
    "sample_qc_notes",
    "passes_expression_qc",
    "recommended_for_absolute_tpm_floor",
    "source_matrix_path",
]

_EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS = [
    "cancer_code",
    "source_cohort",
    "build_source_cohort",
    "source_version",
    "source_matrix_path",
    "sample_qc",
    "sample_qc_policy_version",
    "n_source_samples",
    "n_cohort_samples",
    "n_qc_pass",
    "n_qc_warn",
    "n_qc_fail",
    "partition_policy_version",
    "partition_validation_target",
    "partition_status",
    "n_partition_train",
    "n_partition_validation",
    "n_partition_validation_external",
    "n_partition_audit_only",
]
DEFAULT_MIN_DETECTED_GENES_FOR_QC = 5000
DEFAULT_MIN_HOUSEKEEPING_GENES_FOR_QC = 10
DEFAULT_MAX_ZERO_FRACTION_FOR_QC = 0.70
DEFAULT_MAX_TOP_GENE_FRACTION_FOR_QC = 0.20
DEFAULT_MAX_TOP10_GENE_FRACTION_FOR_QC = 0.50
DEFAULT_HOUSEKEEPING_DETECTION_FLOOR_TPM = 30.0
_BLOCKING_SAMPLE_QC_REASONS = frozenset(
    {
        "low_detected_genes",
        "low_housekeeping_detection",
        "low_housekeeping_floor_fraction",
        "high_zero_fraction",
        "high_top_gene_fraction",
        "high_top10_gene_fraction",
        "broad_background",
        "several_unrelated_tissues",
    }
)


def _clean_tpm_housekeeping_panel_ids() -> frozenset[str]:
    from .gene_families import clean_tpm_biological_housekeeping_gene_ids

    return clean_tpm_biological_housekeeping_gene_ids()


def _min_housekeeping_detected(panel_ids) -> int | None:
    n = len(panel_ids)
    return min(DEFAULT_MIN_HOUSEKEEPING_GENES_FOR_QC, n) if n else None


def cancer_reference_expression_source_metadata(
    cancer_type: str, *, source_cohort: str | None = None
) -> dict[str, str | bool | int | float | None]:
    """Structured provenance and scale metadata for one expression source.

    ``cancer_type`` accepts canonical codes and aliases. With ``source_cohort``
    omitted, the source selected by the source-matrix registry is returned; an
    explicit cohort instead resolves that physical source without falling back to
    another source registered for the same cancer type. The result includes the
    stable source identity, PMID, native unit, scale class, linear-TPM
    comparability, sample count when selected, and build provenance.

    Unknown build-time placeholder codes are retained as given and return
    ``source_scale_class="unknown"`` rather than borrowing unrelated metadata.
    """
    try:
        code = _expr().resolve_cancer_type(cancer_type) or str(cancer_type)
    except ValueError:
        code = str(cancer_type)
    source_type = unit = None
    n_reference_samples = None
    source_was_explicit = source_cohort is not None
    if source_cohort is None:
        try:
            info = _expr().source_matrices.cohort_info(code)
            source_cohort = str(info.get("source_cohort") or "")
            n_reference_samples = info.get("n_samples")
        except _expr().source_matrices.SourceMatrixError:
            source_cohort = ""
    else:
        source_cohort = str(source_cohort)

    from .expression_registry import sources_for_cancer_code

    sources = sources_for_cancer_code(code)
    selected = next(
        (s for s in sources if source_cohort and s.source_cohort == source_cohort),
        None if source_was_explicit or not sources else sources[0],
    )
    if selected is not None:
        source_type = selected.source_type
        unit = selected.unit
        if not source_cohort:
            source_cohort = selected.source_cohort or ""
        special = selected.special_handling or ""
    else:
        special = ""
    cohort_metadata = _expr()._source_cohort_metadata_map().get(source_cohort, {})
    source_type = source_type or cohort_metadata.get("kind")
    source_project = (
        getattr(selected, "source_project", None)
        or cohort_metadata.get("source_project")
        or getattr(selected, "project_id", None)
        or getattr(selected, "accession", None)
    )
    processing_pipeline = getattr(selected, "processing_pipeline", None) or cohort_metadata.get(
        "provenance"
    )
    scale_evidence = (
        source_type,
        unit,
        special,
        cohort_metadata.get("assay"),
        processing_pipeline,
    )
    text = " ".join(str(value or "") for value in scale_evidence).lower()
    inferred_microarray_proxy = "microarray" in text or "tpm-proxy" in text or "tpm proxy" in text
    inferred_scrna_proxy = "scrna" in text or "single-cell" in text
    explicit_scale_class = getattr(selected, "source_scale_class", None)
    explicit_comparability = getattr(selected, "linear_tpm_comparable", None)
    explicit_tpm_proxy = getattr(selected, "tpm_proxy", None)
    tpm_proxy = bool(explicit_tpm_proxy) or inferred_microarray_proxy or inferred_scrna_proxy
    if explicit_scale_class is not None:
        source_scale_class = explicit_scale_class
        linear_tpm_comparable = (
            bool(explicit_comparability) if explicit_comparability is not None else False
        )
    elif inferred_microarray_proxy:
        source_scale_class = "microarray_tpm_proxy"
        linear_tpm_comparable = False
        unit = unit or "TPM proxy"
    elif inferred_scrna_proxy:
        source_scale_class = "scrna_pseudobulk_ntpm"
        linear_tpm_comparable = False
    elif selected is not None or "rna-seq" in text or "rnaseq" in text:
        source_scale_class = "linear_rnaseq_tpm"
        linear_tpm_comparable = True
        unit = unit or "TPM"
    else:
        source_scale_class = "unknown"
        linear_tpm_comparable = False
    return {
        "source_cohort": source_cohort or None,
        "source_project": source_project,
        "source_version": getattr(selected, "source_version", None),
        "source_pmid": getattr(selected, "source_pmid", None),
        "source_type": source_type,
        "unit": unit,
        "source_scale_class": source_scale_class,
        "linear_tpm_comparable": linear_tpm_comparable,
        "tumor_origin": getattr(selected, "tumor_origin", None),
        "metastasis_site": getattr(selected, "metastasis_site", None),
        "n_reference_genes": None,
        "n_reference_samples": n_reference_samples,
        "processing_pipeline": processing_pipeline,
        "notes": getattr(selected, "notes", None) or special or None,
        "tpm_proxy": tpm_proxy,
    }


def _sample_qc_status(reasons: list[str]) -> str:
    if set(reasons) & _BLOCKING_SAMPLE_QC_REASONS:
        return "fail"
    return "warn" if reasons else "pass"


def _validate_sample_qc(sample_qc: str) -> str:
    mode = str(sample_qc).lower()
    if mode not in _SAMPLE_QC_MODES:
        raise ValueError(f"sample_qc must be one of {_SAMPLE_QC_MODES}")
    return mode


def _validate_artifact_sample_qc(sample_qc: str) -> str:
    mode = str(sample_qc).lower()
    if mode not in _ARTIFACT_SAMPLE_QC_MODES:
        raise ValueError(f"sample_qc must be one of {_ARTIFACT_SAMPLE_QC_MODES}")
    return mode


def _apply_sample_qc_filter(
    df: pd.DataFrame, code: str, *, sample_qc: str, auto_fetch: bool
) -> pd.DataFrame:
    mode = _validate_sample_qc(sample_qc)
    if mode == "all":
        return df
    samples = sample_columns(df)
    if not samples:
        return df
    qc = _expr().sample_expression_qc(code, auto_fetch=auto_fetch)
    if qc.empty:
        return df[id_columns(df)].copy()
    if mode == "pass":
        allowed = set(qc.loc[qc["sample_qc_status"] == "pass", "sample_id"].astype(str))
    else:
        allowed = set(
            qc.loc[qc["sample_qc_status"].isin(["pass", "warn"]), "sample_id"].astype(str)
        )
    keep_samples = [s for s in samples if s in allowed]
    return df[[*id_columns(df), *keep_samples]].copy()


def sample_expression_qc_from_matrix(
    raw: pd.DataFrame,
    *,
    cancer_type=None,
    source_metadata: dict[str, str | bool | None] | None = None,
    min_detected_genes: int = DEFAULT_MIN_DETECTED_GENES_FOR_QC,
    min_housekeeping_detected: int | None = None,
    min_housekeeping_fraction_above_floor: float | None = None,
    housekeeping_detection_floor_tpm: float = DEFAULT_HOUSEKEEPING_DETECTION_FLOOR_TPM,
    max_zero_fraction: float = DEFAULT_MAX_ZERO_FRACTION_FOR_QC,
    max_top_gene_fraction: float = DEFAULT_MAX_TOP_GENE_FRACTION_FOR_QC,
    max_top10_gene_fraction: float = DEFAULT_MAX_TOP10_GENE_FRACTION_FOR_QC,
) -> pd.DataFrame:
    """Per-sample QC metrics for an already-loaded raw expression matrix.

    This is the shared policy core behind the public read-path QC and the offline
    artifact rebuild. It canonicalizes gene rows before measuring sparsity so the QC
    contract is evaluated in the same gene-ID space used by expression accessors.

    Raw top-gene fractions remain in the output as assay diagnostics. Concentration
    pass/fail gates use ``top1_fraction_clean`` and ``top10_fraction_clean`` because
    clean TPM first confines protocol-sensitive RNA to its fixed compartment; raw
    fractions would systematically reject valid ribo-depleted libraries.

    Reproductive programs are reported in ``sample_qc_notes``. These annotations
    never change QC status, including for callers selecting only ``sample_qc="pass"``.
    """
    code = (
        _expr().resolve_cancer_type(cancer_type, strict=False) or cancer_type
        if cancer_type is not None
        else None
    )
    raw = (
        _expr()._canonicalize_gene_rows(raw, sample_cols=sample_columns(raw)).reset_index(drop=True)
    )
    samples = sample_columns(raw)
    if not samples:
        return pd.DataFrame()

    id_cols = id_columns(raw)
    ids = raw["Ensembl_Gene_ID"].astype(str).map(unversioned)
    panel_ids = _clean_tpm_housekeeping_panel_ids()
    panel = {str(g).split(".")[0] for g in panel_ids}
    hk_rows = ids.isin(panel)
    from .gene_families import clean_tpm_censored_gene_ids

    biological_rows = ~ids.isin(clean_tpm_censored_gene_ids())
    empty_mode = (
        "record"
        if str(code or "") in _expr().source_matrices.PARTIAL_TECHNICAL_COMPARTMENT_COHORTS
        else "warn"
    )
    clean = clean_tpm(
        raw[samples],
        gene_table=raw[id_cols],
        on_empty_compartment=empty_mode,
    )
    min_hk = (
        min_housekeeping_detected
        if min_housekeeping_detected is not None
        else _min_housekeeping_detected(panel_ids)
    )
    if source_metadata is None:
        meta = (
            cancer_reference_expression_source_metadata(str(code))
            if code is not None
            else {
                "source_cohort": None,
                "source_type": None,
                "unit": None,
                "source_scale_class": "unknown",
                "linear_tpm_comparable": False,
                "tpm_proxy": False,
            }
        )
    else:
        meta = {
            "source_cohort": None,
            "source_type": None,
            "unit": None,
            "source_scale_class": "unknown",
            "linear_tpm_comparable": False,
            "tpm_proxy": False,
            **source_metadata,
        }

    from . import library_qc

    signatures = (
        library_qc.signatures(raw, samples)
        if meta["source_scale_class"] == "linear_rnaseq_tpm"
        else None
    )
    rows: list[dict] = []
    for sample in samples:
        vals = pd.to_numeric(raw[sample], errors="coerce")
        clean_vals = pd.to_numeric(clean[sample], errors="coerce")
        measured = vals.notna()
        positive = vals > 0
        n_measured = int(measured.sum())
        n_detected = int(positive.sum())
        n_zero = int((vals == 0).sum())
        n_missing = int(vals.isna().sum())
        n_detected_clean = int((clean_vals > 0).sum())
        n_detected_clean_biological = int((clean_vals[biological_rows] > 0).sum())
        total = float(vals.sum(skipna=True))
        top_idx = vals.idxmax(skipna=True) if n_measured else None
        top_tpm = float(vals.loc[top_idx]) if top_idx is not None and pd.notna(top_idx) else 0.0
        top10_tpm = float(vals.nlargest(min(10, n_measured)).sum()) if n_measured else 0.0
        clean_total = float(clean_vals.sum(skipna=True))
        clean_top_tpm = float(clean_vals.max(skipna=True)) if n_measured else 0.0
        clean_top10_tpm = (
            float(clean_vals.nlargest(min(10, n_measured)).sum()) if n_measured else 0.0
        )
        hk_vals = vals[hk_rows]
        hk_measured = int(hk_vals.notna().sum())
        hk_detected = int((hk_vals > 0).sum())
        hk_zero = int((hk_vals == 0).sum())
        hk_floor_count = int((hk_vals >= housekeeping_detection_floor_tpm).sum())
        top_fraction = top_tpm / total if total > 0 else np.nan
        top10_fraction = top10_tpm / total if total > 0 else np.nan
        clean_top_fraction = clean_top_tpm / clean_total if clean_total > 0 else np.nan
        clean_top10_fraction = clean_top10_tpm / clean_total if clean_total > 0 else np.nan
        detected_fraction = n_detected / n_measured if n_measured else np.nan
        zero_fraction = n_zero / n_measured if n_measured else np.nan
        parse_missing_fraction = n_missing / len(vals) if len(vals) else np.nan
        hk_zero_fraction = hk_zero / hk_measured if hk_measured else np.nan
        hk_floor_fraction = hk_floor_count / hk_measured if hk_measured else np.nan

        flags: list[str] = []
        if bool(meta["linear_tpm_comparable"]):
            if n_detected < min_detected_genes:
                flags.append("low_detected_genes")
            if min_hk is not None and hk_detected < min_hk:
                flags.append("low_housekeeping_detection")
            if (
                min_housekeeping_fraction_above_floor is not None
                and pd.notna(hk_floor_fraction)
                and hk_floor_fraction < min_housekeeping_fraction_above_floor
            ):
                flags.append("low_housekeeping_floor_fraction")
            if pd.notna(zero_fraction) and zero_fraction > max_zero_fraction:
                flags.append("high_zero_fraction")
            if pd.notna(clean_top_fraction) and clean_top_fraction > max_top_gene_fraction:
                flags.append("high_top_gene_fraction")
            if (
                n_measured > 10
                and pd.notna(clean_top10_fraction)
                and clean_top10_fraction > max_top10_gene_fraction
            ):
                flags.append("high_top10_gene_fraction")
        else:
            flags.append("nonlinear_or_proxy_expression_scale")

        if signatures is not None:
            flags.extend(library_qc.reasons(signatures.loc[sample]))
        status = _sample_qc_status(flags)
        rows.append(
            {
                "cancer_code": code,
                **(signatures.loc[sample].to_dict() if signatures is not None else {}),
                "source_cohort": meta["source_cohort"],
                "source_type": meta["source_type"],
                "unit": meta["unit"],
                "source_scale_class": meta["source_scale_class"],
                "linear_tpm_comparable": bool(meta["linear_tpm_comparable"]),
                "sample_id": sample,
                "n_measured_genes": n_measured,
                "n_detected_genes": n_detected,
                "n_detected_raw": n_detected,
                "n_detected_clean": n_detected_clean,
                "n_detected_clean_biological": n_detected_clean_biological,
                "detected_gene_fraction": detected_fraction,
                "zero_fraction_raw": zero_fraction,
                "parse_missing_fraction": parse_missing_fraction,
                "total_tpm": total,
                "top_gene_id": raw.loc[top_idx, "Ensembl_Gene_ID"] if top_idx is not None else None,
                "top_gene_symbol": raw.loc[top_idx, "Symbol"] if top_idx is not None else None,
                "top_gene_tpm": top_tpm,
                "top_gene_fraction": top_fraction,
                "top10_tpm": top10_tpm,
                "top10_fraction": top10_fraction,
                "top1_fraction_clean": clean_top_fraction,
                "top10_fraction_clean": clean_top10_fraction,
                "housekeeping_genes_present": hk_measured,
                "housekeeping_genes_detected": hk_detected,
                "housekeeping_detection_floor_tpm": float(housekeeping_detection_floor_tpm),
                "housekeeping_genes_above_floor": hk_floor_count,
                "housekeeping_fraction_above_floor": hk_floor_fraction,
                "housekeeping_zero_fraction": hk_zero_fraction,
                "tpm_proxy": bool(meta["tpm_proxy"]),
                "sample_qc_status": status,
                "sample_qc_reasons": ";".join(flags),
                "sample_qc_notes": (
                    ";".join(library_qc.notes(signatures.loc[sample]))
                    if signatures is not None
                    else ""
                ),
                "passes_expression_qc": status != "fail",
                "recommended_for_absolute_tpm_floor": status != "fail"
                and bool(meta["linear_tpm_comparable"]),
            }
        )
    return pd.DataFrame(rows)


def sample_expression_qc(
    cancer_type,
    *,
    auto_fetch: bool = True,
    min_detected_genes: int = DEFAULT_MIN_DETECTED_GENES_FOR_QC,
    min_housekeeping_detected: int | None = None,
    min_housekeeping_fraction_above_floor: float | None = None,
    housekeeping_detection_floor_tpm: float = DEFAULT_HOUSEKEEPING_DETECTION_FLOOR_TPM,
    max_zero_fraction: float = DEFAULT_MAX_ZERO_FRACTION_FOR_QC,
    max_top_gene_fraction: float = DEFAULT_MAX_TOP_GENE_FRACTION_FOR_QC,
    max_top10_gene_fraction: float = DEFAULT_MAX_TOP10_GENE_FRACTION_FOR_QC,
) -> pd.DataFrame:
    """Per-sample QC metrics for a cohort's raw expression matrix.

    This is an audit surface over the raw per-sample matrix before clean-TPM
    normalization. It is designed to catch source/sample artifacts such as literal-zero
    sparsity in otherwise universal genes, while still making source-type caveats
    explicit (for example microarray TPM-proxy sources). It does not exclude samples by
    itself; downstream code can use ``passes_expression_qc`` or inspect
    ``sample_qc_reasons``.
    """
    code = _expr().resolve_cancer_type(cancer_type, strict=False) or cancer_type
    raw = _expr().per_sample_expression(
        code, normalize="tpm_raw", auto_fetch=auto_fetch, sample_qc="all"
    )
    return sample_expression_qc_from_matrix(
        raw,
        cancer_type=code,
        min_detected_genes=min_detected_genes,
        min_housekeeping_detected=min_housekeeping_detected,
        min_housekeeping_fraction_above_floor=min_housekeeping_fraction_above_floor,
        housekeeping_detection_floor_tpm=housekeeping_detection_floor_tpm,
        max_zero_fraction=max_zero_fraction,
        max_top_gene_fraction=max_top_gene_fraction,
        max_top10_gene_fraction=max_top10_gene_fraction,
    )


_HOUSEKEEPING_CANCER_COVERAGE_COLUMNS = [
    "cancer_code",
    "source_cohort",
    "source_type",
    "unit",
    "source_scale_class",
    "linear_tpm_comparable",
    "recommended_for_absolute_tpm_floor",
    "sample_qc",
    "expression_space",
    "housekeeping_detection_floor_tpm",
    "Ensembl_Gene_ID",
    "Symbol",
    "panel_member_present",
    "n_samples",
    "n_measured_samples",
    "n_detected_samples",
    "n_above_floor_samples",
    "fraction_measured",
    "fraction_detected",
    "fraction_above_floor",
    "min_tpm",
    "p1_tpm",
    "p5_tpm",
    "median_tpm",
    "mean_tpm",
    "max_tpm",
    "passes_all_samples_floor",
    "passes_p5_floor",
]

_HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_COLUMNS = [
    "Ensembl_Gene_ID",
    "Symbol",
    "housekeeping_detection_floor_tpm",
    "sample_qc",
    "expression_space",
    "linear_cohorts_audited",
    "linear_samples_audited",
    "noncomparable_cohorts_observed",
    "noncomparable_samples_observed",
    "linear_cohorts_with_measurement",
    "linear_cohorts_missing_measurement",
    "linear_cohorts_passing_p5_floor",
    "linear_measured_cohort_pass_fraction",
    "minimum_linear_p5_tpm",
    "minimum_linear_fraction_above_floor",
    "worst_linear_p5_cancer_code",
    "worst_linear_p5_source_cohort",
    "linear_floor_status",
]

_HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_REQUIRED = frozenset(
    {
        "cancer_code",
        "source_cohort",
        "recommended_for_absolute_tpm_floor",
        "housekeeping_detection_floor_tpm",
        "Ensembl_Gene_ID",
        "Symbol",
        "panel_member_present",
        "n_samples",
        "n_measured_samples",
        "fraction_above_floor",
        "p5_tpm",
        "passes_p5_floor",
        "sample_qc",
        "expression_space",
    }
)


def _housekeeping_panel_rows(panel_ids=None) -> pd.DataFrame:
    from .gene_families import clean_tpm_biological_housekeeping_genes

    if panel_ids is None:
        panel = clean_tpm_biological_housekeeping_genes()[["Ensembl_Gene_ID", "Symbol"]].copy()
    else:
        ids = [unversioned(str(g)) for g in panel_ids]
        full = clean_tpm_biological_housekeeping_genes(primary_only=False)
        symbol_by_id = dict(
            zip(full["Ensembl_Gene_ID"].astype(str).map(unversioned), full["Symbol"].astype(str))
        )
        panel = pd.DataFrame(
            {
                "Ensembl_Gene_ID": ids,
                "Symbol": [symbol_by_id.get(gid, gid) for gid in ids],
            }
        )
    panel["Ensembl_Gene_ID"] = panel["Ensembl_Gene_ID"].astype(str).map(unversioned)
    panel = panel.drop_duplicates("Ensembl_Gene_ID", keep="first").reset_index(drop=True)
    return panel


def _coverage_quantile(values: pd.Series, q: float) -> float:
    measured = pd.to_numeric(values, errors="coerce").dropna()
    if measured.empty:
        return np.nan
    return float(np.nanquantile(measured.to_numpy(dtype=float), q))


def housekeeping_cancer_expression_coverage_from_matrix(
    matrix: pd.DataFrame,
    *,
    cancer_type=None,
    source_metadata: dict[str, str | bool | None] | None = None,
    panel_ids=None,
    value_cols=None,
    housekeeping_detection_floor_tpm: float = DEFAULT_HOUSEKEEPING_DETECTION_FLOOR_TPM,
    sample_qc: str = "provided",
    expression_space: str = "tpm_clean",
) -> pd.DataFrame:
    """Per-housekeeping-gene cancer expression coverage for one cohort matrix.

    ``matrix`` should be in the expression space being audited, normally clean TPM.
    The helper does not select a new housekeeping panel; it reports the empirical
    coverage and low-tail statistics needed to evaluate one. Absolute TPM floor
    decisions should use ``recommended_for_absolute_tpm_floor``: RNA-seq-like linear
    TPM sources are comparable, while microarray/proxy sources should be warnings or
    rank-calibrated audits rather than hard floor vetoes.
    """
    if matrix is None or matrix.empty:
        return pd.DataFrame(columns=_HOUSEKEEPING_CANCER_COVERAGE_COLUMNS)

    samples = [str(c) for c in (value_cols if value_cols is not None else sample_columns(matrix))]
    samples = [c for c in samples if c in matrix.columns]
    if not samples:
        return pd.DataFrame(columns=_HOUSEKEEPING_CANCER_COVERAGE_COLUMNS)

    df = _expr()._canonicalize_gene_rows(matrix, sample_cols=samples).reset_index(drop=True)
    ids = df["Ensembl_Gene_ID"].astype(str).map(unversioned)
    row_by_id = {gid: idx for idx, gid in ids.items()}
    symbols = dict(zip(ids, df["Symbol"].astype(str)))
    panel = _housekeeping_panel_rows(panel_ids)

    code = (
        _expr().resolve_cancer_type(cancer_type, strict=False) or cancer_type
        if cancer_type is not None
        else None
    )
    if source_metadata is None:
        meta = (
            cancer_reference_expression_source_metadata(str(code))
            if code is not None
            else {
                "source_cohort": None,
                "source_type": None,
                "unit": None,
                "source_scale_class": "unknown",
                "linear_tpm_comparable": False,
            }
        )
    else:
        meta = {
            "source_cohort": None,
            "source_type": None,
            "unit": None,
            "source_scale_class": "unknown",
            "linear_tpm_comparable": False,
            **source_metadata,
        }

    raw_linear_comparability = meta.get("linear_tpm_comparable")
    linear_tpm_comparable = _expr()._optional_bool(raw_linear_comparability)
    if linear_tpm_comparable is None:
        if pd.isna(raw_linear_comparability):
            linear_tpm_comparable = False
        else:
            raise ValueError("source_metadata linear_tpm_comparable must be a boolean value")

    floor = float(housekeeping_detection_floor_tpm)
    if not np.isfinite(floor) or floor < 0:
        raise ValueError("housekeeping_detection_floor_tpm must be finite and non-negative")
    n_samples = len(samples)
    rows: list[dict] = []
    for panel_row in panel.itertuples(index=False):
        gid = unversioned(str(panel_row.Ensembl_Gene_ID))
        idx = row_by_id.get(gid)
        present = idx is not None
        values = (
            pd.to_numeric(df.loc[idx, samples], errors="coerce")
            if present
            else pd.Series([np.nan] * n_samples, index=samples, dtype=float)
        )
        measured = values.notna()
        detected = values > 0
        above_floor = values >= floor
        n_measured = int(measured.sum())
        n_detected = int(detected.sum())
        n_above_floor = int(above_floor.sum())
        symbol = symbols.get(gid, str(panel_row.Symbol))
        p5 = _coverage_quantile(values, 0.05)
        rows.append(
            {
                "cancer_code": code,
                "source_cohort": meta.get("source_cohort"),
                "source_type": meta.get("source_type"),
                "unit": meta.get("unit"),
                "source_scale_class": meta.get("source_scale_class"),
                "linear_tpm_comparable": linear_tpm_comparable,
                "recommended_for_absolute_tpm_floor": linear_tpm_comparable,
                "sample_qc": sample_qc,
                "expression_space": expression_space,
                "housekeeping_detection_floor_tpm": floor,
                "Ensembl_Gene_ID": gid,
                "Symbol": symbol,
                "panel_member_present": bool(present),
                "n_samples": n_samples,
                "n_measured_samples": n_measured,
                "n_detected_samples": n_detected,
                "n_above_floor_samples": n_above_floor,
                "fraction_measured": n_measured / n_samples if n_samples else np.nan,
                "fraction_detected": n_detected / n_samples if n_samples else np.nan,
                "fraction_above_floor": n_above_floor / n_samples if n_samples else np.nan,
                "min_tpm": _coverage_quantile(values, 0.0),
                "p1_tpm": _coverage_quantile(values, 0.01),
                "p5_tpm": p5,
                "median_tpm": _coverage_quantile(values, 0.5),
                "mean_tpm": float(values.mean(skipna=True)) if n_measured else np.nan,
                "max_tpm": _coverage_quantile(values, 1.0),
                "passes_all_samples_floor": bool(n_samples and n_above_floor == n_samples),
                "passes_p5_floor": bool(pd.notna(p5) and p5 >= floor),
            }
        )
    out = pd.DataFrame(rows, columns=_HOUSEKEEPING_CANCER_COVERAGE_COLUMNS)
    out.attrs["issue"] = "#202"
    out.attrs["sample_qc"] = sample_qc
    out.attrs["expression_space"] = expression_space
    return out


def housekeeping_cancer_expression_coverage(
    cancer_types=None,
    *,
    sample_qc: str = "pass_or_warn",
    auto_fetch: bool = True,
    panel_ids=None,
    housekeeping_detection_floor_tpm: float = DEFAULT_HOUSEKEEPING_DETECTION_FLOOR_TPM,
    on_missing: str = "empty",
) -> pd.DataFrame:
    """Audit clean-TPM housekeeping-gene coverage across cancer cohorts.

    ``cancer_types=None`` audits every cohort listed in
    :func:`oncoref.source_matrices.available_cohorts`. Rows are returned for every
    requested housekeeping-panel gene in every readable cohort, including source-scale
    metadata. Use ``recommended_for_absolute_tpm_floor`` before treating ``p5_tpm`` or
    ``fraction_above_floor`` as a hard RNA-seq TPM floor; proxy/non-linear sources are
    retained for visibility but are not comparable on an absolute TPM scale.
    """
    mode = str(on_missing).lower()
    if mode not in {"empty", "raise"}:
        raise ValueError("on_missing must be 'empty' or 'raise'")
    if cancer_types is None:
        codes = _expr().source_matrices.available_cohorts()
    elif isinstance(cancer_types, str):
        codes = [cancer_types]
    else:
        codes = list(cancer_types)

    requested_codes: list[str] = []
    seen_codes: set[str] = set()
    audited_codes: list[str] = []
    missing_codes: list[str] = []
    frames: list[pd.DataFrame] = []
    for requested in codes:
        code = str(_expr().resolve_cancer_type(requested, strict=False) or requested)
        if code in seen_codes:
            continue
        seen_codes.add(code)
        requested_codes.append(code)
        try:
            clean = _expr().per_sample_expression(
                code,
                normalize="tpm_clean",
                sample_qc=sample_qc,
                auto_fetch=auto_fetch,
            )
        except (FileNotFoundError, _expr().source_matrices.SourceMatrixError):
            if mode == "raise":
                raise
            missing_codes.append(code)
            continue
        meta = cancer_reference_expression_source_metadata(str(code))
        audited_codes.append(code)
        frames.append(
            housekeeping_cancer_expression_coverage_from_matrix(
                clean,
                cancer_type=code,
                source_metadata=meta,
                panel_ids=panel_ids,
                housekeeping_detection_floor_tpm=housekeeping_detection_floor_tpm,
                sample_qc=sample_qc,
                expression_space="tpm_clean",
            )
        )

    out = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=_HOUSEKEEPING_CANCER_COVERAGE_COLUMNS)
    )
    out.attrs["issue"] = "#202"
    out.attrs["sample_qc"] = sample_qc
    out.attrs["expression_space"] = "tpm_clean"
    out.attrs["housekeeping_detection_floor_tpm"] = float(housekeeping_detection_floor_tpm)
    out.attrs["on_missing"] = mode
    out.attrs["requested_cancer_codes"] = tuple(requested_codes)
    out.attrs["audited_cancer_codes"] = tuple(audited_codes)
    out.attrs["missing_cancer_codes"] = tuple(missing_codes)
    out.attrs["cohort_audit_complete"] = not missing_codes
    return out


def _strict_coverage_bool_column(values: pd.Series, *, column: str) -> pd.Series:
    normalized = values.astype(str).str.strip().str.lower()
    true_values = {"true", "1", "yes"}
    false_values = {"false", "0", "no"}
    invalid = normalized[~normalized.isin(true_values | false_values)]
    if not invalid.empty:
        bad = sorted(invalid.unique().tolist())
        raise ValueError(f"coverage column {column!r} contains non-boolean values: {bad}")
    return normalized.isin(true_values)


def _prepare_housekeeping_coverage_summary_input(
    coverage: pd.DataFrame,
) -> tuple[pd.DataFrame, float, dict[str, str]]:
    work = coverage.loc[:, sorted(_HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_REQUIRED)].copy()
    raw_ids = work["Ensembl_Gene_ID"]
    invalid_ids = raw_ids.isna() | raw_ids.astype(str).str.strip().str.lower().isin(
        {"", "nan", "none"}
    )
    if invalid_ids.any():
        raise ValueError("coverage contains missing Ensembl_Gene_ID values")
    work["Ensembl_Gene_ID"] = raw_ids.astype(str).map(unversioned)

    for column in (
        "housekeeping_detection_floor_tpm",
        "n_samples",
        "n_measured_samples",
        "fraction_above_floor",
        "p5_tpm",
    ):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    for column in (
        "recommended_for_absolute_tpm_floor",
        "panel_member_present",
        "passes_p5_floor",
    ):
        work[column] = _strict_coverage_bool_column(work[column], column=column)

    if work["housekeeping_detection_floor_tpm"].isna().any():
        raise ValueError("coverage has missing housekeeping_detection_floor_tpm values")
    floors = work["housekeeping_detection_floor_tpm"].unique()
    if len(floors) != 1:
        raise ValueError(
            "coverage must use exactly one housekeeping_detection_floor_tpm; "
            f"found {sorted(float(value) for value in floors)}"
        )
    floor = float(floors[0])
    if not np.isfinite(floor) or floor < 0:
        raise ValueError("housekeeping_detection_floor_tpm must be finite and non-negative")

    for column in ("n_samples", "n_measured_samples"):
        values = work[column]
        invalid = values.isna() | values.lt(0) | values.mod(1).ne(0)
        if invalid.any():
            raise ValueError(f"coverage column {column!r} must contain non-negative integers")
    if work["n_measured_samples"].gt(work["n_samples"]).any():
        raise ValueError("coverage n_measured_samples cannot exceed n_samples")
    invalid_fraction = work["fraction_above_floor"].isna() | ~work["fraction_above_floor"].between(
        0.0, 1.0
    )
    if invalid_fraction.any():
        raise ValueError("coverage fraction_above_floor values must be between 0 and 1")
    invalid_p5 = work["p5_tpm"].notna() & ~np.isfinite(work["p5_tpm"])
    if invalid_p5.any():
        raise ValueError("coverage p5_tpm values must be finite or missing")

    dimensions: dict[str, str] = {}
    for column in ("sample_qc", "expression_space"):
        values = work[column].dropna().astype(str).str.strip().unique()
        if len(values) != 1 or work[column].isna().any() or not values[0]:
            raise ValueError(f"coverage must use exactly one non-missing {column}; found {values}")
        dimensions[column] = str(values[0])

    computed_p5_pass = work["p5_tpm"].notna() & work["p5_tpm"].ge(floor)
    inconsistent_pass = work["passes_p5_floor"] != computed_p5_pass
    if inconsistent_pass.any():
        raise ValueError("coverage passes_p5_floor values disagree with p5_tpm and the floor")
    measured_without_member = ~work["panel_member_present"] & work["n_measured_samples"].gt(0)
    if measured_without_member.any():
        raise ValueError("coverage records measured samples for an absent panel member")

    audit_key = ["Ensembl_Gene_ID", "cancer_code", "source_cohort"]
    duplicates = work.duplicated(audit_key, keep=False)
    if duplicates.any():
        repeated = work.loc[duplicates, audit_key].drop_duplicates().to_dict("records")
        raise ValueError(f"coverage contains duplicate gene/cohort audit rows: {repeated[:5]}")
    return work, floor, dimensions


def _summarize_housekeeping_gene_coverage(
    gene_id: str,
    gene_rows: pd.DataFrame,
    *,
    floor: float,
    dimensions: dict[str, str],
) -> dict:
    linear = gene_rows[gene_rows["recommended_for_absolute_tpm_floor"]]
    noncomparable = gene_rows[~gene_rows["recommended_for_absolute_tpm_floor"]]
    measured = linear[
        linear["panel_member_present"]
        & linear["n_measured_samples"].gt(0)
        & linear["p5_tpm"].notna()
    ]
    missing_measurement = len(linear) - len(measured)
    passing = int(measured["passes_p5_floor"].sum())

    if linear.empty:
        status = "not_audited"
    elif missing_measurement:
        status = "missing_linear_measurement"
    elif passing == len(linear):
        status = "passes_all_linear_cohorts"
    else:
        status = "fails_some_linear_cohorts"

    worst = measured.sort_values("p5_tpm", kind="stable", na_position="last").head(1)
    worst_row = worst.iloc[0] if not worst.empty else None
    symbols = gene_rows["Symbol"].dropna().astype(str)
    symbol = symbols.iloc[0] if not symbols.empty else str(gene_id)
    return {
        "Ensembl_Gene_ID": gene_id,
        "Symbol": symbol,
        "housekeeping_detection_floor_tpm": floor,
        "sample_qc": dimensions["sample_qc"],
        "expression_space": dimensions["expression_space"],
        "linear_cohorts_audited": len(linear),
        "linear_samples_audited": int(linear["n_samples"].sum()),
        "noncomparable_cohorts_observed": len(noncomparable),
        "noncomparable_samples_observed": int(noncomparable["n_samples"].sum()),
        "linear_cohorts_with_measurement": len(measured),
        "linear_cohorts_missing_measurement": missing_measurement,
        "linear_cohorts_passing_p5_floor": passing,
        "linear_measured_cohort_pass_fraction": (
            passing / len(measured) if len(measured) else np.nan
        ),
        "minimum_linear_p5_tpm": (
            float(measured["p5_tpm"].min()) if not measured.empty else np.nan
        ),
        "minimum_linear_fraction_above_floor": (
            float(measured["fraction_above_floor"].min()) if not measured.empty else np.nan
        ),
        "worst_linear_p5_cancer_code": (
            worst_row["cancer_code"] if worst_row is not None else None
        ),
        "worst_linear_p5_source_cohort": (
            worst_row["source_cohort"] if worst_row is not None else None
        ),
        "linear_floor_status": status,
    }


def housekeeping_cancer_expression_coverage_summary(
    coverage: pd.DataFrame,
    *,
    require_complete: bool = True,
) -> pd.DataFrame:
    """Summarize the cancer-side evidence for each housekeeping candidate.

    This policy layer asks whether each candidate meets the recorded cohort p5
    floor across sources whose scale supports absolute clean-TPM comparison.
    Non-comparable sources remain counted, but their TPM values cannot pass or veto
    the floor. ``linear_floor_status`` separates complete passes, measured low-tail
    failures, missing measurements, and candidates with no comparable audit.

    The statuses are panel-review evidence, not an automatic promotion rule. A panel
    change alters every HK-normalized value and still requires a consumer-specific
    benchmark. Results from :func:`housekeeping_cancer_expression_coverage` declare
    whether every requested cohort was audited. Partial or undeclared inputs are
    rejected by default; use ``require_complete=False`` only for explicit exploration.
    """
    missing = sorted(_HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_REQUIRED - set(coverage.columns))
    if missing:
        raise ValueError(f"coverage missing required columns: {missing}")

    completeness = _expr()._optional_bool(coverage.attrs.get("cohort_audit_complete"))
    raw_missing_codes = coverage.attrs.get("missing_cancer_codes", ())
    if isinstance(raw_missing_codes, str):
        missing_codes = (raw_missing_codes,)
    else:
        missing_codes = tuple(str(code) for code in raw_missing_codes)
    if completeness is True and missing_codes:
        raise ValueError("coverage completeness metadata conflicts with missing_cancer_codes")
    if require_complete and completeness is not True:
        detail = (
            f"missing requested cancer codes: {list(missing_codes)}"
            if missing_codes
            else "input does not declare whether every requested cohort was audited"
        )
        raise ValueError(
            f"housekeeping coverage audit is not confirmed complete; {detail}. "
            "Use housekeeping_cancer_expression_coverage to produce a complete audit "
            "or pass require_complete=False for an explicitly partial audit."
        )
    if coverage.empty:
        out = pd.DataFrame(columns=_HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_COLUMNS)
        out.attrs.update(coverage.attrs)
        out.attrs["issue"] = "#202"
        out.attrs["panel_selection_policy"] = "consumer_benchmark_required"
        out.attrs["missing_cancer_codes"] = missing_codes
        out.attrs["cohort_audit_complete"] = completeness
        return out

    work, floor, dimensions = _prepare_housekeeping_coverage_summary_input(coverage)
    rows = [
        _summarize_housekeeping_gene_coverage(
            str(gene_id), gene_rows, floor=floor, dimensions=dimensions
        )
        for gene_id, gene_rows in work.groupby("Ensembl_Gene_ID", sort=False)
    ]

    out = pd.DataFrame(rows, columns=_HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_COLUMNS)
    out.attrs["issue"] = "#202"
    out.attrs["housekeeping_detection_floor_tpm"] = floor
    out.attrs["panel_selection_policy"] = "consumer_benchmark_required"
    out.attrs["sample_qc"] = dimensions["sample_qc"]
    out.attrs["expression_space"] = dimensions["expression_space"]
    out.attrs["requested_cancer_codes"] = coverage.attrs.get("requested_cancer_codes")
    out.attrs["audited_cancer_codes"] = coverage.attrs.get("audited_cancer_codes")
    out.attrs["missing_cancer_codes"] = missing_codes
    out.attrs["cohort_audit_complete"] = completeness
    return out


def _validate_metadata_on_missing(on_missing: str) -> str:
    mode = str(on_missing).lower()
    if mode not in {"empty", "raise"}:
        raise ValueError("on_missing must be 'empty' or 'raise'")
    return mode


def _empty_metadata_frame(columns: list[str], *, schema_version: str, missing_reason: str):
    out = pd.DataFrame(columns=columns)
    if schema_version:
        out.attrs["schema_version"] = schema_version
    out.attrs["data_version"] = DATA_VERSION
    out.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    out.attrs["missing_reason"] = missing_reason
    return out


def _optional_bundle_metadata_path(relative_path: str, *, auto_fetch: bool) -> Path | None:
    path = data_bundle.find(relative_path)
    if path is not None:
        return path
    if auto_fetch:
        data_bundle.ensure_local(auto_fetch=True, verbose=False)
        path = data_bundle.find(relative_path)
    return path


def _requested_cancer_codes(cancer_type) -> set[str] | None:
    if cancer_type is None:
        return None
    if isinstance(cancer_type, str):
        values = [cancer_type]
    else:
        values = list(cancer_type)
    return {_expr().resolve_cancer_type(value, strict=False) or str(value) for value in values}


def source_matrix_sample_qc_manifest(
    cancer_type=None,
    *,
    sample_qc: str = "all",
    auto_fetch: bool = True,
    on_missing: str = "empty",
) -> pd.DataFrame:
    """Read the generated source-matrix sample-QC manifest from the data bundle.

    This is the bundle-level companion to :func:`sample_expression_qc`, which computes
    QC live from a single source matrix. The rebuild script emits
    ``source-matrix-sample-qc.csv`` for the exact samples that fed generated
    expression artifacts. Current released bundles may not yet contain it; by default
    this returns an empty schema-stable frame rather than synthesizing rows.
    """
    mode = _validate_metadata_on_missing(on_missing)
    path = _optional_bundle_metadata_path(
        SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH, auto_fetch=auto_fetch
    )
    if path is None:
        if mode == "raise":
            raise FileNotFoundError(
                f"{SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH} is not present in the "
                f"oncoref data bundle at {data_bundle.cache_dir()}"
            )
        return _empty_metadata_frame(
            _SOURCE_MATRIX_SAMPLE_QC_MANIFEST_COLUMNS,
            schema_version="",
            missing_reason="source-matrix sample QC manifest not present in bundle",
        )

    out = pd.read_csv(path, low_memory=False)
    for col in _SOURCE_MATRIX_SAMPLE_QC_MANIFEST_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    extra_cols = [
        col for col in out.columns if col not in _SOURCE_MATRIX_SAMPLE_QC_MANIFEST_COLUMNS
    ]
    out = out[[*_SOURCE_MATRIX_SAMPLE_QC_MANIFEST_COLUMNS, *extra_cols]]

    codes = _requested_cancer_codes(cancer_type)
    if codes is not None:
        if not codes:
            out = out.iloc[0:0].copy()
        else:
            out = out[out["cancer_code"].astype(str).isin(codes)].copy()
    qc_mode = _validate_sample_qc(sample_qc)
    if qc_mode == "pass":
        out = out[out["sample_qc_status"].astype(str) == "pass"].copy()
    elif qc_mode == "pass_or_warn":
        out = out[out["sample_qc_status"].astype(str).isin(["pass", "warn"])].copy()

    out.attrs["data_version"] = DATA_VERSION
    out.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    out.attrs["path"] = str(path)
    return out.reset_index(drop=True)


def expression_artifact_build_metadata(
    cancer_type=None,
    *,
    auto_fetch: bool = True,
    on_missing: str = "empty",
) -> pd.DataFrame:
    """Read per-cohort build metadata for generated expression artifacts.

    ``source_cohort`` is the canonical cohort identity used by the source-matrix and
    expression artifact loaders. ``build_source_cohort`` preserves the cohort label
    physically recorded when the bundle was built, which can differ in older bundles.
    The remaining fields record the selected source matrix, QC policy, selected sample
    count, QC pass/warn/fail counts, and representative partition coverage. Missing
    current-bundle metadata returns a schema-stable empty frame by default.
    """
    mode = _validate_metadata_on_missing(on_missing)
    path = _optional_bundle_metadata_path(
        EXPRESSION_ARTIFACT_BUILD_METADATA_PATH, auto_fetch=auto_fetch
    )
    if path is None:
        if mode == "raise":
            raise FileNotFoundError(
                f"{EXPRESSION_ARTIFACT_BUILD_METADATA_PATH} is not present in the "
                f"oncoref data bundle at {data_bundle.cache_dir()}"
            )
        return _empty_metadata_frame(
            _EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS,
            schema_version=EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
            missing_reason="expression artifact build metadata not present in bundle",
        )

    out = pd.read_csv(path)
    if "source_cohort" not in out.columns:
        out["source_cohort"] = pd.NA
    if "build_source_cohort" not in out.columns:
        out["build_source_cohort"] = out["source_cohort"]
    else:
        build_source = out["build_source_cohort"].astype("string")
        missing_build_source = build_source.isna() | build_source.str.strip().eq("")
        out.loc[missing_build_source, "build_source_cohort"] = out.loc[
            missing_build_source, "source_cohort"
        ]
    for col in _EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA

    registry = _expr().source_matrices.registry()
    loader_source_by_code = dict(
        zip(registry["cancer_code"].astype(str), registry["source_cohort"].astype(str))
    )
    loader_sources = out["cancer_code"].astype(str).map(loader_source_by_code)
    out["source_cohort"] = loader_sources.where(loader_sources.notna(), out["source_cohort"])
    extra_cols = [
        col for col in out.columns if col not in _EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS
    ]
    out = out[[*_EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS, *extra_cols]]
    codes = _requested_cancer_codes(cancer_type)
    if codes is not None:
        if not codes:
            out = out.iloc[0:0].copy()
        else:
            out = out[out["cancer_code"].astype(str).isin(codes)].copy()
    out.attrs["schema_version"] = EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION
    out.attrs["data_version"] = DATA_VERSION
    out.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    out.attrs["path"] = str(path)
    return out.reset_index(drop=True)


def _artifact_build_metadata_qc_value(row) -> str:
    for col in ("sample_qc_effective", "sample_qc"):
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            return str(value).strip().lower()
    return ""


def _artifact_build_metadata_qc_policies(meta: pd.DataFrame) -> dict[str, tuple[str, ...]]:
    """Return the recorded effective QC policies for each artifact cohort."""
    if meta.empty or "cancer_code" not in meta.columns:
        return {}
    policies: dict[str, tuple[str, ...]] = {}
    for code, group in meta.groupby("cancer_code", sort=False):
        values = sorted(
            {_artifact_build_metadata_qc_value(row) for _, row in group.iterrows()} - {""}
        )
        policies[str(code)] = tuple(values)
    return policies


def _artifact_effective_sample_qc(
    code: str,
    policies: dict[str, tuple[str, ...]],
) -> tuple[str | None, str]:
    """Resolve one live-sample QC policy for an artifact-backed cohort."""
    recorded = policies.get(code, ())
    if not recorded:
        return None, "artifact_build_metadata_missing"
    if len(recorded) != 1:
        return None, "artifact_sample_qc_conflict"
    effective = recorded[0]
    if effective not in _SAMPLE_QC_MODES:
        return None, "artifact_sample_qc_invalid"
    return effective, ""


def _artifact_build_metadata_sample_counts(meta: pd.DataFrame) -> dict[str, int]:
    """Return the selected sample count used to build each artifact cohort."""
    if meta.empty or not {"cancer_code", "n_cohort_samples"} <= set(meta.columns):
        return {}
    counts: dict[str, int] = {}
    for code, group in meta.groupby("cancer_code", sort=False):
        values = pd.to_numeric(group["n_cohort_samples"], errors="coerce").dropna().unique()
        if len(values) > 1:
            raise ValueError(f"artifact build metadata has conflicting sample counts for {code}")
        if len(values) == 1:
            value = float(values[0])
            if value < 0 or not value.is_integer():
                raise ValueError(f"artifact build metadata has invalid sample count for {code}")
            counts[str(code)] = int(value)
    return counts


def _artifact_qc_attrs(
    meta: pd.DataFrame,
    *,
    requested_sample_qc: str,
) -> dict[str, object]:
    values = sorted({_artifact_build_metadata_qc_value(row) for _, row in meta.iterrows()} - {""})
    missing_reason = str(meta.attrs.get("missing_reason", ""))
    return {
        "sample_qc": requested_sample_qc,
        "artifact_sample_qc": ",".join(values) if values else "unknown",
        "artifact_sample_qc_verified": bool(values),
        "artifact_build_metadata_n": len(meta),
        "artifact_build_metadata_missing_reason": missing_reason,
    }


def _require_expression_artifact_sample_qc(
    codes: Iterable[str],
    *,
    sample_qc: str,
) -> pd.DataFrame:
    requested = _validate_artifact_sample_qc(sample_qc)
    code_list = list(dict.fromkeys(str(c) for c in codes))
    # Do not fetch as part of read-path validation. The artifact reader already picked a
    # shard root; validation should inspect metadata that is present beside that bundle,
    # not populate an override/test cache with the global release as a side effect.
    meta = _expr().expression_artifact_build_metadata(
        code_list, auto_fetch=False, on_missing="empty"
    )
    if requested == "artifact" or not code_list:
        return meta

    # Legacy/current-development bundles may not yet ship build metadata. Keep those
    # readable, but mark attrs as unverified rather than pretending the QC policy is known.
    if meta.empty and meta.attrs.get("missing_reason"):
        return meta

    policies = _artifact_build_metadata_qc_policies(meta)
    missing = [code for code in code_list if code not in policies]
    if missing:
        raise ValueError(
            "expression artifact build metadata lacks rows for requested cohort(s): "
            + ", ".join(missing)
            + "; pass sample_qc='artifact' only for explicit legacy/audit reads"
        )

    mismatches = [
        f"{code}={','.join(policies[code]) or 'unknown'}"
        for code in code_list
        if policies[code] != (requested,)
    ]
    if mismatches:
        raise ValueError(
            "expression artifact sample_qc mismatch for requested cohort(s): "
            + ", ".join(mismatches)
            + f"; requested {requested!r}. Regenerate the bundle with that policy or "
            "pass sample_qc='artifact' for an explicit legacy/audit read."
        )
    return meta


def expression_artifact_build_summary(
    *,
    auto_fetch: bool = True,
    on_missing: str = "empty",
) -> dict:
    """Read bundle-level expression artifact build metadata.

    Returns the JSON summary emitted beside the per-cohort metadata CSV, including the
    representative partition policy and unique source-group role counts. The summary
    remains optional until a regenerated data bundle ships it; callers that require it
    should pass ``on_missing="raise"``.
    """
    mode = _validate_metadata_on_missing(on_missing)
    path = _optional_bundle_metadata_path(
        EXPRESSION_ARTIFACT_BUILD_METADATA_JSON_PATH, auto_fetch=auto_fetch
    )
    if path is None:
        if mode == "raise":
            raise FileNotFoundError(
                f"{EXPRESSION_ARTIFACT_BUILD_METADATA_JSON_PATH} is not present in the "
                f"oncoref data bundle at {data_bundle.cache_dir()}"
            )
        return {
            "schema_version": EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
            "data_version": DATA_VERSION,
            "source_matrix_version": SOURCE_MATRIX_VERSION,
            "missing_reason": "expression artifact build summary not present in bundle",
        }
    with path.open() as handle:
        out = json.load(handle)
    out.setdefault("schema_version", EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION)
    out.setdefault("data_version", DATA_VERSION)
    out.setdefault("source_matrix_version", SOURCE_MATRIX_VERSION)
    out["path"] = str(path)
    return out
