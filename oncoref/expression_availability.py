# Licensed under the Apache License, Version 2.0.
"""Cohort expression availability, non-comparable scale diagnostics, and umbrella coverage (#756)."""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from functools import cache, lru_cache

import pandas as pd

from .load_dataset import _register_derived_cache
from .version import DATA_VERSION, SOURCE_MATRIX_VERSION


def _expr():
    from . import expression

    return expression


def cancer_reference_expression_availability(
    cancer_types: str | Iterable[str] | None = None,
    normalize: str | Iterable[str] = "tpm_clean",
    *,
    sample_qc: str = "pass",
    reference_source: str = "artifact",
    source_kind: str | Iterable[str] | None = None,
    source_cohort: str | Iterable[str] | None = None,
    exclude_microarray_proxy: bool = False,
    all_sources: bool = False,
    include_aggregates: bool = False,
) -> pd.DataFrame:
    """Availability/provenance table for :func:`cancer_reference_expression`.

    By default the result is one row per requested cancer code (expanding computed
    aggregate requests to member cohorts) and normalization mode. With
    ``all_sources=True`` and ``reference_source="summary_rows_all"``, it instead
    returns one row per source cohort and normalization without loading the
    multi-million-row expression frame. Both modes are gene-independent, so callers
    can distinguish an empty gene filter from unavailable upstream data. For
    percentile artifacts and raw artifact-policy requests, ``artifact_sample_qc``
    reports the build policy recorded in the compact metadata manifest;
    ``sample_qc_effective`` reports the policy applied by that row. ``n_samples`` mirrors
    ``n_reference_samples``. With current build metadata it is the selected sample
    count used to build the artifact, not necessarily the wider raw source-matrix
    count; legacy bundles without build metadata fall back to the registered matrix
    width.
    """
    modes = _reference_normalize_modes(normalize)
    sample_qc = _expr()._validate_artifact_sample_qc(sample_qc)
    reference_source = _validate_reference_source(reference_source)
    source_kinds = _expr()._normalize_source_filter_values(source_kind)
    source_cohorts = _expr()._normalize_source_cohort_filter_values(source_cohort)
    has_source_filter = source_kinds is not None or source_cohorts is not None
    if sample_qc == "artifact" and reference_source != "artifact":
        raise ValueError("sample_qc='artifact' requires reference_source='artifact'")
    if reference_source == "summary_rows_all" and sample_qc != "all":
        raise ValueError('reference_source="summary_rows_all" requires sample_qc="all"')
    if all_sources and reference_source != "summary_rows_all":
        raise ValueError('all_sources=True requires reference_source="summary_rows_all"')
    if (has_source_filter or exclude_microarray_proxy) and (reference_source != "summary_rows_all"):
        raise ValueError(
            "source_kind, source_cohort, and exclude_microarray_proxy are only "
            "supported with reference_source='summary_rows_all'"
        )
    if all_sources:
        if include_aggregates:
            raise ValueError("include_aggregates=True is not supported with all_sources=True")
        return _reference_source_union_availability(
            cancer_types,
            modes,
            source_kinds=source_kinds,
            source_cohorts=source_cohorts,
            exclude_microarray_proxy=exclude_microarray_proxy,
        )
    requests = _reference_expression_requests(
        cancer_types,
        modes,
        reference_source=reference_source,
        sample_qc=sample_qc,
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
        describe_invalid=True,
        include_aggregates=include_aggregates,
    )
    return _reference_expression_availability_for_requests(
        requests,
        modes,
        sample_qc=sample_qc,
        reference_source=reference_source,
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )


_REFERENCE_NORMALIZE_ALIASES = {
    "clean_tpm": "tpm_clean",
    "tpm_clean": "tpm_clean",
    "clean_tpm_biological": "tpm_clean_biological",
    "tpm_clean_biological": "tpm_clean_biological",
    "biological_clean_tpm": "tpm_clean_biological",
    "tpm_clean_log1p": "tpm_clean_log1p",
    "clean_tpm_log1p": "tpm_clean_log1p",
    "tpm": "tpm_raw",
    "raw_tpm": "tpm_raw",
    "tpm_raw": "tpm_raw",
}

_REFERENCE_NORMALIZE_LABELS = {
    "tpm_clean": "tpm_clean",
    "tpm_clean_biological": "tpm_clean_biological",
    "tpm_clean_log1p": "tpm_clean_log1p",
    "tpm_raw": "tpm_raw",
}

_REFERENCE_WIDE_SUFFIXES = {
    "tpm_clean": "TPM_clean",
    "tpm_clean_biological": "TPM_clean_biological",
    "tpm_clean_log1p": "TPM_clean_log1p",
    "tpm_raw": "TPM_raw",
}

_REFERENCE_PROVENANCE_COLUMNS = [
    "source_cohort",
    "source_project",
    "source_version",
    "source_pmid",
    "source_type",
    "source_unit",
    "source_scale_class",
    "linear_tpm_comparable",
    "tumor_origin",
    "metastasis_site",
    "n_reference_genes",
    "n_reference_samples",
    "n_samples",
    "n_detected",
    "processing_pipeline",
    "notes",
    "reference_method",
    "sample_qc",
    "sample_qc_effective",
    "data_version",
    "source_matrix_version",
]

_REFERENCE_SOURCE_UNION_IDENTITY_COLUMNS = [
    "source_cohort",
    "n_reference_samples",
    "n_samples",
    "n_detected",
]

_REFERENCE_SOURCES = {"artifact", "summary_rows", "summary_rows_all"}
_REFERENCE_SUMMARY_DATASET = "cancer-reference-expression"
_REFERENCE_SUMMARY_AVAILABILITY_DATASET = "cancer-reference-expression-availability"
_REFERENCE_SUMMARY_METHOD = "source_summary_rows"
_REFERENCE_SUMMARY_ALL_METHOD = "source_summary_rows_all"
_REFERENCE_RECOMPUTED_METHOD = "source_matrix_stats"
_REFERENCE_TUMOR_ORIGIN_RANK = {
    "primary": 0,
    "mixed": 1,
    "metastasis": 2,
}

_REFERENCE_REQUEST_COLUMNS = [
    "requested_code",
    "request_kind",
    "available",
    "missing_reason",
]

_REFERENCE_AVAILABILITY_COLUMNS = [
    "requested_code",
    "cancer_code",
    "request_kind",
    "normalization",
    "available",
    "missing_reason",
    "source_cohort",
    "source_project",
    "source_version",
    "source_pmid",
    "source_type",
    "source_unit",
    "source_scale_class",
    "linear_tpm_comparable",
    "tumor_origin",
    "metastasis_site",
    "n_reference_genes",
    "n_reference_samples",
    "n_samples",
    "n_detected",
    "processing_pipeline",
    "notes",
    "reference_method",
    "sample_qc",
    "sample_qc_effective",
    "artifact_sample_qc",
    "artifact_schema_version",
    "data_version",
    "source_matrix_version",
]


def _reference_expression_requests(
    cancer_types: str | Iterable[str] | None,
    modes: list[str],
    *,
    reference_source: str,
    sample_qc: str,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
    describe_invalid: bool = False,
    include_aggregates: bool = False,
) -> list[dict[str, str]]:
    """Resolved request rows while preserving aggregate-vs-direct intent.

    The default listing omits cohorts whose selected source has an invalid expression
    scale (with a warning); naming one explicitly is rejected when its data is read.
    ``describe_invalid=True`` keeps them so availability tables can report the problem.
    ``include_aggregates=True`` appends explicit ``aggregate_member`` expansion rows for
    every ``member_union`` umbrella code when ``cancer_types=None``.
    """
    if cancer_types is None:
        invalid = {} if describe_invalid else _expr()._invalid_selected_source_codes()
        if invalid:
            warnings.warn(
                f"Omitting cohort(s) {sorted(invalid)} from the default listing: "
                f"invalid source expression scale (oncoref#604)",
                UserWarning,
                stacklevel=3,
            )
        rows = [
            {"requested_code": code, "cancer_code": code, "request_kind": "default_available"}
            for code in _expr()._reference_available_codes_for_modes(
                modes,
                reference_source=reference_source,
                sample_qc=sample_qc,
                source_kinds=source_kinds,
                source_cohorts=source_cohorts,
                exclude_microarray_proxy=exclude_microarray_proxy,
            )
            if code not in invalid
        ]
        if include_aggregates:
            aggregates = _expr().cohort_aggregates()
            for umbrella in sorted(aggregates):
                if _expr().cancer_type_reference_source(umbrella) != "member_union":
                    continue
                for member in aggregates[umbrella]:
                    if member in invalid:
                        continue
                    rows.append(
                        {
                            "requested_code": umbrella,
                            "cancer_code": member,
                            "request_kind": "aggregate_member",
                        }
                    )
        return rows
    raw_values = [cancer_types] if isinstance(cancer_types, str) else list(cancer_types)
    aggregates = _expr().cohort_aggregates()
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for raw in raw_values:
        resolved = _expr().resolve_cancer_type(raw)
        members = aggregates.get(str(raw)) or aggregates.get(resolved)
        if members and _expr().cancer_type_reference_source(resolved) == "member_union":
            expanded = [(member, "aggregate_member") for member in members]
        else:
            expanded = [(resolved, "direct")]
        for code, kind in expanded:
            key = (resolved, code)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"requested_code": resolved, "cancer_code": code, "request_kind": kind})
    return rows


def _reference_available_codes_for_modes(
    modes: list[str],
    *,
    reference_source: str,
    sample_qc: str,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
) -> list[str]:
    if reference_source in {"summary_rows", "summary_rows_all"}:
        if sample_qc == "all":
            return _expr()._reference_summary_available_codes(
                source_kinds=source_kinds,
                source_cohorts=source_cohorts,
                exclude_microarray_proxy=exclude_microarray_proxy,
            )
        return sorted(_expr().source_matrices.available_cohorts())
    out: set[str] = set()
    if any(mode in {"tpm_clean", "tpm_clean_biological", "tpm_clean_log1p"} for mode in modes):
        out.update(_expr().available_percentile_cohorts())
    if "tpm_raw" in modes:
        out.update(_expr().source_matrices.available_cohorts())
    return sorted(out)


def _reference_expression_availability_for_requests(
    requests: list[dict[str, str]],
    modes: list[str],
    *,
    sample_qc: str,
    reference_source: str,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
) -> pd.DataFrame:
    percentile_available = (
        set(_expr().available_percentile_cohorts()) if reference_source == "artifact" else set()
    )
    needs_source_matrix = reference_source == "artifact" or (
        reference_source == "summary_rows" and sample_qc != "all"
    )
    source_matrix_available = (
        set(_expr().source_matrices.available_cohorts()) if needs_source_matrix else set()
    )
    summary_available = (
        set(
            _expr()._reference_summary_available_codes(
                source_kinds=source_kinds,
                source_cohorts=source_cohorts,
                exclude_microarray_proxy=exclude_microarray_proxy,
            )
        )
        if reference_source in {"summary_rows", "summary_rows_all"}
        else set()
    )
    artifact_modes = {"tpm_clean", "tpm_clean_biological", "tpm_clean_log1p"}
    artifact_qc_policies: dict[str, tuple[str, ...]] = {}
    artifact_sample_counts: dict[str, int] = {}
    artifact_qc_metadata_present = False
    needs_artifact_qc_metadata = any(mode in artifact_modes for mode in modes) or (
        sample_qc == "artifact" and "tpm_raw" in modes
    )
    if reference_source == "artifact" and needs_artifact_qc_metadata:
        artifact_qc_meta = _expr().expression_artifact_build_metadata(
            auto_fetch=False, on_missing="empty"
        )
        artifact_qc_policies = _expr()._artifact_build_metadata_qc_policies(artifact_qc_meta)
        artifact_sample_counts = _expr()._artifact_build_metadata_sample_counts(artifact_qc_meta)
        artifact_qc_metadata_present = not bool(artifact_qc_meta.attrs.get("missing_reason"))
    rows: list[dict] = []
    for request in requests:
        code = request["cancer_code"]
        for mode in modes:
            recorded_policies = artifact_qc_policies.get(code, ())
            reports_artifact_policy = reference_source == "artifact" and (
                mode in artifact_modes or sample_qc == "artifact"
            )
            effective_sample_qc: str | None = sample_qc
            if sample_qc == "artifact":
                effective_sample_qc, _ = _expr()._artifact_effective_sample_qc(
                    code, artifact_qc_policies
                )
            available, missing_reason = _reference_mode_availability(
                code,
                mode,
                percentile_available,
                source_matrix_available,
                summary_available,
                sample_qc=sample_qc,
                reference_source=reference_source,
                artifact_qc_policies=artifact_qc_policies,
                artifact_qc_metadata_present=artifact_qc_metadata_present,
            )
            method = _reference_expected_method(
                mode, sample_qc=sample_qc, reference_source=reference_source
            )
            row = {
                "requested_code": request["requested_code"],
                "cancer_code": code,
                "request_kind": request["request_kind"],
                "normalization": _REFERENCE_NORMALIZE_LABELS[mode],
                "available": bool(available),
                "missing_reason": "" if available else missing_reason,
                "reference_method": method,
                "artifact_schema_version": _expr().REFERENCE_EXPRESSION_SCHEMA_VERSION,
                "data_version": DATA_VERSION,
                "source_matrix_version": SOURCE_MATRIX_VERSION,
                "artifact_sample_qc": (
                    ",".join(recorded_policies)
                    if reports_artifact_policy and recorded_policies
                    else pd.NA
                ),
            }
            row.update(
                _expr()._reference_expression_provenance(
                    code,
                    method,
                    sample_qc=sample_qc,
                    sample_qc_effective=effective_sample_qc,
                    reference_source=reference_source,
                    source_kinds=source_kinds,
                    source_cohorts=source_cohorts,
                    exclude_microarray_proxy=exclude_microarray_proxy,
                )
            )
            if reports_artifact_policy:
                selected_count = artifact_sample_counts.get(code)
                if selected_count is not None:
                    row["n_reference_samples"] = selected_count
                    row["n_samples"] = selected_count
            rows.append(row)
    return _reference_availability_frame(rows)


def _reference_availability_frame(rows: list[dict]) -> pd.DataFrame:
    """Build the stable public availability schema and attach version metadata."""
    out = pd.DataFrame(rows, columns=_REFERENCE_AVAILABILITY_COLUMNS)
    out.attrs["artifact_schema_version"] = _expr().REFERENCE_EXPRESSION_SCHEMA_VERSION
    out.attrs["data_version"] = DATA_VERSION
    out.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    out.attrs["issues"] = ["#207"]
    return out


_NON_COMPARABLE_REPLACEMENT_PLANS: dict[str, str] = {
    "CMN": "Replace GPL96 microarray TPM-proxy (GSE11482) with pediatric renal bulk RNA-seq when publicly available (#737).",
    "CTCL": "Supplement or replace ECCITE-seq donor pseudobulk nTPM (GSE171811) with sorted malignant T-cell bulk RNA-seq (#737).",
    "EPN": "Supplement diagnosis malignant-cell scRNA pseudobulk (GSE141460) with OpenPBTA/CBTN primary bulk RNA-seq (#737).",
    "HCL": "Supplement malignant-cell scRNA pseudobulk nTPM (Zenodo 14917813) with sorted primary HCL bulk RNA-seq (#737).",
    "MTC": "Replace GPL6480 microarray TPM-proxy (GSE32662) with primary medullary thyroid bulk RNA-seq (#737).",
    "NUTM": "Rebuild controlled 3-case paired-FASTQ cohort with scripts/build_local_fastq_salmon_source.py (#604, #665).",
    "SARC_ESS_HG": "Replace GPL22303 microarray TPM-proxy (GSE85383) with fusion-annotated HG-ESS bulk RNA-seq (#737).",
    "SARC_ESS_LG": "Replace GPL22303 microarray TPM-proxy (GSE85383) with fusion-annotated LG-ESS bulk RNA-seq (#737).",
    "SARC_MYXLPS": "Replace GPL96 microarray TPM-proxy (GSE30929) with FUS::DDIT3 bulk RNA-seq (#737).",
}

_NON_COMPARABLE_PUBLIC_RNASEQ_AUDIT: dict[str, str] = {
    "CMN": "no_public_primary_bulk_rnaseq (GEO RNA-seq audit 2026-11-22: 0 cohorts; Treehouse RiboD has n=1 unselected)",
    "CTCL": "sorted_malignant_scrna_preferred_over_unpurified_bulk (skin biopsy bulk mixes reactive T cells and keratinocytes)",
    "EPN": "diagnosis_malignant_scrna_retained_pending_subtype_resolved_openpbta_split",
    "HCL": "sorted_malignant_scrna_only_public_cohort (Zenodo 14917813; no public sorted bulk RNA-seq)",
    "MTC": "controlled_access_only (EGAS00001001473 primary MTC RNA-seq; public GEO GSE114069/GSE114070 are TT/MZ-CRC-1 cell lines)",
    "NUTM": "controlled_local_fastq_builder_ready (scripts/build_local_fastq_salmon_source.py)",
    "SARC_ESS_HG": "no_public_hg_ess_bulk_rnaseq (GEO audit 2026-11-22: 0 public HG-ESS bulk RNA-seq matrices)",
    "SARC_ESS_LG": "underpowered_public_rnaseq_only (GSE196757 has 2 patients with technical replicates vs n=9 in GSE85383)",
    "SARC_MYXLPS": "unsubtyped_3prime_tag_seq_only (GSE352573 uses Lexogen QuantSeq 3' tag-seq without histological subtype labels)",
}


def non_comparable_expression_cohorts(*, include_invalid: bool = True) -> pd.DataFrame:
    """Structured inventory of selected expression cohorts that are not linear-TPM comparable (#737).

    Separates withheld invalid-scale sources (``NUTM`` / ``invalid_legacy_tpm``) from
    rank-preserving proxy cohorts (microarray TPM-proxies, scRNA malignant-cell
    pseudobulk nTPM, and gene-length-unnormalized nCPM proxies) so callers never pool
    them into absolute-TPM thresholds by accident.
    """
    table = _expr()._reference_summary_source_table()
    sub = table.loc[
        table["selected"] & (~table["linear_tpm_comparable"].fillna(False).astype(bool))
    ].copy()
    if not include_invalid:
        sub = sub.loc[
            ~sub["source_scale_class"].astype("string").isin(_expr().INVALID_SOURCE_SCALE_CLASSES)
        ].copy()

    def _scale_family(scale_cls: str) -> str:
        s = str(scale_cls or "").lower()
        if s in _expr().INVALID_SOURCE_SCALE_CLASSES:
            return "invalid_withheld"
        if "microarray" in s:
            return "microarray_tpm_proxy"
        if "scrna" in s or "pseudobulk" in s:
            return "scrna_pseudobulk_ntpm"
        if "cpm" in s:
            return "bulk_rnaseq_cpm_proxy"
        return "other_non_comparable"

    rows = []
    for row in sub.sort_values("cancer_code").itertuples(index=False):
        code = str(row.cancer_code)
        scale_cls = str(row.source_scale_class)
        is_invalid = scale_cls in _expr().INVALID_SOURCE_SCALE_CLASSES
        rows.append(
            {
                "cancer_code": code,
                "source_cohort": str(row.source_cohort),
                "source_type": str(row.source_type),
                "source_scale_class": scale_cls,
                "scale_family": _scale_family(scale_cls),
                "linear_tpm_comparable": False,
                "is_invalid_withheld": is_invalid,
                "n_reference_samples": int(row.n_reference_samples)
                if pd.notna(row.n_reference_samples)
                else 0,
                "allowed_uses": (
                    "none_withheld" if is_invalid else "within_sample_rank_and_pathway_only"
                ),
                "public_rnaseq_audit_status": _NON_COMPARABLE_PUBLIC_RNASEQ_AUDIT.get(
                    code, "pending_audit"
                ),
                "replacement_plan": _NON_COMPARABLE_REPLACEMENT_PLANS.get(
                    code, "Audit linear bulk RNA-seq replacement (#737)."
                ),
            }
        )
    return pd.DataFrame(rows)


def umbrella_expression_coverage() -> pd.DataFrame:
    """Explicit member-level coverage and subtype-gap audit for umbrella cancer codes (#737).

    Covers all ``member_union`` umbrella codes (including the 14 classification-target
    umbrellas such as ``THYM_EPITHELIAL``, ``NSCLC``, ``CRC``, ``SARC_ESS``, ``SARC_LPS``,
    ``NEN``, ``NET``, ``NEC``) and discloses child subtype gaps (such as ``THYMCA``
    under ``THYM_EPITHELIAL``, where ``TCGA-THYM`` represents thymoma rather than a
    dedicated thymic carcinoma RNA-seq cohort) as well as member cohorts on proxy scales.
    """
    from .cancer_types import expression_reference_coverage

    cov = expression_reference_coverage().set_index("code", drop=False)
    aggregates = _expr().cohort_aggregates()
    src_tbl = _expr()._reference_summary_source_table()
    selected = src_tbl.loc[src_tbl["selected"]].set_index("cancer_code")
    avail = cancer_reference_expression_availability(
        reference_source="summary_rows",
        sample_qc="all",
        include_aggregates=True,
    )
    avail_by_pair = (
        avail.loc[avail["request_kind"].eq("aggregate_member")]
        .set_index(["requested_code", "cancer_code"])["available"]
        .to_dict()
    )

    rows = []
    for code in sorted(cov.index):
        rec = cov.loc[code]
        ref_source = str(rec["reference_source"])
        if ref_source != "member_union":
            continue
        members = tuple(str(m) for m in aggregates.get(code, ()))
        direct_children = tuple(
            str(c)
            for c in cov.index[cov["parent_code"].fillna("").astype(str).eq(code)]
            if str(cov.loc[c, "reference_source"]) in {"parent", "none"}
        )
        avail_members = [
            m for m in members if bool(avail_by_pair.get((code, m), m in selected.index))
        ]
        non_comp = [
            m
            for m in avail_members
            if m in selected.index and not bool(selected.loc[m, "linear_tpm_comparable"])
        ]
        comp_members = [m for m in avail_members if m not in non_comp]
        missing_children = sorted((set(direct_children) | set(members)) - set(avail_members))
        total_n = int(
            sum(
                int(selected.loc[m, "n_reference_samples"])
                for m in avail_members
                if m in selected.index and pd.notna(selected.loc[m, "n_reference_samples"])
            )
        )
        if missing_children or len(avail_members) < len(members):
            status = "partial_member_gap"
        elif non_comp:
            status = "includes_proxy_scale_members"
        else:
            status = "complete_linear_tpm"
        rows.append(
            {
                "cancer_code": code,
                "name": str(rec["name"]),
                "reference_source": ref_source,
                "is_classification_target": bool(rec["is_classification_target"]),
                "member_codes": ";".join(members),
                "available_member_codes": ";".join(avail_members),
                "unbacked_child_subtypes": ";".join(missing_children),
                "non_comparable_member_codes": ";".join(non_comp),
                "n_members": len(members),
                "n_available_members": len(avail_members),
                "n_comparable_members": len(comp_members),
                "total_reference_samples": total_n,
                "coverage_status": status,
            }
        )
    return pd.DataFrame(rows)


@lru_cache(maxsize=1)
def _reference_summary_availability_table() -> pd.DataFrame:
    """Read the compact one-row-per-source summary manifest."""
    table = _expr().get_data(_REFERENCE_SUMMARY_AVAILABILITY_DATASET, copy=False)
    required = {
        "cancer_code",
        "source_cohort",
        "n_reference_genes",
        "n_reference_samples",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{_REFERENCE_SUMMARY_AVAILABILITY_DATASET} is missing columns: {missing}")
    return table


_register_derived_cache(_reference_summary_availability_table.cache_clear)


def _reference_source_union_availability(
    cancer_types: str | Iterable[str] | None,
    modes: list[str],
    *,
    source_kinds: set[str] | None,
    source_cohorts: set[str] | None,
    exclude_microarray_proxy: bool,
) -> pd.DataFrame:
    """Return lightweight availability rows for every loadable summary source."""
    sources = _expr()._filter_reference_summary_sources(
        _expr()._reference_summary_availability_table(),
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )
    if cancer_types is None:
        requests = [
            {"requested_code": code, "cancer_code": code, "request_kind": "default_available"}
            for code in sorted(sources["cancer_code"].astype(str).unique())
        ]
    else:
        requests = _reference_expression_requests(
            cancer_types,
            modes,
            reference_source="summary_rows_all",
            sample_qc="all",
            source_kinds=source_kinds,
            source_cohorts=source_cohorts,
            exclude_microarray_proxy=exclude_microarray_proxy,
            describe_invalid=True,
        )

    rows = []
    for request in requests:
        matching_sources = sources.loc[sources["cancer_code"].astype(str) == request["cancer_code"]]
        source_records = [None] if matching_sources.empty else matching_sources.to_dict("records")
        for mode in modes:
            rows.extend(
                _reference_source_availability_row(request, mode, source)
                for source in source_records
            )
    return _reference_availability_frame(rows)


def _reference_source_availability_row(
    request: dict[str, str], mode: str, source: dict | None
) -> dict:
    """Build one source-specific availability row for a resolved request."""
    available = source is not None
    row = {
        "requested_code": request["requested_code"],
        "cancer_code": request["cancer_code"],
        "request_kind": request["request_kind"],
        "normalization": _REFERENCE_NORMALIZE_LABELS[mode],
        "available": available,
        "missing_reason": "" if available else "no_reference_summary_rows",
        "reference_method": _REFERENCE_SUMMARY_ALL_METHOD,
        "sample_qc": "all",
        "sample_qc_effective": "all",
        "artifact_schema_version": _expr().REFERENCE_EXPRESSION_SCHEMA_VERSION,
        "data_version": DATA_VERSION,
        "source_matrix_version": SOURCE_MATRIX_VERSION,
    }
    if source is None:
        return row

    metadata = _expr()._reference_summary_metadata_from_row(request["cancer_code"], source)
    row.update(
        {
            "source_cohort": metadata.get("source_cohort"),
            "source_project": metadata.get("source_project"),
            "source_version": metadata.get("source_version"),
            "source_pmid": metadata.get("source_pmid"),
            "source_type": metadata.get("source_type"),
            "source_unit": metadata.get("unit"),
            "source_scale_class": metadata.get("source_scale_class"),
            "linear_tpm_comparable": bool(metadata.get("linear_tpm_comparable")),
            "tumor_origin": metadata.get("tumor_origin"),
            "metastasis_site": metadata.get("metastasis_site"),
            "n_reference_genes": metadata.get("n_reference_genes"),
            "n_reference_samples": metadata.get("n_reference_samples"),
            "n_samples": metadata.get("n_reference_samples"),
            "processing_pipeline": source.get("processing_pipeline"),
            "notes": source.get("notes"),
        }
    )
    return row


def _reference_mode_availability(
    code: str,
    mode: str,
    percentile_available: set[str],
    source_matrix_available: set[str],
    summary_available: set[str],
    *,
    sample_qc: str,
    reference_source: str,
    artifact_qc_policies: dict[str, tuple[str, ...]],
    artifact_qc_metadata_present: bool,
) -> tuple[bool, str]:
    if reference_source != "summary_rows_all" and code in _expr()._invalid_selected_source_codes():
        # The artifact may exist, but it is built from an invalid source and withheld.
        return False, "invalid_source_scale"
    if reference_source in {"summary_rows", "summary_rows_all"}:
        if sample_qc == "all":
            return (True, "") if code in summary_available else (False, "no_reference_summary_rows")
        if code not in source_matrix_available:
            return False, "no_source_matrix"
        n_samples = _source_matrix_effective_sample_count(code, sample_qc)
        if n_samples == 0:
            return False, f"no_source_matrix_samples_matching_{sample_qc}_qc"
        return True, ""
    if mode in {"tpm_clean", "tpm_clean_biological", "tpm_clean_log1p"}:
        if code not in percentile_available:
            return False, "no_percentile_artifact"
        if sample_qc == "artifact" or not artifact_qc_metadata_present:
            return True, ""
        if code not in artifact_qc_policies:
            return False, "artifact_build_metadata_missing"
        if artifact_qc_policies[code] != (sample_qc,):
            return False, "artifact_sample_qc_mismatch"
        return True, ""
    if mode == "tpm_raw":
        if code not in source_matrix_available:
            return False, "no_source_matrix"
        if sample_qc != "artifact":
            return True, ""
        if not artifact_qc_metadata_present:
            return False, "artifact_build_metadata_missing"
        effective_sample_qc, missing_reason = _expr()._artifact_effective_sample_qc(
            code, artifact_qc_policies
        )
        if missing_reason:
            return False, missing_reason
        assert effective_sample_qc is not None
        n_samples = _expr()._source_matrix_effective_sample_count(code, effective_sample_qc)
        if n_samples == 0:
            return False, f"no_source_matrix_samples_matching_{effective_sample_qc}_qc"
        return True, ""
    raise AssertionError(f"unhandled reference normalize mode: {mode}")


@cache
def _source_matrix_effective_sample_count(code: str, sample_qc: str) -> int | None:
    """Return selected sample count for cached source matrices, or ``None`` if unknown.

    Availability should not fetch source matrices. When a matrix is cached, use the same
    read-time QC contract as :func:`per_sample_expression` so strict reference summaries
    do not claim availability and then fail after every sample is filtered out.
    """
    mode = _expr()._validate_sample_qc(sample_qc)
    if mode == "all":
        return None
    try:
        qc = _expr().sample_expression_qc(code, auto_fetch=False)
    except (FileNotFoundError, _expr().source_matrices.SourceMatrixError):
        return None
    if qc.empty or "sample_qc_status" not in qc.columns:
        return 0
    statuses = qc["sample_qc_status"].astype(str)
    if mode == "pass":
        return int((statuses == "pass").sum())
    return int(statuses.isin(["pass", "warn"]).sum())


def _reference_expected_method(mode: str, *, sample_qc: str, reference_source: str) -> str:
    if reference_source == "summary_rows_all":
        return _REFERENCE_SUMMARY_ALL_METHOD
    if reference_source == "summary_rows":
        return _REFERENCE_SUMMARY_METHOD if sample_qc == "all" else _REFERENCE_RECOMPUTED_METHOD
    if mode in {"tpm_clean", "tpm_clean_biological"}:
        return "percentile_shard"
    if mode == "tpm_clean_log1p":
        return "percentile_shard_log1p"
    if mode == "tpm_raw":
        return "source_matrix_stats"
    raise AssertionError(f"unhandled reference normalize mode: {mode}")


def _attach_reference_expression_attrs(
    df: pd.DataFrame,
    availability: pd.DataFrame | None,
    *,
    gene_id_style: str,
    gene_universe: str,
    include_gene_universe_flags: bool,
    reference_source: str,
) -> None:
    df.attrs["artifact_schema_version"] = _expr().REFERENCE_EXPRESSION_SCHEMA_VERSION
    df.attrs["data_version"] = DATA_VERSION
    df.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    df.attrs["gene_id_style"] = gene_id_style
    df.attrs["gene_universe"] = gene_universe
    df.attrs["include_gene_universe_flags"] = bool(include_gene_universe_flags)
    df.attrs["reference_source"] = reference_source
    if availability is None:
        return
    missing = availability.loc[~availability["available"]]
    df.attrs["availability"] = availability.to_dict("records")
    df.attrs["missing_requests"] = missing.to_dict("records")


def _reference_normalize_modes(normalize: str | Iterable[str]) -> list[str]:
    if isinstance(normalize, str):
        requested = [normalize]
    else:
        requested = list(normalize)
    modes: list[str] = []
    bad: list[str] = []
    for mode in requested:
        key = str(mode).lower()
        resolved = _REFERENCE_NORMALIZE_ALIASES.get(key)
        if resolved is None:
            bad.append(str(mode))
        elif resolved not in modes:
            modes.append(resolved)
    if bad:
        supported = ", ".join(sorted(_REFERENCE_NORMALIZE_ALIASES))
        raise ValueError(f"unsupported reference normalize mode(s): {bad}; supported: {supported}")
    return modes


def _validate_reference_source(reference_source: str) -> str:
    source = str(reference_source).lower()
    if source not in _REFERENCE_SOURCES:
        supported = ", ".join(sorted(_REFERENCE_SOURCES))
        raise ValueError(f"reference_source must be one of {supported}")
    return source
