# Licensed under the Apache License, Version 2.0.
"""Gene-universe delta classification and reporting helpers for expression artifacts (#756)."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from .gene_ids import gene_biotype
from .gene_qc import TECHNICAL_RNA_GROUPS, classify_gene_qc


def _expr():
    from . import expression

    return expression


_ARTIFACT_TECHNICAL_EXTRA_STATUSES = frozenset(
    {
        "technical_or_noncoding_extra",
        "y_linked_extra",
        "immune_receptor_segment_extra",
    }
)
_ARTIFACT_MISSING_BIOLOGICAL_STATUSES = frozenset(
    {
        "canonical_replacement_absent_from_output",
        "canonical_row_absent_from_oncoref_output",
        "unresolved_missing_oncoref_row",
    }
)
_ARTIFACT_CANONICALIZATION_STATUSES = frozenset(
    {
        "intentional_canonicalization",
        "remapped_to_oncoref",
        "sequence_identical_remapped_to_oncoref",
    }
)
_ARTIFACT_NON_SIGNAL_EXTRA_STATUSES = frozenset({"non_signal_oncoref_extra"})
_ARTIFACT_BIOLOGICAL_EXTRA_STATUSES = frozenset({"biological_oncoref_extra"})
_ARTIFACT_GENE_UNIVERSE_MODES = ("artifact", "tumor_signal", "pirlygenes")
_ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS = [
    "artifact_row_class",
    "is_filterable_extra",
    "is_technical_extra",
    "is_missing_biological",
    "recommended_consumer_action",
]
_ARTIFACT_NON_SIGNAL_EXTRA_BIOTYPES = frozenset(
    {
        "misc_RNA",
        "snRNA",
        "snoRNA",
        "scaRNA",
        "rRNA",
        "Mt_rRNA",
        "Mt_tRNA",
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_V_gene",
        "TR_C_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_V_gene",
    }
)
_ARTIFACT_BIOLOGICAL_EXTRA_BIOTYPES = frozenset({"protein_coding", "lncRNA"})
_ARTIFACT_FILTERABLE_EXTRA_CLASSES = frozenset({"technical_extra", "non_signal_extra"})


def _delta_nonempty(value: object) -> str | None:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return None
    return text


def _delta_gene_biotype(row) -> str | None:
    gene_id = _delta_nonempty(row.oncoref_ensembl_gene_id) or _delta_nonempty(
        row.legacy_ensembl_gene_id
    )
    return gene_biotype(gene_id) if gene_id else None


def _artifact_delta_class(row) -> str:
    status = str(row.status)
    if status in _ARTIFACT_TECHNICAL_EXTRA_STATUSES:
        return "technical_extra"
    if status in _ARTIFACT_NON_SIGNAL_EXTRA_STATUSES:
        return "non_signal_extra"
    if status in _ARTIFACT_BIOLOGICAL_EXTRA_STATUSES:
        return "biological_extra"
    if status == "unresolved_oncoref_extra":
        qc = classify_gene_qc(
            _delta_nonempty(row.symbol),
            ensembl_id=_delta_nonempty(row.oncoref_ensembl_gene_id),
        )
        biotype = _delta_nonempty(_delta_gene_biotype(row))
        if qc.group in TECHNICAL_RNA_GROUPS or qc.group in {
            "small_ncrna",
            "immune_receptor",
            "ribosomal_protein_pseudogene",
        }:
            return "technical_extra"
        if biotype in _ARTIFACT_NON_SIGNAL_EXTRA_BIOTYPES or (
            biotype is not None and "pseudogene" in biotype
        ):
            return "non_signal_extra"
        if biotype in _ARTIFACT_BIOLOGICAL_EXTRA_BIOTYPES:
            return "biological_extra"
        return "unresolved_extra"
    if status in _ARTIFACT_MISSING_BIOLOGICAL_STATUSES:
        return "missing_biological"
    if status in _ARTIFACT_CANONICALIZATION_STATUSES:
        return "canonicalized"
    return "unclassified"


def _artifact_delta_consumer_action(row) -> str:
    status = str(row.status)
    if status in _ARTIFACT_TECHNICAL_EXTRA_STATUSES:
        return "filter_from_signal_views"
    if status in _ARTIFACT_NON_SIGNAL_EXTRA_STATUSES:
        return "filter_from_signal_views"
    if status in _ARTIFACT_BIOLOGICAL_EXTRA_STATUSES:
        return "keep_oncoref_biological_row"
    if status == "unresolved_oncoref_extra":
        row_class = _artifact_delta_class(row)
        if row_class in _ARTIFACT_FILTERABLE_EXTRA_CLASSES:
            return "filter_from_signal_views"
        if row_class == "biological_extra":
            return "keep_oncoref_biological_row"
        return "audit_before_filtering"
    if status in _ARTIFACT_MISSING_BIOLOGICAL_STATUSES:
        return "restore_or_remap_in_next_bundle"
    if status in _ARTIFACT_CANONICALIZATION_STATUSES:
        return "accept_canonical_mapping"
    return "audit_status"


def _annotate_expression_artifact_gene_universe_deltas(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    statuses = out["status"].astype(str)
    out["gene_biotype"] = [_delta_gene_biotype(row) for row in out.itertuples(index=False)]
    rows = list(out.itertuples(index=False))
    out["artifact_row_class"] = [_artifact_delta_class(row) for row in rows]
    out["is_filterable_extra"] = out["artifact_row_class"].isin(_ARTIFACT_FILTERABLE_EXTRA_CLASSES)
    out["is_technical_extra"] = out["artifact_row_class"].eq("technical_extra")
    out["is_missing_biological"] = statuses.isin(_ARTIFACT_MISSING_BIOLOGICAL_STATUSES)
    out["recommended_consumer_action"] = [_artifact_delta_consumer_action(row) for row in rows]
    return out


def expression_artifact_gene_universe_deltas(
    *,
    product: str | None = None,
    cancer_type: str | None = None,
    delta_kind: str | None = None,
    status: str | None = None,
) -> pd.DataFrame:
    """Known row-universe deltas between pirlygenes and oncoref expression artifacts.

    This is a provenance/audit table, not a value-transforming compatibility shim. It
    records the known remapped, missing, and extra rows from the pirlygenes 5.23.2 vs
    oncoref 5.23.3 parity run tracked in #191/#193/#278 so downstream migration code can
    distinguish intentional canonicalization, biological rows missing from oncoref
    artifacts, strict technical extras, broader filterable non-signal extras, and
    biological oncoref-only rows that should stay visible until the heavy artifacts are
    regenerated. The explicit ``non_signal_oncoref_extra`` and
    ``biological_oncoref_extra`` statuses encode those policy decisions directly.

    Optional filters are exact for ``product``, ``delta_kind``, and ``status``.
    ``cancer_type`` resolves aliases and matches semicolon-separated cohort lists in the
    table (for deltas shared by PRAD/COAD_MSI/READ_MSI, for example).
    """
    df = _annotate_expression_artifact_gene_universe_deltas(
        _expr().get_data("expression-artifact-gene-universe-deltas")
    )
    if product is not None:
        df = df[df["product"].astype(str) == str(product)]
    if delta_kind is not None:
        df = df[df["delta_kind"].astype(str) == str(delta_kind)]
    if status is not None:
        df = df[df["status"].astype(str) == str(status)]
    if cancer_type is not None:
        code = _expr().resolve_cancer_type(cancer_type)

        def _matches(cell: object) -> bool:
            codes = [c for c in str(cell or "").split(";") if c]
            return code in codes

        df = df[df["cancer_code"].map(_matches)]
    df = df.reset_index(drop=True)
    df.attrs["comparison"] = "pirlygenes_5.23.2_vs_oncoref_5.23.3"
    df.attrs["issues"] = ["#191", "#193", "#278"]
    return df


def expression_artifact_technical_extra_gene_ids(
    *,
    product: str | None = None,
    cancer_type: str | None = None,
) -> list[str]:
    """Oncoref artifact row IDs that are known technical extras for a request.

    This is the machine-readable subset of
    :func:`expression_artifact_gene_universe_deltas` where ``is_technical_extra`` is
    true: small/noncoding RNA, Y-linked rows, and immune-receptor segments that appear
    in oncoref artifacts but are not expected in pirlygenes-compatible tumor-signal
    views. It is an audit/filter list over the current bundle; it does not mutate
    expression values or claim those rows have been removed from shipped shards.
    """
    df = expression_artifact_gene_universe_deltas(
        product=product,
        cancer_type=cancer_type,
        delta_kind="oncoref_only",
    )
    df = df[df["is_technical_extra"]]
    ids = {gid for gid in df["oncoref_ensembl_gene_id"].map(_delta_nonempty) if gid is not None}
    return sorted(ids)


def _validate_artifact_gene_universe(gene_universe: str) -> str:
    mode = str(gene_universe)
    if mode not in _ARTIFACT_GENE_UNIVERSE_MODES:
        allowed = "', '".join(_ARTIFACT_GENE_UNIVERSE_MODES)
        raise ValueError(f"gene_universe must be one of '{allowed}'")
    return mode


def _artifact_row_key(series: pd.Series) -> pd.Series:
    return series.astype(str).str.split(".").str[0]


def _row_universe_annotation_frame(product: str, cancer_codes: Iterable[str]) -> pd.DataFrame:
    """Per-output-row annotations for known oncoref-side artifact deltas."""
    source_product = _delta_source_product(product)
    df = expression_artifact_gene_universe_deltas(
        product=source_product,
        delta_kind="oncoref_only",
    )
    df = _filter_delta_rows_by_codes(df, [str(code) for code in cancer_codes])
    if df.empty:
        return pd.DataFrame(
            columns=[
                "Ensembl_Gene_ID",
                *_ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS,
            ]
        )
    work = df.copy()
    work["Ensembl_Gene_ID"] = work["oncoref_ensembl_gene_id"].map(_delta_nonempty)
    work = work.dropna(subset=["Ensembl_Gene_ID"])

    rows: list[dict[str, object]] = []
    class_priority = (
        "technical_extra",
        "non_signal_extra",
        "unresolved_extra",
        "biological_extra",
        "canonicalized",
        "missing_biological",
        "unclassified",
    )
    action_by_class = {
        "technical_extra": "filter_from_signal_views",
        "non_signal_extra": "filter_from_signal_views",
        "unresolved_extra": "audit_before_filtering",
        "biological_extra": "keep_oncoref_biological_row",
        "canonicalized": "accept_canonical_mapping",
        "missing_biological": "restore_or_remap_in_next_bundle",
        "unclassified": "audit_status",
    }
    for gene_id, group in work.groupby("Ensembl_Gene_ID", sort=False):
        classes = set(group["artifact_row_class"].astype(str))
        row_class = next((name for name in class_priority if name in classes), "artifact")
        rows.append(
            {
                "Ensembl_Gene_ID": str(gene_id),
                "artifact_row_class": row_class,
                "is_technical_extra": row_class == "technical_extra",
                "is_filterable_extra": row_class in _ARTIFACT_FILTERABLE_EXTRA_CLASSES,
                "is_missing_biological": row_class == "missing_biological",
                "recommended_consumer_action": action_by_class.get(row_class, "keep"),
            }
        )
    return pd.DataFrame(rows)


def _pirlygenes_compatible_row_ids(product: str, cancer_codes: Iterable[str]) -> set[str]:
    """Oncoref artifact rows with a documented pirlygenes legacy counterpart."""
    source_product = _delta_source_product(product)
    id_map, _ = _expr()._artifact_legacy_gene_id_map(
        product=source_product,
        cancer_codes=[str(code) for code in cancer_codes],
    )
    return set(id_map)


def _apply_artifact_gene_universe(
    df: pd.DataFrame,
    *,
    product: str,
    cancer_codes: Iterable[str],
    gene_universe: str,
    include_gene_universe_flags: bool,
) -> pd.DataFrame:
    """Filter/annotate known expression-artifact row-universe deltas.

    ``gene_universe="artifact"`` preserves the exact shipped row set.
    ``"tumor_signal"`` drops rows explicitly audited as filterable extras
    (strict technical rows plus biotype-resolved non-signal extras such as
    pseudogene, small-RNA, and immune-receptor rows). It keeps biological
    oncoref-only rows and never synthesizes missing biological rows.
    ``"pirlygenes"`` is a stricter migration-compatibility view: it keeps the
    artifact rows and documented canonical/remap targets needed to reproduce the
    pirlygenes row universe, but drops oncoref-only audited extras that have no
    pirlygenes counterpart. It still never synthesizes rows or values.
    """
    mode = _validate_artifact_gene_universe(gene_universe)
    if "Ensembl_Gene_ID" not in df.columns:
        return df
    attrs = dict(df.attrs)
    out = df.copy()
    annotations = _row_universe_annotation_frame(product, cancer_codes)
    ann = annotations.set_index("Ensembl_Gene_ID") if not annotations.empty else annotations

    if mode in {"tumor_signal", "pirlygenes"} and not annotations.empty:
        drop_mask = annotations["is_filterable_extra"]
        if mode == "pirlygenes":
            compatible_ids = _pirlygenes_compatible_row_ids(product, cancer_codes)
            drop_mask = ~annotations["Ensembl_Gene_ID"].astype(str).isin(compatible_ids)
        drop_ids = set(annotations.loc[drop_mask, "Ensembl_Gene_ID"].astype(str))
        if drop_ids:
            keys = _artifact_row_key(out["Ensembl_Gene_ID"])
            out = out.loc[~keys.isin(drop_ids)].reset_index(drop=True)

    if include_gene_universe_flags:
        keys = _artifact_row_key(out["Ensembl_Gene_ID"])
        if annotations.empty:
            out["artifact_row_class"] = "artifact"
            out["is_filterable_extra"] = False
            out["is_technical_extra"] = False
            out["is_missing_biological"] = False
            out["recommended_consumer_action"] = "keep"
        else:
            out["artifact_row_class"] = (
                keys.map(ann["artifact_row_class"]).fillna("artifact").to_numpy()
            )
            out["is_technical_extra"] = keys.map(ann["is_technical_extra"]).eq(True).to_numpy()
            out["is_filterable_extra"] = keys.map(ann["is_filterable_extra"]).eq(True).to_numpy()
            out["is_missing_biological"] = (
                keys.map(ann["is_missing_biological"]).eq(True).to_numpy()
            )
            out["recommended_consumer_action"] = (
                keys.map(ann["recommended_consumer_action"]).fillna("keep").to_numpy()
            )

    out.attrs.update(attrs)
    out.attrs["gene_universe"] = mode
    out.attrs["include_gene_universe_flags"] = bool(include_gene_universe_flags)
    return out


def expression_artifact_gene_universe_delta_summary() -> pd.DataFrame:
    """Counts of known expression-artifact row-universe deltas by product/status."""
    df = expression_artifact_gene_universe_deltas()
    if df.empty:
        return pd.DataFrame(
            columns=[
                "product",
                "cancer_code",
                "delta_kind",
                "status",
                "artifact_row_class",
                "is_technical_extra",
                "is_missing_biological",
                "recommended_consumer_action",
                "n",
            ]
        )
    return (
        df.groupby(
            [
                "product",
                "cancer_code",
                "delta_kind",
                "status",
                "artifact_row_class",
                "is_technical_extra",
                "is_missing_biological",
                "recommended_consumer_action",
            ],
            dropna=False,
        )
        .size()
        .rename("n")
        .reset_index()
    )


_DELTA_REPORT_COLUMNS = [
    "accessor",
    "product",
    "cancer_code",
    "delta_kind",
    "status",
    "artifact_row_class",
    "is_filterable_extra",
    "is_technical_extra",
    "is_missing_biological",
    "recommended_consumer_action",
    "n",
    "legacy_ensembl_gene_ids",
    "oncoref_ensembl_gene_ids",
    "gene_biotypes",
    "symbols",
    "issues",
]


def _delta_source_product(product: str) -> str:
    aliases = {
        "cancer_reference_expression": "cohort_gene_percentiles",
        "cohort_gene_percentiles": "cohort_gene_percentiles",
        "representative_cohort_samples": "representative_cohort_samples",
    }
    try:
        return aliases[str(product)]
    except KeyError as e:
        allowed = "', '".join(sorted(aliases))
        raise ValueError(f"product must be one of '{allowed}'") from e


def _resolve_delta_report_codes(cancer_types: str | Iterable[str] | None) -> list[str] | None:
    if cancer_types is None:
        return None
    if isinstance(cancer_types, str):
        requested = [cancer_types]
    else:
        requested = list(cancer_types)
    if not requested:
        return []
    return _expr()._resolve_cancer_types(requested, expand_aggregates=True) or []


def _filter_delta_rows_by_codes(df: pd.DataFrame, codes: list[str] | None) -> pd.DataFrame:
    if codes is None:
        return df
    if not codes:
        return df.iloc[0:0].copy()
    code_set = set(codes)

    def _matches(cell: object) -> bool:
        cell_codes = {c for c in str(cell or "").split(";") if c}
        return bool(cell_codes & code_set)

    return df[df["cancer_code"].map(_matches)].copy()


def _join_unique(values: pd.Series) -> str:
    out = sorted({str(v) for v in values if pd.notna(v) and str(v)})
    return ";".join(out)


def _expression_artifact_gene_universe_delta_report_for_codes(
    product: str,
    codes: list[str] | None,
) -> pd.DataFrame:
    source_product = _delta_source_product(product)
    df = expression_artifact_gene_universe_deltas(product=source_product)
    df = _filter_delta_rows_by_codes(df, codes)
    if df.empty:
        out = pd.DataFrame(columns=_DELTA_REPORT_COLUMNS)
    else:
        grouped = df.groupby(
            [
                "product",
                "cancer_code",
                "delta_kind",
                "status",
                "artifact_row_class",
                "is_technical_extra",
                "is_missing_biological",
                "recommended_consumer_action",
            ],
            dropna=False,
        )
        out = grouped.agg(
            n=("status", "size"),
            is_filterable_extra=("is_filterable_extra", "first"),
            legacy_ensembl_gene_ids=("legacy_ensembl_gene_id", _join_unique),
            oncoref_ensembl_gene_ids=("oncoref_ensembl_gene_id", _join_unique),
            gene_biotypes=("gene_biotype", _join_unique),
            symbols=("symbol", _join_unique),
            issues=("issue", _join_unique),
        ).reset_index()
        out.insert(0, "accessor", str(product))
        out = out[_DELTA_REPORT_COLUMNS]
    out.attrs["comparison"] = "pirlygenes_5.23.2_vs_oncoref_5.23.3"
    out.attrs["issues"] = ["#191", "#193", "#278"]
    out.attrs["requested_cancer_codes"] = None if codes is None else list(codes)
    return out


def expression_artifact_gene_universe_delta_report(
    product: str,
    cancer_types: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Request-scoped row-universe delta report for expression accessors.

    This is a compact companion to
    :func:`expression_artifact_gene_universe_deltas`. It filters the known
    pirlygenes/oncoref row-set audit to the accessor/product and requested cancer
    codes, then returns counts by ``delta_kind`` and ``status`` with the affected
    legacy/oncoref identifiers summarized. It does not filter, synthesize, or alter
    expression values.

    ``product="cancer_reference_expression"`` reports against the underlying
    ``cohort_gene_percentiles`` artifact used by clean-TPM reference-expression
    modes. Empty ``cancer_types`` is preserved as an empty report.
    """
    codes = _resolve_delta_report_codes(cancer_types)
    return _expression_artifact_gene_universe_delta_report_for_codes(product, codes)


def _attach_gene_universe_delta_attrs(
    df: pd.DataFrame,
    *,
    product: str,
    cancer_codes: Iterable[str],
) -> None:
    report = _expression_artifact_gene_universe_delta_report_for_codes(
        product, [str(code) for code in cancer_codes]
    )
    records = report.to_dict("records")
    df.attrs["gene_universe_delta_summary"] = records
    df.attrs["gene_universe_delta_n"] = int(report["n"].sum()) if "n" in report else 0
    df.attrs["gene_universe_delta_comparison"] = report.attrs["comparison"]
    df.attrs["gene_universe_delta_issues"] = report.attrs["issues"]
