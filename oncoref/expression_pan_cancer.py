# Licensed under the Apache License, Version 2.0.
"""Wide HPA normal + TCGA tumor pan-cancer reference expression loader (#756)."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

from .gene_ids import unversioned
from .load_dataset import _register_derived_cache
from .normalization import clean_tpm, percentile_rank, tpm_to_housekeeping_normalized
from .version import DATA_VERSION, SOURCE_MATRIX_VERSION


def _expr():
    from . import expression

    return expression


def pan_cancer_expression(
    genes: str | Iterable[str] | None = None,
    *,
    normalize: str | Iterable[str] | None = "tpm_clean",
    to_tpm: bool | None = None,
    column_style: str | None = None,
) -> pd.DataFrame:
    """Wide pan-cancer reference: each gene's expression across **50 HPA normal
    tissues** and **33 TCGA tumor cohorts**, tumor and normal side by side in one
    frame — the combined companion to the per-cohort accessors above.

    Columns are entity-first: ``<tissue>_nTPM_raw`` for HPA normal tissues,
    ``<CODE>_FPKM_raw`` for TCGA source/provenance values, and ``<CODE>_TPM_raw``
    for deterministic FPKM→TPM tumor companions.
    The default ``normalize="tpm_clean"`` (alias ``"clean_tpm"``) appends
    clean TPM companions named ``<entity>_<measure>_clean``.
    ``normalize="housekeeping"`` / ``"hk"`` appends ``*_hk`` columns, and
    ``"percentile"`` appends ``*_percentile`` columns. ``"tpm_log1p"`` and
    ``"tpm_clean_log1p"`` append natural-log companions over the raw and clean
    TPM/nTPM values, respectively. Pass ``normalize=None`` or ``"tpm"`` for the
    raw TPM/nTPM companions only.

    ``column_style="pirlygenes"`` returns the same data with pirlygenes-style
    unsuffixed analysis names (for example ``LUAD_FPKM``, ``LUAD_TPM``, and
    ``liver_nTPM``). The legacy ``to_tpm`` keyword is accepted as a compatibility
    alias; when used without an explicit ``column_style`` it selects the
    pirlygenes-style output and, for the default normalization, maps to
    ``normalize="tpm"``.

    ``genes`` filters to the given Ensembl gene ids (version-insensitive) or
    symbols; ``None`` returns the full reference matrix. FPKM→TPM conversion and
    clean TPM normalization run over **all reference rows** before filtering,
    so a filtered slice retains that scaling and need not sum to 1e6. The source
    reference is not a whole transcriptome: clean compartment budgets are
    distributed over its measured rows, independently for each column. An
    unavailable compartment warns and leaves its budget unfilled."""
    if to_tpm is not None:
        if normalize == "tpm_clean":
            normalize = "tpm" if to_tpm else None
        if column_style is None:
            column_style = "pirlygenes"
    column_style = _pan_cancer_column_style(column_style)

    df = _expr().get_data("pan-cancer-expression")
    # Dense canonical space: sum alt-haplotype/patch copies into their primary gene and
    # relabel retired ids (oncoref#135 item 6). The cohort/tissue columns are per-gene
    # abundances (nTPM / FPKM) — linear-additive — so summing the rows of a fragmented
    # gene is exact, and it precedes the FPKM->TPM rescale (whose per-cohort 1e6 total is
    # conserved under row-summing, so every other gene's conversion is unchanged).
    id_cols = ["Ensembl_Gene_ID", "Symbol"]
    value_cols = [c for c in df.columns if c not in id_cols]
    df = _expr()._canonicalize_gene_rows(df, sample_cols=value_cols)
    out = df[id_cols].copy()

    raw_cols: list[str] = []
    raw_parts: list[pd.DataFrame] = []
    ntpm_cols = [c for c in df.columns if c.startswith("nTPM_")]
    ntpm_data: dict[str, pd.Series] = {}
    for col in ntpm_cols:
        entity = col[len("nTPM_") :]
        target = f"{entity}_nTPM_raw"
        ntpm_data[target] = pd.to_numeric(df[col], errors="coerce")
        raw_cols.append(target)
    if ntpm_data:
        raw_parts.append(pd.DataFrame(ntpm_data, index=df.index))

    fpkm_cols = [c for c in df.columns if c.startswith("FPKM_")]
    if fpkm_cols:
        from .normalization import fpkm_to_tpm

        converted, _ = fpkm_to_tpm(df[id_cols + fpkm_cols], value_cols=fpkm_cols)
        fpkm_data: dict[str, pd.Series] = {}
        for col in fpkm_cols:
            entity = col[len("FPKM_") :]
            provenance_col = f"{entity}_FPKM_raw"
            fpkm_data[provenance_col] = pd.to_numeric(df[col], errors="coerce")
            target = f"{entity}_TPM_raw"
            fpkm_data[target] = converted[col]
            raw_cols.append(target)
        raw_parts.append(pd.DataFrame(fpkm_data, index=df.index))

    tpm_cols = [c for c in df.columns if c.startswith("TPM_")]
    tpm_data: dict[str, pd.Series] = {}
    for col in tpm_cols:
        entity = col[len("TPM_") :]
        target = f"{entity}_TPM_raw"
        tpm_data[target] = pd.to_numeric(df[col], errors="coerce")
        raw_cols.append(target)
    if tpm_data:
        raw_parts.append(pd.DataFrame(tpm_data, index=df.index))
    if raw_parts:
        out = pd.concat([out, *raw_parts], axis=1)

    out, raw_cols, computed_aggregate_cols = _add_computed_pan_cancer_raw_columns(out, raw_cols)

    modes = _pan_cancer_normalize_modes(normalize)
    clean_cols: list[str] = []
    if modes & {"tpm_clean", "housekeeping", "hk", "percentile", "tpm_clean_log1p"}:
        clean = clean_tpm(out[raw_cols], gene_table=out[id_cols])
        clean_data: dict[str, pd.Series] = {}
        for col in raw_cols:
            target = col[: -len("_raw")] + "_clean"
            clean_data[target] = clean[col]
            clean_cols.append(target)
        out = pd.concat([out, pd.DataFrame(clean_data, index=out.index)], axis=1)

    if modes & {"housekeeping", "hk"}:
        hk_input_cols = clean_cols or raw_cols
        hk_input = out[id_cols + hk_input_cols].copy()
        panel_ids = _expr()._clean_tpm_housekeeping_panel_ids()
        hk, _ = tpm_to_housekeeping_normalized(
            hk_input,
            value_cols=hk_input_cols,
            panel_ids=panel_ids,
            panel_name="clean_tpm_biological_housekeeping",
            min_panel_detected=_expr()._min_housekeeping_detected(panel_ids),
            drop_zero_panel_values=True,
            warn_on_unreliable=True,
        )
        hk_data: dict[str, pd.Series] = {}
        for col in hk_input_cols:
            target = col.rsplit("_", 1)[0] + "_hk"
            hk_data[target] = hk[col]
        out = pd.concat([out, pd.DataFrame(hk_data, index=out.index)], axis=1)

    if "percentile" in modes:
        pct_input_cols = clean_cols or raw_cols
        pct = percentile_rank(out[id_cols + pct_input_cols], value_cols=pct_input_cols)
        pct_data: dict[str, pd.Series] = {}
        for col in pct_input_cols:
            target = col.rsplit("_", 1)[0] + "_percentile"
            pct_data[target] = pct[col]
        out = pd.concat([out, pd.DataFrame(pct_data, index=out.index)], axis=1)

    if "tpm_log1p" in modes:
        log_data: dict[str, np.ndarray] = {}
        for col in raw_cols:
            target = col + "_log1p"
            log_data[target] = np.log1p(out[col].to_numpy(dtype=float))
        out = pd.concat([out, pd.DataFrame(log_data, index=out.index)], axis=1)

    if "tpm_clean_log1p" in modes:
        clean_log_data: dict[str, np.ndarray] = {}
        for col in clean_cols:
            target = col + "_log1p"
            clean_log_data[target] = np.log1p(out[col].to_numpy(dtype=float))
        out = pd.concat([out, pd.DataFrame(clean_log_data, index=out.index)], axis=1)

    if genes is not None:
        wanted = {genes} if isinstance(genes, str) else set(genes)
        wanted = {str(g) for g in wanted}
        wanted_unversioned = {unversioned(g) for g in wanted}
        ids = out["Ensembl_Gene_ID"].astype(str)
        mask = (
            ids.isin(wanted)
            | ids.map(unversioned).isin(wanted_unversioned)
            | out["Symbol"].astype(str).isin(wanted)
        )
        out = out[mask].reset_index(drop=True)
    out = _format_pan_cancer_columns(out, column_style=column_style)
    out.attrs["oncoref"] = {
        "dataset": "pan-cancer-expression",
        "data_version": DATA_VERSION,
        "source_matrix_version": SOURCE_MATRIX_VERSION,
        "column_style": column_style,
        "computed_aggregate_columns": tuple(computed_aggregate_cols),
    }
    return out


_PAN_CANCER_COMPUTED_AGGREGATES = ("NET", "CRC", "NSCLC", "BTC", "SGC")


def _add_computed_pan_cancer_raw_columns(
    out: pd.DataFrame, raw_cols: list[str]
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Append TPM columns for grouping codes whose references are computed from members."""
    if out.empty:
        return out, raw_cols, []
    added: list[str] = []
    raw_cols = list(raw_cols)
    gene_key = out["Ensembl_Gene_ID"].astype(str).map(unversioned)
    frames: list[pd.DataFrame] = []
    for aggregate_code in _PAN_CANCER_COMPUTED_AGGREGATES:
        target = f"{aggregate_code}_TPM_raw"
        if target in out.columns:
            continue
        members = _expr()._computed_expression_reference_members(aggregate_code)
        if not members:
            continue
        values = _pooled_reference_expression_series_for_pan_cancer(members)
        if values.empty:
            continue
        frames.append(pd.DataFrame({target: gene_key.map(values)}, index=out.index))
        raw_cols.append(target)
        added.append(target)
    if frames:
        out = pd.concat([out, *frames], axis=1)
    return out, raw_cols, added


_PAN_CANCER_POOLED_SERIES_CACHE: dict[tuple[str, ...], pd.Series] = {}


def _clear_pan_cancer_pooled_series_cache() -> None:
    _PAN_CANCER_POOLED_SERIES_CACHE.clear()


_register_derived_cache(_clear_pan_cancer_pooled_series_cache)


def _pooled_reference_expression_series_for_pan_cancer(member_codes: tuple[str, ...]) -> pd.Series:
    """n-sample-weighted raw TPM medians for a computed pan-cancer aggregate."""
    if member_codes not in _PAN_CANCER_POOLED_SERIES_CACHE:
        _PAN_CANCER_POOLED_SERIES_CACHE[member_codes] = _compute_pooled_reference_expression_series(
            member_codes
        )
    return _PAN_CANCER_POOLED_SERIES_CACHE[member_codes]


def _compute_pooled_reference_expression_series(member_codes: tuple[str, ...]) -> pd.Series:
    try:
        rows = _selected_reference_summary_rows_for_pan_cancer(member_codes)
    except (FileNotFoundError, KeyError, TypeError):
        return pd.Series(dtype=float)
    if rows.empty:
        return pd.Series(dtype=float)
    work = rows.copy()
    work["_gene_key"] = work["Ensembl_Gene_ID"].astype(str).map(unversioned)
    work["_expression"] = pd.to_numeric(work["TPM_median"], errors="coerce")
    weights = pd.to_numeric(work.get("n_samples"), errors="coerce")
    if not isinstance(weights, pd.Series):
        weights = pd.Series(1.0, index=work.index)
    work["_weight"] = weights.reindex(work.index).fillna(0.0)
    rows: dict[str, float] = {}
    for gene_id, group in work.groupby("_gene_key", sort=False):
        valid = group["_expression"].notna() & (group["_weight"] > 0)
        if not valid.any():
            continue
        rows[str(gene_id)] = float(
            np.average(group.loc[valid, "_expression"], weights=group.loc[valid, "_weight"])
        )
    return pd.Series(rows, dtype=float)


def _selected_reference_summary_rows_for_pan_cancer(member_codes: tuple[str, ...]) -> pd.DataFrame:
    source_table = _expr()._reference_summary_source_table()
    if source_table.empty:
        return pd.DataFrame()
    wanted = {str(code) for code in member_codes}
    selected = source_table.loc[
        source_table["cancer_code"].astype(str).isin(wanted) & source_table["selected"]
    ].copy()
    if selected.empty:
        return pd.DataFrame()
    parts = [
        _expr()._reference_summary_rows_for_sources(
            str(code),
            sources,
        )
        for code, sources in selected.groupby("cancer_code", sort=False, observed=True)
    ]
    parts = [part for part in parts if not part.empty]
    if not parts:
        return pd.DataFrame(columns=["Ensembl_Gene_ID", "TPM_median", "n_samples"])
    needed_columns = ["Ensembl_Gene_ID", "TPM_median", "n_samples"]
    return pd.concat([part[needed_columns] for part in parts], ignore_index=True)


def _pan_cancer_column_style(column_style: str | None) -> str:
    if column_style is None:
        return "oncoref"
    style = str(column_style).lower()
    aliases = {"default": "oncoref", "entity_first": "oncoref", "legacy": "pirlygenes"}
    style = aliases.get(style, style)
    if style not in {"oncoref", "pirlygenes"}:
        raise ValueError("column_style must be None, 'oncoref', or 'pirlygenes'")
    return style


def _format_pan_cancer_columns(df: pd.DataFrame, *, column_style: str) -> pd.DataFrame:
    if column_style == "oncoref":
        return df
    rename: dict[str, str] = {}
    for col in df.columns:
        if col.endswith("_raw_log1p"):
            rename[col] = col.replace("_raw_log1p", "_log1p")
        elif col.endswith(("_nTPM_raw", "_FPKM_raw", "_TPM_raw")):
            rename[col] = col[: -len("_raw")]
    return df.rename(columns=rename)


def _pan_cancer_normalize_modes(normalize: str | Iterable[str] | None) -> set[str]:
    if normalize is None:
        return set()
    if isinstance(normalize, str):
        modes = {normalize.lower()}
    else:
        modes = {str(mode).lower() for mode in normalize}
    aliases = {"clean_tpm": "tpm_clean", "housekeeping": "housekeeping", "hk": "hk"}
    modes = {aliases.get(mode, mode) for mode in modes}
    allowed = {
        "tpm",
        "raw",
        "tpm_log1p",
        "tpm_clean",
        "housekeeping",
        "hk",
        "percentile",
        "tpm_clean_log1p",
    }
    unknown = modes - allowed
    if unknown:
        raise ValueError(
            "normalize must be None, 'tpm', 'tpm_clean'/'clean_tpm', "
            "'housekeeping', 'hk', 'percentile', 'tpm_log1p', "
            "or 'tpm_clean_log1p'"
        )
    return modes - {"tpm", "raw"}
