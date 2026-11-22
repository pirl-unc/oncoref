# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Per-patient antigen coverage over the full per-sample matrices.

These answer questions a per-cohort summary cannot, because they need the *joint*
per-sample matrix — which patient expresses which antigen, not just per-gene
prevalence:

  - :func:`cta_patient_fractions` — per gene, the fraction of a cohort's patients
    expressing it above a TPM threshold;
  - :func:`addressable_fraction` — the fraction of patients expressing **≥1** of a
    gene panel (the union — "what share of this cohort a CTA-directed therapy could
    address"), which is NOT the per-gene fractions summed;
  - :func:`greedy_coverage` — a minimal antigen panel by greedy set cover: at each
    step add the gene covering the most still-uncovered patients.

All operate on :func:`oncoref.expression.per_sample_expression` (clean TPM), so
they need the cohort's per-sample matrix fetched (see :mod:`oncoref.source_matrices`).
The default gene panel is the current ``target_selection`` CTA set
(:func:`oncoref.antigen_sets.gene_ids`).
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import antigen_sets
from .cancer_genes import cancer_key_genes_df, cancer_type_genes_df
from .cancer_types import format_cancer_code_label, resolve_cancer_type
from .expression import (
    cancer_reference_expression,
    cancer_reference_expression_availability,
    per_sample_expression,
)
from .expression_engine import ID_COLUMNS as _ANTIGEN_ID_COLS
from .figure_export import save_figure
from .gene_families import gene_families, gene_family_ids, housekeeping_gene_ids
from .gene_ids import (
    canonical_gene_space,
    canonical_gene_symbol,
    gene_id_column,
    is_ensembl_gene_id,
    resolve_symbol,
    symbol_column,
    symbol_synonyms,
)
from .incidence import burden_category, cancer_burden
from .normalization import percentile_rank


def cta_gene_ids() -> frozenset[str]:
    """Canonical member Ensembl gene IDs for current ``target_selection`` CTAs."""
    return frozenset(_clean_id(g) for g in antigen_sets.gene_ids("CTA"))


def cta_gene_id_to_name() -> dict[str, str]:
    """``{unversioned CTA gene id -> canonical display symbol}`` for current CTAs."""
    from . import target_selection

    group_labels: dict[str, str] = {}
    ctas = target_selection.results(category="CTA")
    for row in ctas[["protein_group_label", "member_gene_ids"]].itertuples(index=False):
        for gid in json.loads(row.member_gene_ids):
            group_labels[_clean_id(gid)] = str(row.protein_group_label)
    out: dict[str, str] = {}
    for gid in sorted(cta_gene_ids()):
        clean = _clean_id(gid)
        out[clean] = canonical_gene_symbol(clean) or group_labels.get(clean, clean)
    return out


#: Default clean-TPM threshold for "expressed in a patient". 10 TPM is a common
#: cut for a confidently-expressed transcript.
DEFAULT_EXPRESSED_TPM: float = 10.0

DEFAULT_SOURCE: str = "all"
DEFAULT_TPM_THRESHOLDS: tuple[int, ...] = (25, 50, 100, 200)
DEFAULT_PERCENTILE_THRESHOLDS: tuple[int, ...] = (90, 95)

# Mirrors pirlygenes' patient-coverage tabulation cutoffs so the generated CSVs
# are comparable across packages.
DEFAULT_THRESHOLDS: tuple[int, ...] = DEFAULT_TPM_THRESHOLDS

ACTIONABLE_TPM: float = 30.0
COVERING_SET_STATS: dict[str, tuple[str, str]] = {
    "median": ("expression", "median (≥50% of patients)"),
    "q3": ("q3", "upper quartile (≥25% of patients)"),
}


@dataclass(frozen=True)
class Threshold:
    """One coverage cutoff (``kind`` is ``"tpm"`` or ``"percentile"``).

    ``"pctile"`` is accepted as an alias for ``"percentile"``. Absolute TPM uses
    strict ``>`` comparison; percentile mode uses ``>=`` because a p95 call means
    a gene ranks at or above the 95th percentile within that sample.
    """

    kind: str
    value: int

    def __post_init__(self) -> None:
        kind = str(self.kind).strip().lower()
        if kind == "pctile":
            kind = "percentile"
        if kind not in {"tpm", "percentile"}:
            raise ValueError("threshold kind must be 'tpm' or 'percentile'")
        numeric = float(self.value)
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise ValueError("coverage threshold values must be finite integers")
        value = int(numeric)
        if kind == "tpm" and value < 0:
            raise ValueError("TPM thresholds must be non-negative")
        if kind == "percentile" and not 0 < value <= 100:
            raise ValueError("percentile thresholds must be in (0, 100]")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "value", value)

    @property
    def slug(self) -> str:
        return f"t{self.value}" if self.kind == "tpm" else f"p{self.value}"

    @property
    def count_suffix(self) -> str:
        return f"gt{self.value}" if self.kind == "tpm" else f"p{self.value}"

    @property
    def xlabel(self) -> str:
        return (
            f"> {self.value} clean TPM"
            if self.kind == "tpm"
            else f"≥ {self.value}th within-sample percentile"
        )

    @property
    def label(self) -> str:
        return self.xlabel

    @property
    def count_col(self) -> str:
        return f"n_{self.count_suffix}"

    def cutoff(self, cols, pctile_cutoffs=None):
        """Scalar TPM or a per-sample percentile cutoff vector."""
        if self.kind == "tpm":
            return self.value
        col = f"p{self.value}"
        series = pctile_cutoffs[col] if pctile_cutoffs is not None else None
        return np.array(
            [series.get(c, np.inf) if series is not None else np.inf for c in cols],
            dtype=float,
        )

    def compare(self, values, cutoff=None):
        """Return the boolean hit mask for ``values``."""
        threshold = self.value if cutoff is None else cutoff
        arr = np.asarray(values, dtype=float)
        if self.kind == "tpm":
            return (arr > threshold) & np.isfinite(arr) & np.isfinite(threshold)
        return (arr >= threshold) & np.isfinite(arr) & np.isfinite(threshold)

    def hits(self, values, cols=(), pctile_cutoffs=None):
        """Compare raw values using this threshold's scalar/vector cutoff."""
        return self.compare(values, self.cutoff(cols, pctile_cutoffs))


# The two within-sample percentile cuts used by the pirlygenes CTA plot catalog.
# Values are percentile ranks on [0, 1], so 0.95 means the top 5% of biological
# expression values within each individual tumor.
DEFAULT_WITHIN_SAMPLE_PERCENTILES: tuple[float, ...] = (0.90, 0.95)

_BASE = ["Ensembl_Gene_ID", "Symbol"]

_GREEDY_COVERAGE_COLUMNS = [
    "rank",
    "Ensembl_Gene_ID",
    "Symbol",
    "marginal_patients",
    "marginal_fraction",
    "cumulative_patients",
    "cumulative_fraction",
    "proteoform_key",
    "proteoform_members",
]

_PERCENTILE_COVERAGE_COLUMNS = [
    "cancer_code",
    "within_sample_percentile",
    "n_patients",
    *_GREEDY_COVERAGE_COLUMNS,
]


def _clean_id(value: str) -> str:
    return str(value).split(".", 1)[0].strip()


def _symbol_to_id_map() -> dict[str, str]:
    """Case-insensitive symbol/synonym -> canonical ENSG, using bundled data only."""
    df = canonical_gene_space()
    out = {
        str(sym).strip().upper(): _clean_id(gid)
        for gid, sym in zip(df["ensembl_gene_id"], df["symbol"])
        if not pd.isna(gid) and not pd.isna(sym) and str(gid).strip() and str(sym).strip()
    }
    # Keep canonical symbols authoritative, then add NCBI synonym aliases when the
    # official symbol exists in the canonical space.
    for alias, official in symbol_synonyms().items():
        if official and official.upper() in out:
            out.setdefault(alias.upper(), out[official.upper()])
    return out


def _symbols_to_ids(symbols: Iterable[str]) -> set[str]:
    by_symbol = _symbol_to_id_map()
    out: set[str] = set()
    for symbol in symbols:
        s = str(symbol).strip()
        if not s or s.lower() == "nan":
            continue
        if is_ensembl_gene_id(s):
            out.add(_clean_id(s).upper())
            continue
        gid = by_symbol.get(resolve_symbol(s).upper()) or by_symbol.get(s.upper())
        if gid:
            out.add(gid)
    return out


def _ids_from_frame(df: pd.DataFrame, *, symbol_col: str | None = None) -> set[str]:
    ids: set[str] = set()
    id_col = gene_id_column(df)
    has_id = pd.Series(False, index=df.index)
    if id_col is not None:
        raw_ids = df[id_col].fillna("").astype(str).str.strip()
        has_id = raw_ids.ne("") & raw_ids.str.lower().ne("nan")
        ids |= {_clean_id(v) for v in raw_ids[has_id]}
    resolved_sym_col = (
        symbol_col if symbol_col is not None and symbol_col in df.columns else symbol_column(df)
    )
    if resolved_sym_col is not None:
        ids |= _symbols_to_ids(df.loc[~has_id, resolved_sym_col].dropna().astype(str))
    return ids


def _gene_set_from_file(path: Path) -> tuple[str, set[str]]:
    # Read row 0 as data unless it matches a recognized header column name; bare
    # one-gene-per-line TXT/CSV/TSV panels otherwise silently drop their first gene.
    suffix = path.suffix.lower()
    sep = "\t" if suffix in {".tsv", ".tab"} else ","
    try:
        raw = pd.read_csv(path, header=None, dtype=str, keep_default_na=False, sep=sep)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        try:
            raw = pd.read_csv(path, header=None, dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError as exc:
            raise ValueError(
                f"gene-set file {path} did not contain resolvable Ensembl IDs/symbols"
            ) from exc
    if raw.empty:
        raise ValueError(f"gene-set file {path} did not contain resolvable Ensembl IDs/symbols")
    headers = [str(v).strip() for v in raw.iloc[0].tolist()]
    lower_headers = [h.lower() for h in headers]
    recognized = {"ensembl_gene_id", "symbol", "gene", "gene_id", "gene_symbol"}
    if any(value in recognized for value in lower_headers):
        df = raw.iloc[1:].copy()
        df.columns = headers
        col_map = {}
        for col in df.columns:
            cl = col.lower()
            if cl in {"ensembl_gene_id", "gene_id"} and "Ensembl_Gene_ID" not in col_map.values():
                col_map[col] = "Ensembl_Gene_ID"
            elif cl in {"symbol", "gene", "gene_symbol"} and "Symbol" not in col_map.values():
                col_map[col] = "Symbol"
        df = df.rename(columns=col_map)
        ids = _ids_from_frame(df)
    else:
        tokens = [str(v).strip() for v in raw.iloc[:, 0].tolist() if str(v).strip()]
        ids = _symbols_to_ids(tokens)
    if not ids:
        raise ValueError(f"gene-set file {path} did not contain resolvable Ensembl IDs/symbols")
    return path.stem, ids


_ANTIGEN_SET_LOOKUP = {s.lower(): s for s in antigen_sets.CANDIDATE_SETS}


def resolve_gene_set(name: str) -> tuple[str, set[str]]:
    """Resolve a plot gene-set token to ``(label, unversioned Ensembl IDs)``.

    Supported tokens mirror pirlygenes where oncoref owns the data:
    any current :data:`oncoref.antigen_sets.CANDIDATE_SETS` token (``cta``, ``ca``,
    ``ra``, ``ca|ra``, ``mra``, ``fra``, ``era``, ``sra``, ``cta-m``, ``cta-f``,
    ``cta-e``, ``cta-mf``), ``legacy-cta`` (historical 1.9.15 HPA snapshot),
    ``mito``/``mitochondrial``, ``housekeeping``, ``family:<name>``,
    ``therapy:<agent_class>``, ``lineage:<cancer_code>``, or a CSV/TSV/TXT file
    with ``Ensembl_Gene_ID`` and/or ``Symbol`` (or bare one-gene-per-line tokens).
    """
    token = str(name).strip()
    low = token.lower()
    if low == "cta":
        return "CTA", {_clean_id(g) for g in cta_gene_ids()}
    if low in _ANTIGEN_SET_LOOKUP:
        canonical = _ANTIGEN_SET_LOOKUP[low]
        return canonical, {_clean_id(g) for g in antigen_sets.gene_ids(canonical)}
    if low == "legacy-cta":
        from . import cta as legacy_cta

        return "legacy-cta", {_clean_id(g) for g in legacy_cta.cta_gene_ids()}
    if low in {"mito", "mitochondrial"}:
        return "mitochondrial", set(gene_family_ids("mitochondrial"))
    if low == "housekeeping":
        return "housekeeping", set(housekeeping_gene_ids())
    if low.startswith("family:"):
        family = token.split(":", 1)[1].strip()
        if family not in gene_families():
            raise ValueError(f"unknown gene family {family!r}; one of {sorted(gene_families())}")
        return f"family:{family}", set(gene_family_ids(family))
    if low.startswith("therapy:"):
        agent_class = token.split(":", 1)[1].strip().lower()
        df = cancer_key_genes_df()
        target = df["role"].astype(str).str.lower() == "target"
        cls = df["agent_class"].astype(str).str.lower() == agent_class
        ids = _symbols_to_ids(df.loc[target & cls, "symbol"].dropna().astype(str))
        if not ids:
            raise ValueError(f"no therapy target genes found for agent_class={agent_class!r}")
        return f"therapy:{agent_class}", ids
    if low.startswith("lineage:"):
        code = resolve_cancer_type(token.split(":", 1)[1])
        df = cancer_type_genes_df()
        ids = _ids_from_frame(df[df["Cancer_Type"].astype(str) == code])
        if not ids:
            raise ValueError(f"no lineage gene panel found for {code}")
        return f"lineage:{code}", ids
    path = Path(token).expanduser()
    if path.exists():
        return _gene_set_from_file(path)
    raise ValueError(
        f"unknown gene set {name!r}; use cta, mito, housekeeping, family:<name>, "
        "therapy:<agent_class>, lineage:<code>, or a path to a CSV/TSV/TXT gene panel"
    )


def _panel_ids(gene_ids: Iterable[str] | None) -> set[str]:
    ids = set(gene_ids) if gene_ids is not None else set(cta_gene_ids())
    return {str(g).split(".")[0] for g in ids}


def _hit_matrix(cancer_type, *, threshold_tpm: float, gene_ids, proteoform: bool = True):
    """``(panel rows, sample columns, gene×sample boolean hit matrix, id columns)``
    where a hit is clean TPM > ``threshold_tpm`` for one antigen in one patient.

    With ``proteoform=True`` (default), identical-protein paralogs in the panel
    (CTAG1A+CTAG1B = NY-ESO-1, XAGE1A+XAGE1B, …) are **summed per patient** into one
    antigen row before thresholding — the biologically correct unit for antigen
    coverage: RNA-seq reads multi-map between the loci (so per-gene TPM under-counts
    the proteoform), and a TCR/vaccine targets the shared protein once. With
    ``proteoform=False`` each gene is kept separate."""
    df = per_sample_expression(cancer_type, normalize="tpm_clean")
    unversioned = df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    panel = _panel_ids(gene_ids)
    sub = df[unversioned.isin(panel)].reset_index(drop=True)
    if proteoform and len(sub):
        from .proteoforms import collapse_to_proteoforms

        sub = collapse_to_proteoforms(sub).reset_index(drop=True)
    id_cols = [c for c in _ANTIGEN_ID_COLS if c in sub.columns]
    samples = [c for c in sub.columns if c not in id_cols]
    hits = sub[samples].to_numpy(dtype=float) > float(threshold_tpm)
    return sub, samples, hits, id_cols


def cta_patient_fractions(
    cancer_type,
    *,
    threshold_tpm: float = DEFAULT_EXPRESSED_TPM,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
) -> pd.DataFrame:
    """Per antigen, the fraction of a cohort's patients expressing it above
    ``threshold_tpm`` (clean TPM). Returns ``Ensembl_Gene_ID`` (a real ENSG — the
    canonical member for a collapsed proteoform), ``Symbol`` (the proteoform symbol —
    alias ``NY-ESO-1`` / contracted ``XAGE1A/B`` — when collapsed), ``proteoform_members``
    (the sorted member symbols, present when ``proteoform=True``), ``fraction_expressing``,
    ``n_patients_expressing``, ``n_patients`` — sorted by prevalence. The default
    gene panel is the expressed CTA set; identical-protein paralogs are summed to
    one antigen (``proteoform=True``)."""
    sub, samples, hits, id_cols = _hit_matrix(
        cancer_type, threshold_tpm=threshold_tpm, gene_ids=gene_ids, proteoform=proteoform
    )
    n = len(samples)
    out = sub[id_cols].copy()
    counts = hits.sum(axis=1)
    out["n_patients_expressing"] = counts
    out["fraction_expressing"] = counts / n if n else 0.0
    out["n_patients"] = n
    return out.sort_values("fraction_expressing", ascending=False).reset_index(drop=True)


def addressable_fraction(
    cancer_type,
    *,
    threshold_tpm: float = DEFAULT_EXPRESSED_TPM,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
) -> float:
    """Fraction of a cohort's patients expressing **at least one** antigen in the
    panel above ``threshold_tpm`` — the faithful "addressable" share (the union across
    patients, which the per-gene fractions can't be summed into). Identical-protein
    paralogs are summed to one antigen (``proteoform=True``). 0.0 for an empty
    cohort/panel."""
    _, samples, hits, _ = _hit_matrix(
        cancer_type, threshold_tpm=threshold_tpm, gene_ids=gene_ids, proteoform=proteoform
    )
    if not samples or hits.size == 0:
        return 0.0
    return float(hits.any(axis=0).mean())


def _empty_greedy_coverage() -> pd.DataFrame:
    """Schema-stable empty result for a cohort with no coverable patients."""
    return pd.DataFrame(columns=_GREEDY_COVERAGE_COLUMNS)


def _greedy_coverage_from_hits(
    panel: pd.DataFrame,
    samples: list[str],
    hits: np.ndarray,
    *,
    max_genes: int | None,
) -> pd.DataFrame:
    """Run deterministic greedy set cover over an already-defined hit matrix."""
    n_patients = len(samples)
    if n_patients == 0 or hits.size == 0:
        return _empty_greedy_coverage()

    rows: list[dict] = []
    covered = np.zeros(n_patients, dtype=bool)
    remaining = set(range(len(panel)))
    total_prevalence = hits.sum(axis=1)
    limit = len(panel) if max_genes is None else max_genes

    while remaining and len(rows) < limit:
        best_gene = None
        best_new = 0
        best_prevalence = -1
        for gene_index in sorted(remaining):
            newly_covered = int(np.count_nonzero(hits[gene_index] & ~covered))
            score = (newly_covered, int(total_prevalence[gene_index]))
            if score > (best_new, best_prevalence):
                best_gene = gene_index
                best_new, best_prevalence = score
        if best_gene is None or best_new == 0:
            break

        covered |= hits[best_gene]
        remaining.remove(best_gene)
        cumulative = int(covered.sum())
        row = {
            "rank": len(rows) + 1,
            "Ensembl_Gene_ID": str(panel.at[best_gene, "Ensembl_Gene_ID"]),
            "Symbol": str(panel.at[best_gene, "Symbol"]),
            "marginal_patients": best_new,
            "marginal_fraction": best_new / n_patients,
            "cumulative_patients": cumulative,
            "cumulative_fraction": cumulative / n_patients,
            "proteoform_key": pd.NA,
            "proteoform_members": pd.NA,
        }
        for column in ("proteoform_key", "proteoform_members"):
            if column in panel.columns:
                row[column] = str(panel.at[best_gene, column])
        rows.append(row)
    columns = list(_GREEDY_COVERAGE_COLUMNS)
    if "proteoform_key" not in panel.columns:
        columns = columns[:-2]
    return pd.DataFrame(rows, columns=columns)


def greedy_coverage(
    cancer_type,
    *,
    threshold_tpm: float = DEFAULT_EXPRESSED_TPM,
    gene_ids: Iterable[str] | None = None,
    max_genes: int | None = None,
    proteoform: bool = True,
) -> pd.DataFrame:
    """Greedy set-cover panel: at each step add the antigen covering the most
    patients not yet covered by the chosen panel, until every coverable patient is
    covered (or ``max_genes`` is reached). Identical-protein paralogs are summed to
    one antigen (``proteoform=True``), so e.g. CTAG1A/CTAG1B counts once.

    Returns one row per chosen antigen, in selection order: ``rank``,
    ``Ensembl_Gene_ID``, ``Symbol``, ``marginal_patients`` (newly covered),
    ``marginal_fraction``, ``cumulative_patients``, ``cumulative_fraction``. The
    cumulative fraction is the coverage curve; its last value equals
    :func:`addressable_fraction` once the panel is exhausted (unless ``max_genes``
    truncates it first). Ties are broken by total prevalence (deterministic)."""
    panel, samples, hits, _ = _hit_matrix(
        cancer_type, threshold_tpm=threshold_tpm, gene_ids=gene_ids, proteoform=proteoform
    )
    return _greedy_coverage_from_hits(panel, samples, hits, max_genes=max_genes)


def _validate_within_sample_percentiles(percentiles) -> tuple[float, ...]:
    """Normalize percentile ranks and enforce the shipped within-sample contract."""
    from .expression_builders import WITHIN_SAMPLE_THRESHOLDS

    raw = (percentiles,) if isinstance(percentiles, (int, float)) else tuple(percentiles)
    values = tuple(dict.fromkeys(float(value) for value in raw))
    allowed = tuple(sorted(WITHIN_SAMPLE_THRESHOLDS))
    invalid = [value for value in values if value not in WITHIN_SAMPLE_THRESHOLDS]
    if not values or invalid:
        raise ValueError(f"percentiles must contain one or more of {list(allowed)}")
    return values


def _normalize_cohort_codes(cohorts: Iterable[str] | str) -> list[str]:
    """Resolve one code or an iterable while retaining first-seen order."""
    requested = [cohorts] if isinstance(cohorts, str) else list(cohorts)
    return list(
        dict.fromkeys(resolve_cancer_type(code, strict=False) or str(code) for code in requested)
    )


def _within_sample_percentile_panel_ranks(
    cancer_type,
    *,
    gene_ids: Iterable[str] | None,
    proteoform: bool,
    sample_qc: str,
) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    """Return panel rows and their whole-transcriptome ranks in each patient."""
    from .expression import _biological_per_sample

    biological = _biological_per_sample(
        cancer_type,
        proteoform=proteoform,
        auto_fetch=False,
        scope="cta",
        sample_qc=sample_qc,
    )
    samples = [column for column in biological.columns if column not in _ANTIGEN_ID_COLS]
    if not samples:
        raise ValueError(f"no per-sample columns to rank for {cancer_type!r}")

    panel_ids = _panel_ids(gene_ids)
    if proteoform:
        from .proteoforms import gene_to_proteoform_id

        panel_keys = set(gene_to_proteoform_id(sorted(panel_ids), scope="cta").values())
        in_panel = biological["proteoform_key"].astype(str).isin(panel_keys)
    else:
        row_ids = biological["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        in_panel = row_ids.isin(panel_ids)

    panel = biological.loc[in_panel].reset_index(drop=True)
    if panel.empty:
        return panel, samples, np.empty((0, len(samples)), dtype=float)

    # Use the same rank definition as the shipped within-sample artifacts:
    # average ranks for ties, divided by the number of measured biological rows.
    ranks = biological[samples].rank(axis=0, pct=True)
    panel_ranks = ranks.loc[in_panel].to_numpy(dtype=float)
    return panel, samples, panel_ranks


def _zero_coverage_row(code: str, percentile: float, n_patients: int) -> dict:
    """Represent an audited cohort with no CTA-positive patient at one cutoff."""
    return {
        "cancer_code": code,
        "within_sample_percentile": percentile,
        "n_patients": n_patients,
        "rank": 0,
        "Ensembl_Gene_ID": pd.NA,
        "Symbol": pd.NA,
        "marginal_patients": 0,
        "marginal_fraction": 0.0,
        "cumulative_patients": 0,
        "cumulative_fraction": 0.0,
        "proteoform_key": pd.NA,
        "proteoform_members": pd.NA,
    }


def within_sample_percentile_coverage_sweep(
    cohorts: Iterable[str] | None = None,
    *,
    percentiles: Iterable[float] = DEFAULT_WITHIN_SAMPLE_PERCENTILES,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
    sample_qc: str = "pass",
    on_missing: str = "raise",
) -> pd.DataFrame:
    """Co-occurrence-aware antigen coverage at within-sample percentile cuts.

    Each tumor is ranked across its biological clean-TPM transcriptome.  A panel
    antigen is positive at ``0.95`` when it is in that tumor's top 5% of ranked
    expression values.  The returned long table contains the deterministic greedy
    set-cover order for every cohort and requested percentile.  Identical-protein
    CTA loci are summed before ranking by default.

    ``cohorts=None`` uses only locally cached per-sample matrices and never fetches.
    Explicit cohorts are also read with downloads disabled.  ``on_missing='record'``
    keeps available cohorts and records missing paths in ``attrs['missing_cohorts']``;
    the default raises.  Audited cohorts with zero coverage receive one ``rank=0``
    row so absence is distinguishable from missing data.
    """
    if on_missing not in {"raise", "record"}:
        raise ValueError("on_missing must be 'raise' or 'record'")
    percentile_values = _validate_within_sample_percentiles(percentiles)

    from . import source_matrices

    if cohorts is None:
        requested = [
            code for code in source_matrices.available_cohorts() if source_matrices.is_cached(code)
        ]
    else:
        requested = cohorts
    codes = _normalize_cohort_codes(requested)

    rows: list[dict] = []
    audited: list[str] = []
    missing: dict[str, str] = {}
    for code in codes:
        try:
            panel, samples, panel_ranks = _within_sample_percentile_panel_ranks(
                code,
                gene_ids=gene_ids,
                proteoform=proteoform,
                sample_qc=sample_qc,
            )
        except FileNotFoundError as error:
            if on_missing == "raise":
                raise
            missing[code] = str(error)
            continue

        audited.append(code)
        for percentile in percentile_values:
            hits = panel_ranks >= percentile
            greedy = _greedy_coverage_from_hits(panel, samples, hits, max_genes=None)
            if greedy.empty:
                rows.append(_zero_coverage_row(code, percentile, len(samples)))
                continue
            for record in greedy.to_dict(orient="records"):
                rows.append(
                    {
                        "cancer_code": code,
                        "within_sample_percentile": percentile,
                        "n_patients": len(samples),
                        **record,
                    }
                )

    out = pd.DataFrame(rows, columns=_PERCENTILE_COVERAGE_COLUMNS)
    out.attrs.update(
        {
            "threshold_basis": "biological_clean_tpm_rank_within_sample",
            "requested_cohorts": codes,
            "audited_cohorts": audited,
            "missing_cohorts": missing,
            "percentiles": percentile_values,
            "proteoform": proteoform,
            "sample_qc": sample_qc,
        }
    )
    return out


def _coverage_at_percentile(coverage: pd.DataFrame, percentile: float) -> pd.DataFrame:
    """Validate and select one percentile from a coverage sweep."""
    required = set(_PERCENTILE_COVERAGE_COLUMNS)
    missing = sorted(required - set(coverage.columns))
    if missing:
        raise ValueError(f"coverage table is missing required columns: {missing}")
    target = _validate_within_sample_percentiles((percentile,))[0]
    values = pd.to_numeric(coverage["within_sample_percentile"], errors="coerce")
    if values.isna().any():
        raise ValueError("coverage table contains a nonnumeric within_sample_percentile")
    matches = np.isclose(values.to_numpy(dtype=float), target)
    if not coverage.empty and not matches.any():
        raise ValueError(f"coverage table does not contain percentile {target:g}")
    return coverage.loc[matches].copy()


def within_sample_percentile_addressable_fraction_by_cohort(
    cohorts: Iterable[str] | None = None,
    *,
    percentile: float = 0.95,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
    sample_qc: str = "pass",
    coverage: pd.DataFrame | None = None,
    on_missing: str = "raise",
) -> pd.Series:
    """Fraction of patients with at least one panel antigen above a rank cutoff.

    This is the faithful patient union, derived from the joint per-sample matrix;
    it is not the maximum of independent per-gene prevalence summaries.  Supply a
    precomputed :func:`within_sample_percentile_coverage_sweep` table to reuse one
    matrix pass across plots.
    """
    if coverage is None:
        coverage = within_sample_percentile_coverage_sweep(
            cohorts,
            percentiles=(percentile,),
            gene_ids=gene_ids,
            proteoform=proteoform,
            sample_qc=sample_qc,
            on_missing=on_missing,
        )
    selected = _coverage_at_percentile(coverage, percentile)
    if cohorts is not None:
        requested = set(_normalize_cohort_codes(cohorts))
        selected = selected[selected["cancer_code"].astype(str).isin(requested)]
    if selected.empty:
        return pd.Series(dtype=float, name="addressable_fraction")
    ordered = selected.sort_values(["cancer_code", "rank"], kind="stable")
    final = ordered.groupby("cancer_code", sort=False)["cumulative_fraction"].last()
    return pd.to_numeric(final, errors="raise").rename("addressable_fraction")


def mean_antigens_per_patient(
    cancer_type,
    *,
    threshold_tpm: float = DEFAULT_EXPRESSED_TPM,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
) -> float:
    """Mean number of panel antigens a patient in the cohort expresses above
    ``threshold_tpm`` — the per-patient antigen *load* (how many CTAs the average
    patient presents, not just whether ≥1). Equals the sum over antigens of their
    per-patient prevalence; identical-protein paralogs count once
    (``proteoform=True``). 0.0 for an empty cohort/panel."""
    _, samples, hits, _ = _hit_matrix(
        cancer_type, threshold_tpm=threshold_tpm, gene_ids=gene_ids, proteoform=proteoform
    )
    if not samples or hits.size == 0:
        return 0.0
    return float(hits.sum(axis=0).mean())


def _selected_coverage_codes(
    *,
    source_id: str = DEFAULT_SOURCE,
    cohorts: Iterable[str] | None = None,
    codes: Iterable[str] | None = None,
) -> list[str]:
    """Select cohort codes for patient coverage across ``source_id`` and explicit filters."""
    from . import source_matrices

    explicit = codes if codes is not None else cohorts
    if explicit is not None:
        return [
            resolve_cancer_type(str(c), strict=False) or str(c)
            for c in dict.fromkeys(str(x).strip() for x in explicit if str(x).strip())
        ]

    reg = source_matrices.registry()
    cached = [c for c in source_matrices.available_cohorts() if source_matrices.is_cached(c)]
    if source_id in (None, "", "all"):
        return cached
    target = str(source_id).strip().lower()
    allowed: set[str] = set()
    for _, row in reg.iterrows():
        code = str(row.get("cancer_code", ""))
        candidates = {
            str(row.get("source_id", "")).strip().lower(),
            str(row.get("source_cohort", "")).strip().lower(),
            str(row.get("source_cohort", "")).strip().lower().replace("_", "-"),
        }
        if any(c and (c == target or c.startswith(target)) for c in candidates):
            allowed.add(code)
    return [c for c in cached if c in allowed]


def _coverage_source_metadata(
    codes: Iterable[str],
    *,
    source_id: str = DEFAULT_SOURCE,
    synthetic_default_comparable: bool = False,
) -> dict[str, dict]:
    """Owner scale/provenance metadata for selected matrix codes."""
    from . import source_matrices

    code_list = list(dict.fromkeys(str(c) for c in codes))
    registry = source_matrices.registry()
    registry_rows = {str(row["cancer_code"]): row for _, row in registry.iterrows()}
    known_codes = [code for code in code_list if code in registry_rows]
    owner = pd.DataFrame()
    if known_codes:
        owner = cancer_reference_expression_availability(
            cancer_types=known_codes,
            normalize="tpm_clean",
            sample_qc="all",
            reference_source="summary_rows_all",
            all_sources=True,
        )

    out: dict[str, dict] = {}
    for code in code_list:
        registry_row = registry_rows.get(code)
        source_cohort = (
            str(registry_row["source_cohort"])
            if registry_row is not None
            else str(source_id or "custom")
        )
        if owner.empty:
            match = owner
        else:
            match = owner.loc[
                owner["cancer_code"].astype(str).eq(code)
                & owner["source_cohort"].astype(str).eq(source_cohort)
            ]
        if match.empty and not owner.empty:
            match = owner.loc[owner["cancer_code"].astype(str).eq(code)]
        row = match.iloc[0] if not match.empty else None
        comparable = None
        if row is not None and not pd.isna(row.get("linear_tpm_comparable")):
            comparable = bool(row["linear_tpm_comparable"])
        elif registry_row is None and synthetic_default_comparable:
            comparable = True
        out[code] = {
            "source_cohort": source_cohort,
            "source_type": str(row.get("source_type", "unknown")) if row is not None else "unknown",
            "source_scale_class": (
                str(row.get("source_scale_class", "unknown")) if row is not None else "unknown"
            ),
            "linear_tpm_comparable": comparable,
            "normalization": "tpm_clean",
        }
    return out


def _resolve_threshold_mode(requested: str, metadata: dict[str, dict]) -> str:
    mode = str(requested).strip().lower()
    if mode == "pctile":
        mode = "percentile"
    if mode not in {"auto", "tpm", "percentile"}:
        raise ValueError("threshold_mode must be one of 'auto', 'tpm', or 'percentile'")
    if mode == "auto":
        if metadata and all(item["linear_tpm_comparable"] is True for item in metadata.values()):
            return "tpm"
        return "percentile"
    if mode == "tpm":
        incompatible = sorted(
            code for code, item in metadata.items() if item["linear_tpm_comparable"] is False
        )
        if incompatible:
            joined = ", ".join(incompatible)
            raise ValueError(
                "absolute TPM coverage is invalid for non-comparable source "
                f"scale(s): {joined}; use threshold_mode='percentile'"
            )
    return mode


def _threshold_values(mode: str, thresholds: Iterable[float] | None) -> tuple[int, ...]:
    if thresholds is None:
        thresholds = DEFAULT_TPM_THRESHOLDS if mode == "tpm" else DEFAULT_PERCENTILE_THRESHOLDS
    values: list[int] = []
    for value in thresholds:
        normalized = Threshold(mode, value).value
        if normalized not in values:
            values.append(normalized)
    return tuple(values)


def _panel_coverage_matrix(
    cancer_type: str,
    *,
    gene_ids: Iterable[str],
    proteoform: bool,
    percentile_mode: bool,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Load per-sample clean TPM (or full-transcriptome within-sample percentile ranks),
    optionally collapse identical-protein paralogs, and restrict to ``gene_ids``."""
    df = per_sample_expression(cancer_type, normalize="tpm_clean")
    id_cols_all = [c for c in _ANTIGEN_ID_COLS if c in df.columns]
    samples = [c for c in df.columns if c not in id_cols_all]
    panel = _panel_ids(gene_ids)

    if percentile_mode:
        # Rank the full biological transcriptome BEFORE restricting to the panel so
        # percentiles reflect whole-sample rank rather than rank within the panel.
        work = df.copy()
        if proteoform and len(work):
            from .proteoforms import collapse_to_proteoforms, gene_to_proteoform_id

            work = collapse_to_proteoforms(work).reset_index(drop=True)
            work = percentile_rank(work, value_cols=samples)
            panel_keys = set(gene_to_proteoform_id(sorted(panel), scope="cta").values()) | panel
            in_panel = work["proteoform_key"].astype(str).isin(panel_keys) | (
                work["Ensembl_Gene_ID"].astype(str).str.split(".").str[0].isin(panel)
            )
            sub = work.loc[in_panel].reset_index(drop=True)
        else:
            work = percentile_rank(work, value_cols=samples)
            unversioned = work["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
            sub = work.loc[unversioned.isin(panel)].copy()
            sub["Ensembl_Gene_ID"] = unversioned.loc[unversioned.isin(panel)]
            sub = sub.reset_index(drop=True)
    else:
        unversioned = df["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
        sub = df.loc[unversioned.isin(panel)].copy()
        sub["Ensembl_Gene_ID"] = unversioned.loc[unversioned.isin(panel)]
        sub = sub.reset_index(drop=True)
        if proteoform and len(sub):
            from .proteoforms import collapse_to_proteoforms

            sub = collapse_to_proteoforms(sub).reset_index(drop=True)

    id_cols = [c for c in _ANTIGEN_ID_COLS if c in sub.columns]
    samples = [c for c in sub.columns if c not in id_cols]
    return sub, samples, id_cols


def _coverage_frame(
    gene_ids: Iterable[str],
    selected_codes: list[str],
    metadata: dict[str, dict],
    mode: str,
    threshold_values: tuple[int, ...],
    *,
    proteoform: bool = True,
    greedy_threshold: Threshold | None = None,
) -> tuple[pd.DataFrame, list[tuple[str, int, list[float], list[str]]]]:
    """Compute per-(cohort × gene × threshold) counts and optional greedy curves in one pass."""
    thresholds = [Threshold(mode, value) for value in threshold_values]
    rows: list[dict] = []
    per: list[tuple[str, int, list[float], list[str]]] = []

    for code in selected_codes:
        sub, samples, id_cols = _panel_coverage_matrix(
            code,
            gene_ids=gene_ids,
            proteoform=proteoform,
            percentile_mode=(mode == "percentile"),
        )
        n = len(samples)
        if n == 0 or sub.empty:
            continue
        source = metadata[code]
        values = sub[samples].to_numpy(dtype=float)
        for i, r in sub.iterrows():
            vals = values[i]
            n_available = int(np.isfinite(vals).sum())
            rec = {
                "cancer_code": str(code),
                "label": format_cancer_code_label(str(code)),
                "source_cohort": source["source_cohort"],
                "source_type": source["source_type"],
                "source_scale_class": source["source_scale_class"],
                "linear_tpm_comparable": source["linear_tpm_comparable"],
                "normalization": source["normalization"],
                "threshold_mode": mode,
                "n_samples": n,
                "n_available": n_available,
                "n_patients": n,
            }
            for col in id_cols:
                rec[col] = r[col]
            any_hit = False
            for threshold in thresholds:
                count = int(threshold.compare(vals).sum())
                rec[f"n_{threshold.count_suffix}"] = count
                rec[f"pct_{threshold.count_suffix}"] = (
                    round(100.0 * count / n_available, 2) if n_available else np.nan
                )
                any_hit = any_hit or count > 0
            if any_hit:
                rows.append(rec)

        if greedy_threshold is not None:
            hits = np.vstack([greedy_threshold.compare(values[i]) for i in range(len(sub))])
            gc = _greedy_coverage_from_hits(sub, samples, hits, max_genes=None)
            if not gc.empty:
                per.append(
                    (
                        str(code),
                        n,
                        [float(v) for v in gc["cumulative_fraction"]],
                        [
                            str(sym) if str(sym) else str(gid)
                            for sym, gid in zip(gc["Symbol"], gc["Ensembl_Gene_ID"])
                        ],
                    )
                )

    columns = (
        [
            "cancer_code",
            "label",
            "source_cohort",
            "source_type",
            "source_scale_class",
            "linear_tpm_comparable",
            "normalization",
            "threshold_mode",
            "n_samples",
            "n_available",
            "n_patients",
        ]
        + [
            c
            for c in ("Ensembl_Gene_ID", "Symbol", "proteoform_key", "proteoform_members")
            if any(c in r for r in rows) or c in ("Ensembl_Gene_ID", "Symbol")
        ]
        + [
            f"{prefix}_{threshold.count_suffix}"
            for threshold in thresholds
            for prefix in ("n", "pct")
        ]
    )
    out = pd.DataFrame(rows, columns=columns)
    out.attrs.update(
        {
            "threshold_mode": mode,
            "thresholds": tuple(threshold.value for threshold in thresholds),
            "source_metadata": metadata,
        }
    )
    return out, per


def patient_coverage(
    gene_set: str = "cta",
    cohorts: Iterable[str] | None = None,
    thresholds: Iterable[float] | None = None,
    *,
    source_id: str = DEFAULT_SOURCE,
    codes: Iterable[str] | None = None,
    threshold_mode: str | None = None,
    proteoform: bool = True,
) -> pd.DataFrame:
    """Per-cohort/per-gene patient coverage under an explicit threshold contract.

    ``threshold_mode="auto"`` (default) uses owner ``linear_tpm_comparable`` metadata:
    clean TPM for wholly comparable sources, otherwise within-sample percentile rank.
    Explicit ``"tpm"`` is rejected for sources marked non-comparable; ``"percentile"``
    is always platform-safe. TPM output uses ``n_gt25``/``pct_gt25`` columns;
    percentile output uses ``n_p90``/``pct_p90``. Percentages divide by ``n_available``
    (finite measurements).
    """
    _label, gene_ids = resolve_gene_set(gene_set)
    selected = _selected_coverage_codes(source_id=source_id, cohorts=cohorts, codes=codes)
    legacy_synthetic = threshold_mode is None and cohorts is not None and codes is None
    metadata = _coverage_source_metadata(
        selected,
        source_id=source_id,
        synthetic_default_comparable=legacy_synthetic,
    )
    mode = _resolve_threshold_mode(threshold_mode or "auto", metadata)
    threshold_values = _threshold_values(mode, thresholds)
    return _coverage_frame(
        gene_ids,
        selected,
        metadata,
        mode,
        threshold_values,
        proteoform=proteoform,
    )[0]


def _scalar_by_cohort(
    scalar_fn, name, cohorts, *, threshold_tpm, gene_ids, proteoform
) -> pd.Series:
    """Map a per-cohort scalar coverage function over cohorts → ``Series``. When
    ``cohorts`` is ``None`` every cohort with a *cached* per-sample matrix is used
    (uncached ones are skipped, never fetched implicitly); an explicit ``cohorts``
    list is taken as-is."""
    from . import source_matrices

    codes = list(cohorts) if cohorts is not None else source_matrices.available_cohorts()
    out: dict[str, float] = {}
    for code in codes:
        if cohorts is None and not source_matrices.is_cached(code):
            continue
        out[str(code)] = scalar_fn(
            code, threshold_tpm=threshold_tpm, gene_ids=gene_ids, proteoform=proteoform
        )
    return pd.Series(out, name=name)


def mean_antigens_per_patient_by_cohort(
    cohorts: Iterable[str] | None = None,
    *,
    threshold_tpm: float = DEFAULT_EXPRESSED_TPM,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
) -> pd.Series:
    """``{cohort code -> mean antigens per patient}`` over the cohorts that have a
    cached per-sample matrix (default: all cached ones); uncached cohorts are skipped
    rather than fetched."""
    return _scalar_by_cohort(
        mean_antigens_per_patient,
        "mean_antigens_per_patient",
        cohorts,
        threshold_tpm=threshold_tpm,
        gene_ids=gene_ids,
        proteoform=proteoform,
    )


def addressable_fraction_by_cohort(
    cohorts: Iterable[str] | None = None,
    *,
    threshold_tpm: float = DEFAULT_EXPRESSED_TPM,
    gene_ids: Iterable[str] | None = None,
    proteoform: bool = True,
) -> pd.Series:
    """``{cohort code -> addressable fraction}`` over the cohorts that have a cached
    per-sample matrix (default: all cached ones); uncached cohorts are skipped rather
    than fetched."""
    return _scalar_by_cohort(
        addressable_fraction,
        "addressable_fraction",
        cohorts,
        threshold_tpm=threshold_tpm,
        gene_ids=gene_ids,
        proteoform=proteoform,
    )


def cta_id_to_name() -> dict[str, str]:
    """``{unversioned CTA gene id -> symbol}`` for labelling coverage outputs."""
    return cta_gene_id_to_name()


# --- pirlygenes-style rendering -------------------------------------------------------

_PALETTE = [
    "#e6194B",
    "#3cb44b",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#469990",
    "#9A6324",
    "#800000",
    "#808000",
    "#000075",
    "#e6beff",
    "#aaffc3",
    "#ffd8b1",
    "#a9a9a9",
    "#fabed4",
]


def _slug(label: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(label).lower()).strip("_")


def _gene_color_map(genes_ordered):
    seen, colors = [], {}
    for gene in genes_ordered:
        if gene not in colors:
            colors[gene] = _PALETTE[len(seen) % len(_PALETTE)]
            seen.append(gene)
    return colors


def _as_threshold(threshold: Threshold | float | int) -> Threshold:
    return threshold if isinstance(threshold, Threshold) else Threshold("tpm", int(threshold))


def _stacked_bar(per, label, threshold, path, plt):
    """Horizontal stacked bar: each cohort's greedy plateau split into each
    gene's marginal new-patient contribution (segments sum to the plateau)."""
    from collections import Counter

    threshold_obj = _as_threshold(threshold)
    totals = Counter()
    for _code, _n, cum, names in per:
        prev = 0.0
        for name, c in zip(names, cum):
            totals[name] += (c - prev) * 100
            prev = c
    color = _gene_color_map([g for g, _ in totals.most_common()])

    fig, ax = plt.subplots(figsize=(13, max(6, len(per) * 0.28)))
    labels = []
    for y, (code, n, cum, names) in enumerate(per):
        labels.append(f"{format_cancer_code_label(code)}  (n={n})")
        left, prev = 0.0, 0.0
        for j, (name, c) in enumerate(zip(names, cum)):
            marginal = (c - prev) * 100
            prev = c
            if marginal <= 0:
                continue
            ax.barh(
                y,
                marginal,
                left=left,
                color=color.get(name, "#cccccc"),
                edgecolor="white",
                linewidth=0.3,
            )
            if marginal >= 1.5 and (marginal >= 3.0 or j == 0):
                ax.text(
                    left + marginal / 2,
                    y,
                    name,
                    va="center",
                    ha="center",
                    fontsize=4.5,
                    clip_on=True,
                )
            left += marginal
    ax.set_yticks(range(len(per)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(0, 100)
    ax.set_xlabel(
        f"% of patients with ≥1 {label} gene {threshold_obj.xlabel} "
        "(stacked by each gene's marginal new-patient share, greedy)"
    )
    ax.grid(axis="x", alpha=0.3)
    ax.set_title(
        f"{label} coverage by cancer type, split by gene "
        f"({threshold_obj.xlabel}, {len(per)} cohorts)",
        fontsize=11,
    )
    fig.tight_layout()
    save_figure(fig, path, dpi=300)
    plt.close(fig)
    return path


def _coverage_curves(per, label, threshold, path, plt):
    """Small-multiples of each cohort's greedy coverage curve (sorted by
    plateau, broadest first)."""
    threshold_obj = _as_threshold(threshold)
    ordered = sorted(per, key=lambda t: t[2][-1], reverse=True)
    ncol = 6
    nrow = (len(ordered) + ncol - 1) // ncol
    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(ncol * 2.5, nrow * 1.9),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    axes = axes.ravel()
    for ax, (code, n, cum, names) in zip(axes, ordered):
        xs = range(1, len(cum) + 1)
        ax.plot(xs, [c * 100 for c in cum], color="#b5179e", lw=1.2)
        ax.fill_between(xs, [c * 100 for c in cum], alpha=0.15, color="#b5179e")
        for x, (name, c) in enumerate(zip(names[:3], cum[:3]), start=1):
            ax.annotate(
                name,
                (x, c * 100),
                fontsize=4,
                rotation=45,
                textcoords="offset points",
                xytext=(1, 2),
            )
        ax.set_title(
            f"{format_cancer_code_label(code)} (n={n}) {cum[-1] * 100:.0f}%",
            fontsize=7,
        )
        ax.set_xlim(0, 25)
        ax.set_ylim(0, 100)
        ax.tick_params(labelsize=5)
        ax.grid(alpha=0.25)
    for ax in axes[len(ordered) :]:
        ax.axis("off")
    fig.suptitle(
        f"{label} panel coverage by cancer type — distinct patients "
        f"with ≥1 gene {threshold_obj.xlabel} (sorted by plateau)",
        fontsize=11,
    )
    fig.supxlabel("# genes added (greedy)", fontsize=8)
    fig.supylabel("% patients covered", fontsize=8)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.97))
    save_figure(fig, path, dpi=300)
    plt.close(fig)
    return path


def render_patient_coverage(
    gene_set: str = "cta",
    source_id: str = DEFAULT_SOURCE,
    codes: Iterable[str] | None = None,
    threshold: float | int | None = None,
    thresholds: Iterable[float] | None = None,
    out_dir="coverage_out",
    *,
    cohorts: Iterable[str] | None = None,
    threshold_mode: str | None = None,
    proteoform: bool = True,
) -> dict:
    """Compute patient coverage for ``gene_set`` and write a counts CSV plus
    300-dpi PNG + vector PDF figures (stacked coverage bar and coverage-curve
    small-multiples) into ``out_dir``.
    """
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    label, gene_ids = resolve_gene_set(gene_set)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    slug = _slug(label)

    selected = _selected_coverage_codes(source_id=source_id, cohorts=cohorts, codes=codes)
    legacy_synthetic = threshold_mode is None and cohorts is not None and codes is None
    metadata = _coverage_source_metadata(
        selected,
        source_id=source_id,
        synthetic_default_comparable=legacy_synthetic,
    )
    mode = _resolve_threshold_mode(threshold_mode or "auto", metadata)
    requested_plot_value = threshold if threshold is not None else (25 if mode == "tpm" else 95)
    plot_threshold = Threshold(mode, requested_plot_value)
    plot_value = plot_threshold.value
    table_thresholds = list(_threshold_values(mode, thresholds))
    if plot_value not in table_thresholds:
        table_thresholds.append(plot_value)

    counts, per = _coverage_frame(
        gene_ids,
        selected,
        metadata,
        mode,
        tuple(table_thresholds),
        proteoform=proteoform,
        greedy_threshold=plot_threshold,
    )
    csv_path = out / f"{slug}_patient_counts.csv"
    if plot_threshold.count_col in counts.columns:
        counts = counts.sort_values(
            ["cancer_code", plot_threshold.count_col], ascending=[True, False]
        )
    counts.to_csv(csv_path, index=False)

    per.sort(key=lambda t: t[2][-1])

    paths = {"counts_csv": str(csv_path)}
    if per:
        paths["stacked_bar"] = str(
            _stacked_bar(
                per,
                label,
                plot_threshold,
                out / f"{slug}_stacked_coverage_{plot_threshold.slug}.png",
                plt,
            )
        )
        paths["coverage_curves"] = str(
            _coverage_curves(
                per,
                label,
                plot_threshold,
                out / f"{slug}_coverage_curves_{plot_threshold.slug}.png",
                plt,
            )
        )
    return {
        "paths": paths,
        "counts": counts,
        "label": label,
        "n_cohorts": len(per),
        "threshold_mode": mode,
        "threshold": plot_value,
        "threshold_label": plot_threshold.xlabel,
    }


#: pirlygenes compatibility alias.
render = render_patient_coverage


# --- covering set: one panel across all cancer types -----------------------


def burden_weights(
    codes: Iterable[str], *, metric: str = "us_incidence_pct"
) -> tuple[pd.Series, list[str]]:
    """Patient weight per cancer code, from the curated disease-burden table.

    Returns ``(weights, unmapped)``: a ``code -> weight`` Series and the list of
    codes with no burden category. Each category's share is split evenly across
    the codes from ``codes`` that map to it; unmapped codes retain weight ``0.0``
    and are returned explicitly rather than given a fabricated floor.
    """
    unique_codes = list(dict.fromkeys(str(c) for c in codes))
    shares = cancer_burden(metric=metric)
    by_category: dict[str, list[str]] = {}
    unmapped: list[str] = []
    for code in unique_codes:
        cat = burden_category(code)
        pct = shares.get(cat) if cat else None
        if cat is None or pct is None or pd.isna(pct):
            unmapped.append(code)
            continue
        by_category.setdefault(cat, []).append(code)
    weights = dict.fromkeys(unique_codes, 0.0)
    for cat, members in by_category.items():
        per_member = float(shares[cat]) / len(members)
        for code in members:
            weights[code] = per_member
    return pd.Series(weights, dtype=float), unmapped


def representative_sources(
    df: pd.DataFrame, *, min_samples: int = 10, min_display_samples: int = 5
) -> pd.DataFrame:
    """One ``source_cohort`` per ``cancer_code``: the most gene-rich source that
    still clears the display floor."""
    meta = df.groupby(["cancer_code", "source_cohort"], as_index=False).agg(
        n_samples=("n_samples", "first"),
        n_genes=("Ensembl_Gene_ID", "nunique"),
    )
    eligible = set(meta.loc[meta["n_samples"] >= min_samples, "cancer_code"])
    displayable = meta[
        meta["cancer_code"].isin(eligible) & (meta["n_samples"] >= min_display_samples)
    ]
    return (
        displayable.sort_values(["cancer_code", "n_genes", "n_samples"])
        .groupby("cancer_code", as_index=False)
        .tail(1)
    )


def _fold_to_cdna_canonical_ids(ensgs: Iterable[str]) -> list[str]:
    """Fold member ENSGs to their cDNA-identical group identity (singletons keep ENSG)."""
    from .expression import _identical_locus_identity_maps

    member_to_identity, _ = _identical_locus_identity_maps("cdna", "oncoref")
    out: list[str] = []
    for gid in ensgs:
        clean = _clean_id(gid)
        if not clean:
            continue
        ident = member_to_identity.get(clean, clean)
        if ident not in out:
            out.append(ident)
    return out


def covering_set_matrix(
    stat: str = "q3",
    *,
    gene_set: str = "cta",
    codes: Iterable[str] | None = None,
    df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """``cancer_code`` × gene matrix of clean TPM at one statistic, one representative
    source per code — the input to :func:`greedy_covering_set`."""
    if stat not in COVERING_SET_STATS:
        raise ValueError(f"stat must be one of {sorted(COVERING_SET_STATS)}, got {stat!r}")
    value_col, _ = COVERING_SET_STATS[stat]
    _label, ensgs = resolve_gene_set(gene_set)
    panel_ids = set(_fold_to_cdna_canonical_ids(ensgs))
    if df is None:
        df = cancer_reference_expression(collapse_cdna_identical=True)
    rep = representative_sources(df)
    keep = set(zip(rep["cancer_code"], rep["source_cohort"]))
    group_ids: dict[str, str] = {}
    if "Member_Ensembl_Gene_IDs" in df.columns:
        groups = df.loc[
            ~df["Ensembl_Gene_ID"].astype(str).map(is_ensembl_gene_id),
            ["Ensembl_Gene_ID", "Member_Ensembl_Gene_IDs"],
        ].drop_duplicates()
        for group_id, members in groups.itertuples(index=False, name=None):
            if pd.notna(members):
                folded = _fold_to_cdna_canonical_ids(str(members).split(";"))
                if len(folded) == 1:
                    group_ids[str(group_id)] = folded[0]
    from .expression import _identical_locus_identity_maps

    member_to_identity, _ = _identical_locus_identity_maps("cdna", "oncoref")
    frame_ids = (
        df["Ensembl_Gene_ID"]
        .astype(str)
        .map(lambda g: member_to_identity.get(_clean_id(g), str(g)))
        .replace(group_ids)
    )
    selected = frame_ids.isin(panel_ids)
    sub = df.loc[selected].copy()
    sub["Ensembl_Gene_ID"] = frame_ids.loc[selected]
    sub = sub.loc[[(c, s) in keep for c, s in zip(sub["cancer_code"], sub["source_cohort"])]]
    population = sorted(set(rep["cancer_code"]))
    if codes is not None:
        wanted = {resolve_cancer_type(str(c), strict=False) or str(c) for c in codes}
        sub = sub[sub["cancer_code"].isin(wanted)]
        population = [code for code in population if code in wanted]
    matrix = sub.pivot_table(
        index="cancer_code", columns="Ensembl_Gene_ID", values=value_col, aggfunc="max"
    ).reindex(population)
    fallback = (
        sub[["Ensembl_Gene_ID", "Symbol"]]
        .dropna()
        .sort_values(["Ensembl_Gene_ID", "Symbol"])
        .drop_duplicates("Ensembl_Gene_ID")
        .set_index("Ensembl_Gene_ID")["Symbol"]
    )
    matrix.attrs["symbols"] = {
        gid: (canonical_gene_symbol(gid) or fallback.get(gid, gid))
        if is_ensembl_gene_id(str(gid))
        else str(gid)
        for gid in matrix.columns
    }
    return matrix


@dataclass(frozen=True)
class CoverStep:
    """One pick of the greedy panel."""

    gene: str
    new_codes: tuple[str, ...]
    cum_weight: float
    cum_codes: int


def greedy_covering_set(
    matrix: pd.DataFrame, threshold: float, weights: pd.Series
) -> tuple[list[CoverStep], list[str]]:
    """Greedy weighted set cover over cancer types. Returns ``(steps, coverable)``."""
    hits = matrix > threshold
    coverable = list(hits.index[hits.any(axis=1)])
    w = weights.reindex(matrix.index).fillna(0.0)
    remaining: set[str] = set(coverable)
    steps: list[CoverStep] = []
    while remaining:
        best_gene, best_codes, best_gain = None, set(), 0.0
        for gene in hits.columns:
            covered = {str(c) for c in hits.index[hits[gene]] if str(c) in remaining}
            if not covered:
                continue
            gain = float(w.reindex(sorted(covered)).sum())
            if (gain, len(covered)) > (best_gain, len(best_codes)):
                best_gene, best_codes, best_gain = str(gene), covered, gain
        if best_gene is None:
            break
        remaining -= best_codes
        covered_so_far = sorted(set(coverable) - remaining)
        steps.append(
            CoverStep(
                gene=best_gene,
                new_codes=tuple(sorted(best_codes)),
                cum_weight=float(w.reindex(covered_so_far).sum()),
                cum_codes=len(covered_so_far),
            )
        )
    return steps, coverable


def _covering_set_plot(
    steps,
    coverable,
    weights,
    unmapped,
    *,
    total_types: int,
    label: str,
    metric: str,
    stat: str,
    threshold: float,
    n_genes: int,
    path: Path,
    plt,
) -> Path:
    """Cumulative burden and type coverage over the full input population."""
    from matplotlib.ticker import MaxNLocator, PercentFormatter

    shown = steps[:n_genes] if n_genes else steps
    total_weight = float(weights.sum())
    metric_label = {
        "us_incidence_pct": "US annual incidence",
        "world_incidence_pct": "world annual incidence",
        "us_mortality_pct": "US annual mortality",
        "world_mortality_pct": "world annual mortality",
    }.get(metric, metric.replace("_", " "))
    xs = list(range(1, len(shown) + 1))
    by_patients = [100.0 * s.cum_weight / total_weight if total_weight else 0.0 for s in shown]
    by_type = [100.0 * s.cum_codes / total_types if total_types else 0.0 for s in shown]

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(xs, by_patients, "o-", ms=5, lw=2, color="#3a0ca3", label=f"{metric_label} burden")
    ax.plot(xs, by_type, "s--", ms=5, lw=2, color="#e07a00", label="cancer types")
    ax.axhline(80, color="0.7", ls=":", lw=1)
    ax.set_xlabel(f"# genes in panel (> {threshold:g} clean TPM at {stat})")
    ax.set_ylabel("cumulative coverage")
    ax.set_ylim(0, 102)
    ax.set_xlim(left=1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.set_title(
        f"{label} covering set — how few targets reach most of the burden\n"
        f"{len(coverable)} of {total_types} cancer types coverable; "
        f"weighted by curated {metric_label}"
    )
    ax.legend(title="cumulative % covered")
    ax.grid(alpha=0.3)
    if unmapped:
        fig.text(
            0.01,
            -0.01,
            f"{len(unmapped)} cancer type(s) have no curated burden category and carry zero "
            f"burden weight (they still count toward cancer types): {', '.join(sorted(unmapped))}",
            fontsize=7,
            color="0.35",
            va="top",
            wrap=True,
        )
    fig.tight_layout()
    save_figure(fig, path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def render_covering_set(
    gene_set: str = "cta",
    *,
    stat: str = "q3",
    threshold: float = ACTIONABLE_TPM,
    codes: Iterable[str] | None = None,
    n_genes: int = 25,
    metric: str = "us_incidence_pct",
    out_dir="covering_set_out",
    df: pd.DataFrame | None = None,
) -> dict:
    """Write the covering-set figure (PNG + vector PDF) and its step CSV."""
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    label, _ensgs = resolve_gene_set(gene_set)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    slug = _slug(label)

    matrix = covering_set_matrix(stat, gene_set=gene_set, codes=codes, df=df)
    weights, unmapped = burden_weights(matrix.index, metric=metric)
    steps, coverable = greedy_covering_set(matrix, threshold, weights)

    total_weight = float(weights.sum())
    rows = pd.DataFrame(
        [
            {
                "rank": i,
                "gene": matrix.attrs.get("symbols", {}).get(s.gene, s.gene),
                "Ensembl_Gene_ID": s.gene,
                "cum_pct_patients": (
                    100.0 * s.cum_weight / total_weight if total_weight else float("nan")
                ),
                "cum_cancer_types": s.cum_codes,
                "cum_pct_cancer_types": (
                    100.0 * s.cum_codes / len(matrix.index) if len(matrix.index) else float("nan")
                ),
                "newly_covered": ";".join(s.new_codes),
            }
            for i, s in enumerate(steps, 1)
        ],
        columns=[
            "rank",
            "gene",
            "Ensembl_Gene_ID",
            "cum_pct_patients",
            "cum_cancer_types",
            "cum_pct_cancer_types",
            "newly_covered",
        ],
    )
    csv_path = out / f"{slug}_covering_set_{stat}.csv"
    rows.to_csv(csv_path, index=False)

    png = _covering_set_plot(
        steps,
        coverable,
        weights,
        unmapped,
        stat=stat,
        threshold=threshold,
        total_types=len(matrix.index),
        label=label,
        metric=metric,
        n_genes=n_genes,
        path=out / f"{slug}_covering_set_{stat}.png",
        plt=plt,
    )
    return {
        "paths": {
            "covering_set_csv": str(csv_path),
            "covering_set_png": str(png),
            "covering_set_pdf": str(png.with_suffix(".pdf")),
        },
        "steps": rows,
        "label": label,
        "stat": stat,
        "threshold": threshold,
        "metric": metric,
        "n_cancer_types": len(matrix.index),
        "n_coverable": len(coverable),
        "unmapped_codes": unmapped,
    }
