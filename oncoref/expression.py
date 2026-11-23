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

"""Read accessors over the per-cohort expression data.

This module is a small, regular surface over one underlying object — a cohort's
**per-sample expression matrix** (genes × samples) — plus the pre-computed summary
artifacts derived from it. The design is two orthogonal axes; everything here is a
point on that grid, which is why the accessors look repetitive by construction
rather than by accident.

**Axis 1 — expression level** (``proteoforms.py``). Every value is available at
*gene* level (one row per Ensembl gene) or *proteoform* level (identical-protein
paralogs summed per sample — CTAG1A+CTAG1B → NY-ESO-1 — keyed by ``proteoform_key``;
see :mod:`oncoref.proteoforms`). The base functions take ``proteoform=``/``scope=``
flags; each also has an explicit ``gene_*`` / ``proteoform_*`` wrapper so the level is
legible at the call site and discoverable by name. The wrappers are deliberately
thin — the collapse logic lives in exactly one place (the base function).

**Axis 2 — the dataset / reduction**:
  - ``per_sample_expression`` — the raw matrix itself (one column per sample).
  - ``cohort_mean_expression`` — one across-patient statistic (mean/median).
  - ``cohort_stats`` / ``pooled_cohort_stats`` — the full per-gene statistic suite,
    for one cohort or a heterogeneity-safe pool of several.
  - ``cohort_gene_percentiles`` — the within-cohort percentile vector.
  - ``within_sample_top_fraction`` — within-sample top-expression prevalence.
  - ``representative_cohort_samples`` — bounded medoid per-sample vectors.
  - ``pan_cancer_expression`` — the wide tumor+normal HPA/TCGA reference (a distinct
    data product, not the per-sample-matrix family).

**Normalize spaces** (``per_sample_expression``'s ``normalize=``): raw TPM, clean
two-compartment TPM (the default biological view), its ``log1p``, and the
housekeeping ratio (``tpm_clean_hk``). Every summary computed live inherits the space.

**Where the numbers come from.** The light artifacts (percentiles, within-sample,
representatives) are pre-computed offline by the ``scripts/generate_*`` build cores
and shipped as per-cohort parquet *shards*; the read path just resolves a code and
reads a parquet. When a shard isn't shipped (every proteoform variant, today), the
summary is **recomputed on the fly** from the per-sample matrix via the *same* build
core, so shipped and on-the-fly values agree. Which artifact ships which variant, and
whether a missing shard may be fetched, is recorded once in the :class:`ShardDataset`
registry (:data:`SHARD_DATASETS`) rather than restated per accessor.

Higher-level consumer accessors such as ``cancer_reference_expression`` are thin
views over these same shipped artifacts rather than separate data products.
"""

from __future__ import annotations

import os
import re
import warnings
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from . import data_bundle, source_matrices
from ._reference_sources import (
    TREEHOUSE_TCGA_LEGACY_COHORT,
    TREEHOUSE_TCGA_SAMPLES_COHORT,
    TREEHOUSE_TCGA_SARC_HISTOLOGY_CODES,
    TREEHOUSE_TCGA_SARC_HISTOLOGY_COHORT,
)
from .cancer_types import (
    _computed_expression_reference_members as _computed_expression_reference_members,
)
from .cancer_types import (
    cancer_type_reference_source as cancer_type_reference_source,
)
from .cancer_types import (
    cohort_aggregates,
    cohort_registry_df,
    resolve_cancer_type,
    resolve_cohort_id,
)
from .expression_builders import (
    PERCENTILE_BREAKPOINTS as _PERCENTILE_BREAKPOINTS,
)
from .expression_builders import (
    WITHIN_SAMPLE_THRESHOLDS as _WITHIN_SAMPLE_THRESHOLD_COLS,
)
from .expression_engine import id_columns, sample_columns, validate_nonnegative_expression
from .gene_ids import (
    canonical_gene_symbol,
    cdna_identical_groups,
    ensembl_id_alias_symbols,
    is_ensembl_gene_id,
    proteoform_collapse_overrides,
    resolve_ensembl_id,
    unversioned,
)
from .load_dataset import _register_derived_cache
from .load_dataset import get_data as get_data
from .normalization import clean_tpm, tpm_to_housekeeping_normalized
from .representative_partitions import (
    REPRESENTATIVE_PARTITION_POLICY_VERSION as REPRESENTATIVE_PARTITION_POLICY_VERSION,
)
from .version import DATA_VERSION, SOURCE_MATRIX_VERSION


@dataclass(frozen=True)
class ShardDataset:
    """A per-cohort summary artifact: one parquet shard per cohort code, available at
    gene level and (optionally) proteoform level. The single declarative record of
    everything the read path needs for one artifact, so the accessors stay generic over
    it rather than special-casing each.

    Fields:
      - ``noun`` — human label used in error messages.
      - ``gene_dir`` — bundle subdirectory of the gene-level shards.
      - ``gene_fetches`` — whether the gene-level shards ship in a released bundle (so a
        missing one may be auto-fetched) vs. must be recomputed on the fly.
      - ``proteoform_stem`` — the *stem* of the proteoform shard directory; ``None`` if
        the artifact has no proteoform variant. The on-disk directory is **scope-suffixed**
        (``f"{proteoform_stem}-{scope}"``) because identical-protein members group
        differently under ``"cta"`` vs ``"genome"``, so each scope is its own shard set.
      - ``proteoform_fetch_scopes`` — proteoform scopes whose shards ship in the active
        bundle. Other scopes still recompute on demand from the per-sample matrix
        because proteoform groups are scope-specific.
      - ``build_attr`` — name of the :mod:`oncoref.expression_builders` core that regenerates a
        missing shard from the per-sample matrix (the same core that produced the shipped
        shards, so on-the-fly and shipped values agree).

    The registry of the concrete artifacts is :data:`SHARD_DATASETS`.
    """

    noun: str
    gene_dir: str
    gene_fetches: bool
    proteoform_stem: str | None = None
    proteoform_fetch_scopes: tuple[str, ...] = ()
    build_attr: str | None = None

    def subdir(self, *, proteoform: bool, scope: str = "cta") -> str:
        """Bundle subdirectory holding this artifact's shards at the requested level.
        Gene-level is scope-independent; the proteoform variant is **scope-specific**
        (``f"{proteoform_stem}-{scope}"``)."""
        if not proteoform:
            return self.gene_dir
        if self.proteoform_stem is None:
            raise ValueError(f"{self.noun} has no proteoform variant")
        return f"{self.proteoform_stem}-{scope}"

    def fetches(self, *, proteoform: bool, scope: str = "cta") -> bool:
        """Whether the requested level's shards ship in a released bundle (so a missing
        shard may be auto-fetched rather than only recomputed)."""
        if not proteoform:
            return self.gene_fetches
        return scope in self.proteoform_fetch_scopes


#: The expression summary artifacts, keyed by short name. The QC-policy data bundle ships
#: gene-level percentiles, representatives, within-sample prevalence, and CTA-scope
#: proteoform percentile/within-sample shards. Other proteoform scopes still recompute
#: on demand from the source matrix.
SHARD_DATASETS: dict[str, ShardDataset] = {
    "representatives": ShardDataset(
        noun="representative-samples shard",
        gene_dir="cancer-reference-expression-representatives",
        gene_fetches=True,
    ),
    "percentiles": ShardDataset(
        noun="percentile vector",
        gene_dir="cancer-reference-expression-percentiles",
        gene_fetches=True,
        proteoform_stem="cancer-reference-expression-percentiles-proteoform",
        proteoform_fetch_scopes=("cta",),
        build_attr="cohort_percentile_vectors",
    ),
    "within_sample": ShardDataset(
        noun="within-sample top-fraction vector",
        gene_dir="cancer-reference-expression-within-sample-top5",
        gene_fetches=True,
        proteoform_stem="cancer-reference-expression-within-sample-top5-proteoform",
        proteoform_fetch_scopes=("cta",),
        build_attr="within_sample_top_fractions",
    ),
}
_REPRESENTATIVES = SHARD_DATASETS["representatives"]
_PERCENTILES = SHARD_DATASETS["percentiles"]
_WITHIN_SAMPLE = SHARD_DATASETS["within_sample"]

REPRESENTATIVE_ARTIFACT_SCHEMA_VERSION = "representative_expression_v2"
PERCENTILE_ARTIFACT_SCHEMA_VERSION = "cohort_percentile_expression_v1"
REFERENCE_EXPRESSION_SCHEMA_VERSION = "cancer_reference_expression_v5"
REPRESENTATIVE_SELECTION_METHOD = "central_medoid_then_farthest_first"
REPRESENTATIVE_SELECTION_BASIS = "biological_clean_tpm_log1p_distance"


# How many cleaned per-sample matrices to keep in the in-process LRU. Each frame is
# a full gene x sample matrix (~100MB+), so the default is intentionally small. A
# workflow that pools the same N>2 cohorts repeatedly (e.g. a gene then a proteoform
# pool over one cohort set) can raise CANCERDATA_PER_SAMPLE_CACHE to keep them all
# warm and skip the re-read, trading memory for latency.
def _per_sample_cache_size(default: int = 2) -> int:
    """Parse the cache-size env knob, falling back to ``default`` on any malformed
    value (empty string, non-integer) — a tuning knob must never break ``import``."""
    try:
        return max(1, int(os.environ.get("CANCERDATA_PER_SAMPLE_CACHE", str(default))))
    except (TypeError, ValueError):
        return default


_PER_SAMPLE_CACHE_SIZE = _per_sample_cache_size()


def _bundle_subdir(name: str, *, auto_fetch: bool = True) -> Path:
    """Locate a bundle shard directory: an in-repo checkout (``oncoref/data/…``)
    wins, else the downloaded bundle cache; the bundle is fetched if absent.

    ``auto_fetch=False`` skips the potentially large bundle download; the
    returned path simply won't exist when the shard is absent locally."""
    local = data_bundle.find_local_item(name)
    if local is not None:
        return local
    if auto_fetch:
        data_bundle.ensure_local()
    return data_bundle.cache_dir() / name


def _available_shard_codes(root: Path) -> list[str]:
    """Sorted cohort codes that ship a parquet shard under ``root``."""
    if not root.exists():
        return []
    return sorted(p.stem for p in root.glob("*.parquet") if not p.name.startswith("._"))


def _local_shard_dirs(
    dataset: ShardDataset, *, proteoform: bool = False, scope: str = "cta"
) -> tuple[Path, ...]:
    """Every locally present directory for one shard variant, without fetching."""
    return data_bundle.local_item_paths(dataset.subdir(proteoform=proteoform, scope=scope))


def _available_local_cohorts(
    dataset: ShardDataset, *, proteoform: bool = False, scope: str = "cta"
) -> list[str]:
    """Sorted union of shard codes across package and active-cache roots.

    Cohorts whose selected source has an invalid expression scale are excluded:
    their shards are withheld by :func:`_shard_path`.
    """
    invalid = _invalid_selected_source_codes()
    return sorted(
        {
            code
            for root in _local_shard_dirs(dataset, proteoform=proteoform, scope=scope)
            for code in _available_shard_codes(root)
            if code not in invalid
        }
    )


def _shard_dir(dataset: ShardDataset, *, proteoform: bool = False, scope: str = "cta") -> Path:
    """The bundle directory holding ``dataset``'s per-cohort shards at the requested
    level and ``scope`` (the proteoform variant is scope-specific — see
    :meth:`ShardDataset.subdir`). The per-variant auto-fetch policy lives on the dataset
    record, so it's decided in one place rather than re-stated at each call site."""
    return _bundle_subdir(
        dataset.subdir(proteoform=proteoform, scope=scope),
        auto_fetch=dataset.fetches(proteoform=proteoform, scope=scope),
    )


def _available_cohorts(
    dataset: ShardDataset, *, proteoform: bool = False, scope: str = "cta"
) -> list[str]:
    """Sorted cohort codes with a shipped shard for ``dataset`` (gene, or the
    scope-specific proteoform variant)."""
    local = _available_local_cohorts(dataset, proteoform=proteoform, scope=scope)
    if local:
        return local
    invalid = _invalid_selected_source_codes()
    return [
        code
        for code in _available_shard_codes(_shard_dir(dataset, proteoform=proteoform, scope=scope))
        if code not in invalid
    ]


def _shard_path(
    dataset: ShardDataset,
    code: str,
    *,
    proteoform: bool = False,
    scope: str = "cta",
) -> Path:
    _reject_invalid_source_code(code, context=f"{dataset.noun} shard")
    filename = f"{code}.parquet"
    for root in _local_shard_dirs(dataset, proteoform=proteoform, scope=scope):
        candidate = root / filename
        if candidate.is_file():
            return candidate
    return _shard_dir(dataset, proteoform=proteoform, scope=scope) / filename


def _resolve_cancer_types(
    cancer_types: str | Iterable[str] | None,
    *,
    expand_aggregates: bool = False,
) -> list[str] | None:
    """Resolve a code / alias / iterable to canonical codes. With
    ``expand_aggregates``, a computed-aggregate code (``SARC`` and the
    ``SARC_RMS`` / ``SARC_LPS`` rollups) expands to its member subtypes."""
    if cancer_types is None:
        return None
    if isinstance(cancer_types, str):
        requested = [cancer_types]
    else:
        requested = list(cancer_types)
    if not expand_aggregates:
        return [resolve_cancer_type(code) for code in requested]

    aggregates = cohort_aggregates()
    out: list[str] = []
    for code in requested:
        members = aggregates.get(str(code))
        if members is None:
            resolved = resolve_cancer_type(code)
            members = aggregates.get(resolved)
            if members is None:
                out.append(resolved)
                continue
        out.extend(members)
    return list(dict.fromkeys(out))


def _gene_filter_mask(df: pd.DataFrame, genes: str | Iterable[str] | None) -> pd.Series:
    if genes is None:
        return pd.Series(True, index=df.index)
    wanted = {genes} if isinstance(genes, str) else set(genes)
    wanted = {str(g) for g in wanted}
    wanted_ids = set(wanted)
    wanted_unversioned = {unversioned(g) for g in wanted}
    for gene in wanted:
        base = unversioned(gene)
        if is_ensembl_gene_id(base):
            resolved = resolve_ensembl_id(base)
            wanted_ids.add(resolved)
            wanted_unversioned.add(unversioned(resolved))
    ids = df["Ensembl_Gene_ID"].astype(str)
    return (
        ids.isin(wanted_ids)
        | ids.map(unversioned).isin(wanted_unversioned)
        | df["Symbol"].astype(str).isin(wanted)
    )


_REPRESENTATIVE_PROVENANCE_COLUMNS = [
    "source_cohort",
    "source_version",
    "source_project",
    "source_diagnosis",
    "source_morphology",
    "source_sample",
    "source_group_id",
    "n_cohort_samples",
    "sample_qc",
    "sample_qc_requested",
    "source_sample_qc",
    "sample_qc_effective",
    "sample_qc_fallback_reason",
    "sample_qc_policy_version",
    "source_sample_qc_reasons",
    "n_qc_pass",
    "n_qc_warn",
    "n_qc_fail",
    "source_scale_class",
    "linear_tpm_comparable",
    "recommended_for_absolute_tpm_floor",
    "selection_scale_class",
    "representative_role",
    "benchmark_eligible",
    "partition_role",
    "partition_reason",
    "partition_policy_version",
    "partition_validation_target",
    "partition_status",
    "n_partition_train",
    "n_partition_validation",
    "n_partition_validation_external",
    "n_partition_audit_only",
    "review_source",
    "review_note",
    "selection_rank",
    "selection_method",
    "selection_basis",
    "artifact_schema_version",
    "data_version",
    "source_matrix_version",
]

_REPRESENTATIVE_AVAILABILITY_COLUMNS = [
    "cancer_code",
    "n_representatives",
    "source_cohort",
    "source_scale_class",
    "linear_tpm_comparable",
    "recommended_for_absolute_tpm_floor",
    "sample_qc_requested",
    "sample_qc_effective",
    "sample_qc_fallback_reason",
    "sample_qc_policy_version",
    "n_qc_pass",
    "n_qc_warn",
    "n_qc_fail",
    "representative_role",
    "benchmark_eligible",
    "partition_policy_version",
    "partition_validation_target",
    "partition_status",
    "n_partition_train",
    "n_partition_validation",
    "n_partition_validation_external",
    "n_partition_audit_only",
    "selection_scale_class",
    "selection_method",
    "selection_basis",
    "availability_reason",
]

_REPRESENTATIVE_PARTITION_COLUMNS = [
    "cancer_code",
    "representative_id",
    "source_group_id",
    "source_cohort",
    "source_project",
    "source_sample",
    "representative_role",
    "benchmark_eligible",
    "partition_role",
    "partition_reason",
    "partition_policy_version",
    "partition_validation_target",
    "partition_status",
    "n_partition_train",
    "n_partition_validation",
    "n_partition_validation_external",
    "n_partition_audit_only",
]


_INTERNAL_REPRESENTATIVE_ID_RE = re.compile(r"^(?P<code>.+)__rep(?P<rank>\d+)$")
_PIRLYGENES_REPRESENTATIVE_ID_RE = re.compile(r"^(?P<code>.+)_rep(?P<rank>\d+)$")


def _representative_rank(representative_id: object) -> int | None:
    text = str(representative_id)
    for pattern in (_INTERNAL_REPRESENTATIVE_ID_RE, _PIRLYGENES_REPRESENTATIVE_ID_RE):
        match = pattern.match(text)
        if match:
            return int(match.group("rank"))
    return None


def _representative_id_for_style(representative_id: object, *, style: str) -> str:
    _validate_representative_id_style(style)
    text = str(representative_id)
    if style == "internal":
        return text
    match = _INTERNAL_REPRESENTATIVE_ID_RE.match(text)
    if not match:
        return text
    return f"{match.group('code')}_rep{int(match.group('rank')):02d}"


def _validate_representative_id_style(style: str) -> None:
    if style not in {"pirlygenes", "internal"}:
        raise ValueError("representative_id_style must be 'pirlygenes' or 'internal'")


def _validate_gene_id_style(style: str) -> None:
    if style not in {"oncoref", "pirlygenes"}:
        raise ValueError("gene_id_style must be 'oncoref' or 'pirlygenes'")


def _artifact_legacy_gene_id_map(
    *,
    product: str,
    cancer_codes: Iterable[str],
) -> tuple[dict[str, str], dict[str, str]]:
    code_set = {str(c) for c in cancer_codes}
    if not code_set:
        return {}, {}
    deltas = expression_artifact_gene_universe_deltas(
        product=product,
        delta_kind="pirlygenes_only",
    )
    deltas = deltas[deltas["status"].isin(_ARTIFACT_CANONICALIZATION_STATUSES)]

    def _matches(cell: object) -> bool:
        codes = {c for c in str(cell or "").split(";") if c}
        return bool(codes & code_set)

    deltas = deltas[deltas["cancer_code"].map(_matches)]
    if deltas.empty:
        return {}, {}
    deltas = deltas.drop_duplicates("oncoref_ensembl_gene_id", keep="first")
    id_map = dict(zip(deltas["oncoref_ensembl_gene_id"], deltas["legacy_ensembl_gene_id"]))
    symbol_map = dict(zip(deltas["oncoref_ensembl_gene_id"], deltas["symbol"]))
    return id_map, symbol_map


def _artifact_current_only_gene_ids(
    *,
    product: str,
    cancer_codes: Iterable[str],
) -> set[str]:
    code_set = {str(c) for c in cancer_codes}
    if not code_set:
        return set()
    deltas = expression_artifact_gene_universe_deltas(
        product=product,
        delta_kind="oncoref_only",
    )

    def _matches(cell: object) -> bool:
        codes = {c for c in str(cell or "").split(";") if c}
        return bool(codes & code_set)

    deltas = deltas[deltas["cancer_code"].map(_matches)]
    return {
        gene_id
        for gene_id in deltas["oncoref_ensembl_gene_id"].map(_delta_nonempty)
        if gene_id is not None
    }


def _apply_gene_id_style(
    df: pd.DataFrame,
    *,
    product: str,
    cancer_codes: Iterable[str],
    gene_id_style: str,
    alias_expand_remaps: bool = False,
) -> pd.DataFrame:
    _validate_gene_id_style(gene_id_style)
    if gene_id_style == "oncoref" or "Ensembl_Gene_ID" not in df.columns:
        return df
    id_map, symbol_map = _artifact_legacy_gene_id_map(
        product=product,
        cancer_codes=cancer_codes,
    )
    if not id_map:
        return df
    out = df.copy()
    current = out["Ensembl_Gene_ID"].astype(str)
    current_only_ids = _artifact_current_only_gene_ids(
        product=product,
        cancer_codes=cancer_codes,
    )
    replacement_id_map = {
        gene_id: legacy for gene_id, legacy in id_map.items() if gene_id in current_only_ids
    }
    duplicate_id_map = (
        {gene_id: legacy for gene_id, legacy in id_map.items() if gene_id not in current_only_ids}
        if alias_expand_remaps
        else {}
    )
    replacement_symbol_map = {
        gene_id: symbol for gene_id, symbol in symbol_map.items() if gene_id in replacement_id_map
    }
    duplicate_symbol_map = {
        gene_id: symbol for gene_id, symbol in symbol_map.items() if gene_id in duplicate_id_map
    }

    legacy_ids = current.map(replacement_id_map)
    mapped = legacy_ids.notna()
    out.loc[mapped, "Ensembl_Gene_ID"] = legacy_ids[mapped].to_numpy()
    if "Symbol" in out.columns:
        legacy_symbols = current.map(replacement_symbol_map)
        out.loc[mapped, "Symbol"] = legacy_symbols[mapped].to_numpy()
    duplicate_legacy_ids = current.map(duplicate_id_map)
    duplicate_mask = duplicate_legacy_ids.notna()
    if duplicate_mask.any():
        existing_ids = set(out["Ensembl_Gene_ID"].astype(str))
        duplicate_mask &= ~duplicate_legacy_ids.isin(existing_ids)
    if duplicate_mask.any():
        duplicates = out.loc[duplicate_mask].copy()
        duplicates.loc[:, "Ensembl_Gene_ID"] = duplicate_legacy_ids[duplicate_mask].to_numpy()
        if "Symbol" in duplicates.columns:
            duplicate_symbols = current.map(duplicate_symbol_map)
            duplicates.loc[:, "Symbol"] = duplicate_symbols[duplicate_mask].to_numpy()
        out = pd.concat([out, duplicates], ignore_index=True)
    return out


def _representative_attrs(
    *,
    codes: Iterable[str],
    normalize: str,
    format: str,
    k: int | None,
    representative_id_style: str,
    gene_id_style: str,
    gene_universe: str,
    sample_qc: str,
    artifact_qc_meta: pd.DataFrame | None = None,
) -> dict[str, object]:
    attrs = {
        "artifact": "representative_samples",
        "schema_version": REPRESENTATIVE_ARTIFACT_SCHEMA_VERSION,
        "data_version": DATA_VERSION,
        "source_matrix_version": SOURCE_MATRIX_VERSION,
        "cancer_codes": list(codes),
        "normalize": normalize,
        "format": format,
        "k": k,
        "representative_id_style": representative_id_style,
        "gene_id_style": gene_id_style,
        "gene_universe": gene_universe,
        "selection_method": REPRESENTATIVE_SELECTION_METHOD,
        "selection_basis": REPRESENTATIVE_SELECTION_BASIS,
    }
    if artifact_qc_meta is not None:
        attrs.update(_artifact_qc_attrs(artifact_qc_meta, requested_sample_qc=sample_qc))
    availability = representative_cohort_availability(codes)
    if not availability.empty:
        scale_classes = sorted(
            {str(v) for v in availability["source_scale_class"] if pd.notna(v) and str(v)}
        )
        attrs["source_scale_class"] = (
            scale_classes[0]
            if len(scale_classes) == 1
            else ("mixed" if scale_classes else "unknown")
        )
        for col in (
            "linear_tpm_comparable",
            "recommended_for_absolute_tpm_floor",
            "benchmark_eligible",
        ):
            values = [_optional_bool(v) for v in availability[col]]
            known = [v for v in values if v is not None]
            attrs[col] = bool(known and len(known) == len(values) and all(known))
    return attrs


def _representative_empty_frame(*, include_provenance: bool) -> pd.DataFrame:
    cols = ["Ensembl_Gene_ID", "Symbol", "cancer_code", "representative_id", "expression"]
    if include_provenance:
        cols.extend(_REPRESENTATIVE_PROVENANCE_COLUMNS)
    return pd.DataFrame(columns=cols)


def _optional_bool(value: object) -> bool | None:
    if pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def _representative_provenance(root: Path) -> pd.DataFrame:
    prov_path = root / "_provenance.csv"
    if not prov_path.exists():
        return pd.DataFrame(columns=["representative_id", "cancer_code"])
    prov = pd.read_csv(prov_path)
    if "cancer_code" not in prov.columns:
        prov["cancer_code"] = (
            prov["representative_id"].astype(str).str.replace(r"__rep\d+$", "", regex=True)
        )
    enrich_cols = {
        "sample_qc_status": "source_sample_qc",
        "sample_qc_reasons": "source_sample_qc_reasons",
        "source_scale_class": "source_scale_class",
        "linear_tpm_comparable": "linear_tpm_comparable",
        "recommended_for_absolute_tpm_floor": "recommended_for_absolute_tpm_floor",
    }
    needs_enrichment = any(
        col not in prov.columns or prov[col].isna().any() for col in enrich_cols.values()
    )
    if needs_enrichment and "source_sample" in prov.columns:
        qc = source_matrix_sample_qc_manifest(
            prov["cancer_code"].dropna().astype(str).unique(),
            auto_fetch=False,
        )
        if not qc.empty:
            qc_keep = ["cancer_code", "sample_id", *enrich_cols]
            qc = qc[[col for col in qc_keep if col in qc.columns]].rename(
                columns={"sample_id": "source_sample", **enrich_cols}
            )
            prov = prov.merge(
                qc.drop_duplicates(["cancer_code", "source_sample"]),
                on=["cancer_code", "source_sample"],
                how="left",
                suffixes=("", "__qc"),
            )
            for col in enrich_cols.values():
                qc_col = f"{col}__qc"
                if qc_col in prov.columns:
                    if col in prov.columns:
                        prov[col] = prov[col].where(prov[col].notna(), prov[qc_col])
                    else:
                        prov[col] = prov[qc_col]
                    prov = prov.drop(columns=qc_col)
    if "selection_scale_class" not in prov.columns:
        prov["selection_scale_class"] = prov.get("source_scale_class", pd.NA)
    return prov


def _attach_representative_provenance(long: pd.DataFrame, root: Path) -> pd.DataFrame:
    prov = _representative_provenance(root)
    if not prov.empty:
        keep = ["representative_id", *_REPRESENTATIVE_PROVENANCE_COLUMNS]
        long = long.merge(
            prov[[c for c in keep if c in prov.columns]],
            on="representative_id",
            how="left",
        )
    for col in (
        "source_cohort",
        "source_version",
        "source_project",
        "source_diagnosis",
        "source_morphology",
        "source_sample",
        "source_group_id",
        "representative_role",
        "partition_role",
        "partition_reason",
        "partition_policy_version",
        "partition_status",
        "review_source",
        "review_note",
    ):
        if col not in long.columns:
            long[col] = pd.NA
    for col in (
        "n_cohort_samples",
        "selection_rank",
        "n_qc_pass",
        "n_qc_warn",
        "n_qc_fail",
        "partition_validation_target",
        "n_partition_train",
        "n_partition_validation",
        "n_partition_validation_external",
        "n_partition_audit_only",
    ):
        if col not in long.columns:
            long[col] = pd.NA
    for col in (
        "sample_qc",
        "sample_qc_requested",
        "source_sample_qc",
        "sample_qc_effective",
        "sample_qc_fallback_reason",
        "sample_qc_policy_version",
        "source_sample_qc_reasons",
        "source_scale_class",
        "linear_tpm_comparable",
        "recommended_for_absolute_tpm_floor",
        "selection_scale_class",
        "benchmark_eligible",
    ):
        if col not in long.columns:
            long[col] = pd.NA
    if "selection_rank" in long.columns:
        ranks = pd.Series(long["representative_id"].map(_representative_rank), index=long.index)
        long["selection_rank"] = (
            long["selection_rank"].astype("Int64").fillna(ranks).astype("Int64")
        )
    long["selection_method"] = REPRESENTATIVE_SELECTION_METHOD
    long["selection_basis"] = REPRESENTATIVE_SELECTION_BASIS
    long["artifact_schema_version"] = REPRESENTATIVE_ARTIFACT_SCHEMA_VERSION
    long["data_version"] = DATA_VERSION
    long["source_matrix_version"] = SOURCE_MATRIX_VERSION
    return long


def _one_provenance_value(group: pd.DataFrame, col: str, default=pd.NA):
    if col not in group:
        return default
    values = [v for v in group[col] if pd.notna(v) and str(v) != ""]
    unique = list(dict.fromkeys(values))
    return unique[0] if len(unique) == 1 else ("mixed" if unique else default)


def _all_provenance_values_true(group: pd.DataFrame, col: str):
    if col not in group or group.empty:
        return pd.NA
    values = [_optional_bool(v) for v in group[col]]
    if any(v is None for v in values):
        return pd.NA
    return all(values)


def _percentile_cols() -> list[str]:
    return [f"p{bp}" for bp in _PERCENTILE_BREAKPOINTS]


def _percentile_identity_cols(*, proteoform: bool) -> list[str]:
    if proteoform:
        return ["proteoform_key", "Ensembl_Gene_ID", "Symbol", "proteoform_members"]
    return ["Ensembl_Gene_ID", "Symbol"]


def _percentile_provenance_columns() -> list[str]:
    return [
        "cancer_code",
        "normalization",
        "expression_unit",
        "percentile_basis",
        "artifact_schema_version",
        "data_version",
        "source_matrix_version",
        "sample_qc",
        "sample_qc_policy_version",
    ]


def _percentile_attrs(
    *,
    code: str,
    as_tpm: bool,
    proteoform: bool,
    scope: str,
    source: str,
    gene_id_style: str,
    gene_universe: str,
    sample_qc: str,
    artifact_qc_meta: pd.DataFrame | None = None,
    missing_reason: str | None = None,
) -> dict[str, object]:
    attrs: dict[str, object] = {
        "artifact": "cohort_gene_percentiles",
        "schema_version": PERCENTILE_ARTIFACT_SCHEMA_VERSION,
        "data_version": DATA_VERSION,
        "source_matrix_version": SOURCE_MATRIX_VERSION,
        "cancer_code": code,
        "as_tpm": as_tpm,
        "proteoform": proteoform,
        "scope": scope,
        "source": source,
        "gene_id_style": gene_id_style,
        "gene_universe": gene_universe,
        "normalization": "tpm_clean" if as_tpm else "tpm_clean_log1p",
        "expression_unit": "tpm_clean" if as_tpm else "log1p_tpm_clean",
        "percentile_basis": "biological_clean_tpm_across_samples",
    }
    if artifact_qc_meta is not None:
        attrs.update(_artifact_qc_attrs(artifact_qc_meta, requested_sample_qc=sample_qc))
    if missing_reason is not None:
        attrs["missing_reason"] = missing_reason
    return attrs


def _empty_percentile_frame(
    *,
    code: str,
    as_tpm: bool,
    proteoform: bool,
    scope: str,
    include_provenance: bool,
    gene_id_style: str,
    gene_universe: str,
    sample_qc: str,
    artifact_qc_meta: pd.DataFrame | None = None,
    include_gene_universe_flags: bool,
    missing_reason: str,
) -> pd.DataFrame:
    cols = [*_percentile_identity_cols(proteoform=proteoform), *_percentile_cols()]
    if include_provenance:
        cols.extend(_percentile_provenance_columns())
    if include_gene_universe_flags:
        cols.extend(_ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS)
    out = pd.DataFrame(columns=cols)
    out.attrs.update(
        _percentile_attrs(
            code=code,
            as_tpm=as_tpm,
            proteoform=proteoform,
            scope=scope,
            source="missing",
            gene_id_style=gene_id_style,
            gene_universe=gene_universe,
            sample_qc=sample_qc,
            artifact_qc_meta=artifact_qc_meta,
            missing_reason=missing_reason,
        )
    )
    return out


def _attach_percentile_provenance(
    df: pd.DataFrame,
    *,
    code: str,
    as_tpm: bool,
    sample_qc: str,
    artifact_qc_meta: pd.DataFrame,
) -> pd.DataFrame:
    out = df.copy()
    out["cancer_code"] = code
    out["normalization"] = "tpm_clean" if as_tpm else "tpm_clean_log1p"
    out["expression_unit"] = "tpm_clean" if as_tpm else "log1p_tpm_clean"
    out["percentile_basis"] = "biological_clean_tpm_across_samples"
    out["artifact_schema_version"] = PERCENTILE_ARTIFACT_SCHEMA_VERSION
    out["data_version"] = DATA_VERSION
    out["source_matrix_version"] = SOURCE_MATRIX_VERSION
    qc_values = sorted(
        {_artifact_build_metadata_qc_value(row) for _, row in artifact_qc_meta.iterrows()} - {""}
    )
    out["sample_qc"] = ",".join(qc_values) if qc_values else sample_qc
    policy_versions = sorted(
        {
            str(v)
            for v in artifact_qc_meta.get("sample_qc_policy_version", pd.Series(dtype=object))
            if pd.notna(v) and str(v)
        }
    )
    out["sample_qc_policy_version"] = ",".join(policy_versions) if policy_versions else pd.NA
    return out


def representative_cohort_availability(
    cancer_types: str | Iterable[str] | None = None,
    *,
    linear_tpm_comparable: bool | None = None,
    benchmark_eligible: bool | None = None,
) -> pd.DataFrame:
    """One row per shipped representative cohort with classifier compatibility.

    Proxy-scale cohorts remain available for rank/percentile consumers, while
    ``linear_tpm_comparable`` and ``benchmark_eligible`` let classifiers fail closed
    without privately joining the source-matrix QC manifest. Boolean filters retain
    only cohorts whose released provenance explicitly matches the requested value.
    """
    root = _shard_dir(_REPRESENTATIVES)
    shipped = _available_shard_codes(root)
    if cancer_types is None:
        codes = shipped
    else:
        raw = [cancer_types] if isinstance(cancer_types, str) else list(cancer_types)
        aggregates = cohort_aggregates()
        requested = []
        for value in raw:
            text = str(value)
            if text in shipped:
                requested.append(text)
                continue
            resolved = resolve_cancer_type(value, strict=False)
            if resolved is not None:
                requested.extend(aggregates.get(resolved, [resolved]))
        codes = [code for code in dict.fromkeys(requested) if code in shipped]
    provenance = _representative_provenance(root)
    rows = []
    for code in codes:
        group = provenance[provenance.get("cancer_code", pd.Series(dtype=str)) == code]
        comparable = _all_provenance_values_true(group, "linear_tpm_comparable")
        benchmark = _all_provenance_values_true(group, "benchmark_eligible")
        role = _one_provenance_value(group, "representative_role")
        fallback_reason = _one_provenance_value(group, "sample_qc_fallback_reason")
        if benchmark is False:
            availability_reason = str(fallback_reason if pd.notna(fallback_reason) else role)
        elif comparable is False:
            availability_reason = "nonlinear_or_proxy_expression_scale"
        elif benchmark is True and comparable is True:
            availability_reason = "available"
        else:
            availability_reason = "provenance_unavailable"
        rows.append(
            {
                "cancer_code": code,
                "n_representatives": (
                    int(group["representative_id"].nunique())
                    if not group.empty and "representative_id" in group
                    else pd.NA
                ),
                "source_cohort": _one_provenance_value(group, "source_cohort"),
                "source_scale_class": _one_provenance_value(group, "source_scale_class"),
                "linear_tpm_comparable": comparable,
                "recommended_for_absolute_tpm_floor": _all_provenance_values_true(
                    group, "recommended_for_absolute_tpm_floor"
                ),
                "sample_qc_requested": _one_provenance_value(group, "sample_qc_requested"),
                "sample_qc_effective": _one_provenance_value(group, "sample_qc_effective"),
                "sample_qc_fallback_reason": fallback_reason,
                "sample_qc_policy_version": _one_provenance_value(
                    group, "sample_qc_policy_version"
                ),
                "n_qc_pass": _one_provenance_value(group, "n_qc_pass"),
                "n_qc_warn": _one_provenance_value(group, "n_qc_warn"),
                "n_qc_fail": _one_provenance_value(group, "n_qc_fail"),
                "representative_role": role,
                "benchmark_eligible": benchmark,
                "partition_policy_version": _one_provenance_value(
                    group, "partition_policy_version"
                ),
                "partition_validation_target": _one_provenance_value(
                    group, "partition_validation_target"
                ),
                "partition_status": _one_provenance_value(group, "partition_status"),
                "n_partition_train": _one_provenance_value(group, "n_partition_train"),
                "n_partition_validation": _one_provenance_value(group, "n_partition_validation"),
                "n_partition_validation_external": _one_provenance_value(
                    group, "n_partition_validation_external"
                ),
                "n_partition_audit_only": _one_provenance_value(group, "n_partition_audit_only"),
                "selection_scale_class": _one_provenance_value(group, "selection_scale_class"),
                "selection_method": REPRESENTATIVE_SELECTION_METHOD,
                "selection_basis": REPRESENTATIVE_SELECTION_BASIS,
                "availability_reason": availability_reason,
            }
        )
    out = pd.DataFrame(rows, columns=_REPRESENTATIVE_AVAILABILITY_COLUMNS)
    for col, expected in (
        ("linear_tpm_comparable", linear_tpm_comparable),
        ("benchmark_eligible", benchmark_eligible),
    ):
        if expected is not None:
            out = out[out[col].map(_optional_bool) == bool(expected)].copy()
    out.attrs["schema_version"] = REPRESENTATIVE_ARTIFACT_SCHEMA_VERSION
    out.attrs["data_version"] = DATA_VERSION
    out.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    return out.reset_index(drop=True)


def representative_partition_manifest(
    cancer_types: str | Iterable[str] | None = None,
    *,
    representative_id_style: str = "pirlygenes",
) -> pd.DataFrame:
    """One machine-readable train/validation row per representative vector.

    ``source_group_id`` is the leakage boundary: every parent, subtype, or alias
    row backed by the same physical sample has one shared partition role.
    ``validation_external`` holds out a complete independent source project;
    ``audit_only`` keeps reviewed or QC-fallback vectors out of both model fitting
    and evaluation. Cohorts without two independent eligible groups remain
    ``train`` and report ``partition_status="insufficient_independent_groups"``.
    """

    _validate_representative_id_style(representative_id_style)
    root = _shard_dir(_REPRESENTATIVES)
    provenance = _representative_provenance(root)

    if cancer_types is not None and not provenance.empty:
        requested = set(_resolve_cancer_types(cancer_types, expand_aggregates=True))
        provenance = provenance[provenance["cancer_code"].isin(requested)].copy()

    for column in _REPRESENTATIVE_PARTITION_COLUMNS:
        if column not in provenance:
            provenance[column] = pd.NA
    out = provenance[_REPRESENTATIVE_PARTITION_COLUMNS].copy()
    out["representative_id"] = out["representative_id"].map(
        lambda value: _representative_id_for_style(
            value,
            style=representative_id_style,
        )
    )
    out = out.sort_values(["cancer_code", "representative_id"], kind="stable")
    policy_versions = [
        str(value) for value in out["partition_policy_version"].dropna().unique() if str(value)
    ]
    out.attrs["schema_version"] = REPRESENTATIVE_ARTIFACT_SCHEMA_VERSION
    out.attrs["partition_policy_version"] = (
        policy_versions[0] if len(policy_versions) == 1 else ("mixed" if policy_versions else None)
    )
    out.attrs["data_version"] = DATA_VERSION
    out.attrs["source_matrix_version"] = SOURCE_MATRIX_VERSION
    return out.reset_index(drop=True)


def available_representative_cohorts(
    *,
    linear_tpm_comparable: bool | None = None,
    benchmark_eligible: bool | None = None,
) -> list[str]:
    """Sorted cohort codes with a representative shard.

    Optional compatibility filters are backed by released representative provenance;
    for example, both ``True`` filters return classifier-ready linear-TPM cohorts.
    """
    if linear_tpm_comparable is None and benchmark_eligible is None:
        return _available_cohorts(_REPRESENTATIVES)
    availability = representative_cohort_availability(
        linear_tpm_comparable=linear_tpm_comparable,
        benchmark_eligible=benchmark_eligible,
    )
    return availability["cancer_code"].astype(str).tolist()


def available_percentile_cohorts(*, proteoform: bool = False, scope: str = "cta") -> list[str]:
    """Cohort codes that ship a percentile-vector shard (sorted). With
    ``proteoform=True``, the proteoform-summed variant (one vector per proteoform
    key, identical-protein members collapsed before ranking, in ``scope``)."""
    return _available_cohorts(_PERCENTILES, proteoform=proteoform, scope=scope)


def _locally_available_summary_cohorts(
    dataset: ShardDataset,
    *,
    proteoform: bool,
    scope: str,
    include_recomputable: bool,
) -> list[str]:
    """Local shard union, optionally extended by cached source matrices."""
    available = set(_available_local_cohorts(dataset, proteoform=proteoform, scope=scope))
    if proteoform:
        from .proteoforms import _proteoform_registry_metadata

        expected = _proteoform_registry_metadata(scope)
        current = set()
        for code in available:
            attrs = pd.read_parquet(
                _shard_path(dataset, code, proteoform=True, scope=scope), columns=[]
            ).attrs
            if all(attrs.get(key) == value for key, value in expected.items()):
                current.add(code)
        available = current
    if include_recomputable:
        available.update(
            code for code in source_matrices.available_cohorts() if source_matrices.is_cached(code)
        )
    return sorted(available)


def locally_available_percentile_cohorts(
    *,
    proteoform: bool = False,
    scope: str = "cta",
    include_recomputable: bool = True,
) -> list[str]:
    """Cohorts readable as percentile vectors without a download.

    The result is the union of shards in package/artifact data and the active
    partial bundle cache. With ``include_recomputable=True`` (default), it also
    includes cohorts whose per-sample source matrix is cached locally. This
    function never fetches either data source.
    """
    return _locally_available_summary_cohorts(
        _PERCENTILES,
        proteoform=proteoform,
        scope=scope,
        include_recomputable=include_recomputable,
    )


def per_sample_expression(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    proteoform: bool = False,
    scope: str = "cta",
    sample_qc: str = "all",
) -> pd.DataFrame:
    """Full per-sample expression matrix (genes x **every** sample) for a cohort —
    the raw **TPM values** at gene level (default) or proteoform level.

    The packaged references are summaries — per-gene percentile vectors, bounded
    medoid :func:`representative_cohort_samples`, within-sample top fractions. This
    returns the raw material behind them: one column per individual sample or
    library. Columns are not necessarily independent patients. Join
    :func:`sample_metadata` by sample_id and group resolved participant_id values
    for patient-level analyses. It fetches the cohort's per-sample matrix via
    :mod:`oncoref.source_matrices` (a per-cohort release asset, tens of MB) and
    normalizes it.

    ``normalize``:
      - ``"tpm_clean"`` (default) — two-compartment clean TPM (the comparable
        biological view the summaries are built on);
      - ``"tpm_clean_log1p"`` — clean TPM, ``log1p``-transformed;
      - ``"tpm_clean_hk"`` — clean TPM divided per sample by the biological
        clean-TPM housekeeping-panel median-of-ratios size factor against the
        fixed HPA-derived reference profile (unit-free ratio-to-baseline,
        robust to library-depth drift);
      - ``"tpm_raw"`` — the matrix as shipped (raw TPM), no normalization.

    Every source read requires nonnegative linear TPM (NaN remains unmeasured).
    Invalid values raise before canonicalization or normalization; clipping a
    transformed source field would hide a unit error rather than repair it.

    With ``proteoform=True``, identical-protein paralogs are **summed per sample** to
    proteoform level (:func:`oncoref.proteoforms.collapse_to_proteoforms`, ``scope``
    = ``"cta"``/``"genome"``) — a **proteoform-level** frame carrying ``proteoform_key``
    (see :func:`oncoref.proteoforms.expression_level`). The sum is always taken in
    **linear** TPM and the ``log1p`` transform (if any) applied *after*, so the
    proteoform value is ``log1p(Σ member TPM)``, not the meaningless ``Σ log1p``.
    (``scope`` is ignored when ``proteoform=False``.)

    ``auto_fetch=False`` raises instead of downloading if the matrix isn't cached.
    ``sample_qc`` can be ``"all"`` (default), ``"pass"``, or ``"pass_or_warn"``; the
    latter two use :func:`sample_expression_qc` to drop failing sample columns while
    retaining the gene rows.
    Returns ``Ensembl_Gene_ID``, ``Symbol`` and one column per sample (plus the
    proteoform identity columns when collapsed).
    """
    if normalize not in _PER_SAMPLE_NORMALIZE:
        raise ValueError(f"normalize must be one of {_PER_SAMPLE_NORMALIZE}")
    sample_qc = _validate_sample_qc(sample_qc)
    code = resolve_cancer_type(cancer_type, strict=False) or cancer_type
    if auto_fetch:
        path = source_matrices.ensure(code)
    else:
        path = source_matrices.local_path(code)
        if not path.exists():
            raise FileNotFoundError(
                f"per-sample matrix for {code!r} not cached at {path}. "
                f"Run source_matrices.fetch({code!r}) to download it."
            )
    # The read + clean_tpm of a tens-of-MB matrix is the dominant cost on every
    # coverage / 9-mer call; memoize it (keyed on the matrix path + its mtime +
    # normalize) and hand callers a fresh copy so the shared frame can't be mutated.
    # The mtime in the key self-invalidates the cache if the matrix is re-fetched.
    mtime = os.path.getmtime(path)
    if not proteoform:
        out = _load_per_sample_matrix(str(path), mtime, normalize).copy()
        return _apply_sample_qc_filter(out, str(code), sample_qc=sample_qc, auto_fetch=auto_fetch)
    # Proteoform level: sum members in LINEAR TPM, then apply the requested transform.
    from .proteoforms import collapse_to_proteoforms

    linear = "tpm_raw" if normalize == "tpm_raw" else "tpm_clean"
    out = collapse_to_proteoforms(_load_per_sample_matrix(str(path), mtime, linear), scope=scope)
    samples = sample_columns(out)
    if normalize == "tpm_clean_log1p":
        out[samples] = np.log1p(out[samples].to_numpy(dtype=float))
    elif normalize == "tpm_clean_hk":
        out = _housekeeping_normalize(out, samples)
    return _apply_sample_qc_filter(out, str(code), sample_qc=sample_qc, auto_fetch=auto_fetch)


def sample_metadata(
    cancer_type, *, auto_fetch: bool = True, sample_qc: str = "all"
) -> pd.DataFrame:
    """Patient identity and source metadata aligned to actual expression columns.

    Uses the same cached matrix and QC selection as :func:`per_sample_expression`.
    Returns one row per retained sample in matrix-column order, with sample_id,
    cancer_code, source_cohort, physical sample_key, participant_id and identity
    provenance. Blank participant_id means unresolved; it must not count as an
    independent donor. Unknown identities can still supply single-sample evidence.
    """
    from .samples import sample_identities

    code = resolve_cancer_type(cancer_type, strict=False) or cancer_type
    matrix = per_sample_expression(
        code, normalize="tpm_raw", auto_fetch=auto_fetch, sample_qc=sample_qc
    )
    info = source_matrices.cohort_info(code)
    rows = pd.DataFrame(
        {
            "sample_id": sample_columns(matrix),
            "cancer_code": str(code),
            "source_cohort": info["source_cohort"],
        }
    )
    result = sample_identities(rows)
    result.attrs["sample_qc"] = sample_qc
    return result


@lru_cache(maxsize=_PER_SAMPLE_CACHE_SIZE)
def _load_per_sample_matrix(path: str, mtime: float, normalize: str) -> pd.DataFrame:
    """Read + normalize one cohort's per-sample matrix (the cached canonical frame).
    ``path``/``mtime`` identify the on-disk parquet (mtime keys cache invalidation);
    the matrix must already be present."""
    raw = pd.read_parquet(path)
    validate_nonnegative_expression(raw[sample_columns(raw)], context=f"{path}: linear TPM")
    # Make the matrix dense in the CANONICAL gene-id space before anything else: sum
    # alt-haplotype/patch copies into their primary gene (in LINEAR raw TPM) and relabel
    # retired ids to their successor (resolve_ensembl_id). Doing it pre-clean_tpm means
    # the alt copy inherits the primary's compartment and column totals are unchanged, so
    # clean_tpm's renormalization of every other gene is untouched; every downstream
    # accessor (cohort_stats, coverage, pooled, percentile/within-sample recompute) then
    # shares one canonical key space. (pirlygenes#465 / oncoref#135 item 6.)
    raw = _canonicalize_gene_rows(raw, sample_cols=sample_columns(raw))
    base = id_columns(raw)
    samples = sample_columns(raw)
    if normalize == "tpm_raw":
        return raw
    cohort_stem = Path(path).stem.removesuffix("_per_sample_tpm")
    empty_mode = (
        "record" if cohort_stem in source_matrices.PARTIAL_TECHNICAL_COMPARTMENT_COHORTS else "warn"
    )
    clean = clean_tpm(raw[samples], gene_table=raw[base], on_empty_compartment=empty_mode)
    out = pd.concat([raw[base].reset_index(drop=True), clean.reset_index(drop=True)], axis=1)
    if "unfilled_compartments" in clean.attrs:
        out.attrs["unfilled_compartments"] = dict(clean.attrs["unfilled_compartments"])
    if normalize == "tpm_clean_log1p":
        out[samples] = np.log1p(out[samples].to_numpy(dtype=float))
    elif normalize == "tpm_clean_hk":
        out = _housekeeping_normalize(out, samples)
    return out


def _housekeeping_normalize(df: pd.DataFrame, sample_cols) -> pd.DataFrame:
    """Divide each clean-TPM sample column by its housekeeping median-of-ratios factor.

    The factor is per-column, so this commutes with the proteoform sum and can be
    applied before or after collapse.
    """

    panel_ids = _clean_tpm_housekeeping_panel_ids()
    out, _ = tpm_to_housekeeping_normalized(
        df,
        value_cols=list(sample_cols),
        panel_ids=panel_ids,
        panel_name="clean_tpm_biological_housekeeping",
        min_panel_detected=_min_housekeeping_detected(panel_ids),
        drop_zero_panel_values=True,
        warn_on_unreliable=True,
    )
    return out


_register_derived_cache(_load_per_sample_matrix.cache_clear)


def cohort_mean_expression(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    statistic: str = "mean",
    auto_fetch: bool = True,
    proteoform: bool = False,
    scope: str = "cta",
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Per-gene **across-patient summary** of a cohort's expression (one value per
    gene, collapsed over all patients).

    The continuous cohort-level expression that downstream mechanism/correlation
    analyses need (e.g. cohort-mean TGF-β-signature expression vs aPD1 ORR) — which
    the percentile vectors and n=5 medoids don't give directly. Reduces
    :func:`per_sample_expression` over the sample axis with ``statistic`` (``"mean"``
    / ``"median"``). ``normalize`` is passed through (clean TPM by default; use
    ``"tpm_clean_log1p"`` to average in log space). A **gene-level** frame.

    With ``proteoform=True``, identical-protein paralogs are summed per sample
    (:func:`oncoref.proteoforms.collapse_to_proteoforms`) **before** the
    across-patient reduction, so the summary is over the reduced proteoform key space
    (rows carry ``proteoform_key`` — a **proteoform-level** frame, see
    :func:`oncoref.proteoforms.expression_level`). ``sample_qc`` defaults to
    ``"pass"`` so live summaries exclude source/sample QC failures; pass ``"all"`` for
    forensic parity with the current source matrix. ``scope`` selects the gene
    universe to collapse: ``"cta"`` (focused) or ``"genome"`` (every protein-coding
    gene). Returns ``Ensembl_Gene_ID``, ``Symbol`` (plus the proteoform identity
    columns when collapsed) and one ``expression`` column."""
    if statistic not in ("mean", "median"):
        raise ValueError("statistic must be 'mean' or 'median'")
    df = per_sample_expression(
        cancer_type,
        normalize=normalize,
        auto_fetch=auto_fetch,
        proteoform=proteoform,
        scope=scope,
        sample_qc=sample_qc,
    )
    id_cols = id_columns(df)
    samples = sample_columns(df)
    reducer = df[samples].mean(axis=1) if statistic == "mean" else df[samples].median(axis=1)
    out = df[id_cols].copy()
    out["expression"] = reducer.to_numpy()
    return out


def cancer_reference_expression(
    cancer_types: str | Iterable[str] | None = None,
    genes: str | Iterable[str] | None = None,
    normalize: str | Iterable[str] = "tpm_clean",
    *,
    format: str = "long",
    include_provenance: bool = True,
    include_request_metadata: bool = False,
    on_missing: str = "omit",
    auto_fetch: bool = False,
    sample_qc: str = "pass",
    reference_source: str = "artifact",
    gene_id_style: str = "oncoref",
    gene_universe: str = "artifact",
    include_gene_universe_flags: bool = False,
    source_kind: str | Iterable[str] | None = None,
    source_cohort: str | Iterable[str] | None = None,
    exclude_microarray_proxy: bool = False,
    pool: bool = False,
    collapse_cdna_identical: bool = False,
    collapse_protein_identical: bool = False,
) -> pd.DataFrame:
    """Observed tumor expression references as cohort-level clean TPM summaries.

    The default long form mirrors the downstream reference contract:
    ``Ensembl_Gene_ID``, ``Symbol``, ``cancer_code``, ``source_cohort``,
    ``normalization``, ``expression`` (median clean TPM), ``q1`` and ``q3``.
    Wide form returns one row per gene with ``<CODE>_<normalization>`` columns.

    ``normalize`` accepts one mode or an iterable of modes:
      - ``"tpm_clean"`` / ``"clean_tpm"``: shipped biological clean-TPM percentiles;
      - ``"tpm_clean_biological"``: explicit alias for that biological-only artifact;
      - ``"tpm_clean_log1p"``: stored log1p biological clean-TPM percentiles;
      - ``"tpm_raw"`` / ``"tpm"`` / ``"raw_tpm"``: raw-TPM cohort stats recomputed from
        the source matrix.

    ``reference_source="artifact"`` preserves the historical behavior: shipped
    percentile artifacts supply clean/log clean TPM and raw TPM is recomputed from
    source matrices. ``reference_source="summary_rows"`` reads the shipped
    per-source ``cancer-reference-expression`` sidecars for ``sample_qc="all"`` and
    uses the one source explicitly selected in the availability manifest. Gene and
    sample counts never promote an alternative physical source. For
    ``sample_qc="pass"`` or ``"pass_or_warn"``, it recomputes the same
    reference summaries from source matrices at read time so QC-filtered views are not
    silently backed by all-sample build artifacts.
    ``reference_source="summary_rows_all"`` is the source-union view over those
    sidecars: it preserves one row per ``(gene, cancer_code, source_cohort)`` and
    supports ``source_kind``, ``source_cohort``, ``exclude_microarray_proxy``, and
    long-form ``pool=True``. Source/code/gene filters are applied while reading the
    selected per-source shards, before any concatenation. Because the shipped
    sidecars are all-sample summaries, this mode requires ``sample_qc="all"``. Use
    ``summary_rows`` for the current QC-filtered selected-source recompute path. For
    both source filters, ``None`` means unfiltered; an empty iterable or blank scalar
    explicitly matches no sources.

    ``sample_qc="artifact"`` uses the QC policy recorded when each requested
    artifact was built. Clean/log-clean modes read the shipped percentile shard;
    raw TPM is recomputed from the bounded per-code source matrix with that code's
    recorded effective policy. Long-form provenance keeps ``sample_qc="artifact"``
    and reports the applied policy in ``sample_qc_effective``.
    ``n_reference_samples`` is the selected sample count used to build the returned
    artifact, not the wider raw source-matrix count. This mode requires
    ``reference_source="artifact"``.

    Raw-TPM mode needs the per-sample source matrix available; pass
    ``auto_fetch=True`` to download it. Raw-TPM summaries default to
    ``sample_qc="pass"`` so sparse/source-QC-failed samples do not shape new
    derived summaries; use ``"pass_or_warn"`` or ``"all"`` for audit/parity views.
    ``gene_id_style="oncoref"`` returns canonical oncoref ENSG IDs. Opt into
    ``"pirlygenes"`` only for migration wrappers that need known legacy ENSG IDs
    for rows with one-to-one remaps recorded in
    ``expression-artifact-gene-universe-deltas.csv``.
    ``gene_universe="artifact"`` preserves exact shipped rows.
    ``gene_universe="tumor_signal"`` drops filterable extras audited for the requested
    cohort/product while retaining biological oncoref-only rows.
    ``gene_universe="pirlygenes"`` is stricter: it drops audited oncoref-only
    extras without a pirlygenes counterpart while keeping documented remap targets;
    ``include_gene_universe_flags=True`` appends row-level audit columns in long output.
    Neither option synthesizes missing biological rows.
    Long output always carries ``Proteoform_ID`` and ``Member_Ensembl_Gene_IDs`` so
    gene-level and folded views share one schema. By default this is a bridge over
    the cDNA/read-recovery identity space without folding rows. Set
    ``collapse_cdna_identical=True`` or ``collapse_protein_identical=True`` to sum
    identical-locus rows in linear expression space before output. cDNA collapse is
    the read-recovery view (byte-identical CDS plus curated overrides); protein
    collapse is the genome-wide identical-protein view. Set at most one. Wide output
    keeps the historical ``Ensembl_Gene_ID``/``Symbol`` plus value columns shape.
    Missing requested cohorts are omitted by default to preserve the historical
    behavior; pass ``on_missing="empty"`` to preserve a schema-stable empty result
    with missing-request metadata in ``df.attrs["missing_requests"]`` or
    ``on_missing="raise"`` to fail when any requested cohort/mode is unavailable.
    """
    modes = _reference_normalize_modes(normalize)
    sample_qc = _validate_artifact_sample_qc(sample_qc)
    reference_source = _validate_reference_source(reference_source)
    if sample_qc == "artifact" and reference_source != "artifact":
        raise ValueError("sample_qc='artifact' requires reference_source='artifact'")
    _validate_gene_id_style(gene_id_style)
    gene_universe = _validate_artifact_gene_universe(gene_universe)
    source_kinds = _normalize_source_filter_values(source_kind)
    source_cohorts = _normalize_source_cohort_filter_values(source_cohort)
    has_source_filter = source_kinds is not None or source_cohorts is not None
    if reference_source == "summary_rows_all" and sample_qc != "all":
        raise ValueError('reference_source="summary_rows_all" requires sample_qc="all"')
    if (has_source_filter or exclude_microarray_proxy or pool) and (
        reference_source != "summary_rows_all"
    ):
        raise ValueError(
            "source_kind, source_cohort, exclude_microarray_proxy, and pool are "
            "only supported with reference_source='summary_rows_all'"
        )
    if format not in ("long", "wide"):
        raise ValueError("format must be 'long' or 'wide'")
    if reference_source == "summary_rows_all" and format != "long":
        raise ValueError('reference_source="summary_rows_all" requires format="long"')
    if pool and format != "long":
        raise ValueError("pool=True requires format='long'")
    if collapse_cdna_identical and collapse_protein_identical:
        raise ValueError("set at most one of collapse_cdna_identical or collapse_protein_identical")
    if on_missing not in ("omit", "empty", "raise"):
        raise ValueError("on_missing must be 'omit', 'empty', or 'raise'")
    if include_request_metadata and format != "long":
        raise ValueError("include_request_metadata=True requires format='long'")
    if include_gene_universe_flags and format != "long":
        raise ValueError("include_gene_universe_flags=True requires format='long'")

    requests = _reference_expression_requests(
        cancer_types,
        modes,
        reference_source=reference_source,
        sample_qc=sample_qc,
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )
    if reference_source != "summary_rows_all":
        # A named cohort whose selected source is invalid is rejected, never silently
        # dropped as "unavailable". The source-union view applies its own source filter.
        for request in requests:
            _reject_invalid_source_code(
                request["cancer_code"], context="cancer reference expression"
            )
    availability = _reference_expression_availability_for_requests(
        requests,
        modes,
        sample_qc=sample_qc,
        reference_source=reference_source,
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )
    missing = availability.loc[~availability["available"]].reset_index(drop=True)
    if on_missing == "raise" and not missing.empty:
        detail = ", ".join(
            f"{r.cancer_code}/{r.normalization}: {r.missing_reason}"
            for r in missing.itertuples(index=False)
        )
        raise ValueError(f"missing cancer reference expression artifact(s): {detail}")
    available_keys = {
        (str(r.requested_code), str(r.cancer_code), str(r.normalization))
        for r in availability.loc[availability["available"]].itertuples(index=False)
    }
    request_lookup = {
        (str(r.requested_code), str(r.cancer_code), str(r.normalization)): r
        for r in availability.loc[availability["available"]].itertuples(index=False)
    }

    long_parts: list[pd.DataFrame] = []
    wide_parts: list[pd.DataFrame] = []
    reads_summary_shards = reference_source == "summary_rows_all" or (
        reference_source == "summary_rows" and sample_qc == "all"
    )
    for request in requests:
        code = request["cancer_code"]
        summary_shard_cache: dict[Path, pd.DataFrame] | None = {} if reads_summary_shards else None
        for mode in modes:
            request_key = (request["requested_code"], code, mode)
            if request_key not in available_keys:
                continue
            request_row = request_lookup[request_key]
            effective_sample_qc = request_row.sample_qc_effective
            load_sample_qc = sample_qc
            if sample_qc == "artifact" and mode == "tpm_raw":
                if pd.isna(effective_sample_qc):
                    raise AssertionError(f"available raw artifact request lacks QC policy: {code}")
                load_sample_qc = _validate_sample_qc(str(effective_sample_qc))
            ref, method = _reference_expression_frame(
                code,
                mode,
                auto_fetch=auto_fetch,
                sample_qc=load_sample_qc,
                reference_source=reference_source,
                source_kinds=source_kinds,
                source_cohorts=source_cohorts,
                exclude_microarray_proxy=exclude_microarray_proxy,
                genes=genes,
                summary_shard_cache=summary_shard_cache,
            )
            if method not in {_REFERENCE_SUMMARY_METHOD, _REFERENCE_SUMMARY_ALL_METHOD}:
                ref = ref[_gene_filter_mask(ref, genes)].reset_index(drop=True)
            ref = _apply_artifact_gene_universe(
                ref,
                product="cancer_reference_expression",
                cancer_codes=[code],
                gene_universe=gene_universe,
                include_gene_universe_flags=include_gene_universe_flags,
            )
            collapse_kind = (
                "cdna"
                if collapse_cdna_identical
                else "protein"
                if collapse_protein_identical
                else None
            )
            if collapse_kind is not None:
                ref = _collapse_reference_identical_loci(
                    ref,
                    kind=collapse_kind,
                    group_keys=_reference_collapse_group_keys(ref),
                    identity_style=gene_id_style,
                )
            ref = _apply_gene_id_style(
                ref,
                product="cohort_gene_percentiles",
                cancer_codes=[code],
                gene_id_style=gene_id_style,
                alias_expand_remaps=gene_universe == "pirlygenes",
            )
            if collapse_kind is None:
                ref = _annotate_reference_proteoform_bridge(
                    ref, kind="cdna", identity_style=gene_id_style
                )
            label = _REFERENCE_NORMALIZE_LABELS[mode]
            if format == "long":
                value_cols = [
                    "Ensembl_Gene_ID",
                    "Symbol",
                    "p25",
                    "p50",
                    "p75",
                    "Proteoform_ID",
                    "Member_Ensembl_Gene_IDs",
                ]
                if include_gene_universe_flags:
                    value_cols.extend(_ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS)
                part = ref[value_cols].copy()
                part.insert(2, "cancer_code", code)
                part["normalization"] = label
                if include_request_metadata:
                    part["requested_code"] = request_row.requested_code
                    part["request_kind"] = request_row.request_kind
                    part["available"] = bool(request_row.available)
                    part["missing_reason"] = str(request_row.missing_reason)
                part = part.rename(columns={"p50": "expression", "p25": "q1", "p75": "q3"})
                source_union_identity = method == _REFERENCE_SUMMARY_ALL_METHOD
                if source_union_identity:
                    for col in _REFERENCE_SOURCE_UNION_IDENTITY_COLUMNS:
                        if col in ref.columns:
                            part[col] = ref[col].to_numpy()
                        else:
                            part[col] = pd.NA
                if include_provenance:
                    if method == _REFERENCE_SUMMARY_ALL_METHOD:
                        for col in _REFERENCE_PROVENANCE_COLUMNS:
                            if col in ref.columns:
                                part[col] = ref[col].to_numpy()
                    else:
                        # Availability has already applied the artifact build count
                        # and QC policy. Reuse it so returned provenance cannot drift
                        # back to the wider raw source-matrix metadata.
                        for col in _REFERENCE_PROVENANCE_COLUMNS:
                            part[col] = getattr(request_row, col)
                    for col in _REFERENCE_PROVENANCE_COLUMNS:
                        if col not in part.columns:
                            part[col] = pd.NA
                long_parts.append(
                    part[
                        _reference_long_columns(
                            include_provenance,
                            include_request_metadata,
                            include_gene_universe_flags,
                            source_union_identity=source_union_identity,
                            include_proteoform_columns=True,
                        )
                    ]
                )
            else:
                suffix = _REFERENCE_WIDE_SUFFIXES[mode]
                part = ref[["Ensembl_Gene_ID", "Symbol", "p50"]].rename(
                    columns={"p50": f"{code}_{suffix}"}
                )
                wide_parts.append(part)

    if format == "long":
        cols = _reference_long_columns(
            include_provenance,
            include_request_metadata,
            include_gene_universe_flags,
            source_union_identity=reference_source == "summary_rows_all",
            include_proteoform_columns=True,
        )
        if not long_parts:
            out = pd.DataFrame(columns=cols)
        else:
            out = pd.concat(long_parts, ignore_index=True)
        if pool and not out.empty:
            out = _pool_reference_expression_rows(out)
        _attach_reference_expression_attrs(
            out,
            availability if on_missing != "omit" else None,
            gene_id_style=gene_id_style,
            gene_universe=gene_universe,
            include_gene_universe_flags=include_gene_universe_flags,
            reference_source=reference_source,
        )
        if any(mode in {"tpm_clean", "tpm_clean_biological", "tpm_clean_log1p"} for mode in modes):
            _attach_gene_universe_delta_attrs(
                out,
                product="cancer_reference_expression",
                cancer_codes=[r["cancer_code"] for r in requests],
            )
        return out

    if not wide_parts:
        out = pd.DataFrame(columns=["Ensembl_Gene_ID", "Symbol"])
    else:
        out = _merge_cancer_reference_wide_parts(wide_parts)
    _attach_reference_expression_attrs(
        out,
        availability if on_missing != "omit" else None,
        gene_id_style=gene_id_style,
        gene_universe=gene_universe,
        include_gene_universe_flags=include_gene_universe_flags,
        reference_source=reference_source,
    )
    if any(mode in {"tpm_clean", "tpm_clean_biological", "tpm_clean_log1p"} for mode in modes):
        _attach_gene_universe_delta_attrs(
            out,
            product="cancer_reference_expression",
            cancer_codes=[r["cancer_code"] for r in requests],
        )
    return out


def _reference_collapse_group_keys(ref: pd.DataFrame) -> list[str]:
    """Columns that keep identical-locus summation inside one source context.

    ``n_detected`` varies by gene and is aggregated with ``max`` during a fold;
    treating it as an identity key would split valid groups and then remove the
    key before the aggregate merge.
    """
    return [
        c
        for c in _REFERENCE_SOURCE_UNION_IDENTITY_COLUMNS
        if c in ref.columns and c != "n_detected"
    ]


def _clean_identity_label(value) -> str | None:
    if pd.isna(value):
        return None
    label = str(value).strip()
    if not label or label.lower() == "nan":
        return None
    return label


_PIRLYGENES_PROTEIN_GROUP_EXCLUSIONS = frozenset(
    {
        frozenset({"ENSG00000169789", "ENSG00000169807"}),  # PRY/PRY2
    }
)

_PIRLYGENES_PROTEIN_LABEL_OVERRIDES = {
    frozenset({"ENSG00000183889", "ENSG00000233024"}): "NPIPA6/9",
    frozenset({"ENSG00000163611", "ENSG00000285943"}): "SPICE1/SPICE1-CFAP44",
}


def _pirlygenes_identity_label(symbols: Iterable[object], *, fallback: str) -> str:
    """Recreate pirlygenes' compact identical-locus group identifier."""
    parts: list[str] = []
    for value in symbols:
        label = _clean_identity_label(value)
        if label is None or re.fullmatch(r"ENSG\d+", label) or label in parts:
            continue
        parts.append(label)
    parts.sort(
        key=lambda value: [
            int(token) if token.isdigit() else token for token in re.split(r"(\d+)", value)
        ]
    )
    if not parts:
        return fallback
    if len(parts) == 1:
        return parts[0]
    first, last = min(parts), max(parts)
    prefix_len = 0
    for ch_a, ch_b in zip(first, last):
        if ch_a != ch_b:
            break
        prefix_len += 1
    prefix = first[:prefix_len]
    if (
        prefix
        and prefix[-1].isdigit()
        and any(len(value) > len(prefix) and value[len(prefix)].isdigit() for value in parts)
    ):
        prefix = prefix.rstrip("0123456789")
    suffixes = [value[len(prefix) :] for value in parts]
    if prefix and all(suffixes):
        return prefix + "/".join(suffixes)
    return "/".join(parts)


@lru_cache(maxsize=4)
def _identical_locus_identity_maps(
    kind: str, identity_style: str = "oncoref"
) -> tuple[dict[str, str], dict[str, str]]:
    """Return ``({member ENSG -> identity}, {identity -> member ENSGs})``.

    ``kind="cdna"`` is the read-recovery collapse: byte-identical CDS groups plus
    the curated override table. ``kind="protein"`` is the genome-wide
    identical-protein collapse. The public registry/accessor remains gene-level by
    default; these maps are used only by explicit compatibility collapse requests.
    """
    _validate_gene_id_style(identity_style)
    if kind == "cdna":
        groups = cdna_identical_groups()
        member_to_canonical = {
            unversioned(member): unversioned(canon)
            for member, canon in zip(
                groups["ensembl_gene_id"],
                groups["group_canonical_ensembl_gene_id"],
            )
        }
        if identity_style == "pirlygenes":
            canonical_to_identity = {
                unversioned(canonical): _pirlygenes_identity_label(
                    group["symbol"], fallback=unversioned(canonical)
                )
                for canonical, group in groups.groupby(
                    "group_canonical_ensembl_gene_id", sort=False
                )
            }
        else:
            canonical_to_identity = {
                unversioned(canon): _clean_identity_label(symbol) or unversioned(canon)
                for canon, symbol in zip(
                    groups["group_canonical_ensembl_gene_id"],
                    groups["group_canonical_symbol"],
                )
            }
        overrides = proteoform_collapse_overrides()
        if not overrides.empty:
            from .proteoforms import proteoform_groups

            protein_groups = proteoform_groups(scope="genome")
            override_symbols = {
                unversioned(canon): _clean_identity_label(symbol) or unversioned(canon)
                for canon, symbol in zip(
                    overrides["group_canonical_ensembl_gene_id"],
                    overrides["group_symbol"],
                )
            }
            for canonical_id, identity in override_symbols.items():
                labels = protein_groups.loc[
                    protein_groups["member_gene_id"].astype(str).map(unversioned) == canonical_id,
                    "proteoform_id",
                ]
                if labels.empty:
                    continue
                label = str(labels.iloc[0])
                members = protein_groups.loc[
                    protein_groups["proteoform_id"].astype(str) == label,
                    "member_gene_id",
                ]
                for member in members:
                    member_to_canonical[unversioned(member)] = canonical_id
                canonical_to_identity[canonical_id] = identity
        member_to_identity = {
            member: canonical_to_identity.get(canonical, canonical)
            for member, canonical in member_to_canonical.items()
        }
    elif kind == "protein":
        from .proteoforms import proteoform_groups

        groups = proteoform_groups(scope="genome")
        if identity_style == "pirlygenes":
            member_to_identity = {}
            for _, group in groups.groupby("proteoform_id", sort=False):
                members = frozenset(group["member_gene_id"].astype(str).map(unversioned))
                if members in _PIRLYGENES_PROTEIN_GROUP_EXCLUSIONS:
                    continue
                identity = _PIRLYGENES_PROTEIN_LABEL_OVERRIDES.get(members)
                if identity is None:
                    identity = _pirlygenes_identity_label(
                        group["member_symbol"], fallback=min(members)
                    )
                member_to_identity.update(dict.fromkeys(members, identity))
        else:
            member_to_identity = {
                unversioned(member): str(identity)
                for member, identity in zip(groups["member_gene_id"], groups["proteoform_id"])
            }
    else:
        raise ValueError("kind must be 'cdna' or 'protein'")

    identity_to_members: dict[str, set[str]] = defaultdict(set)
    for member, identity in member_to_identity.items():
        identity_to_members[identity].add(member)
    return member_to_identity, {
        identity: ";".join(sorted(members)) for identity, members in identity_to_members.items()
    }


def _collapse_reference_identical_loci(
    df: pd.DataFrame,
    *,
    kind: str,
    group_keys: list[str],
    identity_style: str = "oncoref",
) -> pd.DataFrame:
    """Collapse identical loci in a reference-expression long-like frame."""
    if df.empty:
        out = df.copy()
        out["Proteoform_ID"] = []
        out["Member_Ensembl_Gene_IDs"] = []
        return out
    member_to_identity, identity_members = _identical_locus_identity_maps(kind, identity_style)
    work = df.reset_index(drop=True).copy()
    work["_ord"] = range(len(work))
    gene_ids = work["Ensembl_Gene_ID"].astype(str).map(unversioned)
    identity = gene_ids.map(member_to_identity).fillna(gene_ids)
    work["_identity"] = identity
    work["_is_singleton"] = ~gene_ids.isin(member_to_identity)
    full_keys = ["_identity", *group_keys]

    sum_cols = [c for c in ("p25", "p50", "p75") if c in work.columns]
    max_cols = [c for c in ("n_detected",) if c in work.columns]
    rep = (
        work.sort_values([*full_keys, "_is_singleton", "_ord"], ascending=True)
        .drop_duplicates(full_keys, keep="first")
        .drop(columns=sum_cols + max_cols, errors="ignore")
    )
    observed_groups = work.groupby(full_keys, as_index=False, observed=True)
    out = rep
    if sum_cols:
        out = out.merge(
            observed_groups[sum_cols].sum(min_count=1),
            on=full_keys,
            how="left",
        )
    if max_cols:
        out = out.merge(
            observed_groups[max_cols].max(),
            on=full_keys,
            how="left",
        )
    folded = out["_identity"].isin(identity_members)
    out.loc[folded, "Ensembl_Gene_ID"] = out.loc[folded, "_identity"]
    out.loc[folded, "Symbol"] = out.loc[folded, "_identity"]
    out["Proteoform_ID"] = out["_identity"]
    out["Member_Ensembl_Gene_IDs"] = [
        identity_members.get(identity_id, identity_id) for identity_id in out["_identity"]
    ]
    keep = list(df.columns)
    for col in ("Proteoform_ID", "Member_Ensembl_Gene_IDs"):
        if col not in keep:
            keep.append(col)
    return out.sort_values("_ord").reset_index(drop=True)[keep]


def _annotate_reference_proteoform_bridge(
    df: pd.DataFrame, *, kind: str, identity_style: str = "oncoref"
) -> pd.DataFrame:
    """Add pirlygenes-compatible gene/proteoform bridge columns without folding rows."""
    out = df.copy()
    if "Proteoform_ID" not in out.columns:
        if out.empty:
            out["Proteoform_ID"] = pd.Series(dtype="object")
        else:
            member_to_identity, _ = _identical_locus_identity_maps(kind, identity_style)
            gene_ids = out["Ensembl_Gene_ID"].astype(str).map(unversioned)
            out["Proteoform_ID"] = gene_ids.map(member_to_identity).fillna(gene_ids).to_numpy()
    if "Member_Ensembl_Gene_IDs" not in out.columns:
        if out.empty:
            out["Member_Ensembl_Gene_IDs"] = pd.Series(dtype="object")
        else:
            out["Member_Ensembl_Gene_IDs"] = out["Ensembl_Gene_ID"].astype(str).to_numpy()
    return out


def _reference_long_columns(
    include_provenance: bool,
    include_request_metadata: bool,
    include_gene_universe_flags: bool,
    *,
    source_union_identity: bool = False,
    include_proteoform_columns: bool = False,
) -> list[str]:
    cols = ["Ensembl_Gene_ID", "Symbol"]
    if include_proteoform_columns:
        cols += ["Proteoform_ID", "Member_Ensembl_Gene_IDs"]
    cols += ["cancer_code", "normalization"]
    if include_request_metadata:
        cols += _REFERENCE_REQUEST_COLUMNS
    if include_provenance:
        cols += _REFERENCE_PROVENANCE_COLUMNS
    elif source_union_identity:
        cols += _REFERENCE_SOURCE_UNION_IDENTITY_COLUMNS
    if include_gene_universe_flags:
        cols += _ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS
    return [*cols, "expression", "q1", "q3"]


def _reference_expression_frame(
    code: str,
    mode: str,
    *,
    auto_fetch: bool,
    sample_qc: str,
    reference_source: str,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
    genes: str | Iterable[str] | None = None,
    summary_shard_cache: dict[Path, pd.DataFrame] | None = None,
) -> tuple[pd.DataFrame, str]:
    if reference_source == "summary_rows_all":
        return (
            _reference_summary_expression_frame(
                code,
                mode,
                all_sources=True,
                source_kinds=source_kinds,
                source_cohorts=source_cohorts,
                exclude_microarray_proxy=exclude_microarray_proxy,
                genes=genes,
                shard_cache=summary_shard_cache,
            ),
            _REFERENCE_SUMMARY_ALL_METHOD,
        )
    if reference_source == "summary_rows":
        if sample_qc == "all":
            return (
                _reference_summary_expression_frame(
                    code,
                    mode,
                    genes=genes,
                    shard_cache=summary_shard_cache,
                ),
                _REFERENCE_SUMMARY_METHOD,
            )
        normalize = "tpm_clean" if mode == "tpm_clean_biological" else mode
        return (
            cohort_stats(code, normalize=normalize, auto_fetch=auto_fetch, sample_qc=sample_qc),
            _REFERENCE_RECOMPUTED_METHOD,
        )
    if mode in {"tpm_clean", "tpm_clean_biological"}:
        return (
            cohort_gene_percentiles(code, as_tpm=True, auto_fetch=auto_fetch, sample_qc=sample_qc),
            "percentile_shard",
        )
    if mode == "tpm_clean_log1p":
        return (
            cohort_gene_percentiles(code, as_tpm=False, auto_fetch=auto_fetch, sample_qc=sample_qc),
            "percentile_shard_log1p",
        )
    if mode == "tpm_raw":
        return (
            cohort_stats(code, normalize="tpm_raw", auto_fetch=auto_fetch, sample_qc=sample_qc),
            "source_matrix_stats",
        )
    raise AssertionError(f"unhandled reference normalize mode: {mode}")


def _canonical_reference_summary_source_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Expose exact Treehouse TCGA source identities in old summary bundles.

    DDLPS, WDLPS, and PLEOLPS are GDC-histology-derived overlays. All other
    legacy rows belong to the generic cohort selected by TCGA provenance.
    """

    if df.empty:
        return df
    codes = df["cancer_code"]
    sources = df["source_cohort"]
    sarc_histology = codes.isin(TREEHOUSE_TCGA_SARC_HISTOLOGY_CODES) & sources.isin(
        [TREEHOUSE_TCGA_LEGACY_COHORT, TREEHOUSE_TCGA_SAMPLES_COHORT]
    )
    generic_tcga = sources.eq(TREEHOUSE_TCGA_LEGACY_COHORT) & ~sarc_histology
    if not (sarc_histology.any() or generic_tcga.any()):
        return df

    out = df.copy(deep=False)
    canonical_sources = sources
    if isinstance(sources.dtype, pd.CategoricalDtype):
        missing_categories = [
            replacement
            for replacement, mask in (
                (TREEHOUSE_TCGA_SARC_HISTOLOGY_COHORT, sarc_histology),
                (TREEHOUSE_TCGA_SAMPLES_COHORT, generic_tcga),
            )
            if mask.any() and replacement not in sources.cat.categories
        ]
        if missing_categories:
            canonical_sources = sources.cat.add_categories(missing_categories)
        canonical_sources = canonical_sources.copy()
    else:
        canonical_sources = sources.copy()
    canonical_sources.loc[sarc_histology] = TREEHOUSE_TCGA_SARC_HISTOLOGY_COHORT
    canonical_sources.loc[generic_tcga] = TREEHOUSE_TCGA_SAMPLES_COHORT
    if isinstance(canonical_sources.dtype, pd.CategoricalDtype):
        canonical_sources = canonical_sources.cat.remove_unused_categories()
    out["source_cohort"] = canonical_sources
    return out


def _reference_summary_available_codes(
    *,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
) -> list[str]:
    table = _filter_reference_summary_sources(
        _reference_summary_availability_table(),
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )
    if table.empty:
        return []
    return sorted(table["cancer_code"].astype(str).unique())


@lru_cache(maxsize=1)
def _reference_summary_source_table() -> pd.DataFrame:
    """Validated source choices from the compact manifest, without expression rows."""

    columns = [
        "cancer_code",
        "source_cohort",
        "source_project",
        "source_version",
        "source_pmid",
        "source_type",
        "source_scale_class",
        "linear_tpm_comparable",
        "tumor_origin",
        "metastasis_site",
        "n_reference_genes",
        "n_reference_samples",
        "processing_pipeline",
        "notes",
        "selected",
    ]
    table = _canonical_reference_summary_source_labels(
        _reference_summary_availability_table()
    ).copy()
    if table.empty:
        return pd.DataFrame(columns=columns)
    has_selected_column = "selected" in table.columns
    for column in columns:
        if column not in table.columns:
            table[column] = pd.NA
    if not has_selected_column:
        source_counts = table.groupby("cancer_code", observed=True)["source_cohort"].size()
        ambiguous = source_counts[source_counts.ne(1)]
        if not ambiguous.empty:
            raise ValueError(
                "reference availability without a selected column must contain exactly "
                f"one source per cancer code: {ambiguous.astype(int).to_dict()}"
            )
        selected = pd.Series(True, index=table.index, dtype=bool)
    else:
        selected = table["selected"].map(_optional_bool)
    if selected.isna().any():
        invalid = sorted(table.loc[selected.isna(), "selected"].astype(str).unique())
        raise ValueError(f"reference availability has invalid selected values: {invalid}")
    table["selected"] = selected.astype(bool)
    selected_counts = table.groupby("cancer_code", observed=True)["selected"].sum()
    invalid_codes = selected_counts[selected_counts.ne(1)]
    if not invalid_codes.empty:
        raise ValueError(
            "reference availability must select exactly one source per cancer code: "
            f"{invalid_codes.astype(int).to_dict()}"
        )
    table["_origin_rank"] = (
        table["tumor_origin"]
        .fillna("")
        .astype(str)
        .str.lower()
        .map(_REFERENCE_TUMOR_ORIGIN_RANK)
        .fillna(99)
        .astype(int)
    )
    table["_source_sort"] = table["source_cohort"].astype("string").fillna("")
    table = table.sort_values(
        [
            "cancer_code",
            "selected",
            "n_reference_genes",
            "n_reference_samples",
            "_origin_rank",
            "_source_sort",
        ],
        ascending=[True, False, False, False, True, True],
        kind="stable",
    )
    return table.drop(columns=["_origin_rank", "_source_sort"])[columns].reset_index(drop=True)


_register_derived_cache(_reference_summary_source_table.cache_clear)


def _reference_summary_selected_source(code: str) -> dict | None:
    table = _reference_summary_source_table()
    if table.empty:
        return None
    matches = table.loc[(table["cancer_code"].astype(str) == code) & table["selected"]]
    if matches.empty:
        return None
    return matches.iloc[0].to_dict()


#: Source scale classes whose expression values are known to be invalid. Everything
#: derived from such a source is withheld, not served with a warning.
INVALID_SOURCE_SCALE_CLASSES = frozenset({"invalid_legacy_tpm"})
_INVALID_SOURCE_GUIDANCE = (
    "rebuild it from verified linear expression "
    "(see docs/audits/nutm-expression-scale.md and oncoref#604)"
)


def _invalid_source_message(row) -> str:
    return (
        f"source {row['source_cohort']} has an invalid expression scale "
        f"({row['source_scale_class']}); its derived artifacts are withheld until "
        f"{_INVALID_SOURCE_GUIDANCE}"
    )


@lru_cache(maxsize=1)
def _invalid_selected_source_codes() -> dict[str, str]:
    """Cancer codes whose *selected* reference source has an invalid expression scale.

    Derived artifacts (percentiles, representatives, within-sample fractions and the
    cohort summaries) are built from the selected source, so a code listed here must
    not be served from any of them.
    """
    table = _reference_summary_source_table()
    if table.empty:
        return {}
    flagged = table["source_scale_class"].astype("string").isin(INVALID_SOURCE_SCALE_CLASSES)
    return {
        str(row["cancer_code"]): _invalid_source_message(row)
        for _, row in table.loc[table["selected"] & flagged].iterrows()
    }


_register_derived_cache(_invalid_selected_source_codes.cache_clear)


def _reject_invalid_source_code(code: str, *, context: str) -> None:
    reason = _invalid_selected_source_codes().get(str(code))
    if reason is not None:
        raise ValueError(f"{context} for {code}: {reason}")


def _reference_summary_shard_path(root: Path, code: str, source_cohort: str) -> Path:
    """Resolve the one physical shard that contains a source/code pair."""

    source_names = [source_cohort]
    if source_cohort == TREEHOUSE_TCGA_SARC_HISTOLOGY_COHORT and code in (
        TREEHOUSE_TCGA_SARC_HISTOLOGY_CODES
    ):
        source_names.extend([TREEHOUSE_TCGA_SAMPLES_COHORT, TREEHOUSE_TCGA_LEGACY_COHORT])
    elif source_cohort == TREEHOUSE_TCGA_SAMPLES_COHORT:
        source_names.append(TREEHOUSE_TCGA_LEGACY_COHORT)
    stems = tuple(stem for source in source_names for stem in (f"{source}__{code}", source))
    for stem in stems:
        for suffix in (".csv.gz", ".csv"):
            path = root / f"{stem}{suffix}"
            if path.is_file():
                return path
    expected = ", ".join(f"{stem}.csv[.gz]" for stem in stems)
    raise FileNotFoundError(
        f"no cancer reference expression shard for {code}/{source_cohort}; expected {expected}"
    )


def _read_reference_summary_shard(
    path: Path,
    *,
    shard_cache: dict[Path, pd.DataFrame] | None,
) -> pd.DataFrame:
    """Read one source shard, optionally reusing it within one public request."""

    if shard_cache is not None and path in shard_cache:
        return shard_cache[path]
    frame = _canonical_reference_summary_source_labels(pd.read_csv(path, low_memory=False))
    if shard_cache is not None:
        shard_cache[path] = frame
    return frame


def _reference_summary_rows_for_sources(
    code: str,
    sources: pd.DataFrame,
    *,
    genes: str | Iterable[str] | None = None,
    shard_cache: dict[Path, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    """Load only the shards named by ``sources`` and filter before concatenation."""

    if sources.empty:
        return pd.DataFrame()
    invalid = sources.loc[
        sources["source_scale_class"].astype("string").isin(INVALID_SOURCE_SCALE_CLASSES)
    ]
    if not invalid.empty:
        raise ValueError(f"{code}: {_invalid_source_message(invalid.iloc[0])}")
    root = _bundle_subdir(_REFERENCE_SUMMARY_DATASET, auto_fetch=True)
    source_ids = set(sources["source_cohort"].dropna().astype(str))
    paths = {_reference_summary_shard_path(root, str(code), source) for source in source_ids}
    parts = []
    for path in sorted(paths):
        frame = _read_reference_summary_shard(path, shard_cache=shard_cache)
        mask = frame["cancer_code"].astype(str).eq(str(code))
        mask &= frame["source_cohort"].astype(str).isin(source_ids)
        selected = frame.loc[mask]
        if genes is not None and not selected.empty:
            selected = selected.loc[_gene_filter_mask(selected, genes)]
        if not selected.empty:
            parts.append(selected.copy())
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def _reference_summary_selected_rows(
    code: str,
    *,
    genes: str | Iterable[str] | None = None,
    shard_cache: dict[Path, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    selected = _reference_summary_selected_source(code)
    if selected is None:
        return pd.DataFrame()
    return _reference_summary_rows_for_sources(
        code,
        pd.DataFrame([selected]),
        genes=genes,
        shard_cache=shard_cache,
    )


def _normalize_source_filter_values(values: str | Iterable[str] | None) -> set[str] | None:
    """Normalize a source filter while preserving ``None`` versus no matches."""
    if values is None:
        return None
    raw = [values] if isinstance(values, str) else list(values)
    return {str(v) for v in raw if str(v)}


def _normalize_source_cohort_filter_values(
    values: str | Iterable[str] | None,
) -> set[str] | None:
    """Normalize exact cohort IDs while preserving unknown filters as no-match values."""
    normalized = _normalize_source_filter_values(values)
    if normalized is None:
        return None
    return {resolve_cohort_id(value, strict=False) or value for value in normalized}


@lru_cache(maxsize=1)
def _source_cohort_kind_map() -> dict[str, str]:
    registry = cohort_registry_df()
    if registry.empty or not {"cohort_id", "kind"} <= set(registry.columns):
        return {}
    return dict(zip(registry["cohort_id"].astype(str), registry["kind"].astype(str)))


@lru_cache(maxsize=1)
def _source_cohort_metadata_map() -> dict[str, dict]:
    """Metadata keyed by exact first-class source cohort ID."""
    registry = cohort_registry_df()
    if registry.empty or "cohort_id" not in registry.columns:
        return {}
    return {
        str(row["cohort_id"]): {
            column: None if pd.isna(row[column]) else row[column]
            for column in ("kind", "source_project", "assay", "provenance")
            if column in registry.columns
        }
        for _, row in registry.iterrows()
    }


_register_derived_cache(_source_cohort_kind_map.cache_clear)
_register_derived_cache(_source_cohort_metadata_map.cache_clear)


def _filter_reference_summary_sources(
    table: pd.DataFrame,
    *,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
) -> pd.DataFrame:
    out = table.copy()
    if out.empty:
        return out
    if source_cohorts is not None:
        out = out.loc[out["source_cohort"].astype(str).isin(source_cohorts)]
    if source_kinds is not None:
        kind_map = _source_cohort_kind_map()
        out = out.loc[out["source_cohort"].astype(str).map(kind_map).isin(source_kinds)]
    if exclude_microarray_proxy:
        scale = out.get("source_scale_class", pd.Series(pd.NA, index=out.index))
        source_type = out.get("source_type", pd.Series(pd.NA, index=out.index))
        scale = scale.astype("string").fillna("").str.strip()
        source_type = source_type.astype("string").fillna("").str.strip()

        missing_scale = scale.str.lower().isin({"", "unknown", "nan", "none"})
        missing_type = source_type.str.lower().isin({"", "unknown", "nan", "none"})
        if missing_scale.any() or missing_type.any():
            inferred = [
                cancer_reference_expression_source_metadata(str(code), source_cohort=str(source))
                for code, source in zip(out["cancer_code"], out["source_cohort"])
            ]
            inferred_scale = pd.Series(
                [row.get("source_scale_class") for row in inferred], index=out.index
            ).astype("string")
            inferred_type = pd.Series(
                [row.get("source_type") for row in inferred], index=out.index
            ).astype("string")
            scale = scale.mask(missing_scale, inferred_scale.fillna(""))
            source_type = source_type.mask(missing_type, inferred_type.fillna(""))

        normalized_scale = scale.str.lower()
        normalized_type = source_type.str.lower()
        scale_is_known = ~normalized_scale.isin({"", "unknown", "nan", "none"})
        is_microarray = scale_is_known & normalized_scale.str.contains("microarray", regex=False)
        is_microarray |= ~scale_is_known & normalized_type.str.contains("microarray", regex=False)

        # Compatibility for older manifests whose structured scale/type fields are
        # both unavailable or generic (for example ``source_type='geo-matrix'``).
        pipeline = out.get("processing_pipeline", pd.Series("", index=out.index))
        pipeline = pipeline.astype("string").fillna("")
        text = pipeline.astype(str).str.lower()
        legacy_proxy = text.str.contains("microarray|tpm_proxy|tpm-proxy", regex=True)
        is_microarray |= ~scale_is_known & legacy_proxy
        out = out.loc[~is_microarray]
    return out


def _reference_summary_all_rows(
    code: str,
    *,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
    genes: str | Iterable[str] | None = None,
    shard_cache: dict[Path, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    sources = _filter_reference_summary_sources(
        _reference_summary_source_table(),
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )
    sources = sources.loc[sources["cancer_code"].astype(str) == str(code)]
    flagged = sources["source_scale_class"].astype("string").isin(INVALID_SOURCE_SCALE_CLASSES)
    if flagged.any():
        asked = {str(c) for c in (source_cohorts or ())}
        if not (flagged & sources["source_cohort"].astype(str).isin(asked)).any():
            omitted = sorted(sources.loc[flagged, "source_cohort"].astype(str))
            warnings.warn(
                f"Omitting {code} source(s) {omitted}: invalid expression scale; "
                f"rebuild from verified linear expression (oncoref#604)",
                UserWarning,
                stacklevel=2,
            )
            sources = sources.loc[~flagged]
    return _reference_summary_rows_for_sources(
        code,
        sources,
        genes=genes,
        shard_cache=shard_cache,
    )


def _reference_summary_expression_frame(
    code: str,
    mode: str,
    *,
    all_sources: bool = False,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
    genes: str | Iterable[str] | None = None,
    shard_cache: dict[Path, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    rows = (
        _reference_summary_all_rows(
            code,
            source_kinds=source_kinds,
            source_cohorts=source_cohorts,
            exclude_microarray_proxy=exclude_microarray_proxy,
            genes=genes,
            shard_cache=shard_cache,
        )
        if all_sources
        else _reference_summary_selected_rows(
            code,
            genes=genes,
            shard_cache=shard_cache,
        )
    )
    if rows.empty:
        return pd.DataFrame(columns=["Ensembl_Gene_ID", "Symbol", "p25", "p50", "p75"])
    if mode == "tpm_raw":
        col_map = {"TPM_q1": "p25", "TPM_median": "p50", "TPM_q3": "p75"}
    elif mode in {"tpm_clean", "tpm_clean_biological", "tpm_clean_log1p"}:
        col_map = {"TPM_clean_q1": "p25", "TPM_clean_median": "p50", "TPM_clean_q3": "p75"}
    else:
        raise AssertionError(f"unhandled reference normalize mode: {mode}")
    provenance_cols = [
        "cancer_code",
        "source_cohort",
        "source_project",
        "source_version",
        "n_samples",
        "n_detected",
        "processing_pipeline",
        "notes",
        "tumor_origin",
        "metastasis_site",
    ]
    keep = ["Ensembl_Gene_ID", "Symbol", *col_map]
    if all_sources:
        keep.extend(c for c in provenance_cols if c in rows.columns)
    out = rows[keep].rename(columns=col_map).copy()
    if mode == "tpm_clean_log1p":
        for col in ("p25", "p50", "p75"):
            out[col] = np.log1p(pd.to_numeric(out[col], errors="coerce"))
    if all_sources:
        out = _attach_summary_row_provenance(out, code=code)
    return out


def _attach_summary_row_provenance(df: pd.DataFrame, *, code: str) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    observed_counts = (
        out.groupby("source_cohort", dropna=False, observed=True)
        .agg(
            n_reference_genes=("Ensembl_Gene_ID", "nunique"),
            n_reference_samples=("n_samples", "max"),
        )
        .reset_index()
        .set_index("source_cohort")
        .to_dict("index")
    )
    source_table = _reference_summary_source_table()
    source_table = source_table.loc[source_table["cancer_code"].astype(str).eq(str(code))]
    manifest_counts = source_table.set_index("source_cohort")[
        ["n_reference_genes", "n_reference_samples"]
    ].to_dict("index")
    source_records = source_table.set_index("source_cohort").to_dict("index")
    meta_by_source = {
        str(source): _reference_summary_metadata_from_row(
            code,
            {
                "source_cohort": str(source),
                **source_records.get(str(source), {}),
            },
        )
        for source in out["source_cohort"].astype("string").fillna("").unique()
    }
    source_key = out["source_cohort"].astype("string").fillna("")
    out["source_type"] = source_key.map(
        lambda source: meta_by_source.get(source, {}).get("source_type")
    )
    out["source_pmid"] = source_key.map(
        lambda source: meta_by_source.get(source, {}).get("source_pmid")
    )
    out["source_unit"] = source_key.map(lambda source: meta_by_source.get(source, {}).get("unit"))
    out["source_scale_class"] = source_key.map(
        lambda source: meta_by_source.get(source, {}).get("source_scale_class")
    )
    out["linear_tpm_comparable"] = source_key.map(
        lambda source: bool(meta_by_source.get(source, {}).get("linear_tpm_comparable"))
    )
    out["n_reference_genes"] = source_key.map(
        lambda source: manifest_counts.get(source, observed_counts.get(source, {})).get(
            "n_reference_genes"
        )
    )
    out["n_reference_samples"] = source_key.map(
        lambda source: manifest_counts.get(source, observed_counts.get(source, {})).get(
            "n_reference_samples"
        )
    )
    out["reference_method"] = _REFERENCE_SUMMARY_ALL_METHOD
    out["sample_qc"] = "all"
    out["sample_qc_effective"] = "all"
    out["data_version"] = DATA_VERSION
    out["source_matrix_version"] = SOURCE_MATRIX_VERSION
    return out


def _first_nonempty(values: pd.Series) -> object:
    for value in values:
        if pd.notna(value) and str(value):
            return value
    return pd.NA


def _join_nonempty(values: pd.Series) -> str:
    vals = sorted({str(v) for v in values if pd.notna(v) and str(v)})
    return ";".join(vals)


def _pool_reference_expression_rows(long: pd.DataFrame) -> pd.DataFrame:
    """Pool source-union reference rows by gene with n-sample weighted means."""
    if long.empty:
        return long
    out_rows: list[dict[str, object]] = []
    group_cols = [
        col
        for col in (
            "Ensembl_Gene_ID",
            "cancer_code",
            "normalization",
        )
        if col in long.columns
    ]
    for keys, group in long.groupby(group_cols, dropna=False, sort=False, observed=True):
        row = dict(zip(group_cols, keys))
        physical = group
        if "source_cohort" in group.columns and group["source_cohort"].duplicated().any():
            for source, duplicate_rows in group.groupby(
                "source_cohort", dropna=False, sort=False, observed=True
            ):
                expression_values = pd.to_numeric(
                    duplicate_rows["expression"], errors="coerce"
                ).dropna()
                sample_values = pd.to_numeric(
                    duplicate_rows.get(
                        "n_reference_samples",
                        duplicate_rows.get("n_samples", pd.Series(1, index=duplicate_rows.index)),
                    ),
                    errors="coerce",
                ).dropna()
                if expression_values.nunique() > 1 or sample_values.nunique() > 1:
                    raise ValueError(
                        "cannot pool duplicate rows for one physical source and output gene: "
                        f"{source}/{row.get('Ensembl_Gene_ID')}"
                    )
            physical = group.drop_duplicates("source_cohort", keep="first")

        gene_id = str(row.get("Ensembl_Gene_ID", ""))
        row["Symbol"] = canonical_gene_symbol(gene_id) or _first_nonempty(group["Symbol"])
        proteoform_ids = {
            str(value)
            for value in group.get("Proteoform_ID", pd.Series(dtype=object))
            if pd.notna(value) and str(value)
        }
        row["Proteoform_ID"] = next(iter(proteoform_ids)) if len(proteoform_ids) == 1 else gene_id
        members = {
            member
            for value in group.get("Member_Ensembl_Gene_IDs", pd.Series(dtype=object))
            if pd.notna(value)
            for member in str(value).split(";")
            if member
        }
        row["Member_Ensembl_Gene_IDs"] = ";".join(sorted(members)) or gene_id
        for column in _REFERENCE_REQUEST_COLUMNS:
            if column not in group.columns:
                continue
            if column == "available":
                row[column] = bool(group[column].fillna(False).all())
            else:
                row[column] = _join_nonempty(group[column])
        for column in _ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS:
            if column not in group.columns:
                continue
            if column in {"is_filterable_extra", "is_technical_extra", "is_missing_biological"}:
                row[column] = any(_optional_bool(value) is True for value in group[column])
            else:
                row[column] = _join_nonempty(group[column])

        expr = pd.to_numeric(physical["expression"], errors="coerce")
        weights = pd.to_numeric(
            physical.get(
                "n_reference_samples",
                physical.get("n_samples", pd.Series(1, index=physical.index)),
            ),
            errors="coerce",
        ).fillna(0.0)
        valid = expr.notna() & (weights > 0)
        if valid.any():
            row["expression"] = float(np.average(expr[valid], weights=weights[valid]))
            pooled_n = float(weights[valid].sum())
        else:
            row["expression"] = np.nan
            pooled_n = 0.0
        row["q1"] = np.nan
        row["q3"] = np.nan
        row["source_cohort"] = "POOLED"
        row["source_project"] = "pooled_source_union"
        row["source_version"] = _join_nonempty(
            physical.get("source_version", pd.Series(dtype=object))
        )
        row["source_pmid"] = _join_nonempty(physical.get("source_pmid", pd.Series(dtype=object)))
        row["source_type"] = "pooled"
        row["source_unit"] = _first_nonempty(physical.get("source_unit", pd.Series(dtype=object)))
        scale_classes = {
            str(v)
            for v in physical.get("source_scale_class", pd.Series(dtype=object))
            if pd.notna(v) and str(v)
        }
        row["source_scale_class"] = scale_classes.pop() if len(scale_classes) == 1 else "mixed"
        comparable = physical.get("linear_tpm_comparable", pd.Series(False, index=physical.index))
        row["linear_tpm_comparable"] = bool(pd.Series(comparable).fillna(False).all())
        row["tumor_origin"] = _join_nonempty(physical.get("tumor_origin", pd.Series(dtype=object)))
        row["metastasis_site"] = _join_nonempty(
            physical.get("metastasis_site", pd.Series(dtype=object))
        )
        row["n_reference_genes"] = 1
        row["n_reference_samples"] = pooled_n
        row["n_samples"] = pooled_n
        if "n_detected" in physical.columns:
            row["n_detected"] = float(pd.to_numeric(physical["n_detected"], errors="coerce").sum())
        else:
            row["n_detected"] = np.nan
        row["processing_pipeline"] = "pooled_n_weighted"
        row["notes"] = _join_nonempty(physical.get("notes", pd.Series(dtype=object)))
        row["reference_method"] = "pooled_source_summary_rows"
        row["sample_qc"] = _first_nonempty(physical.get("sample_qc", pd.Series(dtype=object)))
        row["sample_qc_effective"] = _first_nonempty(
            physical.get("sample_qc_effective", pd.Series(dtype=object))
        )
        row["data_version"] = DATA_VERSION
        row["source_matrix_version"] = SOURCE_MATRIX_VERSION
        out_rows.append(row)
    cols = [c for c in long.columns if c in out_rows[0]]
    return pd.DataFrame(out_rows, columns=cols)


def _unavailable_reference_summary_metadata() -> dict[str, str | bool | int | float | None]:
    return {
        "source_cohort": None,
        "source_project": None,
        "source_version": None,
        "source_pmid": None,
        "source_type": None,
        "unit": None,
        "source_scale_class": "unknown",
        "linear_tpm_comparable": False,
        "tumor_origin": None,
        "metastasis_site": None,
        "n_reference_genes": None,
        "n_reference_samples": None,
        "processing_pipeline": None,
        "notes": None,
        "tpm_proxy": False,
    }


def _reference_summary_source_metadata(
    code: str,
    *,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
) -> dict[str, str | bool | int | float | None] | None:
    """Read one source's provenance from the compact availability manifest."""
    sources = _filter_reference_summary_sources(
        _reference_summary_source_table(),
        source_kinds=source_kinds,
        source_cohorts=source_cohorts,
        exclude_microarray_proxy=exclude_microarray_proxy,
    )
    matches = sources.loc[sources["cancer_code"].astype(str) == str(code)]
    selected_matches = matches.loc[matches["selected"].fillna(False).astype(bool)]
    if not selected_matches.empty:
        matches = selected_matches
    selected = None if matches.empty else matches.iloc[0].to_dict()
    if selected is None:
        return None
    return _reference_summary_metadata_from_row(code, selected)


def _reference_summary_metadata_from_row(
    code: str, source: dict
) -> dict[str, str | bool | int | float | None]:
    """Combine one compact summary-source row with the source registry."""
    meta = cancer_reference_expression_source_metadata(
        code, source_cohort=str(source["source_cohort"])
    )
    manifest_pmid = source.get("source_pmid")
    if pd.isna(manifest_pmid) or not str(manifest_pmid):
        manifest_pmid = meta.get("source_pmid")
    manifest_type = source.get("source_type")
    if pd.isna(manifest_type) or not str(manifest_type):
        manifest_type = meta.get("source_type")
    manifest_scale = source.get("source_scale_class")
    if pd.isna(manifest_scale) or not str(manifest_scale):
        manifest_scale = meta.get("source_scale_class")
    manifest_comparable = source.get("linear_tpm_comparable")
    if pd.isna(manifest_comparable):
        manifest_comparable = meta.get("linear_tpm_comparable")
    text = " ".join(
        str(value) if pd.notna(value) else ""
        for value in (
            source.get("source_cohort"),
            source.get("source_project"),
            source.get("source_version"),
            source.get("tumor_origin"),
        )
    ).lower()
    scale_text = str(manifest_scale or "").lower()
    tpm_proxy = (
        bool(meta.get("tpm_proxy"))
        or "proxy" in scale_text
        or "microarray" in text
        or "tpm-proxy" in text
    )
    meta.update(
        {
            "source_cohort": source.get("source_cohort"),
            "source_project": source.get("source_project"),
            "source_version": source.get("source_version"),
            "source_pmid": manifest_pmid,
            "source_type": manifest_type
            or _source_cohort_kind_map().get(str(source.get("source_cohort"))),
            "tumor_origin": source.get("tumor_origin"),
            "metastasis_site": source.get("metastasis_site"),
            "n_reference_genes": source.get("n_reference_genes"),
            "n_reference_samples": source.get("n_reference_samples"),
            "processing_pipeline": source.get("processing_pipeline"),
            "notes": source.get("notes"),
            "unit": meta.get("unit") or "TPM",
            "source_scale_class": manifest_scale
            or ("microarray_tpm_proxy" if tpm_proxy else "linear_rnaseq_tpm"),
            "linear_tpm_comparable": False
            if tpm_proxy
            else _optional_bool(manifest_comparable) is True,
            "tpm_proxy": tpm_proxy,
        }
    )
    return meta


def _reference_expression_provenance(
    code: str,
    method: str,
    *,
    sample_qc: str,
    sample_qc_effective: str | None,
    reference_source: str,
    source_kinds: set[str] | None = None,
    source_cohorts: set[str] | None = None,
    exclude_microarray_proxy: bool = False,
) -> dict:
    if method == _REFERENCE_SUMMARY_ALL_METHOD:
        meta = (
            _reference_summary_source_metadata(
                code,
                source_kinds=source_kinds,
                source_cohorts=source_cohorts,
                exclude_microarray_proxy=exclude_microarray_proxy,
            )
            or _unavailable_reference_summary_metadata()
        )
    elif method == _REFERENCE_SUMMARY_METHOD:
        meta = _reference_summary_source_metadata(
            code
        ) or cancer_reference_expression_source_metadata(code)
    else:
        meta = cancer_reference_expression_source_metadata(code)
    return {
        "source_cohort": meta["source_cohort"]
        or (None if method == _REFERENCE_SUMMARY_ALL_METHOD else code),
        "source_project": meta.get("source_project"),
        "source_version": meta.get("source_version"),
        "source_pmid": meta.get("source_pmid"),
        "source_type": meta.get("source_type"),
        "source_unit": meta.get("unit"),
        "source_scale_class": meta.get("source_scale_class"),
        "linear_tpm_comparable": bool(meta.get("linear_tpm_comparable")),
        "tumor_origin": meta.get("tumor_origin"),
        "metastasis_site": meta.get("metastasis_site"),
        "n_reference_genes": meta.get("n_reference_genes"),
        "n_reference_samples": meta.get("n_reference_samples"),
        "n_samples": meta.get("n_reference_samples"),
        "processing_pipeline": meta.get("processing_pipeline"),
        "notes": meta.get("notes"),
        "reference_method": method,
        "sample_qc": sample_qc,
        "sample_qc_effective": sample_qc_effective,
        "data_version": DATA_VERSION,
        "source_matrix_version": SOURCE_MATRIX_VERSION,
    }


def _merge_cancer_reference_wide_parts(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge cohort reference vectors by canonical ENSG only.

    Cohorts can carry different release-era symbols for the same canonical gene.
    Joining on ``(Ensembl_Gene_ID, Symbol)`` fragments those loci into sparse rows.
    """
    if not parts:
        return pd.DataFrame(columns=["Ensembl_Gene_ID", "Symbol"])

    value_cols: list[str] = []
    symbols: dict[str, Counter] = defaultdict(Counter)
    keyed_parts: list[pd.DataFrame] = []
    for part in parts:
        value_cols.extend(c for c in part.columns if c not in {"Ensembl_Gene_ID", "Symbol"})
        for gid, symbol in zip(part["Ensembl_Gene_ID"].astype(str), part["Symbol"].astype(str)):
            symbols[gid][symbol] += 1
        keyed = part.drop(columns=["Symbol"], errors="ignore")
        if keyed["Ensembl_Gene_ID"].duplicated().any():
            keyed = keyed.groupby("Ensembl_Gene_ID", sort=False).first().reset_index()
        keyed_parts.append(keyed)

    alias_symbols = ensembl_id_alias_symbols()

    def _symbol(gid: str) -> str:
        auth = alias_symbols.get(gid)
        if auth:
            return auth
        named = Counter({s: c for s, c in symbols.get(gid, {}).items() if s and s != gid})
        return named.most_common(1)[0][0] if named else gid

    out = keyed_parts[0]
    for part in keyed_parts[1:]:
        out = out.merge(part, on="Ensembl_Gene_ID", how="outer")
    out = out.sort_values("Ensembl_Gene_ID").reset_index(drop=True)
    out.insert(1, "Symbol", out["Ensembl_Gene_ID"].astype(str).map(_symbol))
    return out[["Ensembl_Gene_ID", "Symbol", *value_cols]]


#: Per-gene cohort summary statistic -> output column. The percentiles are taken across
#: the cohort's samples (the same axis as :func:`cohort_gene_percentiles`).
_COHORT_STAT_PERCENTILES = {
    0: "min",
    1: "p1",
    5: "p5",
    10: "p10",
    15: "p15",
    20: "p20",
    25: "p25",
    50: "p50",
    75: "p75",
    80: "p80",
    85: "p85",
    90: "p90",
    95: "p95",
    99: "p99",
    100: "max",
}


def _write_cohort_stat_columns(out: pd.DataFrame, mat: np.ndarray) -> None:
    """Write the ``mean``/``std`` + percentile-ladder columns onto ``out`` from a
    ``(n_genes, n_samples)`` value matrix — **availability-aware** (``NaN`` cells are
    skipped, so each gene reduces only over the samples that measured it). ``std`` is
    ``NaN`` for a gene measured by fewer than two samples (a lone observation has no
    spread — the same ``n >= 2`` rule used for ``std_between``). Shared by
    :func:`cohort_stats` and :func:`pooled_cohort_stats` so the suite is defined once."""
    pcts = list(_COHORT_STAT_PERCENTILES)
    n_obs = np.sum(~np.isnan(mat), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN gene rows -> NaN
        out["mean"] = np.nanmean(mat, axis=1)
        out["std"] = np.where(n_obs >= 2, np.nanstd(mat, axis=1), np.nan)
        q = np.nanpercentile(mat, pcts, axis=1)  # (len(pcts), n_genes)
    for i, p in enumerate(pcts):
        out[_COHORT_STAT_PERCENTILES[p]] = q[i]


def cohort_stats(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    proteoform: bool = False,
    scope: str = "cta",
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Per-gene **summary statistics** across a cohort's samples, in one pass —
    ``mean``, ``std`` and a uniform percentile ladder ``min, p1, p5, p10, p15, p20,
    p25, p50, p75, p80, p85, p90, p95, p99, max`` (``p25``/``p50``/``p75`` are
    q1/median/q3).

    The richer companion to :func:`cohort_mean_expression` (a single statistic) and
    :func:`cohort_gene_percentiles` (the dense percentile vector): everything a consumer
    needs to describe a gene's distribution over the cohort in one frame. Computed on
    the per-sample matrix in the ``normalize`` space (linear clean TPM by default; pass
    ``"tpm_clean_log1p"`` for log-space stats).

    ``proteoform=True`` summarizes the reduced proteoform key space (members summed per
    sample first, ``scope`` ``"cta"``/``"genome"``) — a proteoform-level frame carrying
    ``proteoform_key`` (see :func:`oncoref.proteoforms.expression_level`). Returns the
    id columns plus one column per statistic. ``sample_qc`` defaults to ``"pass"``
    for live summaries; use ``"all"`` to include every source-matrix sample."""
    df = per_sample_expression(
        cancer_type,
        normalize=normalize,
        auto_fetch=auto_fetch,
        proteoform=proteoform,
        scope=scope,
        sample_qc=sample_qc,
    )
    id_cols = id_columns(df)
    samples = sample_columns(df)
    if not samples:
        raise ValueError(f"no per-sample columns to summarize for {cancer_type!r}")
    out = df[id_cols].copy()
    _write_cohort_stat_columns(out, df[samples].to_numpy(dtype=float))
    return out


def _canonicalize_gene_rows(df: pd.DataFrame, *, sample_cols=None) -> pd.DataFrame:
    """Re-key one cohort's rows onto the **canonical** gene id so a locus can't be
    fragmented across cohorts (the pirlygenes#465-class bug).

    Each ``Ensembl_Gene_ID`` is unversioned *and* migration-resolved through the shipped
    ensembl-id-aliases map (:func:`resolve_ensembl_id`), so an alt-haplotype / patch /
    retired id collapses onto its canonical primary-assembly id. When a cohort carries
    BOTH an alias id and its canonical sibling (a full-assembly quantification annotates
    the gene on the primary contig *and* its alt-haplotype copy), their per-sample TPMs
    are **summed** under the canonical id: RNA-seq reads multi-map between the copies, so
    each row individually under-counts the gene — the same rationale as proteoform
    summation, one level up (gene rather than protein). A cross-release retired id and its
    successor never co-occur in one sample, so summing them degenerates to the lone value
    (a relabel). All-``NaN`` cells stay ``NaN`` (``min_count=1``) so an unmeasured gene is
    never turned into a measured zero — the canonical symbol is taken from the
    primary-contig row (sorted first) deterministically. The fast path (no collisions)
    only rewrites the id column.

    ``sample_cols`` names the value columns to **sum** (the rest are kept via ``first``);
    summing must be done in **linear** TPM, so the caller transforms (log1p/hk) only
    afterwards. When ``None`` the value columns are inferred as the numeric dtypes — fine
    for the gene×sample frames here, but pass them explicitly at the per-sample chokepoint
    so a stray numeric id column can never be summed."""
    canon = df["Ensembl_Gene_ID"].astype(str).map(resolve_ensembl_id)
    if not canon.duplicated().any():
        return df.assign(Ensembl_Gene_ID=canon.to_numpy())
    orig = df["Ensembl_Gene_ID"].astype(str).map(unversioned)
    is_primary = orig.to_numpy() == canon.to_numpy()
    df = df.assign(Ensembl_Gene_ID=canon.to_numpy(), _primary=is_primary)
    df = df.sort_values("_primary", ascending=False, kind="stable").drop(columns="_primary")
    if sample_cols is None:
        sum_cols = df.select_dtypes("number").columns.tolist()
    else:
        sum_cols = [c for c in sample_cols if c in df.columns]
    keep_cols = [c for c in df.columns if c != "Ensembl_Gene_ID" and c not in sum_cols]
    grouped = df.groupby("Ensembl_Gene_ID", sort=False)
    parts = []
    if keep_cols:  # canonical symbol / id columns: keep the primary-contig row's value
        parts.append(grouped[keep_cols].first())
    if sum_cols:  # per-sample TPM: SUM alt-haplotype reads into the canonical gene
        parts.append(grouped[sum_cols].sum(min_count=1))
    out = pd.concat(parts, axis=1).reset_index()
    return out[list(df.columns)]


def representative_cohort_samples(
    cancer_types: str | Iterable[str] | None = None,
    *,
    k: int | None = None,
    normalize: str = "tpm_clean",
    format: str = "wide",
    include_provenance: bool = False,
    representative_id_style: str = "pirlygenes",
    gene_id_style: str = "oncoref",
    gene_universe: str = "artifact",
    sample_qc: str = "pass",
    include_gene_universe_flags: bool = False,
) -> pd.DataFrame:
    """Representative real per-sample expression vectors per cohort.

    The packaged cohort references are per-cohort aggregates; this returns a
    bounded set of real joint per-sample vectors per cohort — medoids spanning
    the within-cohort variation — in the same ``clean TPM`` basis.

    ``cancer_types`` accepts a code, alias, or iterable; a computed-aggregate
    code expands to its member subtypes; ``None`` returns every cohort that ships
    representatives. ``k`` keeps at most the first ``k`` reps per cohort.
    ``format`` is ``"wide"`` (genes × reps) or ``"long"``. Public representative
    IDs default to pirlygenes-compatible ``CODE_rep01`` columns/values; pass
    ``representative_id_style="internal"`` for the shard/provenance IDs
    (``CODE__rep1``). Long output can attach representative-level provenance:
    source cohort/project/sample and diagnosis, a stable source-group ID shared
    by aliases of the same physical vector, any reviewed role/evidence, selection
    rank/method/basis, and package/data schema versions. Row-level ``sample_qc``
    and ``source_sample_qc`` are the
    selected sample's actual status; ``sample_qc_requested`` and
    ``sample_qc_effective`` preserve the artifact policies separately.
    ``benchmark_eligible=False`` with a QC-fallback ``representative_role``
    keeps an all-fail representative inspectable without presenting it as
    ordinary validation truth.

    ``gene_id_style="oncoref"`` returns canonical oncoref ENSG IDs. Opt into
    ``"pirlygenes"`` only for migration wrappers that need known legacy ENSG IDs
    for rows recorded as ``remapped_to_oncoref`` in the expression artifact delta
    table; missing rows and values are not synthesized.
    ``gene_universe="artifact"`` preserves exact shipped rows.
    ``gene_universe="tumor_signal"`` drops filterable extras while retaining
    biological oncoref-only rows. ``gene_universe="pirlygenes"`` additionally
    drops audited oncoref-only biological/unresolved extras unless they are
    documented remap targets; ``include_gene_universe_flags=True`` appends row-level
    audit columns.
    ``sample_qc="pass"`` requires any shipped build metadata to show the
    representative shard was built from QC-passing samples. Pass
    ``sample_qc="artifact"`` only for explicit legacy/audit reads of whatever
    policy the bundle used.
    """
    _validate_representative_id_style(representative_id_style)
    _validate_gene_id_style(gene_id_style)
    gene_universe = _validate_artifact_gene_universe(gene_universe)
    sample_qc = _validate_artifact_sample_qc(sample_qc)
    if normalize not in ("tpm_clean", "tpm_clean_log1p"):
        raise ValueError(
            "representative_cohort_samples normalize must be 'tpm_clean' or "
            "'tpm_clean_log1p' (the artifact ships only in clean TPM)"
        )
    for code in _resolve_cancer_types(cancer_types, expand_aggregates=True) or []:
        _reject_invalid_source_code(code, context="representative samples")
    if format not in ("wide", "long"):
        raise ValueError("format must be 'wide' or 'long'")
    if include_provenance and format != "long":
        # Provenance is per-representative (one row each); it only attaches to the
        # long form. Fail loudly rather than silently dropping the request.
        raise ValueError("include_provenance=True requires format='long'")

    root = _shard_dir(_REPRESENTATIVES)
    available = set(available_representative_cohorts())
    if cancer_types is None:
        codes = sorted(available)
    else:
        requested = _resolve_cancer_types(cancer_types, expand_aggregates=True)
        codes = [c for c in dict.fromkeys(requested) if c in available]
    artifact_qc_meta = _require_expression_artifact_sample_qc(codes, sample_qc=sample_qc)

    base = ["Ensembl_Gene_ID", "Symbol"]
    # Combine cohorts on a CANONICAL gene id only — never the (Ensembl_Gene_ID, Symbol)
    # pair. Cohorts were quantified against different Ensembl releases, so the same locus
    # carries a different release-alias symbol per cohort; merging on the pair fragments one
    # gene into many mutually-disjoint sparse rows (the pirlygenes#465-class bug). We key
    # each cohort by the canonical gene id — unversioned AND migration-resolved through the
    # shipped ensembl-id-aliases map (resolve_ensembl_id), so an alt-haplotype/archived id
    # collapses onto its primary-contig id instead of standing as a separate row — then
    # resolve one canonical symbol per gene afterwards.
    wide_parts = []
    long_parts = []
    symbols: dict[str, Counter] = defaultdict(Counter)  # canonical id -> symbol counts
    for code in codes:
        shard = pd.read_parquet(root / f"{code}.parquet")
        rep_cols = sample_columns(shard)
        if k is not None:
            rep_cols = rep_cols[:k]
        # Key on the canonical gene id (unversion + migration-resolve), collapsing an
        # alt-haplotype row onto its primary-contig sibling within the cohort first. The
        # shard ships in LINEAR clean TPM, so the alt-copy sum happens in linear space;
        # the log1p transform (if any) is applied AFTER — never sum log-space values
        # (log1p(a)+log1p(b) != log1p(a+b)), mirroring the proteoform linear-then-log1p
        # contract in per_sample_expression.
        shard = _canonicalize_gene_rows(shard)
        if normalize == "tpm_clean_log1p":
            shard[rep_cols] = np.log1p(shard[rep_cols].to_numpy(dtype=float))
        gid = shard["Ensembl_Gene_ID"].astype(str)
        for g, s in zip(gid, shard["Symbol"].astype(str)):
            symbols[g][s] += 1
        display_cols = [
            _representative_id_for_style(c, style=representative_id_style) for c in rep_cols
        ]
        mat = shard[rep_cols].set_axis(gid.to_numpy(), axis=0).set_axis(display_cols, axis=1)
        if format == "wide":
            wide_parts.append(mat)
        else:
            internal_mat = shard[rep_cols].set_axis(gid.to_numpy(), axis=0)
            melted = internal_mat.reset_index(names="Ensembl_Gene_ID").melt(
                id_vars="Ensembl_Gene_ID", var_name="representative_id", value_name="expression"
            )
            melted.insert(1, "cancer_code", code)
            long_parts.append(melted)

    alias_symbols = ensembl_id_alias_symbols()

    def _canonical_symbol(gid: str) -> str:
        # The curated symbol for a migrated locus is authoritative and wins outright;
        # otherwise prefer a real name over the raw-ENSG backfill that release-unaware
        # cohorts carry (the most common alias across cohorts, deterministic), else the id.
        auth = alias_symbols.get(gid)
        if auth:
            return auth
        named = Counter({s: c for s, c in symbols.get(gid, {}).items() if s and s != gid})
        return named.most_common(1)[0][0] if named else gid

    if format == "wide":
        if not wide_parts:
            out = pd.DataFrame(columns=base)
            out.attrs.update(
                _representative_attrs(
                    codes=codes,
                    normalize=normalize,
                    format=format,
                    k=k,
                    representative_id_style=representative_id_style,
                    gene_id_style=gene_id_style,
                    gene_universe=gene_universe,
                    sample_qc=sample_qc,
                    artifact_qc_meta=artifact_qc_meta,
                )
            )
            out = _apply_artifact_gene_universe(
                out,
                product="representative_cohort_samples",
                cancer_codes=codes,
                gene_universe=gene_universe,
                include_gene_universe_flags=include_gene_universe_flags,
            )
            _attach_gene_universe_delta_attrs(
                out, product="representative_cohort_samples", cancer_codes=codes
            )
            return out
        combined = pd.concat(wide_parts, axis=1, join="outer").sort_index()
        if combined.index.has_duplicates:  # belt-and-suspenders: one row per gene id
            combined = combined.groupby(level=0).first()
        out = combined.reset_index(names="Ensembl_Gene_ID")
        out.insert(1, "Symbol", out["Ensembl_Gene_ID"].map(_canonical_symbol))
        out = _apply_artifact_gene_universe(
            out,
            product="representative_cohort_samples",
            cancer_codes=codes,
            gene_universe=gene_universe,
            include_gene_universe_flags=include_gene_universe_flags,
        )
        out = _apply_gene_id_style(
            out,
            product="representative_cohort_samples",
            cancer_codes=codes,
            gene_id_style=gene_id_style,
            alias_expand_remaps=gene_universe == "pirlygenes",
        )
        out.attrs.update(
            _representative_attrs(
                codes=codes,
                normalize=normalize,
                format=format,
                k=k,
                representative_id_style=representative_id_style,
                gene_id_style=gene_id_style,
                gene_universe=gene_universe,
                sample_qc=sample_qc,
                artifact_qc_meta=artifact_qc_meta,
            )
        )
        _attach_gene_universe_delta_attrs(
            out, product="representative_cohort_samples", cancer_codes=codes
        )
        return out

    if not long_parts:
        out = _representative_empty_frame(include_provenance=include_provenance)
        out.attrs.update(
            _representative_attrs(
                codes=codes,
                normalize=normalize,
                format=format,
                k=k,
                representative_id_style=representative_id_style,
                gene_id_style=gene_id_style,
                gene_universe=gene_universe,
                sample_qc=sample_qc,
                artifact_qc_meta=artifact_qc_meta,
            )
        )
        out = _apply_artifact_gene_universe(
            out,
            product="representative_cohort_samples",
            cancer_codes=codes,
            gene_universe=gene_universe,
            include_gene_universe_flags=include_gene_universe_flags,
        )
        _attach_gene_universe_delta_attrs(
            out, product="representative_cohort_samples", cancer_codes=codes
        )
        return out
    long = pd.concat(long_parts, ignore_index=True)
    long.insert(1, "Symbol", long["Ensembl_Gene_ID"].map(_canonical_symbol))
    if include_provenance:
        long = _attach_representative_provenance(long, root)
    long["representative_id"] = long["representative_id"].map(
        lambda x: _representative_id_for_style(x, style=representative_id_style)
    )
    long = _apply_artifact_gene_universe(
        long,
        product="representative_cohort_samples",
        cancer_codes=codes,
        gene_universe=gene_universe,
        include_gene_universe_flags=include_gene_universe_flags,
    )
    long = _apply_gene_id_style(
        long,
        product="representative_cohort_samples",
        cancer_codes=codes,
        gene_id_style=gene_id_style,
        alias_expand_remaps=gene_universe == "pirlygenes",
    )
    long.attrs.update(
        _representative_attrs(
            codes=codes,
            normalize=normalize,
            format=format,
            k=k,
            representative_id_style=representative_id_style,
            gene_id_style=gene_id_style,
            gene_universe=gene_universe,
            sample_qc=sample_qc,
            artifact_qc_meta=artifact_qc_meta,
        )
    )
    _attach_gene_universe_delta_attrs(
        long, product="representative_cohort_samples", cancer_codes=codes
    )
    return long


def _biological_per_sample(
    code, *, proteoform: bool, auto_fetch: bool, scope: str = "cta", sample_qc: str = "pass"
) -> pd.DataFrame:
    """Clean-TPM per-sample matrix with technical/censored genes dropped — the
    biological view the summary artifacts are built on — collapsed to proteoform level
    (in ``scope``) when requested. The runtime input to the percentile / within-sample
    build cores, so a summary can be recomputed on the fly from the per-sample matrix
    (no shard). ``scope`` is ignored when ``proteoform`` is False."""
    from .gene_families import clean_tpm_censored_gene_ids

    clean = per_sample_expression(
        code, normalize="tpm_clean", auto_fetch=auto_fetch, sample_qc=sample_qc
    )
    censored = clean_tpm_censored_gene_ids()
    unversioned = clean["Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    bio = clean[~unversioned.isin(censored)].reset_index(drop=True)
    if proteoform:
        from .proteoforms import collapse_to_proteoforms

        bio = collapse_to_proteoforms(bio, scope=scope, sample_cols=sample_columns(bio))
    return bio


def _read_shard_or_recompute(
    dataset: ShardDataset,
    code: str,
    *,
    proteoform: bool,
    auto_fetch: bool,
    scope: str = "cta",
    sample_qc: str = "pass",
    shard: Path | None = None,
) -> pd.DataFrame:
    """Read ``code``'s shard for ``dataset`` (the ``scope``-specific one at proteoform
    level); if no current shard is present, recompute it from the per-sample matrix
    via the dataset's ``expression_builders`` core (the same core that produced the shipped shards —
    so the on-the-fly and shipped values agree).

    The single home of the shard-or-recompute fallback shared by the percentile and
    within-sample readers. Proteoform shards must identify the current registry;
    legacy/unverified or stale shards cannot be repaired by adding their summary
    values. Raises a clear :class:`ValueError` when neither a usable shard nor its
    per-sample matrix is available. The returned attrs identify the actual source.
    """
    if shard is None:
        shard = _shard_path(dataset, code, proteoform=proteoform, scope=scope)
    stale_registry = False
    if shard.exists():
        df = pd.read_parquet(shard)
        if proteoform:
            from .proteoforms import _proteoform_registry_metadata

            metadata = _proteoform_registry_metadata(scope)
            stale_registry = any(df.attrs.get(key) != value for key, value in metadata.items())
        if df.attrs.get("artifact_status") == "unavailable":
            raise ValueError(
                f"{code}: stored {dataset.noun} unavailable: {df.attrs.get('unavailable_reason', 'source validation failed')}"
            )
        if not stale_registry:
            df.attrs["source"] = "shard"
            return df
    effective_sample_qc = sample_qc
    if sample_qc == "artifact":
        meta = expression_artifact_build_metadata([code], auto_fetch=False, on_missing="empty")
        if meta.empty and meta.attrs.get("missing_reason"):
            # Legacy bundles without build metadata retain the strict fallback.
            effective_sample_qc = "pass"
        else:
            effective_sample_qc, reason = _artifact_effective_sample_qc(
                code, _artifact_build_metadata_qc_policies(meta)
            )
            if reason:
                raise ValueError(f"cannot recompute {dataset.noun} for {code!r}: {reason}")
    try:
        bio = _biological_per_sample(
            code,
            proteoform=proteoform,
            auto_fetch=auto_fetch,
            scope=scope,
            sample_qc=effective_sample_qc,
        )
    except FileNotFoundError as e:
        variant = "proteoform-summed " if proteoform else ""
        reason = (
            " The cached shard uses a stale or unverified proteoform registry and must "
            "be rebuilt from per-sample expression."
            if stale_registry
            else ""
        )
        raise ValueError(
            f"no {variant}{dataset.noun} for {code!r} and its per-sample matrix isn't "
            f"cached — fetch it (source_matrices.fetch / auto_fetch=True).{reason}"
        ) from e
    from importlib import import_module

    build_core = getattr(import_module("oncoref.expression_builders"), dataset.build_attr)
    out = build_core(bio, sample_columns(bio))
    out.attrs["source"] = "recomputed"
    out.attrs["sample_qc_effective"] = effective_sample_qc
    return out


def cohort_gene_percentiles(
    cancer_type,
    *,
    as_tpm: bool = True,
    proteoform: bool = False,
    scope: str = "cta",
    auto_fetch: bool = False,
    include_provenance: bool = False,
    on_missing: str = "raise",
    gene_id_style: str = "oncoref",
    gene_universe: str = "artifact",
    sample_qc: str = "pass",
    include_gene_universe_flags: bool = False,
) -> pd.DataFrame:
    """Tail-weighted per-gene percentile vector for one cohort.

    Returns one row per gene (``Ensembl_Gene_ID`` + ``Symbol``) with 26
    breakpoint columns — ``p0, p1, p5, p10 … p90, p95, p96, p97, p98, p99,
    p100`` — dense in the actionable upper tail. Lets a consumer place a sample's
    gene as a **percentile rank within the cohort** instead of an absolute TPM.

    Computed on the biological clean TPM view (technical genes dropped).
    Stored compactly as ``log1p`` + float16; ``as_tpm=True`` (default)
    ``expm1``-restores clean-TPM values, ``as_tpm=False`` returns the stored
    log1p values.

    ``include_provenance=True`` appends stable artifact/provenance columns. Missing
    cohorts still raise by default; pass ``on_missing="empty"`` to return an
    empty schema-stable frame with ``attrs["missing_reason"]`` instead.

    ``gene_id_style="oncoref"`` returns canonical oncoref ENSG IDs. Opt into
    ``"pirlygenes"`` only for gene-level migration wrappers that need known
    legacy ENSG IDs for rows recorded as ``remapped_to_oncoref`` in the
    expression artifact delta table; missing rows and values are not synthesized.
    ``gene_universe="artifact"`` preserves exact shipped rows.
    ``gene_universe="tumor_signal"`` drops filterable extras while retaining
    biological oncoref-only rows. ``gene_universe="pirlygenes"`` additionally
    drops audited oncoref-only biological/unresolved extras unless they are
    documented remap targets; ``include_gene_universe_flags=True`` appends row-level
    audit columns.
    ``sample_qc="pass"`` requires any shipped build metadata to show the
    percentile shard was built from QC-passing samples. Pass
    ``sample_qc="artifact"`` only for explicit legacy/audit reads of whatever
    policy the bundle used. Recomputed artifact reads apply the recorded effective
    policy; bundles without build metadata fall back to ``"pass"``.

    With ``proteoform=True``, the vector is one row per proteoform key
    (``proteoform_key``/``Symbol`` carry the collapsed identity), identical-protein
    members summed **before** the percentiles are computed (``scope`` ``"cta"``/``"genome"``,
    ignored when ``proteoform`` is False).

    The shipped percentile **shard** can't be converted across proteoform scopes (you
    can't sum already-computed percentiles). When no matching shard verified against
    the current registry is present, the vector is **recomputed on the fly** from
    the per-sample matrix via the same build core. That needs the cohort's per-sample
    matrix cached (pass ``auto_fetch=True`` to download it); otherwise a clear error.
    """
    if on_missing not in ("raise", "empty"):
        raise ValueError("on_missing must be 'raise' or 'empty'")
    _validate_gene_id_style(gene_id_style)
    gene_universe = _validate_artifact_gene_universe(gene_universe)
    sample_qc = _validate_artifact_sample_qc(sample_qc)
    if proteoform and gene_id_style != "oncoref":
        raise ValueError("gene_id_style='pirlygenes' is only supported for gene-level artifacts")
    if proteoform and (gene_universe != "artifact" or include_gene_universe_flags):
        raise ValueError(
            "gene_universe filtering/flags are only supported for gene-level artifacts"
        )

    code = resolve_cancer_type(cancer_type)
    shard = _shard_path(_PERCENTILES, code, proteoform=proteoform, scope=scope)
    artifact_qc_meta = _empty_metadata_frame(
        _EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS,
        schema_version=EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
        missing_reason="no usable percentile artifact shard; recomputing from source matrix",
    )
    try:
        df = _read_shard_or_recompute(
            _PERCENTILES,
            code,
            proteoform=proteoform,
            auto_fetch=auto_fetch,
            scope=scope,
            sample_qc=sample_qc,
            shard=shard,
        )
        source = df.attrs["source"]
        if source == "shard":
            artifact_qc_meta = _require_expression_artifact_sample_qc([code], sample_qc=sample_qc)
    except ValueError as e:
        if on_missing != "empty" or "per-sample matrix isn't cached" not in str(e):
            raise
        out = _empty_percentile_frame(
            code=code,
            as_tpm=as_tpm,
            proteoform=proteoform,
            scope=scope,
            include_provenance=include_provenance,
            gene_id_style=gene_id_style,
            gene_universe=gene_universe,
            sample_qc=sample_qc,
            artifact_qc_meta=artifact_qc_meta,
            include_gene_universe_flags=include_gene_universe_flags,
            missing_reason=str(e),
        )
        if not proteoform:
            _attach_gene_universe_delta_attrs(
                out, product="cohort_gene_percentiles", cancer_codes=[code]
            )
        return out
    bp_cols = sample_columns(df)
    df[bp_cols] = df[bp_cols].astype("float32")
    if as_tpm:
        df[bp_cols] = np.expm1(df[bp_cols])
    if include_provenance:
        df = _attach_percentile_provenance(
            df,
            code=code,
            as_tpm=as_tpm,
            sample_qc=sample_qc,
            artifact_qc_meta=artifact_qc_meta,
        )
    df = _apply_artifact_gene_universe(
        df,
        product="cohort_gene_percentiles",
        cancer_codes=[code],
        gene_universe=gene_universe,
        include_gene_universe_flags=include_gene_universe_flags,
    )
    df = _apply_gene_id_style(
        df,
        product="cohort_gene_percentiles",
        cancer_codes=[code],
        gene_id_style=gene_id_style,
        alias_expand_remaps=gene_universe == "pirlygenes",
    )
    df.attrs.update(
        _percentile_attrs(
            code=code,
            as_tpm=as_tpm,
            proteoform=proteoform,
            scope=scope,
            source=source,
            gene_id_style=gene_id_style,
            gene_universe=gene_universe,
            sample_qc=sample_qc,
            artifact_qc_meta=artifact_qc_meta,
        )
    )
    if not proteoform:
        _attach_gene_universe_delta_attrs(
            df, product="cohort_gene_percentiles", cancer_codes=[code]
        )
    return df


# ---------- within-sample percentile prevalence (signal a) ----------

#: within-sample percentile-rank threshold -> output column: ``_WITHIN_SAMPLE_THRESHOLD_COLS``
#: (imported above from :data:`oncoref.expression_builders.WITHIN_SAMPLE_THRESHOLDS`) is the single
#: source of truth shared with the generator, so the read side and write side can't drift.


def available_within_sample_cohorts(*, proteoform: bool = False, scope: str = "cta") -> list[str]:
    """Cohort codes that ship a within-sample top-fraction shard (sorted).

    With ``proteoform=True``, the proteoform-summed variant (identical-protein
    members collapsed before ranking, in ``scope`` — see :func:`within_sample_top_fraction`)."""
    return _available_cohorts(_WITHIN_SAMPLE, proteoform=proteoform, scope=scope)


def locally_available_within_sample_cohorts(
    *,
    proteoform: bool = False,
    scope: str = "cta",
    include_recomputable: bool = True,
) -> list[str]:
    """Cohorts readable as within-sample prevalence without a download.

    The result is the union of shards in the package/artifact checkout and the
    active partial bundle cache.  With ``include_recomputable=True`` (default),
    it also includes cohorts whose per-sample source matrix is cached locally,
    because :func:`within_sample_top_fraction` can derive the same artifact from
    that matrix.  This function never fetches either data source.
    """
    return _locally_available_summary_cohorts(
        _WITHIN_SAMPLE,
        proteoform=proteoform,
        scope=scope,
        include_recomputable=include_recomputable,
    )


def within_sample_top_fraction(
    cancer_type,
    *,
    threshold: float = 0.95,
    proteoform: bool = False,
    scope: str = "cta",
    auto_fetch: bool = False,
) -> pd.DataFrame:
    """Per-gene fraction of a cohort's samples in which the gene is highly
    expressed *within that sample* — the "top ~5% expressed gene in this tumor"
    prevalence (signal a, the producer side of the within-sample signal).

    Returns one row per gene (``Ensembl_Gene_ID`` + ``Symbol``) with the
    ``frac_samples_top{1,5,10}pct`` column for the requested ``threshold``
    (0.99 / 0.95 / 0.90) plus ``n_samples``.

    With ``proteoform=True``, identical-protein paralogs (CTAG1A+CTAG1B, the CT47A
    family, …) are summed per sample *before* the within-sample ranking, so a
    duplicated antigen is ranked as one proteoform rather than several individually-
    diluted genes (``proteoform_key``/``Symbol`` carry the collapsed identity;
    ``scope`` ``"cta"``/``"genome"``, ignored when ``proteoform`` is False). Note
    collapsing members shrinks the gene axis the within-sample rank is computed over,
    so an ungrouped gene's fraction can shift slightly vs the gene variant.

    Reads the shipped shard when present and, for proteoforms, verified against the
    current registry. Otherwise **recomputes on the fly** from the per-sample matrix
    via the same build core and the bundle's recorded effective sample-QC policy
    (``"pass"`` when build metadata is unavailable). Recompute needs the cohort's per-sample
    matrix cached (pass ``auto_fetch=True`` to download it), else a clear error.
    """
    col = _WITHIN_SAMPLE_THRESHOLD_COLS.get(threshold)
    if col is None:
        raise ValueError(f"threshold must be one of {sorted(_WITHIN_SAMPLE_THRESHOLD_COLS)}")
    code = resolve_cancer_type(cancer_type)
    df = _read_shard_or_recompute(
        _WITHIN_SAMPLE,
        code,
        proteoform=proteoform,
        auto_fetch=auto_fetch,
        scope=scope,
        sample_qc="artifact",
    )
    keep = [*id_columns(df), col]
    if "n_samples" in df.columns:
        keep.append("n_samples")
    return df[keep]


# ---------- proteoform-level summation (identical-protein paralogs) ----------


def proteoform_representative_samples(
    cancer_types: str | Iterable[str] | None = None,
    *,
    k: int | None = None,
) -> pd.DataFrame:
    """Representative per-sample vectors with identical-protein genes summed to
    proteoform level.

    Same per-sample medoid vectors as :func:`representative_cohort_samples`
    (``format="wide"``), but the member genes of each proteoform group
    (CTAG1A+CTAG1B → ``CTAG1A/CTAG1B``, SSX4+SSX4B → ``SSX4/SSX4B``, the 12-member
    CT47A family, …) are *summed* per sample — the multi-mapping-correct unit
    when reads can't be uniquely assigned between identical-protein loci. Genes in
    no group pass through unchanged.

    Always operates on linear ``clean TPM`` (summing log1p values would be
    wrong); ``log1p`` afterward if you need it. This is the runtime,
    every-cohort proteoform view over the shipped medoid samples; the same
    :func:`oncoref.expression_builders.sum_proteoform_tpm` core can run inside the offline
    percentile/within-sample generators to ship proteoform-summed artifacts.
    """
    from .proteoforms import collapse_to_proteoforms

    wide = representative_cohort_samples(cancer_types, k=k, normalize="tpm_clean", format="wide")
    sample_cols = [c for c in wide.columns if c not in ("Ensembl_Gene_ID", "Symbol")]
    if wide.empty or not sample_cols:
        return wide
    return collapse_to_proteoforms(wide, sample_cols=sample_cols)


# ----- Symmetric gene-level / proteoform-level expression accessors -----
#
# Every expression dataset is exposed as a matched pair whose name carries the level —
# no boolean flag to read. ``gene_*`` is one row per Ensembl gene; ``proteoform_*``
# sums identical-protein paralogs to one row per proteoform key (see
# :func:`oncoref.proteoforms.expression_level`). Both are thin wrappers over the one
# base accessor, so the collapse logic lives in exactly one place; the unprefixed base
# names (``per_sample_expression`` etc.) remain the gene-level implementation.
#
#   dataset          gene_*                            proteoform_*
#   per-sample TPM   gene_per_sample_expression        proteoform_per_sample_expression
#   cohort-mean TPM  gene_cohort_mean_expression       proteoform_cohort_mean_expression
#   percentiles      gene_cohort_percentiles           proteoform_cohort_percentiles
#   within-sample    gene_within_sample_top_fraction   proteoform_within_sample_top_fraction
#   representatives  gene_representative_samples       proteoform_representative_samples


def gene_per_sample_expression(
    cancer_type, *, normalize: str = "tpm_clean", auto_fetch: bool = True, sample_qc: str = "all"
) -> pd.DataFrame:
    """Gene-level per-sample **TPM values** (one row per Ensembl gene). Proteoform
    counterpart: :func:`proteoform_per_sample_expression`."""
    return per_sample_expression(
        cancer_type, normalize=normalize, auto_fetch=auto_fetch, sample_qc=sample_qc
    )


def proteoform_per_sample_expression(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    scope: str = "cta",
    sample_qc: str = "all",
) -> pd.DataFrame:
    """Proteoform-level per-sample **TPM values** — identical-protein paralogs summed
    per sample. Gene-level counterpart: :func:`gene_per_sample_expression`."""
    return per_sample_expression(
        cancer_type,
        normalize=normalize,
        auto_fetch=auto_fetch,
        proteoform=True,
        scope=scope,
        sample_qc=sample_qc,
    )


def gene_cohort_mean_expression(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    statistic: str = "mean",
    auto_fetch: bool = True,
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Gene-level across-patient **TPM** summary. Proteoform counterpart:
    :func:`proteoform_cohort_mean_expression`."""
    return cohort_mean_expression(
        cancer_type,
        normalize=normalize,
        statistic=statistic,
        auto_fetch=auto_fetch,
        sample_qc=sample_qc,
    )


def proteoform_cohort_mean_expression(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    statistic: str = "mean",
    auto_fetch: bool = True,
    scope: str = "cta",
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Proteoform-level across-patient **TPM** summary. Gene-level counterpart:
    :func:`gene_cohort_mean_expression`."""
    return cohort_mean_expression(
        cancer_type,
        normalize=normalize,
        statistic=statistic,
        auto_fetch=auto_fetch,
        proteoform=True,
        scope=scope,
        sample_qc=sample_qc,
    )


def gene_cohort_stats(
    cancer_type, *, normalize: str = "tpm_clean", auto_fetch: bool = True, sample_qc: str = "pass"
) -> pd.DataFrame:
    """Gene-level per-gene cohort **summary statistics** (mean/std + the percentile ladder
    min/p1/p5/p10/p15/p20/p25/p50/p75/p80/p85/p90/p95/p99/max). Proteoform counterpart:
    :func:`proteoform_cohort_stats`."""
    return cohort_stats(
        cancer_type, normalize=normalize, auto_fetch=auto_fetch, sample_qc=sample_qc
    )


def proteoform_cohort_stats(
    cancer_type,
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    scope: str = "cta",
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Proteoform-level per-gene cohort **summary statistics**. Gene-level counterpart:
    :func:`gene_cohort_stats`."""
    return cohort_stats(
        cancer_type,
        normalize=normalize,
        auto_fetch=auto_fetch,
        proteoform=True,
        scope=scope,
        sample_qc=sample_qc,
    )


def gene_pooled_cohort_stats(
    cancer_types,
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    min_cohorts: int = 1,
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Gene-level heterogeneity-safe cross-cohort pool. Proteoform counterpart:
    :func:`proteoform_pooled_cohort_stats`. (Alias of :func:`pooled_cohort_stats`.)"""
    return pooled_cohort_stats(
        cancer_types,
        normalize=normalize,
        auto_fetch=auto_fetch,
        min_cohorts=min_cohorts,
        sample_qc=sample_qc,
    )


def proteoform_pooled_cohort_stats(
    cancer_types,
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    scope: str = "cta",
    min_cohorts: int = 1,
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """Proteoform-level heterogeneity-safe cross-cohort pool. Gene-level counterpart:
    :func:`gene_pooled_cohort_stats`."""
    return pooled_cohort_stats(
        cancer_types,
        normalize=normalize,
        auto_fetch=auto_fetch,
        proteoform=True,
        scope=scope,
        min_cohorts=min_cohorts,
        sample_qc=sample_qc,
    )


def gene_cohort_percentiles(
    cancer_type,
    *,
    as_tpm: bool = True,
    auto_fetch: bool = False,
    include_provenance: bool = False,
    on_missing: str = "raise",
    gene_id_style: str = "oncoref",
) -> pd.DataFrame:
    """Gene-level per-cohort **percentile vectors**. Proteoform counterpart:
    :func:`proteoform_cohort_percentiles`. (Alias of :func:`cohort_gene_percentiles`.)"""
    return cohort_gene_percentiles(
        cancer_type,
        as_tpm=as_tpm,
        auto_fetch=auto_fetch,
        include_provenance=include_provenance,
        on_missing=on_missing,
        gene_id_style=gene_id_style,
    )


def proteoform_cohort_percentiles(
    cancer_type,
    *,
    as_tpm: bool = True,
    scope: str = "cta",
    auto_fetch: bool = True,
    include_provenance: bool = False,
    on_missing: str = "raise",
) -> pd.DataFrame:
    """Proteoform-level per-cohort **percentile vectors** (members summed before
    ranking, ``scope`` ``"cta"``/``"genome"``). Gene-level counterpart:
    :func:`gene_cohort_percentiles`. No proteoform shard ships yet, so this **always**
    recomputes from the per-sample matrix — hence ``auto_fetch`` defaults to ``True``
    here (unlike the gene variant, which reads a shipped shard); pass ``False`` to
    require the matrix already be cached."""
    return cohort_gene_percentiles(
        cancer_type,
        as_tpm=as_tpm,
        proteoform=True,
        scope=scope,
        auto_fetch=auto_fetch,
        include_provenance=include_provenance,
        on_missing=on_missing,
    )


def gene_within_sample_top_fraction(
    cancer_type, *, threshold: float = 0.95, auto_fetch: bool = False
) -> pd.DataFrame:
    """Gene-level within-sample top-fraction prevalence. Proteoform counterpart:
    :func:`proteoform_within_sample_top_fraction`."""
    return within_sample_top_fraction(cancer_type, threshold=threshold, auto_fetch=auto_fetch)


def proteoform_within_sample_top_fraction(
    cancer_type, *, threshold: float = 0.95, scope: str = "cta", auto_fetch: bool = True
) -> pd.DataFrame:
    """Proteoform-level within-sample top-fraction prevalence (``scope``
    ``"cta"``/``"genome"``). Gene-level counterpart: :func:`gene_within_sample_top_fraction`.
    No proteoform shard ships yet, so this **always** recomputes from the per-sample
    matrix — hence ``auto_fetch`` defaults to ``True`` here (unlike the gene variant,
    which reads a shipped shard); pass ``False`` to require the matrix already be cached."""
    return within_sample_top_fraction(
        cancer_type, threshold=threshold, proteoform=True, scope=scope, auto_fetch=auto_fetch
    )


def gene_representative_samples(
    cancer_types: str | Iterable[str] | None = None,
    *,
    k: int | None = None,
    normalize: str = "tpm_clean",
    format: str = "wide",
    include_provenance: bool = False,
    gene_id_style: str = "oncoref",
) -> pd.DataFrame:
    """Gene-level representative per-sample vectors. Proteoform counterpart:
    :func:`proteoform_representative_samples`. (Alias of
    :func:`representative_cohort_samples`.)"""
    return representative_cohort_samples(
        cancer_types,
        k=k,
        normalize=normalize,
        format=format,
        include_provenance=include_provenance,
        gene_id_style=gene_id_style,
    )


def pooled_cohort_stats(
    cancer_types: str | Iterable[str],
    *,
    normalize: str = "tpm_clean",
    auto_fetch: bool = True,
    proteoform: bool = False,
    scope: str = "cta",
    min_cohorts: int = 1,
    sample_qc: str = "pass",
) -> pd.DataFrame:
    """**Heterogeneity-safe** per-gene summary pooled across several cohorts.

    The cross-cohort companion to the single-cohort :func:`cohort_stats`. Pools the
    requested cohorts' per-sample matrices at the **sample** level on the shared key
    (``Ensembl_Gene_ID``, or ``proteoform_key`` when ``proteoform=True``) and
    summarizes each gene over the union of samples — **availability-aware**: a cell
    a cohort never measured is ``NaN`` and never treated as a zero, so the per-gene
    denominator ``n_available`` (not the constant ``n_samples``) is the honest one.

    Returns an id-keyed frame with the same statistic suite as :func:`cohort_stats`
    (``mean, std, min, p1, p5, p10, p15, p20, p25, p50, p75, p80, p85, p90, p95, p99, max`` over the
    pooled samples) plus the pooling columns:

    - ``balanced_mean`` — the mean of each cohort's per-gene mean, **equal weight
      per cohort**. The heterogeneity-safe central value: a large cohort can't
      dominate it the way it dominates the sample-pooled ``mean``. The gap between
      ``mean`` and ``balanced_mean`` is itself the cross-cohort imbalance signal.
    - ``std_between`` — std *across* the per-cohort means (between-cohort spread:
      how differently the cancer types express the gene), ``NaN`` for a single cohort.
    - ``n_samples`` — total pooled sample columns (constant across genes).
    - ``n_available`` — per-gene count of samples that **measured** it (non-``NaN``).
    - ``n_detected`` — measured **and** ``> 0``.
    - ``n_cohorts`` — how many of the pooled cohorts measured the gene.

    ``min_cohorts`` drops genes measured by fewer than that many cohorts (default
    ``1`` keeps everything). ``sample_qc`` defaults to ``"pass"`` so QC-failed
    sample columns do not affect live pooled summaries. ``normalize`` selects the
    pooling space (linear clean TPM by default; see :func:`per_sample_expression`).
    ``proteoform=True`` pools the reduced proteoform key space (``scope`` ``"cta"``/
    ``"genome"``). An aggregate code (e.g. ``"SARC"``) expands to its member subtypes
    — pooling them is exactly what a rollup cohort means."""
    codes = _resolve_cancer_types(cancer_types, expand_aggregates=True)
    codes = list(dict.fromkeys(codes or []))
    if not codes:
        raise ValueError("pooled_cohort_stats needs at least one cancer type")

    key = "proteoform_key" if proteoform else "Ensembl_Gene_ID"
    sample_frames: list[pd.DataFrame] = []  # per cohort: key-indexed sample matrix
    cohort_means: list[pd.Series] = []  # per cohort: key -> per-gene mean
    id_rows: list[pd.DataFrame] = []  # per cohort: id columns, key-indexed
    # per_sample_expression already returns the dense CANONICAL space (alt-haplotype copies
    # summed in linear TPM, retired ids relabeled, transform applied after) so an alias id
    # can't stand as a separate sparse row here — no extra canonicalization needed.
    for code in codes:
        df = per_sample_expression(
            code,
            normalize=normalize,
            auto_fetch=auto_fetch,
            proteoform=proteoform,
            scope=scope,
            sample_qc=sample_qc,
        )
        id_cols = id_columns(df)
        samples = sample_columns(df)
        if not samples:
            continue
        indexed = df.set_index(key)
        mat = indexed[samples]
        # Sample labels can collide across cohorts ("s1"); namespace them so the
        # outer concat keeps every sample as its own column.
        mat = mat.add_prefix(f"{code}::")
        sample_frames.append(mat)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN gene -> NaN mean
            cohort_means.append(indexed[samples].mean(axis=1, skipna=True).rename(code))
        id_rows.append(df.set_index(key)[[c for c in id_cols if c != key]])
    if not sample_frames:
        raise ValueError(f"no per-sample columns to pool for {cancer_types!r}")

    pooled = pd.concat(sample_frames, axis=1, join="outer").sort_index()
    per_cohort_mean = pd.concat(cohort_means, axis=1).reindex(pooled.index)
    ids = pd.concat(id_rows).groupby(level=0).first().reindex(pooled.index)

    measured = pooled.notna()
    mat = pooled.to_numpy(dtype=float)
    out = ids.reset_index()
    _write_cohort_stat_columns(out, mat)
    cohort_mean_mat = per_cohort_mean.to_numpy(dtype=float)
    n_cohorts = per_cohort_mean.notna().to_numpy().sum(axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out["balanced_mean"] = np.nanmean(cohort_mean_mat, axis=1)
        # Between-cohort spread is undefined for a single measuring cohort (nanstd
        # would report a misleading 0) -> NaN there, mirroring std's >=2 rule.
        out["std_between"] = np.where(n_cohorts >= 2, np.nanstd(cohort_mean_mat, axis=1), np.nan)
    out["n_samples"] = pooled.shape[1]
    out["n_available"] = measured.to_numpy().sum(axis=1)
    out["n_detected"] = ((mat > 0) & measured.to_numpy()).sum(axis=1)
    out["n_cohorts"] = n_cohorts
    if min_cohorts > 1:
        out = out[out["n_cohorts"] >= min_cohorts].reset_index(drop=True)
    return out


from .expression_availability import (  # noqa: E402, F401
    _NON_COMPARABLE_PUBLIC_RNASEQ_AUDIT,
    _NON_COMPARABLE_REPLACEMENT_PLANS,
    _REFERENCE_AVAILABILITY_COLUMNS,
    _REFERENCE_NORMALIZE_ALIASES,
    _REFERENCE_NORMALIZE_LABELS,
    _REFERENCE_PROVENANCE_COLUMNS,
    _REFERENCE_RECOMPUTED_METHOD,
    _REFERENCE_REQUEST_COLUMNS,
    _REFERENCE_SOURCE_UNION_IDENTITY_COLUMNS,
    _REFERENCE_SOURCES,
    _REFERENCE_SUMMARY_ALL_METHOD,
    _REFERENCE_SUMMARY_AVAILABILITY_DATASET,
    _REFERENCE_SUMMARY_DATASET,
    _REFERENCE_SUMMARY_METHOD,
    _REFERENCE_TUMOR_ORIGIN_RANK,
    _REFERENCE_WIDE_SUFFIXES,
    _attach_reference_expression_attrs,
    _reference_availability_frame,
    _reference_available_codes_for_modes,
    _reference_expected_method,
    _reference_expression_availability_for_requests,
    _reference_expression_requests,
    _reference_mode_availability,
    _reference_normalize_modes,
    _reference_source_availability_row,
    _reference_source_union_availability,
    _reference_summary_availability_table,
    _source_matrix_effective_sample_count,
    _validate_reference_source,
    cancer_reference_expression_availability,
    non_comparable_expression_cohorts,
    umbrella_expression_coverage,
)
from .expression_pan_cancer import (  # noqa: E402, F401
    _PAN_CANCER_COMPUTED_AGGREGATES,
    _PAN_CANCER_POOLED_SERIES_CACHE,
    _add_computed_pan_cancer_raw_columns,
    _clear_pan_cancer_pooled_series_cache,
    _compute_pooled_reference_expression_series,
    _format_pan_cancer_columns,
    _pan_cancer_column_style,
    _pan_cancer_normalize_modes,
    _pooled_reference_expression_series_for_pan_cancer,
    _selected_reference_summary_rows_for_pan_cancer,
    pan_cancer_expression,
)
from .expression_qc import (  # noqa: E402, F401
    _ARTIFACT_SAMPLE_QC_MODES,
    _BLOCKING_SAMPLE_QC_REASONS,
    _EXPRESSION_ARTIFACT_BUILD_METADATA_COLUMNS,
    _HOUSEKEEPING_CANCER_COVERAGE_COLUMNS,
    _HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_COLUMNS,
    _HOUSEKEEPING_CANCER_COVERAGE_SUMMARY_REQUIRED,
    _PER_SAMPLE_NORMALIZE,
    _SAMPLE_QC_MODES,
    _SOURCE_MATRIX_SAMPLE_QC_MANIFEST_COLUMNS,
    DEFAULT_HOUSEKEEPING_DETECTION_FLOOR_TPM,
    DEFAULT_MAX_TOP10_GENE_FRACTION_FOR_QC,
    DEFAULT_MAX_TOP_GENE_FRACTION_FOR_QC,
    DEFAULT_MAX_ZERO_FRACTION_FOR_QC,
    DEFAULT_MIN_DETECTED_GENES_FOR_QC,
    DEFAULT_MIN_HOUSEKEEPING_GENES_FOR_QC,
    EXPRESSION_ARTIFACT_BUILD_METADATA_JSON_PATH,
    EXPRESSION_ARTIFACT_BUILD_METADATA_PATH,
    EXPRESSION_ARTIFACT_BUILD_METADATA_SCHEMA_VERSION,
    SAMPLE_EXPRESSION_QC_POLICY_VERSION,
    SOURCE_MATRIX_SAMPLE_QC_MANIFEST_PATH,
    _apply_sample_qc_filter,
    _artifact_build_metadata_qc_policies,
    _artifact_build_metadata_qc_value,
    _artifact_build_metadata_sample_counts,
    _artifact_effective_sample_qc,
    _artifact_qc_attrs,
    _clean_tpm_housekeeping_panel_ids,
    _coverage_quantile,
    _empty_metadata_frame,
    _housekeeping_panel_rows,
    _min_housekeeping_detected,
    _optional_bundle_metadata_path,
    _prepare_housekeeping_coverage_summary_input,
    _requested_cancer_codes,
    _require_expression_artifact_sample_qc,
    _sample_qc_status,
    _strict_coverage_bool_column,
    _summarize_housekeeping_gene_coverage,
    _validate_artifact_sample_qc,
    _validate_metadata_on_missing,
    _validate_sample_qc,
    cancer_reference_expression_source_metadata,
    expression_artifact_build_metadata,
    expression_artifact_build_summary,
    housekeeping_cancer_expression_coverage,
    housekeeping_cancer_expression_coverage_from_matrix,
    housekeeping_cancer_expression_coverage_summary,
    sample_expression_qc,
    sample_expression_qc_from_matrix,
    source_matrix_sample_qc_manifest,
)
from .expression_universe_deltas import (  # noqa: E402, F401
    _ARTIFACT_BIOLOGICAL_EXTRA_BIOTYPES,
    _ARTIFACT_BIOLOGICAL_EXTRA_STATUSES,
    _ARTIFACT_CANONICALIZATION_STATUSES,
    _ARTIFACT_FILTERABLE_EXTRA_CLASSES,
    _ARTIFACT_GENE_UNIVERSE_FLAG_COLUMNS,
    _ARTIFACT_GENE_UNIVERSE_MODES,
    _ARTIFACT_MISSING_BIOLOGICAL_STATUSES,
    _ARTIFACT_NON_SIGNAL_EXTRA_BIOTYPES,
    _ARTIFACT_NON_SIGNAL_EXTRA_STATUSES,
    _ARTIFACT_TECHNICAL_EXTRA_STATUSES,
    _DELTA_REPORT_COLUMNS,
    _annotate_expression_artifact_gene_universe_deltas,
    _apply_artifact_gene_universe,
    _artifact_delta_class,
    _artifact_delta_consumer_action,
    _artifact_row_key,
    _attach_gene_universe_delta_attrs,
    _delta_gene_biotype,
    _delta_nonempty,
    _delta_source_product,
    _expression_artifact_gene_universe_delta_report_for_codes,
    _filter_delta_rows_by_codes,
    _join_unique,
    _pirlygenes_compatible_row_ids,
    _resolve_delta_report_codes,
    _row_universe_annotation_frame,
    _validate_artifact_gene_universe,
    expression_artifact_gene_universe_delta_report,
    expression_artifact_gene_universe_delta_summary,
    expression_artifact_gene_universe_deltas,
    expression_artifact_technical_extra_gene_ids,
)
