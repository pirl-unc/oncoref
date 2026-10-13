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

"""Registry of every cohort-level expression data source.

The bundled ``expression_sources.yaml`` is the single source of truth for which
upstream datasets feed each cancer cohort — Treehouse, BEAT-AML, GDC/TCGA,
CLL-map, MMRF CoMMpass, TARGET, GEO series, recount3, … — with the provenance a
re-runner needs (project/accession/url, native unit, library prep, builder,
citation, expected size). It drives the legible ``oncoref expression-sources`` listing
and (incrementally) the fetch + regeneration of the per-cohort summaries.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import pandas as pd

from .expression_contract import EXPRESSION_CONTEXTS, ExpressionContract, ExpressionMeasurement
from .species import resolve_species

_DATA_DIR = Path(__file__).parent / "data"
_REGISTRY_PATH = _DATA_DIR / "expression_sources.yaml"


@dataclass(frozen=True)
class ExpressionSource:
    """One upstream expression dataset feeding one or more cancer cohorts."""

    id: str
    category: str
    cancer_codes: tuple[str, ...]
    source_type: str
    builder: str | None = None
    builder_args: tuple[str, ...] = ()
    external_build_exemption: str | None = None
    routed_samples_may_overlap: bool = False
    project_id: str | None = None
    project_ids: tuple[str, ...] = ()
    accession: str | None = None
    url: str | None = None
    unit: str | None = None
    expected_size_gb: float | None = None
    citation: str | None = None
    special_handling: str | None = None
    recount3_srp: str | None = None
    source_cohort: str | None = None
    library_prep: str | None = None
    source_project: str | None = None
    source_version: str | None = None
    tumor_origin: str | None = None
    metastasis_site: str | None = None
    processing_pipeline: str | None = None
    notes: str | None = None
    access_level: str | None = None
    acquisition_status: str | None = None
    data_access_url: str | None = None
    molecular_annotation_status: str | None = None
    source_pmid: str | None = None
    source_scale_class: str | None = None
    linear_tpm_comparable: bool | None = None
    tpm_proxy: bool | None = None
    taxon_id: int = 9606
    context: str = "tumor"
    contract: ExpressionContract | None = None

    def __post_init__(self):
        if type(self.taxon_id) is not int or self.taxon_id < 1:
            raise ValueError("Positive integer taxonomy identifier required")
        if self.context not in EXPRESSION_CONTEXTS:
            raise ValueError("Unknown expression context")
        if self.contract is None:
            if self.taxon_id != 9606:
                raise ValueError("Nonhuman expression requires an explicit contract")
        elif not isinstance(self.contract, ExpressionContract) or (
            self.contract.taxon_id != self.taxon_id
            or self.contract.source_id != self.id
            or self.contract.context != self.context
            or self.contract.source_version != self.source_version
            or self.contract.measurement.unit != self.unit
            or self.contract.measurement.is_proxy != self.tpm_proxy
            or self.contract.measurement.scale != self.source_scale_class
        ):
            raise ValueError("Expression source and contract disagree")

    @property
    def expression_contract(self) -> ExpressionContract:
        """Explicit metadata for new sources; honest unknowns for historical human data."""
        if self.contract is not None:
            return self.contract
        return ExpressionContract(
            taxon_id=9606,
            source_id=self.id,
            source_version=self.source_version,
            context=self.context,
            reference=None,
            reference_status="legacy_human",
            feature_namespace=None,
            measurement=ExpressionMeasurement(
                "unknown",
                self.unit,
                self.source_scale_class,
                self.tpm_proxy,
            ),
            content_provenance={},
            provenance_notes="Legacy human registry; exact quantification reference and feature level were not declared by this schema.",
        )


def _coerce_tuple(value) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    return tuple(str(v) for v in value)


def _coerce_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_bool(value, *, field: str) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    raise ValueError(f"{field} must be a YAML boolean")


def _coerce_optional_bool(value, *, field: str) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    raise ValueError(f"{field} must be a YAML boolean")


def _clean(value):
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def expression_source_registry_path(*, species="human") -> Path:
    """Filesystem path to the bundled ``expression_sources.yaml`` registry.

    This is the stable public path helper for downstream build tooling that
    needs to pass oncoref's source registry to a subprocess or watcher. Prefer
    :func:`expression_source_registry_entries` for ordinary Python access.
    """
    selected = resolve_species(species)
    if selected.taxon_id == 9606:
        return _REGISTRY_PATH
    return _DATA_DIR / "species" / selected.name / "expression_sources.yaml"


def expression_source_registry_text(*, species="human") -> str:
    """Raw text of the bundled ``expression_sources.yaml`` registry."""
    path = expression_source_registry_path(species=species)
    return path.read_text() if path.exists() else "sources: []\n"


def expression_source_registry_entries(
    source_type: str | Iterable[str] | None = None,
    *,
    species="human",
) -> tuple[dict, ...]:
    """Raw source-registry entries from ``expression_sources.yaml``.

    The typed :func:`expression_sources` view intentionally exposes the stable
    inspection fields. This helper preserves the full YAML dictionaries for
    build tools that need source-specific fields such as ``file_name``,
    ``sample_filter``, or ``sample_to_cancer_code``. ``source_type`` optionally
    filters to one or more source families, e.g. ``"geo-matrix"``.
    """
    import yaml

    doc = yaml.safe_load(expression_source_registry_text(species=species)) or {}
    entries = tuple(dict(entry) for entry in doc.get("sources", []))
    if source_type is None:
        return entries
    wanted = {source_type} if isinstance(source_type, str) else set(source_type)
    return tuple(entry for entry in entries if entry.get("source_type") in wanted)


@cache
def load_registry(species="human") -> tuple[ExpressionSource, ...]:
    """Parse the bundled ``expression_sources.yaml`` into ``ExpressionSource``s."""
    out = []
    selected = resolve_species(species)
    for entry in expression_source_registry_entries(species=selected):
        taxon_id = entry.get("taxon_id", 9606)
        if taxon_id != selected.taxon_id:
            raise ValueError(f"Source {entry['id']} is in the wrong species registry")
        out.append(
            ExpressionSource(
                taxon_id=taxon_id,
                context=entry.get("context", "tumor"),
                contract=ExpressionContract.from_dict(entry["contract"])
                if entry.get("contract") is not None
                else None,
                id=str(entry["id"]),
                category=str(entry.get("category", "expression")),
                cancer_codes=_coerce_tuple(entry.get("cancer_codes")),
                source_type=str(entry.get("source_type", "")),
                builder=_clean(entry.get("builder")),
                builder_args=_coerce_tuple(entry.get("builder_args")),
                external_build_exemption=_clean(entry.get("external_build_exemption")),
                routed_samples_may_overlap=_coerce_bool(
                    entry.get("routed_samples_may_overlap"),
                    field=f"{entry['id']}.routed_samples_may_overlap",
                ),
                project_id=_clean(entry.get("project_id")),
                project_ids=_coerce_tuple(entry.get("project_ids")),
                accession=_clean(entry.get("accession")),
                url=_clean(entry.get("url")),
                unit=_clean(entry.get("unit")),
                source_scale_class=_clean(entry.get("source_scale_class")),
                linear_tpm_comparable=_coerce_optional_bool(
                    entry.get("linear_tpm_comparable"),
                    field=f"{entry['id']}.linear_tpm_comparable",
                ),
                tpm_proxy=_coerce_optional_bool(
                    entry.get("tpm_proxy"),
                    field=f"{entry['id']}.tpm_proxy",
                ),
                expected_size_gb=_coerce_float(entry.get("expected_size_gb")),
                citation=_clean(entry.get("citation")),
                source_pmid=_clean(entry.get("source_pmid")),
                special_handling=_clean(entry.get("special_handling")),
                recount3_srp=_clean(entry.get("recount3_srp")),
                source_cohort=_clean(entry.get("source_cohort")),
                library_prep=_clean(entry.get("library_prep")),
                source_project=_clean(entry.get("source_project")),
                source_version=_clean(entry.get("source_version")),
                tumor_origin=_clean(entry.get("tumor_origin")),
                metastasis_site=_clean(entry.get("metastasis_site")),
                processing_pipeline=_clean(entry.get("processing_pipeline")),
                notes=_clean(entry.get("notes")),
                access_level=_clean(entry.get("access_level")),
                acquisition_status=_clean(entry.get("acquisition_status")),
                data_access_url=_clean(entry.get("data_access_url")),
                molecular_annotation_status=_clean(entry.get("molecular_annotation_status")),
            )
        )
    if len({s.id for s in out}) != len(out):
        raise ValueError("Duplicate expression source ID within species")
    return tuple(out)


def expression_sources(*, species="human", context=None) -> tuple[ExpressionSource, ...]:
    """Every registered expression source."""
    if context is not None and context not in EXPRESSION_CONTEXTS:
        raise ValueError(f"Unknown expression context {context!r}")
    return tuple(
        s
        for s in load_registry(resolve_species(species))
        if context is None or s.context == context
    )


def expression_source(source_id: str, *, species="human") -> ExpressionSource | None:
    """Look up one source by id."""
    for s in load_registry(resolve_species(species)):
        if s.id == source_id:
            return s
    return None


def sources_for_cancer_code(code: str, *, species="human") -> list[ExpressionSource]:
    """Every source that feeds a given cancer code."""
    return [s for s in expression_sources(species=species) if code in s.cancer_codes]


def expression_sources_df(*, species="human") -> pd.DataFrame:
    """The registry as a tabular view (one row per source) for inspection."""
    return pd.DataFrame(
        {
            "id": s.id,
            "taxon_id": s.taxon_id,
            "context": s.context,
            "reference_status": s.expression_contract.reference_status,
            "contract_id": s.expression_contract.contract_id,
            "source_type": s.source_type,
            "cancer_codes": ";".join(s.cancer_codes),
            "n_codes": len(s.cancer_codes),
            "source_cohort": s.source_cohort,
            "source_project": s.source_project,
            "source_version": s.source_version,
            "unit": s.unit,
            "source_scale_class": s.source_scale_class,
            "linear_tpm_comparable": s.linear_tpm_comparable,
            "tpm_proxy": s.tpm_proxy,
            "library_prep": s.library_prep,
            "project_id": ";".join(s.project_ids) or s.project_id or s.accession or s.recount3_srp,
            "tumor_origin": s.tumor_origin,
            "metastasis_site": s.metastasis_site,
            "processing_pipeline": s.processing_pipeline,
            "notes": s.notes,
            "access_level": s.access_level,
            "acquisition_status": s.acquisition_status,
            "data_access_url": s.data_access_url,
            "molecular_annotation_status": s.molecular_annotation_status,
            "expected_size_gb": s.expected_size_gb,
            "citation": s.citation,
            "source_pmid": s.source_pmid,
            "builder": s.builder,
            "builder_args": " ".join(s.builder_args),
            "external_build_exemption": s.external_build_exemption,
            "routed_samples_may_overlap": s.routed_samples_may_overlap,
        }
        for s in expression_sources(species=species)
    )


def expression_source_candidates(cancer_code: str | None = None) -> pd.DataFrame:
    """Candidate expression sources per cancer type (``cancer_code``,
    ``source_status``, ``reference_code``, ``source_project``, ``accession``,
    ``assay``, ``processing_plan``, …) — the build-planning view of which cohort
    could back a type's reference. Filtered to one code when given. Defensive copy."""
    from .load_dataset import get_data

    df = get_data("cancer-expression-source-candidates").copy()
    if cancer_code is not None:
        from .cancer_types import resolve_cancer_type

        df = df[df["cancer_code"].astype(str) == resolve_cancer_type(cancer_code)].reset_index(
            drop=True
        )
    return df
