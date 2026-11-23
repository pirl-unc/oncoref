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

"""Per-sample curation manifest for the cohort expression references.

``cancer-reference-expression-samples.csv.gz`` records, for every individual
sample that went into a cohort summary, its provenance and the curation decision
(which cancer code it was assigned, which cohort/project it came from, the source
file, the raw unit + processing pipeline, and whether it was included or excluded
— with the exclusion reason). This is the sample-level curation needed to
reproduce / audit the per-cohort summaries.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache

import pandas as pd

from .cancer_types import resolve_cancer_type
from .load_dataset import get_data

PATIENT_IDENTITY_POLICY = "oncoref.patient-identity.2026-11-10.1"


def _patient_key(source: str, identifier: str) -> str:
    """Use a patient domain shared by aliquots and derived cohort views."""
    for pattern, prefix in [
        (r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})(?:-|$)", "TCGA"),
        (r"^(TARGET-\d{2}-[A-Z0-9]{6})(?:-|$)", "TARGET"),
        (r"^MMRF[_-](\d+)(?:[_-]|$)", "MMRF"),
    ]:
        match = re.match(pattern, identifier)
        if match:
            return prefix + ":" + match[1]
    domain = source
    if source.startswith("TREEHOUSE_"):
        domain = "TREEHOUSE"
    elif source.startswith("OPENPBTA_"):
        domain = "OPENPBTA"
    return domain + ":" + identifier


def sample_identities(rows: pd.DataFrame) -> pd.DataFrame:
    """Join sample rows to documented patient identities without inventing donors.

    Requires ``source_cohort`` and ``sample_id``; preserves row order and all
    other input columns. Adds a physical ``namespace``/``sample_key``, globally
    namespaced ``participant_id``, status, basis, references and candidate IDs.
    An empty participant ID means unresolved or conflicting evidence, never a
    distinct patient. Count resolved participant IDs, not matrix columns, when
    testing independent-patient support.

    The packaged case manifest is not universally a patient manifest: several
    GEO importers copy sample labels into case_id. Only reviewed case mappings
    (BeatAML, CGCI and the curated nonredundant CLL-map cohort), explicit molecular
    donor mappings, documented project barcodes and Treehouse donor families
    establish identity here. A conflict between supported mappings stays visible.
    This function annotates identity; it does not filter samples or apply QC.
    """
    from .source_matrices import source_sample_namespace

    required = {"source_cohort", "sample_id"}
    if not required <= set(rows):
        raise ValueError("sample_identities requires source_cohort and sample_id")
    if rows[list(required)].isna().any().any() or any(
        rows[c].astype(str).str.strip().eq("").any() for c in required
    ):
        raise ValueError("Sample and source identifiers must be nonempty")

    mappings = {}
    case_labels = set()

    def clean(value):
        return "" if pd.isna(value) else str(value).strip()

    def add(source, sample, donor, basis, reference):
        source, sample, donor = clean(source), clean(sample), clean(donor)
        if source and sample and donor:
            key = (source_sample_namespace(source), sample)
            mappings.setdefault(key, set()).add(
                (_patient_key(source, donor), basis, clean(reference))
            )

    cases = sample_manifest()
    for row in cases.to_dict("records"):
        source, sample = clean(row["source_cohort"]), clean(row["sample_id"])
        if clean(row.get("case_id")):
            case_labels.add((source_sample_namespace(source), sample))
        if clean(row.get("included")).lower() not in {"true", "1"}:
            continue
        if source in {"BEATAML_OHSU_2022", "CGCI_BLGSP"}:
            add(source, sample, row.get("case_id"), "curated_case_manifest", row.get("source_url"))
        elif source == "CLLMAP_2022" and sample not in {
            "CRC-0007",
            "CRC-0011",
            "CRC-0028",
            "CRC-0033",
            "DFCI-5053",
            "JB-0010",
            "GCLL-0136",  # suspected MCL, separate from the six duplicate exclusions
        }:
            add(
                source,
                sample,
                row.get("case_id"),
                "curated_nonredundant_cohort",
                "https://www.cllmap.org/downloads.html",
            )
    for row in sample_molecular_provenance().to_dict("records"):
        for column in ("sample_id", "library_id"):
            add(
                row["source_cohort"],
                row.get(column),
                row.get("donor_id"),
                "curated_donor_manifest",
                row.get("evidence_source"),
            )

    output = rows.copy()
    annotations = []
    for row in rows.itertuples(index=False):
        source, sample = str(row.source_cohort), str(row.sample_id)
        namespace = source_sample_namespace(source)
        key = namespace, sample
        candidates = set(mappings.get(key, set()))
        # An empty fallback domain prevents a source literally named TCGA/TARGET/
        # MMRF from making an arbitrary sample label look like a project barcode.
        barcode = _patient_key("", sample)
        if barcode.startswith(("TCGA:", "TARGET:", "MMRF:")):
            candidates.add(
                (
                    barcode,
                    "project_patient_barcode",
                    "https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/"
                    if barcode.startswith("TCGA:")
                    else "https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=case",
                )
            )
        # Prefer curated donor IDs; the family fallback is the existing Treehouse
        # builder's rule for releases lacking study_donor_id, not a generic regex.
        if not candidates and source.startswith("TREEHOUSE_"):
            match = re.fullmatch(r"(THR?\d+_[^_]+)_S\d+", sample)
            if match:
                candidates.add(
                    (
                        _patient_key(source, match[1]),
                        "treehouse_donor_family",
                        "https://treehousegenomics.ucsc.edu/explore-our-data/",
                    )
                )
        ids = sorted({c[0] for c in candidates})
        status = "resolved" if len(ids) == 1 else "conflict" if ids else "unresolved"
        reason = (
            "documented_patient_mapping"
            if len(ids) == 1
            else "conflicting_patient_mappings"
            if ids
            else "case_field_not_verified_as_patient"
            if key in case_labels
            else "no_documented_patient_mapping"
        )
        annotations.append(
            {
                "namespace": namespace,
                "sample_key": namespace + ":" + sample,
                "participant_id": ids[0] if len(ids) == 1 else "",
                "patient_identity_status": status,
                "patient_identity_reason": reason,
                "patient_identity_basis": ";".join(sorted({c[1] for c in candidates})),
                "patient_identity_references": json.dumps(
                    sorted({c[2] for c in candidates if c[2]})
                ),
                "patient_identity_candidates": json.dumps(ids),
                "patient_identity_policy": PATIENT_IDENTITY_POLICY,
            }
        )
    columns = [
        "namespace",
        "sample_key",
        "participant_id",
        "patient_identity_status",
        "patient_identity_reason",
        "patient_identity_basis",
        "patient_identity_references",
        "patient_identity_candidates",
        "patient_identity_policy",
    ]
    metadata = pd.DataFrame(annotations, index=rows.index, columns=columns)
    for column in columns:
        output[column] = metadata[column]
    return output


@lru_cache(maxsize=1)
def _manifest() -> pd.DataFrame:
    return get_data("cancer-reference-expression-samples", copy=False)


def sample_manifest() -> pd.DataFrame:
    """The full per-sample curation manifest (one row per sample). Defensive copy."""
    return _manifest().copy()


def samples_for_cancer_code(code: str, *, included_only: bool = True) -> pd.DataFrame:
    """Manifest rows assigned to a cancer code (included samples by default).

    The code is alias-resolved, so ``"prostate"`` / ``"lung_adeno"`` / a pre-rename
    code all work (an unknown code passes through and simply matches nothing)."""
    resolved = resolve_cancer_type(code, strict=False) or code
    df = _manifest()
    sub = df[df["cancer_code"].astype(str) == resolved]
    if included_only and "included" in sub.columns:
        sub = sub[sub["included"].astype(str).str.lower().isin(("true", "1"))]
    return sub.copy()


def samples_for_cohort(source_cohort: str, *, included_only: bool = True) -> pd.DataFrame:
    """Manifest rows from a source cohort (included samples by default)."""
    df = _manifest()
    sub = df[df["source_cohort"].astype(str) == source_cohort]
    if included_only and "included" in sub.columns:
        sub = sub[sub["included"].astype(str).str.lower().isin(("true", "1"))]
    return sub.copy()


def sample_counts_by_cancer_code(*, included_only: bool = True) -> pd.Series:
    """Number of curated samples per cancer code (a quick coverage view)."""
    df = _manifest()
    if included_only and "included" in df.columns:
        df = df[df["included"].astype(str).str.lower().isin(("true", "1"))]
    return df["cancer_code"].astype(str).value_counts()


@lru_cache(maxsize=1)
def _molecular_provenance() -> pd.DataFrame:
    return get_data("cancer-reference-sample-molecular-provenance", copy=False)


def sample_molecular_provenance() -> pd.DataFrame:
    """Public diagnosis and molecular evidence for reference/acquisition samples.

    Rows are libraries. ``donor_id`` makes technical or longitudinal libraries
    auditable without treating them as independent donors. A diagnosis-only row
    has ``molecular_status="unknown"`` and ``orthogonal_confirmed=False``;
    diagnosis never supplies a fusion call.
    """
    return _molecular_provenance().copy()


def molecular_provenance_for_cancer_code(code: str) -> pd.DataFrame:
    """Sample-level molecular provenance for one alias-resolved cancer code."""
    resolved = resolve_cancer_type(code, strict=False) or code
    rows = _molecular_provenance()
    return rows.loc[rows["cancer_code"].astype(str).eq(resolved)].reset_index(drop=True).copy()


def molecular_provenance_for_sample(sample_id: str) -> pd.DataFrame:
    """All provenance rows matching a sample, library, or donor identifier."""
    query = str(sample_id)
    rows = _molecular_provenance()
    mask = (
        rows["sample_id"].astype(str).eq(query)
        | rows["library_id"].astype(str).eq(query)
        | rows["donor_id"].astype(str).eq(query)
    )
    return rows.loc[mask].reset_index(drop=True).copy()


def molecular_sample_counts(code: str | None = None) -> pd.DataFrame:
    """Library and distinct-donor counts by physical source cohort.

    The optional cancer code is alias-resolved. Counts include public expression
    and controlled acquisition rows; ``n_expression_libraries`` states how many
    of the listed libraries currently have loadable expression.
    """
    rows = _molecular_provenance() if code is None else molecular_provenance_for_cancer_code(code)
    group_cols = ["cancer_code", "source_id", "source_cohort", "access_level"]
    if rows.empty:
        return pd.DataFrame(
            columns=[*group_cols, "n_libraries", "n_donors", "n_expression_libraries"]
        )
    counts = (
        rows.groupby(group_cols, dropna=False, sort=True)
        .agg(
            n_libraries=("library_id", "nunique"),
            n_donors=("donor_id", "nunique"),
            n_expression_libraries=("expression_available", "sum"),
        )
        .reset_index()
    )
    return counts
