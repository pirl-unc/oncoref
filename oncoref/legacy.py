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

"""Machine-readable disposition of frozen legacy compatibility datasets."""

from __future__ import annotations

import threading
import warnings

import pandas as pd

from .load_dataset import get_data

LEGACY_DATASET_STATUSES = frozenset(
    {
        "frozen_replacement_available",
        "frozen_migration_pending",
        "frozen_downstream_owned",
    }
)

LEGACY_COMPATIBILITY_POLICIES = frozenset({"ship_read_only_no_new_rows"})

LEGACY_TYPED_ACCESSORS = {
    "cancer-driver-genes": "oncoref.cancer_genes.cancer_driver_genes_df",
    "cancer-driver-variants": "oncoref.cancer_genes.cancer_driver_variants_df",
    "cancer-key-genes": "oncoref.cancer_genes.cancer_key_genes_df",
    "cancer-type-genes": "oncoref.cancer_genes.cancer_type_genes_df",
    "cancer-viral-antigens": "oncoref.cancer_genes.cancer_viral_antigens_df",
    "disease-state-rules": "oncoref.cancer_genes.disease_state_rules_df",
    "fusion-expression-effects": "oncoref.fusions.fusion_expression_effect_rules_df",
    "fusion-surrogate-expression": "oncoref.fusions.fusion_surrogate_expression_df",
    "narrative-gene-sets": "oncoref.cancer_genes.narrative_gene_sets_df",
    "cancer-response-signatures": "oncoref.response_signatures.response_signatures_df",
    "rare-cancer-fusion-rules": "oncoref.fusions.rare_cancer_fusion_rules_df",
}

LEGACY_CTA_TOP_LEVEL_EXPORTS: tuple[str, ...] = (
    "CTA_by_axes",
    "CTA_excluded_gene_ids",
    "CTA_excluded_gene_names",
    "CTA_filtered_gene_ids",
    "CTA_filtered_gene_names",
    "CTA_gene_ids",
    "CTA_gene_names",
    "CTA_never_expressed_gene_ids",
    "CTA_never_expressed_gene_names",
    "CTA_placental_restricted_gene_ids",
    "CTA_placental_restricted_gene_names",
    "CTA_relaxed_reproductive_gene_ids",
    "CTA_relaxed_reproductive_gene_names",
    "CTA_testis_restricted_gene_ids",
    "CTA_testis_restricted_gene_names",
    "CTA_unfiltered_gene_ids",
    "CTA_unfiltered_gene_names",
    "cta_by_axes",
    "cta_candidate_references",
    "cta_clinical_target_evidence",
    "cta_clinical_target_gene_ids",
    "cta_clinical_target_gene_names",
    "cta_clinical_target_references",
    "cta_df",
    "cta_evidence",
    "cta_excluded_clinical_target_gene_ids",
    "cta_excluded_clinical_target_gene_names",
    "cta_excluded_gene_ids",
    "cta_excluded_gene_names",
    "cta_extended_df",
    "cta_extended_gene_ids",
    "cta_extended_gene_names",
    "cta_filtered_gene_ids",
    "cta_filtered_gene_names",
    "cta_gene_id_to_name",
    "cta_gene_ids",
    "cta_gene_names",
    "cta_never_expressed_gene_ids",
    "cta_never_expressed_gene_names",
    "cta_placental_restricted_gene_ids",
    "cta_placental_restricted_gene_names",
    "cta_relaxed_reproductive_gene_ids",
    "cta_relaxed_reproductive_gene_names",
    "cta_specificity_audit",
    "cta_specificity_audit_references",
    "cta_symbol_for_alias",
    "cta_testis_restricted_gene_ids",
    "cta_testis_restricted_gene_names",
    "cta_unfiltered_gene_ids",
    "cta_unfiltered_gene_names",
    "cta_warning_gene_ids",
    "cta_warning_gene_names",
    "cta_warning_references",
)

LEGACY_API_MIGRATIONS: tuple[dict[str, str], ...] = (
    {
        "legacy_surface": "oncoref.cta / oncoref.cta_* / oncoref.CTA_*",
        "generation": "v1_hpa_nominated_cta",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref.antigen_sets.gene_ids('CTA') / oncoref.target_selection.results()",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": (
            "Historical HPA v23 + literature CTA list (~574 expressed genes, including genes now "
            "classified Unsafe under rare-cell and literature exclusions); import oncoref.cta "
            "explicitly only for historical reproduction."
        ),
    },
    {
        "legacy_surface": "oncoref.target_sets",
        "generation": "v1_20_target_sets_213",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref.antigen_sets.select('CTA') / oncoref.target_selection.results()",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": "Frozen 213-gene expression snapshot superseded by set-based CA/RA/CTA catalog.",
    },
    {
        "legacy_surface": "oncoref.cta_screen",
        "generation": "v1_22_genome_wide_screen",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref.target_selection.results()",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": "Early genome-wide HPA screen superseded by multi-atlas + single-cell + literature catalog.",
    },
    {
        "legacy_surface": "oncoref.normal_expression_categories",
        "generation": "v1_21_normal_expression_categories",
        "current_status": "frozen_replacement_available",
        "replacement_surface": (
            "oncoref.normal_expression_categories.reconcile_with_target_selection() / "
            "oncoref.target_selection.results()"
        ),
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": (
            "HPA v23 Normal/Reproductive/Other(Disposable,Silent)/Unassessed screen reconciled "
            "against active Unsafe/CA/RA/CTA/Other subtypes without weakening Unsafe precedence."
        ),
    },
    {
        "legacy_surface": "oncoref.cta_membership",
        "generation": "v1_26_biological_membership_1416",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref.antigen_sets.proteins('CA|RA') / oncoref.target_selection.results()",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": "Literature-inclusive 1,416-gene biological nomination snapshot.",
    },
    {
        "legacy_surface": "oncoref.cta_quantitative / oncoref.cta_research_policy",
        "generation": "v1_31_quantitative_cta",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref.target_selection.results()",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": "Intermediate quantitative screen superseded by rebuilt target_selection catalog.",
    },
    {
        "legacy_surface": "oncoref.cta_selection",
        "generation": "v1_37_cta_selection_catalog",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref.target_selection.results() / oncoref.antigen_sets",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": "Versioned pre-CA/RA unified catalog snapshot retained for release-to-release audits.",
    },
    {
        "legacy_surface": "oncoref.peptides / oncoref.cta_peptides",
        "generation": "v1_cta_9mer_counts",
        "current_status": "active_facade_over_current_catalog",
        "replacement_surface": (
            "oncoref.peptide_subtraction.PeptideIndex / oncoref.antigen_sets.select('CTA').background_index()"
        ),
        "active_policy_id": "oncoref.target-selection-sequence-roles.2026-11-21.1",
        "notes": "Precomputed 9-mer summary helpers use active CTA defaults; use PeptideIndex for custom k-mers/edit distance.",
    },
    {
        "legacy_surface": "oncoref.cta (CLI)",
        "generation": "v1_cli_cta",
        "current_status": "frozen_replacement_available",
        "replacement_surface": "oncoref antigens list --set CTA",
        "active_policy_id": "oncoref.target-selection.2026-11-20.2",
        "notes": "Legacy CLI prints historical snapshot warning to stderr; use 'oncoref antigens'.",
    },
)

_WARNED_LEGACY_DATASETS: set[str] = set()
_WARNED_LEGACY_APIS: set[str] = set()
_WARNING_LOCK = threading.Lock()


def legacy_dataset_dispositions() -> pd.DataFrame:
    """Return the reviewed migration/disposition row for every legacy table.

    Legacy tables remain importable for compatibility, but are frozen. The
    returned frame names the owner and replacement surface for new work, so
    callers do not need to infer the stack boundary from comments in the data
    manifest.
    """
    return get_data("legacy-dataset-dispositions").copy()


def legacy_dataset_disposition(dataset: str) -> dict | None:
    """Disposition record for one legacy dataset, or ``None`` if it is not legacy."""
    name = str(dataset).strip()
    rows = legacy_dataset_dispositions()
    hit = rows.loc[rows["dataset"].astype(str).eq(name)]
    return None if hit.empty else hit.iloc[0].to_dict()


def warn_legacy_dataset_access(dataset: str, *, stacklevel: int = 2) -> None:
    """Warn once per process when a typed API reads a frozen snapshot.

    The generic :func:`oncoref.load_dataset.get_data` loader intentionally stays
    warning-free as a low-level compatibility escape hatch. Typed accessors call
    this helper so their owner, replacement surface, and freeze policy are not
    mistaken for current curation.
    """
    name = str(dataset).strip()
    disposition = legacy_dataset_disposition(name)
    if disposition is None:
        raise ValueError(f"{name!r} is not a registered legacy dataset")
    with _WARNING_LOCK:
        if name in _WARNED_LEGACY_DATASETS:
            return
        warnings.warn(
            f"{name!r} is a frozen compatibility snapshot "
            f"({disposition['current_status']}); new work belongs to "
            f"{disposition['replacement_owner']} via "
            f"{disposition['replacement_surface']}. "
            f"Policy: {disposition['compatibility_policy']}.",
            DeprecationWarning,
            stacklevel=stacklevel + 1,
        )
        _WARNED_LEGACY_DATASETS.add(name)


def legacy_api_migrations() -> pd.DataFrame:
    """Return the structured migration table across historical CTA/target API generations."""
    from . import sequence_bundle, target_selection

    rows = []
    for entry in LEGACY_API_MIGRATIONS:
        row = dict(entry)
        if row["active_policy_id"].startswith("oncoref.target-selection."):
            row["active_policy_id"] = target_selection.POLICY_ID
        elif row["active_policy_id"].startswith("oncoref.target-selection-sequence-roles."):
            row["active_policy_id"] = sequence_bundle.ROLE_POLICY_ID
        rows.append(row)
    return pd.DataFrame(rows)


def legacy_api_migration(surface: str) -> dict[str, str] | None:
    """Lookup the migration record for a legacy module, top-level accessor, or CLI surface."""
    token = str(surface).strip()
    if (
        token in LEGACY_CTA_TOP_LEVEL_EXPORTS
        or token.removeprefix("oncoref.") in LEGACY_CTA_TOP_LEVEL_EXPORTS
    ):
        return legacy_api_migrations().iloc[0].to_dict()
    frame = legacy_api_migrations()
    for row in frame.to_dict(orient="records"):
        parts = [p.strip() for p in row["legacy_surface"].split("/")]
        if token in parts or f"oncoref.{token}" in parts:
            return row
    return None


def warn_legacy_api_access(
    name: str,
    *,
    replacement: str = "oncoref.antigen_sets.gene_ids('CTA') / oncoref.target_selection.results()",
    stacklevel: int = 2,
) -> None:
    """Emit a once-per-process ``DeprecationWarning`` for top-level legacy CTA accessors."""
    from . import target_selection

    qualified = name if name.startswith("oncoref.") else f"oncoref.{name}"
    with _WARNING_LOCK:
        if qualified in _WARNED_LEGACY_APIS:
            return
        warnings.warn(
            f"{qualified!r} is a frozen compatibility snapshot "
            f"(frozen_replacement_available); new work belongs to oncoref via "
            f"{replacement} (active policy {target_selection.POLICY_ID}). "
            "Import 'oncoref.cta' explicitly only when reproducing historical v1 CTA counts.",
            DeprecationWarning,
            stacklevel=stacklevel + 1,
        )
        _WARNED_LEGACY_APIS.add(qualified)
