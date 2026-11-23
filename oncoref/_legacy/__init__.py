"""Quarantined pre-target_selection legacy CTA modules (#756).

Active antigen selection and multi-species target screening live in
:mod:`oncoref.target_selection`, :mod:`oncoref.species_target_selection`,
and :mod:`oncoref.antigen_sets`. Top-level ``oncoref.cta*`` modules remain
thin compatibility shims forwarding to this package.
"""

from __future__ import annotations

QUARANTINED_CTA_MODULES: tuple[str, ...] = (
    "cta",
    "cta_cell_tiers",
    "cta_continuous_screen",
    "cta_contract",
    "cta_curation_plots",
    "cta_gate_audit",
    "cta_landscape",
    "cta_membership",
    "cta_normal_screen",
    "cta_policy",
    "cta_provenance",
    "cta_quantitative",
    "cta_rare_cell",
    "cta_regen",
    "cta_representatives",
    "cta_research_policy",
    "cta_review",
    "cta_screen",
    "cta_selection",
    "cta_selection_audit",
    "cta_sources",
    "cta_tissues",
)

__all__ = ["QUARANTINED_CTA_MODULES"]
