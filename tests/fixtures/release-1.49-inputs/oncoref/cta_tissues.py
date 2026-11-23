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

"""Tissue definitions and HPA confidence thresholds for CTA restriction analysis.

Ported (HPA-only) from ``tsarina.tissues`` so oncoref can regenerate the
bundled ``cancer-testis-antigens.csv`` restriction/filter columns from HPA data
alone. The MS-related constants and the CTA-universe membership rules (the
gene-family exclusion and the never-expressed rescue) are intentionally NOT
here — they live in :mod:`oncoref.cta`.

These HPA label sets are projections of the shared, cited normal-tissue taxonomy
in :mod:`oncoref.normal_tissues`. Its species-aware scopes also cover developmental
contexts and conditionally disposable normal compartments.

- **Core**: testis, ovary, placenta -- the classic CT restriction definition.
- **Extended**: all reproductive tissues represented in HPA, including mammary
  tissue. Bulk breast does not establish glandular-cell specificity.
- **Permissive**: compatibility alias of extended.

Thymus is somatic. The explicitly named ``*_and_thymus_frac`` columns remain
available as separate historical summaries. Tissue scopes do not establish
immune privilege or target safety.
"""

from __future__ import annotations

from .normal_tissues import normal_tissue_mappings, normal_tissue_scope

# ── Tissue sets ─────────────────────────────────────────────────────────────


def _hpa_tissues(scope: str, source: str) -> frozenset[str]:
    members = normal_tissue_scope(scope, species="human")
    return frozenset(
        row["source_tissue"]
        for row in normal_tissue_mappings(species="human", source_id=source, source_version="v23")
        if not row["source_cell_type"] and row["concept_id"] in members
    )


CORE_REPRODUCTIVE_TISSUES = _hpa_tissues("reproductive_core", "hpa_rna_consensus")
EXTENDED_REPRODUCTIVE_TISSUES = _hpa_tissues("reproductive", "hpa_rna_consensus")
PERMISSIVE_REPRODUCTIVE_TISSUES = EXTENDED_REPRODUCTIVE_TISSUES


def reproductive_tissues(scope: str = "core") -> frozenset[str]:
    """Tissues counted in the named CTA panel's RNA numerator."""
    if scope == "core":
        return CORE_REPRODUCTIVE_TISSUES
    if scope == "extended":
        return EXTENDED_REPRODUCTIVE_TISSUES
    raise ValueError(f"Unknown CTA tissue scope {scope!r}; expected 'core' or 'extended'")


def cta_dataset_name(scope: str = "core") -> str:
    """Bundled evidence snapshot for a named tissue scope."""
    reproductive_tissues(scope)  # Validate before choosing a dataset.
    return "cancer-testis-antigens" + ("-extended" if scope == "extended" else "")


#: Reproductive RNA tissues; thymus and disposable compartments remain somatic.
NON_SOMATIC_TISSUES: frozenset[str] = EXTENDED_REPRODUCTIVE_TISSUES

#: All tissues considered reproductive for the protein (IHC) restriction call.
#: HPA uses "endometrium 1"/"endometrium 2" in the normal_tissue table, so those
#: split labels are included alongside the consensus "endometrium".
ALL_REPRODUCTIVE_TISSUES = _hpa_tissues("reproductive", "hpa_normal_tissue")


# ── Adaptive HPA protein → RNA thresholds ──────────────────────────────────

#: Mapping from HPA antibody reliability tier to the minimum deflated RNA
#: reproductive fraction required for a gene to pass the adaptive confidence
#: filter. Higher-confidence protein data allows a more relaxed RNA threshold.
HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS: dict[str, float] = {
    "Enhanced": 0.80,
    "Supported": 0.90,
    "Approved": 0.95,
    "Uncertain": 0.97,
    "Missing": 0.97,
}

#: RNA expression floor (nTPM) for the ``never_expressed`` flag. A gene with no
#: HPA protein (IHC) data and a maximum RNA nTPM below this value across all
#: tissues is flagged ``never_expressed`` -- it passes the reproductive-
#: restriction filter only because the ``+1`` deflation pseudocount yields a 1.0
#: fraction when every tissue is below 1 nTPM, but HPA lacks the signal to
#: confirm restriction.
HPA_EXPRESSION_FLOOR_NTPM: float = 2.0


# ── Protein reliability ordering ────────────────────────────────────────────

#: HPA antibody reliability tiers, ordered from strongest to weakest.
PROTEIN_RELIABILITY_ORDER: list[str] = [
    "Enhanced",
    "Supported",
    "Approved",
    "Uncertain",
]

#: HPA IHC levels that count as "protein detected" in a tissue.
PROTEIN_DETECTED_LEVELS: frozenset[str] = frozenset({"Low", "Medium", "High"})


# ── Safety-critical tissue groups ──────────────────────────────────────────

#: Safety-critical tissue groups with nTPM threshold. Genes with max nTPM >=
#: threshold in any tissue in the group get flagged.
#:
#: The ``brain`` group spans two HPA vocabularies on purpose. Ten names are
#: present in the bulk ``rna_tissue_consensus`` and so drive the RNA safety flag
#: (``rna_brain_max_ntpm``): amygdala, basal ganglia, cerebellum, cerebral
#: cortex, choroid plexus, hippocampal formation, hypothalamus, midbrain,
#: retina, spinal cord. Four names (medulla oblongata, pons, thalamus, white
#: matter) come from HPA's finer brain-specific dataset; they are absent from
#: the consensus -- inert for the RNA flag -- but are kept for vocabulary parity
#: with tsarina. These are conceptual/RNA safety labels, not an IHC selection
#: vocabulary; use ``hpa.resolve_safety_tissue_group`` before matching an HPA
#: normal-tissue source.
SAFETY_TISSUE_GROUPS: dict[str, set[str]] = {
    "brain": {
        "amygdala",
        "basal ganglia",
        "cerebellum",
        "cerebral cortex",
        "choroid plexus",
        "hippocampal formation",
        "hypothalamus",
        "medulla oblongata",
        "midbrain",
        "pons",
        "retina",
        "spinal cord",
        "thalamus",
        "white matter",
    },
    "heart": {"heart muscle"},
    "lung": {"lung"},
    "liver": {"liver"},
    "pancreas": {"pancreas"},
}

#: Default nTPM threshold for safety flags.
SAFETY_NTPM_THRESHOLD: float = 5.0


def adaptive_rna_threshold(protein_reliability: str) -> float:
    """Return the minimum deflated RNA reproductive fraction for a protein tier.

    Unknown labels (incl. "no data" / "Missing") use the strict "Missing" floor.
    """
    normalized = str(protein_reliability).strip().casefold()
    by_casefold = {k.casefold(): k for k in HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS}
    key = by_casefold.get(normalized, "Missing")
    return HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS[key]
