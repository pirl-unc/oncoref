#!/usr/bin/env python3
"""Explain proposed Other subtypes and audit omitted Reproductive evidence.

This is a diagnostic replay, not a released policy or a change to packaged
assignments. The counterfactual changes only the stale bulk Reproductive input;
cell evidence is reported separately and does not automatically admit candidates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cta_profile_data import hpa_profiles  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import cell_rna_policy, normal_cell_screen, reference_data  # noqa: E402
from oncoref import cta_selection_policy as policy  # noqa: E402
from oncoref.cta_funnel import EvidenceLayer  # noqa: E402
from oncoref.expression_diffuseness import HPA_REPRODUCTIVE_BRANCHES  # noqa: E402

CURRENT = ROOT / "docs/audits/cta-selection-20261028"
OUT = ROOT / "docs/audits/cta-other-20261030"
ID = "protein_equivalence_group_id"


def bulk_counterfactual(facts, scores):
    """Use the existing full Reproductive reference, preserving all normal gates."""
    if scores.entity_id.duplicated().any() or set(scores.entity_id) != set(facts[ID]):
        raise ValueError("Require one source-matched profile per protein")
    d = facts.copy()
    d["reproductive_ntpm"] = d[ID].map(scores.set_index("entity_id").reference_value)
    d = policy.concentration_evidence(policy.association_evidence(d), scores)
    d = policy.tumor_restricted_evidence(policy.neither_profiles(d))
    d["hold:Reproductive / embryonic association"] = ~d.positive_context
    d["hold:Reproductive concentration"] = d.reproductive_concentration_hold
    layers = [EvidenceLayer(c[5:], d[c], d["unsafe:" + c[5:]]) for c in d if c.startswith("hold:")]
    result = policy.classify(d, layers, entity_column=ID).classifications
    # Positive association cannot manufacture normal-safety evidence.
    for c in facts.filter(regex="^unsafe:"):
        pd.testing.assert_series_equal(facts[c], result[c])
    assert (
        facts.category.eq("Unsafe").tolist() == result.category.eq("Unsafe").fillna(False).tolist()
    )
    assert result.loc[facts.category.eq("CTA"), "category"].eq("CTA").all()
    return result


def cell_peaks(facts):
    """Retain supported, mapped cell RNA using the existing identity/QC rules."""
    source = reference_data.provenance("hpa_single_cell_tissue", "v23", verify_content=True)
    identity_path = CURRENT / "hpa_single_cell_tissue-accepted-identities.csv.gz"
    grouped_path = ROOT / "docs/audits/cta-harmonized-20261027/group-normal/cluster-contexts.csv.gz"
    inputs = input_hashes([Path(__file__), identity_path, grouped_path])
    receipt_path = OUT / "cell-peaks-receipt.json"
    output = OUT / "supported-cell-peaks.csv.gz"
    request = {
        "inputs": inputs,
        "proteins": sorted(facts[ID]),
        "source": {k: source[k] for k in ["name", "version", "url", "sha256"]},
    }
    if receipt_path.exists() and output.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt.get("request") == request and receipt.get("output_sha256") == sha256(output):
            return pd.read_csv(output, low_memory=False)
    singleton = set(facts.loc[facts.member_count.eq(1), "representative_gene_id"])
    lookup = pd.read_csv(identity_path).set_index("source_gene_id").canonical_gene_id
    rows = []
    for number, chunk in enumerate(pd.read_csv(source["path"], sep="\t", chunksize=500000)):
        chunk["source_gene_id"] = chunk.Gene
        chunk["Gene"] = chunk.Gene.map(lookup)
        d = (
            chunk[chunk.Gene.isin(singleton) & chunk["Read count"].ge(20)]
            .rename(
                columns={
                    "Gene": "gene_id",
                    "Tissue": "tissue",
                    "Cell type": "cell_type",
                    "nTPM": "value",
                }
            )
            .copy()
        )
        if not d.empty:
            d = d.assign(
                source="hpa_single_cell_tissue",
                version="v23",
                taxonomy_id=9606,
                state="unspecified",
                stage="adult",
                unit="nTPM",
            )
            d = cell_rna_policy.protected_immune_contexts(
                cell_rna_policy.reproductive_immune_contexts(normal_cell_screen.contexts(d))
            )
            rows.append(
                d.sort_values("value", ascending=False).drop_duplicates(["gene_id", "category"])
            )
        if number % 10 == 0:
            print(f"Audited cell chunk {number}", flush=True)
    grouped = pd.read_csv(grouped_path, low_memory=False)
    rows.append(grouped[grouped.gene_id.isin(set(facts.representative_gene_id) - singleton)])
    result = pd.concat(rows, ignore_index=True).sort_values("value", ascending=False)
    result = result[result.mapping_status.eq("mapped")].drop_duplicates(["gene_id", "category"])
    result["protein"] = result.gene_id.map(
        facts.set_index("representative_gene_id").protein_group_label
    )
    write_csv(result, output)
    receipt_path.write_text(
        json.dumps({"request": request, "output_sha256": sha256(output)}, indent=2) + "\n"
    )
    return result


def table(frame):
    def value(v):
        if pd.isna(v):
            return "unmeasured"
        return str(round(v, 3) if isinstance(v, float) else v).replace("|", ";")

    rows = [frame.columns, ["---"] * len(frame.columns), *frame.itertuples(index=False, name=None)]
    return "\n".join("| " + " | ".join(value(v) for v in row) + " |" for row in rows)


def build():
    if policy.POLICY_ID != "oncoref.cta-selection.2026-10-30.1":
        raise RuntimeError(
            "Historical diagnostic: reproduce from commit fcfb7a1. "
            "Use report_cta_125_changes.py for the current categorical release. "
            "Original inputs are preserved in cta-other-20261030/frozen-inputs.zip."
        )
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [CURRENT / "protein-group-classifications.csv.gz", CURRENT / "hpa-diffuseness.csv.gz"]
    facts, scores = [pd.read_csv(p, low_memory=False) for p in paths]
    corrected = bulk_counterfactual(facts, scores)
    original = facts[facts.category.eq("Neither")].copy()
    cells = cell_peaks(original)
    reproductive_cells = cells[cells.category.eq("reproductive") & cells.value.ge(5)]
    matrix = hpa_profiles(facts)
    reference_tissues = sorted(t for ts in HPA_REPRODUCTIVE_BRANCHES.values() for t in ts)
    refs = matrix[reference_tissues]
    top_tissue = refs.loc[refs.notna().any(axis=1)].idxmax(axis=1)
    original["full_reproductive_ntpm"] = original[ID].map(
        scores.set_index("entity_id").reference_value
    )
    original["reproductive_top_tissue"] = original[ID].map(top_tissue)
    original["category_after_bulk_correction"] = original[ID].map(corrected.set_index(ID).category)
    indexed_cells = reproductive_cells.set_index("gene_id")
    for source, destination in [
        ("value", "reproductive_cell_ntpm"),
        ("tissue", "reproductive_cell_tissue"),
        ("cell_type", "reproductive_cell_type"),
        ("Read count", "reproductive_cell_reads"),
    ]:
        original[destination] = original.representative_gene_id.map(indexed_cells[source])
    original["proposed_top_level"] = "Other"
    original["HPA_v23"] = "https://v23.proteinatlas.org/" + original.representative_gene_id
    original["why_not_CTA"] = original.exclusion_reasons.replace(
        {
            "Reproductive / embryonic association": "No qualifying bulk/developmental RNA or reliable reproductive IHC under the current admission rule",
            "Tumor sample ≥2 TPM/pTPM": "No recorded tumor RNA reaches 2 TPM/pTPM",
            "Reproductive concentration": "Measured reproductive/adult profile is too diffuse",
        }
    )
    original.loc[original.category_after_bulk_correction.eq("CTA"), "why_not_CTA"] = (
        "Stale Reproductive bulk input; correcting the tissue set admits this protein"
    )
    original.loc[
        original.category_after_bulk_correction.eq("Neither")
        & original.reproductive_cell_ntpm.notna(),
        "why_not_CTA",
    ] = "Supported reproductive-cell RNA is not an admission route in the current rule"
    write_csv(original, OUT / "all-other-proteins.csv")
    catch_all = original[original.neither_profile.eq("normal_pass_noncandidate")]
    write_csv(catch_all, OUT / "current-generic-bucket.csv")
    restored = original[original.category_after_bulk_correction.eq("CTA")]
    write_csv(restored, OUT / "bulk-correction-CTA.csv")
    remaining = original[original.category_after_bulk_correction.eq("Neither")]
    write_csv(
        remaining[remaining.reproductive_cell_ntpm.notna()],
        OUT / "cell-evidence-reconsideration.csv",
    )
    write_csv(reproductive_cells, OUT / "reproductive-cell-observations.csv")
    tumor = original[original.tumor_restricted_lead]
    write_csv(tumor, OUT / "tumor-supported.csv")
    positive = original[original.positive_context & ~original.reproductive_concentration_hold]
    assert len(positive) == 1 and positive.protein_group_label.iloc[0] == "CDY1/B"
    assert len(catch_all) == 400 and len(restored) == 25 and len(remaining) == 700
    assert len(reproductive_cells) == 113 and remaining.reproductive_cell_ntpm.notna().sum() == 107
    assert len(tumor) == 66 and not tumor.positive_context.any()
    assert tumor.full_reproductive_ntpm.max() == 0.9
    summary = {
        "status": "diagnostic, not released assignments",
        "current_counts": facts.category.fillna("Unassigned").value_counts().to_dict(),
        "bulk_correction_counts": corrected.category.fillna("Unassigned").value_counts().to_dict(),
        "current_reproductive_only": len(positive),
        "tumor_supported": len(tumor),
        "current_generic_bucket": len(catch_all),
        "bulk_correction_promotions": len(restored),
        "supported_reproductive_cell_profiles": len(reproductive_cells),
        "additional_cell_profiles_after_bulk_correction": int(
            remaining.reproductive_cell_ntpm.notna().sum()
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    display = original.set_index("protein_group_label")
    sample = display.loc[
        ["FOXR2", "CAGE1", "PRDM9", "KRT82", "KRT28", "IL1F10", "IFNB1", "OR2T7"],
        [
            "full_reproductive_ntpm",
            "reproductive_cell_ntpm",
            "oncoref_maximum_tpm",
            "healthy_maximum_ntpm",
            "healthy_maximum_tpm",
        ],
    ].reset_index()
    sample.columns = [
        "Protein",
        "Reproductive bulk nTPM",
        "Reproductive cell nTPM (≥5; ≥20 reads)",
        "Tumor TPM",
        "Non-permitted normal nTPM",
        "Non-permitted normal TPM",
    ]
    report = f"""# Other proteins: what the proposed subcategories actually contain

This is a diagnostic audit of the pending 1.25.0 data, **not a released taxonomy change**. The package still calls this top-level category Neither. The proposed user-facing name is **Other**. Counts refer to identical-reference-sequence proteins, not gene loci. Missing required normal evidence remains Unassigned; Other requires a normal-screen pass.

[All 725 current proteins and reasons](all-other-proteins.csv) · [All 400 current generic-bucket proteins](current-generic-bucket.csv) · [66 tumor-supported proteins](tumor-supported.csv) · [25 bulk corrections](bulk-correction-CTA.csv) · [107 further cell-evidence cases](cell-evidence-reconsideration.csv)

## The three names in question

| Proposed description | Count under current rules | Meaning |
|---|---:|---|
| Other:Reproductive | 1 | CDY1/B: qualifying reproductive IHC but no qualifying recorded tumor RNA. This is not the count of all proteins with any biological reproductive expression. |
| Other:Tumor-supported | 66 | Two distinct tumors ≥10 native TPM, low measured normals and no qualifying reproductive evidence; rare-normal-cell/state restriction is not yet established. |
| Existing generic bucket | 400 | This is a legacy fallback, not a biological class. It includes CDY1/B and 25 proteins affected by the stale Reproductive input. It must not simply become Other:Other. |

These three rows are an audit of existing flags, **not a disjoint partition**. The current 725 rows comprise 400 generic, 321 low-normal, three legacy Manageable-expression profiles and one diffuse profile. The 66 tumor-supported leads are contained in the 321 low-normal rows. The three legacy Manageable profiles are OXGR1, KRT9 and TMDD1; they must not be relabelled Disposable.

## Reproductive-only: CDY1/B

The identical-reference-sequence protein is supported by High/Supported testicular IHC at both loci. Its full HPA v23 reproductive bulk maximum is **3.1 nTPM**, compared with the stale input 2.7. OncoRef tumor RNA reaches **1.597 TPM** and the separate HPA tumor lower bound reaches **1.208 pTPM**, below the current one-tumor ≥2 rule. It therefore remains Other under the present cancer-admission rule.

[HPA v23 CDY1B pathology evidence](https://v23.proteinatlas.org/ENSG00000172352-CDY1B/pathology) also describes liver-cancer staining. That is a tumor-protein follow-up, not an automatic RNA-cutoff override: antibody family attribution and the staining localization need review. Reliable tumor-protein admission is currently implemented only in the Tumor-restricted prototype, another asymmetry to resolve before a broad admission change.

## Tumor-supported does not mean reproductive or tumor-exclusive

All 66 fail current Reproductive admission. Their full adult/placental bulk maximum is at most **0.9 nTPM** and available developmental maxima are below **2 TPM**. None has qualifying reproductive IHC or a ≥5 nTPM/≥20-read reproductive cell observation in this audit. Missing optional measurements remain missing, and early developmental coverage is not exhaustive fetal-organ coverage.

Forty labels are olfactory receptors. Other examples are **GFRAL, SCN10A, CNGA2, IFNA4, IFNW1, OTOR and OC90**. The relevant normal contexts include specialized sensory cells, brainstem nuclei and induced immune states. These are tumor-supported leads, not established tumor-exclusive targets. See the [full normal-coverage follow-up](https://github.com/pirl-unc/oncoref/issues/653).

## The generic bucket exposes two admission gaps

**1. Stale bulk tissue definition.** Reproductive admission inherited an older, narrower maximum although the current full HPA profile already includes epididymis and the other mapped branches. MAGEA11 was assessed at 4.3 nTPM even though epididymis is **20.5**. Correcting only this input, while retaining every normal hold, tumor rule and concentration gate, gives **598 CTA / 700 Other / 18,031 Unsafe / 534 Unassigned**. Every existing CTA remains CTA and Unsafe membership is identical.

{table(restored[["protein_group_label", "reproductive_ntpm", "full_reproductive_ntpm", "reproductive_top_tissue"]].rename(columns={"protein_group_label": "Protein", "reproductive_ntpm": "Old input", "full_reproductive_ntpm": "Full Reproductive nTPM", "reproductive_top_tissue": "Tissue"}))}

**2. Supported reproductive-cell RNA is not an admission route.** The current rule uses bulk RNA, supported developmental RNA or reliable reproductive IHC. The existing source cache nevertheless contains **113 current Other proteins** with mapped reproductive-cell RNA ≥5 nTPM and the existing ≥20-read support requirement. Six overlap the 25 bulk corrections, leaving **107 further cases**. This diagnostic evidence flag does not automatically promote them: it needs a consistent positive-association route, cell identity/attribution review and an explicit treatment of concentration when bulk is diluted. Endothelial or other generic cells found inside a reproductive organ must not automatically establish reproductive-specific identity.

{table(sample)}

After removing the 25 bulk-input errors and separating CDY1/B, **374 of the original 400 generic rows remain**. The three legacy Manageable profiles sit outside this count. The residual includes weak bulk reproductive signals, overlooked cell evidence and expression in other normal tissues. 373 of the 374 pass the current cancer gate; **OR2T7** fails both cancer and Reproductive admission. Full row-level explanations and native units are in the linked tables.

## Naming recommendation

Use **CTA / Unsafe / Other**, with Unassigned retained as an assessment status for incomplete or unresolved normal evidence. Within Other:

- **Low** means low measured expression on a stated panel, not never expressed. The existing 321 flag is low *normal* expression and says nothing about cancer abundance.
- **Disposable** requires affirmative Disposable expression and the agreed below-threshold conditions elsewhere. Skin/Manageable expression is not automatically Disposable.
- **Diffuse** means a measured dispersed profile, not a synonym for missing reproductive evidence.
- **Reproductive** means the Reproductive route qualifies but cancer evidence does not.
- **Tumor-supported** means strong cancer evidence without qualifying Reproductive identity or established tumor restriction.

Do not invent an Other:Other biological subtype. Leave Other without a suffix when the evidence does not establish one of these profiles, and show the specific failed or missing admission evidence. Display secondary evidence flags as well as a disjoint primary subtype; a low-normal profile and tumor support can coexist. Final subtype counts require the two admission gaps and the Disposable definition to be resolved first.

## Sources and reproduction

[Human Protein Atlas v23 downloadable datasets](https://v23.proteinatlas.org/about/download) defines separate consensus bulk and tissue-cluster RNA tables. This audit uses the same pinned cache, accepted gene identities, anatomy mappings and supported cell rows as the existing screen; it does not pool bulk and single-cell values into one numerical assay. The ≥5/≥20 diagnostic flag reuses project conventions, not an HPA clinical-safety cutoff.

Run `python scripts/audit_cta_other_subtypes.py`. The bulk counterfactual reruns association, concentration and classification while preserving every original normal-safety mask. Source hashes and output hashes are recorded in [manifest.json](manifest.json). No packaged selection table, report PDF or released artifact is rewritten by this script.
"""
    (OUT / "index.md").write_text(report)
    input_paths = [
        *paths,
        Path(__file__),
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "scripts/cta_profile_data.py",
        CURRENT / "hpa_rna_consensus-accepted-identities.csv.gz",
    ]
    (OUT / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(input_paths),
                "outputs": {
                    p.name: sha256(p) for p in sorted(OUT.iterdir()) if p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    build()
