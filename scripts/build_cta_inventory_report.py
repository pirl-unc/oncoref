#!/usr/bin/env python3
"""Account for the complete candidate inventory and bound expression-only claims."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_cta_contract_report import digest, write_table  # noqa: E402

from oncoref import cta_selection, gene_knowledge  # noqa: E402

OUT = ROOT / "docs/audits/cta-inventory-20261105"
REPLAY = ROOT / "docs/audits/cta-selection-20261028"


def build():
    facts = pd.read_csv(REPLAY / "protein-group-classifications.csv.gz", low_memory=False)
    old = pd.read_csv(OUT / "baseline-1.30.0.csv.gz", low_memory=False)
    former = pd.read_csv(OUT / "formerly-out-of-universe.csv.gz")
    restored = facts[
        facts.protein_equivalence_group_id.isin(former.protein_equivalence_group_id)
    ].copy()
    historical_path = ROOT / "docs/audits/cta-unassigned-20261029/before-diffuseness.csv.gz"
    historical = pd.read_csv(historical_path, low_memory=False).set_index(
        "protein_equivalence_group_id"
    )
    for column in ["category", "Unsafe_reasons"]:
        restored["historical_pre_diffuseness_" + column] = restored[
            "protein_equivalence_group_id"
        ].map(historical[column])
    assert len(restored) == len(former) == 259
    assert not restored.category.eq("CTA").any()
    assert set(restored.classification) <= {"Unsafe", "Other:Translation"}
    write_table(restored, OUT / "restored-records.csv.gz")
    historical_reassessments = restored[
        restored.historical_pre_diffuseness_category.eq("Unsafe") & ~restored.category.eq("Unsafe")
    ]
    current = facts.set_index("protein_equivalence_group_id")
    previous = old.set_index("protein_equivalence_group_id")
    changed = current.loc[previous.index].classification.ne(previous.classification)
    write_table(
        current.loc[changed.index[changed]].reset_index(),
        OUT / "existing-classification-changes.csv.gz",
    )

    # A stored zero is a property of the supplied summaries, not absence in every
    # donor/cell. Include reproductive/placental bulk; missing never equals zero.
    rna = facts.required_coverage & facts.cancer_association
    adult_zero = rna & facts.normal_maximum_floor_ratio.eq(0) & facts.reproductive_ntpm.eq(0)
    development_known = adult_zero & facts.developmental_measured
    zero_rna = development_known & facts.developmental_all_lineage_max_tpm.eq(0)
    zero_with_no_protein = zero_rna & ~facts.any_detected_protein
    write_table(facts.loc[adult_zero], OUT / "adult-zero-rna-records.csv.gz")
    profiles = cta_selection.profiles()
    leads = profiles[profiles.preliminary_tumor_restriction_profile]
    write_table(leads, OUT / "preliminary-expression-profile-matches.csv.gz")
    conflicts = leads[leads.category.eq("Unsafe")]
    reasons = []
    findings = gene_knowledge.findings()
    for row in conflicts.itertuples():
        members = set(row.gene_ids)
        relevant = findings[
            findings.gene_ids.map(lambda ids, members=members: bool(set(ids) & members))
        ]
        for finding in relevant.itertuples():
            reasons.append(
                {
                    "label": row.label,
                    "classification": row.classification,
                    "finding_id": finding.id,
                    "effect": finding.effect,
                    "finding": finding.finding,
                    "decision_reason": finding.decision_reason,
                    "source_url": finding.source_url,
                    "limitations": finding.limitations,
                }
            )
    write_table(pd.DataFrame(reasons), OUT / "profile-counterevidence.csv")
    summary = {
        "policy_id": cta_selection.coverage()["policy_id"],
        "input_records": len(facts),
        "counts": facts.category.value_counts().to_dict(),
        "restored_records": len(restored),
        "restored_classifications": restored.classification.value_counts().to_dict(),
        "existing_classification_changes": int(changed.sum()),
        "restored_historical_exclusions_reassessed": historical_reassessments[
            ["protein_group_label", "historical_pre_diffuseness_Unsafe_reasons", "classification"]
        ].to_dict("records"),
        "stored_zero_rna_screen": {
            "cancer_present_required_bulk_complete": int(rna.sum()),
            "adult_and_reproductive_bulk_zero": int(adult_zero.sum()),
            "also_developmental_measured": int(development_known.sum()),
            "developmental_measurement_missing": int(
                (adult_zero & ~facts.developmental_measured).sum()
            ),
            "also_developmental_zero": int(zero_rna.sum()),
            "also_no_recorded_adult_protein": int(zero_with_no_protein.sum()),
            "interpretation": "Stored RNA summary values equal zero; no claim of whole-protein absence or exhaustive normal-cell coverage.",
        },
        "preliminary_expression_profile": {
            "matches": len(leads),
            "classification_counts": leads.classification.value_counts().to_dict(),
            "disposition_counts": leads.tumor_restriction_status.value_counts().to_dict(),
            "interpretation": "Existing numerical screen: replicated tumor RNA at least 10 TPM in two samples (or reviewed tumor protein), measured adult RNA below nominal 2-unit assay floors, no positive adult protein call in the profile's measurement inputs, measured developmental lineage means below 2 TPM, and no qualifying reproductive/developmental association. Curated knowledge and rare-cell counterevidence is assessed separately and overrides the screen.",
        },
        "limitations": [
            "Required bulk coverage does not establish coverage of every cell type, age, physiological state or protein assay.",
            "Routine tissue summaries and developmental lineage means are not maxima across all donors or individual cells.",
            "Small positive developmental RNA values can reflect assay/mapping uncertainty; stored-zero screening is not a biological detection limit.",
            "Missing measurements are reported as gaps, never imputed to zero.",
            "No numerical RNA profile establishes translation, HLA presentation, receptor specificity or clinical safety.",
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    zero = summary["stored_zero_rna_screen"]
    restored_counts = summary["restored_classifications"]
    labels = {
        "normal_coverage_review_required": "Normal-context coverage review required",
        "normal_evidence_hold": "Unresolved normal evidence",
        "excluded_by_normal_evidence": "Unsafe: normal evidence",
        "translation_review": "Translation review",
        "reviewed_panel_restriction": "Restriction reviewed within named panels",
        "incomplete_evidence": "Incomplete required evidence",
    }
    dispositions = "\n".join(
        f"| {labels[name]} | {n} |"
        for name, n in summary["preliminary_expression_profile"]["disposition_counts"].items()
    )
    counter = "\n".join(
        f"| {r['label']} | {r['decision_reason']} [Source]({r['source_url']}) |" for r in reasons
    )
    text = f"""# Complete candidate inventory and cancer-restriction accounting

All **{len(facts):,} candidate protein records** enter the main assessment under
`{summary["policy_id"]}`. Identical selected reference sequences count together;
member loci and original sequences remain traceable.

## Translation is an admission gate

The 259 formerly omitted records now receive **{restored_counts.get("Other:Translation", 0)} Other:Translation**
and **{restored_counts.get("Unsafe", 0)} Unsafe** calls. None becomes CTA. Normal exclusion takes
precedence; lack of translation support cannot erase it. Other:Translation requires
sequence review and supplies no direct-target allowance. All {len(old):,} previously
classified records retain their classifications: **{int(changed.sum())} changes**.

The older pre-diffuseness snapshot called {len(historical_reassessments)} of these restored
records Unsafe through DICE immune evidence: {", ".join(historical_reassessments.protein_group_label)}.
They are Other:Translation under the current policy. The export preserves their
historical categories and reasons alongside the current calls; restoring records
does not reinstate superseded thresholds.

[Every restored record and all reasons](restored-records.csv.gz) ·
[Historical omission ledger](formerly-out-of-universe.csv.gz) ·
[Full current classification](../cta-selection-20261028/protein-group-classifications.csv.gz).

The retained translation rule accepts qualifying reference coding annotation or
ncORF support from mass-spectrometry peptides or at least two independent Ribo-seq
studies, followed by exact-record sequence qualification. Annotation is not
experimental proof of production. These two questions remain explicit in the API.

## What can we count as cancer-only with the available data?

The data permit a reproducible **observed-expression screen**, with coverage limits:

| Stored-zero RNA screen | Records |
|---|---:|
| Cancer expression and complete required normal bulk coverage | {zero["cancer_present_required_bulk_complete"]:,} |
| Also zero recorded adult RNA and reproductive/placental bulk RNA | {zero["adult_and_reproductive_bulk_zero"]} |
| Also have developmental measurements | {zero["also_developmental_measured"]} |
| Also zero recorded developmental RNA | {zero["also_developmental_zero"]} |
| Also no positive adult protein call in the measurement inputs | {zero["also_no_recorded_adult_protein"]} |

**{zero["developmental_measurement_missing"]}** adult-zero records lack developmental measurements;
they remain unresolved for that comparison. [All adult-zero records](adult-zero-rna-records.csv.gz).
A zero here is a stored tissue/lineage summary, not absence from every individual
cell or state. Tiny positive values are also not automatically validated protein
production. RNA and protein claims must remain separate.

## Why an expression-profile match can still be Unsafe

The existing, less literal screen returns **{len(leads)} preliminary matches**. It uses
replicated tumor expression (at least 10 TPM in two samples, or reviewed tumor
protein), low measured adult RNA, no positive adult protein call in its measurement inputs, and low measured
developmental RNA. Its nominal adult floor is 2 TPM/nTPM, with the separately
specified 0.02 CP10K corneal floor; these scales are not pooled or converted.

These are numerical matches before all curated and rare-cell evidence is applied.
In particular, the profile's protein-detection flag does not incorporate the
separate curated normal-cell knowledge ledger; those findings can establish
normal production even when the profile flag is negative.
They are not {len(leads)} established cancer-only proteins:

| Disposition after independent evidence | Records |
|---|---:|
{dispositions}

[Every match, identity, classification and review status](preliminary-expression-profile-matches.csv.gz).
The public API exposes `preliminary_tumor_restriction_profile` separately from
`tumor_restriction_status`; an Unsafe record is explicitly `excluded_by_normal_evidence`.
The former undifferentiated `tumor_restricted_lead` evidence flag is no longer exported.

The following reviewed findings demonstrate why a low routine-panel maximum
cannot establish absence from a rare normal population or an activated state:

| Protein | Counterevidence |
|---|---|
{counter}

[Counterevidence with source links and limitations](profile-counterevidence.csv).
None of these numerical screens establishes antigen presentation or therapeutic safety.

[Machine-readable summary](summary.json) · [Reproduction receipt](manifest.json).
"""
    (OUT / "index.md").write_text(text)
    inputs = [
        Path(__file__),
        REPLAY / "protein-group-classifications.csv.gz",
        OUT / "baseline-1.30.0.csv.gz",
        OUT / "formerly-out-of-universe.csv.gz",
        historical_path,
        ROOT / "oncoref/data/cta-selection.json",
        ROOT / "oncoref/data/gene-knowledge.json",
        ROOT / "scripts/build_cta_selection_release.py",
    ]
    receipt = {
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "outputs": {
            p.name: digest(p)
            for p in sorted(OUT.iterdir())
            if p.name
            not in {"manifest.json", "baseline-1.30.0.csv.gz", "formerly-out-of-universe.csv.gz"}
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    build()
