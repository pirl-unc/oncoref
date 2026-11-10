#!/usr/bin/env python3
"""Publish policy thresholds, named CTA placements and the patient-mapping correction."""

import gzip
import hashlib
import io
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_selection_policy, samples, target_selection  # noqa: E402

BASELINE = "2e24b126"
BASE = ROOT / "docs/audits/target-selection-20261109"
OUT = BASE / "policy-reference"


def old_csv(path):
    data = subprocess.check_output(["git", "show", f"{BASELINE}:{path}"], cwd=ROOT)
    return pd.read_csv(io.BytesIO(gzip.decompress(data)), low_memory=False).fillna("")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    catalog = target_selection.results().fillna("")
    old = old_csv("oncoref/data/target-selection.csv.gz")
    compare = catalog.merge(
        old[["entity_id", "category"]], on="entity_id", suffixes=("", "_before")
    )
    changes = compare[compare.category.ne(compare.category_before)]
    changes.to_csv(OUT / "patient-mapping-category-changes.csv", index=False)
    identities = pd.read_csv(BASE / "rebuild/tumor-samples.csv.gz").fillna("")
    old_ids = old_csv("docs/audits/target-selection-20261109/rebuild/tumor-samples.csv.gz")
    joined = identities.merge(
        old_ids[["sample_key", "participant_id"]], on="sample_key", suffixes=("", "_before")
    )
    gained = joined[joined.participant_id_before.eq("") & joined.participant_id.ne("")]
    coverage = identities.groupby("namespace").agg(
        samples=("sample_id", "size"),
        resolved_samples=("participant_id", lambda s: s.ne("").sum()),
        known_patients=("participant_id", lambda s: s[s.ne("")].nunique()),
    )
    coverage["newly_resolved_samples"] = (
        gained.groupby("namespace").size().reindex(coverage.index, fill_value=0)
    )
    coverage.to_csv(OUT / "patient-coverage.csv")
    old_proposal = old_csv(
        "docs/audits/target-selection-20261109/thresholds/recommended-eligible-inventory.csv.gz"
    )
    proposed = pd.read_csv(BASE / "thresholds/recommended-eligible-inventory.csv.gz").fillna("")
    proposal_changes = proposed.merge(
        old_proposal[["entity_id", "proposed_category", "proposed_classification"]],
        on="entity_id",
        suffixes=("", "_before"),
    )
    proposal_changes = proposal_changes[
        proposal_changes.proposed_category.ne(proposal_changes.proposed_category_before)
    ]
    proposal_changes.to_csv(OUT / "patient-mapping-proposed-changes.csv", index=False)

    # Explicit panel, not a claim that every member of a named family is a CTA.
    named = {
        "NY-ESO-1",
        "CTAG2",
        "CTAGE1",
        "PRAME",
        "CT83",
        "SAGE1",
        "CAGE1",
        "AKAP4",
        "ACRBP",
        "CTCFL",
        "FATE1",
        "HORMAD1",
        "IGF2BP3",
        "LY6K",
        "PLAC1",
        "PBK",
        "SPA17",
        "SPAG9",
        "ROPN1",
        "ROPN1B",
        "ROPN1L",
        "PIWIL1",
        "PIWIL2",
        "PIWIL4",
        "TSGA10",
        "TEX101",
        "DPEP3",
        "HSPB9",
        "SPACA1",
        "LIPI",
        "OIP5",
        "SPANXN1",
        "SPANXN2",
        "SPANXN3",
        "SPANXN4",
        "SPANXN5",
        "LDHC",
        "SYCP2",
        "SYCP1",
        "SYCP3",
        "TPTE",
        "DDX43",
        "CEP55",
        "SPAG1",
        "PAGE4",
    }
    pattern = re.compile(
        r"^(MAGE[ABC]\d|SSX\d+(?:/|$)|GAGE\d|PAGE\d|XAGE\d|BAGE\d|CT45A\d|CT47A\d|SPANX[ABCDN]|TSPY\d|CSAG\d)"
    )
    panel = catalog[
        catalog.protein_group_label.map(lambda s: s in named or bool(pattern.match(s)))
    ].copy()
    with gzip.open(BASE / "rebuild/paper-attestations.json.gz", "rt") as handle:
        papers = json.load(handle)
    proposed_by_id = proposed.set_index("entity_id")

    def links(row):
        genes = set(json.loads(row.member_gene_ids))
        records = [papers[i] for i in json.loads(row.paper_attestation_indices)]
        return sorted(
            {
                (p["citation"], p["source_url"])
                for p in records
                if genes & set(p["nomination_gene_ids"])
            }
        )

    panel["proposed_classification"] = [
        proposed_by_id.loc[r.entity_id, "proposed_classification"]
        if r.entity_id in proposed_by_id.index
        else r.classification
        for r in panel.itertuples()
    ]
    panel["nomination_links"] = [json.dumps(links(r)) for r in panel.itertuples()]
    panel["ensembl_url"] = panel.representative_gene_id.map(
        lambda g: f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={g}"
    )
    panel["hpa_url"] = panel.representative_gene_id.map(
        lambda g: f"https://www.proteinatlas.org/{g}"
    )
    panel = panel.sort_values("protein_group_label")
    panel.to_csv(OUT / "named-cta-proteins.csv", index=False)
    names_table = "\n".join(
        f"| [{r.protein_group_label}]({r.ensembl_url}) · [HPA]({r.hpa_url}) | {r.classification} | {r.proposed_classification} | "
        + (r.Unsafe_reasons or r.decision_reason).replace("|", ";")
        + " | "
        + ("; ".join(f"[{name}]({url})" for name, url in links(r)) or "—")
        + " |"
        for r in panel.itertuples()
    )
    (OUT / "named-ctas.md").write_text(f"""# Named CTA families and studied targets

This explicit panel contains {len(panel)} protein records from the named MAGE-A/B/C, SSX, GAGE, PAGE, XAGE, BAGE, CT45/47, SPANX, TSPY and CSAG families, plus individually listed targets. SSX2IP is not an SSX antigen-family member and is excluded from this panel. BAGE has no resolved record in the current reference identity catalog; no category is inferred for that missing identity. Family membership and historical nomination do not override a normal-expression exclusion. The table uses current protein records, so identical-reference loci can share a row. The CSV retains every member gene, transcript and paper link.

Current means the corrected fivefold-rule draft; proposed means the tenfold single-tumor alternative. The two-patient route is unchanged. Unsafe records retain the existing normal/clinical screen. Gate names identify the recorded cause; they do not assert that clinical injury has been proved for every protein. The current [full inventory](../rebuild/report.html) includes proteins outside this named panel.

[Download all panel data](named-cta-proteins.csv).

| Protein | Current classification | Proposed classification | Recorded reason | Exact-locus published nominations |
|---|---|---|---|---|
{names_table}
""")
    thresholds = [
        (
            "Cancer: repeated",
            "≥2 documented patients; each tumor ≥1 TPM and strictly above N",
            "Current draft",
            "N is the largest supported nonreproductive, non-Disposable normal context statistic; repeated libraries count once.",
        ),
        (
            "Cancer: single tumor",
            "≥5 × max(N, 1) TPM in ≥1 sample",
            "Current draft",
            "One strong sample suffices, including one without a resolved patient ID.",
        ),
        (
            "Reproductive bulk RNA",
            "≥1 native nTPM",
            "Current draft",
            "Includes reproductive, placental and developmental scope; normal exclusions still apply.",
        ),
        (
            "Reproductive-cell RNA",
            "≥1 nTPM and ≥20 reads",
            "Current draft",
            "Indirect sites also require complete outside-cluster coverage, every outside value <1, and no curated normal-induction conflict.",
        ),
        (
            "Developmental RNA",
            "Mean ≥1 TPM; ≥20 cells, ≥5 positive and ≥10% positive; positive means ≥1 TPM",
            "Current draft",
            "Qualifying developmental lineage/stage; bulk fetal/reproductive scope is retained separately.",
        ),
        (
            "Reproductive protein",
            "Reliable mapped positive IHC or reviewed human reproductive protein for every member",
            "Current draft",
            "An alternative to RNA association; antibody/identity conflicts remain unresolved.",
        ),
        (
            "Normal RNA ceilings",
            "Vital 4/8/16/32; Protected 8/16/32/64; Manageable 16/32/64/128",
            "Retained normal screen",
            "Columns are missing/uncertain, Approved, Supported, Enhanced matched negative IHC. Units remain native TPM or nTPM.",
        ),
        (
            "Reviewed clinical RNA credit",
            "One IHC rank increment, capped at rank 3",
            "Retained normal screen",
            "Only qualifying product-specific evidence; no stacking or override of injury/protein exclusions.",
        ),
        (
            "Normal RNA consequence",
            "Above Vital/Protected ceiling → Unsafe; above Manageable ceiling → hold; Manageable >128 → Unsafe",
            "Retained normal screen",
            "Equality to the ceiling passes. Disposable and reproductive scope have no current numerical ceiling in this gate.",
        ),
        (
            "HPA tissue-cell support",
            "≥20 reads",
            "Retained normal screen",
            "Cell-specific tier; only exact/declared compartment-matched IHC can supply credit.",
        ),
        (
            "HPA/Monaco immune",
            "Cell-specific 4/8/16 bases on native nTPM/TPM",
            "Retained normal screen",
            "No matched negative-IHC credit; a reproductive location does not exempt shared immune cells.",
        ),
        (
            "DICE immune",
            "Third-highest reported TPM; ≥3 observations",
            "Retained normal screen",
            "Apply cell-specific ceiling; unpaired distributions keep their bounds.",
        ),
        (
            "Adult cornea",
            "≥5 positive cells and ≥20 counts in ≥2 sample groups; RNA ceilings divided by 100",
            "Retained normal screen",
            "Native CP10K screening scale, not a TPM conversion.",
        ),
        (
            "Rare-cell detection",
            "≥1% detected, ≥3,000 cells and ≥10 expressing cells",
            "Retained normal screen",
            "Existing brain/spinal cord/heart/lung panel: Vital → Unsafe; Protected → review hold.",
        ),
        (
            "Normal protein",
            "Reliable non-Low positive in Vital/Protected → Unsafe; Manageable → hold",
            "Retained normal screen",
            "Low is flagged; uncertain is unresolved. Reviewed/external findings retain reagent, cell and identity scope.",
        ),
        (
            "Required RNA coverage",
            "Every nonreproductive combined-bulk tissue concept for every member",
            "Retained normal screen",
            "The comparison uses 36 concepts. Missing optional cell/immune/cornea panels alone does not block; observed exclusions still apply.",
        ),
        (
            "Translation entry",
            "Reference coding CDS, or supported ncORF: MS peptide or ≥2 independent Ribo-seq studies",
            "Current draft",
            "All proteins still undergo the normal screen.",
        ),
        (
            "Representative translation",
            "Complete primary-locus coding sequence; no LoF/NMD, incomplete tags, internal stops or non-AA20 residues",
            "Current draft",
            "MANE Select first, Ensembl canonical second; conflicting evidence requires review. X is not allowed.",
        ),
        (
            "Nested coding sequence",
            "<50 amino acids inside a same-strand coding gene without independent ncORF evidence",
            "Current draft",
            "Other:Sequence.",
        ),
        (
            "Unsafe sequence overlap",
            "≥95% of 9-mers occur in Unsafe proteins",
            "Current draft",
            "Other:Sequence; this does not create expression-based Unsafe evidence.",
        ),
        (
            "Reproductive concentration",
            "Recorded normal/(normal + reproductive maximum); previous >10% veto removed",
            "Descriptive only",
            "Developmental and adult reproductive association now share the same admission flow.",
        ),
        (
            "Published nomination",
            "Exact-member claim in a curated CA/CTA/oncofetal target source",
            "Annotation",
            "Retain source, table and claim; neither automatic admission nor exemption from normal findings.",
        ),
        (
            "Proposed single-tumor rule",
            "≥10 × max(N, 1) TPM",
            "Counterfactual; not adopted",
            "Repeated route remains ≥2 patients at ≥1 TPM and above N.",
        ),
        (
            "Proposed target maximum",
            "max(cancer, reproductive/developmental) ≥1 and strictly above every measured nonreproductive value",
            "Counterfactual; not adopted",
            "Source-compatible values only; missing compartments are not zeros.",
        ),
        (
            "Proposed Disposable caps",
            "≤32 TPM and ≤32 nTPM, assessed separately",
            "Counterfactual; not adopted",
            "No unit conversion or IHC-derived numerical allowance for these proposed caps.",
        ),
        (
            "Percentile alternatives",
            "Within-tumor p80/p90/p95; across-patient p80/p90/p95 with ≥20 measured patients",
            "Sensitivity analysis",
            "Within-tumor ranks use the full measured catalog; across-patient values use the lower observed order statistic.",
        ),
    ]
    table = pd.DataFrame(thresholds, columns=["rule", "threshold", "status", "meaning"])
    table.to_csv(OUT / "thresholds.csv", index=False)
    threshold_md = "\n".join("| " + " | ".join(row) + " |" for row in thresholds)
    counts = catalog.category.value_counts()
    changes_md = "\n".join(
        f"| {r.protein_group_label} | {r.category_before} | {r.category} | {r.independent_patients_above_normal:g} |"
        for r in changes.itertuples()
    )
    coverage_md = "\n".join(
        f"| {name} | {r.newly_resolved_samples} |"
        for name, r in coverage[coverage.newly_resolved_samples.gt(0)].iterrows()
    )
    (OUT / "index.md").write_text(f"""# Thresholds, patient identity and named CTA placements

**Corrected draft: {counts["CTA"]} CTA, {counts["Cancer Antigen only"]} Cancer Antigen only, {counts["Reproductive Antigen only"]} Reproductive Antigen only, {counts["Other"]} Other and {counts["Unsafe"]:,} Unsafe.** These use the existing fivefold single-tumor rule. [Every protein](../rebuild/report.html) · [Named CTA families and studied targets](named-ctas.md) · [All threshold values as CSV](thresholds.csv).

## Samples and patients

OncoRef's cancer matrices have one column per sample/library. The previous comparison resolved only TCGA, TARGET and MMRF patient barcodes, omitting patient metadata already available in OncoRef. `expression.sample_metadata` now joins actual matrix columns to the central `samples.sample_identities` policy. The rebuild and threshold analysis use that same mapping. Source-specific patient IDs, physical sample keys, status, basis and references stay separate.

The correction resolves **{len(gained):,} additional samples**. **{identities.participant_id.ne("").sum():,} of {len(identities):,} sample columns** have documented identities, representing **{identities.loc[identities.participant_id.ne(""), "participant_id"].nunique():,} patient IDs**; **{identities.participant_id.eq("").sum():,} samples** remain unresolved. Identities use each source's documented patient namespace; unavailable cross-study linkages remain unknown. A copied sample label in a GEO case field is not accepted as a patient mapping.

| Source namespace | Additional resolved samples |
|---|---:|
{coverage_md}

The [CLL-map source](https://www.cllmap.org/downloads.html) specifies six duplicate exclusions plus one suspected MCL exclusion, already applied to OncoRef's 708-sample cohort. Treehouse uses curated donor IDs where present and the existing builder's donor-family rule for older identifiers. Multiple S01/S02 libraries from one family count once. [Complete source coverage](patient-coverage.csv).

**Multiple samples are not mandatory.** One qualified tumor ≥5 × max(N,1) still qualifies. Patient identity is needed only to establish independent support for the repeated route and patient-based percentiles. Unknown IDs do not erase measured expression.

| Protein | Before correction | Corrected category | Patients above normal |
|---|---|---|---:|
{changes_md}

The expression values and normal exclusions did not change. [Current-rule changes](patient-mapping-category-changes.csv) and [tenfold-alternative changes](patient-mapping-proposed-changes.csv) retain complete evidence. The baseline is [commit {BASELINE}](https://github.com/pirl-unc/oncoref/commit/{BASELINE}); the [threshold comparison](../thresholds/index.md) has been regenerated with the corrected identities.

## Threshold inventory

These are the rules used in the current draft or inherited normal screen. Proposed alternatives are labeled separately. Historical 1.36.0 used a ≥5 reproductive RNA threshold, tumor presence >0 and the >10% concentration veto; those are not the current draft's association rules. Older compatibility taxonomy parameters do not replace the active `cta_selection_policy` ceilings.

| Rule | Threshold | Status | Meaning |
|---|---|---|---|
{threshold_md}

The base normal ceilings are **4 for Vital, 8 for Protected and 16 for Manageable**, before qualifying IHC or clinical adjustments. They are screening parameters, not measured injury thresholds. [Normal-screen implementation](https://github.com/pirl-unc/oncoref/blob/feat/antigen-set-policy/oncoref/cta_selection_policy.py), [antigen membership implementation](https://github.com/pirl-unc/oncoref/blob/feat/antigen-set-policy/oncoref/target_selection.py) and the [tissue registry with cited rationales](../../../tissue-hazards.md) identify the distinct policy layers.

## Hair follicles

**Hair follicles are already classified Disposable; whole skin is Protected.** The missing component in the numerical cap comparison is isolated hair-follicle RNA, not a missing ontology bucket. HPA v23's tissue-cell RNA file has no hair-follicle cluster. It does have cutaneous melanocytes on its separate nTPM scale; the [cell audit](../expression-maxima/index.md) retains read support and missing loci.

HPA v23 also contains hair protein staining: cortex/medulla, cuticle, internal/external root sheath, and hair-follicle cells in skin sections. Thus “no hair data” was too broad. A separate mapping defect was found: the skin cell-label regex `hair follic` does not match `hair follicles` and falls back to Protected skin. [Issue #687](https://github.com/pirl-unc/oncoref/issues/687) tracks correction together with a replay of the inherited normal-protein gates and affected sequence containment; the patient-only correction here does not revise those gates.

Human hair-follicle RNA datasets are available, including [GSE129611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129611) and [the human scalp single-cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9126928/). They have not been added to the pinned cap panel. Missing hair RNA cannot count as zero Disposable expression.
""")
    assert cta_selection_policy.TIER_BASES == {"vital": 4, "protected": 8, "manageable": 16}
    manifest = {
        "baseline_commit": BASELINE,
        "patient_identity_policy": samples.PATIENT_IDENTITY_POLICY,
        "inputs_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [
                ROOT / "oncoref/data/target-selection.csv.gz",
                BASE / "rebuild/tumor-samples.csv.gz",
                BASE / "thresholds/recommended-eligible-inventory.csv.gz",
                BASE / "rebuild/paper-attestations.json.gz",
            ]
        },
        "implementation_sha256": {"scripts/build_target_policy_reference.py": sha(Path(__file__))},
        "outputs_sha256": {p.name: sha(p) for p in OUT.iterdir() if p.name != "manifest.json"},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "named_proteins": len(panel),
                "category_changes": len(changes),
                "newly_resolved_samples": len(gained),
            }
        )
    )


if __name__ == "__main__":
    build()
