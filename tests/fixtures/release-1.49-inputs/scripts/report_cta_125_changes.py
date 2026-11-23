#!/usr/bin/env python3
"""Generate the 1.25 categorical release comparison from pinned decisions."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from cta_report_common import input_hashes, sha256, write_csv

from oncoref.cta_selection_policy import OTHER_SUBTYPE_RULES, POLICY_ID

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-neither-20261030"
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
ID = "protein_equivalence_group_id"


def table(frame):
    rows = [frame.columns.tolist(), ["---"] * len(frame.columns)]
    rows.extend(frame.fillna("unmeasured").itertuples(index=False, name=None))
    return "\n".join("| " + " | ".join(str(v).replace("|", ";") for v in r) + " |" for r in rows)


def build():
    paths = [
        OUT / "baseline-1.24.0.csv.gz",
        CURRENT / "protein-group-classifications.csv.gz",
        CURRENT / "reproductive-protein-evidence.csv.gz",
        CURRENT / "reproductive-cell-evidence.csv.gz",
        CURRENT / "tumor-restricted-observations.csv.gz",
    ]
    before, after, ihc, cells, tumors = [pd.read_csv(p, low_memory=False) for p in paths]
    before, after = before.set_index(ID), after.set_index(ID)
    assert set(before.index) == set(after.index)
    after = after.reindex(before.index)
    old, new = before.category.fillna("Unassigned"), after.category.fillna("Unassigned")
    assert after.loc[old.eq("CTA"), "category"].eq("CTA").all()
    assert set(before.index[old.eq("Unsafe")]) == set(after.index[new.eq("Unsafe")])
    assert after.loc[new.eq("Other"), "subtype"].notna().all()
    assert (
        after.loc[new.eq("Other"), ["healthy_maximum_ntpm", "healthy_maximum_tpm"]]
        .max()
        .le(30)
        .all()
    )
    counts = pd.DataFrame({"Category": ["CTA", "Other", "Unsafe", "Unassigned"]})
    counts["1.24.0"] = (
        counts.Category.map(old.replace("Neither", "Other").value_counts()).fillna(0).astype(int)
    )
    counts["1.25.0"] = counts.Category.map(new.value_counts()).fillna(0).astype(int)
    transitions = pd.crosstab(old.rename("1.24"), after.classification.rename("1.25")).reset_index()
    write_csv(transitions, OUT / "transitions-1.25.0.csv")
    changes = after.loc[old.replace("Neither", "Other").ne(new)].copy()
    changes.insert(0, "category_1_24", old.loc[changes.index])
    write_csv(changes.reset_index(), OUT / "changes-1.25.0.csv.gz")
    columns = [
        "protein_group_label",
        "classification",
        "subtype",
        "member_gene_ids",
        "reproductive_ntpm",
        "reproductive_cell_ntpm",
        "reproductive_cell_indirect_ntpm",
        "reproductive_cell_indirect_association",
        "outside_cluster_max_ntpm",
        "outside_cluster_complete",
        "normal_induction_evidence",
        "reproductive_protein_association",
        "developmental_mean_tpm",
        "oncoref_maximum_tpm",
        "hpa_tumor_group_lower_bound_ptpm",
        "hpa_diffuseness",
        "adult_profile_status",
        "healthy_maximum_ntpm",
        "healthy_maximum_tpm",
        "healthy_maximum_cp10k",
        "exclusion_reasons",
        "optional_DICE_complete",
        "optional_cornea_complete",
    ]
    exported = after[columns].sort_values(["classification", "protein_group_label"]).copy()
    exported["category_definition"] = exported.subtype.map(OTHER_SUBTYPE_RULES)
    exported["HPA_v23"] = "https://v23.proteinatlas.org/" + after.representative_gene_id
    write_csv(exported.reset_index(), OUT / "protein-categories-1.25.0.csv")
    restored = after.loc[new.eq("CTA") & ~old.eq("CTA")].copy().sort_values("protein_group_label")
    write_csv(restored.reset_index(), OUT / "restored-CTA-1.25.0.csv")
    restored_ids = {g for ids in restored.member_gene_ids for g in json.loads(ids)}
    write_csv(
        ihc[ihc.qualifying_reproductive_positive & ihc.gene_id.isin(restored_ids)],
        OUT / "restored-reproductive-IHC-1.25.0.csv",
    )
    write_csv(
        cells[cells.gene_id.isin(restored.representative_gene_id)],
        OUT / "restored-reproductive-cell-1.25.0.csv",
    )
    other = after[new.eq("Other")]
    category_rows = []
    for name, description in OTHER_SUBTYPE_RULES.items():
        part = other[other.subtype.eq(name)].sort_values("protein_group_label")
        write_csv(exported.loc[part.index].reset_index(), OUT / f"other-{name.lower()}-1.25.0.csv")
        category_rows.append(
            {
                "Category": "Other:" + name,
                "Proteins": len(part),
                "Definition": description,
                "Examples": ", ".join(part.protein_group_label.head(5)),
            }
        )
    categories = pd.DataFrame(category_rows)
    write_csv(categories, OUT / "other-categories-1.25.0.csv")
    leads = after[after.subtype.eq("Tumor-supported")]
    obs = tumors[tumors[ID].isin(leads.index)].copy()
    assert obs.groupby(ID).sample_key.nunique().eq(2).all()
    obs.insert(1, "protein", obs[ID].map(after.protein_group_label))
    write_csv(obs, OUT / "tumor-restricted-lead-observations-1.25.0.csv")
    write_csv(exported.loc[leads.index].reset_index(), OUT / "tumor-restricted-leads-1.25.0.csv")
    held = after[old.eq("Neither") & new.eq("Unassigned")]
    write_csv(held.reset_index(), OUT / "previous-neither-normal-holds-1.25.0.csv")
    named = after[
        after.protein_group_label.isin(
            [
                "MAGEA4",
                "NY-ESO-1",
                "PRAME",
                "CTAG2",
                "XAGE1A/B",
                "MAGEA11",
                "FOXR2",
                "CAGE1",
                "PRDM9",
                "CDY1/B",
                "OR2T7",
                "IFNA2",
                "IFNA14",
                "XKR5",
                "MTCP1",
                "USP17L3",
                "FAM90A9",
            ]
        )
    ]
    named = exported.loc[
        named.index,
        [
            "protein_group_label",
            "classification",
            "reproductive_ntpm",
            "reproductive_cell_ntpm",
            "reproductive_cell_indirect_ntpm",
            "oncoref_maximum_tpm",
        ],
    ].round(3)
    report = f"""# OncoRef 1.25: categorical CTA, Unsafe and Other decisions

**{int(new.eq("CTA").sum())} CTA, {len(other)} Other, {int(new.eq("Unsafe").sum()):,} Unsafe and {int(new.eq("Unassigned").sum())} unassessed proteins.** All {int(old.eq("CTA").sum())} previous CTAs remain included; Unsafe membership is unchanged. This release adds {len(restored)} CTA proteins. Counts collapse identical selected reference protein sequences; they are not gene-locus counts.

{table(counts)}

[Full report and funnel](../cta-source-walkthrough-20261028/index.md) · [PDF logic slides](../cta-source-walkthrough-20261028/cta-selection-1.25.0.pdf) · [Every protein, category and measurement](protein-categories-1.25.0.csv) · [Every transition](transitions-1.25.0.csv) · [Added CTAs](restored-CTA-1.25.0.csv)

## One category per assessed protein

The categorical label is `classification`: CTA, Unsafe, or one of the named Other categories below. `category` retains the broad family and `subtype` the Other name. Every Other protein receives one subtype; there is no Other:Other or unnamed remainder. Insufficient normal evidence is explicitly Unassessed, not a normal-pass category.

{table(categories)}

These categories are assigned in the table's order after safety and CTA admission. Secondary evidence flags remain available. For example, a protein can have both low normal expression and strong tumor evidence, but receives Other:Tumor-supported. Descriptive 2-unit floors mean native TPM/nTPM, with the pre-existing separate nominal CP10K floor; the assays are not pooled or converted into equivalent biological measurements. Low does not mean never expressed. Normal and Inducible do not mean demonstrated clinical injury. All Other proteins pass every recorded normal gate; supported non-permitted maxima remain below 30 TPM/nTPM.

Raw maxima are not always supported expression. **OR2B6 and SCGB1C2** each have a 47.6 nTPM platelet-cluster value supported by only 5 reads (the rule requires 20). **OR6C2** has a 30.92 TPM DICE maximum, but the source-specific third-highest value is 1.14 TPM. These remain **Other:Weak**, with the raw observations retained, rather than being described as never expressed. This support distinction applies consistently across proteins.

## The repaired admission rules

1. **Normal safety:** unchanged tissue/IHC RNA ceilings, protein rules and clinical exclusions. Required normal coverage and evidence resolution remain mandatory. No reproductive or cancer signal overrides a normal exclusion.
2. **Cancer presence:** one sample ≥2 native TPM or HPA pTPM, or reviewed reliable tumor protein. The protein route is now shared by CTA admission and tumor-support assessment. Its current review ledger has no accepted rows; aggregate cancer-IHC counts alone do not qualify.
3. **Reproductive presence:** full mapped HPA v23 reproductive bulk maximum ≥5 nTPM, supported reproductive-cell RNA ≥5 nTPM with ≥20 reads, reliable positive reproductive IHC, or supported developmental RNA ≥5 TPM. Bulk, cells and developmental values retain separate scales and provenance. All matched protein-member loci must be present before summing a normal cell context.
4. **Reproductive concentration:** the existing complete HPA bulk contrast remains active, including Disposable and excluding Permissive thymus. Cell evidence cannot override a measured concentration failure. A missing/zero bulk reference remains explicitly unassessed rather than being replaced by cell RNA.

The full bulk tissue definition fixes the MAGEA11 example: **20.5 nTPM in epididymis**, not the stale 4.3 input. The cell route recovers direct signals such as **CAGE1 208.9 nTPM in early spermatids**, **PRDM9 61.8 in spermatocytes**, and **FOXR2 5.2 in spermatogonia**. [Supporting cell rows and read counts](restored-reproductive-cell-1.25.0.csv) · [Supporting IHC rows](restored-reproductive-IHC-1.25.0.csv).

### Generic cells and macrophages

A generic, immune or mixed-cell signal inside reproductive tissue may provide **indirect sample-site support** when it reaches ≥5 nTPM with ≥20 reads, the same-source outside-reproductive cluster panel is complete, and every outside cluster is <5 nTPM. This includes all measured outside cell types, not only the matching cell label. All normal-source and concentration gates still apply. It establishes the sampling site, not which cell produced the transcript. [Single-cell RNA sequencing uncovers heterogenous transcriptional signatures in macrophages during efferocytosis (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7459098/) documents passenger transcripts in phagocytes; this is an interpretation limit, not proof of the mechanism for any gene here.

Published normal induction prevents treating an indirect-only atlas signal as reproductive restriction. **IFNA2 and IFNA14** have high breast-macrophage atlas signals, but [Plasmacytoid dendritic cells respond to Epstein-Barr virus infection with a distinct type I interferon subtype profile (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6457222/) reports expression in activated human nonreproductive immune cells. They therefore become **Other:Inducible**, not CTA. The curated source table also records IFNA5, IFNA7 and IFNA17. This qualitative counterevidence is not fabricated TPM and does not itself supply Unsafe evidence. Unreviewed family members do not inherit it by name or sequence similarity.

## Named targets and boundary cases

{table(named)}

MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B remain CTA. CDY1/B has reliable reproductive IHC but fails the cancer-RNA floor; OR2T7 now has reproductive-cell support but also lacks qualifying cancer RNA. HPA tumor staining for CDY1/B remains a protein-evidence follow-up because localization and family attribution need adjudication; the ledger does not silently accept it. [HPA v23 CDY1B pathology](https://v23.proteinatlas.org/ENSG00000172352-CDY1B/pathology).

## Tumor-supported and normal-coverage caveats

{len(leads)} Other:Tumor-supported proteins have two qualifying ≥10 TPM tumors. They lack qualifying Reproductive identity and are not declared tumor-exclusive. Rare sensory cells and activated immune states remain relevant coverage gaps. [Every qualifying tumor observation](tumor-restricted-lead-observations-1.25.0.csv) · [Coverage follow-up #653](https://github.com/pirl-unc/oncoref/issues/653). Missing optional matrices are not measured zeros.

The seven previous Neither proteins with normal holds remain unassessed: {", ".join(held.protein_group_label)}. Their older narrow profile cannot override the normal gates. Historical audits and the separate exact-protein/k-mer policy remain frozen.

## API and reproducibility

```python
from oncoref import cta_selection

cta_selection.results(category="CTA")
cta_selection.results(category="Other", subtype="Inducible")
cta_selection.results()[["label", "classification", "reasons"]]
cta_selection.candidates()  # CTA only
cta_selection.audit().explain("FOXR2")
```

The former `category="Neither"` query is a deprecated alias for Other; returned labels use Other. The proposed but unreleased Tumor-restricted top-level category is removed. Dog and mouse keep the same schema with explicit Unassessed status; human admission is not transferred by orthology. Pin policy **{POLICY_ID}** and the source snapshot. The existing legacy CTA and exact-sequence exclusion APIs remain separate products, so downstream code must not silently change k-mer backgrounds.

Run `scripts/regenerate_cta_walkthrough.py --verify-replay`, build the PDF/PPTX slides, finalize the receipts, and run this script. Inputs and output checksums are in [update-manifest.json](update-manifest.json). The [pre-fix diagnostic audit](../cta-other-20261030/index.md) is historical; its counts are not this release.
"""
    (OUT / "update-1.25.0.md").write_text(report)
    outputs = [p for p in OUT.iterdir() if "1.25.0" in p.name]
    (OUT / "update-manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(
                    [
                        *paths,
                        Path(__file__),
                        ROOT / "oncoref/cta_selection_policy.py",
                        ROOT / "oncoref/data/normal-induction-evidence.csv",
                    ]
                ),
                "outputs": {p.name: sha256(p) for p in sorted(outputs)},
            },
            indent=2,
        )
        + "\n"
    )
    print(counts.to_string(index=False))
    print(categories[["Category", "Proteins"]].to_string(index=False))


if __name__ == "__main__":
    build()
