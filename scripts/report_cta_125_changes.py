#!/usr/bin/env python3
"""Compare the frozen 1.24 selection with the 1.25 normal-pass/IHC update."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from cta_report_common import input_hashes, sha256, write_csv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-neither-20261030"
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
ID = "protein_equivalence_group_id"


def table(frame):
    rows = [frame.columns.tolist(), ["---"] * len(frame.columns)]
    rows.extend(frame.fillna("unmeasured").itertuples(index=False, name=None))
    return "\n".join("| " + " | ".join(str(v).replace("|", ";") for v in r) + " |" for r in rows)


def build():
    inputs = [
        OUT / "baseline-1.24.0.csv.gz",
        CURRENT / "protein-group-classifications.csv.gz",
        CURRENT / "reproductive-protein-evidence.csv.gz",
        CURRENT / "tumor-restricted-observations.csv.gz",
    ]
    before, after, ihc, tumors = [pd.read_csv(p, low_memory=False) for p in inputs]
    before = before.set_index(ID)
    after = after.set_index(ID).reindex(before.index)
    old, new = before.category.fillna("Unassigned"), after.category.fillna("Unassigned")
    assert set(before.index) == set(after.index)
    assert after.loc[old.eq("CTA"), "category"].eq("CTA").all()
    assert set(before.index[old.eq("Unsafe")]) == set(after.index[new.eq("Unsafe")])
    changes = after.loc[old.ne(new)].copy()
    changes.insert(0, "category_1_24", old.loc[changes.index])
    write_csv(changes.reset_index(), OUT / "changes-1.25.0.csv.gz")
    transitions = pd.crosstab(old.rename("1.24"), new.rename("1.25")).reset_index()
    write_csv(transitions, OUT / "transitions-1.25.0.csv")
    columns = [
        "protein_group_label",
        "member_gene_ids",
        "reproductive_ntpm",
        "developmental_mean_tpm",
        "oncoref_maximum_tpm",
        "hpa_tumor_group_lower_bound_ptpm",
        "hpa_diffuseness",
        "adult_profile_status",
    ]
    restored = after.loc[new.eq("CTA") & ~old.eq("CTA"), columns].sort_values("protein_group_label")
    positive = ihc[ihc.qualifying_reproductive_positive]
    restored["reproductive_IHC"] = restored.member_gene_ids.map(
        lambda ids: "; ".join(
            sorted(
                {
                    f"{r.gene_id}: {r.tissue}/{r.cell_type}: {r.Level} ({r.Reliability})"
                    for r in positive[positive.gene_id.isin(json.loads(ids))].itertuples(
                        index=False
                    )
                }
            )
        )
    )
    write_csv(restored.reset_index(), OUT / "restored-CTA-1.25.0.csv")
    evidence = positive[
        positive.gene_id.isin({g for ids in restored.member_gene_ids for g in json.loads(ids)})
    ]
    write_csv(evidence, OUT / "restored-reproductive-IHC-1.25.0.csv")
    leads = after[after.tumor_restricted_lead].copy().sort_values("protein_group_label")
    lead_obs = tumors[tumors[ID].isin(leads.index)].copy()
    assert lead_obs.groupby(ID).sample_key.nunique().eq(2).all()
    lead_obs.insert(1, "protein", lead_obs[ID].map(after.protein_group_label))
    write_csv(lead_obs, OUT / "tumor-restricted-lead-observations-1.25.0.csv")
    leads["coverage_review"] = "Relevant normal cells and physiological states require review"
    leads.loc[leads.protein_group_label.str.match(r"OR\d"), "coverage_review"] = (
        "Olfactory sensory epithelium; locus-specific expression/sequence attribution"
    )
    special = {
        "SCN10A": "Human peripheral sensory neurons / dorsal root ganglia",
        "GFRAL": "Area postrema and nucleus tractus solitarius",
        "OTOR": "Specialized inner-ear cells",
        "OC90": "Specialized inner-ear cells",
        "CNGA2": "Olfactory sensory neurons",
        "IFNA4": "Induced antiviral immune states",
        "IFNW1": "Induced antiviral immune states",
    }
    for label, gap in special.items():
        leads.loc[leads.protein_group_label.eq(label), "coverage_review"] = gap
    write_csv(
        leads.reset_index()[
            [
                ID,
                "protein_group_label",
                "tumor_samples_ge10",
                "normal_maximum_floor_ratio",
                "developmental_all_lineage_max_tpm",
                "coverage_review",
            ]
        ],
        OUT / "tumor-restricted-leads-1.25.0.csv",
    )
    held = after[old.eq("Neither") & new.eq("Unassigned")]
    write_csv(held.reset_index(), OUT / "previous-neither-normal-holds-1.25.0.csv")
    counts = pd.DataFrame(
        {"Category": ["CTA", "Tumor-restricted", "Neither", "Unsafe", "Unassigned"]}
    )
    counts["1.24.0"] = counts.Category.map(old.value_counts()).fillna(0).astype(int)
    counts["1.25.0"] = counts.Category.map(new.value_counts()).fillna(0).astype(int)
    restored_display = restored[
        [
            "protein_group_label",
            "reproductive_ntpm",
            "oncoref_maximum_tpm",
            "hpa_tumor_group_lower_bound_ptpm",
        ]
    ].round(3)
    restored_display.columns = [
        "Protein",
        "Reproductive nTPM",
        "OncoRef tumor TPM",
        "HPA tumor pTPM",
    ]
    report = f"""# Neither, reproductive IHC and Tumor-restricted candidates

This update restores **{len(restored)} CTA proteins**, expands Neither to **{int(new.eq("Neither").sum())}**, and preserves every previous CTA and Unsafe assignment. Identical selected reference protein sequences remain collapsed; these are protein counts, not locus counts.

{table(counts)}

[Full current report and figures](../cta-source-walkthrough-20261028/index.md) · [1.25 logic slides](../cta-source-walkthrough-20261028/cta-selection-1.25.0.pdf) · [Original 20-Neither review](index.md) · [Every category change](changes-1.25.0.csv.gz) · [Transition table](transitions-1.25.0.csv)

## The corrected rules

1. **Normal risk comes first.** Existing Vital/Protected/clinical exclusions establish Unsafe. Missing required normal coverage or unresolved normal evidence remains unassigned. Manageable ceiling excess is a normal hold.
2. **CTA identity requires cancer and Reproductive evidence.** The existing tumor rule stays ≥2 native TPM or HPA pTPM in one sample. Reproductive RNA ≥5 nTPM, supported developmental RNA ≥5 TPM, or reliable positive reproductive IHC can establish association. The existing profile-concentration assessment still applies where RNA is comparable.
3. **Tumor-restricted is a separate candidate subtype.** It needs ≥10 native TPM in two distinct tumor samples, or reviewed reliable tumor protein evidence. It also needs low measured normal/developmental expression, no detected normal protein, resolved attribution, all normal gates passing and a review addressing relevant rare normal cells and physiological states.
4. **Neither is the fallback after a complete normal pass.** Cancer or Reproductive evidence may be high, low or missing. Failure to qualify as a candidate does not establish Unsafe.

All 725 Neither proteins remain below 30 in the measured non-permitted normal TPM/nTPM summaries: the maxima are 16.0 nTPM and 15.3163 TPM, kept as separate assays. Missing optional measurements remain missing.

Neither describes a research screen. It is not evidence that targeting a protein is safe. Existing exact-protein background/k-mer APIs keep their independent policy and are not silently rewritten by this selection update.

## Why the 52 IHC cases were overlooked

The newer selection-policy rewrite kept reproductive IHC as an antibody control but omitted its independent positive-association route. Admission used only reproductive or developmental RNA ≥5. The earlier screen had allowed protein evidence below that RNA floor. The migration tests did not cover that route.

The restored rule requires human HPA v23, a mapped reproductive compartment, Low/Medium/High staining with Approved/Supported/Enhanced reliability, positive evidence for every member locus, and no unresolved antibody concern. It cannot override normal exclusions, cancer admission or a measured concentration failure. A zero reproductive RNA reference leaves RNA concentration unassessed; an IHC category does not supply a fabricated RNA abundance.

**51 of the 52 become CTA. CDY1/B remains Neither because it lacks qualifying cancer expression.** All 522 earlier CTAs, including MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B, remain CTA. The full HPA observations retain cell type, staining, reliability, mapped anatomy and source version.

[Restored proteins and measurements](restored-CTA-1.25.0.csv) · [Individual positive IHC rows](restored-reproductive-IHC-1.25.0.csv)

{table(restored_display)}

## Tumor-restricted leads need normal-context evidence

**{len(leads)} proteins meet the numerical profile, but zero currently have an accepted normal-restriction review.** They remain Neither with a Tumor-restricted lead flag. The normal-low convention is <2 native TPM/nTPM, with a separate nominal CP10K floor; it is not a detection limit or a conversion between assays. Developmental coverage is required for this restrictive subtype. Missing optional sources remain explicitly missing.

Forty leads are olfactory-receptor genes. A general bulk normal panel cannot establish absence from rare sensory cells. [Durante et al., *Single cell analysis of olfactory neurogenesis and differentiation in adult humans* (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7065961/) demonstrates this coverage problem at the receptor-family level; it is not a gene-by-gene verdict on all forty leads. [Human dorsal-root-ganglion single-nucleus profiling](https://pmc.ncbi.nlm.nih.gov/articles/PMC8626086/) supplies relevant normal-cell evidence for SCN10A. [The GDF15–GFRAL receptor study (2017)](https://pubmed.ncbi.nlm.nih.gov/28846097/) reports GFRAL in the human area postrema and nucleus tractus solitarius. These findings motivate coverage review without converting missing bulk measurements into either absence or clinical toxicity.

[All 66 leads and coverage questions](tumor-restricted-leads-1.25.0.csv) · [Their two qualifying tumor observations](tumor-restricted-lead-observations-1.25.0.csv)

The tumor reader uses the existing checksum-pinned source-matrix cache and QC-selected samples. It collapses duplicate source views and TCGA participant/aliquot identities. Both observations must resolve to the protein's loci, and all identical-sequence member loci must be measured before summing within a sample. Counts stop at two and are lower bounds, not prevalence. HPA pTPM is not substituted for native TPM. The reliable-tumor-protein route is implemented, but its review ledger currently has no accepted entries; aggregate HPA cancer-IHC counts alone are insufficient reagent evidence.

USP17L3 and FAM90A9 are now Neither. Their low measured normals no longer become unassigned merely because cancer RNA exceeds the CTA cancer floor. Neither has two qualifying ≥10 TPM observations in this screen.

## Seven old Neither assignments now retain their normal holds

**{", ".join(held.protein_group_label)}** previously entered a narrow multiple-Manageable-tissue profile despite a normal-screen hold. The fallback now requires every normal gate to pass. They remain unassigned with a normal-evidence hold, rather than gaining an unsupported normal-pass label. Unsafe membership is unchanged.

[Measurements and all holds for these seven proteins](previous-neither-normal-holds-1.25.0.csv)

## API and reproducibility

```python
from oncoref import cta_selection

cta_selection.results(category="CTA")
cta_selection.results(category="Tumor-restricted")  # empty until a review qualifies
cta_selection.candidates()  # union of both candidate subtypes
cta_selection.results(category="Neither")
cta_selection.audit().explain("SAGE1")
```

Policy: `oncoref.cta-selection.2026-10-30.1`. Run `scripts/build_tumor_restricted_support.py` with the pinned matrices cached, then `scripts/regenerate_cta_walkthrough.py --verify-replay` and `scripts/report_cta_125_changes.py`. `--download` on the tumor support builder fills missing cache artifacts through the existing uniform source API. Full and cached decision replay are byte-identical. The [current rules](../../cta-selection.md) document thresholds, exact species scope and clinical-credit limitations.
"""
    (OUT / "update-1.25.0.md").write_text(report)
    outputs = [p for p in OUT.iterdir() if "1.25.0" in p.name]
    (OUT / "update-manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes([*inputs, Path(__file__)]),
                "outputs": {p.name: sha256(p) for p in sorted(outputs)},
            },
            indent=2,
        )
        + "\n"
    )
    print(counts.to_string(index=False))


if __name__ == "__main__":
    build()
