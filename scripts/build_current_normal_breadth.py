#!/usr/bin/env python3
"""Compare named historical/current bulk panels without changing CTA admission."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_normal_expression_breadth import digest, write  # noqa: E402

from oncoref import cta_policy, cta_selection, hpa, reference_data  # noqa: E402
from oncoref import normal_expression_breadth as breadth  # noqa: E402

OUTPUT = ROOT / "docs/audits/normal-breadth-20261104"


def build():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    policy = breadth.policy()
    records = cta_selection.results()
    wanted = sorted({gene for genes in records.gene_ids for gene in genes})
    source = reference_data.provenance("hpa_rna_consensus", "v23", verify_content=True)
    native = hpa._read_hpa("hpa_rna_consensus", "v23").rename(
        columns={"Gene": "gene_id", "Tissue": "compartment", "nTPM": "value"}
    )
    panel = breadth.anatomical_panel()
    grid = pd.MultiIndex.from_product(
        [wanted, panel.compartment], names=["gene_id", "compartment"]
    ).to_frame(index=False)
    profiles = grid.merge(
        native[["gene_id", "compartment", "value"]], how="left", validate="one_to_one"
    )
    profiles = profiles.assign(
        source="hpa_rna_consensus",
        version="v23",
        taxonomy_id=9606,
        unit="nTPM",
        cohort="HPA consensus",
        summary="consensus",
    )
    historical_ceilings = {
        tier: values[0] for tier, values in cta_policy.manifest()["rna_ceiling_ntpm"].items()
    }
    results, panels = [], []
    for scenario, resolver, ceilings in [
        ("historical", breadth.historical_panel, historical_ceilings),
        ("current_tiers_historical_ceilings", breadth.anatomical_panel, historical_ceilings),
        ("current", breadth.anatomical_panel, policy["ceilings"]),
    ]:
        for weight in [0.5, 1.0, 2.0]:
            p = resolver(brain_weight=weight)
            profiles["ceiling"] = profiles.compartment.map(p.set_index("compartment").category).map(
                ceilings
            )
            parts, scores = breadth.burden(profiles, p)
            results.append(scores.assign(scenario=scenario, brain_weight=weight))
            panels.append(p.assign(scenario=scenario, brain_weight=weight))
            if weight == 1 and scenario in {"historical", "current"}:
                write(parts, OUTPUT / f"{scenario}-profiles.csv.gz")
    scores = pd.concat(results, ignore_index=True)
    write(scores, OUTPUT / "breadth.csv.gz")
    write(pd.concat(panels, ignore_index=True), OUTPUT / "panels.csv")
    comparisons = scores[scores.brain_weight.eq(1)].pivot(
        index=["gene_id", "category"],
        columns="scenario",
        values=["lower", "upper", "observed_weight"],
    )
    comparisons.columns = ["_".join(x) for x in comparisons.columns]
    write(comparisons.reset_index(), OUTPUT / "before-after.csv.gz")
    # The donor source is independent of consensus; never manufacture matching donors.
    donors = pd.read_csv(ROOT / "oncoref/data/cta-donor-rna.csv").rename(
        columns={
            "Ensembl_Gene_ID": "gene_id",
            "tissue": "compartment",
            "source_version": "version",
            "rna_ntpm": "value",
        }
    )
    if not donors.unit.eq("HPA-v23-nTPM").all():
        raise ValueError("Unexpected donor units")
    donors = donors.assign(source="hpa_tissue_donors", unit="nTPM")
    summary = breadth.donor_summaries(donors)
    summary["ceiling"] = summary.compartment.map(panel.set_index("compartment").category).map(
        policy["ceilings"]
    )
    _, donor_scores = breadth.burden(summary, panel)
    write(summary, OUTPUT / "donor-summaries.csv")
    write(donor_scores, OUTPUT / "donor-breadth.csv")
    (OUTPUT / "policy.json").write_text(json.dumps(policy, indent=2) + "\n")
    text = f"""# Current normal-expression breadth

Descriptive policy `{policy["policy_id"]}` uses bulk tiers from
`{policy["hazard_registry_id"]}`. **No admission cutoff or membership change is introduced.**
The report covers {len(wanted):,} distinct member loci in the qualified reference
universe, each once; these are gene profiles, not protein-group counts.

Skin, spleen and thymus now resolve to **Protected**. The historical October panel
remains available as `historical_panel()`; `anatomical_panel()` uses the registry.
The non-retinal eye residual uses the bulk eye rule. Cell exceptions never lower
a whole-organ observation. Each organ has equal weight within its tier; brain
regions share one organ's weight. Missing compartments retain their fixed weight.

For each tier, **B = sum(w × x / (x + L))**. Reproductive expression does not enter
this calculation and cannot improve it. The separate reproductive concentration
score `mean(N / (N + R))` answers a different question.

Current descriptive baseline ceilings are {policy["ceilings"]} nTPM, without
gene-specific IHC or clinical-credit adjustments. These anchors are research
choices, not injury thresholds. The [before/after table](before-after.csv.gz)
separates historical tiers/ceilings, current tiers with historical ceilings,
and current tiers/ceilings. [All scores](breadth.csv.gz) include brain-weight
sensitivity 0.5, 1 and 2. [Panels](panels.csv) expose every weight and assessment ID.

The [current](current-profiles.csv.gz) and [historical](historical-profiles.csv.gz)
profiles retain native source, version, unit, measurement and missing coverage.
Lower/upper bounds describe missing anatomical weight, not confidence intervals.
An entirely missing profile remains [0, 1]. HPA consensus is a published aggregate,
not an independent donor cohort. Duplicate observations cannot increase weight.

[Donor summaries](donor-summaries.csv) and [donor burdens](donor-breadth.csv) remain
separate. The available four-donor CTAG2-heart series supports a median but not
the prespecified ten-donor q90; missing q90 is not zero. Marginal bootstrap
intervals do not provide a joint tissue interval or population-risk estimate.

Rebuild with `python scripts/build_current_normal_breadth.py`. The
[receipt](manifest.json) pins source content, ontology, ceilings, code and outputs.
The [October report](../normal-breadth-20261026/index.md) stays frozen and reproducible.
"""
    (OUTPUT / "index.md").write_text(text)
    inputs = [
        Path(__file__),
        ROOT / "oncoref/normal_expression_breadth.py",
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "oncoref/data/normal-tissue-hazards.yaml",
        ROOT / "oncoref/tissue_hazards.py",
        ROOT / "oncoref/data/normal-tissue-taxonomy.yaml",
        ROOT / "oncoref/data/cta-policy.json",
        ROOT / "oncoref/data/cta-donor-rna.csv",
        ROOT / "oncoref/data/cta-selection-human.csv.gz",
    ]
    manifest = {
        "policy": policy,
        "eligibility_changes": 0,
        "source": {k: source[k] for k in ["url", "sha256", "version"]},
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "outputs": {
            p.name: digest(p) for p in sorted(OUTPUT.iterdir()) if p.name != "manifest.json"
        },
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    build()
