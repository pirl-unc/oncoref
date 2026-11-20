"""Re-evaluate normal gates whose recorded peak can change under an organ revision."""

import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd
from cta_report_common import sha256, write_csv
from replay_cta_selection import ALTERNATIVE, BASE, CROSSWALK, FROZEN, OUTPUT, ROOT, replay_contexts

from oncoref import cta_selection_policy as policy
from oncoref import reference_data, tissue_tiers

OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"


class GateCollector:
    """Reduce supported RNA and conservative bounds across all organs for a selected gene."""

    def __init__(self, facts):
        self.credit = facts.set_index("representative_gene_id")[
            ["clinical_rna_credit", "clinical_product_ids"]
        ]
        self.parts = []

    def collect(self, frame, assay):
        if assay not in {"Normal bulk", "HPA tissue clusters"}:
            return
        frame = frame.loc[frame.gene_id.isin(self.credit.index)].copy()
        if frame.empty:
            return
        for column in self.credit:
            frame[column] = frame.gene_id.map(self.credit[column])
        d = policy.assess_contexts(frame)
        if "value_upper" in d:
            d["candidate_exclusion"] |= d.measured_members.eq(d.member_count) & d.value_upper.gt(
                d.candidate_ceiling
            )
        if "upper_supported_value" in d:
            d["candidate_exclusion"] |= d.upper_supported_value.gt(d.candidate_ceiling)
        self.parts.append(
            d.groupby("gene_id")[["candidate_exclusion", "unsafe_evidence"]]
            .any()
            .assign(assay=assay)
        )

    def finish(self):
        return pd.concat(self.parts).groupby(["gene_id", "assay"]).any().reset_index()


def main():
    facts = pd.read_csv(OUTPUT / "protein-group-classifications.csv.gz", low_memory=False)
    peaks = pd.read_csv(OUTPUT / "normal-peak-evidence.csv.gz", usecols=["gene_id", "tissue"])
    revisions = tissue_tiers.assessments().dropna(subset=["previous_tier"])
    if not (
        revisions.level.eq("bulk")
        & revisions.previous_tier.eq("protected")
        & revisions.tier.eq("disposable")
    ).all():
        raise ValueError("This replay covers whole-organ Protected to Disposable revisions")
    organs = {t for organs in revisions.tissues for t in organs}
    # Separate worst positive/negative peaks are recorded per assay. An unchanged
    # positive peak independently supports the gate. Replay all multi-locus groups
    # as their conservative upper bounds need not coincide with a measured peak.
    candidates = facts.representative_gene_id.isin(
        peaks.loc[peaks.tissue.isin(organs), "gene_id"]
    ) | facts.member_count.gt(1)
    selected = facts[candidates].copy()
    with tempfile.TemporaryDirectory() as directory:
        stage = Path(directory)
        for path in OUTPUT.glob("*-accepted-identities.csv.gz"):
            shutil.copyfile(path, stage / path.name)
        rna = replay_contexts(selected, stage, collector=GateCollector(selected))
    # Every positive IHC record is retained, including positives outside the
    # revised organ. Preserve the source's eligibility and reagent adjudications.
    ihc = pd.read_csv(OUTPUT / "hpa-protein-evidence.csv.gz", low_memory=False)
    affected = ihc.tissue.isin(organs)
    tier = tissue_tiers.context_tiers(ihc.loc[affected])
    ihc.loc[affected, "candidate_exclusion"] &= tier.isin(policy.TIER_BASES)
    ihc.loc[affected, "unsafe_evidence"] &= tier.isin(["vital", "protected"])
    protein = (
        ihc.groupby("representative_gene_id")[["candidate_exclusion", "unsafe_evidence"]]
        .any()
        .rename_axis("gene_id")
        .reset_index()
        .assign(assay="HPA normal protein")
    )
    rows = pd.concat([rna, protein], ignore_index=True)
    base = facts.set_index("representative_gene_id")
    changed = []
    for assay, group in rows.groupby("assay"):
        group = group.set_index("gene_id")
        original = base.loc[group.index, ["hold:" + assay, "unsafe:" + assay]]
        original.columns = ["candidate_exclusion", "unsafe_evidence"]
        diff = group[original.columns].ne(original).any(axis=1)
        changed.append(group[diff].reset_index())
    write_csv(
        pd.concat(changed).sort_values(["gene_id", "assay"]), OUT / "revised-tissue-gates.csv"
    )
    inputs = [
        OUTPUT / f
        for f in [
            "protein-group-classifications.csv.gz",
            "normal-peak-evidence.csv.gz",
            "hpa-protein-evidence.csv.gz",
        ]
    ]
    inputs.extend(OUTPUT.glob("*-accepted-identities.csv.gz"))
    inputs.extend(
        [
            CROSSWALK,
            FROZEN / "bulk-contexts.csv.gz",
            ALTERNATIVE / "normal-bulk.csv.gz",
            BASE / "group-normal/cluster-contexts.csv.gz",
            ROOT / "oncoref/data/normal-tissue-tier-updates.yaml",
        ]
    )
    receipt = {str(p.relative_to(ROOT)): sha256(p) for p in inputs}
    for source in ["hpa_rna_consensus", "hpa_single_cell_tissue", "hpa_normal_tissue"]:
        receipt[f"reference:{source}:v23"] = reference_data.provenance(source, verify_content=True)[
            "sha256"
        ]
    (OUT / "revised-tissue-gate-inputs.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
