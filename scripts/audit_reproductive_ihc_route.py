#!/usr/bin/env python3
"""Measure an alternative reproductive IHC route without changing catalog policy."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_screen, curation, reference_data, target_selection  # noqa: E402
from oncoref.normal_tissues import resolve_normal_tissue  # noqa: E402

OUT = ROOT / "docs/audits/antigen-evidence-review-20261117"


def assess_group(rows, members, full_panel):
    """Any reliable reproductive High; every reported other cell must be negative.

    Missing/uncertain positive controls never qualify. Missing cell/tissue levels
    are not negatives; an absent member locus cannot borrow another's negatives.
    Complete coverage of every tissue in the source is reported separately.
    """
    high = rows.reproductive & rows.Level.eq("High") & rows.Reliability.isin(cta_screen.RELIABLE)
    other = rows[~rows.reproductive]
    tested = other[other.Tissue.ne("") & other.Tissue.ne("N/A")]
    measured_members = set(tested.canonical_gene_id)
    negative = (
        bool(len(other))
        and set(members) <= measured_members
        and other.Level.eq("Not detected").all()
        and len(tested) == len(other)
        and other.Reliability.isin(cta_screen.RELIABLE).all()
    )
    covered = {gene: set(group.Tissue) for gene, group in tested.groupby("canonical_gene_id")}
    complete = all(full_panel <= covered.get(gene, set()) for gene in members)
    return {
        "high_reproductive": bool(high.any()),
        "high_reliabilities": "; ".join(sorted(rows.loc[high, "Reliability"].unique())),
        "high_locations": "; ".join(
            sorted(set(rows.loc[high, "Tissue"] + ": " + rows.loc[high, "Cell type"]))
        ),
        "nonreproductive_tissues_measured": len(set(tested.Tissue)),
        "nonreproductive_rows_measured": len(tested),
        "nonreproductive_tissues": json.dumps(sorted(set(tested.Tissue))),
        "all_members_measured": set(members) <= measured_members,
        "observed_panel_negative": bool(negative),
        "full_source_tissue_panel_complete": complete,
        "missing_source_tissues": json.dumps(
            {gene: sorted(full_panel - covered.get(gene, set())) for gene in members},
            sort_keys=True,
        ),
        "ihc_route": bool(high.any() and negative),
    }


def main():
    declared = curation.manifest()["reference_inputs"]["host_protein"]
    source = reference_data.local_path(declared["source"], declared["version"])
    pin = reference_data.REFERENCE_SOURCES[declared["source"]]["pins"][declared["version"]]
    if curation.digest(source) != pin["sha256"] or source.stat().st_size != pin["bytes"]:
        raise ValueError("HPA IHC source does not match its immutable pin")
    mapping = (
        ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz"
    )
    catalog = target_selection.results()
    members = {r.entity_id: json.loads(r.member_gene_ids) for r in catalog.itertuples()}
    gene_entity = {gene: entity for entity, genes in members.items() for gene in genes}
    lookup = pd.read_csv(mapping).set_index("source_gene_id").canonical_gene_id
    data = pd.read_csv(source, sep="\t", keep_default_na=False)
    data["canonical_gene_id"] = data.Gene.map(lookup)
    data["entity_id"] = data.canonical_gene_id.map(gene_entity)
    locations = {}
    for tissue, cell in data[["Tissue", "Cell type"]].drop_duplicates().itertuples(index=False):
        annotation = resolve_normal_tissue(
            tissue,
            cell_type=cell,
            species="human",
            source_id="hpa_normal_tissue",
            source_version="v23",
        )
        locations[tissue, cell] = "reproductive" in annotation.scopes
    data["reproductive"] = [locations[t, c] for t, c in zip(data.Tissue, data["Cell type"])]
    full_panel = set(data.loc[~data.reproductive, "Tissue"]) - {"", "N/A"}
    high_ids = set(
        data.loc[
            data.reproductive & data.Level.eq("High") & data.Reliability.isin(cta_screen.RELIABLE),
            "entity_id",
        ].dropna()
    )
    rows = [
        {"entity_id": entity, **assess_group(group, members[entity], full_panel)}
        for entity, group in data[data.entity_id.isin(high_ids)].groupby("entity_id", sort=True)
    ]
    result = catalog.merge(pd.DataFrame(rows), on="entity_id", validate="one_to_one")
    result["shared_eligible"] = result.category.ne("Unsafe") & result[
        list(target_selection.SHARED_CHECKS)
    ].all(axis=1)
    result["additional_ra"] = (
        result.ihc_route & result.shared_eligible & ~result.is_reproductive_antigen
    )
    result["additional_cta"] = result.additional_ra & result.cancer_evidence.eq("qualifies")
    result["additional_union"] = result.additional_ra & ~result.is_target_candidate
    result["failed_shared_checks"] = result.apply(
        lambda row: "; ".join(k for k in target_selection.SHARED_CHECKS if not row[k]), axis=1
    )
    fields = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "classification",
        "reproductive_fold",
        "tumor_fold",
        "rare_cell_status",
        "failed_shared_checks",
        *[key for key in rows[0] if key != "entity_id"],
        "shared_eligible",
        "additional_ra",
        "additional_cta",
        "additional_union",
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    result[fields].to_csv(
        OUT / "ihc-counterfactual.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
    )
    summary = {
        "catalog_policy": target_selection.POLICY_ID,
        "catalog_sha256": target_selection.metadata()["data_sha256"],
        "source": "hpa_normal_tissue:v23",
        "source_sha256": curation.digest(source),
        "identity_mapping_sha256": curation.digest(mapping),
        "implementation_sha256": curation.digest(Path(__file__)),
        "full_source_nonreproductive_tissues": sorted(full_panel),
        "definition": "Reliable reproductive High plus reliable Not detected in every reported non-reproductive cell/tissue row, with every member locus represented. All shared gates retained. Missing tissues remain missing.",
        "limitations": [
            "Counterfactual only; no catalog membership or tenfold policy is changed.",
            "High/Not detected is not a measured tenfold protein ratio or clinical safety margin.",
            "HPA consensus rows do not expose matched antibody/run/donor pairs or prove locus-specific staining.",
            "All reported negative tissues is not all tissues in the body. Complete source-panel coverage is reported separately.",
        ],
        "baseline": {
            name: int(catalog[column].sum())
            for name, column in [
                ("CTA", "is_cta"),
                ("RA", "is_reproductive_antigen"),
                ("union", "is_target_candidate"),
            ]
        },
        "additional": {
            name: int(result[column].sum())
            for name, column in [
                ("CTA", "additional_cta"),
                ("RA", "additional_ra"),
                ("union", "additional_union"),
            ]
        },
        "additional_cta_labels": result.loc[result.additional_cta, "protein_group_label"].tolist(),
        "additional_ctas_with_complete_source_panel": int(
            (result.additional_cta & result.full_source_tissue_panel_complete).sum()
        ),
    }
    (OUT / "ihc-counterfactual.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {key: summary[key] for key in ("baseline", "additional", "additional_cta_labels")}
        )
    )


if __name__ == "__main__":
    main()
