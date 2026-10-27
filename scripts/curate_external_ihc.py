#!/usr/bin/env python3
"""Build structured external staining evidence and a reviewable evidence report.

Supplier extraction, reviewed publication findings and policy eligibility remain
separate columns. Counts from adjacent tumors, RNA assays or flow cytometry are
never relabeled as healthy normal IHC donor counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cta_report_common import (  # noqa: E402
    input_hashes,
    sha256,
    write_csv,
)

from oncoref.external_ihc import assess_evidence  # noqa: E402
from oncoref.gene_matching import gene_identity_catalog  # noqa: E402


def publication_observations() -> pd.DataFrame:
    """Explicit findings checked against the cited primary texts and figure tables."""
    symbols = gene_identity_catalog().set_index("canonical_symbol").canonical_gene_id
    rows: list[dict] = []

    def add(
        symbol: str,
        tissue: str,
        cells: str,
        antibody: str,
        result: str,
        url: str,
        locator: str,
        **fields: object,
    ) -> None:
        rows.append(
            {
                "gene_id": symbols[symbol],
                "tissue": tissue,
                "cell_type": cells,
                "antibody_id": antibody,
                "reagent_lineage": antibody,
                "result": result,
                "source_url": url,
                "source_locator": locator,
                "source_kind": "peer_reviewed",
                "taxonomy_id": 9606,
                "stage": "adult",
                "specimen_state": "normal",
                "assay": "IHC",
                "review_status": "reviewed",
                "positive_control": False,
                "specificity_validation": "not_established_in_ihc",
                "validation_locator": "",
                "exact_cell_coverage": cells != "unspecified",
                "normal_donors": None,
                "independent_of_hpa": True,
                "notes": "",
                "intensity": "unspecified",
                "localization": "unspecified",
                "target_scope": "protein_or_recognized_isoform",
                **fields,
            }
        )

    add(
        "TBXT",
        "lung",
        "unspecified",
        "54-1",
        "not_detected",
        "https://www.oncotarget.com/article/3086/text/",
        "Results; Figure 4C",
        positive_control=True,
        specificity_validation="genetic_western_blot",
        validation_locator="Figures 1C, 2A, 4; WB knockdown, TBX19 ELISA and IHC isotype control",
        notes="Normal lung negative; normal donor count and rare-cell coverage not enumerated. WB knockdown is not IHC knockdown.",
    )
    add(
        "TBXT",
        "lung",
        "macrophages",
        "HPA003322",
        "reported_cross_reactivity",
        "https://www.oncotarget.com/article/3086/text/",
        "Results preceding Figure 4",
        independent_of_hpa=False,
        notes="Authors describe Prestige/HPA macrophage signal as background. This is an attribution, not a new quantified negative panel.",
    )
    add(
        "TBXT",
        "intervertebral disc",
        "nucleus pulposus cells",
        "ab57480",
        "detected",
        "https://www.nature.com/articles/s41598-017-01567-w",
        "Figures 2–5; Methods / immunohistochemistry",
        specimen_state="surgical_disc_degeneration",
        specimens_assayed=56,
        normal_donors=None,
        positive_control=True,
        validation_locator="Isotype controls; companion RNA N=60, flow N=5",
        notes="Adult surgical samples across degeneration grades; 56 IHC samples are not 56 healthy donors. A separate adolescent image is not pooled into the adult evidence.",
    )
    add(
        "TBXT",
        "intervertebral disc",
        "nucleus pulposus cells",
        "IC20851A",
        "detected",
        "https://www.nature.com/articles/s41598-017-01567-w",
        "Figure 6; flow cytometry",
        assay="flow_cytometry",
        specimen_state="surgical_disc_degeneration",
        specimens_assayed=5,
        positive_fraction_mean=76.2,
        positive_fraction_low=59,
        positive_fraction_high=97,
        positive_fraction_sd=16.3,
        notes="Percent positive cells in five freshly isolated surgical NP samples; retained as flow, not IHC.",
    )
    add(
        "CT83",
        "stomach",
        "epithelial cells",
        "CL4762",
        "detected",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC7953676/",
        "Tables 2–3; Methods / immunohistochemistry",
        specimen_state="tumor_adjacent",
        specimens_assayed=84,
        positive_specimens=76,
        score_name="H-score",
        score_median=50,
        score_low=0,
        score_high=150,
        score_scale_min=0,
        score_scale_max=300,
        localization="cytoplasmic",
        intensity="low;moderate",
        validation_locator="Two pathologists; specificity validation not established",
        notes="Adjacent gastric tissue from cancer resections is not a healthy normal panel. Conflicting evidence retained for a named CTA control; not a gene-specific waiver.",
    )
    add(
        "ADAM30",
        "cerebral cortex",
        "neurons",
        "Interchim-custom-ADAM30",
        "detected",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC4972530/",
        "Figure 1C, S4; Methods / human brain immunohistochemistry",
        specimen_state="autopsy_controls_and_alzheimer",
        specificity_validation="preimmune_control_only",
        validation_locator="Custom Interchim polyclonal and preimmune serum; IHC donor count unreported",
        notes="Cortical neuronal staining. RNA cases=51 and controls=42 are a separate sample set, not IHC counts. GeneTex reagent from other experiments is not assigned to this IHC.",
    )
    add(
        "CDX4",
        "lung",
        "unspecified",
        "PA5-112380",
        "detected",
        "https://www.thermofisher.com/antibody/product/CDX4-Antibody-Polyclonal/PA5-112380",
        "Human lung IHC product image",
        source_kind="manufacturer_image",
        independent_of_hpa=False,
        notes="Staining documented; reagent origin, specificity and donor count unresolved. Potential reseller copies are not independent replicates.",
    )
    add(
        "RNASE11",
        "testis",
        "seminiferous and Leydig cells",
        "PA5-51991",
        "detected",
        "https://www.thermofisher.com/antibody/product/RNase11-Antibody-Polyclonal/PA5-51991",
        "Human testis IHC image",
        source_kind="manufacturer_image",
        intensity="moderate",
        localization="cytoplasmic",
        independent_of_hpa=False,
        notes="Reproductive positive context; Atlas-origin independence not established. Does not provide adult nonreproductive negative coverage.",
    )
    add(
        "GOT1L1",
        "duodenum",
        "glandular cells",
        "PA5-55977",
        "detected",
        "https://www.thermofisher.com/antibody/product/GOT1L1-Antibody-Polyclonal/PA5-55977",
        "Human duodenum IHC image",
        source_kind="manufacturer_image",
        intensity="high",
        localization="cytoplasmic",
        independent_of_hpa=False,
        notes="Potential Atlas-origin reagent; source independence and specificity remain unresolved.",
    )
    add(
        "FOXR2",
        "placenta",
        "trophoblast cells",
        "2C1",
        "detected",
        "https://www.scbt.com/p/foxr2-antibody-2c1",
        "Normal placenta IHC image",
        stage="placental",
        source_kind="manufacturer_image",
        localization="nuclear",
        notes="Same clone marketed by other suppliers is one reagent lineage; placental staining does not establish adult absence.",
    )
    result = pd.DataFrame(rows)
    result["evidence_id"] = [
        "curated:"
        + hashlib.sha256((r.source_url + r.source_locator + r.antibody_id).encode()).hexdigest()[
            :20
        ]
        for r in result.itertuples()
    ]
    return result


def build(output: Path) -> None:
    """Publish parsed categories, reviewed findings and explicit decision limitations."""
    observations = pd.read_csv(output / "external-panel-observations.csv.gz", keep_default_na=False)
    inventory = pd.read_csv(output / "external-panel-inventory.csv")
    d = observations[observations.identity_status.eq("matched")].copy()
    d = d.rename(columns={"canonical_gene_id": "gene_id", "source_class": "source_kind"})
    d["antibody_id"] = d.reagent_title.str.extract(r"\(([^()]+)\)")[0].fillna("")
    d["reagent_lineage"] = d.antibody_id
    d["stage"] = "adult"
    d.loc[d.tissue.str.contains("placent|amnion|chorion", case=False), "stage"] = "placental"
    d.loc[d.source_statement.str.contains("fetal|children|pediatric", case=False), "stage"] = (
        "mixed_or_developmental"
    )
    d["specimen_state"] = "normal"
    d["taxonomy_id"] = 9606
    d["positive_control"] = False
    d["specificity_validation"] = "supplier_claim_unadjudicated"
    d["validation_locator"] = ""
    d["normal_donors"] = None
    d["independent_of_hpa"] = False
    d["notes"] = (
        "Reported manufacturer observation. Exact cell coverage, donor counts and independent IHC validation not established."
    )
    d["evidence_id"] = [
        "catalog:"
        + hashlib.sha256(
            (r.source_url + r.source_locator + r.source_statement).encode()
        ).hexdigest()[:20]
        for r in d.itertuples()
    ]
    d = d.drop_duplicates("evidence_id")
    reviewed = publication_observations()
    combined = pd.concat([d, reviewed], ignore_index=True)
    assessed = assess_evidence(combined)
    write_csv(combined, ROOT / "oncoref/data/external-ihc-evidence.csv.gz")
    write_csv(assessed, output / "assessed-evidence.csv.gz")
    write_csv(
        assessed[assessed.review_status.eq("reviewed")],
        output / "reviewed-publication-and-image-evidence.csv",
    )
    gaps = inventory[inventory["Antibody concern"].eq(True) | inventory.ihc_missing.eq(True)]
    summary = {
        "observations": len(assessed),
        "genes": assessed.gene_id.nunique(),
        "reviewed_observations": len(reviewed),
        "result_counts": assessed.result.value_counts().to_dict(),
        "context_validated_observations": int(assessed.eligible.sum()),
        "catalog_genes_with_missing_or_uncertain_ihc": sorted(gaps.canonical_symbol.unique()),
    }
    (output / "evidence-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "index.md").write_text(f"""# Structured external staining evidence

{summary["observations"]:,} observations across {summary["genes"]} matched genes are represented as
**detected, not detected, mixed, reported cross-reactivity or not assessed**.
Intensity, localization, explicitly reported percentages, H-scores, specimen counts,
age/stage and specimen state have separate fields. Missing numbers remain missing.
Source statements and table/figure locations accompany each observation.

The 201-product MSVA collection supplies additional panels for **{", ".join(summary["catalog_genes_with_missing_or_uncertain_ihc"])}**
among proteins with missing or uncertain HPA evidence. Paragraph observations have
categorical results even when anatomical mapping still requires review.
{len(reviewed)} additional findings were reviewed against publication text or product images.

## How an external result can supersede an Uncertain HPA observation

Use the same rule for every gene: reviewed human adult normal-tissue IHC, independent
identified reagent, a positive control, specificity validation in IHC, the exact
tissue and cell compartment, and at least two normal donors for a negative result.
Two donors is a provisional review-policy minimum, not a validated absence cutoff.
Documented positives require at least one normal donor and the same context and
specificity checks. Validated positive/negative disagreement remains a conflict.
Western-blot specificity does not automatically validate tissue IHC. Supplier
panels do not inherit HPA's Approved, Supported or Enhanced grades.

**{summary["context_validated_observations"]} observations currently meet every requirement.**
Consequently this evidence collection produces no automatic negative-IHC overrides.
Each failure is machine-readable in `ineligible_reasons`; collection is not clearance.

## Findings that matter to the current review

- TBXT: 54-1 reports negative normal lung, while HMV321 reports negative organ panels.
  Missing donor/cell-level documentation prevents a context-specific override.
  Adult disc staining remains recorded separately, including surgical specimen state
  and the companion flow assay (5 specimens; mean 76.2% positive cells).
- CT83: CL4762 staining was reported in 76/84 tumor-adjacent gastric specimens;
  median H-score 50, range 0–150 on a 0–300 scale. This is a review flag for a named
  control, not a healthy-donor panel and not grounds for pretending the signal is zero.
- ADAM30: cortical neuronal staining with a custom Interchim reagent is a positive
  finding; the separate RNA donor counts are not borrowed for IHC.
- CDX4 and GOT1L1 have additional positive tissue images, with specificity and reagent
  independence unresolved. RNASE11 and FOXR2 findings are reproductive/placental.

The detailed rows link each finding to its original source.
[Reviewed findings](reviewed-publication-and-image-evidence.csv),
[all assessed evidence](assessed-evidence.csv.gz),
[catalog and gene identity mapping](external-panel-inventory.csv),
[original categorical extraction](external-panel-observations.csv.gz),
[source receipts](source-receipts.json), [summary](evidence-summary.json).

Reproduce with `python scripts/collect_external_ihc.py` followed by
`python scripts/curate_external_ihc.py`. These are research evidence annotations;
released classifications and sequence references are unchanged.
""")
    (output / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(
                    [
                        Path(__file__),
                        ROOT / "oncoref/external_ihc.py",
                        ROOT / "oncoref/ihc_observations.py",
                        output / "source-receipts.json",
                    ]
                ),
                "outputs": {
                    p.name: sha256(p)
                    for p in output.iterdir()
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/audits/external-ihc-20261027")
    build(parser.parse_args().output)
