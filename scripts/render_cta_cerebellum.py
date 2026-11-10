#!/usr/bin/env python3
"""Render the Sepp cerebellum CTA comparison from its saved audit tables."""

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/plac1-brain-20261110/cta-comparison"


def protein_link(row):
    ids = json.loads(row.member_gene_ids)
    return f"[{row.protein}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ids[0]})"


def table(frame):
    lines = [
        "| Protein | Current category | RNA detections / 16,833 nuclei | Detection % | Nuclei with ≥2 counts | Feature coverage |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in frame.itertuples():
        if row.coverage == "unmeasured":
            values = "— | — | —"
        else:
            values = f"{int(row.n_positive):,} | {row.pct_positive:.3f}% | {int(row.n_ge2):,}"
        lines.append(f"| {protein_link(row)} | {row.classification} | {values} | {row.coverage} |")
    return "\n".join(lines)


def main():
    adult = pd.read_csv(OUT / "adult-granule-all-proteins.csv")
    all_scopes = pd.read_csv(OUT / "protein-detection.csv.gz")
    cta = adult[adult.category.eq("CTA")]
    hits = cta[cta.pct_positive.ge(1)].sort_values("pct_positive", ascending=False).copy()
    for scope, label in [
        ("adult_granule_donor_#105", "donor_42"),
        ("adult_granule_donor_#138", "donor_52"),
    ]:
        lookup = all_scopes[all_scopes.scope.eq(scope)].set_index("entity_id")
        for key in ["n_positive", "n_cells", "pct_positive"]:
            hits[f"{label}_{key}"] = hits.entity_id.map(lookup[key])
    hits.to_csv(OUT / "current-cta-granule-hits.csv", index=False)
    both = int((hits.donor_42_pct_positive.ge(1) & hits.donor_52_pct_positive.ge(1)).sum())
    singleton_pct = 100 * (hits.n_positive.sum() - hits.n_ge2.sum()) / hits.n_positive.sum()
    cell_hits = all_scopes[
        all_scopes.category.eq("CTA")
        & all_scopes.scope.str.startswith("adult_type:")
        & all_scopes.detection_screen_hit.eq(True)
    ]
    extra = cell_hits[~cell_hits.entity_id.isin(hits.entity_id)]
    extra.to_csv(OUT / "additional-adult-cell-hits.csv", index=False)
    familiar = adult[adult.familiar_target].sort_values("pct_positive", ascending=False)
    (OUT / "familiar-targets.md").write_text(
        "# Familiar targets in adult cerebellar granule nuclei\n\n"
        "[Methods and interpretation](index.md). Highest-detection available member per protein; "
        "absent features are unmeasured. Current category reflects the existing catalog, not a reclassification.\n\n"
        + table(familiar)
        + "\n"
    )
    (OUT / "cta-hits.md").write_text(
        f"# {len(hits)} current CTAs at or above 1% adult granule detection\n\n"
        "[Methods and interpretation](index.md). These cross the local RNA detection criterion. "
        "This audit does not change their catalog category.\n\n" + table(hits) + "\n"
    )

    names = [
        "ENSG00000279983",
        "CBY3",
        "PLAC1",
        "DMRT1",
        "CAGE1",
        "TMEM269",
        "NUTM1",
        "SYCP1",
        "PRAME",
    ]
    chosen = adult.set_index("protein").loc[names]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1, 1.6]}
    )
    labels = [
        "No feature in matrix",
        "Zero adult detections",
        "Positive, below 1%",
        "At least 1% detection",
    ]
    values = [
        int(cta.coverage.eq("unmeasured").sum()),
        int(cta.n_positive.eq(0).sum()),
        int((cta.n_positive.gt(0) & cta.pct_positive.lt(1)).sum()),
        len(hits),
    ]
    left.barh(labels, values, color=["#9ca3af", "#568b7f", "#6288ac", "#c08341"])
    for y, value in enumerate(values):
        left.text(value + 6, y, str(value), va="center")
    left.invert_yaxis()
    left.set_xlim(0, 500)
    left.set_xlabel("Current CTA proteins")
    left.set_title("Coverage and detection: 638 CTAs", loc="left", weight="bold")
    y = np.arange(len(chosen))
    right.scatter(chosen.pct_positive, y, s=45, label="At least one count", color="#476d91")
    right.scatter(chosen.pct_ge2, y, s=45, marker="^", label="At least two counts", color="#bf8040")
    right.axvline(1, color="#777777", linestyle="--", linewidth=1)
    right.set_yticks(y, names)
    right.invert_yaxis()
    right.set_xlim(-0.1, 5.5)
    right.set_xlabel("Positive adult granule nuclei (%)")
    right.set_title("Most detections have only one count", loc="left", weight="bold")
    right.legend(loc="lower right", frameon=False)
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.02,
        0.025,
        "Same two normal-labelled adults; 16,833 granule nuclei. PLAC1 is an Unsafe comparator.\nThe two-count calculation is a sensitivity analysis, not a replacement policy. RNA detection does not establish protein production.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.1, 1, 1), w_pad=3)
    fig.savefig(OUT / "cta-comparison.png", dpi=170)
    fig.savefig(OUT / "cta-comparison.svg")
    plt.close(fig)

    examples = (
        adult.set_index("protein")
        .loc[
            [
                "PLAC1",
                "CBY3",
                "DMRT1",
                "CAGE1",
                "TMEM269",
                "NUTM1",
                "SYCP1",
                "SSX2/B",
                "MAGEB10",
                "PRAME",
                "SPAG9",
                "TSGA10",
            ]
        ]
        .reset_index()
    )
    report = f"""# Other CTAs in the PLAC1 cerebellum dataset

**{len(hits)} of the 638 current CTAs reach 1% RNA detection in the same adult granule nuclei used for the PLAC1 audit.** {both} reach 1% in both adult donors. These are newly calculated observations; the existing catalog categories have not changed.

The data are from [Sepp et al.](https://doi.org/10.1038/s41586-023-06884-x), using the [same version and normal-adult filters as the PLAC1 audit](../index.md). There are 14,102 granule nuclei from donor #105 (42 years) and 2,731 from #138 (52 years), totaling 16,833.

## Coverage comes first

| Result among current CTAs | Proteins |
|---|---:|
| No member feature in the supplied matrix | {values[0]} |
| Feature present; zero adult granule detections | {values[1]} |
| Positive; below 1% adult granule detection | {values[2]} |
| At least 1% adult granule detection | {values[3]} |
| Total | {len(cta)} |

Of the 215 proteins with at least one available gene feature, 213 have complete member coverage and two have partial coverage: SSX2/B and TP53TG3/B/C/D/E/F. Their missing members remain unmeasured. All 67 zero-detection proteins have complete member coverage in this file.

**MAGE-A4, NY-ESO-1, XAGE-1 and CT83 have no member feature in the supplied matrix.** They cannot be called brain-negative from this file. PRAME is present and has zero counts in the adult granule subset. These are different results.

The later [normal-cell atlas audit](../../normal-cell-atlas-20261110/index.md) found that the authors' workflow also filters genes within each library. Consequently, feature presence in the merged matrix does not guarantee an uncensored measurement in every donor. Zero counts here do not establish absence of adult brain RNA.

![CTA coverage and RNA-count sensitivity](cta-comparison.png)

## Named comparisons

{table(examples)}

For a protein represented by several genes, the displayed value is the member gene with the highest detection rate. Frequencies are never added across paralogs. Gene-level rows and the selected member are retained in the downloads. These are RNA observations attached to protein records, not measurements of protein abundance.

[All 46 CTA hits](cta-hits.md) · [113 familiar targets](familiar-targets.md) · [Every current CTA](adult-granule-current-ctas.csv) · [All 966 selected proteins](adult-granule-all-proteins.csv).

## The one-count result matters

Across the 46 CTA hits, **{singleton_pct:.1f}% of positive protein–nucleus observations contain exactly one count**. None reaches 1% detection if a positive nucleus must contain at least two counts. PLAC1 also drops from 3.05% to 0.048% under that calculation. A two-count threshold is a sensitivity analysis, not an adopted rule and not a validated safety boundary.

Some familiar targets already classified Unsafe have much stronger signals: SPAG9 is detected in 69.74% of granule nuclei, including 6,871 nuclei with at least two counts; TSGA10 is detected in 38.30%, including 2,193 with at least two counts. The assay therefore distinguishes widespread, multi-count expression from the sparse signals seen for PLAC1 and the current CTA hits.

Within adult oligodendrocytes (5,309 nuclei), another five current CTAs reach the same local detection criterion: ASZ1, STARD6, RIPPLY3, CPHXL and TP53TG3/B/C/D/E/F. That makes **51 current CTAs** with a local hit in either supported adult cell population. [Additional cell-type results](additional-adult-cell-hits.csv). Smaller populations are reported but lack the 3,000-cell support required for this comparison. Fetal and childhood granule results are separate in the full download.

## Why these are still CTA in the current catalog

The catalog uses the older WMG aggregate across datasets and ages. The packaged cerebellar-granule values include CAGE1 at 0.8346%, DMRT1 at 0.9851%, CBY3 at 0.8294% and NUTM1 at 0.5968%—below its 1% line. This audit uses the current downloadable matrix and adult donors, where those values are 2.002%, 2.620%, 3.333% and 1.277%. Both the source representation and the selected population differ; this comparison does not isolate which change caused the discrepancy.

The local calculation uses ≥1% detection, ≥10 positive nuclei and ≥3,000 total nuclei. The complete existing rare-cell policy additionally uses a population-share filter against all sampled brain cells. This single-study comparison does not reconstruct that global denominator or automatically rerun the classification.

**The result exposes a consistency problem for the normal-expression screen.** PLAC1 is not unique in crossing the local detection line. All such observations need the same age, source-version, mapping and background review. PLAC1 additionally has independent cerebellar PCR and adult-brain atlas support, described in the [parent audit](../index.md). The [expanded audit](../../normal-cell-atlas-20261110/index.md) now checks all 46 against the independent Siletti granule-cell matrix: 41 have available features, 40 have positive counts, and only TMEM269 reaches 1% pooled detection. Repeated sparse RNA assignments do not establish protein expression or exclude shared technical artifacts.

Sparse counts alone do not decide whether neural cells produce a recognizable antigen. The immediate follow-up is to reconcile source versions and apply a uniform adult-tissue analysis across the catalog. This audit does not silently change the threshold or relabel the proteins.

## Name of the cancer–reproductive union

The current definition calls it **Target candidates**:

**Target candidates = Cancer Antigens ∪ Reproductive Antigens.**

| Membership | Current draft proteins |
|---|---:|
| CTA: cancer and reproductive/developmental | 638 |
| Cancer Antigen only | 260 |
| Reproductive Antigen only | 22 |
| Target candidates: union | 920 |

Cancer Antigens contain 898 proteins; Reproductive Antigens contain 660. Their shared 638 CTAs are counted once in the union. Membership means the protein passed the catalog's shared checks and at least one expression route. It does not assert clinical safety. [Existing formal definition](../../target-selection-20261109/definition.json).

## Data and reproduction

- [Protein results across all analysed populations](protein-detection.csv.gz)
- [Gene-level counts across all analysed populations](gene-detection.csv.gz)
- [The 46 hits with both donor-specific results](current-cta-granule-hits.csv)
- [Category counts and coverage](category-summary.csv)
- [Source hashes, mapping rules and population sizes](manifest.json)

Run `scripts/audit_cta_cerebellum.py --h5ad tmp/plac1-brain/sepp-human.h5ad --output docs/audits/plac1-brain-20261110/cta-comparison`, then `scripts/render_cta_cerebellum.py`. The input matrix is pinned in the parent audit. The analysis requires h5py, scipy, numpy, pandas and matplotlib; no shared-environment packages were changed. Extracted adult counts for PLAC1, CAGE1, NUTM1, PRAME, CBY3 and ENSG00000279983 were independently checked with AnnData slicing.
"""
    (OUT / "index.md").write_text(report)


if __name__ == "__main__":
    main()
