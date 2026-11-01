#!/usr/bin/env python3
"""Explain the simplified CTA decision and preserve native example measurements."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from oncoref import cta_selection
from oncoref.version import __version__

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-simplification-20261031"
CURRENT = ROOT / "docs/audits/cta-selection-20261028"
EXAMPLES = ["DBX1", "OR2B2", "KRTAP20-2", "OR10A4"]


def markdown(frame):
    def cell(value):
        if isinstance(value, float):
            return "Unmeasured" if np.isnan(value) else f"{value:.4g}"
        return str(value).replace("|", "/").replace("\n", " ")

    return "\n".join(
        [
            "| " + " | ".join(frame.columns) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |",
            *[
                "| " + " | ".join(map(cell, row)) + " |"
                for row in frame.itertuples(index=False, name=None)
            ],
        ]
    )


def build():
    OUT.mkdir(exist_ok=True, parents=True)
    facts = pd.read_csv(CURRENT / "protein-group-classifications.csv.gz", low_memory=False)
    previous = pd.read_csv(OUT / "baseline-1.26.0.csv.gz", low_memory=False)
    merged = previous[["protein_equivalence_group_id", "classification", "category"]].merge(
        facts[
            ["protein_equivalence_group_id", "protein_group_label", "classification", "category"]
        ],
        on="protein_equivalence_group_id",
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    assert merged.loc[merged.category_before.eq("CTA"), "category_after"].eq("CTA").all()
    assert merged.category_before.eq("Unsafe").equals(merged.category_after.eq("Unsafe"))
    merged.to_csv(OUT / "classification-transitions.csv", index=False)
    profiles = cta_selection.profiles()
    profiles.to_csv(
        OUT / "protein-evidence-profiles.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    counts = (
        facts.category.fillna("Unassessed")
        .value_counts()
        .rename_axis("Category")
        .reset_index(name="Proteins")
    )
    reasons = (
        facts[facts.category.eq("Other")]
        .subtype.value_counts()
        .rename_axis("Other reason")
        .reset_index(name="Proteins")
    )
    attributes = pd.DataFrame(
        [
            {
                "Independent attribute": "Replicated tumor RNA (two identified tumors >=10 TPM)",
                "Other proteins": int(
                    profiles.loc[profiles.category.eq("Other"), "replicated_tumor_rna"].sum()
                ),
            },
            {
                "Independent attribute": "Adult RNA below the descriptive floor",
                "Other proteins": int(
                    profiles.loc[profiles.category.eq("Other"), "adult_rna_below_floor"].sum()
                ),
            },
            {
                "Independent attribute": "Adult Disposable-only observed pattern",
                "Other proteins": int(
                    profiles.loc[profiles.category.eq("Other"), "adult_disposable_only"].sum()
                ),
            },
        ]
    )
    selected = facts.set_index("protein_group_label").loc[EXAMPLES]
    peaks = pd.read_csv(CURRENT / "normal-peak-evidence.csv.gz", low_memory=False)
    peaks = peaks[peaks.gene_id.isin(selected.representative_gene_id)].copy()
    peaks["Protein"] = peaks.gene_id.map(
        selected["representative_gene_id"]
        .reset_index()
        .set_index("representative_gene_id")
        .protein_group_label
    )
    peak_columns = [
        "Protein",
        "source",
        "version",
        "tissue",
        "cell_type",
        "value",
        "raw_maximum",
        "unit",
        "Read count",
        "candidate_ceiling",
    ]
    peaks[peak_columns].to_csv(OUT / "example-normal-evidence.csv", index=False)
    rows = []
    for name in EXAMPLES:
        row = selected.loc[name]
        ntpm = (
            peaks[peaks.Protein.eq(name) & peaks.unit.eq("nTPM")]
            .sort_values("value", ascending=False)
            .iloc[0]
        )
        tpm = (
            peaks[peaks.Protein.eq(name) & peaks.unit.eq("TPM")]
            .sort_values("value", ascending=False)
            .iloc[0]
        )
        rows.append(
            {
                "Protein": name,
                "Cancer TPM": row.oncoref_maximum_tpm,
                "Adult atlas nTPM": row.healthy_maximum_ntpm,
                "Atlas tissue / cells": ntpm.tissue + " / " + str(ntpm.cell_type),
                "Separate adult TPM": row.healthy_maximum_tpm,
                "TPM tissue / cells": tpm.tissue + " / " + str(tpm.cell_type),
                "Reproductive bulk nTPM": row.reproductive_ntpm,
                "Developmental lineage mean TPM": row.developmental_all_lineage_max_tpm,
                "Decision": row.classification,
            }
        )
    examples = pd.DataFrame(rows)
    examples.to_csv(OUT / "four-protein-review.csv", index=False)
    disposable = profiles[profiles.category.eq("Other") & profiles.adult_disposable_only]
    changed = merged[
        merged.category_before.fillna("Unassessed").ne(merged.category_after.fillna("Unassessed"))
    ]
    (
        OUT / "index.md"
    ).write_text(f"""# Simpler CTA decisions: safety, cancer and reproductive identity

**OncoRef {__version__}.** The main decision is CTA, Unsafe, Other or Unassessed.
Other reports the unmet identity requirement. Normal abundance, cancer support,
adult tissue disposition and evidence confidence are independent attributes.

[Rules and API](../../cta-selection.md) · [Full source report](../cta-source-walkthrough-20261028/index.md) · [Evidence profiles for every protein](protein-evidence-profiles.csv.gz) · [Every classification transition](classification-transitions.csv)

{markdown(counts)}

All 701 prior CTAs and all 18,062 Unsafe proteins keep their top-level categories.
The only top-level change is OR10A4: Other → Unassessed after newly reviewed adult
human nasal RNA evidence. This does not assert protein production or clinical injury.

## Five questions

1. **Normal exclusion evidence?** Assign Unsafe. Every normal layer is evaluated.
2. **Normal evidence complete and resolved?** Otherwise assign Unassessed.
3. **Qualifying reproductive presence?** Otherwise assign Other:Nonreproductive.
4. **Measured cancer expression?** Otherwise assign Other:Reproductive.
5. **Reproductive concentration passes?** Otherwise assign Other:Diffuse. All pass → CTA.

Reproductive includes adult reproductive organs, placenta, embryo and fetus.
Only normal/clinical exclusion evidence establishes Unsafe. Missing identity
evidence can establish Other after a resolved normal-screen pass.

{markdown(reasons)}

![Coarse category funnel](../cta-source-walkthrough-20261028/figures/category-funnel.png)

This view displays normal safety and coverage first. The final categories exactly
match the classifier; only the display order changes. [Stage counts](../cta-source-walkthrough-20261028/figures/category-funnel.csv).

## Why Low, Weak, Tumor-supported and Normal were confusing

The former labels answered different questions: Low described raw adult RNA;
Weak mixed subthreshold reproductive RNA with unsupported raw observations;
Normal described supported normal expression below safety ceilings;
Tumor-supported combined cancer replication with a low measured normal profile.
Their precedence made mutually compatible observations look like competing outcomes.

They are replaced by the primary Other reasons above and quantitative/boolean
columns in `cta_selection.profiles()`. A protein can simultaneously have high
tumor RNA, low recorded adult RNA and incomplete rare-cell evidence. A support
flag never erases a normal exclusion or supplies reproductive identity.

{markdown(attributes)}

Attribute counts overlap and must not be added as categories. The descriptive
2 native TPM/nTPM floor (separate 0.02 CP10K floor) describes abundance, not safety.
Normal safety still uses the existing tissue × usable-IHC ceiling and independent
protein/clinical evidence. Raw observations and source support remain inspectable.

## Disposable means an allowed adult compartment

Disposable is an adult tissue disposition: the policy permits targeting that
mapped compartment, subject to its stated limitations. It is independent of
embryonic, fetal or placental expression. The adult Disposable-only protein
attribute requires positive Disposable observations and low/no detected expression
elsewhere in the measured adult panel. It never requires developmental measurements.

Current Other proteins with this adult pattern: **{", ".join(disposable.label)}**.
LCE3B and LCE3C now receive the attribute despite missing developmental measurements.
Their developmental coverage stays explicitly unknown. Whole skin remains
Manageable; mapped follicles and cutaneous melanocytes have conditional Disposable
allowances. Bulk skin cannot be silently assigned to a particular skin compartment.

## Should the four high-tumor proteins be CTAs?

None currently meets the reproductive-presence requirement. DBX1, OR2B2 and
KRTAP20-2 remain Other:Nonreproductive, with strong cancer evidence exposed as an
attribute. They are potential cancer-target leads; low bulk atlas values alone
do not establish restriction across relevant rare normal cells. OR10A4 requires
normal nasal review before direct targeting or shared-sequence allowance.

{markdown(examples)}

These adult values are **source-supported atlas summaries**, not the largest
possible donor observation. DICE uses a supported third-highest value; sparse
HPA clusters require at least 20 reads. TPM, nTPM and qPCR copy estimates remain
different measurement scales. The developmental column is the highest recorded
lineage mean, including lineages that lack admission support; none passes the
current ≥5 TPM presence floor. No convincing reproductive IHC route is established
for these four. [Per-source adult evidence](example-normal-evidence.csv).

### OR10A4 has additional normal human nasal evidence

[Verbeurgt et al., *Profiling of Olfactory Receptor Gene Expression in Whole Human Olfactory Mucosa* (2014)](https://doi.org/10.1371/journal.pone.0096333)
reports OR10A4 in all 26 adult donors: mean estimated **167.6 RNA copies per 20 ng
total RNA**, range **33.1–548.9**, with Cq 27.21–31.61. The original Tables S1/S3 are
pinned as `olfactory_verbeurgt_cq` and `olfactory_verbeurgt_copies`, version
`2014-05-06`, through the shared `reference_data` / `catalog` cache API.
OR2B2 was not found under the retrieved exact/legacy labels; that is unassessed,
not zero expression or permission to substitute a paralog.

[*Novel Nasal Epithelial Cell Markers of Parkinson’s Disease…* (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7408990/)
independently reports OR10A4 RNA in nasal-fluid samples, including 12 healthy
controls. Both are RNA evidence in a mixed sampling compartment. Neither supplies
normal protein abundance, a TPM conversion or demonstrated T-cell injury.

## API migration

```python
from oncoref import cta_selection

decisions = cta_selection.results()
profiles = cta_selection.profiles()
targets = decisions[decisions.category.eq("CTA")]
tumor_leads = profiles[
    profiles.category.eq("Other") & profiles.replicated_tumor_rna
]
adult_disposable = profiles[profiles.adult_disposable_only.fillna(False)]
```

The top-level API and identity fields remain stable. Other subtypes are now
Nonreproductive, Reproductive and Diffuse. Replace filters for Low, Weak, Normal,
Disposable or Tumor-supported with the relevant independent profile fields.
Human evidence never transfers to dog or mouse; their profile measurements remain
missing. Existing historical sequence bundles retain their own policy/reference
identity and require deliberate consumer migration.
""")
    source_files = [
        CURRENT / "protein-group-classifications.csv.gz",
        ROOT / "oncoref/data/gene-knowledge.json",
        Path(__file__),
    ]
    manifest = {
        "policy": facts.policy.unique().tolist(),
        "category_counts": counts.set_index("Category").Proteins.to_dict(),
        "other_reason_counts": reasons.set_index("Other reason").Proteins.to_dict(),
        "top_level_changes": changed.protein_group_label.tolist(),
        "inputs": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source_files
        },
        "outputs": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(OUT.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(counts.to_string(index=False))
    print(reasons.to_string(index=False))


if __name__ == "__main__":
    build()
