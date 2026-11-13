#!/usr/bin/env python3
"""Audit abundance, source dependence and publication claims without relabeling proteins."""

import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import target_selection  # noqa: E402

OUT = ROOT / "docs/audits/target-selection-20261109/scrutiny"
REBUILD = OUT.parent / "rebuild"


def family(label):
    for name, pattern in [
        ("Olfactory receptors", r"^OR\d"),
        ("Keratin-associated proteins", r"^KRTAP"),
        ("Other sensory receptors", r"^(TAS[12]R|TAAR|VN1R)"),
        ("SCYGR", r"^SCYGR"),
        ("USP17", r"^USP17"),
        ("Unnamed Ensembl loci", r"^ENSG"),
    ]:
        if re.match(pattern, label):
            return name
    return "Other families"


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    catalog = target_selection.results()
    matrix_path = ROOT / "tmp/target-catalog-cache/tumor-comparison.parquet"
    matrix_hash = hashlib.sha256(matrix_path.read_bytes()).hexdigest()
    assert matrix_hash == target_selection.metadata()["comparison"]["tumor_matrix_cache_sha256"]
    matrix = pd.read_parquet(matrix_path)
    samples = pd.read_csv(REBUILD / "tumor-samples.csv.gz", keep_default_na=False)
    samples = samples.set_index("sample_key").loc[matrix.index]
    x = catalog[catalog.category.eq("Cancer Antigen only")].copy()
    x["family"] = x.protein_group_label.map(family)
    x["tumor_peak_over_normal"] = x.tumor_max_tpm / x.normal_reference_tpm.where(
        x.normal_reference_tpm.gt(0)
    )
    scenarios = {
        "without_hl": {"HL"},
        "without_merkel": {"NEC_MERKEL"},
        "without_both": {"HL", "NEC_MERKEL"},
    }
    sensitivity = []
    for row in catalog[catalog.is_cancer_antigen].itertuples():
        result = {
            "entity_id": row.entity_id,
            "protein": row.protein_group_label,
            "category": row.category,
        }
        for name, omitted in scenarios.items():
            keep = ~samples.cancer_code.isin(omitted)
            values = matrix.loc[keep, row.entity_id]
            above = values.ge(1) & values.gt(row.normal_reference_tpm)
            people = samples.loc[keep, "participant_id"]
            repeated = people[above & people.ne("")].nunique() >= 2
            single = values.ge(5 * max(row.normal_reference_tpm, 1)).any()
            result[name] = bool(repeated or single)
        sensitivity.append(result)
    sensitivity = pd.DataFrame(sensitivity)
    x = x.merge(
        sensitivity.drop(columns=["protein", "category"]), on="entity_id", validate="one_to_one"
    )
    x.to_csv(OUT / "cancer-antigen-only-support.csv", index=False)
    sensitivity.to_csv(OUT / "cohort-sensitivity-all-cancer-antigens.csv", index=False)
    summary = {
        "policy_id": target_selection.POLICY_ID,
        "screen_matches": len(x),
        "published_catalog_nominations": int(x.published_nomination.sum()),
        "any_paper_link": int(x.paper_source_count.gt(0).sum()),
        "families": x.family.value_counts().to_dict(),
        "peak_below_5_tpm": int(x.tumor_max_tpm.lt(5).sum()),
        "peak_below_10_tpm": int(x.tumor_max_tpm.lt(10).sum()),
        "single_sample_route_only": int(
            x.cancer_routes.eq("one_tumor_fivefold_above_normal").sum()
        ),
        "peak_cohorts": x.tumor_top_cancer_code.value_counts().to_dict(),
        "peak_sample_counts": x.tumor_top_sample.value_counts().to_dict(),
        "lose_cancer_membership_on_removal": {
            name: sensitivity[~sensitivity[name]].category.value_counts().to_dict()
            for name in scenarios
        },
        "tumor_matrix_sha256": matrix_hash,
        "catalog_sha256": target_selection.metadata()["data_sha256"],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    table = "\n".join(f"| {name} | {count} |" for name, count in summary["families"].items())
    known = ", ".join(x.loc[x.published_nomination, "protein_group_label"])
    text = f"""# Scrutiny of the 250 Cancer Antigens only

**250 proteins pass the current expression rules. This is not a count of experimentally validated cancer antigens.** The category and its comparison methods require further review before release. PR #686 is draft.

## What is in the set?

| Family | Proteins |
|---|---:|
{table}

The family counts use explicit symbol patterns, retained in the audit script. They are a descriptive breakdown, not a new annotation ontology.

## Tumor-expression concerns

- **22** have no measured tumor above 5 TPM; **69** have none at 10 TPM or higher.
- **53** qualify through the single-sample route without meeting the confirmed two-patient route. This does not mean only one measured tumor expresses them: patient identities are unresolved in several cohorts.
- **73** have their maximum in Hodgkin lymphoma (HL), including **69 in SRS3811196**. Another **64** peak in the Merkel cohort. These are concentrations of maxima, not proof of technical artifacts.
- Removing HL alone makes **21** lose qualification; removing Merkel alone makes **7** lose it. Removing both makes **31** lose it, leaving **219** under the otherwise unchanged rules. This is a sensitivity analysis, not a corrected final count.
- The same two-cohort removal affects **9 CTAs**. The concern applies to the shared cancer gate, not only the nonreproductive branch.

HL uses recount3 coverage-derived gene TPM; Merkel uses symbol-keyed gene counts divided by reference gene lengths. The source registry calls both `linear_rnaseq_tpm`, but that does not establish calibrated comparability with the normal reference. The comparison also takes tumor-sample maxima against tissue/cell summaries, including donor medians. It does not compare against the highest individual healthy donor.

The [HL study](https://pubmed.ncbi.nlm.nih.gov/30546079/) investigates LTR-driven transcription and includes laser-captured HRS cells. The [Merkel source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235092) is a bulk-RNA study. Bulk expression does not by itself identify which cell produced the RNA. None of these source descriptions proves that the suspicious values are false.

Olfactory receptors particularly need transcript-level scrutiny. A [primary RNA-seq study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3566163/) reports cases, including OR5K2, where broad read coverage could be misinterpreted as receptor transcription. OR5K2 is in this set. That is a reason to inspect its tumor alignments, not evidence that this tumor measurement is itself an artifact. Exact gene-label mapping does not establish uniquely assigned reads or an intact coding transcript.

## Which are previously nominated?

The earlier **47** flag was inflated by assigning a source-wide nomination role to every member of a mixed paper intake. The da Silva paper has a 1,103-gene testis-biased background and a distinct 745-gene tumor-positive list. These are now separated at the exact source-table and gene level.

**{summary["published_catalog_nominations"]} of the 250 have exact catalog-nomination support** in the retained sources: 22 in [da Silva et al.](https://www.oncotarget.com/article/21715/text/) and five in [Wang et al.](https://www.nature.com/articles/ncomms10499), with three overlapping. These are computational expression-based nominations, not 24 demonstrations of tumor protein production or antigen presentation.

{known}.

All **51** proteins with any retained paper link keep those links, including the 23 background-only records whose nomination flag was removed. A missing link in this curated collection does not prove that no relevant paper exists.

## Normal and developmental coverage

The 141 olfactory receptors, 31 keratin-associated proteins and other sensory receptors concentrate the question in specialized normal compartments. Ordinary tissue summaries can miss rare sensory cells. The reproductive/developmental gate records qualifying evidence in the imported sources; it does not establish absence from every fetal organ or developmental stage. Disposable allowances also deliberately permit some normal expression.

## Required follow-up

1. Review source-specific comparability, raw-count support, transcript coverage and mapping uniqueness for the shared cancer gate. Keep valid rare high expressors when their measurements are supported.
2. Audit rare normal sensory cells, skin appendages, and developmental compartments with exact-gene data. Apply the same normal-evidence rules to both cancer branches.
3. Resolve patient identities from original sample metadata, without treating aliquots as independent donors.
4. Keep source-table nominations separate from tumor protein, immune recognition and therapeutic evidence. Missing or conflicting required evidence belongs in Other:Review; a family name or a suspicious cohort alone does not establish Unsafe.

[All 250 support records](cancer-antigen-only-support.csv) · [Cohort-removal analysis for all Cancer Antigens](cohort-sensitivity-all-cancer-antigens.csv) · [Audit summary](summary.json) · [Catalog table](../rebuild/report.html).
"""
    (OUT / "index.md").write_text(text)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k not in {"peak_cohorts", "peak_sample_counts"}},
            indent=2,
        )
    )


if __name__ == "__main__":
    build()
