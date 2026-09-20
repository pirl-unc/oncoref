#!/usr/bin/env python3
"""Build a guide from verified current selections, with a semantic evidence index."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd
from cta_report_common import deterministic_zip, seal_stage, verify_analysis, verify_stage
from cta_report_render import Report
from pypdf import PdfReader, PdfWriter


def replace_start_section(content, replacement):
    pattern = r"(?ms)^## Start here\n.*?(?=^## |\Z)"
    section = "## Start here\n\n" + replacement.strip() + "\n\n"
    return (
        re.sub(pattern, lambda match: section, content, count=1)
        if re.search(pattern, content)
        else content.rstrip() + "\n\n" + section
    )


def supporting_sections(index, page_count):
    if len(index) != page_count or sorted(index.pdf_page) != list(range(1, page_count + 1)):
        raise ValueError("PDF page index does not match the evidence PDF")
    if index.category.isna().any():
        raise ValueError("Every evidence page requires a category")
    return [
        (category, rows.sort_values("pdf_page"))
        for category, rows in index.groupby("category", sort=True)
    ]


def build_guide(out):
    out = Path(out)
    verify_analysis(out)
    for stage in ("primary", "burden", "render", "hpa"):
        verify_stage(out, stage)
    dest = out / "primary_panel"
    primary = pd.read_csv(dest / "primary_proteoforms.csv")
    burden = pd.read_csv(dest / "primary_burden_coverage.csv")
    if set(primary.proteoform_key) != set(burden.proteoform_key):
        raise ValueError("Primary selection and burden rows disagree")
    primary = primary.merge(
        burden.drop(columns="Symbol"), on="proteoform_key", validate="one_to_one"
    )
    primary = primary.sort_values(
        ["p90_n_cancer_type_groups_gt10", "Symbol"], ascending=[False, True]
    )
    annotations = pd.read_csv(out / "proteoform_member_annotations.csv").set_index("Symbol")
    handoff = []
    for row in primary.itertuples():
        for symbol in row.member_symbols.split(";"):
            gene = annotations.loc[symbol]
            if gene.Ensembl_Gene_ID not in row.member_gene_ids.split(";"):
                raise ValueError(f"Mismatched member identity: {symbol}")
            handoff.append(
                {
                    "protein_group": row.Symbol,
                    "gene_symbol": symbol,
                    "ensembl_gene_id": gene.Ensembl_Gene_ID,
                    "proteoform_key": row.proteoform_key,
                    "protein_length_aa": row.protein_length_aa,
                    "cancer_type_groups_gt10_p90": row.p90_n_cancer_type_groups_gt10,
                    "cancer_type_groups_denominator": row.p90_n_eligible_cancer_type_groups,
                    "world_incidence_pct_represented": row.world_incidence_pct_represented,
                    "world_mortality_pct_represented": row.world_mortality_pct_represented,
                }
            )
    frame = pd.DataFrame(handoff)
    frame.to_csv(dest / "primary_gene_handoff.csv", index=False)
    for column, filename in [
        ("gene_symbol", "primary_gene_symbols.txt"),
        ("ensembl_gene_id", "primary_ensembl_gene_ids.txt"),
    ]:
        (dest / filename).write_text("\n".join(sorted(set(frame[column]))) + "\n")
    guide_path = out / "cta-working-list-guide.pdf"
    guide = Report(guide_path)
    guide.start(f"Current working list: {len(primary)} protein groups")
    guide.paragraph(
        f"The >50% prevalence / p90 selection contains {len(primary)} protein identities and {len(frame)} member genes. Each identity qualifies in at least one cohort with at least 20 patient groups. This is a threshold-selected candidate pool; unique patient coverage has not been optimized."
    )
    guide.paragraph(
        "Broadest observed p90 breadth: "
        + "; ".join(
            f"{r.Symbol}: {int(r.p90_n_cancer_type_groups_gt10)}/{int(r.p90_n_eligible_cancer_type_groups)} cancer-type groups"
            for r in primary.head(4).itertuples()
        )
        + "."
    )
    guide.paragraph(
        "Genes with identical longest Ensembl 112 protein sequences share one entry. Import member-gene symbols from the handoff CSV. Transcriptome cutoffs use a fixed common background, collapsed consistently; missing members cannot qualify. Normal RNA uses within-tissue sums of the same loci. Rounded zero and unavailable measurements remain distinct."
    )
    method = json.loads((dest / "burden_coverage_method.json").read_text())
    guide.paragraph(
        f"Breadth audit: {method['n_cohort_views']} cohort views; {method['n_eligible_cohort_views']} with at least 20 patients; {method['n_eligible_cancer_type_groups']} cancer-type groups. Burden shares sum each represented category once. A subtype hit represents its parent category, so these are not percentages of worldwide patients expressing a protein."
    )
    guide.finish("Current selection and interpretation", "guide")
    display = primary[
        [
            "Symbol",
            "member_symbols",
            "protein_length_aa",
            "p90_n_cancer_type_groups_gt10",
            "world_incidence_pct_represented",
            "world_mortality_pct_represented",
        ]
    ].copy()
    display.columns = [
        "Protein",
        "Member genes",
        "Length (aa)",
        "Cancer-type groups >10%",
        "World incidence represented (%)",
        "World mortality represented (%)",
    ]
    for col in display.columns[-2:]:
        display[col] = display[col].map(lambda value: f"{value:.1f}")
    guide.table(display, "Importable identifiers and evidence", "guide")
    for name, title in [
        ("primary_common_cancer_heatmap", "Observed coverage in broad cancer cohorts"),
        ("list_size_coverage_comparison", "Compare selections at matched expression thresholds"),
        ("primary_normal_tissue_scatter", "Normal-tissue RNA and cancer breadth"),
    ]:
        guide.figure(dest / f"{name}.png", title, "guide", name)
    guide.save()
    source = PdfReader(out / "cta-analysis-all-figures.pdf")
    index = pd.read_csv(out / "pdf_page_index.csv")
    writer = PdfWriter()
    for page in PdfReader(guide_path).pages:
        writer.add_page(page)
    writer.add_outline_item("Current selection guide", 0)
    mapping = []
    for category, rows in supporting_sections(index, len(source.pages)):
        parent = writer.add_outline_item(str(category), len(writer.pages))
        for row in rows.itertuples():
            writer.add_page(source.pages[int(row.pdf_page) - 1])
            writer.add_outline_item(row.title, len(writer.pages) - 1, parent=parent)
            mapping.append(
                {
                    "original_page": row.pdf_page,
                    "new_page": len(writer.pages),
                    "section": category,
                    "title": row.title,
                }
            )
    writer.add_metadata({"/Title": "CTA working list and complete evidence library"})
    writer.write(out / "cta-analysis-primary-panel.pdf")
    pd.DataFrame(mapping).to_csv(dest / "supporting_page_map.csv", index=False)
    readme = dest / "README.md"
    content = readme.read_text() if readme.exists() else "# Primary CTA panel\n"
    readme.write_text(
        replace_start_section(
            content,
            f"Use ../cta-working-list-guide.pdf for the current guide and ../cta-analysis-primary-panel.pdf for all evidence. primary_gene_handoff.csv maps {len(frame)} genes to {len(primary)} protein groups. See burden_coverage_method.json for category-scope caveats and coverage_denominator_audit.csv for the full cohort audit.",
        )
    )
    seal_stage(
        out,
        "guide",
        [
            Path(__file__),
            *[out / f"{stage}_receipt.json" for stage in ("primary", "burden", "render", "hpa")],
        ],
        [
            guide_path,
            out / "cta-analysis-primary-panel.pdf",
            dest / "primary_gene_handoff.csv",
            dest / "primary_gene_symbols.txt",
            dest / "primary_ensembl_gene_ids.txt",
            dest / "supporting_page_map.csv",
            readme,
        ],
    )
    deterministic_zip(out, out / "cta-analysis-bundle.zip")
    print(
        f"Built {len(PdfReader(guide_path).pages)} guide pages and {len(writer.pages)} full-report pages"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs/cta_proteoform_report_20260917",
    )
    build_guide(parser.parse_args().out.resolve())


if __name__ == "__main__":
    main()
