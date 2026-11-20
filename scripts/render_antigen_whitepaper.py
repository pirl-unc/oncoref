"""Render the antigen whitepaper and its evidence tables from a source checkout."""

import argparse
import json
import sys
from pathlib import Path
from string import Template

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_target_filter_audit import FAMILIAR, lookup  # noqa: E402

from oncoref import (  # noqa: E402
    antigen_data,
    catalog,
    cta_selection_policy,
    target_filters,
    target_selection,
)

DOCS = ROOT / "docs"
BUILD = DOCS / "audits/target-selection-20261110/rebuild"
ASSETS = "antigen-whitepaper"
NORMAL_PEAKS = DOCS / "audits/cta-selection-20261028/normal-peak-evidence.csv.gz"
TEMPLATE = ROOT / "scripts/cta_evidence_assets/antigen-whitepaper.md"
PEAK = ["tumor_max_tpm", "protein_group_label", "entity_id"]


def table(headers, rows):
    def cell(x):
        if pd.isna(x) or x == "":
            return "—"
        return str(x).replace("|", " / ").replace("\n", " ")

    return "\n".join(
        "| " + " | ".join(map(cell, line)) + " |"
        for line in [headers, ["---"] * len(headers), *rows]
    )


def number(x):
    return "—" if pd.isna(x) else f"{x:,.2f}".rstrip("0").rstrip(".")


def contexts(raw, field):
    observations = json.loads(raw) if pd.notna(raw) else []
    return " | ".join(dict.fromkeys(str(o[field]) for o in observations if o.get(field)))


def disposition(r):
    if r.category == "Unsafe":
        result = r.Unsafe_reasons
    elif r.subtype == "Unscreened tissue":
        result = f"Coverage gap in {r.home_tissue}"
    else:
        result = r.decision_reason
    if r.rare_cell_status == "established" or (
        r.subtype == "Review" and r.rare_cell_status == "unresolved"
    ):
        result += f" — {r.finding_tissue}, {r.finding_top_cell_type} ({r.rare_cell_basis})"
    return result


def study_table(proteins):
    counts = target_filters.paper_counts()
    data = proteins.join(counts, on="entity_id")
    data[counts.columns] = data[counts.columns].fillna(0).astype(int)
    familiar = {}
    for members in FAMILIAR.values():
        for name in members.split():
            row = lookup(data, name)
            if row is None:
                raise ValueError(f"Familiar protein missing from catalog: {name}")
            familiar.setdefault(row.entity_id, []).append(name)
    data["familiar_names"] = data.entity_id.map(
        {key: " | ".join(names) for key, names in familiar.items()}
    )
    trials = pd.read_csv(ROOT / "oncoref/data/cta-trial-evidence.csv")
    by_gene = data.set_index("entity_id").member_gene_ids.map(json.loads).explode()
    trials = trials.merge(
        by_gene.rename("Ensembl_Gene_ID").reset_index(),
        on="Ensembl_Gene_ID",
        how="left",
        validate="many_to_one",
    )
    if trials.entity_id.isna().any():
        raise ValueError("Clinical evidence with no catalog locus")
    for source, name in [("source_url", "clinical_sources"), ("evidence_kind", "clinical_scope")]:
        data[name] = data.entity_id.map(
            trials.groupby("entity_id")[source].agg(lambda x: " | ".join(sorted(set(x))))
        )
    data = data[data.papers_cta.gt(0) | data.familiar_names.notna() | data.clinical_sources.notna()]
    data = data.copy()
    data["disposition"] = [disposition(r) for r in data.itertuples()]
    # Quote the strongest excluding normal-expression measurement whose specific
    # gate and current injury tier still exclude the group.
    measurements = pd.read_csv(
        NORMAL_PEAKS,
        usecols=[
            "gene_id",
            "tissue",
            "cell_type",
            "source",
            "unit",
            "value",
            "candidate_ceiling",
            "selection_tier",
            "category",
            "taxonomy_id",
            "stage",
            "mapping_status",
            "unsafe_evidence",
            "assay",
        ],
        low_memory=False,
    )
    flags = data.loc[data.category.eq("Unsafe"), ["entity_id", "Unsafe_reasons"]]
    measured = (
        measurements.loc[
            measurements.unsafe_evidence
            & measurements.selection_tier.eq(cta_selection_policy.compartment_tiers(measurements))
        ]
        .merge(by_gene.rename("gene_id").reset_index(), on="gene_id", validate="many_to_one")
        .merge(flags, on="entity_id", validate="many_to_one")
    )
    matched = measured.loc[
        [r.assay in r.Unsafe_reasons.split(" | ") for r in measured.itertuples()]
    ].copy()
    matched["ceiling_ratio"] = matched.value / matched.candidate_ceiling
    strongest = matched.sort_values(
        ["ceiling_ratio", "source", "tissue", "cell_type"], ascending=[False, True, True, True]
    ).drop_duplicates("entity_id")
    locations = {
        r.entity_id: f"{r.tissue}, {r.cell_type} · {number(r.value)} {r.unit} "
        f"(ceiling {number(r.candidate_ceiling)})"
        for r in strongest.itertuples()
    }
    data["matched_normal_assay"] = data.entity_id.map(locations)
    return data.sort_values(
        ["papers_cta", "protein_group_label", "entity_id"], ascending=[False, True, True]
    )


def ctag2_evidence(data):
    row = data.loc[data.protein_group_label.eq("CTAG2")].iloc[0]
    path = DOCS / "audits/normal-cell-atlas-20261110/rare-cell"
    findings = pd.read_csv(path / "findings.csv")
    findings = findings[findings.entity_id.eq(row.entity_id)]
    risk = catalog.read_table("rare-cell-evidence").set_index("entity_id")
    return {
        "ctag2_cells": table(
            [
                "Organ",
                "Strongest population",
                "Independent studies",
                "Donor IDs with signal (linked where known)",
                "Cells with ≥2 molecules in top donor",
                "Grading",
            ],
            [
                [
                    r.tissue,
                    r.top_cell_type,
                    r.studies_multi_molecule,
                    r.donors_multi_molecule,
                    f"{r.top_n_multi:,}/{r.top_n_cells:,} ({r.top_multi_molecule_pct:.2f}%)",
                    r.status,
                ]
                for r in findings.itertuples()
            ],
        ),
        "ctag2_shared": f"{100 * risk.at[row.entity_id, 'shared_kmer_fraction']:.1f}%",
        "ctag2_classification": row.classification,
        "ctag2_site": contexts(row.normal_reference_contexts, "tissue"),
        "ctag2_reference": number(row.normal_reference_tpm),
        "ctag2_bar": number(row.single_tumor_threshold_tpm),
        "ctag2_repro": number(row.reproductive_supported_tpm),
        "ctag2_fold": number(row.reproductive_fold),
        "ctag2_tumor": number(row.tumor_max_tpm),
        "ctag2_patients": int(row.independent_patients_at_tenfold),
    }


def tissue_hierarchy():
    from oncoref import tissue_tiers

    sections = [
        "# Human tissue injury tiers",
        "The adult human rules describe loss of an organ or a named population. Expression, antigen presentation and exposure are assessed separately. Each cell rule applies to the organs listed in its assessment.",
    ]
    rules = tissue_tiers.assessments()
    references = tissue_tiers.references().set_index("reference_id")
    for tier in antigen_data.tiers().itertuples():
        sections.extend([f"## {tier.tier.title()}", tier.definition])
        records = []
        for rule in rules[rules.tier.eq(tier.tier)].sort_values(["level", "label"]).itertuples():
            links = ", ".join(f"[{key}]({references.at[key, 'url']})" for key in rule.references)
            records.append(
                [rule.label, rule.level, ", ".join(rule.tissues), rule.persistent_loss, links]
            )
        sections.append(
            table(
                [
                    "Structure or population",
                    "Resolution",
                    "Organ IDs",
                    "Persistent loss",
                    "Evidence",
                ],
                records,
            )
        )
    return "\n\n".join(sections) + "\n"


def build():
    data = target_selection.results()
    meta = target_selection.metadata()["summary"]
    data["loci"] = data.member_gene_ids.map(json.loads).map(len)
    loci = data.member_gene_ids.map(json.loads).explode()
    if not loci.is_unique or not data.entity_id.is_unique:
        raise ValueError("Catalog is not a partition of source loci")
    ca, ra = data.is_cancer_antigen, data.is_reproductive_antigen
    if not ((ca & ra).equals(data.is_cta) and (ca | ra).equals(data.is_target_candidate)):
        raise ValueError("Catalog membership flags disagree")
    for source, prefix in [
        ("normal_reference_contexts", "normal"),
        ("reproductive_supported_observations", "reproductive"),
    ]:
        for field in ["source_id", "tissue", "source_context"]:
            data[f"{prefix}_{field}"] = data[source].map(lambda x, field=field: contexts(x, field))
    funnel = pd.read_csv(BUILD / "summary/funnel-all.csv")
    if not funnel.remaining.eq(len(data) - funnel.removed.cumsum()).all() or (
        funnel.remaining.iloc[-1] != data.is_target_candidate.sum()
    ):
        raise ValueError("Funnel does not reconcile with protein membership")

    sets = {"CA": ca, "RA": ra, "CTA": ca & ra, "CA only": ca & ~ra, "RA only": ra & ~ca}
    candidates = data[ca | ra].sort_values(PEAK, ascending=[False, True, True]).copy()
    for name, mask in sets.items():
        ranked = candidates.index[
            candidates.index.isin(data.index[mask]) & candidates.tumor_max_tpm.notna()
        ]
        candidates.loc[ranked, "rank_" + name.lower().replace(" ", "_")] = pd.array(
            range(1, len(ranked) + 1), dtype="Int64"
        )
    witnesses = pd.read_csv(
        BUILD / "tumor-witnesses.csv.gz",
        usecols=["entity_id", "source_id", "sample_id", "value", "qc_pass", "identity_resolved"],
    )
    in_range = witnesses.loc[
        witnesses.entity_id.isin(candidates.entity_id)
        & witnesses.qc_pass
        & witnesses.identity_resolved
        & witnesses.value.between(136, 1300, inclusive="neither")
    ]
    peak_sections = []
    for label, mask in [
        ("CTA", ca & ra),
        ("CA only", ca & ~ra),
        ("RA only", ra & ~ca),
        ("Other", data.category.eq("Other")),
        ("Unsafe", data.category.eq("Unsafe")),
    ]:
        group = data[mask]
        measured = group.dropna(subset=["tumor_max_tpm"]).sort_values(
            PEAK, ascending=[False, True, True]
        )
        peak_sections.append(
            f"### {label}\n\nTumor comparison available for {len(measured):,} of {len(group):,} groups."
        )
        if measured.empty:
            continue
        peak_sections.append(
            table(
                [
                    "Rank",
                    "Protein",
                    "Peak cohort",
                    "Tumor TPM",
                    "N TPM",
                    "Reproductive peak TPM",
                    "Tumor fold",
                    "Known patients at 10×",
                ],
                [
                    [
                        i,
                        r.protein_group_label,
                        r.tumor_top_cancer_code,
                        number(r.tumor_max_tpm),
                        number(r.normal_reference_tpm),
                        number(r.reproductive_supported_tpm),
                        number(r.tumor_fold),
                        number(r.independent_patients_at_tenfold),
                    ]
                    for i, r in enumerate(measured.head(10).itertuples(), 1)
                ],
            )
        )
    study = study_table(data)
    inventory = antigen_data.sources().copy()
    inventory["contents"] = inventory.dataset.map(lambda name: catalog.dataset(name).description)

    # Expose the actual schema for shipped tables, including evidence fields
    # that a source title alone cannot describe.
    def columns(source):
        path = ROOT / source
        if source.startswith("oncoref/data/") and source.endswith((".csv", ".csv.gz")):
            return " | ".join(pd.read_csv(path, nrows=0).columns)
        return ""

    inventory["table_columns"] = inventory.source_file.map(columns)
    coverage = pd.read_csv(BUILD / "tumor-source-coverage.csv")
    fields = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "classification",
        "normal_reference_tpm",
        "normal_source_id",
        "normal_tissue",
        "normal_source_context",
        "disposable_max_tpm",
        "single_tumor_threshold_tpm",
        "tumor_max_tpm",
        "tumor_top_cancer_code",
        "tumor_top_sample",
        "tumor_fold",
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "anonymous_samples_at_tenfold",
        "cancer_types_at_tenfold",
        "reproductive_supported_tpm",
        "reproductive_source_id",
        "reproductive_tissue",
        "reproductive_source_context",
        "reproductive_fold",
    ]
    studied_fields = [
        "familiar_names",
        "papers_cta",
        "papers_counted",
        "nomination_source_ids",
        "clinical_sources",
        "clinical_scope",
        "disposition",
        "rare_cell_status",
        "finding_tissue",
        "finding_top_cell_type",
        "healthy_maximum_ntpm",
        "matched_normal_assay",
    ]
    outputs = {
        "tissue-hierarchy.md": tissue_hierarchy(),
        f"{ASSETS}/sources.csv": inventory.to_csv(index=False, lineterminator="\n"),
        f"{ASSETS}/tumor-ranks.csv": candidates[
            fields + [c for c in candidates if c.startswith("rank_")]
        ].to_csv(index=False, lineterminator="\n"),
        f"{ASSETS}/studied-proteins.csv": study[fields + studied_fields].to_csv(
            index=False, lineterminator="\n"
        ),
    }
    examples = study[study.familiar_names.notna() | study.clinical_scope.notna()]
    held = examples[examples.category.isin(["Unsafe", "Other"])].copy()
    other_examples = study[
        study.protein_group_label.isin(["CAGE1", "DAZL", "POTEC"]) & study.category.eq("Other")
    ].copy()
    highlighted_holds = pd.concat([held, other_examples]).drop_duplicates("entity_id").copy()

    def concise_gate(r):
        if r.category == "Unsafe":
            reasons = r.Unsafe_reasons
            if reasons == "Established rare normal-cell RNA" and pd.notna(r.rare_cell_basis):
                return f"Established rare-cell RNA ({r.rare_cell_basis})"
            return reasons
        if r.subtype == "Review" and r.rare_cell_status == "unresolved":
            return f"Unresolved rare-cell RNA ({r.rare_cell_basis})"
        return r.disposition

    def strongest_signal(r):
        if pd.notna(r.matched_normal_assay) and r.matched_normal_assay:
            return r.matched_normal_assay
        if pd.notna(r.finding_tissue) and pd.notna(r.finding_top_multi_molecule_pct):
            return (
                f"{r.finding_tissue}, {r.finding_top_cell_type} · "
                f"{r.finding_top_multi_molecule_pct:.2f}% of cells (≥2 molecules)"
            )
        return "—"

    highlighted_holds["concise_gate"] = [concise_gate(r) for r in highlighted_holds.itertuples()]
    highlighted_holds["strongest_signal"] = [
        strongest_signal(r) for r in highlighted_holds.itertuples()
    ]
    most = study.head(10)
    # Include familiar reference points even when they have fewer curated CTA papers.
    anchors = study[study.protein_group_label.isin(["NY-ESO-1", "PRAME", "CT83", "CTAG2", "CAGE1"])]
    highlighted = pd.concat([most, anchors]).drop_duplicates("entity_id")
    comparison_fields = ["normal_reference_tpm", "reproductive_supported_tpm", "tumor_max_tpm"]
    skipped = (
        highlighted.category.eq("Unsafe")
        & highlighted.cancer_evidence.eq("unassessed")
        & highlighted[comparison_fields].isna().all(axis=1)
    )
    skipped_ids = set(highlighted.loc[skipped, "entity_id"])

    numeric_cols = {"normal": "normal_reference_tpm", "reproductive": "reproductive_supported_tpm"}

    def measure(r, prefix):
        if r.entity_id in skipped_ids:
            return "NR"
        v = getattr(r, numeric_cols[prefix])
        where = getattr(r, prefix + "_source_context")
        return f"{number(v)} · {where}" if pd.notna(v) else "—"

    partition = data.groupby("classification").loci.agg(["size", "sum"])
    context_path = BUILD / "normal-rna/source-tissues.csv"
    source_coverage = pd.read_csv(context_path)
    lit_sources = pd.read_csv(BUILD / "summary/literature-sources.csv")
    published = data.published_nomination.fillna(False).astype(bool)
    other_counts = {
        key: int(data.classification.eq(f"Other:{key}").sum())
        for key in ("Review", "Sequence", "Expression", "Unscreened tissue")
    }
    profiles = target_filters.profile_summary(data).set_index("profile")
    testis_profile = profiles.loc["male-patient (protect testis)"]
    pregnant_profile = profiles.loc["female-patient (pregnant)"]
    headline = table(
        ["Stage or set", "Protein groups", "Notes"],
        [
            [
                "Human coding loci tracked",
                f"{len(loci):,}",
                f"Collapsed by identical representative amino-acid sequence into {len(data):,} protein groups",
            ],
            [
                "Excluded by normal-tissue screens (`Unsafe`)",
                f"{int(data.category.eq('Unsafe').sum()):,}",
                "Bulk RNA, single-cell clusters, immune, cornea, rare-cell RNA, normal IHC, or curated normal findings",
            ],
            [
                "Held or non-qualifying (`Other`)",
                f"{int(data.category.eq('Other').sum()):,}",
                (
                    f"Coverage/translation/rare-cell review ({other_counts['Review']:,}), "
                    f"sequence containment ({other_counts['Sequence']:,}), "
                    f"expression bounds ({other_counts['Expression']:,}), "
                    f"unscreened home tissue ({other_counts['Unscreened tissue']:,})"
                ),
            ],
            [
                "Shared-eligible before enrichment",
                f"{meta['shared_eligible']:,}",
                f"Native TPM compared for {meta['assessed_for_cancer_evidence']:,} groups (including counterfactuals for held proteins)",
            ],
            [
                "Cancer Antigens (`CA`)",
                f"{int(ca.sum()):,}",
                f"{int((ca & ra).sum()):,} CTA + {int((ca & ~ra).sum()):,} CA only (≥10× non-thymic healthy reference)",
            ],
            [
                "Reproductive Antigens (`RA`)",
                f"{int(ra.sum()):,}",
                f"{int((ca & ra).sum()):,} CTA + {int((ra & ~ca).sum()):,} RA only (≥10× on the same reference)",
            ],
            [
                "Candidate antigen union (`CA ∪ RA`)",
                f"{int((ca | ra).sum()):,}",
                f"{int(((ca | ra) & published).sum()):,} previously nominated in curated CTA papers + {int(((ca | ra) & ~published).sum()):,} newly identified",
            ],
        ],
    )
    values = {
        "policy": target_selection.POLICY_ID,
        "n": f"{len(data):,}",
        "loci": f"{len(loci):,}",
        "missing_seq": f"{data.representative_sequence_sha256s.map(json.loads).map(len).eq(0).sum():,}",
        "ca": f"{int(ca.sum()):,}",
        "ca_only": f"{int((ca & ~ra).sum()):,}",
        "ra": f"{int(ra.sum()):,}",
        "ra_only": f"{int((ra & ~ca).sum()):,}",
        "cta": f"{int((ca & ra).sum()):,}",
        "cta_published": f"{int((ca & ra & published).sum()):,}",
        "cta_new": f"{int((ca & ra & ~published).sum()):,}",
        "published_noms": f"{int(published.sum()):,}",
        "published_surviving": f"{int(((ca | ra) & published).sum()):,}",
        "published_outside": f"{int((~ca & ~ra & published).sum()):,}",
        "union": f"{int((ca | ra).sum()):,}",
        "unsafe_n": f"{int(data.category.eq('Unsafe').sum()):,}",
        "other_n": f"{int(data.category.eq('Other').sum()):,}",
        "testis_protected_cta": f"{int(testis_profile.CTA):,}",
        "pregnant_retained": f"{int(pregnant_profile.CTA):,}",
        "pregnant_unresolved": f"{int(pregnant_profile.CTA_unresolved):,}",
        "pregnant_excluded": f"{int(pregnant_profile.CTA_excluded):,}",
        "fold": f"{target_selection.ENRICHMENT_FOLD:g}",
        "floor": f"{target_selection.DETECTION_FLOOR:g}",
        "cap": f"{target_selection.DISPOSABLE_CAP:g}",
        "cells": target_selection.RARE_CELL_MIN_CELLS,
        "pct": f"{target_selection.RARE_CELL_DETECTION_PCT:g}",
        "frequent": f"{target_selection.RARE_CELL_SINGLE_STUDY_PCT:g}",
        "rna_bases": ", ".join(
            f"{v:g} for {k.title()}" for k, v in cta_selection_policy.TIER_BASES.items()
        ),
        "manageable_ceiling": f"{cta_selection_policy.MANAGEABLE_UNSAFE_FLOOR:g}",
        "shared": f"{meta['shared_eligible']:,}",
        "assessed": f"{meta['assessed_for_cancer_evidence']:,}",
        "required": len(meta["comparison_required_contexts"]),
        "files": f"{len(inventory):,}",
        "cohorts": f"{int(coverage.status.eq('assessed').sum()):,}",
        "tissue_rows": f"{len(source_coverage):,}",
        "lit_sources": f"{len(lit_sources):,}",
        "studied": f"{len(study):,}",
        "study_extra": f"{len(study) - int(published.sum()):,}",
        "held": f"{len(held):,}",
        "anchors": f"{len(examples):,}",
        "nr_proteins": ", ".join(highlighted.loc[skipped, "protein_group_label"]),
        "headline": headline,
        "sources": table(
            ["Inventory role", "Versioned files"],
            [[role, f"{count:,}"] for role, count in sorted(inventory.role.value_counts().items())],
        ),
        "partition": table(
            ["Exclusive category", "Protein records", "Member loci"],
            [[k, f"{r['size']:,}", f"{r['sum']:,}"] for k, r in partition.iterrows()],
        ),
        "funnel": table(
            ["First attributed exclusion", "Removed", "Remaining"],
            [[r.step, f"{r.removed:,}", f"{r.remaining:,}"] for r in funnel.itertuples()],
        ),
        "peaks": "\n\n".join(peak_sections),
        "union_ge_100": f"{int(candidates.tumor_max_tpm.ge(100).sum()):,}",
        "union_ge_1000": f"{int(candidates.tumor_max_tpm.ge(1000).sum()):,}",
        "union_patients_ge_10": f"{int(candidates.independent_patients_at_tenfold.ge(10).sum()):,}",
        "cta_median_tpm": number(data.loc[ca & ra, "tumor_max_tpm"].median()),
        "ca_only_median_tpm": number(data.loc[ca & ~ra, "tumor_max_tpm"].median()),
        "range_observations": f"{len(in_range):,}",
        "range_samples": f"{len(in_range[['source_id', 'sample_id']].drop_duplicates()):,}",
        "midrange": int(candidates.tumor_max_tpm.between(136, 1300, inclusive="neither").sum()),
        "publications": table(
            [
                "Protein",
                "CTA papers",
                "Catalog",
                "N TPM · site",
                "RA peak TPM · context",
                "Tumor TPM · cohort",
            ],
            [
                [
                    r.protein_group_label,
                    r.papers_cta,
                    r.classification,
                    measure(r, "normal"),
                    measure(r, "reproductive"),
                    f"{number(r.tumor_max_tpm)} · {r.tumor_top_cancer_code}"
                    if pd.notna(r.tumor_max_tpm)
                    else ("NR" if r.entity_id in skipped_ids else "—"),
                ]
                for r in highlighted.itertuples()
            ],
        ),
        "excluded": table(
            [
                "Studied protein",
                "Category",
                "Excluding or holding gate",
                "Strongest measured normal signal",
            ],
            highlighted_holds.loc[
                highlighted_holds.protein_group_label.isin(
                    [
                        "MAGEA3",
                        "MAGEA12",
                        "CTCFL",
                        "CTAG2",
                        "PLAC1",
                        "SPA17",
                        "DDX43",
                        "CLDN6",
                        "CAGE1",
                        "DAZL",
                        "POTEC",
                    ]
                ),
                ["protein_group_label", "classification", "concise_gate", "strongest_signal"],
            ].itertuples(index=False, name=None),
        ),
    }
    values.update(ctag2_evidence(data))
    outputs["antigen-whitepaper.md"] = Template(TEMPLATE.read_text()).substitute(values)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Compare to the committed artifacts")
    args = parser.parse_args()
    for name, content in build().items():
        path = DOCS / name
        if args.check:
            if not path.exists() or path.read_bytes() != content.encode("utf-8"):
                parser.exit(1, f"Stale whitepaper artifact: {name}\n")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content.encode("utf-8"))


if __name__ == "__main__":
    main()
