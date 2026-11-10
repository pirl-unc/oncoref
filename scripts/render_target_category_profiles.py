#!/usr/bin/env python3
"""List all published nominations, exclusive antigen sets and normal RNA peaks."""

import gzip
import html
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from audit_target_category_profiles import BASE, CACHE, OUT, sha, write  # noqa: E402

from oncoref import reference_data, target_selection, tissue_hazards  # noqa: E402
from oncoref.gene_matching import gene_identity_catalog  # noqa: E402

SOURCES = {
    "rna_tissue_hpa": "HPA tissue mean",
    "rna_tissue_gtex": "HPA/GTEx tissue mean",
    "rna_brain_gtex": "HPA/GTEx brain mean",
    "rna_brain_hpa": "HPA brain mean",
    "gtex_tissue_rna": "GTEx v11 median",
    "hpa_brain": "HPA brain",
    "hpa_rna_consensus": "HPA consensus",
    "hpa_single_cell_tissue": "HPA tissue cells",
    "hpa_immune_cell": "HPA immune",
    "hpa_immune_monaco": "Monaco immune",
    "GSE155683_adult": "Adult cornea",
    "cornea_collin_adult_counts": "Adult cornea",
}

CELL_LABELS = {
    "B_CELL_NAIVE": "naive B cells",
    "CD4_NAIVE": "naive CD4 T cells",
    "CD4_STIM": "stimulated CD4 T cells",
    "CD8_NAIVE": "naive CD8 T cells",
    "CD8_STIM": "stimulated CD8 T cells",
    "TREG_MEM": "memory regulatory T cells",
    "TREG_NAIVE": "naive regulatory T cells",
    "TFH": "follicular helper T cells",
    "NK": "natural killer cells",
    "MONOCYTES": "monocytes",
}


def source_label(source):
    source = source.split(":")[0]
    return "DICE third-highest" if source.startswith("dice_") else SOURCES.get(source, source)


def location(row):
    if pd.notna(row.value) and row.value == 0:
        return "No supported positive measurement"
    context = (
        row.source_context if pd.notna(row.source_context) and row.source_context else row.tissue
    )
    cell = row.cell_type if pd.notna(row.cell_type) else ""
    label = CELL_LABELS.get(cell, cell)
    context = str(context).replace("_", " ")
    return context if not cell or cell == "bulk" or cell == context else f"{context} / {label}"


def url(source, gene):
    source = source.split(":")[0]
    if source.startswith("dice_"):
        return "https://dice-database.org/"
    if source == "gtex_tissue_rna":
        return f"https://gtexportal.org/home/gene/{gene}"
    if source in {"GSE155683_adult", "cornea_collin_adult_counts"}:
        return "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155683"
    return f"https://v23.proteinatlas.org/{gene}"


def value_text(value):
    if pd.isna(value):
        return "—"
    if value == 0:
        return "0"
    return f"{value:.3g}" if value < 0.01 else f"{value:,.3f}".rstrip("0").rstrip(".")


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    receipt = json.loads((CACHE / "receipt.json").read_text())
    for path, expected in receipt["inputs"].items():
        assert sha(ROOT / path) == expected, path
    for path, expected in receipt["outputs"].items():
        assert sha(CACHE / path) == expected, path
    catalog = target_selection.results()
    screened = pd.read_csv(CACHE / "screen-peaks.csv.gz", low_memory=False)
    native = pd.read_csv(CACHE / "comparison-peaks.csv.gz", low_memory=False)
    native["source"] = native.source_id.str.split(":").str[0]
    native["version"] = native.source_id.str.split(":").str[-1]
    context_keys = ["tissue", "cell_type", "source"]
    tiers = {}
    for t, c, source in native[context_keys].drop_duplicates().itertuples(index=False, name=None):
        finding = tissue_hazards.assess(
            t, c, resolution="cell" if source in tissue_hazards.CELL_SOURCES else "bulk"
        )
        tiers[(t, c, source)] = finding.get("fallback_tier", finding["tier"])
    native["selection_tier"] = [
        tiers[key] for key in native[context_keys].itertuples(index=False, name=None)
    ]
    native["assay"] = "Additional native-TPM comparison"
    native["gene_id"] = native.entity_id.map(catalog.set_index("entity_id").representative_gene_id)
    all_peaks = pd.concat([screened, native], ignore_index=True)
    assert set(all_peaks.unit) <= {"TPM", "nTPM", "CP10K"}
    all_peaks = all_peaks.sort_values(
        ["value", "source", "tissue", "cell_type"], ascending=[False, True, True, True]
    )
    peaks = all_peaks.drop_duplicates(["entity_id", "unit"]).copy()
    peaks["location"] = peaks.apply(location, axis=1)
    peaks["source_label"] = peaks.source.map(source_label)
    peaks["source_url"] = [url(r.source, r.gene_id) for r in peaks.itertuples()]
    peaks = peaks.merge(
        catalog[["entity_id", "protein_group_label", "category"]],
        on="entity_id",
        validate="many_to_one",
    )
    write(peaks, OUT / "normal-peaks.csv.gz")
    write(all_peaks, OUT / "normal-and-disposable-source-peaks.csv.gz")
    shutil_inputs = [
        "duplicate-normal-observations.csv.gz",
        "incomplete-normal-contexts.csv.gz",
        "receipt.json",
    ]
    for name in shutil_inputs:
        (OUT / ("extraction-receipt.json" if name == "receipt.json" else name)).write_bytes(
            (CACHE / name).read_bytes()
        )

    # The complete audit must reproduce the native-TPM comparator for every
    # shared-eligible record. Its new Disposable-inclusive maximum can be higher.
    non_disposable = native[native.scope.eq("normal")].set_index("entity_id").value
    expected = catalog.set_index("entity_id").normal_reference_tpm.dropna()
    assert np.allclose(non_disposable.reindex(expected.index), expected)

    result = catalog.copy()
    result = result.rename(columns={"protein_group_label": "protein"})
    identity = gene_identity_catalog().set_index("canonical_gene_id").canonical_symbol
    result["member_symbols"] = result.member_gene_ids.map(
        lambda ids: json.dumps([identity.get(g, g) for g in json.loads(ids)])
    )
    for unit, key in [("TPM", "tpm"), ("nTPM", "ntpm"), ("CP10K", "cp10k")]:
        subset = peaks[peaks.unit.eq(unit)].set_index("entity_id")
        for source_col, target in [
            ("value", "normal_peak"),
            ("location", "normal_location"),
            ("source_label", "normal_source"),
            ("source_url", "normal_source_url"),
            ("selection_tier", "normal_tier"),
            ("scope", "normal_scope"),
            ("candidate_ceiling", "normal_ceiling"),
            ("negative_ihc_grade", "normal_ihc"),
        ]:
            result[f"{target}_{key}"] = result.entity_id.map(subset[source_col])
    with gzip.open(BASE / "rebuild/paper-attestations.json.gz", "rt") as f:
        attestations = json.load(f)
    claims = []
    paper_links = {}
    for row in result.itertuples():
        members = set(json.loads(row.member_gene_ids))
        links = set()
        for i in json.loads(row.paper_attestation_indices):
            paper = attestations[i]
            covered = members & set(paper["nomination_gene_ids"])
            if not covered:
                continue
            links.add((paper["citation"], paper["source_url"]))
            claims.append(
                {
                    "entity_id": row.entity_id,
                    "protein": row.protein,
                    "category": row.category,
                    "source_id": paper["source_id"],
                    "citation": paper["citation"],
                    "source_url": paper["source_url"],
                    "attested_members": json.dumps(sorted(covered)),
                    "claim": paper["claim_scope"],
                    "table_receipts": json.dumps(paper["publication_claims"]),
                }
            )
        paper_links[row.entity_id] = sorted(links)
        assert bool(links) == row.published_nomination
    result["published_target_links"] = result.entity_id.map(lambda e: json.dumps(paper_links[e]))
    fields = [
        "entity_id",
        "protein",
        "representative_gene_id",
        "member_gene_ids",
        "member_symbols",
        "category",
        "classification",
        "decision_reason",
        "Unsafe_reasons",
        "published_nomination",
        "nomination_scope",
        "published_target_links",
        "cancer_evidence",
        "reproductive_evidence",
        "normal_reference_tpm",
        "tumor_max_tpm",
        "tumor_top_cancer_code",
        "tumor_top_sample",
        "independent_patients_above_normal",
        "unresolved_patient_samples_above_normal",
        "reproductive_ntpm",
        "reproductive_cell_ntpm",
        "developmental_mean_tpm",
        *[
            c
            for c in result
            if c.startswith(
                (
                    "normal_peak_",
                    "normal_location_",
                    "normal_source_",
                    "normal_tier_",
                    "normal_scope_",
                    "normal_ceiling_",
                    "normal_ihc_",
                )
            )
        ],
    ]
    result = result[fields].copy()
    write(pd.DataFrame(claims).drop_duplicates(), OUT / "publication-claims.csv.gz")
    published = result[result.published_nomination]
    write(published, OUT / "published-targets.csv.gz")
    write(result, OUT / "all-protein-profiles.csv.gz")
    tops = []
    for category in target_selection.CATEGORIES:
        subset = result[result.category.eq(category)].sort_values("protein")
        slug = category.lower().replace(" ", "-")
        suffix = ".csv.gz" if category in {"Other", "Unsafe"} else ".csv"
        write(subset, OUT / f"{slug}{suffix}")
        if suffix == ".csv.gz":
            (OUT / f"{slug}.csv").unlink(missing_ok=True)
        genes = subset.assign(gene_id=subset.member_gene_ids.map(json.loads)).explode("gene_id")
        genes["gene_symbol"] = genes.gene_id.map(identity).fillna(genes.gene_id)
        write(genes, OUT / f"{slug}-genes.csv.gz")
        for unit, key in [("TPM", "tpm"), ("nTPM", "ntpm"), ("CP10K", "cp10k")]:
            top = (
                subset[subset[f"normal_peak_{key}"].notna()]
                .sort_values([f"normal_peak_{key}", "protein"], ascending=[False, True])
                .head(20)
                .copy()
            )
            top["ranking_unit"] = unit
            top["rank"] = range(1, len(top) + 1)
            write(top, OUT / f"top20-{slug}-{key}.csv")
            tops.append(top)
    write(pd.concat(tops, ignore_index=True), OUT / "top20-all-categories.csv.gz")
    counts = published.category.value_counts()
    counts_md = "\n".join(f"| {c} | {counts.get(c, 0)} |" for c in target_selection.CATEGORIES)
    source_rows = (
        pd.DataFrame(claims).groupby(["citation", "source_url"]).entity_id.nunique().sort_index()
    )
    sources_md = "\n".join(
        f"| [{citation}]({source_url}) | {count:,} |"
        for (citation, source_url), count in source_rows.items()
    )
    familiar = pd.read_csv(BASE / "policy-reference/named-cta-proteins.csv")
    familiar = result[result.entity_id.isin(familiar.entity_id)].sort_values("protein")
    write(familiar, OUT / "familiar-cta-profiles.csv")
    familiar_rows = []
    for r in familiar.itertuples():
        reason = r.Unsafe_reasons if pd.notna(r.Unsafe_reasons) else r.decision_reason
        papers = "; ".join(f"[{n}]({u})" for n, u in paper_links[r.entity_id]) or "—"
        familiar_rows.append(
            f"| [{r.protein}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r.representative_gene_id}) | {r.classification} | {value_text(r.normal_peak_tpm)} | {str(r.normal_location_tpm).replace('|', ';')} | {str(reason).replace('|', ';')} | {papers} |"
        )
    (OUT / "familiar-ctas.md").write_text(
        "# Familiar published CTA families and studied targets\n\nThe named panel has 113 protein records. The complete curated published-target list has "
        + str(len(published))
        + " records and is available in the [searchable report](report.html). Current categories use the fivefold rule.\n\n| Protein | Current classification | Highest normal TPM | Context | Recorded reason | Published nominations |\n|---|---|---:|---|---|---|\n"
        + "\n".join(familiar_rows)
        + "\n"
    )
    rank_sections = []
    for category in target_selection.CATEGORIES:
        slug = category.lower().replace(" ", "-")
        rank_sections.append(f"## {category}\n\n")
        for unit, key in [("TPM", "tpm"), ("nTPM", "ntpm"), ("CP10K", "cp10k")]:
            rows = (
                result[result.category.eq(category) & result[f"normal_peak_{key}"].notna()]
                .sort_values([f"normal_peak_{key}", "protein"], ascending=[False, True])
                .head(20)
            )
            lines = []
            for r in rows.itertuples():
                loc = getattr(r, f"normal_location_{key}")
                src = getattr(r, f"normal_source_{key}")
                link = getattr(r, f"normal_source_url_{key}")
                tier = getattr(r, f"normal_tier_{key}")
                ceiling = value_text(getattr(r, f"normal_ceiling_{key}"))
                ihc = getattr(r, f"normal_ihc_{key}")
                ihc = "—" if pd.isna(ihc) else ihc
                lines.append(
                    f"| [{r.protein}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r.representative_gene_id}) | {value_text(getattr(r, f'normal_peak_{key}'))} | {loc} | {tier} | {ceiling} | {ihc} | [{src}]({link}) |"
                )
            rank_sections.append(
                f"### {unit}\n\n[Download ranking](top20-{slug}-{key}.csv). Equal values are ordered by protein label. A zero is not evidence of localization to the arbitrarily selected zero-valued context. A missing ceiling means no ceiling was assigned on that source scale (including additional comparison-only sources and Disposable contexts).\n\n| Protein | Highest supported {unit} | Normal tissue/cell | Tier | Screen ceiling | Matched negative IHC | Source |\n|---|---:|---|---|---:|---|---|\n"
                + "\n".join(lines)
                + "\n\n"
            )
    (OUT / "normal-top20.md").write_text(
        "# Highest noncancer, nonreproductive RNA expression\n\nDisposable tissues are included. Rankings are by supported RNA abundance, not predicted injury. TPM, nTPM and CP10K are separate rankings. Values are source-reported means, donor medians, supported cell summaries or DICE third-highest observations; they are not maxima over every healthy individual. Native-TPM comparison sources supplement the normal-screen sources. The table shows one context when maxima tie.\n\n"
        + "".join(rank_sections)
    )
    lists = []
    for category in ["Cancer Antigen only", "Reproductive Antigen only"]:
        rows = result[result.category.eq(category)].sort_values("protein")
        slug = category.lower().replace(" ", "-")
        loci = sum(len(json.loads(v)) for v in rows.member_gene_ids)
        lines = [
            f"| [{r.protein}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r.representative_gene_id}) | {', '.join(json.loads(r.member_symbols))} | {value_text(r.normal_peak_tpm)} | {r.normal_location_tpm} | {value_text(r.tumor_max_tpm)} | {r.cancer_evidence} | {r.reproductive_evidence} |"
            for r in rows.itertuples()
        ]
        lists.append(
            f"## {category}: {len(rows)} proteins / {loci} genes\n\n[Protein CSV]({slug}.csv) · [One row per gene]({slug}-genes.csv.gz). “Only” describes membership in one qualifying set; it does not mean zero expression in the other context.\n\n| Protein record | Member genes | Highest normal TPM | Normal context | Highest tumor TPM | Cancer evidence | Reproductive evidence |\n|---|---|---:|---|---:|---|---|\n"
            + "\n".join(lines)
            + "\n"
        )
    (OUT / "exclusive-sets.md").write_text(
        "# Every exclusive antigen-set member\n\nCurrent draft, with corrected patient identities and the fivefold single-tumor rule. Normal peaks include Disposable tissues.\n\n"
        + "\n".join(lists)
    )
    main_peaks = []
    for c in ["CTA", "Cancer Antigen only", "Reproductive Antigen only"]:
        r = result[result.category.eq(c)].sort_values("normal_peak_tpm", ascending=False).iloc[0]
        main_peaks.append(
            f"| {c} | {r.protein} | {value_text(r.normal_peak_tpm)} | {r.normal_location_tpm} | {r.normal_source_tpm} |"
        )
    summary = {
        "categories": catalog.category.value_counts().to_dict(),
        "curated_published_nomination_proteins": len(published),
        "published_nomination_categories": counts.to_dict(),
        "published_sources": len({c["source_id"] for c in claims}),
        "cancer_only_proteins": 260,
        "cancer_only_loci": 261,
        "reproductive_only_proteins": 22,
        "reproductive_only_cancer_status": result[result.category.eq("Reproductive Antigen only")]
        .cancer_evidence.value_counts()
        .to_dict(),
        "ranking_scope": "All supported nonreproductive normal contexts, including Disposable; each unit separate; classifications unchanged",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (
        OUT / "index.md"
    ).write_text(f"""# Published CTAs, exclusive antigen sets and normal-expression peaks

**MAGE-A4, PRAME, NY-ESO-1, XAGE-1 and CT83 remain CTA.** These are the current draft's classifications, after the patient-identity correction, using the fivefold single-tumor rule.

[Search and filter every protein](report.html) · [113 familiar CTA-family proteins and studied targets](familiar-ctas.md) · [All exclusive-set members](exclusive-sets.md) · [The highest 20 in every category and unit](normal-top20.md).

## All curated published target nominations

The OncoRef ledger contains **{len(published):,} protein records with a curated published CA/CTA/oncofetal target nomination from {summary["published_sources"]} sources**. This includes broad catalog nominations and experimental target papers. It does not label every entry as an experimentally validated CTA. Reproductive-only background lists are excluded from this nomination flag, while their citations remain in the main catalog.

| Current category | Published nominations |
|---|---:|
{counts_md}

[Every published protein](published-targets.csv.gz) · [Every exact-member paper claim and table receipt](publication-claims.csv.gz). The [familiar-name table](familiar-ctas.md) shows why MAGEC2 is Other:Review, CT45A3 is Other:Sequence, and MAGEA3/MAGEA12, CTCFL/BORIS and HORMAD1 are Unsafe. Published nomination does not override normal-expression evidence or an identity/translation hold.

The 113 familiar CTA-family proteins and studied targets divide into **65 CTA, 43 Unsafe, three Other, one Cancer Antigen only and one Reproductive Antigen only**.

These sources supply the nomination flag. Counts are current protein records with an exact-member claim, not the original papers' gene totals; sources overlap. The claim download retains each source's role and table locations.

| Paper | Matching protein records |
|---|---:|
{sources_md}

## Cancer-only and reproductive-only

**260 Cancer Antigen only protein records cover 261 genes.** They include 147 olfactory receptors, 32 keratin-associated proteins and 81 other records. PPIAL4E/F is one protein record with two member genes. [Complete list and values](exclusive-sets.md#cancer-antigen-only-260-proteins-261-genes).

**22 proteins are Reproductive Antigen only.** Fifteen have unresolved patient independence for potentially qualifying repeated tumor measurements; seven do not meet either cancer RNA rule. Their category therefore does not mean “never expressed in cancer.” [Complete list](exclusive-sets.md#reproductive-antigen-only-22-proteins-22-genes).

## Largest normal values in the three antigen sets

These peaks **include Disposable tissues** and exclude reproductive/developmental contexts. The primary table uses native TPM. The complete ranking has separate nTPM and CP10K tables, plus Other and Unsafe. RNA values are attributed to protein records through their member genes; they are not protein-abundance measurements.

| Category | Protein | Highest normal TPM | Context | Source |
|---|---|---:|---|---|
{chr(10).join(main_peaks)}

The normal screen permits different RNA ceilings by cell/tissue tier and matched IHC evidence, so a high normal value need not itself create Unsafe. “Highest” here means the highest supported source statistic: tissue means, GTEx donor medians, cell summaries or DICE third-highest TPM. It is not the largest value in any individual healthy donor. DICE unpaired-member upper bounds retain that label in the download. HPA cell RNA requires ≥20 reads; low-read raw signals do not enter this supported ranking. No conversion or pooled ranking is made between TPM, nTPM and CP10K.

PRAME has **20.3 nTPM in undifferentiated colon cells**, and CTAG2/LAGE-1 has **24.9 nTPM in bronchial basal cells**. Both retain CTA status under a **32-nTPM ceiling with Enhanced matched negative IHC**. These are concrete consequences of the current IHC-adjusted rule. The rankings display the source-scale ceiling and IHC grade where assigned.

The full-inventory expansion quarantines conflicting duplicate normal locus/context vectors and collapses identical duplicates; [every affected source row](duplicate-normal-observations.csv.gz) is retained. [Incomplete protein-member observations](incomplete-normal-contexts.csv.gz), including unbounded DICE group estimates, do not enter a quantitative maximum. Infinite upper bounds represent missing coverage, not measured expression. The resulting non-Disposable native-TPM maxima reproduce the existing comparator for all shared-eligible proteins. Classification outputs are unchanged. The hair-follicle mapping and coverage work remains tracked in [#687](https://github.com/pirl-unc/oncoref/issues/687).

Reproduce with `scripts/audit_target_category_profiles.py`, then `scripts/render_target_category_profiles.py`. [Provenance](manifest.json) · [All normal peaks](normal-peaks.csv.gz) · [Separate normal/Disposable source peaks](normal-and-disposable-source-peaks.csv.gz).
""")
    render_html(result)
    inputs = [
        ROOT / "oncoref/data/target-selection.csv.gz",
        BASE / "rebuild/paper-attestations.json.gz",
        BASE / "policy-reference/named-cta-proteins.csv",
        CACHE / "receipt.json",
        CACHE / "screen-peaks.csv.gz",
        CACHE / "comparison-peaks.csv.gz",
    ]
    inputs += [
        ROOT / "scripts/replay_cta_selection.py",
        ROOT / "scripts/rebuild_target_catalog.py",
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "oncoref/tissue_hazards.py",
        ROOT / "oncoref/data/normal-tissue-hazards.yaml",
        ROOT / "docs/audits/cta-selection-20261028/alternative-evidence/normal-bulk.csv.gz",
        ROOT
        / "docs/audits/cta-transcript-selection-20261108/inputs/regrouped/group-normal/cluster-contexts.csv.gz",
    ]
    frozen = ROOT / "docs/audits/cta-transcript-selection-20261108/inputs/frozen"
    inputs += [
        frozen / n
        for n in [
            "bulk-contexts.csv.gz",
            "immune-cornea-contexts.csv.gz",
            "dice-group-contexts.csv.gz",
            "cornea-group-contexts.csv.gz",
        ]
    ]
    manifest = {
        "classification_policy": target_selection.POLICY_ID,
        "classifications_changed": False,
        "external_source_sha256": {
            str(reference_data.local_path(source, "v23")): sha(
                reference_data.local_path(source, "v23")
            )
            for source in [
                "hpa_rna_consensus",
                "hpa_normal_tissue",
                "hpa_single_cell_tissue",
                "hpa_immune_cell",
                "hpa_immune_monaco",
            ]
        },
        "inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        "implementation_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [Path(__file__), ROOT / "scripts/audit_target_category_profiles.py"]
        },
        "outputs_sha256": {
            p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def render_html(rows):
    columns = [
        "protein",
        "representative_gene_id",
        "member_symbols",
        "member_gene_ids",
        "category",
        "classification",
        "published_nomination",
        "published_target_links",
        "normal_peak_tpm",
        "normal_location_tpm",
        "normal_source_tpm",
        "normal_source_url_tpm",
        "normal_peak_ntpm",
        "normal_location_ntpm",
        "normal_source_ntpm",
        "normal_source_url_ntpm",
        "normal_peak_cp10k",
        "normal_location_cp10k",
        "normal_source_cp10k",
        "normal_source_url_cp10k",
        "tumor_max_tpm",
        "cancer_evidence",
        "reproductive_evidence",
        "decision_reason",
    ]
    payload = rows[columns].astype(object).where(rows[columns].notna(), None).to_dict("records")
    for row in payload:
        for c in ["member_symbols", "member_gene_ids", "published_target_links"]:
            row[c] = json.loads(row[c])
    data = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace(
        "<", "\\u003c"
    )
    options = "".join(f"<option>{html.escape(c)}</option>" for c in target_selection.CATEGORIES)
    page = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Published targets and normal-expression peaks</title><style>body{font:15px system-ui;color:#173042;margin:2rem;line-height:1.5}h1{font-size:1.8rem}table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:.6rem;border-bottom:1px solid #d6dfe4;text-align:left;vertical-align:top}th{position:sticky;top:0;background:#eaf0f3}input,select,button{padding:.5rem;margin:.4rem}a{color:#08677b}.scroll{overflow:auto;max-height:75vh}details{max-width:30em}.note{max-width:95ch}</style><h1>Published targets and normal-expression peaks</h1><p class="note">Current draft: fivefold single-tumor rule. Normal peaks include Disposable tissues and exclude reproductive/developmental contexts. TPM, nTPM and CP10K stay separate. Missing values appear as —. <a href="../category-profiles/">Methods and downloads</a> · <a href="normal-top20/">Top 20 per category</a> · <a href="exclusive-sets/">Complete exclusive sets</a></p><label>Category <select id="category"><option value="">All</option>__OPTIONS__</select></label><label>Published <select id="published"><option value="yes" selected>Published target nominations</option><option value="">All proteins</option><option value="no">No target nomination</option></select></label><label>Search <input id="search" placeholder="Protein or Ensembl gene"></label><label>Sort normal peak <select id="unit"><option value="tpm">TPM</option><option value="ntpm">nTPM</option><option value="cp10k">CP10K</option></select></label><label>Limit <select id="limit"><option value="100">100</option><option value="20">20</option><option value="0">All matches</option></select></label><p id="count"></p><div class="scroll"><table><thead><tr><th>Protein / genes</th><th>Classification</th><th>Highest normal TPM</th><th>Highest normal nTPM</th><th>Highest normal CP10K</th><th>Highest tumor TPM</th><th>Evidence / papers</th></tr></thead><tbody id="body"></tbody></table></div><script>const rows=__DATA__;const el=id=>document.getElementById(id);function link(label,url){const a=document.createElement('a');a.textContent=label;a.href=url;return a}function cell(value){const td=document.createElement('td');td.textContent=value;return td}function fmt(v){return v===null?'—':Number(v.toPrecision(5)).toLocaleString()}function peak(r,u){const td=cell(fmt(r['normal_peak_'+u]));if(r['normal_location_'+u]){td.append(document.createElement('br'));td.append(document.createTextNode(r['normal_location_'+u]))}if(r['normal_source_'+u]){td.append(document.createElement('br'));td.append(link(r['normal_source_'+u],r['normal_source_url_'+u]))}return td}function render(){const q=el('search').value.toLowerCase().replace(/-/g,'');let hits=rows.filter(r=>(!el('category').value||r.category===el('category').value)&&(!el('published').value||r.published_nomination===(el('published').value==='yes'))&&(!q||[r.protein,...r.member_symbols,...r.member_gene_ids].join(' ').toLowerCase().replace(/-/g,'').includes(q)));const k='normal_peak_'+el('unit').value;hits.sort((a,b)=>(b[k]??-1)-(a[k]??-1)||a.protein.localeCompare(b.protein));const n=Number(el('limit').value);const shown=n?hits.slice(0,n):hits;el('count').textContent=hits.length+' matching proteins; showing '+shown.length;el('body').replaceChildren(...shown.map(r=>{const tr=document.createElement('tr'),name=cell('');name.append(link(r.protein,'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+r.representative_gene_id));name.append(document.createElement('br'));name.append(link('HPA','https://www.proteinatlas.org/'+r.representative_gene_id));name.append(document.createElement('br'));name.append(document.createTextNode(r.member_symbols.join(', ')));const info=cell(r.cancer_evidence+' / '+r.reproductive_evidence);info.append(document.createElement('br'));info.append(document.createTextNode(r.decision_reason));for(const [title,url] of r.published_target_links){info.append(document.createElement('br'));info.append(link(title,url))}tr.append(name,cell(r.classification),peak(r,'tpm'),peak(r,'ntpm'),peak(r,'cp10k'),cell(fmt(r.tumor_max_tpm)),info);return tr}));}for(const id of ['category','published','search','unit','limit'])el(id).addEventListener('input',render);render();</script></html>""".replace(
        "__OPTIONS__", options
    ).replace("__DATA__", data)
    (OUT / "report.html").write_text(page)


if __name__ == "__main__":
    build()
