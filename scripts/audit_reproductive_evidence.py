#!/usr/bin/env python3
"""Compare assay-specific expression and audit the requested 78/31 former CTAs.

This report never changes membership. TPM, HPA nTPM, donor UMI frequencies and
staining remain separate. Compact codes refer to a single report legend.
"""

from __future__ import annotations

import hashlib
import html
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import cta_screen, normal_tissues, reference_data, target_selection  # noqa: E402

OUT = ROOT / "docs/audits/antigen-evidence-review-20261118"
REBUILD = ROOT / "docs/audits/target-selection-20261110/rebuild"
RARE = ROOT / "docs/audits/normal-cell-atlas-20261110/rare-cell"
LEGEND = {
    "H": "Same-strand host overlap: resolve spliced/exonic locus attribution.",
    "P": "Paralog ambiguity: validate gene-unique sequence or a specific protein assay.",
    "S": "Selected finding qualifies in one catalog collection; verify independent study support.",
    "F": "Accepted reproductive enrichment remains below 10×.",
    "E": "Eye bulk comparison covers retina only; no matched HPA bulk or reliable eye staining.",
    "C": "Integrated retina atlas contributes; audit original-study/donor overlap before counting replication.",
    "I": "Reliable positive IHC in the implicated organ; inspect cell and antibody specificity.",
    "—": "Not measured or not attributable; never a measured zero.",
}
BRAIN = {
    "amygdala",
    "basal ganglia",
    "caudate",
    "cerebellum",
    "cerebral cortex",
    "choroid plexus",
    "hippocampal formation",
    "hippocampus",
    "hypothalamus",
    "medulla oblongata",
    "midbrain",
    "pons",
    "substantia nigra",
    "thalamus",
    "white matter",
}


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def organ(tissue):
    t = str(tissue).lower()
    if t in BRAIN:
        return "brain"
    if t in {"retina", "eye", "cornea", "conjunctiva"}:
        return "eye"
    return {"heart muscle": "heart"}.get(t, t)


def scope(tissue):
    # HPA's PBMC panel is peripheral blood mononuclear cells, not an organ alias.
    concept = "blood" if tissue == "pbmc" else tissue
    a = normal_tissues.resolve_normal_tissue(concept, species="human", stage="adult")
    category = normal_tissues.normal_tissue_category(a)
    if category == "unresolved":
        raise ValueError(f"Unresolved source tissue: {tissue}")
    return "reproductive" if category == "reproductive" else "normal"


def group_values(raw, lookup, genes, sizes, contexts):
    """Sum loci only within one source context, requiring every group member."""
    d = raw.copy()
    d["gene_id"] = d.Gene.map(lookup)
    d["entity_id"] = d.gene_id.map(genes)
    d = d[d.entity_id.notna()]
    if d.duplicated(["gene_id", *contexts]).any():
        raise ValueError("Duplicate source measurements cannot be added")
    good = d[d.nTPM.notna()]
    g = (
        good.groupby(["entity_id", *contexts], dropna=False)
        .agg(value=("nTPM", "sum"), measured=("gene_id", "nunique"))
        .reset_index()
    )
    g = g[g.measured.eq(g.entity_id.map(sizes))].copy()
    if not np.isfinite(g.value).all() or g.value.lt(0).any():
        raise ValueError("Invalid source expression")
    return g


def peak(rows):
    if rows.empty:
        return np.nan, "", ""
    r = rows.sort_values(
        ["value", "source_id", "source_context"], ascending=[False, True, True]
    ).iloc[0]
    return float(r.value), str(r.source_context), str(r.source_id)


def write(frame, name):
    frame.to_csv(
        OUT / name,
        index=False,
        compression={"method": "gzip", "mtime": 0} if name.endswith(".gz") else None,
    )


def main():
    inputs = {}

    def read(path, **kwargs):
        p = Path(path)
        inputs[str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)] = sha(p)
        return pd.read_csv(p, **kwargs)

    cohort = read(OUT / "requested-cohorts.csv")
    assert cohort.baseline_classification.value_counts().to_dict() == {
        "Cancer Antigen only": 78,
        "Other:Review": 31,
        "Other:Tissue": 15,
    }
    catalog = read(ROOT / "oncoref/data/target-selection.csv.gz").set_index("entity_id")
    members = {r.entity_id: json.loads(r.member_gene_ids) for r in cohort.itertuples()}
    genes = {g: e for e, ids in members.items() for g in ids}
    sizes = {e: len(ids) for e, ids in members.items()}
    assert all(json.loads(catalog.loc[e, "member_gene_ids"]) == ids for e, ids in members.items())
    mapping = read(
        ROOT / "docs/audits/cta-selection-20261028/hpa_rna_consensus-accepted-identities.csv.gz"
    )
    lookup = mapping.set_index("source_gene_id").canonical_gene_id
    source_genes = set(lookup[lookup.isin(genes)].index)
    normal = read(REBUILD / "normal-comparison.csv.gz", low_memory=False)
    normal = normal[
        normal.entity_id.isin(members)
        & normal.identity_resolved
        & normal.qc_pass
        & normal.value.notna()
    ].copy()
    repro = read(REBUILD / "reproductive-comparison.csv.gz", low_memory=False)
    repro = repro[
        repro.entity_id.isin(members) & repro.identity_resolved & repro.value.notna()
    ].copy()
    assert normal.unit.eq("TPM").all() and repro.unit.eq("TPM").all()
    normal["organ"] = normal.tissue.map(organ)
    panels = []
    for source, kind, contexts in [
        ("hpa_single_cell_tissue", "cell", ["Tissue", "Cluster", "Cell type"]),
        ("hpa_rna_consensus", "bulk", ["Tissue"]),
    ]:
        p = reference_data.local_path(source, "v23")
        pin = reference_data.REFERENCE_SOURCES[source]["pins"].get("v23")
        if pin is None and source == "hpa_rna_consensus":
            # Exact v23 bytes already bound by the released catalog input receipt.
            pin = {
                "sha256": "c49b5b33a076e3b9c1eb4f7d15d063ee936e07045e192ceccfeb9edcf1d90ddb",
                "bytes": 37949495,
            }
        digest = sha(p)
        assert digest == pin["sha256"] and p.stat().st_size == pin["bytes"]
        inputs[f"reference:{source}:v23"] = digest
        chunks = [d[d.Gene.isin(source_genes)] for d in pd.read_csv(p, sep="\t", chunksize=500000)]
        g = group_values(pd.concat(chunks, ignore_index=True), lookup, genes, sizes, contexts)
        g["source_id"] = source + ":v23"
        g["source_context"] = g.Tissue
        if kind == "cell":
            g["source_context"] += " / " + g["Cell type"] + " / " + g.Cluster
        g["scope"] = g.Tissue.map({t: scope(t) for t in g.Tissue.unique()})
        g["organ"] = g.Tissue.map(organ)
        g["kind"] = kind
        g["unit"] = "nTPM"
        panels.append(g)
    panels = pd.concat(panels, ignore_index=True)
    write(panels, "hpa-assay-measurements.csv.gz")

    findings = read(RARE / "findings.csv")
    findings = findings[findings.entity_id.isin(members)].copy()
    risk = read(ROOT / "oncoref/data/gene-mapping-risk.csv.gz").set_index("gene_id")
    donors = read(RARE / "donor-populations.csv.gz")
    donors = donors[donors.gene_id.isin(genes)].copy()
    donors["entity_id"] = donors.gene_id.map(genes)
    donors["multi_pct"] = 100 * donors.n_multi / donors.n_cells
    donors["qualifies"] = (
        donors.n_cells.ge(1000)
        & donors.n_multi.ge(10)
        & donors.multi_pct.ge(donors.tier.map({"vital": 1, "protected": 2}))
    )
    write(donors, "donor-evidence.csv.gz")
    write(findings, "all-organ-findings.csv")
    sources = read(ROOT / "docs/audits/normal-cell-atlas-20261110/source-registry.csv")
    write(
        sources[sources.collection_id.isin(donors.loc[donors.qualifies, "source_id"])],
        "census-sources.csv",
    )
    ihc_path = reference_data.local_path("hpa_normal_tissue", "v23")
    pin = reference_data.REFERENCE_SOURCES["hpa_normal_tissue"]["pins"]["v23"]
    assert sha(ihc_path) == pin["sha256"]
    inputs["reference:hpa_normal_tissue:v23"] = pin["sha256"]
    ihc = pd.read_csv(ihc_path, sep="\t", keep_default_na=False)
    ihc["gene_id"] = ihc.Gene.map(lookup)
    ihc["entity_id"] = ihc.gene_id.map(genes)
    ihc = ihc[ihc.entity_id.notna()].copy()
    ihc["organ"] = ihc.Tissue.map(organ)
    ihc["reliable"] = ihc.Reliability.isin(cta_screen.RELIABLE)
    write(ihc, "protein-staining.csv.gz")

    rows = []
    for r in cohort.itertuples():
        c = catalog.loc[r.entity_id]
        n = normal[normal.entity_id.eq(r.entity_id)]
        dev = repro[repro.entity_id.eq(r.entity_id) & repro.scope.eq("developmental")]
        p = panels[panels.entity_id.eq(r.entity_id)]
        row = {
            "entity_id": r.entity_id,
            "protein": r.protein_group_label,
            "baseline_classification": r.baseline_classification,
            "classification": c.classification,
            "N_reference_TPM": c.normal_reference_tpm,
            "D_TPM": c.disposable_max_tpm,
            "R_accepted_TPM": c.reproductive_supported_tpm,
            "R_fold": c.reproductive_fold,
            "T_fold": c.tumor_fold,
        }
        for name, data in {
            "N_bulk_TPM": n[n.scope.isin(["normal", "disposable"]) & n.cell_type.eq("bulk")],
            "R_bulk_TPM": n[n.scope.eq("reproductive") & n.cell_type.eq("bulk")],
            "embryo_observed_TPM": dev,
            "N_cell_nTPM": p[p.scope.eq("normal") & p.kind.eq("cell")],
            "R_cell_nTPM": p[p.scope.eq("reproductive") & p.kind.eq("cell")],
            "N_bulk_nTPM": p[p.scope.eq("normal") & p.kind.eq("bulk")],
            "R_bulk_nTPM": p[p.scope.eq("reproductive") & p.kind.eq("bulk")],
        }.items():
            row[name], row[name + "_context"], row[name + "_source"] = peak(data)
        d = donors[donors.entity_id.eq(r.entity_id) & donors.n_cells.ge(1000)]
        row["normal_UMI_peak_pct"] = float(d.multi_pct.max()) if len(d) else np.nan
        rows.append(row)
    comparison = pd.DataFrame(rows)
    write(comparison, "expression-comparison.csv")
    ca = comparison[comparison.baseline_classification.eq("Cancer Antigen only")]
    assert len(ca) == 78
    write(ca, "ca-only-78.csv")

    reviews = []
    for r in cohort[cohort.baseline_classification.eq("Other:Review")].itertuples():
        c = catalog.loc[r.entity_id]
        f = findings[findings.entity_id.eq(r.entity_id) & findings.status.eq("unresolved")].copy()
        f = f.sort_values(["replicated", "top_multi_molecule_pct"], ascending=False)
        assert len(f)
        top = f.iloc[0]
        k = risk.loc[top.gene_id]
        n = normal[
            normal.entity_id.eq(r.entity_id)
            & normal.organ.eq(top.tissue)
            & normal.cell_type.eq("bulk")
        ]
        p = panels[panels.entity_id.eq(r.entity_id) & panels.organ.eq(top.tissue)]
        protein = ihc[ihc.entity_id.eq(r.entity_id) & ihc.organ.eq(top.tissue) & ihc.reliable]
        positive = protein[protein.Level.isin(["Low", "Medium", "High"])]
        codes = []
        if k.host_overlapped:
            codes.append("H")
        if k.paralog_ambiguous:
            codes.append("P")
        if not top.replicated:
            codes.append("S")
        if c.reproductive_fold < 10:
            codes.append("F")
        if top.tissue == "eye":
            codes.append("E")
        implicated = donors[
            donors.entity_id.eq(r.entity_id) & donors.tissue.eq(top.tissue) & donors.qualifies
        ]
        if implicated.source_id.eq("4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e").any():
            codes.append("C")
        if len(positive):
            codes.append("I")
        row = {
            "entity_id": r.entity_id,
            "protein": r.protein_group_label,
            "gene_id": top.gene_id,
            "organ": top.tissue,
            "cell": top.top_cell_type,
            "tier": top.tier,
            "collections": int(top.studies_multi_molecule),
            "source_donor_labels": int(top.donors_multi_molecule),
            "peak_pct": top.top_multi_molecule_pct,
            "peak_cells": int(top.top_n_cells),
            "peak_multi_cells": int(top.top_n_multi),
            "replicated": bool(top.replicated),
            "codes": " ".join(codes),
            "shared_31mer_pct": 100 * k.shared_kmer_fraction,
            "shared_3prime_31mer_pct": 100 * k.shared_3prime_kmer_fraction,
            "host": k.same_strand_hosts,
            "partners": k.sequence_partners,
            "unresolved_organs": ";".join(sorted(f.tissue.unique())),
            "IHC_reliable_rows": len(protein),
            "IHC_positive_rows": len(positive),
            "R_fold": c.reproductive_fold,
            "T_fold": c.tumor_fold,
            "decision": "hold",
        }
        for name, data in {
            "organ_HPA_TPM": n[n.source_id.str.contains("hpa")],
            "organ_GTEx_TPM": n[n.source_id.str.contains("gtex")],
            "organ_cell_nTPM": p[p.kind.eq("cell")],
            "organ_bulk_nTPM": p[p.kind.eq("bulk")],
        }.items():
            row[name], row[name + "_context"], row[name + "_source"] = peak(data)
        reviews.append(row)
    review = pd.DataFrame(reviews)
    assert len(review) == 31
    write(review, "review-31.csv")
    # Publication adjudications live in the shared ledger, not a report-only file.
    ledger = ROOT / "oncoref/data/cta-literature-followup.json"
    inputs[str(ledger.relative_to(ROOT))] = sha(ledger)
    render(ca, review)
    receipt = {
        "schema": "oncoref.antigen-evidence-review.v1",
        "policy": target_selection.POLICY_ID,
        "cohort_release": "1.42.0",
        "cohort_commit": "d3f58e94f5b18d195085af3f8e54071b1a65ccf3",
        "legend": LEGEND,
        "inputs": inputs,
        "implementation_sha256": sha(__file__),
        "outputs": {
            p.name: sha(p)
            for p in sorted(OUT.iterdir())
            if p.is_file() and p.name not in {"manifest.json", "requested-cohorts.csv"}
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {
                "CA_only_compared": len(ca),
                "holds_reviewed": len(review),
                "review_codes": review.codes.value_counts().to_dict(),
            }
        )
    )


def render(ca, review):
    def number(value):
        return "—" if pd.isna(value) else f"{value:.3g}"

    def measured(row, key):
        v = number(row[key])
        return v if v == "—" else v + " · " + str(row[key + "_context"]).replace("|", "/")

    lines = [
        "# Reproductive evidence and normal-cell holds",
        "",
        "The requested cohorts are frozen at release 1.42.0: 78 CA-only proteins and 31 Other:Review proteins with tumor ≥10×. Current measurements are regenerated from pinned sources; this report does not change admission.",
        "",
        "Bulk TPM, HPA cell-cluster nTPM and donor UMI frequencies are separate assays. Maxima include all measured non-reproductive normal tissues on the normal side and all reproductive tissues on the reproductive side. The normal maximum includes Disposable tissues; their separate D value is also retained in the CSV. A cell-cluster maximum is descriptive and cannot be divided by bulk TPM to qualify RA. Whole-organ bulk does not measure every cell population. Embryo TPM is a separate within-embryo/lineage mean; its observed maximum can fail cell/donor support.",
        "",
        "Sources: HPA v23 tissue-resolved single-cell clusters and bulk consensus, native-TPM HPA individual tissues and bulk/brain summaries, GTEx v11 and HPA-hosted GTEx, Petropoulos/Tyser embryo observations, and the pinned Census donor UMI panel. The native normal reference also includes sorted immune panels. No single-cell value is silently substituted for missing bulk coverage.",
        "",
        "## Codes",
        "",
        *[f"- **{code}**: {meaning}" for code, meaning in LEGEND.items()],
        "",
        "## The 78 CA-only proteins",
        "",
        "[Full numeric table, exact contexts and sources](ca-only-78.csv). The accepted RA fold uses the existing native-TPM witness/support rules, not the largest value across these assays.",
        "",
        "| Protein | Normal bulk TPM | Reproductive bulk TPM | Normal cells nTPM | Reproductive cells nTPM | Observed embryo TPM | Accepted R fold |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in ca.to_dict("records"):
        fields = [
            row["protein"],
            *[
                measured(row, k)
                for k in (
                    "N_bulk_TPM",
                    "R_bulk_TPM",
                    "N_cell_nTPM",
                    "R_cell_nTPM",
                    "embryo_observed_TPM",
                )
            ],
            number(row["R_fold"]) + "×",
        ]
        lines.append("| " + " | ".join(fields) + " |")
    lines += [
        "",
        "## Review of the 31 holds",
        "",
        "Eye and retinal findings remain Vital-tissue concerns. All 31 remain excluded from antigen sets. None is rescued by low bulk RNA or a mapping-risk flag. Established gene-specific normal-cell expression would move a hold to Unsafe; the current review does not establish that attribution.",
        "",
        "### What the deeper review found",
        "",
        "- **Ocular evidence:** HPA's retina values under its GTEx label are actually [EyeGEx RNA processed by HPA](https://v23.proteinatlas.org/about/assays%2Bannotation). The older rare-cell confirmation summary omitted eye; the table below includes retinal bulk, including RAD21L1 (1.2 TPM), RGPD3 (1.3) and RGPD4 (1.4). These are additional concerns, not proof that the same cell or exact locus produced both signals. There is no separate HPA retinal bulk measurement here. Retinal RNA cannot assess corneal, conjunctival or optic-nerve populations, and duplicated HPA exports of EyeGEx are one source.",
        "- **Replication:** collection IDs are not guaranteed independent studies. The [integrated retinal atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC13050529/) combines public studies and new samples. Its final publication postdates the pinned Census snapshot, so original-study and donor mapping must be checked against that snapshot. Counts below are collections and source-specific donor labels, not verified unique patients. Possible overlap does not dismiss the signals.",
        "- **PSG8:** the strongest ocular observation reaches 46% multi-molecule cells; the paralog attribution question remains consequential. [Kurji et al.](https://doi.org/10.1167/iovs.09-3622) also report a PSG8-labelled microarray signal in cultured RPE (Table 2). The probe has not been checked for PSG-family specificity and the reported validation panel did not validate PSG8. This is additional RNA evidence, not exact-protein confirmation.",
        "- **Host-overlap hypotheses:** CAGE1, PPP3R2, RAD21L1 and the other H-coded genes need locus-specific attribution. Same-strand overlap does not demonstrate host pre-mRNA contamination. The previous bulk confirmation score was the weaker source maximum, not the highest normal expression: for example, brain HPA reaches 0.5 TPM for CAGE1 and 2.7 for PPP3R2, despite much lower GTEx values. Negative bulk cannot rule out a rare essential population.",
        "- **Single-collection findings:** DCAF8L2 (16% retinal bipolar cells), FBXO47 (6.8% keratocytes) and NKX2-6 (37% atrial cardiomyocytes) lack a mapping-risk flag but still need independent confirmation. Their frequency warrants holding them; it does not justify calling them artifacts.",
        "- **DAZL and CTAG2:** [Zhang et al.](https://doi.org/10.1186/s12885-020-07155-y) provide DAZL glioblastoma-expression support, but their normal astrocyte comparator does not resolve the ocular astrocyte signal. CTAG2 has low but nonzero lung bulk measurements and remains below reproductive 10× even if its review hold were resolved. Both newly reviewed papers are recorded centrally with exact gene IDs and permitted uses.",
        "",
        "No selected organ has reliable positive HPA staining in this extraction. Usually the assay is missing; only measured negative rows count as negative staining. Neither missing staining nor a negative bulk measurement clears a cell-specific concern. Use Unsafe|Other for conservative peptide subtraction so unresolved proteins remain in the exclusion background.",
        "",
        "| Protein | Organ / strongest population | Collections / donor labels | Peak ≥2-UMI cells | Same-organ HPA / GTEx–EyeGEx TPM | Same-organ cell nTPM | Codes |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in review.to_dict("records"):
        fields = [
            row["protein"],
            row["organ"] + " / " + row["cell"],
            f"{row['collections']} / {row['source_donor_labels']}",
            number(row["peak_pct"]) + "%",
            number(row["organ_HPA_TPM"]) + " / " + number(row["organ_GTEx_TPM"]),
            number(row["organ_cell_nTPM"]),
            row["codes"],
        ]
        lines.append("| " + " | ".join(fields) + " |")
    lines += [
        "",
        "[31-gene review details](review-31.csv) · [Every organ finding](all-organ-findings.csv) · [Donor observations](donor-evidence.csv.gz) · [Primary source registry](census-sources.csv) · [HPA assay measurements](hpa-assay-measurements.csv.gz) · [Protein staining](protein-staining.csv.gz) · [All 124 comparisons](expression-comparison.csv) · [Receipt](manifest.json)",
        "",
    ]
    (OUT / "index.md").write_text("\n".join(lines))
    import markdown

    body = markdown.markdown("\n".join(lines), extensions=["tables"])
    title = "Reproductive evidence and normal-cell holds"
    (OUT / "report.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>body{{font:15px/1.5 system-ui;margin:2rem;color:#21343d}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:.6rem;text-align:left;border-bottom:1px solid #ddd;vertical-align:top}}th{{background:#eef3f5;position:sticky;top:0}}p{{max-width:95ch}}a{{color:#126487}}</style><main>"
        + body
        + "</main></html>"
    )


if __name__ == "__main__":
    main()
