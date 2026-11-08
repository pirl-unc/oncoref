#!/usr/bin/env python3
"""Rebuild the named-candidate evidence audit from pinned, native-unit inputs."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import (  # noqa: E402
    cancer_type_registry,
    cta_decision_flow,
    cta_representatives,
    gene_knowledge,
)
from oncoref import cta_selection_policy as policy  # noqa: E402
from oncoref.cta_funnel import EvidenceLayer  # noqa: E402

DEFAULT = ROOT / "docs/audits/cta-candidate-evidence-20261108"
THRESHOLDS = (0.1, 1, 2, 5, 10, 100)


def threshold_key(value):
    return "ge_" + str(value).replace(".", "_")


def count_samples(values):
    """Native TPM thresholds; zero stays measured, missing stays unmeasured."""
    finite = values.notna() & np.isfinite(values) & values.ge(0)
    if not finite.all():
        raise ValueError("Eligible tumor measurements must be finite, nonnegative and measured")
    return {
        "measured_samples": len(values),
        "gt_0": int(values.gt(0).sum()),
        **{threshold_key(t): int(values.ge(t).sum()) for t in THRESHOLDS},
        "maximum_tpm": values.max() if len(values) else None,
    }


def tumor_tables(raw, qc, candidates, coverage):
    """Count exact loci, retaining aliases and known sample namespaces.

    Cohort views overlap and must not be summed. Quality failures and non-TPM
    sources remain in the input ledger but cannot meet an absolute TPM floor.
    """
    if raw.duplicated(["Ensembl_Gene_ID", "cohort", "sample"]).any():
        raise ValueError("Duplicate locus/sample within a cohort")
    joined = raw.merge(
        qc,
        left_on=["cohort", "sample"],
        right_on=["cancer_code", "sample_id"],
        how="left",
        validate="many_to_one",
    )
    joined["eligible"] = joined.recommended_for_absolute_tpm_floor.eq(True)
    if (joined.eligible & ~joined.linear_tpm_comparable.eq(True)).any():
        raise ValueError("An eligible sample lacks a linear RNA-seq TPM scale")
    keys = ["Ensembl_Gene_ID", "namespace", "sample"]
    # An alias may change across matrices, but a shared sample vector must not.
    if joined.groupby(keys).tpm.nunique(dropna=False).gt(1).any():
        raise ValueError("Overlapping cohort views disagree on the same locus/sample")
    if joined.groupby(keys).eligible.nunique().gt(1).any():
        raise ValueError("Overlapping cohort views disagree on sample quality")
    distinct = joined.drop_duplicates(keys).copy()
    by_cohort = dict(iter(joined.groupby(["Ensembl_Gene_ID", "cohort"])))
    counts, cohort_rows = [], []
    for item in candidates:
        gid = item["gene_id"]
        all_rows = distinct[distinct.Ensembl_Gene_ID.eq(gid)]
        eligible = all_rows[all_rows.eligible]
        counts.append(
            {
                "gene_id": gid,
                "protein": item["protein"],
                "observed_sample_records": len(all_rows),
                "excluded_sample_records": int((~all_rows.eligible).sum()),
                **count_samples(eligible.tpm),
            }
        )
        for row in coverage[coverage.gene_id.eq(gid)].itertuples():
            measured = by_cohort.get((gid, row.cohort), joined.iloc[:0])
            passing = measured[measured.eligible]
            if row.status == "gene_absent" and len(measured):
                raise ValueError("Gene-coverage record contradicts observed rows")
            cohort_rows.append(
                {
                    "gene_id": gid,
                    "protein": item["protein"],
                    "cohort": row.cohort,
                    "coverage": row.status,
                    "matrix_asset": row.matrix_asset,
                    "observed_samples": len(measured),
                    "excluded_samples": len(measured) - len(passing),
                    **(
                        count_samples(passing.tpm)
                        if len(passing)
                        else {
                            "measured_samples": 0,
                            "gt_0": None,
                            **{threshold_key(t): None for t in THRESHOLDS},
                            "maximum_tpm": None,
                        }
                    ),
                }
            )
    names = {r["gene_id"]: r["protein"] for r in candidates}
    distinct["protein"] = distinct.Ensembl_Gene_ID.map(names)
    return pd.DataFrame(counts), pd.DataFrame(cohort_rows), distinct


def replay_updates(baseline):
    """Replay accepted curation without rebuilding or editing released assets."""
    layers = [
        EvidenceLayer(c[5:], baseline[c], baseline["unsafe:" + c[5:]])
        for c in baseline
        if c.startswith("hold:")
    ]
    expected_unsafe = baseline.category.eq("Unsafe")

    def unchanged_overlap(mask):
        if not mask.equals(expected_unsafe):
            raise ValueError("Unsafe membership changed; recompute peptide containment")
        return baseline.unsafe_sequence_overlap_pct / 100

    reproduced = policy.classify(
        baseline,
        layers,
        entity_column="protein_equivalence_group_id",
        unsafe_overlap=unchanged_overlap,
        legacy_subtypes=True,
    ).classifications
    if not reproduced.classification.equals(baseline.classification):
        raise ValueError("Stored gates and containment do not reproduce the released decisions")
    updated = gene_knowledge.annotate(
        baseline, evidence=gene_knowledge.findings(include_updates=True)
    )
    updated_layers = [
        EvidenceLayer(
            name,
            hold | updated.knowledge_normal_hold | updated.knowledge_normal_exclusion,
            unsafe | updated.knowledge_normal_exclusion,
        )
        if name == "Reviewed normal-cell knowledge"
        else EvidenceLayer(name, hold, unsafe)
        for name, hold, unsafe in layers
    ]
    result = policy.classify(
        updated,
        updated_layers,
        entity_column="protein_equivalence_group_id",
        unsafe_overlap=unchanged_overlap,
    )
    return result


def assay_coverage(inputs, candidates):
    rows = []
    uniprot = {
        e["primaryAccession"]: e
        for e in json.loads(gzip.decompress((inputs / "uniprot-20261108.json.gz").read_bytes()))
    }
    for item in candidates:
        if uniprot[item["uniprot"]]["proteinExistence"] != item["protein_existence"]:
            raise ValueError("Protein-existence annotation disagrees with its source response")
        entry = ET.fromstring(
            gzip.decompress((inputs / f"{item['protein']}.xml.gz").read_bytes())
        ).find("entry")
        antibodies = [a.attrib["id"] for a in entry.findall("antibody")]
        ms = entry.find("proteinEvidence/evidence[@source='MS']")
        rows.append(
            {
                **item,
                "hpa_ihc": "Uncertain shared antibody" if antibodies else "Not assayed",
                "hpa_antibodies": ";".join(antibodies),
                "hpa_ms": ms.attrib["evidence"] if ms is not None else "Not reported",
            }
        )
    return pd.DataFrame(rows)


def uniform_decisions(out, revised):
    """Export every assessed protein and its member loci using the same flow."""
    contract = cta_decision_flow.definition()
    (out / "decision-flow.json").write_text(json.dumps(contract, indent=2) + "\n")
    columns = [
        "protein_equivalence_group_id",
        "protein_group_label",
        "member_gene_ids",
        "category",
        "subtype",
        "classification",
        "assessment_status",
        "source_classification",
        "decision_flow_id",
        "policy",
        "exclusion_reasons",
        "Unsafe_reasons",
        "required_coverage",
        "unresolved_evidence",
        "positive_context",
        "cancer_association",
        "translation_support_eligible",
        "translation_qualified",
        "representative_selection_approved",
        "representative_transcripts",
        "knowledge_finding_ids",
        "unsafe_sequence_overlap_pct",
    ]
    proteins = revised[columns].copy()
    write_csv(proteins, out / "uniform-protein-decisions.csv.gz")
    genes = proteins.assign(gene_id=proteins.member_gene_ids.map(json.loads)).explode("gene_id")
    if genes.gene_id.isna().any() or not genes.gene_id.is_unique:
        raise ValueError("Every member locus must have exactly one uniform decision")
    write_csv(genes, out / "uniform-gene-decisions.csv.gz")
    selected = revised[revised.protein_group_label.isin([f"SSU72L{i}" for i in range(1, 6)])]
    proof = selected[
        [
            "protein_group_label",
            "representative_gene_id",
            "classification",
            "developmental_mean_tpm",
            "cancer_association",
            "required_coverage",
            "translation_support_eligible",
            "translation_qualified",
            "representative_selection_approved",
            "assessment_status",
            "exclusion_reasons",
        ]
    ].merge(
        cta_representatives.records(),
        left_on="representative_gene_id",
        right_on="gene_id",
        validate="one_to_one",
    )
    write_csv(proof, out / "ssu72-admission.csv")
    fig, ax = plt.subplots(figsize=(12, 7.5), layout="constrained")
    ax.set(xlim=(0, 12), ylim=(0, 7.5))
    ax.axis("off")
    ax.text(0.2, 7.1, "One decision for every assessed human protein", fontsize=18, weight="bold")
    colors = {"Unsafe": "#943c41", "CTA": "#126c63"}
    for step in contract["steps"]:
        y = 6.2 - (step["order"] - 1) * 1.15
        question = step["question"].replace(
            "Cancer, reproductive/developmental or concentration",
            "Cancer, reproductive/developmental\nor concentration",
        )
        ax.text(0.3, y, f"{step['order']}. {question}", va="center", fontsize=12)
        ax.annotate(
            "", xy=(8.9, y), xytext=(8.0, y), arrowprops={"arrowstyle": "->", "color": "#68737b"}
        )
        ax.text(
            9.1,
            y,
            step["if_yes"],
            va="center",
            fontsize=13,
            weight="bold",
            color=colors.get(step["if_yes"], "#76542c"),
        )
        if step["order"] < 5:
            ax.annotate(
                "",
                xy=(0.45, y - 0.83),
                xytext=(0.45, y - 0.34),
                arrowprops={"arrowstyle": "->", "color": "#68737b"},
            )
            ax.text(0.65, y - 0.6, "Otherwise", fontsize=9, color="#68737b")
    ax.text(
        0.3,
        0.25,
        "Paper coverage, protein detection, biological program and tumor prevalence remain separate evidence fields.",
        fontsize=10,
    )
    for suffix in ("svg", "png"):
        fig.savefig(
            out / f"decision-flow.{suffix}",
            dpi=150,
            metadata={"Date": None} if suffix == "svg" else {},
        )
    plt.close(fig)
    return {
        "proteins": len(proteins),
        "member_loci": len(genes),
        "other_subtypes": proteins[proteins.category.eq("Other")].subtype.value_counts().to_dict(),
        "flow_id": contract["flow_id"],
    }


def figures(out, counts, samples):
    plt.rcParams.update({"font.size": 11, "svg.hashsalt": "cta-candidate-evidence"})
    order = counts.sort_values("ge_1", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 5.4), layout="constrained")
    pos = np.arange(len(order))
    for offset, threshold, color in [
        (-0.22, 1, "#167d8d"),
        (0, 5, "#d77a31"),
        (0.22, 10, "#5a466b"),
    ]:
        values = order[threshold_key(threshold)]
        ax.barh(pos + offset, values, height=0.21, label=f"≥{threshold} TPM", color=color)
    ax.set_yticks(pos, order.protein)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlabel("Samples passing quality checks (symlog scale; linear from 0–10)")
    ax.set_title("Tumor RNA is usually sparse in these candidates", loc="left", weight="bold")
    ax.legend(frameon=False)
    ax.spines[["right", "top"]].set_visible(False)
    for ext in ["svg", "png"]:
        fig.savefig(
            out / f"tumor-prevalence.{ext}",
            dpi=160,
            metadata={"Date": None} if ext == "svg" else None,
        )
    plt.close(fig)
    x = samples[samples.protein.eq("SSU72L2") & samples.eligible & samples.tpm.ge(1)]
    x = x.sort_values("tpm", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 4.8), layout="constrained")
    ax.scatter(x.index + 1, x.tpm, color="#167d8d", s=32)
    ax.set_yscale("log")
    ax.set_xlabel("Sample rank among the 50 with ≥1 TPM")
    ax.set_ylabel("Native RNA-seq TPM (log scale)")
    ax.set_title("SSU72L2: one very high sample, few strong positives", loc="left", weight="bold")
    ax.annotate(
        "Embryonal rhabdomyosarcoma\nTHR24_2154_S01 · 990.91 TPM",
        (1, x.tpm.iloc[0]),
        xytext=(10, 250),
        arrowprops={"arrowstyle": "-", "color": "#777"},
    )
    for threshold in [5, 10]:
        ax.axhline(threshold, color="#bbb", linewidth=0.7, linestyle="--")
    ax.spines[["right", "top"]].set_visible(False)
    for ext in ["svg", "png"]:
        fig.savefig(
            out / f"ssu72l2-samples.{ext}",
            dpi=160,
            metadata={"Date": None} if ext == "svg" else None,
        )
    plt.close(fig)


def build(out=DEFAULT):
    inputs = out / "inputs"
    candidates = json.loads((inputs / "candidates.json").read_text())

    def read(name):
        return pd.read_csv(inputs / name, low_memory=False)

    counts, cohorts, samples = tumor_tables(
        read("tumor-samples.csv.gz"),
        read("sample-qc.csv.gz"),
        candidates,
        read("tumor-coverage.csv"),
    )
    cancer_names = cancer_type_registry().set_index("code")["name"]
    cohorts["cancer_type"] = cohorts.cohort.map(cancer_names)
    samples["cancer_type"] = samples.cohort.map(cancer_names)
    if cohorts.cancer_type.isna().any():
        raise ValueError("A cohort lacks its cancer-type name")
    coverage = assay_coverage(inputs, candidates)
    candidate_table = coverage.merge(counts, on=["gene_id", "protein"], validate="one_to_one")
    baseline_path = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
    baseline = pd.read_csv(baseline_path, low_memory=False)
    result = replay_updates(baseline)
    revised = result.classifications
    uniform_summary = uniform_decisions(out, revised)
    lookup = revised.set_index("representative_gene_id")
    candidate_table["released_category"] = candidate_table.gene_id.map(
        baseline.set_index("representative_gene_id").category
    )
    candidate_table["next_classification"] = candidate_table.gene_id.map(lookup.classification)
    candidate_table["next_assessment_status"] = candidate_table.gene_id.map(
        lookup.assessment_status
    )
    changed = revised[revised.category.ne(baseline.category)].copy()
    changed["released_classification"] = baseline.loc[changed.index, "classification"]
    change_cols = [
        "protein_group_label",
        "representative_gene_id",
        "released_classification",
        "classification",
        "assessment_status",
        "knowledge_finding_ids",
    ]
    write_csv(changed[change_cols], out / "classification-changes.csv")
    write_csv(candidate_table, out / "candidates.csv")
    write_csv(cohorts, out / "tumor-cohort-prevalence.csv")
    sample_cols = [
        "protein",
        "Ensembl_Gene_ID",
        "Symbol",
        "namespace",
        "sample",
        "cohort",
        "cancer_type",
        "tpm",
        "eligible",
        "sample_qc_status",
        "sample_qc_reasons",
        "matrix_version",
    ]
    write_csv(samples[sample_cols], out / "tumor-samples.csv.gz")
    ssu = samples[samples.protein.eq("SSU72L2") & samples.eligible]
    write_csv(
        ssu[ssu.tpm.ge(1)][sample_cols].sort_values("tpm", ascending=False),
        out / "ssu72l2-positive-samples.csv",
    )
    write_csv(cohorts[cohorts.protein.eq("SSU72L2")], out / "ssu72l2-cohorts.csv")
    # Missing IHC and MS are explicit strings, never numeric negatives.
    summary = {
        "released_categories": baseline.category.value_counts().to_dict(),
        "proposed_categories": revised.category.value_counts().to_dict(),
        "uniform_decisions": uniform_summary,
        "status": "curation_replay_not_a_release",
        "ssu72l2": counts[counts.protein.eq("SSU72L2")].iloc[0].to_dict(),
        "matrix_cohorts": int(cohorts.cohort.nunique()),
        "sample_count_unit": "source namespace/sample; overlapping views deduplicated; not patients",
        "translation_rule": "Qualifying reference coding annotation remains sufficient; no experimental-protein requirement was added.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures(out, counts, samples)
    template = (ROOT / "scripts/cta_evidence_assets/named-candidate-report.md").read_text()
    top = ssu.nlargest(6, "tpm")
    names = {
        "SARC_RMS_ERMS": "Embryonal rhabdomyosarcoma",
        "SARC_ANGIO": "Angiosarcoma",
        "SKCM": "Melanoma",
        "NEC_LUNG_LARGECELL": "Large-cell neuroendocrine lung carcinoma",
        "WILMS": "Wilms tumor",
        "CESC": "Cervical cancer",
    }
    top_rows = [
        "| Cancer type | Sample | TPM | Samples ≥5 TPM / measured in cohort |",
        "|---|---|---:|---:|",
    ]
    for r in top.itertuples():
        c = cohorts[cohorts.protein.eq("SSU72L2") & cohorts.cohort.eq(r.cohort)].iloc[0]
        top_rows.append(
            f"| {names[r.cohort]} | {r.sample} | {r.tpm:.2f} | {int(c.ge_5)}/{c.measured_samples} |"
        )
    summary_rows = [
        "| Protein | Measured samples | ≥1 TPM | ≥5 TPM | ≥10 TPM | Maximum TPM | Next decision |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in candidate_table.itertuples():
        summary_rows.append(
            f"| [{r.protein}]({r.ensembl_url}) | {r.measured_samples:,} | {r.ge_1} | {r.ge_5} | {r.ge_10} | {r.maximum_tpm:.2f} | {r.next_classification} |"
        )
    report = template.replace("{{TOP_SAMPLES}}", "\n".join(top_rows)).replace(
        "{{COUNTS}}", "\n".join(summary_rows)
    )
    (out / "index.md").write_text(report)
    import markdown

    style = """body{font:16px/1.55 system-ui;color:#22313a;margin:30px auto;max-width:1180px;padding:0 20px}a{color:#126d80}h1{font-size:30px}h2{margin-top:2em}table{border-collapse:collapse;font-size:14px;display:block;overflow:auto}td,th{padding:9px 12px;border-bottom:1px solid #dce2e6;min-width:100px;text-align:left}th{background:#edf4f6}img{max-width:100%;height:auto}pre{overflow:auto;background:#f3f5f7;padding:16px}"""
    (out / "report.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Ten CTA candidates: evidence review</title><style>'
        + style
        + "</style><main>"
        + markdown.markdown(report, extensions=["tables", "fenced_code"])
        + "</main></html>"
    )
    # Standalone searchable tables; source links remain usable offline.
    candidate_view = candidate_table[
        [
            "protein",
            "next_classification",
            "protein_existence",
            "hpa_ihc",
            "hpa_ms",
            "measured_samples",
            "ge_1",
            "ge_5",
            "ge_10",
            "maximum_tpm",
        ]
    ].copy()
    candidate_view["protein"] = [
        f'<a href="{html.escape(r.ensembl_url)}">{html.escape(r.protein)}</a>'
        f' · <a href="{html.escape(r.hpa_url)}">HPA</a>'
        f' · <a href="{html.escape(r.uniprot_url)}">UniProt</a>'
        for r in candidate_table.itertuples()
    ]
    candidate_view.columns = [
        "Protein / sources",
        "Next decision",
        "Protein annotation",
        "HPA IHC",
        "HPA MS",
        "Measured samples",
        "≥1 TPM",
        "≥5 TPM",
        "≥10 TPM",
        "Maximum TPM",
    ]
    cohort_view = cohorts[cohorts.protein.eq("SSU72L2")][
        [
            "cancer_type",
            "cohort",
            "coverage",
            "measured_samples",
            "excluded_samples",
            "ge_1",
            "ge_2",
            "ge_5",
            "ge_10",
            "maximum_tpm",
        ]
    ].copy()
    cohort_view["coverage"] = cohort_view.coverage.replace(
        {"gene_absent": "Gene not measured", "measured": "Gene measured"}
    )
    cohort_view.columns = [
        "Cancer type",
        "Code",
        "Source coverage",
        "Measured samples",
        "Excluded samples",
        "≥1 TPM",
        "≥2 TPM",
        "≥5 TPM",
        "≥10 TPM",
        "Maximum TPM",
    ]
    tables = []
    for title, frame in [
        ("Candidate measurements", candidate_view),
        ("SSU72L2 cohorts", cohort_view),
    ]:
        tables.append(
            f"<h2>{html.escape(title)}</h2>"
            + frame.to_html(
                index=False, na_rep="—", escape=False, border=0, float_format=lambda v: f"{v:,.2f}"
            )
        )
    (out / "table.html").write_text(
        """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CTA candidate evidence</title><style>
body{font:15px/1.5 system-ui;margin:28px;color:#22313a}a{color:#126d80}table{border-collapse:collapse;font-size:13px}td,th{padding:8px 12px;border-bottom:1px solid #dce2e6;max-width:330px;overflow-wrap:anywhere}th{position:sticky;top:0;background:#edf4f6;text-align:left}input{font:inherit;padding:9px;width:330px}main{overflow:auto}h1{font-size:28px}
</style><h1>Ten CTA candidates: evidence and tumor prevalence</h1><p>Release 1.35.0 still labels all ten CTA. The proposed curation places PRAMEF18 and TRIM75 on review hold. A dash means no eligible measurement, not zero. Cohort views overlap; do not add their denominators.</p><p><a href="report.html">Full rationale and linked papers</a> · <a href="candidates.csv">Candidate CSV</a> · <a href="ssu72l2-cohorts.csv">SSU72L2 cohort CSV</a></p><label>Filter rows <input id="filter" type="search" placeholder="Protein, cancer type, source or status"></label><main>"""
        + "\n".join(tables)
        + """</main><script>document.getElementById('filter').addEventListener('input',e=>{const q=e.target.value.toLowerCase();document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></html>"""
    )
    source_paths = [
        Path(__file__),
        ROOT / "scripts/collect_cta_candidate_tumors.py",
        ROOT / "oncoref/gene_knowledge.py",
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "oncoref/cta_decision_flow.py",
        ROOT / "oncoref/cta_representatives.py",
        ROOT / "oncoref/data/gene-knowledge.json",
        ROOT / "oncoref/data/gene-knowledge-updates.json",
        ROOT / "oncoref/data/cancer-type-registry.csv",
        ROOT / "scripts/cta_evidence_assets/named-candidate-report.md",
        baseline_path,
        *sorted(inputs.glob("*")),
    ]
    manifest = {
        "inputs": input_hashes(source_paths),
        "outputs": {
            p.name: sha256(p)
            for p in sorted(out.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    build(parser.parse_args().output)
