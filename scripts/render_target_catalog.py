#!/usr/bin/env python3
"""Create the filterable rebuilt catalog, migrations and count figure."""

from __future__ import annotations

import hashlib
import html
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import render_curation_summary  # noqa: E402
import render_normal_rna_sources  # noqa: E402

from oncoref import cta_literature, gene_knowledge, target_selection, tissue_tiers  # noqa: E402

OUT = ROOT / "docs/audits/target-selection-20261110/rebuild"
DRAFT = ROOT / "docs/audits/target-selection-20261109/rebuild/protein-decisions.csv.gz"


def funnel(data):
    """Sequential removal by the shared screen; every protein lands in one step."""
    new_reasons = {"Established rare normal-cell RNA", "Published normal-tissue expression"}
    only_new = data.Unsafe_reasons.fillna("").map(
        lambda v: bool(v) and set(v.split(" | ")) <= new_reasons
    )
    unsafe = data.category.eq("Unsafe")
    steps = [
        (
            "Bulk, cluster, immune, protein or reviewed normal exclusion",
            unsafe & ~only_new,
            "Unsafe",
        ),
        (
            "Established rare-cell RNA or published normal expression only",
            unsafe & only_new,
            "Unsafe",
        ),
        (
            render_curation_summary.OTHER[0][0],
            data.classification.eq("Other:Review"),
            "Other:Review",
        ),
        (
            render_curation_summary.OTHER[1][0],
            data.classification.eq("Other:Sequence"),
            "Other:Sequence",
        ),
        (
            "No qualifying expression, target maximum or Disposable cap",
            data.classification.eq("Other:Expression"),
            "Other:Expression",
        ),
        (
            "Cancer enrichment with insufficient home-tissue coverage and no qualifying RA witness",
            data.classification.eq("Other:Unscreened tissue"),
            "Other:Unscreened tissue",
        ),
    ]
    rows, remaining = [], len(data)
    for label, mask, result in steps:
        remaining -= int(mask.sum())
        rows.append(
            {"step": label, "result": result, "removed": int(mask.sum()), "remaining": remaining}
        )
    if remaining != int(data.is_target_candidate.sum()):
        raise ValueError("Funnel does not end at the Candidate antigen proteins")
    return pd.DataFrame(rows)


def funnel_figure(table, data, counts):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    labels = ["All assessed proteins", *table.step, "Candidate antigen proteins"]
    values = [len(data), *table.remaining, int(data.is_target_candidate.sum())]
    colors = [
        "#6f7f89",
        "#ac5552",
        "#ac5552",
        "#b88a3c",
        "#8d7aa8",
        "#9aabb5",
        "#9aabb5",
        "#3f8f6b",
    ]
    ax.barh(labels[::-1], values[::-1], color=colors[::-1])
    removed = table.set_index("step").removed
    for i, (label, value) in enumerate(zip(labels[::-1], values[::-1])):
        note = f"{value:,}" + (f"  (−{removed[label]:,})" if label in removed.index else "")
        ax.text(value + len(data) * 0.01, i, note, va="center", fontsize=10)
    ax.set_xlim(0, len(data) * 1.2)
    ax.set_xlabel("Proteins remaining after each filter")
    ax.set_title(
        f"Curation funnel: {counts['CTA']} CTA, {counts['Cancer Antigen only']} Cancer Antigen "
        f"only, {counts['Reproductive Antigen only']} Reproductive Antigen only",
        loc="left",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "funnel.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def build():
    data = target_selection.results()
    meta = target_selection.metadata()
    summary = meta["summary"]
    prior_receipt = json.loads((OUT / "release-1.39-input.json").read_text())
    prior_path = OUT / prior_receipt["projection"]
    if hashlib.sha256(prior_path.read_bytes()).hexdigest() != prior_receipt["sha256"]:
        raise ValueError("Changed 1.39 membership baseline")
    previous = pd.read_csv(prior_path).set_index("entity_id")
    data["previous_classification"] = data.entity_id.map(previous.classification)
    changes = data[data.classification.ne(data.previous_classification)]
    changes.to_csv(OUT / "changes-from-1.39.csv", index=False)
    changes_summary = (
        changes.groupby(["previous_classification", "classification"])
        .size()
        .reset_index(name="proteins")
    )
    changes_summary.to_csv(OUT / "changes-from-1.39-summary.csv", index=False)
    sources = {}
    for row in cta_literature.attestations().to_dict("records"):
        sources.setdefault(
            row["source_id"], {k: row[k] for k in ["source_url", "citation", "title"]}
        )
    draft = pd.read_csv(DRAFT, usecols=["entity_id", "classification"], low_memory=False)
    data["draft_classification"] = data.entity_id.map(draft.set_index("entity_id").classification)
    moved = data[data.draft_classification.ne(data.classification)]
    moved[
        [
            "entity_id",
            "protein_group_label",
            "member_gene_ids",
            "draft_classification",
            "classification",
            "decision_reason",
            "published_nomination",
            "rare_cell_status",
            "rare_cell_basis",
            "rare_cell_mapping_notes",
            "literature_effect",
            "literature_decision_reason",
            "finding_tissue",
            "finding_top_cell_type",
            "finding_studies_multi_molecule",
            "finding_top_multi_molecule_pct",
            "finding_bulk_tpm",
            "max_detection_pct",
            "max_multi_molecule_pct",
            "sequence_partners",
            "same_strand_hosts",
            "cancer_evidence",
            "cancer_routes",
            "tumor_max_tpm",
            "single_tumor_threshold_tpm",
        ]
    ].sort_values(["draft_classification", "classification", "protein_group_label"]).to_csv(
        OUT / "changes-from-draft.csv", index=False
    )
    steps = funnel(data)
    steps.to_csv(OUT / "funnel.csv", index=False)
    columns = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "representative_gene_id",
        "category",
        "cta_tier",
        "classification",
        "released_category",
        "previous_classification",
        "disposable_max_tpm",
        "decision_reason",
        "rare_cell_status",
        "rare_cell_basis",
        "rare_cell_mapping_notes",
        "literature_effect",
        "literature_decision_reason",
        "finding_tissue",
        "finding_top_cell_type",
        "finding_studies_multi_molecule",
        "finding_top_multi_molecule_pct",
        "max_detection_pct",
        "max_multi_molecule_pct",
        "normal_reference_tpm",
        "normal_reference_contexts",
        "tumor_max_tpm",
        "tumor_top_sample",
        "tumor_top_cancer_code",
        "independent_patients_above_normal",
        "unresolved_patient_samples_above_normal",
        "cancer_evidence",
        "cancer_routes",
        "reproductive_routes",
        "reproductive_fold",
        "reproductive_max_sample_tpm",
        "reproductive_supported_tpm",
        "reproductive_supported_context",
        "reproductive_qualifying_compartments",
        "reproductive_expression_contexts",
        "reproductive_top_sample",
        "reproductive_top_context",
        "reproductive_samples_at_tenfold",
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "anonymous_samples_at_tenfold",
        "cancer_types_at_tenfold",
        "reproductive_ntpm",
        "reproductive_cell_ntpm",
        "developmental_mean_tpm",
        "published_nomination",
        "nomination_source_ids",
        "paper_source_ids",
        "representative_transcripts",
    ]
    rows = data[columns].astype(object).where(data[columns].notna(), None).to_dict("records")
    knowledge = gene_knowledge.findings(include_updates=True).set_index("id").to_dict("index")
    superseded = {i for r in knowledge.values() for i in r["supersedes"]}
    finding_ids = data.set_index("entity_id").knowledge_finding_ids
    for row in rows:
        row["loss_consequences"] = []
        for context in json.loads(row.pop("reproductive_expression_contexts") or "[]"):
            loss = tissue_tiers.consequences(
                context["tissue"],
                context["cell_type"],
                stage=context["stage"],
                resolution="bulk" if context["cell_type"] == "bulk" else "cell",
            )
            row["loss_consequences"].append(
                {
                    **context,
                    "tier": loss["tier"],
                    "matched_level": loss["matched_level"],
                    "persistent_loss": loss.get(
                        "persistent_loss", "No applicable tissue-loss assessment"
                    ),
                    "management": loss.get("management", ""),
                }
            )
        row["reviewed_findings"] = [
            {
                "id": i,
                "finding": knowledge[i]["finding"],
                "reason": knowledge[i]["decision_reason"],
                "limitations": knowledge[i]["limitations"],
                "source_url": knowledge[i]["source_url"],
                "superseded": i in superseded,
            }
            for i in json.loads(finding_ids[row["entity_id"]])
            if i in knowledge
        ]
        for key in [
            "member_gene_ids",
            "nomination_source_ids",
            "paper_source_ids",
            "normal_reference_contexts",
        ]:
            row[key] = json.loads(row[key]) if row[key] else []
    payload = json.dumps(rows, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
    source_payload = json.dumps(sources, allow_nan=False).replace("<", "\\u003c")
    counts = summary["categories"]
    count_rows = "".join(
        f"<tr><td>{html.escape(c)}</td><td>{counts.get(c, 0):,}</td><td>{summary['published_nomination_by_category'].get(c, 0):,}</td></tr>"
        for c in target_selection.CATEGORIES
    )
    migrations = "".join(
        f"<tr><td>{html.escape(r.previous_classification)}</td><td>{html.escape(r.classification)}</td><td>{r.proteins:,}</td></tr>"
        for r in changes_summary.itertuples()
    )
    options = "".join(f"<option>{html.escape(c)}</option>" for c in target_selection.CATEGORIES)
    downloads = " · ".join(
        f'<a href="{c.lower().replace(" ", "-")}.csv.gz">{html.escape(c)}</a>'
        for c in target_selection.CATEGORIES
    )
    members = summary["memberships"]
    reproductive_only = data[data.category.eq("Reproductive Antigen only")]
    unresolved_reproductive = int(reproductive_only.cancer_evidence.eq("unassessed").sum())
    other_counts = ", ".join(
        f"{summary['classifications'].get('Other:' + name, 0):,} {name.lower()}"
        for name in ["Review", "Sequence", "Expression", "Unscreened tissue"]
    )
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Rebuilt antigen catalog — OncoRef</title>
<style>body{font:16px/1.5 system-ui,sans-serif;color:#20333d;margin:0}main{max-width:1550px;margin:auto;padding:28px}
h1{font-size:34px}h2{margin-top:2em}a{color:#126487}.summary{max-width:850px}.note{padding:14px;background:#f4f6f7;border-left:4px solid #b17a31}
table{border-collapse:collapse;width:100%}th,td{padding:9px;text-align:left;border-bottom:1px solid #dbe1e4;vertical-align:top}th{background:#eef3f5;position:sticky;top:0}
.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin:20px 0}input,select,button{font:inherit;padding:7px}input[type=search]{min-width:330px}
.scroll{overflow:auto;max-height:78vh}.catalog{min-width:1250px;font-size:14px}small{color:#576771}details{max-width:480px}details ul{padding-left:18px}
.badge{font-weight:600;white-space:nowrap}.Unsafe{color:#a33532}.CTA{color:#20613c}#count{font-weight:bold}.pager{margin-top:12px;display:flex;gap:16px}
</style><main><h1>Rebuilt antigen catalog</h1>
<p>__TOTAL__ proteins assessed under <code>__POLICY__</code>. Identical selected reference proteins are counted together; every member locus remains listed.</p>
<p class="note"><strong>These counts are outputs of screening rules, not validated tumor antigens.</strong> Rare normal-cell RNA is graded per adult donor; established findings exclude, strong unresolved signals require review, and weak or inconsistent findings retain their reasons. <a href="literature/">Literature review of the dismissed and released proteins</a>. <a href="../../normal-cell-atlas-20261110/rare-cell/">Donor-level evidence</a> · <a href="changes-from-draft.csv">Every change from the 2026-11-09 draft</a> · <a href="../../target-selection-20261109/scrutiny/">Earlier 250-protein audit</a>.</p>
<div class="summary"><table><thead><tr><th>Exclusive category</th><th>Proteins</th><th>With a curated target nomination</th></tr></thead><tbody>__COUNTS__</tbody></table>
<p><strong>Cancer Antigens: __CA__ · Reproductive Antigens: __RA__ · CTA intersection: __CTA__ · Candidate antigen proteins union: __UNION__.</strong></p></div>
<figure><img src="summary/source-to-tiers-architecture.png" alt="End-to-end curation architecture from data sources to CTA tiers" style="width:100%;max-width:1350px"></figure>
<figure><img src="summary/cta-tiers-and-profiles.png" alt="CTA compartment tiers and host toxicity profiles" style="width:100%;max-width:1300px"></figure>
<figure><img src="summary/antigen-sets.png" alt="Cancer and reproductive antigen membership" style="width:100%;max-width:900px"></figure>
<figure><img src="summary/funnel-all.png" alt="Sequential curation funnel with removed and retained counts" style="width:100%;max-width:1300px"></figure>
<p><a href="summary/">All curation figures and count tables</a> · <a href="library-qc.csv.gz">Library QC measurements</a> · <a href="../../../antigen-sequences/">Set and peptide subtraction APIs</a></p>
<p>Other: __OTHER_COUNTS__. The __RA_ONLY__ Reproductive Antigens only have no qualifying cancer evidence under this screen.</p>
<p>The shared eligibility screen leaves __ELIGIBLE__ proteins. Cancer RNA needs one tumor sample, from any comparable TPM cohort, at ≥10 × max(N, D, 1), where N is the highest non-reproductive normal value and D the highest Disposable value. All qualifying 10× samples, known patients, anonymous samples and cancer types are counted separately. The highest tumor, reproductive or developmental value must exceed every measured nonreproductive and Disposable value, with Disposable expression at most 32 TPM and 32 nTPM. RA independently requires ≥10 × max(N, D, 1) in a reproductive tissue sample or supported embryo-specific lineage mean. Protein staining and incompatible RNA scales remain independent evidence layers; they cannot bypass the RNA fold.</p>
<p class="note">The comparison uses native gene TPM with resolved identities and source quality checks. N includes within-site GTEx v11 linear donor P95, HPA v23 within-tissue biopsy P95 and supported tissue/cell summaries; it does not cover every healthy donor. Cross-study batch effects are not corrected. Other assay scales retain their independent normal-exclusion rules. <a href="normal-rna/">Normal RNA source coverage</a> · <a href="comparison.json">Exact comparison scope</a> · <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5903355/">Cross-study methods</a>.</p>
<p><strong>NUTM coverage:</strong> the stored UNC three-sample matrix fails the nonnegative-TPM check and contributes no qualifying samples. Its failed source remains in the <a href="tumor-source-coverage.csv">coverage ledger</a>; <a href="https://github.com/pirl-unc/oncoref/issues/604">#604</a> and <a href="https://github.com/pirl-unc/oncoref/issues/665">#665</a> track the repair. No raw reads were re-quantified for this rebuild.</p>
<h2>Changes from the released 1.39.0 catalog</h2><div class="summary"><table><thead><tr><th>Released category</th><th>New category</th><th>Proteins</th></tr></thead><tbody>__MIGRATIONS__</tbody></table></div>
<p><a href="changes-from-1.39.csv">All changes since 1.39 and their evidence</a> · <a href="summary.json">Counts</a> · <a href="gene-decisions.csv.gz">All member genes</a> · <a href="protein-decisions.csv.gz">Full protein decisions</a></p>
<p>Category downloads: __DOWNLOADS__.</p>
<h2>Find a protein</h2><p>“Only” describes qualified set membership. Missing evidence is shown as —, never zero. Measured cancer and reproductive counterfactuals remain visible for newly held rare-cell records. A review hold still blocks membership.</p>
<div class="controls"><input id="search" type="search" placeholder="Protein, Ensembl ID, reason, paper or program" aria-label="Search proteins">
<select id="category" aria-label="Category"><option value="">All categories</option>__OPTIONS__</select>
<label><input type="checkbox" id="changed"> Changed category</label><label><input type="checkbox" id="published"> Published target nomination</label>
<span id="count"></span></div>
<div class="scroll"><table class="catalog"><thead><tr><th>Protein and source links</th><th>New category</th><th>1.39 category</th><th>Normal N / Disposable D<br>TPM</th><th>Tumor maximum<br>TPM</th><th>Tumor 10× samples /<br>known patients</th><th>Reproductive fold</th><th>Reproductive maximum<br>TPM</th><th>Reproductive context</th><th>Evidence and papers</th></tr></thead><tbody id="body"></tbody></table></div>
<div class="pager"><button id="previous">Previous</button><span id="page"></span><button id="next">Next</button></div>
<p><a href="normal-comparison.csv.gz">All comparison normal measurements</a> · <a href="tumor-witnesses.csv.gz">Tumor witnesses</a> · <a href="paper-attestations.json.gz">Every retained protein–paper claim</a> · <a href="inputs.json">Input hashes</a> · <a href="../report.html">Ontology figures and rules</a></p>
<script id="records" type="application/json">__DATA__</script><script id="sources" type="application/json">__SOURCES__</script>
<script>
const records=JSON.parse(document.getElementById('records').textContent),sources=JSON.parse(document.getElementById('sources').textContent);
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=x=>x===null?'—':Number(x).toLocaleString(undefined,{maximumFractionDigits:3});
const link=(url,label)=>`<a href="${esc(url)}" target="_blank" rel="noopener">${esc(label)}</a>`;
let page=0,selected=[];const size=100;
function papers(ids){return ids.map(id=>{let s=sources[id];return s?link(s.source_url,s.citation)+(s.title?' — '+esc(s.title):''):esc(id)}).join('<br>')}
function render(reset=true){if(reset)page=0;let q=document.getElementById('search').value.toLowerCase(),cat=document.getElementById('category').value;
selected=records.filter(r=>(!cat||r.category===cat)&&(!document.getElementById('changed').checked||r.classification!==r.previous_classification)&&(!document.getElementById('published').checked||r.published_nomination)&&(!q||JSON.stringify(r).toLowerCase().includes(q)||r.paper_source_ids.some(id=>JSON.stringify(sources[id]||{}).toLowerCase().includes(q))));
const pages=Math.max(1,Math.ceil(selected.length/size));page=Math.min(page,pages-1);
document.getElementById('count').textContent=selected.length.toLocaleString()+' proteins';
document.getElementById('body').innerHTML=selected.slice(page*size,(page+1)*size).map(r=>`<tr><td><strong>${esc(r.protein_group_label)}</strong><br><small>${r.member_gene_ids.map(g=>link('https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+g,g)).join('<br>')}<br>${link('https://www.proteinatlas.org/'+r.representative_gene_id,'HPA')}</small></td>
<td class="badge ${r.category==='Unsafe'?'Unsafe':r.category==='CTA'?'CTA':''}">${esc(r.classification)}</td><td>${esc(r.previous_classification)}</td><td>${fmt(r.normal_reference_tpm)} / ${fmt(r.disposable_max_tpm)}</td><td>${fmt(r.tumor_max_tpm)}<br><small>${esc(r.tumor_top_cancer_code)}</small></td><td>${fmt(r.samples_at_tenfold)} / ${fmt(r.independent_patients_at_tenfold)}</td><td>${fmt(r.reproductive_fold)}</td><td>${fmt(r.reproductive_supported_tpm)}</td><td>${esc(r.reproductive_supported_context)}</td>
<td><details><summary>${r.published_nomination?'Published nomination':'No catalog nomination'} · ${r.paper_source_ids.length} paper sources</summary><p>${esc(r.decision_reason)}</p><p>Cancer assessment: ${esc(r.cancer_evidence)}<br>${esc(r.cancer_routes)}<br>Tenfold witnesses: ${fmt(r.samples_at_tenfold)} samples, ${fmt(r.independent_patients_at_tenfold)} known patients, ${fmt(r.anonymous_samples_at_tenfold)} anonymous samples, ${fmt(r.cancer_types_at_tenfold)} cancer types</p><p>Rare normal-cell RNA: ${esc(r.rare_cell_status)}${r.rare_cell_basis?' — '+esc(r.rare_cell_basis):''}${r.rare_cell_mapping_notes?' ('+esc(r.rare_cell_mapping_notes)+')':''}${r.literature_effect&&r.literature_effect!=='not_reviewed'?'<br>Literature: '+esc(r.literature_effect)+' — '+esc(r.literature_decision_reason):''}${r.finding_tissue?'<br>'+esc(r.finding_tissue)+': '+esc(r.finding_top_cell_type)+', '+fmt(r.finding_studies_multi_molecule)+' studies, '+fmt(r.finding_top_multi_molecule_pct)+'% of the strongest donor’s cells with ≥2 molecules':''}<br>Highest single-donor detection ${fmt(r.max_detection_pct)}%; with ≥2 molecules ${fmt(r.max_multi_molecule_pct)}%</p><p>Reproductive/developmental routes: ${esc(r.reproductive_routes)}${r.reproductive_qualifying_compartments?'<br>Qualifying compartments: '+esc(r.reproductive_qualifying_compartments):''}<br>Highest cell-supported observation: ${fmt(r.reproductive_max_sample_tpm)} TPM; replication-supported ceiling: ${fmt(r.reproductive_supported_tpm)} TPM / ${fmt(r.reproductive_fold)}× in ${esc(r.reproductive_supported_context)}; ${esc(r.reproductive_top_sample)}<br>Context: ${esc(r.reproductive_top_context)}; ${fmt(r.reproductive_samples_at_tenfold)} observations at 10×</p><p><strong>Consequences of tissue loss</strong> · <a href="../../../tissue-tiers/">Reviewed tissue policy</a><br>Observed reproductive RNA ≥1 TPM; not restricted to tenfold RA witnesses. These observed maxima can include unsupported cell means; QC indicates whether any observation in that context passes the cell-support screen. Organ fallbacks remain explicit.</p>${r.loss_consequences.map(c=>`<p>${esc(c.compartment&&c.compartment!==c.tissue?c.tissue+' ('+c.compartment+')':c.tissue)} / ${esc(c.cell_type)}: ${fmt(c.native_tpm_max)} TPM${c.qc_supported?"":" (no observation passes cell-support screen)"}<br>${esc(c.tier)} (${esc(c.matched_level||"unresolved scope")}) — ${esc(c.persistent_loss)}<br>${esc(c.management)}</p>`).join('')}<p><strong>Reviewed gene findings</strong></p>${r.reviewed_findings.map(f=>`<p><a href="${esc(f.source_url)}">${esc(f.id)}</a>${f.superseded?" (superseded)":""}<br>${esc(f.finding)}<br>${esc(f.reason)}<br><em>${esc(f.limitations)}</em></p>`).join('')}<p>Peak tumor sample: ${esc(r.tumor_top_sample)}</p><p>Selected transcript(s): ${esc(r.representative_transcripts)}</p><p>Normal maximum witnesses (up to 10):</p><ul>${r.normal_reference_contexts.map(x=>`<li>${esc(x.tissue)} / ${esc(x.cell_type)}: ${fmt(x.value)} TPM · ${esc(x.source_id)} · ${esc(x.statistic)} · ${esc(x.sample_id||x.source_context)}</li>`).join('')}</ul><p><strong>Target nominations</strong><br>${papers(r.nomination_source_ids)||'None recorded'}</p><p><strong>All paper roles, including counterevidence</strong><br>${papers(r.paper_source_ids)||'None recorded'}</p></details></td></tr>`).join('');
document.getElementById('page').textContent=`Page ${page+1} of ${pages}`;document.getElementById('previous').disabled=page===0;document.getElementById('next').disabled=page+1>=pages;}
for(const id of ['search','category','changed','published'])document.getElementById(id).addEventListener(id==='search'?'input':'change',()=>render());
document.getElementById('previous').onclick=()=>{page--;render(false)};document.getElementById('next').onclick=()=>{page++;render(false)};render();
</script></main></html>"""
    replace = {
        "TOTAL": f"{len(data):,}",
        "POLICY": target_selection.POLICY_ID,
        "COUNTS": count_rows,
        "CA": f"{members['is_cancer_antigen']:,}",
        "RA": f"{members['is_reproductive_antigen']:,}",
        "CTA": f"{members['is_cta']:,}",
        "UNION": f"{members['is_target_candidate']:,}",
        "ELIGIBLE": str(summary["shared_eligible"]),
        "OTHER_COUNTS": other_counts,
        "RA_ONLY": str(len(reproductive_only)),
        "RA_UNRESOLVED": str(unresolved_reproductive),
        "RA_INSUFFICIENT": str(len(reproductive_only) - unresolved_reproductive),
        "MIGRATIONS": migrations,
        "DOWNLOADS": downloads,
        "OPTIONS": options,
        "DATA": payload,
        "SOURCES": source_payload,
    }
    for key, value in replace.items():
        page = page.replace("__" + key + "__", value)
    (OUT / "report.html").write_text(page)
    funnel_figure(steps, data, counts)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.2, 1]})
    labels = ["CTA", "Cancer Antigen only", "Reproductive Antigen only"]
    axes[0].barh(
        labels[::-1], [counts[x] for x in labels[::-1]], color=["#c99a46", "#3e8bba", "#438a69"]
    )
    for i, label in enumerate(labels[::-1]):
        axes[0].text(counts[label] + 5, i, f"{counts[label]:,}", va="center")
    axes[0].set_xlim(0, max(counts[x] for x in labels) * 1.17)
    axes[0].set_title(f"{members['is_target_candidate']:,} Candidate antigen proteins", loc="left")
    axes[1].barh(
        ["Other", "Unsafe"], [counts["Other"], counts["Unsafe"]], color=["#9aabb5", "#ac5552"]
    )
    for i, label in enumerate(["Other", "Unsafe"]):
        axes[1].text(counts[label] + 130, i, f"{counts[label]:,}", va="center")
    axes[1].set_xlim(0, counts["Unsafe"] * 1.22)
    axes[1].set_title("Remaining inventory", loc="left")
    for ax in axes:
        ax.set_xlabel("Proteins")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "counts.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    new_cta = data[data.category.eq("CTA") & data.released_category.ne("CTA")]
    held_cta = data[data.category.eq("Other") & data.released_category.eq("CTA")]
    migration_table = "\n".join(
        f"| {r['released_category']} | {r['category']} | {r['proteins']} |"
        for r in summary["changes"]
    )
    table_rows = "\n".join(
        f"| {r.protein_group_label} | {r.tumor_max_tpm:.3g} | {r.normal_reference_tpm:.3g} | "
        f"{'Patient identity unresolved' if r.cancer_evidence == 'unassessed' else 'Below cancer rule'} |"
        for r in reproductive_only.itertuples()
    )
    funnel_rows = "\n".join(
        f"| {r.step} | {r.result} | {r.removed:,} | {r.remaining:,} |" for r in steps.itertuples()
    )
    transitions = (
        moved.groupby(["draft_classification", "classification"]).size().reset_index(name="n")
    )
    transition_rows = "\n".join(
        f"| {r.draft_classification} | {r.classification} | {r.n} |"
        for r in transitions.itertuples()
    )
    released = data[data.draft_classification.eq("Unsafe") & data.is_target_candidate]
    new_unsafe = data[data.draft_classification.ne("Unsafe") & data.category.eq("Unsafe")]
    dismissed = data[data.is_target_candidate & data.rare_cell_status.eq("dismissed")]
    document = json.loads((ROOT / "oncoref/data/normal-expression-literature.json").read_text())
    excluded = [e for e in document["entries"] if e["effect"] == "exclude"]
    rules = target_selection.definition()["rare_cell_rna"]
    (OUT / "index.md").write_text(
        "# Rebuilt antigen catalog\n\n"
        f"Policy `{target_selection.POLICY_ID}`. **{counts['CTA']} CTA, "
        f"{counts['Cancer Antigen only']} Cancer Antigen only, "
        f"{counts['Reproductive Antigen only']} Reproductive Antigen only, {counts['Other']} Other "
        f"and {counts['Unsafe']:,} Unsafe.** These are outputs of screening rules, not validated "
        "tumor antigens.\n\n"
        "[Filterable protein table](report.html) · [Changes from 1.39.0](changes-from-1.39.csv) · [Transition counts](changes-from-1.39-summary.csv) · "
        "[Changes from release 1.36.0](category-changes.csv) · [Exact comparison scope](comparison.json)\n\n"
        "## End-to-end architecture and CTA tiers\n\n"
        "![End-to-end curation architecture](summary/source-to-tiers-architecture.png)\n\n"
        "![CTA tiers and host profiles](summary/cta-tiers-and-profiles.png)\n\n"
        "## Curation funnel\n\n![Curation funnel](funnel.png)\n\n"
        "| Filter | Result | Proteins removed | Remaining |\n|---|---|---:|---:|\n"
        + funnel_rows
        + f"\n\nThe {members['is_target_candidate']} remaining Candidate antigen proteins are "
        f"{counts['CTA']} CTA, {counts['Cancer Antigen only']} Cancer Antigen only and "
        f"{counts['Reproductive Antigen only']} Reproductive Antigen only. [Funnel table](funnel.csv) · [Curation summary with figures](summary.md).\n\n"
        "![Rebuilt counts](counts.png)\n\n"
        "## The filters\n\n"
        "| Filter | Rule | Failure |\n|---|---|---|\n"
        "| Normal expression | Existing bulk, cluster, immune, cornea, protein and reviewed-knowledge rules, unchanged | Unsafe for accepted exclusions; otherwise Other:Review |\n"
        f"| Rare normal-cell RNA | {rules['signal']} {rules['replication']} | Established: Unsafe. Unresolved: Other:Review. Dismissed: no effect |\n"
        "| Published normal expression | Gene-specific report in normal adult Vital or Protected tissue, by a fixed rubric | Unsafe |\n"
        "| Coverage, identity and translation | Unchanged | Other:Review |\n"
        "| Sequence | Fewer than 95% of 9-mers inside Unsafe proteins; no unsupported nested protein or read-through locus | Other:Sequence |\n"
        "| Expression bounds | Target maximum above every measured nonreproductive and Disposable value; Disposable at most 32 TPM and 32 nTPM | Other:Expression |\n"
        "| Cancer evidence | One tumor sample ≥10 × max(N, Disposable, 1), N the highest non-reproductive normal value | No Cancer Antigen membership |\n"
        "| Reproductive evidence | One reproductive sample or supported embryo-lineage mean ≥10 × max(N, Disposable, 1), on the same scale and reference as CA | No Reproductive Antigen membership |\n\n"
        "## Rare normal-cell RNA\n\n"
        "A count matrix cannot show where a read aligned, and one mis-assigned read per cell from a "
        "paralog, pseudogene or overlapping gene looks like sparse expression. "
        f"{rules['single_molecules']} Strong unresolved signals are held for review; mapping flags do not prove artifact.\n\n"
        f"- **Established, Unsafe:** {rules['established']}\n"
        f"- **Unresolved, Other:Review:** {rules['unresolved']}\n"
        f"- **Dismissed, no effect:** {rules['dismissed']}\n"
        f"- **Mapping risk:** {rules['mapping_risk']}\n"
        f"- **Bulk confirmation:** {rules['confirmation']}\n"
        f"- **Literature:** {rules['literature']}\n\n"
        f"{len(new_unsafe)} draft non-Unsafe proteins become Unsafe, among them "
        + ", ".join(
            new_unsafe[new_unsafe.draft_classification.isin(["CTA", "Cancer Antigen only"])]
            .sort_values("finding_top_multi_molecule_pct", ascending=False)
            .protein_group_label.head(14)
        )
        + f". {len(dismissed)} candidates carry a dismissed signal, among them "
        + ", ".join(
            dismissed.sort_values(
                "finding_top_multi_molecule_pct", ascending=False
            ).protein_group_label.head(10)
        )
        + ".\n\n"
        f"**{len(released)} proteins that were Unsafe only through the pooled 1% detection flag are "
        "candidates again**: "
        + ", ".join(sorted(released.protein_group_label))
        + ". Each was taken through the literature review below and none has a gene-specific report "
        "of expression in vital normal tissue. That is a search result, not proof of absence. "
        "[Donor-level evidence](../../normal-cell-atlas-20261110/rare-cell/index.md).\n\n"
        "## Literature review\n\n"
        f"{len(document['entries'])} proteins were reviewed: {document['scope'][0].lower()}{document['scope'][1:]} "
        f"{len(excluded)} are Unsafe on published evidence: "
        + "; ".join(f"{e['protein']} ({e['decision_tissue']}, {e['rule']})" for e in excluded)
        + ". [Every protein, rubric, sources and quotes](literature.md).\n\n"
        "## Changes from the 2026-11-09 draft\n\n| Draft | Now | Proteins |\n|---|---|---:|\n"
        + transition_rows
        + "\n\nThe tenfold single-tumor rule accounts for the moves to Reproductive Antigen only and for "
        "Cancer Antigens only whose reproductive compartment does not reach the shared threshold. "
        "Known patient identity is a replication descriptor, not a single-sample admission gate.\n\n"
        f"Cancer Antigens include {members['is_cancer_antigen']} proteins; Reproductive Antigens include "
        f"{members['is_reproductive_antigen']}. Their intersection is CTA and their union contains "
        f"{members['is_target_candidate']} Candidate antigen proteins.\n\n"
        f"Other contains {other_counts}. A curated published target nomination is retained for "
        f"{summary['published_nomination_by_category']['CTA']} CTAs and "
        f"{summary['published_nomination_by_category']['Cancer Antigen only']} Cancer Antigens only. "
        "Broader paper coverage also includes biology, normal-expression findings and other evidence.\n\n"
        "## Changes from release 1.36.0\n\n| Released category | Rebuilt category | Proteins |\n|---|---|---:|\n"
        + migration_table
        + "\n\nNew CTAs: "
        + ", ".join(new_cta.protein_group_label)
        + ".\n\nFormer CTAs now requiring review: "
        + ", ".join(held_cta.protein_group_label)
        + ".\n\n## Reproductive Antigens only\n\n"
        + "No tumor sample reaches 10 × the protein's non-reproductive normal maximum"
        + (
            f"; {unresolved_reproductive} lack a complete normal comparison.\n\n"
            if unresolved_reproductive
            else ".\n\n"
        )
        + "| Protein | Tumor maximum TPM | Normal comparison TPM | Cancer assessment |\n"
        "|---|---:|---:|---|\n" + table_rows + "\n\n## Evidence and reproduction\n\n"
        "The comparison uses source-native TPM, within-site GTEx v11 linear donor P95, HPA v23 within-tissue biopsy P95 and supported tissue/cell summaries. It is not a maximum over every healthy donor. "
        "It does not correct cross-study batch effects. NUTM's stored three-sample matrix fails "
        "the nonnegative-TPM check and is listed in the [coverage ledger](tumor-source-coverage.csv).\n\n"
        "[Full protein decisions](protein-decisions.csv.gz) · [All member genes](gene-decisions.csv.gz) · "
        "[HPA sample QC](normal-library-qc.csv) · [GTEx sample QC and donors](gtex-library-qc.csv.gz) · [Normal RNA source coverage](normal-rna/index.md) · [Normal comparison values](normal-comparison.csv.gz) · [Tumor witnesses](tumor-witnesses.csv.gz) · [Reproductive observations](reproductive-comparison.csv.gz) · [Reproductive witnesses](reproductive-witnesses.csv.gz) · "
        "[Paper claims](paper-attestations.json.gz) · [Input hashes](inputs.json).\n\n"
        "The [pinned HPA download receipt](downloaded-hpa-sources.json) supplies exact URLs, "
        "paths and SHA-256 hashes for the newly retrieved native-TPM tables. With the pinned "
        "source files in the paths listed in `inputs.json`, run:\n\n"
        "```sh\npython scripts/build_rare_cell_evidence.py --regrade-pinned\n"
        "python scripts/rebuild_target_catalog.py\npython scripts/render_target_catalog.py\n"
        "python scripts/build_target_selection_report.py\n```\n\n"
        "The builder verifies the pinned tumor matrices, re-resolves source identities and assesses "
        "all proteins. Temporary tumor matrices stay outside the published report. "
        "The new catalog is available through `target_selection.results()`; the `cta_selection` "
        "API and existing sequence downloads retain the 1.36.0 snapshot.\n"
    )

    def link(e):
        return "; ".join(
            f"[{(x['source_title'] or 'source')[:60]}]({x['url']})"
            for x in e["evidence"][:3]
            if x["url"]
        )

    order = {"exclude": 0, "confirm_detection": 0, "none": 1}
    table = "\n".join(
        f"| {e['protein']} | {e['verdict'].replace('_', ' ')} | "
        f"{'Unsafe (' + e['rule'] + ', ' + e['decision_tissue'] + ')' if e['effect'] != 'none' else 'No effect'} | "
        f"{e['decision_reason']} | {link(e)} |"
        for e in sorted(document["entries"], key=lambda e: (order[e["effect"]], e["protein"]))
    )
    (OUT / "literature.md").write_text(
        "# Literature review of normal expression\n\n"
        f"**{len(document['entries'])} proteins reviewed on {document['review_date']}; {len(excluded)} are "
        "Unsafe on published evidence.** Scope: "
        + document["scope"]
        + "\n\n## Rubric\n\n| Rule | Requirement |\n|---|---|\n"
        + "\n".join(f"| {k} | {v} |" for k, v in document["rubric"].items())
        + "\n\n## Limits\n\n"
        + document["method"]
        + "\n\n"
        + "\n".join(f"- {x}." for x in document["limitations"])
        + "\n\n## Decisions\n\n| Protein | Search verdict | Effect | Reason | Sources |\n|---|---|---|---|---|\n"
        + table
        + "\n\nSummaries, every evidence item and its quote are in "
        "`oncoref/data/normal-expression-literature.json`.\n"
    )
    render_normal_rna_sources.build()
    render_curation_summary.build()
    implementation = [
        "scripts/rebuild_target_catalog.py",
        "scripts/render_target_catalog.py",
        "scripts/render_curation_summary.py",
        "scripts/render_normal_rna_sources.py",
        "scripts/audit_normal_library_qc.py",
        "scripts/gtex_sample_evidence.py",
        "scripts/bulk_normal_percentiles.py",
        "oncoref/normal_library_qc.py",
        "scripts/build_target_selection_report.py",
        "scripts/build_rare_cell_evidence.py",
        "scripts/build_gene_mapping_risk.py",
        "scripts/build_cta_candidate_evidence.py",
        "scripts/build_cta_alternative_evidence.py",
        "scripts/antigen_sample_evidence.py",
        "scripts/cta_report_common.py",
        "oncoref/target_selection.py",
        "oncoref/gene_matching.py",
        "oncoref/gene_knowledge.py",
        "oncoref/library_qc.py",
        "scripts/build_tissue_marker_panels.py",
        "oncoref/cta_literature.py",
        "oncoref/target_filters.py",
        "oncoref/host_context.py",
        "scripts/build_host_exposure.py",
        "scripts/cta_evidence_assets/target-selection-guide.md",
        "scripts/cta_evidence_assets/target-filter-guide.md",
        "oncoref/tissue_tiers.py",
        "oncoref/samples.py",
    ]
    manifest = {
        "policy_id": target_selection.POLICY_ID,
        "implementation_sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in implementation
        },
        "outputs_sha256": {
            str(path.relative_to(OUT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(OUT.rglob("*"))
            if path.is_file() and path.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    build()
