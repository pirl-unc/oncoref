#!/usr/bin/env python3
"""Publish the current reference-selection method, gene changes and review queue."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cta_report_common import input_hashes, sha256, write_csv
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from migrate_transcript_selection import BASELINE_COMMIT, baseline

from oncoref import cta_literature, transcript_selection, translation_qualification

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/cta-transcript-selection-20261108"
CURRENT = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"


def by_gene(facts):
    return {gene: r for r in facts.to_dict("records") for gene in json.loads(r["member_gene_ids"])}


def figure():
    matplotlib.rcParams["svg.hashsalt"] = "oncoref-transcript-selection-20261108"
    fig, ax = plt.subplots(figsize=(11, 5.4))
    fig.patch.set_facecolor("#faf9f5")
    ax.set(xlim=(0, 11), ylim=(0, 5.4))
    ax.axis("off")
    boxes = [
        (0.25, 3.6, "Choose a reference", "MANE Select\notherwise Ensembl canonical"),
        (
            3.9,
            3.6,
            "Check the exact CDS",
            "Primary locus · complete codons\nqualifying annotation · AA20",
        ),
        (
            7.55,
            3.6,
            "Review conflicting evidence",
            "Exact reference / sequence\nlinked sources · explicit decision",
        ),
        (
            0.25,
            0.7,
            "Keep unresolved candidates",
            "No silent fallback by length\nretain alternatives and reasons",
        ),
        (
            3.9,
            0.7,
            "Approve a representative",
            "One qualified sequence per locus\nidentical defaults share an identity",
        ),
        (
            7.55,
            0.7,
            "Assess biological evidence",
            "Cancer · reproductive / fetal\nnormal cells and coverage\nsequence containment",
        ),
    ]
    for x, y, title, body in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                3.15,
                1.25,
                boxstyle="round,pad=0.10",
                facecolor="#e4eee9",
                edgecolor="#467365",
            )
        )
        ax.text(x + 0.16, y + 0.92, title, fontsize=10.4, weight="bold", color="#173e35")
        ax.text(x + 0.16, y + 0.45, body, fontsize=9, linespacing=1.6, va="center")
    for start, end in [
        ((3.45, 4.2), (3.75, 4.2)),
        ((7.1, 4.2), (7.4, 4.2)),
        ((9.15, 3.45), (5.65, 2.15)),
        ((5.55, 3.45), (1.9, 2.15)),
        ((7.1, 1.35), (7.4, 1.35)),
    ]:
        ax.add_patch(
            FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color="#467365")
        )
    ax.text(
        6.8,
        2.6,
        "Complete and resolved",
        fontsize=9,
        color="#173e35",
        bbox={"facecolor": "#faf9f5", "edgecolor": "none", "pad": 3},
    )
    ax.text(
        1.45,
        2.6,
        "Incomplete or unresolved",
        fontsize=9,
        color="#173e35",
        bbox={"facecolor": "#faf9f5", "edgecolor": "none", "pad": 3},
    )
    ax.text(
        0.25,
        0.12,
        "Reference selection does not establish protein production, tumor proteoform abundance or safety.",
        fontsize=10,
    )
    fig.savefig(OUT / "method.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "method.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def browser(rows):
    columns = [
        "protein",
        "gene_id",
        "classification",
        "category",
        "selection_status",
        "selection_basis",
        "candidate_transcripts",
        "selection_reasons",
        "support_flags",
        "candidate_support",
        "candidate_lengths",
        "evidence_ids",
        "default_changed",
    ]
    display = [{key: row[key] for key in columns} for row in rows]
    payload = json.dumps(display, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    html = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>OncoRef transcript choices</title>
<style>body{font:16px/1.5 system-ui;max-width:1450px;margin:2rem auto;padding:0 1rem;color:#163e35;background:#faf9f5}h1{font-size:2rem}input,select,button{font:inherit;padding:.5rem;margin:.25rem;border:1px solid #81998f;border-radius:4px}table{width:100%;border-collapse:collapse;font-size:13px}td,th{padding:.6rem;text-align:left;border-bottom:1px solid #d6dfd9;vertical-align:top}th{background:#e4eee9;position:sticky;top:0}a{color:#16605e}td:nth-child(5){max-width:420px}code{overflow-wrap:anywhere}</style>
<h1>Which transcript represents each gene?</h1><p>MANE Select first, Ensembl canonical second. A complete CDS and resolved coding evidence are required for approval. Cancer and normal-tissue measurements are assessed separately.</p>
<p><a href="./">Method and evidence</a> · <a href="gene-decisions.csv.gz">All decisions (CSV)</a> · <a href="gene-default-changes.csv.gz">Changed defaults (CSV)</a></p>
<label>Search <input id="search" placeholder="Gene, transcript, reason, study ID" size="38"></label><label>Status <select id="status"><option value="">All</option><option value="selected">Approved</option><option value="review">Needs review</option></select></label><label>Category <select id="category"><option value="">All</option><option>CTA</option><option>Unsafe</option><option>Other</option></select></label><label><input type="checkbox" id="changed">Changed default only</label><p id="count"></p>
<button id="prev">Previous</button><button id="next">Next</button><table><thead><tr><th>Protein / gene</th><th>Category</th><th>Decision</th><th>Candidate transcripts</th><th>Reasons and support</th><th>Lengths (aa)</th><th>Reviewed evidence</th></tr></thead><tbody id="rows"></tbody></table>
<script>const data=PAYLOAD;let page=0;const el=id=>document.getElementById(id);const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));function draw(){const query=el('search').value.toLowerCase();const rows=data.filter(r=>(!query||JSON.stringify(r).toLowerCase().includes(query))&&(!el('status').value||r.selection_status===el('status').value)&&(!el('category').value||r.category===el('category').value)&&(!el('changed').checked||r.default_changed));const pages=Math.max(1,Math.ceil(rows.length/100));page=Math.min(page,pages-1);el('count').textContent=`${rows.length.toLocaleString()} genes · page ${page+1} of ${pages} · blank means unavailable`;el('prev').disabled=page===0;el('next').disabled=page===pages-1;el('rows').innerHTML=rows.slice(page*100,(page+1)*100).map(r=>`<tr><td>${esc(r.protein)}<br><a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${encodeURIComponent(r.gene_id)}">${esc(r.gene_id)}</a></td><td>${esc(r.classification)}</td><td>${r.selection_status==='selected'?'Approved':'Needs review'}<br>${esc(r.selection_basis)}</td><td>${r.candidate_transcripts.map(t=>`<a href="https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=${encodeURIComponent(t)}">${esc(t)}</a>`).join('<br>')}</td><td>${esc([...r.selection_reasons,...r.support_flags].join('; ').replaceAll('_',' '))}<br>${r.candidate_support.map(s=>esc(s.transcript+': TSL '+(s.transcript_support_level||'missing'))).join('<br>')}</td><td>${esc(r.candidate_lengths.join(', '))}</td><td>${r.evidence_ids.map(id=>`<a href="./#findings-applied-in-this-release">${esc(id.replaceAll('_',' '))}</a>`).join('<br>')}</td></tr>`).join('')}for(const id of ['search','status','category','changed'])el(id).addEventListener('input',()=>{page=0;draw()});el('prev').onclick=()=>{page--;draw()};el('next').onclick=()=>{page++;draw()};draw();</script></html>"""
    (OUT / "table.html").write_text(html.replace("PAYLOAD", payload))


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    facts = pd.read_csv(CURRENT, keep_default_na=False, low_memory=False)
    oldpath = baseline("docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz")
    old = pd.read_csv(oldpath, keep_default_na=False, low_memory=False)
    now_genes, old_genes = by_gene(facts), by_gene(old)
    selected = transcript_selection.records().set_index("gene_id")
    rows = []
    for gene, current in sorted(now_genes.items()):
        previous = old_genes[gene]
        if gene in selected.index:
            r = selected.loc[gene].to_dict()
        else:
            r = {
                "selection_status": "review",
                "selection_basis": "No reference translation",
                "selection_reasons": ["missing_primary_translation"],
                "candidate_sequence_sha256s": [],
                "selected_sequence_sha256": "",
                "candidate_transcripts": [],
                "candidate_lengths": [],
                "candidate_support": [],
                "support_flags": [],
                "evidence_ids": [],
            }
        before = json.loads(previous["representative_sequence_sha256s"])
        after = [r["selected_sequence_sha256"]] if r["selected_sequence_sha256"] else []
        r.update(
            gene_id=gene,
            protein=current["protein_group_label"],
            category=current["category"],
            classification=current["classification"],
            previous_classification=previous["classification"],
            previous_protein=previous["protein_group_label"],
            previous_sequence_sha256s=before,
            default_changed=sorted(before) != sorted(after),
            protein_identity_changed=previous["protein_equivalence_group_id"]
            != current["protein_equivalence_group_id"],
        )
        rows.append(r)
    decisions = pd.DataFrame(rows)
    encoded = decisions.copy()
    for column in encoded:
        if encoded[column].map(lambda v: isinstance(v, (list, dict))).any():
            encoded[column] = encoded[column].map(json.dumps)
    write_csv(encoded, OUT / "gene-decisions.csv.gz")
    write_csv(encoded.loc[decisions.default_changed], OUT / "gene-default-changes.csv.gz")
    write_csv(
        encoded.loc[decisions.classification.ne(decisions.previous_classification)],
        OUT / "gene-classification-changes.csv.gz",
    )
    write_csv(encoded.loc[decisions.selection_status.eq("review")], OUT / "review-queue.csv.gz")
    surface = []
    for r in cta_literature.surface_domains().to_dict("records"):
        choices = selected.reindex(r["gene_ids"])
        approved = choices.selection_status.eq("selected").all()
        hashes = sorted(set(choices.selected_sequence_sha256)) if approved else []
        exact = hashes == [r["uniprot_sequence_sha256"]]
        categories = sorted({now_genes[g]["classification"] for g in r["gene_ids"]})
        surface.append(
            {
                "protein": r["protein"],
                "gene_ids": json.dumps(r["gene_ids"]),
                "classification": " | ".join(categories),
                "uniprot_url": r["uniprot_url"],
                "uniprot_sequence_sha256": r["uniprot_sequence_sha256"],
                "approved_sequence_sha256s": json.dumps(hashes),
                "exact_sequence_match": exact,
                "coordinate_mapping": "exact_sequence"
                if exact
                else "no_approved_default"
                if not approved
                else "sequence_mismatch",
                "selected_extracellular_total_aa": r["extracellular_total_aa"] if exact else None,
                "selected_extracellular_intervals": json.dumps(
                    r["extracellular_intervals"] if exact else []
                ),
                "feature_evidence_codes": json.dumps(r["extracellular_evidence_codes"]),
            }
        )
    write_csv(pd.DataFrame(surface), OUT / "surface-coordinate-review.csv")
    browser(rows)
    figure()
    counts = facts.category.value_counts().to_dict()
    old_cta = set(decisions.loc[decisions.previous_classification.eq("CTA"), "gene_id"])
    new_cta = set(decisions.loc[decisions.classification.eq("CTA"), "gene_id"])
    summary = {
        "policy_id": transcript_selection.POLICY_ID,
        "baseline_commit": BASELINE_COMMIT,
        "previous_categories": old.category.value_counts().to_dict(),
        "categories": counts,
        "protein_records": len(facts),
        "classified_loci": len(rows),
        "loci_with_changed_default": int(decisions.default_changed.sum()),
        "previous_cta_loci": len(old_cta),
        "current_cta_loci": len(new_cta),
        "cta_loci_added": sorted(new_cta - old_cta),
        "cta_loci_removed": sorted(old_cta - new_cta),
        "representative_review_loci": int(decisions.selection_status.eq("review").sum()),
        "surface_exact_matches": sum(r["exact_sequence_match"] for r in surface),
    }
    summary["previous_sequence_bundle_id"] = json.loads(
        baseline("oncoref/data/sequence-bundle.json").read_text()
    )["bundle_id"]
    summary["sequence_bundle_id"] = json.loads(
        (ROOT / "oncoref/data/sequence-bundle.json").read_text()
    )["bundle_id"]
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    claims = transcript_selection.evidence()["findings"]
    case_rows = []
    for claim in claims:
        gene = claim["gene_id"]
        sources = "; ".join(f"[{s['title']}]({s['url']})" for s in claim["sources"])
        case_rows.append(
            f"| {now_genes[gene]['protein_group_label']} | {claim['rationale']} | {now_genes[gene]['classification']} | {sources} |"
        )
    added = decisions.loc[decisions.gene_id.isin(new_cta - old_cta), "protein"].unique()
    removed = decisions.loc[decisions.gene_id.isin(old_cta - new_cta), "protein"].unique()
    new_rows = []
    for r in decisions.loc[decisions.gene_id.isin(new_cta - old_cta)].to_dict("records"):
        previous = old_genes[r["gene_id"]]
        current = now_genes[r["gene_id"]]
        new_rows.append(
            {
                "protein": r["protein"],
                "gene_id": r["gene_id"],
                "previous_failure": previous["translation_qualification_reasons"],
                "selected_transcripts": json.dumps(r["selected_transcripts"]),
                "lengths": json.dumps(r["candidate_lengths"]),
                "support": " | ".join(
                    s["transcript_support_level"] or "missing" for s in r["candidate_support"]
                ),
                "selection_basis": r["selection_basis"],
                "cancer_maximum_tpm": current["oncoref_maximum_tpm"],
                "cancer_cohort_at_maximum": current["oncoref_cancer_code"],
                "current_classification": current["classification"],
            }
        )
    write_csv(pd.DataFrame(new_rows), OUT / "new-cta-review.csv")
    added_table = "\n".join(
        f"| {r['protein']} · [{r['gene_id']}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r['gene_id']}) | {r['previous_failure'].replace(' | ', '; ')} | {r['lengths']} | {r['support'].replace(' | ', '; ')} |"
        for r in new_rows
    )
    text = f"""# Reference proteins: selection and unresolved evidence

OncoRef chooses **MANE Select first, Ensembl canonical second**, then requires a
complete primary-locus coding sequence and resolution of conflicting coding evidence.
The current catalog contains **{counts.get("CTA", 0):,} CTA, {counts.get("Unsafe", 0):,} Unsafe
and {counts.get("Other", 0):,} Other proteins**. Policy: `{transcript_selection.POLICY_ID}`.

[Search every gene and transcript choice](table.html) · [All decisions](gene-decisions.csv.gz) ·
[Changed defaults](gene-default-changes.csv.gz) · [Review queue](review-queue.csv.gz) ·
[Classification changes by gene](gene-classification-changes.csv.gz) ·
[Surface-coordinate review](surface-coordinate-review.csv)

![Reference selection, qualification and biological evidence](method.png)

## What each decision establishes

**Reference selection** chooses a representative sequence. MANE Select is a joint
NCBI/Ensembl representative with matching RefSeq and Ensembl sequences. Ensembl
canonical provides the fallback when no primary-locus MANE Select translation is
available. MANE Plus Clinical identifies additional clinically useful transcripts;
it does not replace Select in this ordering. [NCBI methodology](https://www.ncbi.nlm.nih.gov/refseq/MANE/),
[MANE paper](https://doi.org/10.1038/s41586-022-04558-8),
[Ensembl canonical method](https://www.ensembl.org/info/genome/genebuild/canonical.html).

**Sequence qualification** checks the exact annotated CDS: primary chromosome,
qualifying coding biotype, complete start and stop codons, no incomplete or LoF/NMD
annotation, methionine start and uninterrupted AA20 sequence. The full source
inventory retains failed records and alternative translations.

**Production evidence** asks whether the protein is made, and where. Coding annotation
continues to satisfy the existing reference translation rule. Experimental proof is
not newly required for every CTA. Gene RNA, shared peptides and a reference designation
do not establish production of a particular tumor proteoform. Normal-tissue findings
and cancer/reproductive evidence remain separate classification inputs.

## A common method for incomplete or additional information

| Situation | Decision | Information needed to change it |
| --- | --- | --- |
| Complete primary MANE Select, no unresolved conflict | Approve that sequence. Retain all alternatives. | New evidence is assessed against the exact locus and sequence. |
| No primary MANE Select; complete Ensembl canonical | Approve canonical and record the missing MANE designation. | A later MANE release or reviewed sequence-specific evidence. |
| Designated transcript is incomplete, LoF/NMD, interrupted or alternate-only | Keep the candidate; approve no default. Do not choose a shorter or patch sequence silently. | Corrected source annotation or a reviewed qualified alternative in the pinned reference. |
| Complete alternatives but no designated representative, or conflicting designated proteins | Keep the alternatives under review. | A curated designation or explicit review selecting one qualified sequence. |
| Weak, missing or unassessed transcript support | Display the raw support value and flags alongside the designation. | Exact-version transcript evidence; missing and unassessed are different states. |
| Published protein evidence disagrees with the representative | Attach a source-linked hold to the gene or exact sequence. | A review of sequence identity, construct versus endogenous detection, assay specificity and biological context. |
| Additional RNA or shared-peptide evidence | Attach support; retain the current selection and any unresolved holds. | Isoform-resolving evidence is needed to choose an alternative sequence. |
| Strong evidence favors a different complete reference protein | Permit an explicit reviewed preference for its exact sequence hash. | Record rationale, source, reference version and scope; address each prior hold explicitly. |
| A proposed protein is absent from the reference inventory | Keep the external finding; no automatic sequence admission. | Add a versioned sequence/occurrence record and apply the same qualification and evidence review. |

TSL describes full-transcript alignment support: 1 is strongest and 5 lacks a single
supporting transcript. `NA` means unassessed; an absent value means not supplied.
Previous-version assignments remain visible. TSL is not a measure of protein abundance.
MANE does not guarantee the best numeric TSL: the CIMIP3 change illustrates this.
[GENCODE tags and TSL definitions](https://www.gencodegenes.org/pages/data_format.html).
APPRIS and CCDS can add supporting annotation; the pinned Ensembl GTF supplies CCDS
and TSL, but no APPRIS assignments. Unfetched annotations are not counted as support.

Every reviewed finding stores a source URL, review date, reference identity, scope,
limitations and disposition: **support**, **hold**, **resolved**, or **prefer**.
A positive paper does not automatically erase a hold. Sequence preferences must
identify one exact reference protein and cannot bypass completeness checks.

## Findings applied in this release

| Protein | Reason for review | Current classification | Sources |
| --- | --- | --- | --- |
{chr(10).join(case_rows)}

GDPD4 needs reconciliation of the MANE 520-aa and UniProt 623-aa proteins. The
SLCO1B3–SLCO1B7 readthrough needs evidence connecting the selected 748-aa sequence
to the studied 687-aa protein; the native liver ER finding and transfected-cell
surface localization have different scopes. FRG2/FRG2B require resolution of the
experimental coding challenge in the tested systems. These translation findings
do not establish normal-tissue Unsafe status.

## Changes from the previous catalog

The previous release contained {len(old):,} protein identities; this release contains
{len(facts):,}. Identical approved defaults define shared protein identities.
CDY1 and CDY1B now have different representatives and separate expression rows;
RBMY1F and RBMY1J still share the same sequence. PSG11 now has one selected protein
instead of two equal-length alternatives. RNA sums are recomputed within each
measured sample after membership changes.

{int(decisions.default_changed.sum()):,} classified loci change approved default or
move into/out of representative review. At the gene level, {len(new_cta - old_cta)}
enter CTA and {len(old_cta - new_cta)} leave it. These counts distinguish changed
membership from changed admission decisions.

New CTA loci: {", ".join(added) or "none"}.

These newly admitted loci previously failed the longest-sequence rule. The new
representatives pass the annotation, normal-evidence, coverage and CTA association
gates. ARL13A and C16orf90 have TSL 5 assigned to an earlier transcript version;
that support limitation remains recorded. ENSG00000274944 uses a 110-aa canonical
protein with TSL 2. These changes do not add experimental protein-production evidence.
[Full new-CTA review](new-cta-review.csv).

| Protein / locus | Previous translation failure | New length (aa) | Raw TSL |
| --- | --- | --- | --- |
{added_table}

Former CTA loci now outside CTA: {", ".join(removed) or "none"}.

The original 602-protein literature and membrane reports remain historical snapshots.
The coordinate review rechecks their 37 membrane proteins against current approved
sequences: {summary["surface_exact_matches"]} match UniProt exactly. Coordinates are
withheld for unapproved or mismatched sequences. RNA cancer frequency in those
reports remains gene-level evidence.

Protein expression summaries are rebuilt for 142 cancer cohorts against the new
membership registry. NUTM's stored input has negative or infinite values in a field
declared as linear TPM. Its historical protein summaries are marked unavailable
and rejected by the reader; no clipping or transformed-scale summation is used.
Its source repair is tracked in [#604](https://github.com/pirl-unc/oncoref/issues/604)
and [#665](https://github.com/pirl-unc/oncoref/issues/665).

## Reproduce and inspect

The source inventory has {len(translation_qualification.records()):,} Ensembl 115
peptide records. Selection is versioned separately from sequence completeness and
CTA admission. Rebuild qualification, selection, protein membership and expression
before generating the classifications and sequence bundle. The migration driver
pins historical inputs at `{BASELINE_COMMIT}` and reuses the original measured
samples; it does not rerun read quantification.

Previous sequence bundle: `{summary["previous_sequence_bundle_id"]}`.
Current sequence bundle: `{summary["sequence_bundle_id"]}`.

```python
from oncoref import transcript_selection, translation_qualification, cta_selection
choices = transcript_selection.records()
findings = transcript_selection.evidence()
alternatives = translation_qualification.records()
profiles = cta_selection.profiles()
```

`representative_sequence_sha256s` contains approved defaults only.
`representative_candidate_sha256s` retains designated sequences under review.
`translation_qualified` describes sequence quality; `representative_selection_approved`
also requires resolved selection evidence. The sequence bundle retains every source
occurrence, ambiguous sequence and coordinate-bearing valid segment. Unsafe evidence
still takes precedence across occurrences.
"""
    (OUT / "index.md").write_text(text)
    inputs = [
        Path(__file__),
        CURRENT,
        ROOT / "oncoref/transcript_selection.py",
        ROOT / "oncoref/data/sequence-bundle.json",
        *[
            ROOT / "oncoref/data" / name
            for name in [
                "transcript-selection.json",
                "transcript-selection.csv.gz",
                "transcript-evidence.json",
                "translation-qualification.json",
                "cta-literature.json.gz",
            ]
        ],
    ]
    (OUT / "manifest.json").write_text(
        json.dumps(
            {
                "baseline_commit": BASELINE_COMMIT,
                "baseline_sha256": sha256(oldpath),
                "inputs": input_hashes(inputs),
                "outputs": {
                    p.name: sha256(p)
                    for p in sorted(OUT.iterdir())
                    if p.is_file() and p.name != "manifest.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    build()
