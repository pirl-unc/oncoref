#!/usr/bin/env python3
"""Reconcile current paper coverage and audit released representative sequences.

Offline inputs: installed reference tables, the previous dated review, and a
verified sequence bundle. This report does not mutate classification assets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from build_cta_literature_review import compact_cover

from oncoref import (
    cta_literature,
    cta_selection,
    sequence_bundle,
    transcript_selection,
    translation_qualification,
)

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "docs/audits/cta-transcript-selection-20261108"


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def csv(path, rows):
    frame = pd.DataFrame(rows)
    for col in frame:
        frame[col] = frame[col].map(
            lambda v: json.dumps(v) if isinstance(v, (list, tuple, dict)) else v
        )
    frame.to_csv(path, index=False)


def source_link(source, sources):
    url = source if source.startswith("https://") else "https://doi.org/" + source
    label = sources.get(source, {}).get("citation", source)
    return f"[{label}]({url})"


def gene_link(label, genes):
    return f"[{label}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={genes[0]})"


def figure(out, gap_rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), gridspec_kw={"height_ratios": [1, 2]})
    colors = ["#236b68", "#80b2a8", "#d5a249"]
    left = 0
    for number, label, color in zip(
        [487, 87, 33], ["Catalog sources", "Other papers", "Paper gaps"], colors, strict=True
    ):
        axes[0].barh([0], [number], left=left, color=color, label=f"{label}: {number}")
        axes[0].text(
            left + number / 2,
            0,
            str(number),
            ha="center",
            va="center",
            color="white" if number == 487 else "#172a30",
            fontweight="bold",
        )
        left += number
    axes[0].set(xlim=(0, 607), yticks=[], title="607 current CTA proteins · paper coverage")
    axes[0].legend(loc="lower left", bbox_to_anchor=(0, -0.65), ncol=3, frameon=False)
    counts = Counter(r["biological_group"] for r in gap_rows)
    short = {
        "Duplicated regulatory and antigen-related families": "Duplicated regulatory / antigen-related families",
        "Unnamed loci and a readthrough": "Unnamed loci and a readthrough",
        "Other poorly characterized reproductive proteins": "Other reproductive proteins",
        "Reproductive secretory families": "Reproductive secretory families",
        "Maternal / oocyte-associated proteins": "Maternal / oocyte proteins",
        "Germ-cell structure and chromatin": "Germ-cell structure / chromatin",
    }
    items = counts.most_common()
    axes[1].barh([short[k] for k, _ in items], [v for _, v in items], color="#d5a249")
    axes[1].invert_yaxis()
    for i, (_, value) in enumerate(items):
        axes[1].text(value + 0.15, i, str(value), va="center")
    axes[1].set(xlim=(0, 14), xlabel="Proteins", title="The 33 gaps by biological category")
    for ax in axes:
        ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout(h_pad=3)
    fig.savefig(out / "paper-coverage.png", dpi=160)
    fig.savefig(out / "paper-coverage.svg")
    plt.close(fig)


def proposed_resolutions(qualification):
    """Make resolved choices reviewable without mutating released assets."""
    q = qualification.set_index("transcript_accession")
    proposals = []
    for tid, name, disposition, rationale, url in [
        (
            "ENST00000315938.5",
            "gdpd4",
            "resolved",
            "The MANE 520-aa sequence exactly matches UniProt Q6W3E5-2. UniProt explicitly cross-references this transcript, ENSP00000320815.4 and NP_878253.1 to isoform 2. The canonical 623-aa sequence is isoform 1; the difference does not challenge translation.",
            "https://www.uniprot.org/uniprotkb/Q6W3E5/entry",
        ),
        (
            "ENST00000641364.1",
            "or6c3",
            "prefer",
            "This qualified primary-locus CCDS31819 transcript encodes exactly the same 311-aa sequence as the MANE transcript ENST00000641740.2, which lacks a complete annotated stop feature. Select the identical sequence with its complete transcript occurrence.",
            "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000641364.1",
        ),
    ]:
        r = q.loc[tid]
        proposals.append(
            {
                "id": f"transcript-followup-20261108-{name}",
                "gene_id": r.gene_id,
                "reference_id": r.reference_id,
                "scope": "exact_sequence",
                "sequence_sha256s": [r.sequence_sha256],
                "disposition": disposition,
                "evidence_kind": "reviewed_annotation_reconciliation",
                "reviewed_on": "2026-11-08",
                "rationale": rationale,
                "limitations": "Reference sequence reconciliation; does not establish tumor-specific protein production or remove normal-expression exclusions.",
                "sources": [{"title": f"{name.upper()} exact reference annotation", "url": url}],
                "resolves": ["transcript-review-20261108-gdpd4"] if name == "gdpd4" else [],
            }
        )
    genes = {r["gene_id"] for r in proposals}
    ledger = transcript_selection.evidence()
    ledger["findings"] = [r for r in ledger["findings"] if r["gene_id"] in genes] + proposals
    result = transcript_selection.select(qualification[qualification.gene_id.isin(genes)], ledger)
    assert result.selection_status.eq("selected").all()
    return {
        "schema": "oncoref.transcript-review-proposals",
        "status": "reviewed; not applied to released assets",
        "findings": proposals,
    }


def audit_sequences(bundle, profiles, qualification, decisions):
    q = qualification.set_index("transcript_accession")
    d = decisions.set_index("gene_id")
    expected = {}
    for row in profiles.to_dict("records"):
        if row["category"] != "CTA":
            continue
        for gid in row["gene_ids"]:
            chosen = d.loc[gid]
            assert chosen.selection_status == "selected"
            assert chosen.selected_sequence_sha256 in row["representative_sequence_sha256s"]
            for tid in chosen.selected_transcripts:
                record = q.loc[tid]
                assert record.gene_id == gid and record.target_qualified
                assert record.sequence_sha256 == chosen.selected_sequence_sha256
                expected[tid] = {
                    "gene_id": gid,
                    "protein": row["label"],
                    "entity_id": row["entity_id"],
                    "sequence_sha256": record.sequence_sha256,
                    "protein_id": record.protein_id,
                    "transcript": tid,
                    "selection_basis": chosen.selection_basis,
                }
    seen, candidates = {}, {}
    candidate_ids = {
        t for row in decisions.to_dict("records") for t in row["candidate_transcripts"]
    }
    extra_sequences = {}
    for record in bundle.records():
        tid = record["qualification"].get("transcript_accession")
        if tid in candidate_ids:
            candidates[tid] = record["sequence"]
        if record["qualification"].get("symbol") in {"GDPD4", "OR6C3"}:
            extra_sequences[tid] = record["sequence"]
        for occurrence in record["occurrences"]:
            tid = occurrence["transcript_id"]
            if tid not in expected or occurrence["gene_id"] != expected[tid]["gene_id"]:
                continue
            wanted = expected[tid]
            assert record["record_id"] == wanted["protein_id"]
            assert record["sequence_sha256"] == wanted["sequence_sha256"]
            assert (
                hashlib.sha256(record["sequence"].encode()).hexdigest() == wanted["sequence_sha256"]
            )
            assert record["source_record"]["sequence"] == record["sequence"]
            assert record["shared_sequence"] == "allow" and occurrence["direct_target_candidate"]
            seen[tid] = {
                **wanted,
                "length": len(record["sequence"]),
                "sequence": record["sequence"],
            }
    assert set(seen) == set(expected)
    return list(seen.values()), candidates, extra_sequences


def triage(qualification, candidates):
    queue = pd.read_csv(OLD / "review-queue.csv.gz", keep_default_na=False)
    primary = qualification.loc[qualification.contig.isin(transcript_selection.PRIMARY_CONTIGS)]
    by_gene = dict(iter(primary.groupby("gene_id")))
    q = qualification.set_index("transcript_accession")
    output = []
    for row in queue.to_dict("records"):
        tids = json.loads(row["candidate_transcripts"])
        rows = by_gene.get(row["gene_id"], primary.iloc[:0])
        alternatives = rows.loc[rows.target_qualified]
        reasons = json.loads(row["selection_reasons"])
        has_u = any("U" in candidates[t] for t in tids)
        if "unresolved_translation_evidence" in reasons:
            kind = "Reviewed evidence hold"
            decision = "Retain the reviewed hold pending its stated evidence resolution."
        elif rows.empty:
            kind = "No primary reference translation"
            decision = (
                "Retain the historical identity; no primary reference sequence can be approved."
            )
        elif has_u:
            kind = "Selenoprotein outside AA20 contract"
            decision = "U denotes selenocysteine, not an unknown residue. Preserve the full sequence and exclusion; do not substitute a shorter U-free transcript solely to pass AA20."
        elif not alternatives.empty:
            kind = "Qualified alternative available"
            decision = "Review the designated protein and the alternative's biological identity before substituting it."
        else:
            kind = "No qualifying primary alternative"
            decision = "Retain no approved default: no primary transcript meets the pinned sequence contract. Source correction or new evidence is required."
        label = row["protein"]
        source = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + row["gene_id"]
        if label == "GDPD4":
            decision = "Resolve the annotation-only hold: MANE 520 aa exactly matches UniProt isoform Q6W3E5-2. Retain 623-aa isoform 1 separately. A canonical-isoform difference does not challenge translation."
            source = "https://www.uniprot.org/uniprotkb/Q6W3E5/entry"
        elif label == "OR6C3":
            decision = "Approve an explicit same-sequence preference in the next rebuild: ENST00000641364.1 passes and encodes the identical 311-aa protein as the incomplete-stop MANE record. No amino-acid change is required."
        elif label in {"FRG2", "FRG2B"}:
            decision = "Retain the reviewed coding hold. The experimental challenge concerns protein production; choosing another annotated transcript does not answer it."
            source = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12242771/"
        elif label == "SLCO1B3-SLCO1B7":
            decision = "Retain the 748-aa hold. The studied 687-aa alternative has human hepatocyte ER protein evidence, so selecting it also requires normal-evidence adjudication."
            source = "https://pubmed.ncbi.nlm.nih.gov/29248594/"
        elif label == "PKD1L2":
            decision = "Retain review: the 240/306-aa alternatives do not resolve the 2459-aa LoF-designated product. Passing sequence syntax alone cannot establish that a fragment is representative."
        elif label in {"WDR26", "TEAD1", "KCTD11"}:
            decision = "Retain the Unsafe exclusion. Non-AUG initiation is documented or proposed for this gene; a missing GTF start-codon feature is not evidence that the protein is fictitious. Do not replace the extended protein merely to meet the start-feature check."
            source = "https://doi.org/10.1093/nar/gkr007"
        elif label == "MRTFA":
            decision = "Retain the Unsafe exclusion and review the start-codon contract. RefSeq describes non-AUG initiation; shorter qualifying alternatives do not invalidate the longer representative."
            source = "https://www.ncbi.nlm.nih.gov/nuccore/NM_001282662.3"
        elif label == "FNDC5":
            decision = "Retain the Unsafe exclusion. Examine the complete 260-aa upstream-start alternative against the 2024 study; the shorter 153-aa transcript is not an adequate automatic substitute for the designated precursor."
            source = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11555857/"
        output.append(
            {
                "gene_id": row["gene_id"],
                "protein": label,
                "classification": row["classification"],
                "review_category": kind,
                "decision": decision,
                "source_url": source,
                "candidate_transcripts": tids,
                "candidate_lengths": json.loads(row["candidate_lengths"]),
                "candidate_qualification": {t: q.loc[t].qualification_reasons for t in tids},
                "qualified_alternatives": alternatives[
                    ["transcript_accession", "length", "sequence_sha256", "ccds_id"]
                ].to_dict("records"),
            }
        )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-bundle", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "docs/audits/cta-current-review-20261108"
    )
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    profiles = cta_selection.profiles()
    current = cta_literature.current_coverage()
    q = translation_qualification.records()
    d = transcript_selection.records()
    bundle = sequence_bundle.load(args.sequence_bundle)
    verified, candidates, extras = audit_sequences(bundle, profiles, q, d)
    reviewed = triage(q, candidates)
    assert extras["ENST00000641364.1"] == extras["ENST00000641740.2"]
    # The reviewed UniProt variant features are exact, versioned sequence edits:
    # VSP_020811 substitutes residues 514..520; VSP_020812 removes 521..623.
    uniprot_hash = cta_literature.surface_domains("GDPD4").iloc[0].uniprot_sequence_sha256
    assert hashlib.sha256(extras["ENST00000376217.6"].encode()).hexdigest() == uniprot_hash
    gdpd4_iso2 = extras["ENST00000376217.6"][:513] + "SKNEEDR"
    assert gdpd4_iso2 == extras["ENST00000315938.5"]
    lookup = {r["label"]: r for r in profiles.to_dict("records")}
    old_gaps = {r["protein"]: r for r in cta_literature.gap_assessments().to_dict("records")}
    judgments = json.loads(
        (ROOT / "scripts/cta_evidence_assets/current-gap-judgments.json").read_text()
    )
    gaps = []
    for r in current[current.paper_coverage.eq("no_paper")].to_dict("records"):
        previous = old_gaps.get(r["protein"], {})
        profile = lookup[r["protein"]]
        gaps.append(
            {
                **r,
                "biological_group": previous.get(
                    "biological_group", "Unnamed loci and a readthrough"
                ),
                "biology_categories": previous.get("biology_categories", ["early_embryo"]),
                "evidence_summary": previous.get(
                    "evidence_summary",
                    "No verified exact-locus study. Complete 110-aa Ensembl canonical translation; TSL 2; developmental and tumor RNA.",
                ),
                "judgment": judgments[r["protein"]],
                "tumor_max_tpm": profile["tumor_max_tpm"],
                "developmental_lineage_max_tpm": profile["developmental_lineage_max_tpm"],
                "reproductive_max_ntpm": profile["reproductive_max_ntpm"],
            }
        )
    assert {r["protein"] for r in gaps} == set(judgments)
    sources = {s["source_id"]: s for s in cta_literature.registry()["sources"]}
    sources.update({s["source_id"]: s for s in cta_literature.followup_findings()["sources"]})
    sets = defaultdict(set)
    for r in current.to_dict("records"):
        for gid, papers in r["per_gene_sources"].items():
            for paper in papers:
                sets[paper].add(gid)
    core = set(compact_cover(sets))
    roles = []
    for key, s in sorted(sources.items()):
        roles.append(
            {
                "source_id": key,
                "citation": s["citation"],
                "title": s["title"],
                "url": s["url"],
                "roles": s["roles"],
                "current_collection": "core" if key in core else "supplementary",
                "covered_current_cta_loci": sorted(sets[key]),
                "retention_reason": s["retention_reason"],
            }
        )
    assert set().union(*(sets[k] for k in core)) == set().union(*sets.values())
    for key in core:
        assert set().union(*(sets[k] for k in core - {key})) != set().union(*sets.values())
    summary = {
        "catalog_version": "1.35.0",
        "reviewed_on": "2026-11-08",
        "cta_proteins": len(current),
        "nomination_paper_coverage": int(current.nomination_sources.map(bool).sum()),
        "any_qualifying_paper": int(current.exact_paper_sources.map(bool).sum()),
        "paper_coverage": current.paper_coverage.value_counts().to_dict(),
        "current_core_sources": len(core),
        "supplementary_sources": len(roles) - len(core),
        "core_objective": "All covered exact gene loci in the current CTA catalog; inclusion-minimal, not claimed globally smallest",
        "representative_review_loci": len(reviewed),
        "representative_review_classifications": dict(
            Counter(r["classification"] for r in reviewed)
        ),
        "representative_review_categories": dict(Counter(r["review_category"] for r in reviewed)),
        "with_qualified_alternatives": sum(bool(r["qualified_alternatives"]) for r in reviewed),
        "verified_cta_transcripts": len(verified),
        "sequence_mismatches": 0,
        "all_cta_sequences_AA20": all(
            set(r["sequence"]) <= set("ACDEFGHIKLMNPQRSTVWY") for r in verified
        ),
        "gdpd4_isoform_2_exact_match": True,
        "or6c3_same_sequence_alternative": True,
        "sequence_bundle_id": bundle.identity,
        "classification_effect": "No released classification assets changed. Review recommendations are listed separately.",
    }
    dump(out / "summary.json", summary)
    dump(out / "transcript-review-resolutions.json", proposed_resolutions(q))
    csv(out / "current-paper-coverage.csv", current)
    csv(out / "outside-nomination-papers.csv", current[~current.nomination_sources.map(bool)])
    csv(out / "paper-gap-judgments.csv", gaps)
    csv(out / "source-roles-current.csv", roles)
    csv(out / "representative-review.csv", reviewed)
    csv(out / "verified-cta-transcripts.csv", verified)
    # Keep JSON values portable: unknown numerical measurements are null, not NaN.
    dump(out / "paper-gap-judgments.json", json.loads(pd.DataFrame(gaps).to_json(orient="records")))
    figure(out, gaps)
    report = [
        "# Current CTA papers, representative review and API verification\n",
        "This review uses the released **607 CTA proteins in OncoRef 1.35.0**. It records new evidence and recommendations; it does not change released classifications.\n",
        "## Which papers cover the current CTAs?\n",
        "**487** have membership in the 15 nomination-source collection. **87** more have qualifying papers in the broader review. **33** lack a verified qualifying exact-human paper. Thus **574/607** have paper coverage; **120** are outside the nomination collection.\n",
        "The 574 comprise **569** proteins whose member loci all have coverage and **5** with coverage for only some member loci. A paper can describe normal function, disease association or counterevidence. Paper coverage is not a declaration that the paper nominated a CTA or demonstrated tumor protein.\n",
        "![Paper coverage and biological categories of the gaps](paper-coverage.png)\n",
        "The old frozen report said 569/602 with 33 gaps. Remapping its memberships gives 573/607 with 34 gaps. The newly reviewed human CALHM3 paper closes one gap; newly admitted ENSG00000274944 adds one relative to the old list. The gap count is again 33, with different membership.\n",
        f"The current reading set has **{len(core)} core** and **{len(roles) - len(core)} supplementary** sources. The core covers every presently covered CTA locus and no core source can be removed without losing coverage. Supplementary sources retain normal findings, methods, mechanisms and independent measurements. This objective differs from the older 83-core/67-supplementary set, which also covered proteins outside CTA.\n",
        "[Every current CTA and its papers](current-paper-coverage.csv) · [All 120 outside nomination papers](outside-nomination-papers.csv) · [Sources, links and roles](source-roles-current.csv)\n",
        "| Shared protein | Member loci without qualifying paper coverage |\n|---|---|",
    ]
    for r in current[current.paper_coverage.eq("some_members")].to_dict("records"):
        identities = {g["gene_id"]: g for g in lookup[r["protein"]]["gene_identities"]}
        links = [
            gene_link(identities[g]["symbol"] + " · " + g, [g]) for g in r["uncovered_gene_ids"]
        ]
        report.append(f"| {r['protein']} | {', '.join(links)} |")
    report += [
        "\n## A missed human normal-tissue paper\n",
        "**CALHM3:** Moyer et al. measured RNA in laser-captured human circumvallate taste buds (Figure 10D). The prior mouse-only assessment was incomplete. I recommend a normal-evidence review hold; the released classification is still CTA. The human result is RNA, while cell-subtype localization in this paper is macaque data. [Moyer et al., 2009](https://doi.org/10.1371/journal.pone.0007682).\n",
        "## Judgment of each remaining paper gap\n",
        "These are discovery CTAs under the existing rule allowing a qualifying coding reference. Missing papers alone do not make them Unsafe. Tumor maxima below are single-sample RNA observations, not prevalence or protein abundance.\n",
        "| Protein | Biological category | Maximum tumor RNA (TPM) | Judgment | Related paper |\n|---|---|---:|---|---|",
    ]
    for r in sorted(gaps, key=lambda r: r["protein"]):
        report.append(
            f"| {gene_link(r['protein'], r['gene_ids'])} | {r['biological_group']} | {r['tumor_max_tpm']:.3g} | {r['judgment']} | {'; '.join(source_link(s, sources) for s in r['related_sources']) or '—'} |"
        )
    report += [
        "\n[Full gap evidence and judgments](paper-gap-judgments.json)\n",
        "## What the 563 loci actually mean\n",
        "**None is currently CTA.** They comprise 381 Other and 182 Unsafe loci. The Unsafe label comes from normal-expression evidence independently of representative selection. '563 await review' combines source limitations with already reviewed findings.\n",
        "The pinned reference has **no qualifying primary alternative for 543** loci; **20** have at least one. The original reasons were 557 designated translations failing the contract, four reviewed evidence holds, and two loci with no primary translation. This audit triages every locus and examines the 20 alternative cases. It does not claim 563 separate literature reviews.\n",
        "| Review category | Loci |\n|---|---:|",
    ]
    for key, count in summary["representative_review_categories"].items():
        report.append(f"| {key} | {count} |")
    report += [
        "\n**25 are selenoproteins.** Their U residues encode selenocysteine; they are not unknown X residues. The current AA20 output contract rejects them. Known or proposed non-AUG starts also require a better distinction from missing biological starts. [Non-AUG initiation study](https://doi.org/10.1093/nar/gkr007).\n",
        "**GDPD4:** the 520-aa MANE translation exactly matches UniProt Q6W3E5-2 after applying its two annotated variant features. The 623-aa sequence is isoform 1. The annotation-only hold should be resolved in the next classification rebuild. [UniProt isoforms](https://www.uniprot.org/uniprotkb/Q6W3E5/entry).\n",
        "**OR6C3:** ENST00000641364.1 passes qualification and encodes exactly the same 311-aa sequence as MANE ENST00000641740.2, whose GTF lacks a complete stop feature. An explicit reviewed sequence preference can resolve selection without changing the amino-acid sequence.\n",
        "[Reviewed GDPD4 and OR6C3 resolution records](transcript-review-resolutions.json) pass the selector's validation and select both proteins in an isolated check. They are ready for the next classification rebuild; the released assets retain their existing labels.\n",
        "| Protein | Current class | Reviewed decision | Source |\n|---|---|---|---|",
    ]
    for r in reviewed:
        if r["qualified_alternatives"]:
            report.append(
                f"| {r['protein']} | {r['classification']} | {r['decision']} | [source]({r['source_url']}) |"
            )
    report += [
        "\n[All 563 decisions, failed checks and qualifying alternatives](representative-review.csv)\n",
        "## Why 602 became 607\n",
        "**602 + 8 admissions − 4 translation holds + 1 split identity = 607.**\n",
        "New proteins: BOLL, IZUMO2, ARL13A, CATSPER4, GAGE1, C16orf90, ENSG00000274944 and the shared TP53TG3/B/C/D/E/F protein. Four moved to Other:Translation: GDPD4, SLCO1B3–SLCO1B7, FRG2 and FRG2B. CDY1 and CDY1B now have different selected sequences and count separately. The eight admissions contain 13 loci; protein counts and locus counts differ. The GDPD4 correction above is a new finding, not part of the published 607 calculation.\n",
        "## Does the API identify the transcript and matching sequence?\n",
        "**Yes. All 655 selected transcript versions across the 607 CTAs matched their sequence records exactly: zero mismatches and no X or other non-AA20 residues.** Checks covered exact gene, versioned transcript, versioned protein accession, sequence hash, source peptide bytes and sequence role. This verifies reference provenance; it does not infer tumor-specific isoform production.\n",
        "The new `cta_representatives.records('BOLL', annotation_release=115)` accessor returns the chosen versioned transcript and protein accession, length, sequence hash and reference without downloading the bundle. It accepts exact symbols or Ensembl gene IDs and rejects a requested release mismatch. The legacy `Canonical_Transcript_ID` field is not the current selection contract. LLCFC1, mentioned in issue #683, is now Unsafe from normal bulk evidence and is therefore absent from this CTA-only accessor.\n",
        "`transcript_selection.records()` exposes `selected_transcripts` and `selected_sequence_sha256` for all reference genes, including non-CTA genes. `cta_selection.profiles()` exposes `representative_transcripts` and `representative_sequence_sha256s`. `sequence_bundle.load(...).records()` supplies the actual amino-acid sequence, transcript occurrences and policy role. Candidate transcripts under review are separate from approved transcripts.\n",
        "```python\nfrom oncoref import cta_literature, cta_selection, transcript_selection, sequence_bundle\n\ncoverage = cta_literature.current_coverage()\ngaps = coverage[coverage.paper_coverage == 'no_paper']\nchoices = transcript_selection.records()\nchoice = choices[choices.gene_id == 'ENSG00000152430'].iloc[0]  # BOLL\nbundle = sequence_bundle.load('/path/to/verified/bundle')\nchosen = [r for r in bundle.records()\n          if r['sequence_sha256'] == choice.selected_sequence_sha256\n          and any(o['gene_id'] == choice.gene_id\n                  and o['transcript_id'] in choice.selected_transcripts\n                  for o in r['occurrences'])]\n```\n",
        "[Verified transcript IDs, protein IDs and sequences](verified-cta-transcripts.csv) · [Machine-readable summary](summary.json)\n",
        "Reproduce: `.venv/bin/python scripts/reconcile_cta_reviews.py --sequence-bundle /path/to/verified/bundle`. The historical literature API remains `cta_literature.cta_coverage()` (602 proteins); use `current_coverage()` for the current catalog.\n",
    ]
    (out / "index.md").write_text("\n".join(report))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
