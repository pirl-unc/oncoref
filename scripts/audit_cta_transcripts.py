#!/usr/bin/env python3
"""Audit released CTA defaults against pinned MANE/Ensembl representatives.

This is a migration preview. It does not assign sequence roles, change protein
membership, or replace the released translation/classification policies.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from oncoref import cta_literature, cta_selection
from oncoref.translation_qualification import PRIMARY_CONTIGS, records

ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "scripts/cta_evidence_assets/cta-transcript-support.json.gz"
OUTPUT = ROOT / "docs/audits/cta-transcripts-20261108"
PREVIEW_POLICY = "oncoref.transcript-selection-preview.2026-11-08.1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    data = (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode()
    path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)


def import_support(gtf: Path, output: Path = SUPPORT):
    """Extract exact transcript-version support from the qualified reference GTF."""
    meta = json.loads((ROOT / "oncoref/data/translation-qualification.json").read_text())
    expected = meta["reference"]["content_provenance"]["gtf"]
    if digest(gtf) != expected:
        raise ValueError("GTF hash differs from the pinned qualification reference")
    genes = {g for ids in cta_selection.candidates().gene_ids for g in ids}
    table = records()
    wanted = table.loc[table.gene_id.isin(genes)]
    tids = set(wanted.transcript_id)
    rows = {}
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "transcript":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)";', fields[8]))
            tid = attrs["transcript_id"]
            if tid not in tids:
                continue
            if tid in rows:
                raise ValueError(f"Repeated transcript: {tid}")
            rows[tid] = {
                key: attrs.get(key, "")
                for key in (
                    "gene_id",
                    "transcript_version",
                    "transcript_source",
                    "transcript_support_level",
                    "ccds_id",
                )
            }
    if set(rows) != tids:
        raise ValueError("GTF does not cover every CTA-locus protein transcript")
    write_json(
        output,
        {
            "schema": "oncoref.cta-transcript-support",
            "schema_version": 1,
            "reference_id": meta["reference"]["reference_id"],
            "gtf_sha256": expected,
            "gtf_url": "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz",
            "records": dict(sorted(rows.items())),
        },
    )


def proposal(transcripts: pd.DataFrame) -> dict:
    """Prefer an explicit representative; never resolve different proteins by length.

    A failed MANE/canonical record remains a review case. A lower-priority or
    alternate-locus transcript cannot silently replace it. All source records
    remain available in the separate transcript inventory.
    """
    primary = transcripts.loc[transcripts.contig.isin(PRIMARY_CONTIGS)]
    for tag, basis in (("MANE_Select", "MANE Select"), ("Ensembl_canonical", "Ensembl canonical")):
        designated = primary.loc[primary.tags.map(lambda value, tag=tag: tag in value)]
        if designated.empty:
            continue
        qualified = designated.loc[designated.target_qualified]
        hashes = sorted(set(qualified.sequence_sha256))
        status = (
            "selected"
            if len(hashes) == 1
            else (
                "review: representative fails qualification"
                if not hashes
                else "review: conflicting representative proteins"
            )
        )
        return {
            "basis": basis,
            "status": status,
            "sequence_sha256": hashes[0] if len(hashes) == 1 else None,
            "designated_transcripts": sorted(designated.transcript_accession),
            "qualified_transcripts": sorted(qualified.transcript_accession),
        }
    return {
        "basis": None,
        "status": "review: no designated primary representative",
        "sequence_sha256": None,
        "designated_transcripts": [],
        "qualified_transcripts": [],
    }


def support_summary(transcripts, support):
    """Best numeric support among transcripts encoding these exact protein bytes.

    Missing and NA stay distinct. Prior-version qualifications remain in the
    transcript inventory; this summary is not an experimental translation score.
    """
    raw = [support[t]["transcript_support_level"] for t in transcripts.transcript_id]
    numeric = [int(v[0]) for v in raw if v and v[0] in "12345"]
    return (
        str(min(numeric))
        if numeric
        else ("NA" if any(v.startswith("NA") for v in raw) else "missing")
    )


def csv(path, rows):
    frame = pd.DataFrame(rows)
    for column in frame:
        frame[column] = frame[column].map(
            lambda v: json.dumps(v, separators=(",", ":")) if isinstance(v, (dict, list)) else v
        )
    frame.to_csv(path, index=False, lineterminator="\n")


def build(output=OUTPUT, support_path=SUPPORT):
    output.mkdir(parents=True, exist_ok=True)
    profiles = cta_selection.profiles()
    ctas = profiles.loc[profiles.category.eq("CTA")].sort_values("label")
    table = records()
    support_input = json.loads(gzip.decompress(support_path.read_bytes()))
    qualification = json.loads((ROOT / "oncoref/data/translation-qualification.json").read_text())
    if (
        support_input["reference_id"] != qualification["reference"]["reference_id"]
        or support_input["gtf_sha256"] != qualification["reference"]["content_provenance"]["gtf"]
    ):
        raise ValueError("Transcript support/reference identity mismatch")
    support = support_input["records"]
    genes = {g for ids in ctas.gene_ids for g in ids}
    locus_table = table.loc[table.gene_id.isin(genes)].copy()
    if set(support) != set(locus_table.transcript_id):
        raise ValueError("Transcript support coverage differs from CTA-locus inventory")
    for r in locus_table.itertuples():
        extra = support[r.transcript_id]
        if extra["gene_id"] != r.gene_id or extra["transcript_version"] != str(
            r.transcript_version
        ):
            raise ValueError(f"Transcript support/version mismatch: {r.transcript_accession}")
    by_gene = dict(iter(locus_table.groupby("gene_id", sort=False)))
    previews = {g: proposal(t) for g, t in by_gene.items()}
    surface = cta_literature.surface_domains().set_index("entity_id").to_dict("index")
    counts = Counter()
    protein_rows, locus_rows, transcript_rows = [], [], []
    new_owners = defaultdict(set)
    for r in ctas.itertuples():
        current_hashes = set(r.representative_sequence_sha256s)
        subset = locus_table.loc[locus_table.gene_id.isin(r.gene_ids)]
        current = subset.loc[subset.sequence_sha256.isin(current_hashes)]
        proposed_hashes = {previews[g]["sequence_sha256"] for g in r.gene_ids}
        proposed_hashes.discard(None)
        chosen = subset.loc[subset.sequence_sha256.isin(proposed_hashes)]
        missing = any(previews[g]["sequence_sha256"] is None for g in r.gene_ids)
        split = not missing and len(proposed_hashes) > 1
        changed = current_hashes != proposed_hashes or missing
        basis = sorted({previews[g]["basis"] or "none" for g in r.gene_ids})
        selected_tags = {tag for tags in current.tags for tag in tags}
        flags = []
        if any(previews[g]["sequence_sha256"] not in current_hashes for g in r.gene_ids):
            flags.append("Replace longest-only default")
        if current_hashes - proposed_hashes and proposed_hashes <= current_hashes:
            flags.append("Remove extra tied default")
        if split:
            flags.append("Rebuild shared-protein membership")
        if missing:
            flags.append("No qualified designated default")
        if "Ensembl canonical" in basis:
            flags.append("No MANE for one or more loci")
        best_support = support_summary(current, support)
        if best_support in ("4", "5"):
            flags.append("Limited full-transcript alignment support")
        if best_support == "missing":
            flags.append("Transcript support not supplied")
        proposed_support = support_summary(chosen, support)
        if (
            best_support.isdigit()
            and proposed_support.isdigit()
            and int(proposed_support) > int(best_support)
        ):
            flags.append("Proposed transcript has weaker alignment support")
        no_strong_tag = not selected_tags.intersection(
            {"MANE_Select", "Ensembl_canonical", "CCDS", "gencode_primary"}
        )
        if no_strong_tag:
            flags.append("Basic-only selected annotation")
        uniprot = surface.get(r.entity_id)
        uniprot_match = proposed_hashes == {uniprot["uniprot_sequence_sha256"]} if uniprot else None
        if uniprot_match is False:
            flags.append("Proposed sequence differs from surface-reference UniProt")
        counts.update(flags)
        counts[f"current_TSL_{best_support}"] += 1
        counts[f"proposed_TSL_{support_summary(chosen, support)}"] += 1
        counts["changed_defaults" if changed else "unchanged_defaults"] += 1
        row = {
            "protein": r.label,
            "entity_id": r.entity_id,
            "gene_ids": r.gene_ids,
            "default_change": changed,
            "review_flags": flags,
            "current_lengths": sorted(set(map(int, current.length))),
            "proposed_lengths": sorted(set(map(int, chosen.length))),
            "current_sequence_sha256s": sorted(current_hashes),
            "proposed_sequence_sha256s": sorted(proposed_hashes),
            "proposal_basis": basis,
            "shared_identity_splits": split,
            "current_best_tsl": best_support,
            "proposed_best_tsl": proposed_support,
            "surface_uniprot_accession": uniprot["accession"] if uniprot else None,
            "surface_uniprot_url": uniprot["uniprot_url"] if uniprot else None,
            "surface_uniprot_length": int(uniprot["uniprot_length"]) if uniprot else None,
            "surface_uniprot_sequence_sha256": uniprot["uniprot_sequence_sha256"]
            if uniprot
            else None,
            "proposed_matches_surface_uniprot": uniprot_match,
            "proposed_transcripts": sorted(
                {t for g in r.gene_ids for t in previews[g]["qualified_transcripts"]}
            ),
        }
        protein_rows.append(row)
        for g in r.gene_ids:
            p = previews[g]
            old = current.loc[current.gene_id.eq(g)]
            new = subset.loc[subset.gene_id.eq(g) & subset.sequence_sha256.eq(p["sequence_sha256"])]
            if p["sequence_sha256"]:
                new_owners[p["sequence_sha256"]].add(r.entity_id)
            locus_rows.append(
                {
                    "protein": r.label,
                    "gene_id": g,
                    **p,
                    "changes_protein": p["sequence_sha256"] not in current_hashes,
                    "current_transcripts": sorted(old.transcript_accession),
                    "current_lengths": sorted(set(map(int, old.length))),
                    "proposed_lengths": sorted(set(map(int, new.length))),
                    "ensembl_url": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={g}",
                }
            )
        for t in subset.itertuples():
            preview = previews[t.gene_id]
            transcript_rows.append(
                {
                    "protein": r.label,
                    "gene_id": t.gene_id,
                    "transcript": t.transcript_accession,
                    "protein_accession": t.protein_id,
                    "length": int(t.length),
                    "sequence_sha256": t.sequence_sha256,
                    "contig": t.contig,
                    "biotype": t.transcript_biotype,
                    "target_qualified": bool(t.target_qualified),
                    "qualification_reasons": t.qualification_reasons,
                    "tags": t.tags,
                    "current_default_sequence": t.sequence_sha256 in current_hashes,
                    "proposed_default_transcript": t.transcript_accession
                    in preview["qualified_transcripts"],
                    **{
                        k: v
                        for k, v in support[t.transcript_id].items()
                        if k not in ("gene_id", "transcript_version")
                    },
                    "ensembl_url": f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={t.transcript_accession}",
                }
            )
    counts["cta_proteins"] = len(ctas)
    counts["cta_loci"] = len(genes)
    counts["cta_locus_translations"] = len(transcript_rows)
    counts["changed_loci"] = sum(r["changes_protein"] for r in locus_rows)
    counts["proposed_merges_between_current_ctas"] = sum(len(v) > 1 for v in new_owners.values())
    csv(output / "protein-defaults.csv", protein_rows)
    csv(output / "locus-defaults.csv", locus_rows)
    csv(output / "transcript-inventory.csv", transcript_rows)
    write_json(
        output / "audit.json.gz",
        {
            "preview_policy": PREVIEW_POLICY,
            "reference_id": qualification["reference"]["reference_id"],
            "counts": dict(counts),
            "proteins": protein_rows,
            "loci": locus_rows,
            "transcripts": transcript_rows,
        },
    )
    build_report(output, counts, protein_rows)
    plot_support(output, counts)
    inputs = [
        Path(__file__),
        support_path,
        ROOT / "oncoref/translation_qualification.py",
        ROOT / "oncoref/cta_selection.py",
        ROOT / "oncoref/data/translation-qualification.json",
        ROOT / "oncoref/data/translation-qualification.csv.gz",
        ROOT / "oncoref/data/cta-selection.json",
        ROOT / "oncoref/data/cta-selection-human.csv.gz",
        ROOT / "oncoref/data/cta-selection-profiles-human.csv.gz",
        ROOT / "oncoref/data/cta-literature.json.gz",
    ]
    write_json(
        output / "manifest.json",
        {
            "preview_policy": PREVIEW_POLICY,
            "counts": dict(counts),
            "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
            "outputs": {
                p.name: digest(p)
                for p in sorted(output.iterdir())
                if p.name != "manifest.json" and p.is_file()
            },
        },
    )
    print(json.dumps(dict(counts), indent=2))
    return counts


def plot_support(output, counts):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    with plt.rc_context({"svg.hashsalt": PREVIEW_POLICY, "font.size": 10}):
        fig, ax = plt.subplots(figsize=(10, 5.7))
        labels = [
            "1\nFull mRNA",
            "2\nMixed support",
            "3\nSingle EST",
            "4\nSuspect EST",
            "5\nNo full match",
            "NA\nNot assessed",
            "Missing\nNot supplied",
        ]
        keys = ["1", "2", "3", "4", "5", "NA", "missing"]
        x = np.arange(len(keys))
        for shift, prefix, color, label in [
            (-0.2, "current", "#a85f2a", "Released longest-first defaults"),
            (0.2, "proposed", "#256b80", "MANE-first preview"),
        ]:
            bars = ax.bar(
                x + shift,
                [counts[f"{prefix}_TSL_{k}"] for k in keys],
                0.36,
                color=color,
                label=label,
            )
            ax.bar_label(bars, padding=3, fontsize=9)
        ax.set_xticks(x, labels)
        ax.set_ylabel("CTA proteins")
        ax.set_ylim(
            0, max(counts[f"{p}_TSL_{k}"] for p in ["current", "proposed"] for k in keys) * 1.2
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, loc="upper right")
        fig.suptitle(
            "Transcript-alignment support for the selected protein",
            x=0.09,
            ha="left",
            weight="bold",
        )
        fig.text(
            0.09,
            0.89,
            f"{counts['cta_proteins']} released CTA records · Ensembl 115 · alternatives retained",
        )
        fig.text(
            0.09,
            0.055,
            "Best numeric TSL across exact-sequence transcripts and member loci; prior-version tags retained in the CSV.\nTSL measures RNA-model support. It does not measure protein detection or tumor abundance.",
            fontsize=9,
            color="#444444",
        )
        fig.subplots_adjust(left=0.09, right=0.98, top=0.83, bottom=0.2)
        fig.savefig(output / "transcript-support.svg", metadata={"Date": None})
        fig.savefig(output / "transcript-support.png", dpi=160)
        plt.close(fig)


def build_report(output, counts, rows):
    """Generate the inventory summary; interpretation lives in the linked review."""
    changes = [r for r in rows if r["default_change"]]
    lines = [
        "# CTA transcript defaults: complete comparison",
        "",
        f"The {counts['cta_proteins']} released CTAs contain {counts['cta_loci']} gene loci. "
        f"A MANE-first comparison changes the selected sequence set for **{len(changes)} proteins**: "
        f"{counts['Replace longest-only default']} replace a longest-only choice and "
        f"{counts['Remove extra tied default']} loses an additional tied sequence.",
        "",
        "[Read the findings and selection rule](index.md). "
        "[All proteins](protein-defaults.csv) · [Every locus](locus-defaults.csv) · "
        "[Every translated transcript, including alternatives](transcript-inventory.csv) · "
        "[Full JSON](audit.json.gz) · [Reproduction manifest](manifest.json)",
        "",
        "This comparison previews a new choice of representative. It does not change the released CTA list or sequence permissions.",
        "",
        "| Protein | Current length (aa) | Proposed length (aa) | Proposed transcript(s) |",
        "|---|---:|---:|---|",
    ]
    for r in changes:
        link = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r['gene_ids'][0]}"
        transcripts = ", ".join(
            f"[{t}](https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={t})"
            for t in r["proposed_transcripts"]
        )
        label = r["protein"] + (" (split)" if r["shared_identity_splits"] else "")
        if "Remove extra tied default" in r["review_flags"]:
            label += " (remove tie)"
        lines.append(
            f"| [{label}]({link}) | {', '.join(map(str, r['current_lengths']))} | {', '.join(map(str, r['proposed_lengths']))} | {transcripts} |"
        )
    lines.extend(
        [
            "",
            "## Transcript support",
            "",
            "Best numeric support among transcripts encoding the selected protein sequence, pooled across member loci. "
            "A strong value at one member does not establish support at every other member; use the transcript inventory for locus-specific review. "
            "Prior-version annotations remain verbatim in that inventory. Missing and NA are separate.",
            "",
            "| Support | Current proteins | Proposed proteins |",
            "|---|---:|---:|",
        ]
    )
    labels = {
        "1": "1: supported by a nonsuspect mRNA",
        "2": "2: suspect mRNA or multiple ESTs",
        "3": "3: single EST",
        "4": "4: suspect EST",
        "5": "5: no single transcript supports the model",
        "NA": "NA: not assessed",
        "missing": "Not supplied",
    }
    for key, label in labels.items():
        lines.append(
            f"| {label} | {counts[f'current_TSL_{key}']} | {counts[f'proposed_TSL_{key}']} |"
        )
    lines.extend(
        [
            "",
            "Definitions: [GENCODE GTF documentation](https://www.gencodegenes.org/pages/data_format.html). "
            "TSL measures transcript-alignment support, not protein detection. NA is not a failed assay.",
            "",
        ]
    )
    (output / "comparison.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--import-gtf", type=Path, help="Extract the pinned support input before building"
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.import_gtf:
        import_support(args.import_gtf)
    build(args.output)
