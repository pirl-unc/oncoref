#!/usr/bin/env python3
"""Publish a verified, source-preserving sequence bundle from a pinned ORFdata snapshot."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oncoref import target_selection, transcript_selection, translation_qualification  # noqa: E402
from oncoref.sequence_bundle import (  # noqa: E402
    ROLE_ORDER,
    ROLE_POLICY_ID,
    SCHEMA,
    content_id,
    digest,
    load,
    target_selection_dispositions,
)

AUDIT_DIR = ROOT / "docs/audits/sequence-bundle-20261121"


def _resolve_source_inputs(source: Path) -> tuple[dict, str, Path]:
    if (source / "reference.json").is_file():
        meta = json.loads((source / "reference.json").read_text())
        proteome_path = source / "proteome.json.gz"
    elif (source / "manifest.json").is_file():
        meta = json.loads((source / "manifest.json").read_text())
        proteome_path = source / "orfdata-proteome.json.gz"
    else:
        raise FileNotFoundError(
            f"{source} must contain reference.json+proteome.json.gz or manifest.json+orfdata-proteome.json.gz"
        )
    return meta["reference"], meta["orfdata_identity"], proteome_path


def build(source: Path, output: Path, release: str, *, asset_tag: str | None = None):
    from orfdata import Proteome

    output.mkdir(parents=True, exist_ok=True)
    reference, expected_orfdata_identity, proteome_path = _resolve_source_inputs(source)
    qual_meta = json.loads((ROOT / "oncoref/data/translation-qualification.json").read_text())
    if reference != qual_meta["reference"]:
        raise ValueError("Translation qualification and sequence references differ")
    proteome = Proteome.load(proteome_path)
    if proteome.identity != expected_orfdata_identity or (
        proteome.identity != qual_meta["orfdata_identity"]
    ):
        raise ValueError("Source ORFdata identity mismatch")
    qualification = translation_qualification.records().set_index("protein_id")
    if set(qualification.index) != set(proteome.by_accession()):
        raise ValueError("Qualification does not cover exactly the source protein inventory")

    decisions = target_selection_dispositions()
    per_gene = defaultdict(list)
    for row in decisions.to_dict("records"):
        row["representative_sequence_sha256s"] = set(row["representative_sequence_sha256s"])
        for gene in row["gene_ids"]:
            per_gene[gene].append(row)

    roles = Counter()
    seen_genes = set()
    no_segments = []
    occurrence_count = 0
    target_records = 0
    cancer_antigen_records = 0
    reproductive_antigen_records = 0

    with (
        (output / "records.jsonl.gz").open("wb") as raw,
        gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped,
    ):
        for accession, protein in sorted(proteome.by_accession().items()):
            original = protein.to_dict()
            sequence = protein.sequence
            seqhash = digest(sequence.encode())
            quality = qualification.loc[accession].to_dict()
            if seqhash != quality["sequence_sha256"]:
                raise ValueError("Quality sequence mismatch")
            occurrences = []
            for i, native_occurrence in enumerate(protein.occurrences):
                occ = asdict(native_occurrence)
                gene = occ["gene_id"]
                seen_genes.add(gene)
                matches = per_gene.get(gene, [])
                reasons = []
                role = "review"
                candidate = False
                cancer_candidate = False
                reproductive_candidate = False
                if any(m["category"] == "Unsafe" for m in matches):
                    role = "exclude"
                    reasons.append("applicable_gene_level_exclusion")
                elif not matches:
                    reasons.append("no_current_decision_for_exact_locus")
                elif not quality["target_qualified"]:
                    reasons.extend(quality["qualification_reasons"])
                elif all(
                    m["shared_sequence"] == "allow"
                    and seqhash in m["representative_sequence_sha256s"]
                    for m in matches
                ):
                    role = "allow"
                    reasons.append("qualified_assessed_representative")
                    candidate = any(bool(m["direct_target_candidate"]) for m in matches)
                    cancer_candidate = any(bool(m["cancer_antigen_candidate"]) for m in matches)
                    reproductive_candidate = any(
                        bool(m["reproductive_antigen_candidate"]) for m in matches
                    )
                else:
                    if any(m["classification"] == "Other:Unscreened tissue" for m in matches):
                        reasons.append("unscreened_home_tissue_requires_review")
                    else:
                        reasons.append("unassessed_isoform_or_decision_requires_review")
                occurrences.append(
                    {
                        **occ,
                        "occurrence_id": f"{accession}:{i}",
                        "shared_sequence": role,
                        "direct_target_candidate": candidate,
                        "cancer_antigen_candidate": cancer_candidate,
                        "reproductive_antigen_candidate": reproductive_candidate,
                        "role_reasons": reasons,
                        "decision_ids": [m["entity_id"] for m in matches],
                        "categories": [m["category"] for m in matches],
                        "classifications": [m["classification"] for m in matches],
                    }
                )
            role = max((o["shared_sequence"] for o in occurrences), key=ROLE_ORDER.get)
            segments = translation_qualification.valid_segments(sequence)
            if not segments:
                no_segments.append(accession)
            direct = role == "allow" and any(o["direct_target_candidate"] for o in occurrences)
            cancer_direct = role == "allow" and any(
                o["cancer_antigen_candidate"] for o in occurrences
            )
            reproductive_direct = role == "allow" and any(
                o["reproductive_antigen_candidate"] for o in occurrences
            )
            record = {
                "record_id": accession,
                "reference_id": reference["reference_id"],
                "sequence": sequence,
                "sequence_sha256": seqhash,
                "source_record": original,
                "qualification": quality,
                "occurrences": occurrences,
                "valid_segments": segments,
                "shared_sequence": role,
                "direct_target_candidate": direct,
                "cancer_antigen_candidate": cancer_direct,
                "reproductive_antigen_candidate": reproductive_direct,
            }
            zipped.write(
                (
                    json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False)
                    + "\n"
                ).encode()
            )
            roles[role] += 1
            occurrence_count += len(occurrences)
            target_records += direct
            cancer_antigen_records += cancer_direct
            reproductive_antigen_records += reproductive_direct

    universe = pd.read_csv(ROOT / "oncoref/data/cta-protein-universe.csv", keep_default_na=False)
    out = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/out-of-universe.csv.gz", keep_default_na=False
    )
    coverage = {
        "reference_id": reference["reference_id"],
        "scope": "all Ensembl 115 peptide FASTA records, including alternate loci and NMD/LoF/incomplete records",
        "classified_gene_loci_without_translation": sorted(set(per_gene) - seen_genes),
        "source_gene_loci_without_current_decision": sorted(seen_genes - set(per_gene)),
        "out_of_universe_groups": out.to_dict("records"),
        "nc_orf_nominations_without_sequence_in_this_reference_bundle": universe[
            universe.protein_source.eq("nc_orf")
        ].to_dict("records"),
        "records_without_aa20_segments": no_segments,
        "limitations": [
            "No empirical ncORF peptide inventory is supplied by this Ensembl peptide bundle; gene-level ncORF evidence cannot validate an unrelated reference peptide.",
            "Nonrepresentative isoforms require review unless applicable gene-level evidence excludes them.",
            "Proteins with unscreened recorded home tissues (olfactory epithelium, nasal mucosa, taste bud, hair follicle, dental pulp, bone) remain review.",
            "Alternate loci receive no primary-locus allowance through aliases.",
            "Gene-level normal evidence is conservatively applied to its annotated isoforms; it does not measure isoform-specific production.",
            "No patient-HLA presentation, approximate similarity threshold or clinical safety is established.",
        ],
    }
    (output / "coverage.json").write_text(json.dumps(coverage, indent=2) + "\n")
    shutil.copyfile(ROOT / "oncoref/data/target-selection.csv.gz", output / "decisions.csv.gz")
    if proteome_path.resolve() != (output / "orfdata-proteome.json.gz").resolve():
        shutil.copyfile(proteome_path, output / "orfdata-proteome.json.gz")

    dependencies = [
        "target-selection.json",
        "target-selection.csv.gz",
        "gene-knowledge.json",
        "gene-knowledge-updates.json",
        "rare-cell-attribution-reviews.json",
        "normal-expression-literature.json",
        "readthrough-loci.csv",
        "normal-tissue-tiers.yaml",
        "normal-tissue-tier-updates.yaml",
        "normal-tissue-taxonomy.yaml",
        "translation-qualification.json",
        "translation-qualification.csv.gz",
        "transcript-selection.json",
        "transcript-selection.csv.gz",
        "transcript-evidence.json",
    ]
    manifest = {
        "schema": SCHEMA,
        "schema_version": 1,
        "release": release,
        "reference": reference,
        "selection_policy_id": target_selection.POLICY_ID,
        "role_policy_id": ROLE_POLICY_ID,
        "translation_policy_id": translation_qualification.POLICY_ID,
        "representative_selection_policy_id": transcript_selection.POLICY_ID,
        "orfdata_identity": proteome.identity,
        "record_count": sum(roles.values()),
        "occurrence_count": occurrence_count,
        "record_roles": dict(roles),
        "direct_target_records": target_records,
        "cancer_antigen_records": cancer_antigen_records,
        "reproductive_antigen_records": reproductive_antigen_records,
        "package_dependencies": {
            n: digest((ROOT / "oncoref/data" / n).read_bytes()) for n in dependencies
        },
        "builder_sha256": digest(Path(__file__).read_bytes()),
        "payloads": {
            p.name: {"bytes": p.stat().st_size, "sha256": digest(p.read_bytes())}
            for p in sorted(output.iterdir())
            if p.name
            in ["records.jsonl.gz", "coverage.json", "decisions.csv.gz", "orfdata-proteome.json.gz"]
        },
    }
    manifest["bundle_id"] = content_id(manifest)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    checked = load(output, expected_id=manifest["bundle_id"], require_current=False)
    if checked.as_orfdata().identity != proteome.identity:
        raise ValueError("Offline ORFdata replay changed identity")
    archive = output.parent / f"oncoref-sequences-{release}.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as zipped:
        for p in sorted(output.iterdir()):
            if p.name in manifest["payloads"] or p.name == "manifest.json":
                info = zipfile.ZipInfo(p.name, date_time=(1980, 1, 1, 0, 0, 0))
                zipped.writestr(info, p.read_bytes())
    release_tag = asset_tag or f"v{release}"
    pin = {
        "schema": "oncoref.sequence-bundle-asset",
        "schema_version": 1,
        "bundle_id": manifest["bundle_id"],
        "url": f"https://github.com/pirl-unc/oncoref/releases/download/{release_tag}/{archive.name}",
        "archive_bytes": archive.stat().st_size,
        "archive_sha256": digest(archive.read_bytes()),
        "unpacked_bytes": sum(
            (output / n).stat().st_size for n in [*manifest["payloads"], "manifest.json"]
        ),
    }
    (ROOT / "oncoref/data/sequence-bundle.json").write_text(json.dumps(pin, indent=2) + "\n")
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    for name in ["manifest.json", "coverage.json"]:
        shutil.copyfile(output / name, AUDIT_DIR / name)
    print(
        json.dumps(
            {
                "bundle_id": checked.identity,
                "records": manifest["record_count"],
                "occurrences": occurrence_count,
                "roles": dict(roles),
                "direct_target_records": target_records,
                "cancer_antigen_records": cancer_antigen_records,
                "reproductive_antigen_records": reproductive_antigen_records,
                "archive": str(archive),
                "archive_bytes": archive.stat().st_size,
            }
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--release", required=True)
    p.add_argument("--asset-tag", default=None)
    a = p.parse_args()
    build(a.source, a.output, a.release, asset_tag=a.asset_tag)
