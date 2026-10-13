#!/usr/bin/env python3
"""Replay the October 2026 canine intake from acquired primary-source files.

Pass --source-root to the existing acquisition snapshot (data/ and docs/ receipts).
No consumer package is imported. Matrix bytes remain outside the wheel. Extracted
sample facts retain original labels, primary-source URLs and exact byte digests.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import yaml

from oncoref.expression_contract import ExpressionContract, ExpressionMeasurement


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def build(root, output):
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((output / "assets.json").read_text())
    manifest["assets"] = [a for a in manifest["assets"] if a["group"] != "expression"]
    sources = []
    overrides = {}
    quantitative = json.loads((root / "docs/bgee-quantitative-assets.json").read_text())["assets"]
    records = []
    for q in quantitative:
        records.append(("bgee-" + q["study"].lower(), q, "TPM", "Ensembl gene", "bgee", q["study"]))
    bark = next(
        r
        for r in json.loads((root / "docs/source-registry.json").read_text())["assets"]
        if r["asset_id"] == "barkbase-tmm-gene.tsv"
    )
    records.append(
        ("barkbase", bark, "TMM-CPM", "BarkBase StringTie gene", "barkbase", "PRJNA396033")
    )
    audit = []
    for source_id, receipt, unit, namespace, adapter, study in records:
        path = root / "data/processed" / receipt["asset_id"]
        if digest(path) != receipt["sha256"] or path.stat().st_size != receipt["size_bytes"]:
            raise ValueError(f"Source receipt mismatch: {path}")
        sample_rows = []
        if adapter == "bgee":
            condition_columns = [
                "Library ID",
                "Anatomical entity ID",
                "Anatomical entity name",
                "Stage ID",
                "Stage name",
                "Sex",
                "Strain",
            ]
            chunks = pd.read_csv(
                path,
                sep="\t",
                usecols=condition_columns,
                dtype=str,
                keep_default_na=False,
                chunksize=100000,
            )
            conditions = pd.concat([chunk.drop_duplicates() for chunk in chunks]).drop_duplicates()
            if conditions["Library ID"].duplicated().any():
                raise ValueError("Conflicting library conditions")
            runs_path = root / f"data/report-reconciliation/{study.lower()}-runs.tsv"
            runs = pd.read_csv(runs_path, sep="\t", dtype=str, keep_default_na=False)
            runs_digest = digest(runs_path)
            runs_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={study}&result=read_run&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,sample_alias,sample_title,library_name,library_layout,read_count,base_count&format=tsv"
            for row in conditions.sort_values("Library ID").to_dict("records"):
                library = row["Library ID"]
                matching = runs[runs.experiment_accession == library]
                if matching.empty or matching.sample_accession.nunique() != 1:
                    raise ValueError(f"Unresolved library-to-sample mapping {library}")
                sample_id = matching.iloc[0].sample_accession
                xml = root / "data/report-reconciliation" / f"{sample_id}.xml"
                tree = ET.parse(xml)
                if tree.findtext(".//TAXON_ID") != "9615":
                    raise ValueError(f"Noncanine BioSample {sample_id}")
                attrs = {
                    item.findtext("TAG"): item.findtext("VALUE")
                    for item in tree.findall(".//SAMPLE_ATTRIBUTE")
                }
                tissue = row["Anatomical entity name"]
                original = tissue
                source_tissue = attrs.get("source_name", attrs.get("isolation source"))
                corrected = study == "SRP181218" and source_tissue == "Heart" and tissue == "colon"
                if corrected:
                    tissue = "heart"
                individual = attrs.get("individual")
                sample_rows.append(
                    {
                        "sample_id": library,
                        "taxon_id": 9615,
                        "study": study,
                        "library_id": library,
                        "run_ids": sorted(matching.run_accession.unique()),
                        "biosample_id": sample_id,
                        "source_donor_label": individual,
                        "donor_id": f"{study}:individual:{individual}" if individual else None,
                        "donor_scope": study,
                        "pooling": "single" if individual else "unknown",
                        "tissue": tissue,
                        "original_tissue": original,
                        "source_tissue": source_tissue,
                        "tissue_correction": "BioSample source_name identifies Heart; Bgee labels colon"
                        if corrected
                        else None,
                        "age": attrs.get("age"),
                        "sex": row["Sex"],
                        "stage": row["Stage name"],
                        "health": attrs.get("health status"),
                        "treatment": attrs.get("treatment"),
                        "source_attributes": attrs,
                        "source_condition": row,
                        "evidence": [
                            {
                                "url": f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_id}",
                                "sha256": digest(xml),
                            },
                            {"url": runs_url, "sha256": runs_digest},
                        ],
                    }
                )
        else:
            with path.open() as stream:
                columns = next(csv.reader(stream, delimiter="\t"))[3:]
            paper = root / "data/report-reconciliation/barkbase-paper.xml"
            for column in columns:
                match = re.fullmatch(r"(Adult|Emb)_(\d+)(?:_day_(\d+))?(?:_(.+))?", column)
                if not match:
                    raise ValueError(f"Unexpected BarkBase column {column}")
                stage, individual, day, tissue = match.groups()
                donor_label = f"{stage}_{individual}"
                sample_rows.append(
                    {
                        "sample_id": column,
                        "taxon_id": 9615,
                        "study": study,
                        "library_id": column,
                        "run_ids": [],
                        "biosample_id": None,
                        "source_donor_label": donor_label,
                        # Only the adult donor-label interpretation is established by the paper.
                        "donor_id": f"{study}:{donor_label}" if stage == "Adult" else None,
                        "donor_scope": study,
                        "pooling": "single" if stage == "Adult" else "unknown",
                        "tissue": tissue,
                        "original_tissue": tissue,
                        "source_tissue": tissue,
                        "tissue_correction": None,
                        "age": f"gestational day {day}" if day else None,
                        "sex": None,
                        "stage": "adult" if stage == "Adult" else "embryonic",
                        "health": None,
                        "treatment": None,
                        "source_attributes": {"matrix_column": column},
                        "source_condition": {
                            "donor_label": donor_label,
                            "run_mapping": "unresolved",
                        },
                        "evidence": [
                            {
                                "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6627511/fullTextXML",
                                "sha256": digest(paper),
                            },
                            {"url": receipt["url"], "sha256": receipt["sha256"]},
                        ],
                    }
                )
        sample_file = output / f"{source_id}-samples.json"
        write_json(
            sample_file,
            {
                "schema": "oncoref.expression-samples",
                "schema_version": 1,
                "taxon_id": 9615,
                "source_id": source_id,
                "samples": sample_rows,
            },
        )
        version = receipt["source_version"]
        contract = ExpressionContract(
            taxon_id=9615,
            source_id=source_id,
            source_version=version,
            context="normal",
            reference=None,
            reference_status="unresolved",
            feature_namespace=namespace,
            measurement=ExpressionMeasurement("gene", unit, "linear", False),
            content_provenance={"matrix": receipt["sha256"], "samples": digest(sample_file)},
            provenance_notes="Native source feature space. Exact quantification assembly/annotation remains unresolved; available canine GTFs and symbol/orthology matches do not establish equivalence. Original sample metadata and corrections are preserved. No isoform allocation or target admission.",
        )
        asset_id = source_id + "-matrix"
        manifest["assets"].append(
            {
                "id": asset_id,
                "taxon_id": 9615,
                "group": "expression",
                "filename": path.name,
                "url": receipt["url"],
                "sha256": receipt["sha256"],
                "size_bytes": receipt["size_bytes"],
                "source_version": version,
                "description": f"{source_id} native gene {unit} matrix",
                "license": receipt["license"],
            }
        )
        overrides["9615:" + asset_id] = str(path.resolve())
        sources.append(
            {
                "id": source_id,
                "taxon_id": 9615,
                "category": "expression",
                "context": "normal",
                "cancer_codes": [],
                "source_type": adapter,
                "accession": study,
                "source_version": version,
                "unit": unit,
                "source_scale_class": "linear",
                "tpm_proxy": False,
                "linear_tpm_comparable": unit == "TPM",
                "contract": contract.as_dict(),
                "matrix_asset": asset_id,
                "sample_metadata": sample_file.name,
                "adapter": adapter,
            }
        )
        audit.append(
            {
                "source": source_id,
                "libraries": len(sample_rows),
                "known_donors": len({r["donor_id"] for r in sample_rows if r["donor_id"]}),
                "tissue_corrections": sum(r["tissue_correction"] is not None for r in sample_rows),
            }
        )
    # This is a source-inventory overlap, not overlap of the two imported Bgee matrices.
    library_file = root / "data/inventory-2026-10-08/Canis_lupus_familiaris_RNA-Seq_libraries.tsv"
    ena_file = root / "data/inventory-2026-10-08/barkbase-ena"
    with library_file.open() as stream:
        libraries = list(csv.DictReader(stream, delimiter="\t"))
    with ena_file.open() as stream:
        ena = list(csv.DictReader(stream, delimiter="\t"))
    bgee_runs = {run.strip() for row in libraries for run in row["Run IDs"].split(",")}
    bark_runs = {
        r["run_accession"]
        for r in ena
        if r["tax_id"] == "9615" and r["library_strategy"] == "RNA-Seq"
    }
    write_json(
        output / "lineage_overlaps.json",
        {
            "taxon_id": 9615,
            "scope": "Bgee global library inventory and BarkBase ENA RNA-seq inventory; does not assert overlap of the two imported Bgee matrices or map BarkBase matrix columns to runs",
            "source_a": "Bgee library inventory retrieved 2026-10-08",
            "source_b": "PRJNA396033",
            "shared_run_ids": sorted(bgee_runs & bark_runs),
            "independent_donors": None,
            "content_provenance": {
                "bgee_library_table": digest(library_file),
                "barkbase_ena": digest(ena_file),
            },
        },
    )
    write_json(output / "assets.json", manifest)
    (output / "expression_sources.yaml").write_text(
        yaml.safe_dump({"sources": sources}, sort_keys=False, allow_unicode=True)
    )
    return {"audit": audit, "local_overrides": overrides}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("oncoref/data/species/dog"))
    args = parser.parse_args()
    print(json.dumps(build(args.source_root, args.output), indent=2))
