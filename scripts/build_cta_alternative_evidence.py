#!/usr/bin/env python3
"""Pin compatible alternative measurements without treating missing rows as zero.

Source adapters retain native units, original identifiers, anatomical subregions
and developmental stages. Only matched canonical loci are collapsed into protein
groups, within one source/context. Aliases never become extra expression copies.
"""

from __future__ import annotations

import json
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import normal_tissues, reference_data  # noqa: E402
from oncoref.gene_matching import (  # noqa: E402
    match_reference_genes,
    resolve_measurement_rows,
)

DOWNLOADS = ROOT / "tmp/cta-source-consolidation"
OUT = ROOT / "docs/audits/cta-selection-20261028/alternative-evidence"
GTEX_SOURCE = "gtex_tissue_rna"
GTEX_VERSION = "v11"  # This saved analysis pins its release; public loaders default to latest.
GTEX_PIN = reference_data.REFERENCE_SOURCES[GTEX_SOURCE]["pins"][GTEX_VERSION]
GTEX_URL = GTEX_PIN["url"]
PETRO_URL = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-3929/download/zip?fileType=quantification-filtered"
# Exact anatomical scope, with the native tissue name retained in source_context.
# Cultured fibroblasts and EBV-transformed lymphocytes are not normal organ bulk.
GTEX_TISSUES = {
    "Adipose - Subcutaneous": "adipose tissue",
    "Adipose - Visceral (Omentum)": "adipose tissue",
    "Adrenal Gland": "adrenal gland",
    "Artery - Aorta": "blood vessel",
    "Artery - Coronary": "blood vessel",
    "Artery - Tibial": "blood vessel",
    "Bladder": "urinary bladder",
    "Brain - Amygdala": "amygdala",
    "Brain - Anterior cingulate cortex (BA24)": "cerebral cortex",
    "Brain - Caudate (basal ganglia)": "basal ganglia",
    "Brain - Cerebellar Hemisphere": "cerebellum",
    "Brain - Cerebellum": "cerebellum",
    "Brain - Cortex": "cerebral cortex",
    "Brain - Frontal Cortex (BA9)": "cerebral cortex",
    "Brain - Hippocampus": "hippocampal formation",
    "Brain - Hypothalamus": "hypothalamus",
    "Brain - Nucleus accumbens (basal ganglia)": "basal ganglia",
    "Brain - Putamen (basal ganglia)": "basal ganglia",
    "Brain - Spinal cord (cervical c-1)": "spinal cord",
    "Brain - Substantia nigra": "midbrain",
    "Breast - Mammary Tissue": "breast",
    "Cervix - Ectocervix": "cervix",
    "Cervix - Endocervix": "cervix",
    "Colon - Sigmoid": "colon",
    "Colon - Transverse": "colon",
    "Esophagus - Gastroesophageal Junction": "esophagus",
    "Esophagus - Mucosa": "esophagus",
    "Esophagus - Muscularis": "esophagus",
    "Fallopian Tube": "fallopian tube",
    "Heart - Atrial Appendage": "heart muscle",
    "Heart - Left Ventricle": "heart muscle",
    "Kidney - Cortex": "kidney",
    "Kidney - Medulla": "kidney",
    "Liver": "liver",
    "Lung": "lung",
    "Minor Salivary Gland": "salivary gland",
    "Muscle - Skeletal": "skeletal muscle",
    "Nerve - Tibial": "peripheral nerve",
    "Ovary": "ovary",
    "Pancreas": "pancreas",
    "Pituitary": "pituitary gland",
    "Prostate": "prostate",
    "Skin - Not Sun Exposed (Suprapubic)": "skin",
    "Skin - Sun Exposed (Lower leg)": "skin",
    "Small Intestine - Terminal Ileum": "small intestine",
    "Spleen": "spleen",
    "Stomach": "stomach",
    "Testis": "testis",
    "Thyroid": "thyroid gland",
    "Uterus": "endometrium",
    "Vagina": "vagina",
    "Whole Blood": "blood",
}


def anatomy(frame: pd.DataFrame) -> pd.DataFrame:
    """Apply the common anatomical policy to explicitly mapped tissue concepts."""
    rows = []
    for tissue in frame.tissue.unique():
        # These four atlas parents are explicit brain regions, even though the
        # older shared vocabulary has no separate node for them.
        concept = (
            "brain"
            if tissue in {"pons", "white matter", "thalamus", "medulla oblongata"}
            else tissue
        )
        match = normal_tissues.resolve_normal_tissue(concept, species="human", stage="adult")
        rows.append(
            {
                "tissue": tissue,
                "category": normal_tissues.normal_tissue_category(match),
                "mapping_status": match.mapping_status,
                "concept_id": match.concept_id,
            }
        )
    result = frame.merge(pd.DataFrame(rows), on="tissue", validate="many_to_one")
    if result.category.eq("unresolved").any():
        raise ValueError(
            f"Unmapped alternative anatomy: {result[result.category.eq('unresolved')].tissue.unique()}"
        )
    return result


def canonical_matrix(
    rows: pd.DataFrame, matrix, *, source: str, version: str, sparse: bool = False
) -> tuple[pd.DataFrame, object]:
    """Resolve complete source vectors through the central identity contract."""
    identity = match_reference_genes(
        rows, id_column="gene_id", symbol_column="label", source=source, source_version=version
    )
    decisions = resolve_measurement_rows(
        identity,
        equivalent=(lambda a, b: (matrix[a] != matrix[b]).nnz == 0)
        if sparse
        else (lambda a, b: np.array_equal(matrix[a], matrix[b], equal_nan=True)),
    )
    write_csv(identity, OUT / f"{source}-identity.csv.gz")
    write_csv(decisions, OUT / f"{source}-row-decisions.csv.gz")
    accepted = decisions[decisions.row_position.notna()].reset_index(drop=True)
    return accepted, matrix[accepted.row_position.astype(int).to_numpy()]


def aggregate_bulk(
    chosen, values, native_contexts, tissue_map, *, source, version, unit, statistic
):
    """Sum protein-identical loci within native contexts before tissue maxima.

    The output retains the measured peak subregion. Member coverage is assessed
    across every source subregion, separately from the group-expression maximum.
    """
    facts = pd.read_csv(
        ROOT / "docs/audits/cta-release-1.21.0/protein-group-classifications.csv.gz",
        usecols=["member_gene_ids", "representative_gene_id", "member_count"],
    )
    member_to_rep = {
        g: r.representative_gene_id
        for r in facts.itertuples()
        for g in json.loads(r.member_gene_ids)
    }
    keep = chosen.gene_id.isin(member_to_rep).to_numpy()
    chosen, values = chosen[keep], values[keep]
    frame = pd.DataFrame(values, index=chosen.gene_id, columns=native_contexts)
    tissue_names = sorted(set(tissue_map.values()))
    complete = pd.DataFrame(
        {
            t: frame[[c for c in native_contexts if tissue_map[c] == t]].notna().all(axis=1)
            for t in tissue_names
        }
    )
    coverage = [
        {
            "gene_id": gid,
            "source": source,
            "available_tissues": json.dumps(
                [t for t, present in zip(tissue_names, row, strict=True) if present]
            ),
        }
        for gid, row in zip(frame.index, complete.to_numpy(), strict=True)
    ]
    write_csv(pd.DataFrame(coverage), OUT / f"{source}-member-coverage.csv.gz")
    frame.index = frame.index.map(member_to_rep)
    values = frame.groupby(level=0).sum(min_count=1)
    counts = frame.groupby(level=0).count()
    member_counts = facts.set_index("representative_gene_id").member_count
    pieces = []
    for tissue in sorted(set(tissue_map.values())):
        columns = [c for c in native_contexts if tissue_map[c] == tissue]
        part = values[columns]
        pieces.append(
            pd.DataFrame(
                {
                    "gene_id": part.index,
                    "tissue": tissue,
                    "source_context": part.idxmax(axis=1).fillna(""),
                    "value": part.max(axis=1),
                    "measured_members": counts[columns].min(axis=1),
                    "member_count": part.index.map(member_counts),
                    "sampled_contexts": json.dumps(columns),
                }
            ).reset_index(drop=True)
        )
    return pd.concat(pieces, ignore_index=True).assign(
        source=source, version=version, unit=unit, statistic=statistic
    )


def bulk_sources() -> pd.DataFrame:
    """Import HPA brain nTPM and direct GTEx median TPM as separate observations."""
    brain = pd.read_csv(DOWNLOADS / "rna_brain_hpa.tsv.zip", sep="\t")
    crosswalk = pd.read_csv(DOWNLOADS / "rna_brain_hpa_subregions-current.tsv.zip", sep="\t")
    # The current anatomical crosswalk is usable only if its full vocabulary
    # exactly matches the pinned v23 measurements; its newer date is recorded.
    if set(brain.Subregion) != set(crosswalk.Subregion) or crosswalk.Subregion.duplicated().any():
        raise ValueError("HPA brain subregion crosswalk does not match v23")
    crosswalk.to_csv(OUT / "hpa-brain-anatomy.csv", index=False)
    wide = brain.pivot(index=["Gene", "Gene name"], columns="Subregion", values="nTPM")
    rows = wide.index.to_frame(index=False).rename(
        columns={"Gene": "gene_id", "Gene name": "label"}
    )
    chosen, values = canonical_matrix(rows, wide.to_numpy(), source="hpa_brain", version="v23")
    b = aggregate_bulk(
        chosen,
        values,
        wide.columns,
        crosswalk.set_index("Subregion")["Brain region"].to_dict(),
        source="hpa_brain",
        version="v23",
        unit="nTPM",
        statistic="maximum of subregion means",
    )
    gtex_path = DOWNLOADS / f"gtex-{GTEX_VERSION}-median-tpm.gct.gz"
    if sha256(gtex_path) != GTEX_PIN["sha256"]:
        raise ValueError("GTEx bytes do not match the registered concrete release")
    gtex = pd.read_csv(gtex_path, sep="\t", skiprows=2)
    tissues = (
        GTEX_TISSUES
        if GTEX_VERSION == "v8"
        else {
            re.sub(r"[ ()-]+", "_", native.replace("c-1", "cMINUS1"))
            .strip("_")
            .replace("cMINUS1", "c-1"): tissue
            for native, tissue in GTEX_TISSUES.items()
        }
    )
    if not set(tissues) <= set(gtex.columns):
        raise ValueError(f"Missing native bulk contexts in GTEx {GTEX_VERSION}")
    excluded = sorted(set(gtex.columns[2:]) - set(tissues))
    (OUT / "gtex-scope.json").write_text(
        json.dumps(
            {
                "source": GTEX_SOURCE,
                "version": GTEX_VERSION,
                "included_bulk_contexts": list(tissues),
                "excluded_contexts": excluded,
                "exclusion_scope": "Cultured/transformed cells and newer microdissected compartments are outside this whole-organ bulk adapter. Native reader retains all columns.",
            },
            indent=2,
        )
        + "\n"
    )
    chosen, values = canonical_matrix(
        gtex[["Name", "Description"]].rename(columns={"Name": "gene_id", "Description": "label"}),
        gtex[list(tissues)].to_numpy(),
        source=GTEX_SOURCE,
        version=GTEX_VERSION,
    )
    g = aggregate_bulk(
        chosen,
        values,
        list(tissues),
        tissues,
        source=GTEX_SOURCE,
        version=GTEX_VERSION,
        unit="TPM",
        statistic="maximum of tissue donor medians",
    )
    d = anatomy(pd.concat([b, g], ignore_index=True))
    return d.assign(taxonomy_id=9606, stage="adult", cell_type="bulk", negative_ihc_grade="Missing")


def embryonic_source() -> pd.DataFrame:
    """Extract stage-by-lineage TPM from 1,519 measured preimplantation cells."""
    with zipfile.ZipFile(DOWNLOADS / "petropoulos-tpm.zip") as archive:
        name = "E-MTAB-3929.expression_tpm.mtx"
        rows = pd.read_csv(
            archive.open(name + "_rows"), sep="\t", header=None, names=["gene_id", "label"]
        )
        cells = pd.read_csv(archive.open(name + "_cols"), sep="\t", header=None)[0]
        matrix = mmread(archive.open(name), spmatrix=True).tocsr()
    chosen, matrix = canonical_matrix(
        rows,
        matrix,
        source="petropoulos_E-MTAB-3929",
        version="Ensembl106_filtered_TPM",
        sparse=True,
    )
    meta = pd.read_csv(DOWNLOADS / "petropoulos-design.tsv", sep="\t").set_index("Assay").loc[cells]
    if len(meta) != 1519 or len(set(cells)) != len(cells):
        raise ValueError("Unexpected or duplicate measured embryo cells")
    meta["context"] = (
        meta["Sample Characteristic[age]"]
        + " / "
        + meta["Factor Value[inferred lineage]"].fillna("unspecified lineage")
    )
    rows = []
    for lineage, indices in meta.groupby("context").indices.items():
        values = matrix[:, indices]
        counts = np.asarray((values >= 1).sum(axis=1)).ravel()
        rows.append(
            pd.DataFrame(
                {
                    "Gene ID": chosen.gene_id,
                    "Lineage": lineage,
                    "Cells": len(indices),
                    "Positive cells TPM>=1": counts,
                    "Fraction TPM>=1": counts / len(indices),
                    "Mean TPM": np.asarray(values.mean(axis=1)).ravel(),
                    "Embryos": meta.iloc[indices]["Sample Characteristic[individual]"].nunique(),
                    "source": "petropoulos_E-MTAB-3929",
                    "stage": "preimplantation, 3–7 days",
                    "unit": "TPM",
                }
            )
        )
    meta.reset_index().to_csv(
        OUT / "petropoulos-measured-cell-metadata.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    return pd.concat(rows, ignore_index=True)


def build() -> None:
    """Generate auditable, source-separated alternative-evidence artifacts."""
    OUT.mkdir(parents=True, exist_ok=True)
    selections = [
        ("hpa_brain", "v23", "rna_brain_hpa.tsv.zip"),
        ("hpa_brain_anatomy", "v25.1-20261028", "rna_brain_hpa_subregions-current.tsv.zip"),
        (GTEX_SOURCE, GTEX_VERSION, f"gtex-{GTEX_VERSION}-median-tpm.gct.gz"),
        ("petropoulos_embryonic_tpm", "E-MTAB-3929-20261028", "petropoulos-tpm.zip"),
        ("petropoulos_embryonic_design", "E-MTAB-3929-20261028", "petropoulos-design.tsv"),
    ]
    sources = []
    for source, version, file in selections:
        concrete = reference_data.resolve_version(source, version)
        pin = reference_data.REFERENCE_SOURCES[source]["pins"][concrete]
        if sha256(DOWNLOADS / file) != pin["sha256"]:
            raise ValueError(f"{source}@{concrete} input checksum mismatch")
        sources.append(
            {
                "name": source,
                "version": concrete,
                "file": file,
                "url": pin["url"],
                "sha256": pin["sha256"],
            }
        )
    write_csv(bulk_sources(), OUT / "normal-bulk.csv.gz")
    write_csv(embryonic_source(), OUT / "developmental-lineages.csv.gz")
    (OUT / "manifest.json").write_text(
        json.dumps(
            {
                "inputs": input_hashes(
                    [
                        Path(__file__),
                        ROOT / "oncoref/gene_matching.py",
                        ROOT / "oncoref/data/reference-releases.json",
                    ]
                ),
                "sources": sources,
                "outputs": {
                    p.name: sha256(p) for p in sorted(OUT.iterdir()) if p.name != "manifest.json"
                },
                "limitations": [
                    "GTEx donor medians do not exclude high expression in individual donors.",
                    "HPA anatomical crosswalk downloaded 2026-10-28; all 193 labels verified against v23.",
                    "Measured preimplantation stages complement one gastrulating embryo; they do not cover all fetal or pediatric stages.",
                ],
            },
            indent=2,
        )
        + "\n"
    )
    print("Alternative evidence imported", flush=True)


if __name__ == "__main__":
    build()
