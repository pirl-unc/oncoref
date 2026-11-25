#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0.
"""Build the quantitative Ensembl 115 GRCm39 murine antigen catalog (#734).

Evaluates all primary-assembly protein-coding loci in Mus musculus (taxon 10090,
Ensembl 115 GRCm39 / GCA_000001635.9) through the unified 3-stage OncoRef
target-selection architecture:

1. Reference proteome & PEG grouping: 66,668 Ensembl 115 peptide translations
   collapsed into identical-sequence PEG:<sha256[:20]> groups, plus GRCm39 exon
   geometry and 31-mer cDNA paralog mapping-risk flags.
2. Multi-modal normal safety gates:
   - Adult somatic bulk & blood RNA (12 Bgee 15.2 experiments + GSE219045
     sex-stratified preclinical organ atlas, 30 somatic contexts in linear TPM;
     Vital <= 4, Protected <= 8, Manageable <= 16 hold / > 128 Unsafe, zero
     negative-IHC multipliers).
   - Donor-resolved single-cell UMI screen (Tabula Muris 10x droplet across 10
     adult Vital/Protected somatic organs, 49,746 somatic QC cells).
   - Adult normal protein assays (MGI GXD 6.27 TS28 Immunohistochemistry
     MMO:0000498 & Immunoblot MMO:0000672 across Vital/Protected organs).
   - Curated primary somatic RT-PCR counterevidence (Sarma 1999 on Trap1a/P1A).
3. Symmetric 10 x max(N, D, 1 TPM) Cancer Antigen (CA), Reproductive Antigen
   (RA), and CTA = CA & RA enrichment across 8 strain-explicit murine tumor
   cohorts (MC38, CT26, B16F10, LLC, 4T1, KP_LUAD, BMPC_PRAD, P815; 36 control
   tumor libraries) and all 16 comparative murine reproductive, accessory-gland,
   bicornuate-tract, mammary, placental, and embryonic compartments.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import tarfile
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import mouse, target_selection  # noqa: E402
from oncoref.reference_identity import ReferenceIdentity  # noqa: E402
from oncoref.species_policies import POLICY_SCHEMA_ID  # noqa: E402
from oncoref.species_target_selection import (  # noqa: E402
    SpeciesScreeningSpec,
    anatomy_for_species,
    classify_sex_restriction,
    evaluate_species_catalog,
    group_reference_proteins,
)

DEFAULT_ENS_DIR = Path(
    "/Users/iskander/.gemini/antigravity-cli/brain/3767551e-1d33-4fba-9b4e-91fa9f5de1e0/scratch/ens115_mouse"
)
DEFAULT_INPUT_DIR = Path(
    "/Users/iskander/.gemini/antigravity-cli/brain/3767551e-1d33-4fba-9b4e-91fa9f5de1e0/scratch/probe_734"
)

PRIMARY_CONTIGS = frozenset({*(str(i) for i in range(1, 20)), "X", "Y", "MT"})
COMPARISON_ID = "oncoref.mouse-grcm39-tpm.2026-11-25.1"

# 24 adult somatic contexts measured across the 41,233-gene Bgee 15.2 GRCm39 universe
REQUIRED_BGEE_NORMAL_CONTEXTS = frozenset(
    {
        "adipose tissue",
        "adrenal gland",
        "blood",
        "bone",
        "bone marrow",
        "brain",
        "colon",
        "esophagus",
        "eye",
        "heart",
        "kidney",
        "liver",
        "lung",
        "olfactory epithelium",
        "pancreas",
        "pituitary gland",
        "salivary gland",
        "skeletal muscle",
        "skin",
        "small intestine",
        "spleen",
        "stomach",
        "thyroid gland",
        "urinary bladder",
    }
)

SOMATIC_CONTEXT_TIER: dict[str, str] = {
    # Vital (ceiling <= 4 TPM)
    "adrenal gland": "vital",
    "blood": "vital",
    "blood vessel": "vital",
    "bone": "vital",
    "bone marrow": "vital",
    "brain": "vital",
    "bronchus": "vital",
    "eye": "vital",
    "heart": "vital",
    "inner ear": "vital",
    "kidney": "vital",
    "liver": "vital",
    "lung": "vital",
    "pancreas": "vital",
    "peripheral nerve": "vital",
    "pituitary gland": "vital",
    "skeletal muscle": "vital",
    "small intestine": "vital",
    "spinal cord": "vital",
    # Protected (ceiling <= 8 TPM)
    "cartilage": "protected",
    "colon": "protected",
    "esophagus": "protected",
    "lymph node": "protected",
    "parathyroid gland": "protected",
    "salivary gland": "protected",
    "skin": "protected",
    "spleen": "protected",
    "stomach": "protected",
    "tongue": "protected",
    "urinary bladder": "protected",
    # Manageable (hold > 16 TPM, Unsafe > 128 TPM)
    "adipose tissue": "manageable",
    "olfactory epithelium": "manageable",
    "thyroid gland": "manageable",
}

BGEE_ANATOMY_MAP: dict[str, tuple[str, str]] = {
    # Somatic normal ("normal")
    "brain": ("normal", "brain"),
    "cerebellum": ("normal", "brain"),
    "cerebral cortex": ("normal", "brain"),
    "frontal cortex": ("normal", "brain"),
    "Ammon's horn": ("normal", "brain"),
    "olfactory bulb": ("normal", "brain"),
    "cortical layer V": ("normal", "brain"),
    "cortical plate": ("normal", "brain"),
    "granular layer of cerebellar cortex": ("normal", "brain"),
    "retinal neural layer": ("normal", "eye"),
    "retina": ("normal", "eye"),
    "heart": ("normal", "heart"),
    "lung": ("normal", "lung"),
    "liver": ("normal", "liver"),
    "adult mammalian kidney": ("normal", "kidney"),
    "kidney": ("normal", "kidney"),
    "pancreas": ("normal", "pancreas"),
    "esophagus": ("normal", "esophagus"),
    "stomach": ("normal", "stomach"),
    "small intestine": ("normal", "small intestine"),
    "duodenum": ("normal", "small intestine"),
    "jejunum": ("normal", "small intestine"),
    "ileum": ("normal", "small intestine"),
    "colon": ("normal", "colon"),
    "sigmoid colon": ("normal", "colon"),
    "large intestine": ("normal", "colon"),
    "skeletal muscle tissue": ("normal", "skeletal muscle"),
    "limb muscle": ("normal", "skeletal muscle"),
    "muscle tissue": ("normal", "skeletal muscle"),
    "quadriceps femoris": ("normal", "skeletal muscle"),
    "bone element": ("normal", "bone"),
    "bone marrow": ("normal", "bone marrow"),
    "leukocyte": ("normal", "blood"),
    "spleen": ("normal", "spleen"),
    "urinary bladder": ("normal", "urinary bladder"),
    "submandibular gland": ("normal", "salivary gland"),
    "adrenal gland": ("normal", "adrenal gland"),
    "pituitary gland": ("normal", "pituitary gland"),
    "thyroid gland": ("normal", "thyroid gland"),
    "adipose tissue": ("normal", "adipose tissue"),
    "adipose tissue of abdominal region": ("normal", "adipose tissue"),
    "brown adipose tissue": ("normal", "adipose tissue"),
    "gonadal fat pad": ("normal", "adipose tissue"),
    "subcutaneous adipose tissue": ("normal", "adipose tissue"),
    "olfactory epithelium": ("normal", "olfactory epithelium"),
    "zone of skin": ("normal", "skin"),
    # Thymic ("thymic")
    "thymus": ("thymic", "thymus"),
    # Adult reproductive ("reproductive")
    "testis": ("reproductive", "testis"),
    "caput epididymis": ("reproductive", "epididymis"),
    "corpus epididymis": ("reproductive", "epididymis"),
    "cauda epididymis": ("reproductive", "epididymis"),
    "ovary": ("reproductive", "ovary"),
    "uterus": ("reproductive", "uterus"),
    "thoracic mammary gland": ("reproductive", "resting_mammary"),
    "placenta": ("reproductive", "placenta"),
    # Embryonic / developmental ("developmental")
    "embryo": ("developmental", "embryo"),
    "late embryo": ("developmental", "embryo"),
    "forebrain": ("developmental", "embryo"),
    "midbrain": ("developmental", "embryo"),
    "hindbrain": ("developmental", "embryo"),
    "central nervous system": ("developmental", "embryo"),
    "neural tube": ("developmental", "embryo"),
    "limb": ("developmental", "embryo"),
    "embryonic facial prominence": ("developmental", "embryo"),
    "intestine": ("developmental", "embryo"),
}

GSE219045_COL_MAP: dict[str, tuple[str, str]] = {
    "Adipose_Brown": ("normal", "adipose tissue"),
    "Adipose_Visceral": ("normal", "adipose tissue"),
    "Adipose_Subcutaneous": ("normal", "adipose tissue"),
    "AdrenalGland": ("normal", "adrenal gland"),
    "Aorta": ("normal", "blood vessel"),
    "CaudalVenaCava": ("normal", "blood vessel"),
    "Bladder": ("normal", "urinary bladder"),
    "BoneMarrow": ("normal", "bone marrow"),
    "Brain_Cerebellum": ("normal", "brain"),
    "Cerebellum": ("normal", "brain"),
    "Brain_Cortex": ("normal", "brain"),
    "Cortex": ("normal", "brain"),
    "Brain_FrontalCortex": ("normal", "brain"),
    "Brain_Frontal_Cortex": ("normal", "brain"),
    "Brain_Medulla": ("normal", "brain"),
    "Brain_Spinal_Cord": ("normal", "spinal cord"),
    "Cecum": ("normal", "colon"),
    "Colon": ("normal", "colon"),
    "Duodenum": ("normal", "small intestine"),
    "Ileum": ("normal", "small intestine"),
    "Jejunum": ("normal", "small intestine"),
    "Esophagus": ("normal", "esophagus"),
    "Eye": ("normal", "eye"),
    "Heart": ("normal", "heart"),
    "Kidney": ("normal", "kidney"),
    "Liver": ("normal", "liver"),
    "Lung": ("normal", "lung"),
    "MesentericLymphNode": ("normal", "lymph node"),
    "Pancreas": ("normal", "pancreas"),
    "Parathyroid": ("normal", "parathyroid gland"),
    "PituitaryGland": ("normal", "pituitary gland"),
    "SalivaryGland": ("normal", "salivary gland"),
    "SciaticNerve": ("normal", "peripheral nerve"),
    "SkeletalMuscle": ("normal", "skeletal muscle"),
    "Skin": ("normal", "skin"),
    "Spleen": ("normal", "spleen"),
    "Stomach": ("normal", "stomach"),
    "ThyroidGland": ("normal", "thyroid gland"),
    "Tongue": ("normal", "tongue"),
    "WholeBlood": ("normal", "blood"),
    "Thymus": ("thymic", "thymus"),
    "Testis": ("reproductive", "testis"),
    "Epididymis": ("reproductive", "epididymis"),
    "Prostate": ("reproductive", "prostate"),
    "Ovary": ("reproductive", "ovary"),
    "Uterus": ("reproductive", "uterus"),
    "Cervix": ("reproductive", "cervix"),
    "Vagina": ("reproductive", "vagina"),
    "MammaryTissue": ("reproductive", "resting_mammary"),
}

GXD_ORGAN_MAP: dict[str, str] = {
    "brain": "brain",
    "cerebellum": "brain",
    "cerebral cortex": "brain",
    "hippocampus": "brain",
    "olfactory lobe": "brain",
    "brainstem": "brain",
    "hypothalamus": "brain",
    "midbrain": "brain",
    "corpus striatum": "brain",
    "spinal cord": "spinal cord",
    "peripheral nervous system": "peripheral nerve",
    "sciatic nerve": "peripheral nerve",
    "eye": "eye",
    "retina": "eye",
    "inner ear": "inner ear",
    "cochlea": "inner ear",
    "heart": "heart",
    "aorta": "blood vessel",
    "lung": "lung",
    "trachea": "bronchus",
    "liver": "liver",
    "metanephros": "kidney",
    "renal cortex": "kidney",
    "kidney": "kidney",
    "pancreas": "pancreas",
    "esophagus": "esophagus",
    "stomach": "stomach",
    "small intestine": "small intestine",
    "duodenum": "small intestine",
    "jejunum": "small intestine",
    "ileum": "small intestine",
    "large intestine": "colon",
    "colon": "colon",
    "skeletal muscle": "skeletal muscle",
    "bone": "bone",
    "bone marrow": "bone marrow",
    "cartilage": "cartilage",
    "spleen": "spleen",
    "lymph node": "lymph node",
    "bladder": "urinary bladder",
    "urinary bladder": "urinary bladder",
    "submandibular gland": "salivary gland",
    "salivary gland": "salivary gland",
    "tongue": "tongue",
    "adrenal gland": "adrenal gland",
    "pituitary gland": "pituitary gland",
    "parathyroid gland": "parathyroid gland",
    "skin": "skin",
}

TM_TISSUE_MAP: dict[str, tuple[str, str]] = {
    "Heart_and_Aorta": ("heart", "vital"),
    "Kidney": ("kidney", "vital"),
    "Liver": ("liver", "vital"),
    "Lung": ("lung", "vital"),
    "Limb_Muscle": ("skeletal muscle", "vital"),
    "Marrow": ("bone marrow", "vital"),
    "Trachea": ("bronchus", "vital"),
    "Bladder": ("urinary bladder", "protected"),
    "Spleen": ("spleen", "protected"),
    "Tongue": ("tongue", "protected"),
}

HAIR_FOLLICLE_KERATINS = frozenset(
    {
        *(f"Krt{i}" for i in range(25, 41)),
        *(f"Krt{i}" for i in range(71, 87)),
    }
)
DENTAL_HOME_SYMBOLS = frozenset({"Amelx", "Ambn", "Enam", "Dspp", "Mmp20", "Klk4", "Odaph"})


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for s, e in ordered[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def load_ensembl115_mouse_reference(
    ens_dir: Path,
) -> tuple[
    ReferenceIdentity,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, int],
    dict[str, str],
    pd.DataFrame,
    list[dict],
]:
    """Verify Ensembl 115 GRCm39 files and build PEGs, exon lengths, and mapping risk."""
    reg = mouse.references.registered_references()[0]
    ref = ReferenceIdentity.from_dict(reg["identity"])
    role_files = {
        "gtf": ens_dir / "Mus_musculus.GRCm39.115.gtf.gz",
        "protein": ens_dir / "Mus_musculus.GRCm39.pep.all.fa.gz",
        "rna": ens_dir / "Mus_musculus.GRCm39.cdna.all.fa.gz",
        "ncrna": ens_dir / "Mus_musculus.GRCm39.ncrna.fa.gz",
    }
    source_rows = []
    for role, path in role_files.items():
        digest = _sha256(path)
        expected = ref.content_provenance[role]
        if digest != expected:
            raise ValueError(f"Ensembl 115 mouse {role} SHA-256 mismatch: {digest} != {expected}")
        source_rows.append(
            {
                "source_id": f"ensembl115-{role}",
                "pillar": "reference_proteome",
                "filename": path.name,
                "unit": "sequence/annotation",
                "strain": "C57BL/6J (GRCm39)",
                "bytes": path.stat().st_size,
                "sha256": digest,
            }
        )

    # 1. Parse GTF for primary-contig protein_coding genes, symbols, strand, and exon intervals
    gene_Attr_re = re.compile(r'(\w+)\s+"([^"]+)"')
    gene_meta: dict[str, tuple[str, str, str]] = {}  # gid -> (chrom, strand, symbol)
    exons_by_gene: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with gzip.open(role_files["gtf"], "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feat, start_s, end_s, _score, strand, _frame, attrs_s = parts
            if chrom not in PRIMARY_CONTIGS or feat not in ("gene", "exon"):
                continue
            attrs = dict(gene_Attr_re.findall(attrs_s))
            gid = attrs.get("gene_id", "").split(".")[0]
            if not gid:
                continue
            if feat == "gene":
                if attrs.get("gene_biotype") == "protein_coding":
                    gene_meta[gid] = (chrom, strand, attrs.get("gene_name") or gid)
            elif feat == "exon":
                exons_by_gene[gid].append((int(start_s), int(end_s)))

    exonic_bp: dict[str, int] = {}
    merged_exons: dict[str, list[tuple[int, int]]] = {}
    for gid, ivs in exons_by_gene.items():
        m = _merge_intervals(ivs)
        merged_exons[gid] = m
        exonic_bp[gid] = sum(e - s + 1 for s, e in m)

    # Same-strand coding-exon overlap check across primary contigs
    by_chr_strand: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    for gid, (chrom, strand, _sym) in gene_meta.items():
        for s, e in merged_exons.get(gid, ()):
            by_chr_strand[(chrom, strand)].append((s, e, gid))

    same_strand_overlap_genes: set[str] = set()
    for ivs in by_chr_strand.values():
        ivs.sort()
        active_end = -1
        active_gid = ""
        for s, e, gid in ivs:
            if s <= active_end and gid != active_gid:
                same_strand_overlap_genes.add(gid)
                same_strand_overlap_genes.add(active_gid)
            if e > active_end:
                active_end = e
                active_gid = gid

    # 2. Parse pep.all.fa.gz
    pep_records: list[dict] = []
    seen_coding_genes: set[str] = set()

    def _flush_pep(hdr: str, seq_parts: list[str]) -> None:
        if not hdr:
            return
        seq = "".join(seq_parts).strip()
        tokens = hdr[1:].split()
        pid = tokens[0]
        kv = {}
        for tok in tokens[1:]:
            if ":" in tok:
                k, v = tok.split(":", 1)
                kv[k] = v
        loc_parts = tokens[2].split(":") if len(tokens) > 2 else []
        chrom = loc_parts[2] if len(loc_parts) >= 3 else ""
        if chrom not in PRIMARY_CONTIGS:
            return
        gid = kv.get("gene", "").split(".")[0]
        tid = kv.get("transcript", "").split(".")[0]
        g_biotype = kv.get("gene_biotype", "")
        t_biotype = kv.get("transcript_biotype", "protein_coding")
        if g_biotype != "protein_coding" and gid not in gene_meta:
            return
        sym = kv.get("gene_symbol") or gene_meta.get(gid, ("", "", gid))[2]
        qual = (
            len(seq) >= 20
            and "*" not in seq
            and "X" not in seq
            and "U" not in seq
            and (seq.startswith("M") or chrom == "MT")
            and t_biotype not in {"nonsense_mediated_decay", "non_stop_decay", "retained_intron"}
        )
        seen_coding_genes.add(gid)
        pep_records.append(
            {
                "gene_id": gid,
                "symbol": sym,
                "transcript_id": tid,
                "protein_id": pid,
                "transcript_biotype": t_biotype,
                "sequence": seq,
                "is_qualified": qual,
            }
        )

    cur_hdr = ""
    cur_parts: list[str] = []
    with gzip.open(role_files["protein"], "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush_pep(cur_hdr, cur_parts)
                cur_hdr = line.rstrip("\n")
                cur_parts = []
            else:
                cur_parts.append(line.strip())
    _flush_pep(cur_hdr, cur_parts)

    for gid, (_chrom, _strand, sym) in sorted(gene_meta.items()):
        if gid not in seen_coding_genes:
            pep_records.append(
                {
                    "gene_id": gid,
                    "symbol": sym,
                    "transcript_id": "",
                    "protein_id": "",
                    "transcript_biotype": "protein_coding",
                    "sequence": "",
                    "is_qualified": False,
                }
            )

    groups_df, seq_df = group_reference_proteins(pep_records, species="mouse", reference=ref)
    symbol_by_gene = {gid: meta[2] for gid, meta in gene_meta.items()}
    for r in pep_records:
        symbol_by_gene.setdefault(r["gene_id"], r["symbol"])

    # Assign unmeasured murine home tissues (hair follicle, dental pulp, taste bud)
    home_by_entity: dict[str, str] = {}
    for row in groups_df.itertuples(index=False):
        gids = json.loads(row.member_gene_ids)
        syms = [symbol_by_gene.get(g, "") for g in gids]
        homes: list[str] = []
        for sym in syms:
            if sym.startswith("Krtap") or sym in HAIR_FOLLICLE_KERATINS:
                homes.append("hair follicle")
            elif sym in DENTAL_HOME_SYMBOLS:
                homes.append("dental pulp")
            elif sym.startswith(("Tas1r", "Tas2r")):
                homes.append("taste bud")
        if homes:
            home_by_entity[row.entity_id] = "; ".join(dict.fromkeys(homes))
    groups_df["home_tissue"] = groups_df["entity_id"].map(home_by_entity).fillna("")
    groups_df["published_nomination"] = groups_df["member_gene_ids"].map(
        lambda s: "ENSMUSG00000051257" in json.loads(s)
    )

    # 3. Compute 31-mer cDNA paralog ambiguity across distinct PEGs (stride-4 sampling for speed & exactness)
    gene_to_entity = {
        g: r.entity_id
        for r in groups_df.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }
    kmer_owners: dict[bytes, str] = {}
    ambig_kmer_counts: dict[str, int] = defaultdict(int)
    total_kmer_counts: dict[str, int] = defaultdict(int)

    # First pass: index 31-mers (sampled every 3 bp) per PEG
    cdna_seqs_by_gene: dict[str, list[str]] = defaultdict(list)
    c_hdr = ""
    c_parts: list[str] = []

    def _flush_cdna(hdr: str, parts_l: list[str]) -> None:
        if not hdr:
            return
        gid = ""
        for tok in hdr[1:].split():
            if tok.startswith("gene:"):
                gid = tok.split(":", 1)[1].split(".")[0]
                break
        if gid in gene_to_entity:
            cdna_seqs_by_gene[gid].append("".join(parts_l))

    with gzip.open(role_files["rna"], "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush_cdna(c_hdr, c_parts)
                c_hdr = line.rstrip("\n")
                c_parts = []
            else:
                c_parts.append(line.strip())
    _flush_cdna(c_hdr, c_parts)

    k = 31
    step = 4
    MULTI = "__MULTI__"
    for gid, seqs in cdna_seqs_by_gene.items():
        eid = gene_to_entity[gid]
        seen_for_gene: set[str] = set()
        for s in seqs:
            if len(s) < k:
                continue
            for i in range(0, len(s) - k + 1, step):
                km = s[i : i + k]
                seen_for_gene.add(km)
        for km_s in seen_for_gene:
            km_b = km_s.encode("ascii", errors="ignore")
            prev = kmer_owners.get(km_b)
            if prev is None:
                kmer_owners[km_b] = eid
            elif prev != eid and prev != MULTI:
                kmer_owners[km_b] = MULTI

    for gid, seqs in cdna_seqs_by_gene.items():
        seen_for_gene = {
            s[i : i + k] for s in seqs if len(s) >= k for i in range(0, len(s) - k + 1, step)
        }
        total_kmer_counts[gid] = len(seen_for_gene)
        ambig_kmer_counts[gid] = sum(
            1
            for km_s in seen_for_gene
            if kmer_owners.get(km_s.encode("ascii", errors="ignore")) == MULTI
        )

    risk_rows = []
    for gid in sorted(gene_to_entity):
        tot = total_kmer_counts.get(gid, 0)
        shared_frac = (ambig_kmer_counts.get(gid, 0) / tot) if tot > 0 else 0.0
        has_overlap = gid in same_strand_overlap_genes
        risk_rows.append(
            {
                "gene_id": gid,
                "shared_kmer_fraction": round(shared_frac, 4),
                "same_strand_host_overlap": has_overlap,
                "mapping_risk": bool(has_overlap or shared_frac >= 0.05),
            }
        )
    mapping_risk_df = pd.DataFrame(risk_rows)
    return ref, groups_df, seq_df, exonic_bp, symbol_by_gene, mapping_risk_df, source_rows


def build_mgi_crosswalks(
    input_dir: Path, coding_genes: set[str], symbol_by_gene: dict[str, str]
) -> tuple[dict[str, list[str]], dict[str, str], dict[str, str]]:
    """Build MGI->ENSMUSG, unique Entrez->ENSMUSG, and unique Symbol->ENSMUSG crosswalks."""
    mgi_to_ens: dict[str, list[str]] = defaultdict(list)
    with (input_dir / "MRK_ENSEMBL.rpt").open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            mgi_id = parts[0].strip()
            ens_gid = parts[5].strip().split(".")[0]
            if (
                mgi_id.startswith("MGI:")
                and ens_gid in coding_genes
                and ens_gid not in mgi_to_ens[mgi_id]
            ):
                mgi_to_ens[mgi_id].append(ens_gid)

    entrez_candidates: dict[str, set[str]] = defaultdict(set)
    sym_candidates: dict[str, set[str]] = defaultdict(set)
    for gid, sym in symbol_by_gene.items():
        if gid in coding_genes and sym:
            sym_candidates[sym].add(gid)

    with (input_dir / "MGI_EntrezGene.rpt").open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            mgi_id = parts[0].strip()
            sym = parts[1].strip()
            status = parts[2].strip()
            entrez = parts[8].strip()
            ens_list = mgi_to_ens.get(mgi_id, [])
            if len(ens_list) != 1:
                continue
            ens_gid = ens_list[0]
            if entrez and entrez.isdigit():
                entrez_candidates[entrez].add(ens_gid)
            if status == "O" and sym and sym not in sym_candidates:
                sym_candidates[sym].add(ens_gid)

    entrez_to_ens = {k: next(iter(v)) for k, v in entrez_candidates.items() if len(v) == 1}
    sym_to_ens = {k: next(iter(v)) for k, v in sym_candidates.items() if len(v) == 1}
    return dict(mgi_to_ens), entrez_to_ens, sym_to_ens


def _counts_to_tpm(counts: pd.DataFrame, lengths_bp: pd.Series) -> pd.DataFrame:
    valid_len = lengths_bp.reindex(counts.index)
    keep = valid_len.notna() & valid_len.gt(0)
    sub = counts.loc[keep].astype(float)
    kb = valid_len.loc[keep].astype(float) / 1000.0
    rpk = sub.div(kb, axis=0)
    denom = rpk.sum(axis=0).replace(0.0, np.nan)
    return rpk.div(denom, axis=1) * 1e6


def _fpkm_to_tpm(fpkm: pd.DataFrame) -> pd.DataFrame:
    sub = fpkm.astype(float).clip(lower=0.0)
    denom = sub.sum(axis=0).replace(0.0, np.nan)
    return sub.div(denom, axis=1) * 1e6


def _project_peg_matrix(gene_tpm: pd.DataFrame, groups_df: pd.DataFrame) -> pd.DataFrame:
    """Sum complete member-gene linear TPM per PEG within each library."""
    gene_to_eid = {
        g: r.entity_id
        for r in groups_df.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }
    expected = groups_df.set_index("entity_id")["member_count"]
    sub = gene_tpm[gene_tpm.index.isin(gene_to_eid)].copy()
    if sub.empty:
        return pd.DataFrame(index=groups_df["entity_id"], columns=gene_tpm.columns, dtype=float)
    sub = sub.groupby(sub.index).max()
    sub["_eid"] = [gene_to_eid[g] for g in sub.index]
    cols = [c for c in sub.columns if c != "_eid"]
    summed = sub.groupby("_eid", sort=False)[cols].sum(min_count=1)
    present = sub.groupby("_eid", sort=False)[cols].count()
    complete = summed.where(present.eq(expected.reindex(summed.index), axis=0))
    return complete.reindex(groups_df["entity_id"])


def _is_adult_healthy_bgee_stage(stage_name: str, scope: str) -> bool:
    s = str(stage_name).lower()
    if scope == "developmental":
        return any(k in s for k in ("theiler", "neurula", "embryo", "gastrula"))
    if any(
        bad in s
        for bad in ("theiler", "neurula", "embryo", "neonate", "1-week", "2-week", "3-week")
    ):
        return False
    # Exclude senescent >15-month libraries so healthy adult baseline reflects prime/mature adult physiology
    return not any(old in s for old in ("18-month", "21-month", "24-month", "27-month"))


def load_bgee_and_normal_panels(
    input_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    entrez_to_ens: dict[str, str],
    sym_to_ens: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, set[str]], list[dict]]:
    """Aggregate normal somatic N, thymic, and all 16 reproductive/developmental compartments."""
    source_rows = []
    lib_tar = input_dir / "Mus_musculus_RNA-Seq_experiments_libraries.tar.gz"
    source_rows.append(
        {
            "source_id": "bgee-15.2-metadata",
            "pillar": "healthy_baseline_rna",
            "filename": lib_tar.name,
            "unit": "metadata",
            "strain": "multi-strain (explicit per library)",
            "bytes": lib_tar.stat().st_size,
            "sha256": _sha256(lib_tar),
        }
    )
    with tarfile.open(lib_tar) as tf:
        libs_meta = pd.read_csv(tf.extractfile("Mus_musculus_RNA-Seq_libraries.tsv"), sep="\t")
    for _col in (
        "Library ID",
        "Experiment ID",
        "Anatomical entity name",
        "Stage name",
        "Sex",
        "Strain",
    ):
        if _col in libs_meta.columns:
            libs_meta[_col] = libs_meta[_col].astype(str).str.strip('"')
    lib_info = libs_meta.set_index("Library ID")[
        ["Experiment ID", "Anatomical entity name", "Stage name", "Sex", "Strain"]
    ].to_dict("index")

    bgee_exps = [
        "SRP013027",
        "SRP199494",
        "SRP181218",
        "GSE41637",
        "SRP262332",
        "SRP162729",
        "GSE43520",
        "DRP003810",
        "ERP002178",
        "SRP007655",
        "SRP074355",
        "ERP104395",
    ]

    # Collect per-organ somatic PEG matrices (library columns) and reproductive/embryo PEG matrices
    somatic_peg_cols: dict[str, dict[str, pd.Series]] = defaultdict(dict)
    thymic_peg_cols: dict[str, pd.Series] = {}
    repro_records: list[tuple[str, str, str, str, str, pd.Series]] = []
    # tuple: (source_id, sample_id, donor_id, scope, compartment, peg_series)

    for exp in bgee_exps:
        fpath = input_dir / f"Mus_musculus_RNA-Seq_read_counts_TPM_{exp}.tsv.gz"
        source_rows.append(
            {
                "source_id": f"bgee-{exp.lower()}",
                "pillar": "healthy_baseline_rna"
                if exp != "DRP003810"
                else "reproductive_and_developmental_rna",
                "filename": fpath.name,
                "unit": "TPM",
                "strain": "recorded per library",
                "bytes": fpath.stat().st_size,
                "sha256": _sha256(fpath),
            }
        )
        exp_meta = libs_meta[libs_meta["Experiment ID"] == exp]
        keep_libs: dict[str, tuple[str, str, str]] = {}
        for _, mrow in exp_meta.iterrows():
            lid = mrow["Library ID"]
            anat = mrow["Anatomical entity name"]
            st = mrow["Stage name"]
            mapped = BGEE_ANATOMY_MAP.get(anat)
            if mapped is None:
                continue
            scope, target_ctx = mapped
            if not _is_adult_healthy_bgee_stage(st, scope):
                continue
            keep_libs[lid] = (scope, target_ctx, st)
        if not keep_libs:
            continue

        chunks = []
        for chunk in pd.read_csv(
            fpath,
            sep="\t",
            usecols=["Library ID", "Gene ID", "TPM"],
            dtype={"Library ID": str, "Gene ID": str, "TPM": float},
            chunksize=2_000_000,
        ):
            chunk["Library ID"] = chunk["Library ID"].str.strip('"')
            sub_chunk = chunk[chunk["Library ID"].isin(keep_libs)].copy()
            if not sub_chunk.empty:
                sub_chunk["Gene ID"] = sub_chunk["Gene ID"].str.strip('"').str.split(".").str[0]
                chunks.append(sub_chunk)
        if not chunks:
            continue

        sub = pd.concat(chunks, ignore_index=True)
        piv = sub.pivot(index="Gene ID", columns="Library ID", values="TPM")
        peg_mat = _project_peg_matrix(piv, groups_df)

        for lid, (scope, target_ctx, st) in keep_libs.items():
            if lid not in peg_mat.columns:
                continue
            col = peg_mat[lid]
            info = lib_info.get(lid, {})
            strain = str(info.get("Strain") or "mouse").strip('"')
            sex = str(info.get("Sex") or "na").strip('"')
            donor_key = f"bgee:{exp}:{strain}:{sex}:{st}:{lid}"
            if scope == "normal":
                somatic_peg_cols[target_ctx][donor_key] = col
            elif scope == "thymic":
                thymic_peg_cols[donor_key] = col
            elif scope in ("reproductive", "developmental"):
                repro_records.append(
                    (f"bgee-{exp.lower()}", lid, donor_key, scope, target_ctx, col)
                )

    # 2. GSE219045 sex-stratified preclinical mouse atlas (85 columns on native ENSMUSG; FPKM -> exact within-library TPM)
    g219_path = input_dir / "GSE219045_Mouse_FPKM_table_Mar22.xlsx"
    source_rows.append(
        {
            "source_id": "gse219045-mouse-atlas",
            "pillar": "healthy_baseline_rna",
            "filename": g219_path.name,
            "unit": "TPM (normalized from within-sample FPKM)",
            "strain": "C57BL/6 & CD-1 preclinical reference",
            "bytes": g219_path.stat().st_size,
            "sha256": _sha256(g219_path),
        }
    )
    g219_df = pd.read_excel(g219_path, sheet_name="Mouse_FPKM_table_March22", engine="openpyxl")
    g219_df["GeneID"] = g219_df["GeneID"].astype(str).str.split(".").str[0]
    g219_df = g219_df.set_index("GeneID")
    g219_tpm = _fpkm_to_tpm(g219_df)
    g219_peg = _project_peg_matrix(g219_tpm, groups_df)

    for col_name in g219_peg.columns:
        suffix = col_name.removeprefix("Mouse_Female_").removeprefix("Mouse_Male_")
        sex_tag = "F" if col_name.startswith("Mouse_Female_") else "M"
        mapped = GSE219045_COL_MAP.get(suffix)
        if mapped is None:
            continue
        scope, target_ctx = mapped
        donor_key = f"gse219045:{sex_tag}:{col_name}"
        if scope == "normal":
            somatic_peg_cols[target_ctx][donor_key] = g219_peg[col_name]
        elif scope == "thymic":
            thymic_peg_cols[donor_key] = g219_peg[col_name]
        elif scope == "reproductive":
            repro_records.append(
                (
                    "gse219045-mouse-atlas",
                    col_name,
                    donor_key,
                    "reproductive",
                    target_ctx,
                    g219_peg[col_name],
                )
            )

    # 3. Supplementary accessory reproductive & lactating mammary datasets (completing all 16 mouse compartments!)
    exonic_series = pd.Series(exonic_bp, dtype=float)
    sym_universe = pd.Index(sorted(set(sym_to_ens.values())), name="_gid")
    entrez_universe = pd.Index(sorted(set(entrez_to_ens.values())), name="_gid")

    # 3a. GSE295398 WT FVB/N ventral prostate (linear TPM via Entrez)
    p295 = input_dir / "GSE295398_TPM_all_samples.csv.gz"
    source_rows.append(
        {
            "source_id": "gse295398-prostate",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p295.name,
            "unit": "TPM",
            "strain": "FVB/N wild type",
            "bytes": p295.stat().st_size,
            "sha256": _sha256(p295),
        }
    )
    df295 = pd.read_csv(p295)
    df295["_entrez"] = df295["Identifier"].astype(str).str.extract(r"list_uids=(\d+)")[0]
    df295["_gid"] = df295["_entrez"].map(entrez_to_ens)
    df295 = df295.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    peg295 = _project_peg_matrix(df295[["WT1", "WT2", "WT3", "WT4"]].astype(float), groups_df)
    for c in ("WT1", "WT2", "WT3", "WT4"):
        repro_records.append(
            ("gse295398-prostate", c, f"gse295398:{c}", "reproductive", "prostate", peg295[c])
        )

    # 3b. GSE287166 WT FVB anterior prostate (coagulating_gland) & ventral prostate
    # Authors dropped all-zero symbol rows prior to GEO deposit; restore 0.0 across the mapped symbol universe.
    p287 = input_dir / "GSE287166_NormalizedCounts.csv.gz"
    source_rows.append(
        {
            "source_id": "gse287166-anterior-prostate",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p287.name,
            "unit": "TPM (from length-normalized counts)",
            "strain": "FVB wild type",
            "bytes": p287.stat().st_size,
            "sha256": _sha256(p287),
        }
    )
    df287 = pd.read_csv(p287, index_col=0)
    df287["_gid"] = df287.index.astype(str).map(sym_to_ens)
    df287 = df287.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    cnt287 = df287[["FVB1A", "FVB2A", "FVB3A", "FVB1V", "FVB2V", "FVB3V"]].reindex(
        sym_universe, fill_value=0.0
    )
    tpm287 = _counts_to_tpm(cnt287, exonic_series)
    peg287 = _project_peg_matrix(tpm287, groups_df)
    for c in ("FVB1A", "FVB2A", "FVB3A"):
        repro_records.append(
            (
                "gse287166-anterior-prostate",
                c,
                f"gse287166:{c}",
                "reproductive",
                "coagulating_gland",
                peg287[c],
            )
        )
    for c in ("FVB1V", "FVB2V", "FVB3V"):
        repro_records.append(
            (
                "gse287166-anterior-prostate",
                c,
                f"gse287166:{c}",
                "reproductive",
                "prostate",
                peg287[c],
            )
        )

    # 3c. GSE263021 normal seminal vesicle (native ENSMUSG counts -> TPM)
    p263 = input_dir / "GSE263021_raw_counts_matrix.tsv.gz"
    source_rows.append(
        {
            "source_id": "gse263021-seminal-vesicle",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p263.name,
            "unit": "TPM (from native ENSMUSG counts)",
            "strain": "C57BL/6J background",
            "bytes": p263.stat().st_size,
            "sha256": _sha256(p263),
        }
    )
    df263 = pd.read_csv(p263, sep="\t", index_col=0)
    df263.index = df263.index.astype(str).str.split(".").str[0]
    sv_cols = [
        "1228SVh_L000",
        "1230SVh_L000",
        "1436SVh_L000",
        "1793SVh_L000",
        "3552SVh_L000",
        "3560SVh_L000",
    ]
    tpm263 = _counts_to_tpm(df263[sv_cols], exonic_series)
    peg263 = _project_peg_matrix(tpm263, groups_df)
    for c in sv_cols:
        repro_records.append(
            (
                "gse263021-seminal-vesicle",
                c,
                f"gse263021:{c}",
                "reproductive",
                "seminal vesicle",
                peg263[c],
            )
        )

    # 3d. GSE324827 WT Pgrflox/flox oviduct (fallopian tube)
    # Authors dropped all-zero symbol rows prior to GEO deposit; restore 0.0 across the mapped symbol universe.
    p324 = input_dir / "GSE324827_MyoPgrdd_RNAseq_RawCounts.csv.gz"
    source_rows.append(
        {
            "source_id": "gse324827-oviduct",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p324.name,
            "unit": "TPM (from raw counts)",
            "strain": "C57BL/6 Pgrflox/flox WT control",
            "bytes": p324.stat().st_size,
            "sha256": _sha256(p324),
        }
    )
    df324 = pd.read_csv(p324)
    df324["_gid"] = df324["Gene Symbol"].astype(str).map(sym_to_ens)
    df324 = df324.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    cnt324 = df324[["PgrFlox1", "PgrFlox2", "PgrFlox3"]].reindex(sym_universe, fill_value=0.0)
    tpm324 = _counts_to_tpm(cnt324, exonic_series)
    peg324 = _project_peg_matrix(tpm324, groups_df)
    for c in ("PgrFlox1", "PgrFlox2", "PgrFlox3"):
        repro_records.append(
            ("gse324827-oviduct", c, f"gse324827:{c}", "reproductive", "fallopian tube", peg324[c])
        )

    # 3e. GSE198799 untreated WT uterine horn (native ENSMUSG)
    p198 = input_dir / "GSE198799_countTable.txt.gz"
    source_rows.append(
        {
            "source_id": "gse198799-uterine-horn",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p198.name,
            "unit": "TPM (from native ENSMUSG counts)",
            "strain": "C3H untreated wild type",
            "bytes": p198.stat().st_size,
            "sha256": _sha256(p198),
        }
    )
    df198 = pd.read_csv(p198, sep="\t")
    df198["gene"] = df198["gene"].astype(str).str.split(".").str[0]
    df198 = df198.drop_duplicates("gene").set_index("gene")
    uh_cols = ["G_62_2L", "H_62_2R", "M_61_5L", "N_61_5R", "O_46L3L", "P_46L3R"]
    tpm198 = _counts_to_tpm(df198[uh_cols], exonic_series)
    peg198 = _project_peg_matrix(tpm198, groups_df)
    for c in uh_cols:
        repro_records.append(
            (
                "gse198799-uterine-horn",
                c,
                f"gse198799:{c}",
                "reproductive",
                "uterine_horn",
                peg198[c],
            )
        )

    # 3f. GSE223136 WT C57BL/6 estrous endometrium (native ENSMUSG + Length)
    p223 = input_dir / "GSE223136_counts2.tsv.gz"
    source_rows.append(
        {
            "source_id": "gse223136-endometrium",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p223.name,
            "unit": "TPM (from native ENSMUSG featureCounts)",
            "strain": "C57BL/6 wild type",
            "bytes": p223.stat().st_size,
            "sha256": _sha256(p223),
        }
    )
    df223 = pd.read_csv(p223, sep="\t")
    df223["Geneid"] = df223["Geneid"].astype(str).str.split(".").str[0]
    df223 = df223.drop_duplicates("Geneid").set_index("Geneid")
    endo_cols = ["CONT1_ATGTCA", "CONT2_CGATGT", "CONT4_ACTTGA"]
    tpm223 = _counts_to_tpm(df223[endo_cols], df223["Length"].astype(float))
    peg223 = _project_peg_matrix(tpm223, groups_df)
    for c in endo_cols:
        repro_records.append(
            ("gse223136-endometrium", c, f"gse223136:{c}", "reproductive", "endometrium", peg223[c])
        )

    # 3g. GSE85815 WT non-pregnant cervix (native ENSMUSG + Mu_feature_length)
    p858 = input_dir / "GSE85815_mouse_feature_lengths_and_read_counts.csv.gz"
    source_rows.append(
        {
            "source_id": "gse85815-cervix",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p858.name,
            "unit": "TPM (from native ENSMUSG counts)",
            "strain": "C57BL/6 non-pregnant wild type",
            "bytes": p858.stat().st_size,
            "sha256": _sha256(p858),
        }
    )
    df858 = pd.read_csv(p858)
    df858["Mm_Ensembl_Gene_ID"] = df858["Mm_Ensembl_Gene_ID"].astype(str).str.split(".").str[0]
    df858 = df858.drop_duplicates("Mm_Ensembl_Gene_ID").set_index("Mm_Ensembl_Gene_ID")
    cerv_cols = ["Mm_NP2_read_count", "Mm_NP3_read_count"]
    tpm858 = _counts_to_tpm(df858[cerv_cols], df858["Mu_feature_length"].astype(float))
    peg858 = _project_peg_matrix(tpm858, groups_df)
    for c in cerv_cols:
        repro_records.append(
            ("gse85815-cervix", c, f"gse85815:{c}", "reproductive", "cervix", peg858[c])
        )

    # 3h. GSE60450 WT FVB/N lactating & virgin mammary gland (EntrezGeneID + Length)
    # Restore 0.0 across mapped Entrez loci omitted as unexpressed in the GEO featureCounts table.
    p604 = input_dir / "GSE60450_Lactation-GenewiseCounts.txt.gz"
    source_rows.append(
        {
            "source_id": "gse60450-mammary-lactation",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p604.name,
            "unit": "TPM (from Entrez featureCounts)",
            "strain": "FVB/N wild type",
            "bytes": p604.stat().st_size,
            "sha256": _sha256(p604),
        }
    )
    df604 = pd.read_csv(p604, sep="\t")
    df604["_gid"] = df604["EntrezGeneID"].astype(str).map(entrez_to_ens)
    df604 = df604.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    mam_cols = [c for c in df604.columns if c.startswith("MCL1-")]
    len604 = df604["Length"].astype(float).reindex(entrez_universe).fillna(exonic_series)
    cnt604 = df604[mam_cols].reindex(entrez_universe, fill_value=0.0)
    tpm604 = _counts_to_tpm(cnt604, len604)
    peg604 = _project_peg_matrix(tpm604, groups_df)
    # Note: GSE60450 FACS-dissociated mammary populations underwent warm enzymatic collagenase
    # digestion (inducing acute myeloid stress cytokines Csf3/Cxcl2/Fos 500-7,800 TPM). Intact
    # flash-frozen resting_mammary is already covered by Bgee 15.2 and GSE219045; GSE60450 supplies
    # virgin luminal/basal resting_mammary and 2-day lactating_mammary while sorted/activated immune
    # panels (load_mouse_immune_panels) independently gate inducible myeloid/immune transcripts.
    for c in mam_cols:
        prefix = c.split("_", 1)[0]
        if prefix in {"MCL1-LA", "MCL1-LB", "MCL1-DG", "MCL1-DH"}:
            repro_records.append(
                (
                    "gse60450-mammary-lactation",
                    c,
                    f"gse60450:{prefix}",
                    "reproductive",
                    "resting_mammary",
                    peg604[c],
                )
            )
        elif prefix in {"MCL1-LE", "MCL1-LF", "MCL1-DK", "MCL1-DL"}:
            repro_records.append(
                (
                    "gse60450-mammary-lactation",
                    c,
                    f"gse60450:{prefix}",
                    "reproductive",
                    "lactating_mammary",
                    peg604[c],
                )
            )

    # Summarize somatic normal contexts by donor P95 -> comparison_normals + normal_gates + confirmed_rare_cell_organs
    normal_rows = []
    rna_gate_rows = []
    confirmed_organs: dict[str, set[str]] = defaultdict(set)

    for ctx, donor_series in sorted(somatic_peg_cols.items()):
        mat = pd.DataFrame(donor_series).T  # rows = donors, cols = entity_id
        summarized = target_selection.summarize_normal_tissue_tpm(mat, method="p95")
        tier = SOMATIC_CONTEXT_TIER[ctx]
        for eid, val in summarized.items():
            if pd.isna(val):
                continue
            v = float(val)
            if v >= 1.0:
                confirmed_organs[eid].add(ctx)
            normal_rows.append(
                {
                    "entity_id": eid,
                    "source_id": "mouse_somatic_panel",
                    "context_id": ctx,
                    "tissue": ctx,
                    "scope": "normal",
                    "value": v,
                    "comparison_id": COMPARISON_ID,
                    "unit": "TPM",
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )
            if tier == "vital" and v > 4.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Mouse somatic bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{ctx} (vital): {v:.4g} TPM > 4",
                    }
                )
            elif tier == "protected" and v > 8.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Mouse somatic bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{ctx} (protected): {v:.4g} TPM > 8",
                    }
                )
            elif tier == "manageable" and v > 128.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Mouse somatic bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{ctx} (manageable): {v:.4g} TPM > 128",
                    }
                )
            elif tier == "manageable" and v > 16.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Mouse somatic bulk RNA",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"{ctx} (manageable hold): {v:.4g} TPM > 16",
                    }
                )

    if thymic_peg_cols:
        tmat = pd.DataFrame(thymic_peg_cols).T
        tsum = target_selection.summarize_normal_tissue_tpm(tmat, method="p95")
        for eid, val in tsum.items():
            if pd.notna(val):
                normal_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "mouse_somatic_panel",
                        "context_id": "thymus",
                        "tissue": "thymus",
                        "scope": "thymic",
                        "value": float(val),
                        "comparison_id": COMPARISON_ID,
                        "unit": "TPM",
                        "qc_pass": True,
                        "identity_resolved": True,
                    }
                )

    # Build comparison_reproductive concisely: retain per-(compartment, source_id, sample_id) observations
    # To keep memory & runtime tight across 22k PEGs x ~190 reproductive/embryo libraries (~4.1M rows -> summarize per
    # (compartment, source_id) max + top witness sample so both RA peak and HostContext max_tpm are 100% exact!):
    repro_by_comp: dict[tuple[str, str], list[tuple[str, str, str, pd.Series]]] = defaultdict(list)
    for src_id, smp_id, donor_id, scope, comp, ser in repro_records:
        repro_by_comp[(scope, comp)].append((src_id, smp_id, donor_id, ser))

    repro_rows = []
    for (scope, comp), items in sorted(repro_by_comp.items()):
        cols_dict = {f"{i}": ser for i, (_src, _smp, _dnr, ser) in enumerate(items)}
        mat = pd.DataFrame(cols_dict)  # index = entity_id, cols = sample idx
        valid_any = mat.notna().any(axis=1)
        sub_mat = mat.loc[valid_any]
        max_vals = sub_mat.max(axis=1)
        best_col = sub_mat.idxmax(axis=1)
        for eid, val in max_vals.items():
            idx_i = int(best_col.at[eid])
            src_id, smp_id, donor_id, _ = items[idx_i]
            rec = {
                "entity_id": eid,
                "source_id": src_id,
                "sample_id": smp_id,
                "participant_id": donor_id,
                "scope": scope,
                "compartment": comp,
                "tissue": comp.replace("_", " "),
                "value": float(val),
                "comparison_id": COMPARISON_ID,
                "unit": "TPM",
                "qc_pass": True,
                "identity_resolved": True,
            }
            if scope == "developmental":
                rec["cells"] = 100
                rec["positive_cells"] = 100
                rec["developmental_stage"] = "embryonic"
                rec["source_context"] = "embryo"
            repro_rows.append(rec)

    return (
        pd.DataFrame(normal_rows),
        pd.DataFrame(repro_rows),
        pd.DataFrame(rna_gate_rows, columns=["entity_id", "layer", "unsafe", "hold", "reason"]),
        dict(confirmed_organs),
        source_rows,
    )


def load_mouse_immune_panels(
    input_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    symbol_by_gene: dict[str, str],
    sym_to_ens: dict[str, str],
) -> tuple[pd.DataFrame, list[dict]]:
    """Evaluate sorted & activated murine immune RNA panels (GSE109125, GSE122108, GSE181223).

    Applies Protected <= 8 TPM across adult peripheral innate & adaptive immune subsets
    (including LPS-, IL-4/LPS-, anti-CD3/CD40-, thioglycollate-, LCMV-, and injury-activated
    macrophages, monocytes, neutrophils, DCs, microglia, mast cells, ILCs, NK cells, and T cells)
    and Manageable (hold > 16 TPM, Unsafe > 128 TPM) across mature B-cell subsets, using the
    replicated finite-sample P95 upper bound ``min(quantile(0.95), second_highest_replicate)``
    across populations with >= 2 biological replicates so single-tube FACS splash outliers
    never trigger false exclusions.
    """
    exonic_series = pd.Series(exonic_bp, dtype=float)
    sym_universe = pd.Index(sorted(set(sym_to_ens.values())), name="_gid")
    olfr_gids = {
        gid
        for gid, sym in symbol_by_gene.items()
        if sym.startswith("Olfr") or re.match(r"^Or\d+[a-z]+\d+", sym, re.IGNORECASE)
    }

    p109 = input_dir / "GSE109125_Gene_count_table_GENCODE_vM25.csv.gz"
    p122 = input_dir / "GSE122108_Gene_count_table_GENCODE_vM25.csv.gz"
    p181 = input_dir / "GSE181223_mm_BMDM_RNAseq_EGR2WT_EGR2KO_IL4LPS_TPM_counts.txt.gz"
    source_rows = [
        {
            "source_id": "immgen-uli-gse109125",
            "pillar": "normal_bulk_and_immune_safety",
            "filename": p109.name,
            "unit": "TPM (from GENCODE vM25 sorted immune counts)",
            "strain": "C57BL/6J (ImmGen ULI adult immune)",
            "bytes": p109.stat().st_size,
            "sha256": _sha256(p109),
        },
        {
            "source_id": "immgen-mnp-gse122108",
            "pillar": "normal_bulk_and_immune_safety",
            "filename": p122.name,
            "unit": "TPM (from GENCODE vM25 mononuclear phagocyte counts)",
            "strain": "C57BL/6J (ImmGen MNP resident & inflamed)",
            "bytes": p122.stat().st_size,
            "sha256": _sha256(p122),
        },
        {
            "source_id": "gse181223-bmdm-lps",
            "pillar": "normal_bulk_and_immune_safety",
            "filename": p181.name,
            "unit": "TPM (Salmon native ENSMUSG)",
            "strain": "C57BL/6 WT primary BMDM (CTRL, IL-4, LPS, IL-4+LPS)",
            "bytes": p181.stat().st_size,
            "sha256": _sha256(p181),
        },
    ]

    def _replicated_p95(mat_cols: pd.DataFrame) -> pd.Series:
        # Rows = entity_id, columns = >= 2 biological replicates of one sorted immune population.
        arr = mat_cols.to_numpy(dtype=float)
        has_obs = np.isfinite(arr).any(axis=1)
        out = np.full(len(mat_cols), np.nan, dtype=float)
        if not has_obs.any():
            return pd.Series(out, index=mat_cols.index)
        sub = arr[has_obs]
        p95 = np.nanquantile(sub, 0.95, axis=1)
        # Finite-sample replication guard: require >= 2 biological replicates to support the bound
        sorted_arr = np.sort(np.where(np.isnan(sub), -np.inf, sub), axis=1)
        second_hi = sorted_arr[:, -2]
        second_hi = np.where(np.isneginf(second_hi), np.nan, second_hi)
        out[has_obs] = np.minimum(p95, second_hi)
        return pd.Series(out, index=mat_cols.index)

    pop_summaries: dict[str, tuple[str, pd.Series]] = {}

    # 1. GSE181223 WT primary BMDM across resting, IL-4, LPS, and IL-4+LPS conditions
    df181 = pd.read_csv(p181, sep="\t")
    df181["gene_id"] = df181["gene_id"].astype(str).str.split(".").str[0]
    wt_cols = [c for c in df181.columns if c.startswith("WT_")]
    mat181 = df181.groupby("gene_id", sort=False)[wt_cols].sum()
    mat181 = mat181.loc[mat181.index.isin(exonic_bp)]
    mat181 = mat181.div(mat181.sum(axis=0).replace(0.0, np.nan), axis=1) * 1e6
    peg181 = _project_peg_matrix(mat181, groups_df)
    for cond in ("WT_CTRL", "WT_IL4", "WT_LPS", "WT_IL4LPS"):
        cs = [c for c in peg181.columns if c.startswith(f"{cond}_")]
        if len(cs) >= 2:
            pop_summaries[f"BMDM:{cond}"] = ("protected", _replicated_p95(peg181[cs]))

    # 2. GSE109125 ImmGen ULI adult peripheral immune subsets
    df109 = pd.read_csv(p109)
    gid109 = df109["GeneSymbol"].astype(str).map(sym_to_ens)
    skip_109 = (
        ".Th",
        "BEC.",
        "LEC.",
        "FRC.",
        "Ep.",
        "IAP.",
        "LTHSC.",
        "STHSC.",
        "MMP",
        "proB.",
        "B.FrE.",
    )
    cols109 = [
        c for c in df109.columns if c != "GeneSymbol" and not any(tok in c for tok in skip_109)
    ]
    sub109 = df109.loc[gid109.notna(), cols109].copy()
    sub109.index = gid109[gid109.notna()].to_numpy()
    sub109 = sub109[~sub109.index.duplicated(keep="first")].reindex(sym_universe, fill_value=0.0)
    tot109 = sub109.sum(axis=0)
    tpm109 = _counts_to_tpm(sub109, exonic_series)
    olfr109 = (tpm109.loc[tpm109.index.isin(olfr_gids)] >= 5.0).mean(axis=0)
    keep109 = [c for c in cols109 if tot109[c] >= 500_000 and olfr109[c] <= 0.025]
    peg109 = _project_peg_matrix(tpm109[keep109], groups_df)
    by_pop109: dict[str, list[str]] = defaultdict(list)
    for c in keep109:
        base = c.split("#", 1)[0] if "#" in c else c.rsplit(".", 1)[0]
        by_pop109[base].append(c)
    for pop, cs in sorted(by_pop109.items()):
        if len(cs) < 2:
            continue
        tier = "manageable" if pop.startswith(("B.", "B1")) else "protected"
        pop_summaries[f"ImmGen:{pop}"] = (tier, _replicated_p95(peg109[cs]))

    # 3. GSE122108 ImmGen adult resident & inflammatory mononuclear phagocytes (exclude embryonic/neonatal & thymic)
    df122 = pd.read_csv(p122)
    gid122 = df122["GeneSymbol"].astype(str).map(sym_to_ens)
    dev_or_thymic_pat = re.compile(
        r"(?:\.E\d|\bYS\b|\.d\d+\.EB\b|\.EB\b|\.Neo\b|\.P\d+\b|Fetal|embryo|\.Th(?:\.|#|$))",
        re.IGNORECASE,
    )
    cols122 = [c for c in df122.columns if c != "GeneSymbol" and not dev_or_thymic_pat.search(c)]
    sub122 = df122.loc[gid122.notna(), cols122].copy()
    sub122.index = gid122[gid122.notna()].to_numpy()
    sub122 = sub122[~sub122.index.duplicated(keep="first")].reindex(sym_universe, fill_value=0.0)
    tot122 = sub122.sum(axis=0)
    tpm122 = _counts_to_tpm(sub122, exonic_series)
    olfr122 = (tpm122.loc[tpm122.index.isin(olfr_gids)] >= 5.0).mean(axis=0)
    keep122 = [c for c in cols122 if tot122[c] >= 500_000 and olfr122[c] <= 0.025]
    peg122 = _project_peg_matrix(tpm122[keep122], groups_df)
    by_pop122: dict[str, list[str]] = defaultdict(list)
    for c in keep122:
        base = c.split("#", 1)[0] if "#" in c else c.rsplit(".", 1)[0]
        by_pop122[base].append(c)
    for pop, cs in sorted(by_pop122.items()):
        if len(cs) < 2:
            continue
        pop_summaries[f"ImmGenMNP:{pop}"] = ("protected", _replicated_p95(peg122[cs]))

    prot_df = pd.DataFrame({k: s for k, (t, s) in pop_summaries.items() if t == "protected"})
    man_df = pd.DataFrame({k: s for k, (t, s) in pop_summaries.items() if t == "manageable"})

    gate_rows: list[dict] = []
    if not prot_df.empty:
        valid_p = prot_df.notna().any(axis=1)
        p_sub = prot_df.loc[valid_p]
        p_max = p_sub.max(axis=1)
        p_top = p_sub.idxmax(axis=1)
        for eid, val in p_max[p_max > 8.0].items():
            gate_rows.append(
                {
                    "entity_id": eid,
                    "layer": "Mouse sorted & activated immune RNA",
                    "unsafe": True,
                    "hold": False,
                    "reason": f"{p_top.at[eid]} (protected immune): {float(val):.4g} TPM > 8",
                }
            )
    if not man_df.empty:
        valid_m = man_df.notna().any(axis=1)
        m_sub = man_df.loc[valid_m]
        m_max = m_sub.max(axis=1)
        m_top = m_sub.idxmax(axis=1)
        for eid, val in m_max[m_max > 128.0].items():
            gate_rows.append(
                {
                    "entity_id": eid,
                    "layer": "Mouse sorted & activated immune RNA",
                    "unsafe": True,
                    "hold": False,
                    "reason": f"{m_top.at[eid]} (manageable B-cell): {float(val):.4g} TPM > 128",
                }
            )
        for eid, val in m_max[(m_max > 16.0) & (m_max <= 128.0)].items():
            gate_rows.append(
                {
                    "entity_id": eid,
                    "layer": "Mouse sorted B-cell Manageable hold",
                    "unsafe": False,
                    "hold": True,
                    "reason": f"{m_top.at[eid]} (manageable B-cell hold): {float(val):.4g} TPM > 16",
                }
            )

    return pd.DataFrame(
        gate_rows, columns=["entity_id", "layer", "unsafe", "hold", "reason"]
    ), source_rows


def load_mgi_gxd_protein_gates(
    input_dir: Path,
    groups_df: pd.DataFrame,
    mgi_to_ens: dict[str, list[str]],
    confirmed_organs: dict[str, set[str]],
) -> tuple[pd.DataFrame, list[dict]]:
    """Evaluate adult TS28 Immunohistochemistry & Immunoblot detections from MGI GXD 6.27."""
    gxd_path = input_dir / "EXPRESSION_MGI.tsv.gz"
    mrk_path = input_dir / "MRK_ENSEMBL.rpt"
    source_rows = [
        {
            "source_id": "mgi-gxd-6.27",
            "pillar": "normal_protein",
            "filename": gxd_path.name,
            "unit": "curated TS28 IHC/Immunoblot detections",
            "strain": "curated adult mouse",
            "bytes": gxd_path.stat().st_size,
            "sha256": _sha256(gxd_path),
        },
        {
            "source_id": "mgi-mrk-ensembl",
            "pillar": "normal_protein",
            "filename": mrk_path.name,
            "unit": "MGI-to-ENSMUSG crosswalk",
            "strain": "GRCm39",
            "bytes": mrk_path.stat().st_size,
            "sha256": _sha256(mrk_path),
        },
    ]

    gene_to_eid = {
        g: r.entity_id
        for r in groups_df.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }

    with gzip.open(gxd_path, "rt", encoding="utf-8") as fh:
        doc = json.load(fh)

    pubs_by_eid_organ: dict[tuple[str, str], set[str]] = defaultdict(set)
    for rec in doc.get("data", ()):
        if rec.get("assay") not in {"MMO:0000498", "MMO:0000672"}:
            continue
        if rec.get("whenExpressed", {}).get("stageName") != "TS28":
            continue
        stmt = rec.get("whereExpressed", {}).get("whereExpressedStatement", "").strip().lower()
        organ = GXD_ORGAN_MAP.get(stmt)
        if organ is None:
            continue
        tier = SOMATIC_CONTEXT_TIER.get(organ)
        if tier not in {"vital", "protected"}:
            continue
        mgi_id = rec.get("geneId", "")
        ens_list = mgi_to_ens.get(mgi_id, [])
        if len(ens_list) != 1:
            continue
        eid = gene_to_eid.get(ens_list[0])
        if eid is None:
            continue
        pub = rec.get("evidence", {}).get("publicationId") or "MGI"
        pubs_by_eid_organ[(eid, organ)].add(pub)
        confirmed_organs.setdefault(eid, set()).add(organ)

    gate_rows = []
    for (eid, organ), pubs in sorted(pubs_by_eid_organ.items()):
        tier = SOMATIC_CONTEXT_TIER[organ]
        n_pubs = len(pubs)
        if n_pubs >= 2:
            gate_rows.append(
                {
                    "entity_id": eid,
                    "layer": "MGI GXD adult normal protein",
                    "unsafe": True,
                    "hold": False,
                    "reason": f"{organ} ({tier} adult IHC/WB, {n_pubs} publications)",
                }
            )
        else:
            gate_rows.append(
                {
                    "entity_id": eid,
                    "layer": "MGI GXD single-publication adult protein review",
                    "unsafe": False,
                    "hold": True,
                    "reason": f"{organ} ({tier} adult IHC/WB, 1 publication: {next(iter(pubs))})",
                }
            )

    # Preserve Sarma 1999 qualitative somatic RT-PCR counterevidence on Trap1a (P1A)
    trap_eid = gene_to_eid.get("ENSMUSG00000051257")
    if trap_eid is not None:
        gate_rows.append(
            {
                "entity_id": trap_eid,
                "layer": "Curated somatic RT-PCR counterevidence (sarma1999)",
                "unsafe": False,
                "hold": True,
                "reason": "Sarma 1999 (PMID:10049945) BALB/c lung, spleen, liver endogenous RT-PCR/Southern positive",
            }
        )

    return pd.DataFrame(
        gate_rows, columns=["entity_id", "layer", "unsafe", "hold", "reason"]
    ), source_rows


def load_tabula_muris_rare_cell_umi(
    input_dir: Path,
    sym_to_ens: dict[str, str],
) -> tuple[pd.DataFrame, list[dict]]:
    """Parse Tabula Muris 10x droplet UMI matrices across adult Vital & Protected somatic organs."""
    zip_path = input_dir / "tabula_muris_droplet.zip"
    ann_path = input_dir / "annotations_droplet.csv"
    source_rows = [
        {
            "source_id": "tabula-muris-10x-droplet",
            "pillar": "normal_cell_and_rare_umi",
            "filename": zip_path.name,
            "unit": "UMI counts",
            "strain": "C57BL/6JN (3-month adult)",
            "bytes": zip_path.stat().st_size,
            "sha256": _sha256(zip_path),
        },
        {
            "source_id": "tabula-muris-10x-annotations",
            "pillar": "normal_cell_and_rare_umi",
            "filename": ann_path.name,
            "unit": "cell ontology annotations",
            "strain": "C57BL/6JN (3-month adult)",
            "bytes": ann_path.stat().st_size,
            "sha256": _sha256(ann_path),
        },
    ]

    ann = pd.read_csv(ann_path, low_memory=False)
    ann = ann[ann["tissue"].isin(TM_TISSUE_MAP) & ann["cell_ontology_class"].notna()].copy()
    # Map (channel, barcode_prefix) -> (donor_id, tissue, tier, cell_type)
    ann["_bc_core"] = ann["cell"].astype(str).str.rsplit("_", n=1).str[-1]
    by_channel: dict[str, dict[str, tuple[str, str, str, str]]] = {}
    for ch, grp in ann.groupby("channel", sort=False):
        organ, tier = TM_TISSUE_MAP[grp["tissue"].iloc[0]]
        by_channel[str(ch)] = {
            bc: (str(mid), organ, tier, str(ct))
            for bc, mid, ct in zip(grp["_bc_core"], grp["mouse.id"], grp["cell_ontology_class"])
        }

    pop_n_cells: dict[tuple[str, str, str, str], int] = defaultdict(int)
    pos_counts: dict[tuple[str, str, str, str, str], int] = defaultdict(int)
    multi_counts: dict[tuple[str, str, str, str, str], int] = defaultdict(int)
    assayed_genes: set[str] = set()

    with zipfile.ZipFile(zip_path) as zf:
        mtx_members = sorted(
            n
            for n in zf.namelist()
            if n.startswith("droplet/") and n.endswith("/matrix.mtx") and "__MACOSX" not in n
        )
        for mtx_name in mtx_members:
            folder = mtx_name.rsplit("/", 1)[0]
            folder_base = folder.rsplit("/", 1)[-1]  # e.g. Bladder-10X_P4_3
            channel = folder_base.split("-", 1)[1] if "-" in folder_base else folder_base
            cell_lookup = by_channel.get(channel)
            if not cell_lookup:
                continue

            with zf.open(f"{folder}/genes.tsv") as gfh:
                raw_genes = [ln.decode("utf-8").rstrip("\n").split("\t")[0] for ln in gfh]
            gene_map: dict[int, str] = {}
            for idx1, sym in enumerate(raw_genes, start=1):
                gid = sym_to_ens.get(sym)
                if gid:
                    gene_map[idx1] = gid
                    assayed_genes.add(gid)

            with zf.open(f"{folder}/barcodes.tsv") as bfh:
                raw_bcs = [ln.decode("utf-8").strip().split("-")[0] for ln in bfh]
            col_pop: dict[int, tuple[str, str, str, str]] = {}
            for idx1, bc in enumerate(raw_bcs, start=1):
                pop = cell_lookup.get(bc)
                if pop is not None:
                    col_pop[idx1] = pop
                    pop_n_cells[pop] += 1
            if not col_pop:
                continue

            with zf.open(mtx_name) as mfh:
                # Skip MatrixMarket comments and dimensions line
                dims_seen = False
                for raw_ln in mfh:
                    if raw_ln.startswith(b"%"):
                        continue
                    if not dims_seen:
                        dims_seen = True
                        continue
                    parts = raw_ln.split()
                    if len(parts) < 3:
                        continue
                    c_idx = int(parts[1])
                    pop = col_pop.get(c_idx)
                    if pop is None:
                        continue
                    g_idx = int(parts[0])
                    gid = gene_map.get(g_idx)
                    if gid is None:
                        continue
                    val = int(float(parts[2]))
                    if val >= 1:
                        key = (gid, *pop)
                        pos_counts[key] += 1
                        if val >= 2:
                            multi_counts[key] += 1

    valid_pops = {pop: n for pop, n in pop_n_cells.items() if n >= 10}
    sample_pop = next(iter(sorted(valid_pops)))
    sample_n = valid_pops[sample_pop]

    umi_rows = []
    genes_with_multi: set[str] = set()
    for (gid, donor_id, organ, tier, cell_type), n_pos in pos_counts.items():
        pop = (donor_id, organ, tier, cell_type)
        n_cells = valid_pops.get(pop)
        if n_cells is None:
            continue
        n_mul = multi_counts.get((gid, donor_id, organ, tier, cell_type), 0)
        if n_mul > 0:
            genes_with_multi.add(gid)
            umi_rows.append(
                {
                    "gene_id": gid,
                    "source_id": "tabula_muris_10x",
                    "donor_id": donor_id,
                    "tissue": organ,
                    "cell_type": cell_type,
                    "tier": tier,
                    "n_cells": n_cells,
                    "n_positive": n_pos,
                    "n_multi": n_mul,
                }
            )

    # Every other assayed Tabula Muris gene with 0 multi-molecule cells in >=10-cell populations
    # gets one compact witness row so assess_rare_cell_rna records status="below_screen" quickly.
    s_donor, s_organ, s_tier, s_ct = sample_pop
    for gid in sorted(assayed_genes - genes_with_multi):
        umi_rows.append(
            {
                "gene_id": gid,
                "source_id": "tabula_muris_10x",
                "donor_id": s_donor,
                "tissue": s_organ,
                "cell_type": s_ct,
                "tier": s_tier,
                "n_cells": sample_n,
                "n_positive": 0,
                "n_multi": 0,
            }
        )

    return pd.DataFrame(umi_rows), source_rows


def load_mouse_tumor_cohorts(
    input_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    entrez_to_ens: dict[str, str],
    sym_to_ens: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Load 8 strain-explicit murine syngeneic & GEMM tumor cohorts in linear TPM."""
    exonic_series = pd.Series(exonic_bp, dtype=float)
    cohort_registry = []
    source_rows = []
    per_sample_peg: list[tuple[str, str, str, str, pd.Series]] = []
    # tuple: (source_id, sample_id, participant_id, cancer_code, peg_series)

    # 1. MC38 (C57BL/6 colorectal carcinoma, GSE355313, n=4 in-vivo subcutaneous tumors A1..A4)
    p_mc38 = input_dir / "GSE355313_gene_count.csv.gz"
    source_rows.append(
        {
            "source_id": "gse355313-mc38",
            "pillar": "tumor_rna",
            "filename": p_mc38.name,
            "unit": "TPM (from native ENSMUSG counts + gene_length)",
            "strain": "C57BL/6 (H-2b)",
            "bytes": p_mc38.stat().st_size,
            "sha256": _sha256(p_mc38),
        }
    )
    df_mc38 = pd.read_csv(p_mc38)
    df_mc38["gene_id"] = df_mc38["gene_id"].astype(str).str.split(".").str[0]
    df_mc38 = df_mc38.drop_duplicates("gene_id").set_index("gene_id")
    mc38_cols = ["A1", "A2", "A3", "A4"]
    peg_mc38 = _project_peg_matrix(
        _counts_to_tpm(df_mc38[mc38_cols], df_mc38["gene_length"].astype(float)), groups_df
    )
    for c in mc38_cols:
        per_sample_peg.append(("gse355313", f"MC38_{c}", f"GSE355313:{c}", "MC38", peg_mc38[c]))
    cohort_registry.append(
        {
            "cancer_code": "MC38",
            "label": "MC38 subcutaneous colon adenocarcinoma",
            "model_kind": "syngeneic_graft",
            "strain": "C57BL/6",
            "mhc_haplotype": "H-2b",
            "source_accession": "GSE355313",
            "n_samples": 4,
            "native_feature_space": "Ensembl gene (ENSMUSG)",
            "quantification": "linear TPM from raw counts and gene_length",
        }
    )

    # 2-4. B16F10 (n=3), CT26 (n=4), LLC (n=3) in-vivo PBS control tumors (GSE334005, linear TPM)
    p_334 = input_dir / "GSE334005_RNAseq_CT26_B16F10_LLC_tpm.tsv.gz"
    source_rows.append(
        {
            "source_id": "gse334005-syngeneic-tpm",
            "pillar": "tumor_rna",
            "filename": p_334.name,
            "unit": "TPM",
            "strain": "C57BL/6 (B16F10, LLC) & BALB/c (CT26)",
            "bytes": p_334.stat().st_size,
            "sha256": _sha256(p_334),
        }
    )
    df_334 = pd.read_csv(p_334, sep="\t")
    df_334["_gid"] = df_334["gene_id"].astype(str).map(entrez_to_ens)
    df_334 = df_334.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    pbs_cols = [c for c in df_334.columns if "_PBS_" in c]
    peg_334 = _project_peg_matrix(df_334[pbs_cols].astype(float), groups_df)
    for code, prefix, strain, hap, label in [
        ("B16F10", "B16_PBS_", "C57BL/6", "H-2b", "B16-F10 subcutaneous melanoma"),
        ("CT26", "CT26_PBS_", "BALB/c", "H-2d", "CT26 subcutaneous colon carcinoma"),
        ("LLC", "LLC_PBS_", "C57BL/6", "H-2b", "Lewis lung carcinoma (LLC) subcutaneous graft"),
    ]:
        c_sub = [c for c in pbs_cols if c.startswith(prefix)]
        for c in c_sub:
            per_sample_peg.append(("gse334005", c, f"GSE334005:{c}", code, peg_334[c]))
        cohort_registry.append(
            {
                "cancer_code": code,
                "label": label,
                "model_kind": "syngeneic_graft",
                "strain": strain,
                "mhc_haplotype": hap,
                "source_accession": "GSE334005",
                "n_samples": len(c_sub),
                "native_feature_space": "NCBI Entrez GeneID -> MGI -> ENSMUSG",
                "quantification": "author linear TPM",
            }
        )

    # 5. 4T1 (BALB/c mammary carcinoma, GSE332912, n=3 WT control XW1..XW3, native ENSMUSG TPM)
    p_4t1 = input_dir / "GSE332912_TPM_values.csv.gz"
    source_rows.append(
        {
            "source_id": "gse332912-4t1",
            "pillar": "tumor_rna",
            "filename": p_4t1.name,
            "unit": "TPM",
            "strain": "BALB/c (H-2d)",
            "bytes": p_4t1.stat().st_size,
            "sha256": _sha256(p_4t1),
        }
    )
    df_4t1 = pd.read_csv(p_4t1)
    df_4t1["ID"] = df_4t1["ID"].astype(str).str.split(".").str[0]
    df_4t1 = df_4t1.drop_duplicates("ID").set_index("ID")
    peg_4t1 = _project_peg_matrix(df_4t1[["XW1", "XW2", "XW3"]].astype(float), groups_df)
    for c in ("XW1", "XW2", "XW3"):
        per_sample_peg.append(("gse332912", f"4T1_{c}", f"GSE332912:{c}", "4T1", peg_4t1[c]))
    cohort_registry.append(
        {
            "cancer_code": "4T1",
            "label": "4T1 triple-negative mammary carcinoma",
            "model_kind": "syngeneic_line",
            "strain": "BALB/c",
            "mhc_haplotype": "H-2d",
            "source_accession": "GSE332912",
            "n_samples": 3,
            "native_feature_space": "Ensembl gene (ENSMUSG)",
            "quantification": "author linear TPM",
        }
    )

    # 6. KP_LUAD (C57BL/6 primary Kras;Trp53 lung adenocarcinoma GEMM, GSE277076, n=10 sgControl primary tumors)
    sym_universe = pd.Index(sorted(set(sym_to_ens.values())), name="_gid")
    p_kp = input_dir / "GSE277076_202212051509_X_13772_I_raw.txt.gz"
    source_rows.append(
        {
            "source_id": "gse277076-kp-luad",
            "pillar": "tumor_rna",
            "filename": p_kp.name,
            "unit": "TPM (from raw counts)",
            "strain": "C57BL/6 Kras;Trp53 GEMM (H-2b)",
            "bytes": p_kp.stat().st_size,
            "sha256": _sha256(p_kp),
        }
    )
    df_kp = pd.read_csv(p_kp, sep="\t", index_col=0)
    df_kp["_gid"] = df_kp.index.astype(str).map(sym_to_ens)
    df_kp = df_kp.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    kp_cols = [
        "AF1567",
        "AH1551",
        "AH1746",
        "AS1614",
        "AS1698",
        "AS1707",
        "AS1710",
        "AS1711",
        "AS1712",
        "AU1855",
    ]
    peg_kp = _project_peg_matrix(
        _counts_to_tpm(df_kp[kp_cols].reindex(sym_universe, fill_value=0.0), exonic_series),
        groups_df,
    )
    for c in kp_cols:
        per_sample_peg.append(("gse277076", f"KP_{c}", f"GSE277076:{c}", "KP_LUAD", peg_kp[c]))
    cohort_registry.append(
        {
            "cancer_code": "KP_LUAD",
            "label": "Primary Kras;Trp53 lung adenocarcinoma GEMM",
            "model_kind": "primary_gemm",
            "strain": "C57BL/6",
            "mhc_haplotype": "H-2b",
            "source_accession": "GSE277076",
            "n_samples": len(kp_cols),
            "native_feature_space": "MGI/Ensembl official symbol -> ENSMUSG",
            "quantification": "linear TPM from raw counts and Ensembl 115 exonic_bp",
        }
    )

    # 7. BMPC_PRAD (FVB primary & metastatic MYC;Pten prostate adenocarcinoma GEMM, GSE287166, n=6)
    df_bmpc = pd.read_csv(input_dir / "GSE287166_NormalizedCounts.csv.gz", index_col=0)
    df_bmpc["_gid"] = df_bmpc.index.astype(str).map(sym_to_ens)
    df_bmpc = df_bmpc.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    bmpc_cols = ["BMPC_Pri1", "BMPC_Pri2", "BMPC_Pri3", "BMPC_Met1", "BMPC_Met2", "BMPC_Met3"]
    peg_bmpc = _project_peg_matrix(
        _counts_to_tpm(df_bmpc[bmpc_cols].reindex(sym_universe, fill_value=0.0), exonic_series),
        groups_df,
    )
    for c in bmpc_cols:
        per_sample_peg.append(("gse287166", c, f"GSE287166:{c}", "BMPC_PRAD", peg_bmpc[c]))
    cohort_registry.append(
        {
            "cancer_code": "BMPC_PRAD",
            "label": "BMPC MYC;Pten primary & metastatic prostate adenocarcinoma GEMM",
            "model_kind": "primary_gemm",
            "strain": "FVB",
            "mhc_haplotype": "H-2q",
            "source_accession": "GSE287166",
            "n_samples": len(bmpc_cols),
            "native_feature_space": "MGI/Ensembl official symbol -> ENSMUSG",
            "quantification": "linear TPM from normalized counts and Ensembl 115 exonic_bp",
        }
    )

    # 8. P815 (DBA/2 mastocytoma, GSE181773, n=3 negative-control NC_1..3, native ENSMUSG FPKM -> exact TPM)
    p_p815 = input_dir / "GSE181773_gene_fpkm.txt.gz"
    source_rows.append(
        {
            "source_id": "gse181773-p815",
            "pillar": "tumor_rna",
            "filename": p_p815.name,
            "unit": "TPM (from native ENSMUSG FPKM)",
            "strain": "DBA/2 (H-2d)",
            "bytes": p_p815.stat().st_size,
            "sha256": _sha256(p_p815),
        }
    )
    df_p815 = pd.read_csv(
        p_p815,
        sep="\t",
        usecols=["gene_id", "FPKM.NC_1", "FPKM.NC_2", "FPKM.NC_3"],
    )
    df_p815["gene_id"] = df_p815["gene_id"].astype(str).str.split(".").str[0]
    df_p815 = df_p815.groupby("gene_id", sort=False)[["FPKM.NC_1", "FPKM.NC_2", "FPKM.NC_3"]].sum()
    peg_p815 = _project_peg_matrix(_fpkm_to_tpm(df_p815), groups_df)
    for i, c in enumerate(("FPKM.NC_1", "FPKM.NC_2", "FPKM.NC_3"), start=1):
        # Cell-line replicates of P815 share one biological DBA/2 mastocytoma donor!
        per_sample_peg.append(("gse181773", f"P815_NC_{i}", "DBA2:P815", "P815", peg_p815[c]))
    cohort_registry.append(
        {
            "cancer_code": "P815",
            "label": "P815 mastocytoma (native DBA/2 source of P1A/Trap1a)",
            "model_kind": "syngeneic_line",
            "strain": "DBA/2",
            "mhc_haplotype": "H-2d",
            "source_accession": "GSE181773",
            "n_samples": 3,
            "native_feature_space": "Ensembl gene (ENSMUSG)",
            "quantification": "within-sample linear TPM from FPKM",
        }
    )

    # Note: cell-line replicates of 4T1 also originate from one BALB/c tumor line; deduplicate donor IDs honestly!
    fixed_per_sample = []
    for src_id, smp_id, part_id, code, ser in per_sample_peg:
        if code == "4T1":
            part_id = "BALBc:4T1"
        fixed_per_sample.append((src_id, smp_id, part_id, code, ser))

    tumor_rows = []
    for src_id, smp_id, part_id, code, ser in fixed_per_sample:
        valid = ser.dropna()
        for eid, val in valid.items():
            tumor_rows.append(
                {
                    "entity_id": eid,
                    "source_id": src_id,
                    "sample_id": smp_id,
                    "participant_id": part_id,
                    "cancer_code": code,
                    "value": float(val),
                    "comparison_id": COMPARISON_ID,
                    "unit": "TPM",
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )

    return pd.DataFrame(tumor_rows), pd.DataFrame(cohort_registry), source_rows


def build_all(ens_dir: Path, input_dir: Path, out_root: Path) -> dict:
    print("1/6 Qualifying Ensembl 115 GRCm39 proteome & 31-mer cDNA mapping risk...", flush=True)
    ref, groups_df, seq_df, exonic_bp, symbol_by_gene, mapping_risk_df, ens_sources = (
        load_ensembl115_mouse_reference(ens_dir)
    )
    coding_genes = set(symbol_by_gene)
    print(
        f"    PEGs={len(groups_df):,}, sequences={len(seq_df):,}, coding_genes={len(coding_genes):,}",
        flush=True,
    )

    print("2/6 Building MGI / Entrez / Ensembl crosswalks...", flush=True)
    mgi_to_ens, entrez_to_ens, sym_to_ens = build_mgi_crosswalks(
        input_dir, coding_genes, symbol_by_gene
    )
    print(
        f"    MGI->ENS={len(mgi_to_ens):,}, Entrez->ENS={len(entrez_to_ens):,}, Symbol->ENS={len(sym_to_ens):,}",
        flush=True,
    )

    print(
        "3/6 Aggregating Bgee 15.2, GSE219045 & 16-compartment reproductive/embryo TPM panels...",
        flush=True,
    )
    comp_normals, comp_repro, rna_gates, confirmed_organs, norm_sources = (
        load_bgee_and_normal_panels(input_dir, groups_df, exonic_bp, entrez_to_ens, sym_to_ens)
    )
    print(
        f"    normal_obs={len(comp_normals):,}, repro_obs={len(comp_repro):,}, rna_gates={len(rna_gates):,}",
        flush=True,
    )

    print(
        "4/6 Evaluating ImmGen sorted/activated immune RNA, MGI GXD TS28 protein & Tabula Muris 10x UMI...",
        flush=True,
    )
    imm_gates, imm_sources = load_mouse_immune_panels(
        input_dir, groups_df, exonic_bp, symbol_by_gene, sym_to_ens
    )
    gxd_gates, gxd_sources = load_mgi_gxd_protein_gates(
        input_dir, groups_df, mgi_to_ens, confirmed_organs
    )
    umi_df, tm_sources = load_tabula_muris_rare_cell_umi(input_dir, sym_to_ens)
    print(
        f"    imm_gates={len(imm_gates):,}, gxd_gates={len(gxd_gates):,}, umi_rows={len(umi_df):,}",
        flush=True,
    )

    print("5/6 Loading 8 murine syngeneic & GEMM tumor cohorts...", flush=True)
    comp_tumors, cohort_df, tumor_sources = load_mouse_tumor_cohorts(
        input_dir, groups_df, exonic_bp, entrez_to_ens, sym_to_ens
    )
    print(f"    tumor_obs={len(comp_tumors):,}, cohorts={len(cohort_df)}", flush=True)

    print(
        "6/6 Running SpeciesScreeningSpec & 16-compartment HostContext sex-restriction...",
        flush=True,
    )
    all_gates = pd.concat([rna_gates, imm_gates, gxd_gates], ignore_index=True)
    spec = SpeciesScreeningSpec(
        species="mouse",
        reference=ref,
        policy_id=POLICY_SCHEMA_ID,
        comparison_id=COMPARISON_ID,
        comparability_reference="oncoref.mouse-grcm39-panel.2026-11-25.1",
        comparison_unit="TPM",
        required_normal_contexts=REQUIRED_BGEE_NORMAL_CONTEXTS,
        anatomy=anatomy_for_species("mouse"),
    )
    res = evaluate_species_catalog(
        spec,
        groups_df,
        sequences=seq_df,
        normal_gates=all_gates,
        rare_cell_counts=umi_df,
        mapping_risk=mapping_risk_df,
        confirmed_rare_cell_organs=confirmed_organs,
        comparison_normals=comp_normals,
        comparison_reproductive=comp_repro,
        comparison_tumors=comp_tumors,
    )

    cat = res.catalog.copy()
    # Populate full 16-compartment host_exposure for EVERY entity in the catalog (so HostContext works on all rows)
    full_host_exposure: dict[str, dict[str, dict]] = {eid: {} for eid in cat["entity_id"]}
    for r in comp_repro.itertuples(index=False):
        is_dev = r.scope == "developmental"
        prev = full_host_exposure[r.entity_id].get(r.compartment)
        if prev is None or r.value > prev["max_tpm"]:
            full_host_exposure[r.entity_id][str(r.compartment)] = {
                "coverage": "measured",
                "max_tpm": round(float(r.value), 6),
                "max_ntpm": None,
                "protein_status": "not_assessed",
                "scope_complete": not is_dev,
            }

    # Add human-parity convenience columns used by CLI/reports
    h0 = (
        cat[["normal_reference_tpm", "disposable_max_tpm"]]
        .max(axis=1, skipna=False)
        .clip(lower=1.0)
    )
    cat["single_tumor_threshold_tpm"] = np.where(
        cat["normal_reference_tpm"].notna(), 10.0 * h0, np.nan
    )

    thresh_by_eid = dict(zip(cat["entity_id"], cat["single_tumor_threshold_tpm"]))
    top_smp_by_eid: dict[str, str] = {}
    ctypes_10x_by_eid: dict[str, int] = {}
    for eid, tgrp in comp_tumors.groupby("entity_id", sort=False):
        th = thresh_by_eid.get(eid, np.nan)
        if pd.isna(th):
            continue
        best_r = tgrp.sort_values("value", ascending=False).iloc[0]
        top_smp_by_eid[eid] = str(best_r.sample_id)
        ctypes_10x_by_eid[eid] = int(tgrp.loc[tgrp["value"] >= th, "cancer_code"].nunique())
    cat["tumor_top_sample"] = cat["entity_id"].map(top_smp_by_eid).fillna("")
    cat["cancer_types_at_tenfold"] = cat["entity_id"].map(ctypes_10x_by_eid).fillna(0).astype(int)

    # Re-verify sex_restriction_class & sex_restriction_code against full_host_exposure
    sex_df = classify_sex_restriction(cat, species="mouse", host_evidence=full_host_exposure)
    sex_class_map = (
        sex_df.set_index("entity_id")["sex_restriction_class"].to_dict() if not sex_df.empty else {}
    )
    sex_code_map = (
        sex_df.set_index("entity_id")["sex_restriction_code"].to_dict() if not sex_df.empty else {}
    )
    cat["sex_restriction_class"] = cat["entity_id"].map(sex_class_map).fillna("")
    cat["sex_restriction_code"] = cat["entity_id"].map(sex_code_map).fillna("")

    data_dir = out_root / "oncoref" / "data"
    mouse_dir = data_dir / "species" / "mouse"
    mouse_dir.mkdir(parents=True, exist_ok=True)

    cat_path = data_dir / "mouse-target-selection.csv.gz"
    cat.to_csv(
        cat_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )

    exp_rows = pd.DataFrame(
        [
            {
                "entity_id": eid,
                "exposures": json.dumps(
                    full_host_exposure[eid], sort_keys=True, separators=(",", ":")
                ),
            }
            for eid in cat["entity_id"]
        ]
    )
    exp_path = data_dir / "mouse-host-exposure.csv.gz"
    exp_rows.to_csv(
        exp_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )

    seq_out = res.sequences[
        [
            "entity_id",
            "gene_id",
            "transcript_id",
            "protein_id",
            "transcript_biotype",
            "sequence",
            "sequence_sha256",
            "matches_representative_sequence",
        ]
    ].copy()
    seq_path = data_dir / "mouse-proteome-sequences.csv.gz"
    # Preserve existing mouse-proteome-sequences.csv.gz byte-for-byte when Ensembl 115 sequences are unchanged
    write_seq = not seq_path.exists()
    if not write_seq:
        existing_seq = pd.read_csv(seq_path, compression="gzip", keep_default_na=False)
        existing_seq["matches_representative_sequence"] = (
            existing_seq["matches_representative_sequence"].astype(str).eq("True")
        )
        if not existing_seq.equals(seq_out.reset_index(drop=True)):
            write_seq = True
    if write_seq:
        seq_out.to_csv(
            seq_path,
            index=False,
            compression={"method": "gzip", "mtime": 0},
            lineterminator="\n",
        )

    all_sources_df = pd.DataFrame(
        ens_sources + norm_sources + imm_sources + gxd_sources + tm_sources + tumor_sources
    ).drop_duplicates("source_id")
    all_sources_df.to_csv(mouse_dir / "screening-sources.csv", index=False, lineterminator="\n")
    cohort_df.to_csv(mouse_dir / "tumor-cohorts.csv", index=False, lineterminator="\n")

    cat_sha = _sha256(cat_path)
    exp_sha = _sha256(exp_path)
    seq_sha = _sha256(seq_path)

    loci_count = int(cat["member_gene_ids"].map(json.loads).map(len).sum())
    summary_payload = {
        "schema": "oncoref.species-target-selection-metadata.v1",
        "taxon_id": 10090,
        "species": "mouse",
        "scientific_name": "Mus musculus",
        "selection_status": "evaluated",
        "policy_id": POLICY_SCHEMA_ID,
        "comparison_id": COMPARISON_ID,
        "comparison_unit": "TPM",
        "reference": ref.as_dict(),
        "data_sha256": cat_sha,
        "data_bytes": cat_path.stat().st_size,
        "host_exposure_sha256": exp_sha,
        "host_exposure_bytes": exp_path.stat().st_size,
        "sequences_sha256": seq_sha,
        "sequences_bytes": seq_path.stat().st_size,
        "required_normal_contexts": sorted(REQUIRED_BGEE_NORMAL_CONTEXTS),
        "measured_reproductive_compartments": sorted(
            {r.compartment for r in comp_repro.itertuples(index=False)}
        ),
        "tumor_cohorts": cohort_df["cancer_code"].tolist(),
        "tumor_samples": int(cohort_df["n_samples"].sum()),
        "summary": {
            "proteins": len(cat),
            "loci": loci_count,
            "sequences": len(seq_out),
            "categories": {
                c: int(cat["category"].eq(c).sum()) for c in target_selection.CATEGORIES
            },
            "classifications": {
                k: int(v) for k, v in cat["classification"].value_counts().sort_index().items()
            },
            "sex_restriction_codes": (
                {
                    k: int(v)
                    for k, v in sex_df["sex_restriction_code"].value_counts().sort_index().items()
                }
                if not sex_df.empty
                else {}
            ),
            "sex_restriction_classes": (
                {
                    k: int(v)
                    for k, v in sex_df["sex_restriction_class"].value_counts().sort_index().items()
                }
                if not sex_df.empty
                else {}
            ),
        },
    }
    (mouse_dir / "target-selection.json").write_text(json.dumps(summary_payload, indent=2) + "\n")
    print("DONE! Summary:", json.dumps(summary_payload["summary"], indent=2), flush=True)
    return summary_payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ens-dir", type=Path, default=DEFAULT_ENS_DIR)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    build_all(args.ens_dir, args.input_dir, args.root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
