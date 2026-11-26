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
import os
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
    os.environ.get("ONCOREF_MOUSE_ENS115_DIR", str(ROOT / "tmp" / "ens115_mouse"))
)
DEFAULT_INPUT_DIR = Path(
    os.environ.get("ONCOREF_MOUSE_INPUTS_DIR", str(ROOT / "tmp" / "mouse_catalog_inputs"))
)
DEFAULT_EXTRA_DIR = Path(
    os.environ.get("ONCOREF_MOUSE_EXTRA_DIR", str(DEFAULT_INPUT_DIR.parent / "mouse_153_extra"))
)

PRIMARY_CONTIGS = frozenset({*(str(i) for i in range(1, 20)), "X", "Y", "MT"})
COMPARISON_ID = "oncoref.mouse-grcm39-tpm.2026-11-26.2"

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
DENTAL_HOME_SYMBOLS = frozenset(
    {"Amelx", "Ambn", "Enam", "Dspp", "Mmp20", "Klk4", "Odaph", "Odam", "Amtn"}
)
TASTE_HOME_SYMBOLS = frozenset({"Otop1"})


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
            elif sym.startswith(("Tas1r", "Tas2r")) or sym in TASTE_HOME_SYMBOLS:
                homes.append("taste bud")
            elif sym.startswith(("Olfr", "Vmn1r", "Vmn2r")) or re.match(
                r"^(?:Or\d+[a-z]+\d+|Taar\d+[a-z]*)$", sym, re.IGNORECASE
            ):
                homes.append("olfactory epithelium")
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

    # Two-pass phase-invariant 31-mer paralog check across all Ensembl 115 cDNA + ncRNA transcripts:
    # 1) Sample probe 31-mers (step=4) across protein-coding transcripts per PEG.
    # 2) Scan EVERY transcript in cdna.all.fa.gz + ncrna.fa.gz at step=1 (every base offset)
    #    so out-of-phase 5'-UTR shifts and processed pseudogenes/ncRNAs are caught identically.
    cdna_seqs_by_gene: dict[str, list[bytes]] = defaultdict(list)
    c_hdr = ""
    c_parts: list[str] = []

    def _extract_gid(hdr: str) -> str:
        for tok in hdr[1:].split():
            if tok.startswith("gene:"):
                return tok.split(":", 1)[1].split(".")[0]
        return ""

    def _flush_coding_cdna(hdr: str, parts_l: list[str]) -> None:
        if not hdr:
            return
        gid = _extract_gid(hdr)
        if gid in gene_to_entity:
            cdna_seqs_by_gene[gid].append("".join(parts_l).encode("ascii", errors="ignore"))

    with gzip.open(role_files["rna"], "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush_coding_cdna(c_hdr, c_parts)
                c_hdr = line.rstrip("\n")
                c_parts = []
            else:
                c_parts.append(line.strip())
    _flush_coding_cdna(c_hdr, c_parts)

    k = 31
    probe_step = 4
    MULTI = "__MULTI__"
    probe_by_gene: dict[str, set[bytes]] = {}
    for gid, seqs in cdna_seqs_by_gene.items():
        eid = gene_to_entity[gid]
        seen_for_gene = {
            s[i : i + k] for s in seqs if len(s) >= k for i in range(0, len(s) - k + 1, probe_step)
        }
        probe_by_gene[gid] = seen_for_gene
        for km_b in seen_for_gene:
            prev = kmer_owners.get(km_b)
            if prev is None:
                kmer_owners[km_b] = eid
            elif prev != eid and prev != MULTI:
                kmer_owners[km_b] = MULTI

    def _scan_tx_step1(hdr: str, parts_l: list[str]) -> None:
        if not hdr:
            return
        gid = _extract_gid(hdr)
        if not gid:
            return
        owner = gene_to_entity.get(gid, f"NC:{gid}")
        s_b = "".join(parts_l).encode("ascii", errors="ignore")
        if len(s_b) < k:
            return
        for i in range(len(s_b) - k + 1):
            km_b = s_b[i : i + k]
            prev = kmer_owners.get(km_b)
            if prev is not None and prev != owner and prev != MULTI:
                kmer_owners[km_b] = MULTI

    for role_key in ("rna", "ncrna"):
        tx_path = role_files.get(role_key)
        if tx_path is None or not tx_path.exists():
            continue
        t_hdr = ""
        t_parts: list[str] = []
        with gzip.open(tx_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    _scan_tx_step1(t_hdr, t_parts)
                    t_hdr = line.rstrip("\n")
                    t_parts = []
                else:
                    t_parts.append(line.strip())
        _scan_tx_step1(t_hdr, t_parts)

    for gid, seen_for_gene in probe_by_gene.items():
        total_kmer_counts[gid] = len(seen_for_gene)
        ambig_kmer_counts[gid] = sum(1 for km_b in seen_for_gene if kmer_owners.get(km_b) == MULTI)

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


def _mouse_olfr_gene_ids(symbol_by_gene: dict[str, str]) -> set[str]:
    return {
        gid
        for gid, sym in symbol_by_gene.items()
        if sym.startswith("Olfr") or re.match(r"^Or\d+[a-z]+\d+", sym, re.IGNORECASE)
    }


def _olfr_share_ge5(gene_tpm: pd.DataFrame, olfr_gids: set[str]) -> pd.Series:
    olfr_sub = gene_tpm.loc[gene_tpm.index.isin(olfr_gids)]
    if len(olfr_sub) < 200:
        return pd.Series(0.0, index=gene_tpm.columns, dtype=float)
    return (olfr_sub >= 5.0).mean(axis=0)


def load_bgee_and_normal_panels(
    input_dir: Path,
    extra_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    symbol_by_gene: dict[str, str],
    entrez_to_ens: dict[str, str],
    sym_to_ens: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, set[str]], list[dict]]:
    """Aggregate normal somatic N, thymic, and all 16 reproductive/developmental compartments."""
    olfr_gids = _mouse_olfr_gene_ids(symbol_by_gene)
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
        olfr_sh = _olfr_share_ge5(piv, olfr_gids)
        qc_keep_cols = [
            lid
            for lid, (_sc, t_ctx, _st) in keep_libs.items()
            if lid in piv.columns
            and (
                t_ctx in {"olfactory epithelium", "testis"} or float(olfr_sh.get(lid, 0.0)) <= 0.025
            )
        ]
        if not qc_keep_cols:
            continue
        peg_mat = _project_peg_matrix(piv[qc_keep_cols], groups_df)

        for lid in qc_keep_cols:
            scope, target_ctx, st = keep_libs[lid]
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
    g219_olfr = _olfr_share_ge5(g219_tpm, olfr_gids)
    g219_peg = _project_peg_matrix(g219_tpm, groups_df)

    for col_name in g219_peg.columns:
        suffix = col_name.removeprefix("Mouse_Female_").removeprefix("Mouse_Male_")
        sex_tag = "F" if col_name.startswith("Mouse_Female_") else "M"
        mapped = GSE219045_COL_MAP.get(suffix)
        if mapped is None:
            continue
        scope, target_ctx = mapped
        if (
            target_ctx not in {"olfactory epithelium", "testis"}
            and float(g219_olfr.get(col_name, 0.0)) > 0.025
        ):
            continue
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

    # 2b. GSE322844 untreated WT murine liver controls (Control1..5, ENSMUSG FPKM -> exact linear TPM)
    p322_ctrl = extra_dir / "GSE322844_Control_FPKM.csv.gz"
    if p322_ctrl.is_file():
        source_rows.append(
            {
                "source_id": "gse322844-wt-liver",
                "pillar": "healthy_baseline_rna",
                "filename": p322_ctrl.name,
                "unit": "TPM (from native ENSMUSG FPKM)",
                "strain": "untreated wild-type C57BL/6 liver (n=5)",
                "bytes": p322_ctrl.stat().st_size,
                "sha256": _sha256(p322_ctrl),
            }
        )
        df322_c = pd.read_csv(p322_ctrl, sep=";")
        df322_c["gene_id"] = df322_c["gene_id"].astype(str).str.split(".").str[0]
        ctrl_cols_322 = [c for c in df322_c.columns if c.startswith("Control")]
        df322_c = df322_c.groupby("gene_id", sort=False)[ctrl_cols_322].sum()
        tpm322_c = _fpkm_to_tpm(df322_c)
        olfr322_c = _olfr_share_ge5(tpm322_c, olfr_gids)
        qc322_c = [c for c in ctrl_cols_322 if float(olfr322_c.get(c, 0.0)) <= 0.025]
        if qc322_c:
            peg322_c = _project_peg_matrix(tpm322_c[qc322_c], groups_df)
            for c in qc322_c:
                somatic_peg_cols["liver"][f"gse322844:{c}"] = peg322_c[c]

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
    extra_dir: Path,
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

    if extra_dir.is_dir():
        # 4. GSE287306 WT splenic naive & IL-15-stimulated NK cells (ENSMUSG expected counts -> TPM)
        p287_nk = extra_dir / "GSE287306_Khanna_-_Rajiv_RNASeq.genes.ExpectedCounts.csv.gz"
        if p287_nk.is_file():
            source_rows.append(
                {
                    "source_id": "gse287306-splenic-il15-nk",
                    "pillar": "normal_bulk_and_immune_safety",
                    "filename": p287_nk.name,
                    "unit": "TPM (from native ENSMUSG expected counts)",
                    "strain": "C57BL/6 WT splenic naive & IL-15-treated NK cells",
                    "bytes": p287_nk.stat().st_size,
                    "sha256": _sha256(p287_nk),
                }
            )
            df287_nk = pd.read_csv(p287_nk)
            df287_nk["Ensembl Gene ID"] = (
                df287_nk["Ensembl Gene ID"].astype(str).str.split(".").str[0]
            )
            nk_cols = [c for c in df287_nk.columns if c.startswith(("naiveNK", "IL15NK"))]
            df287_nk = df287_nk.drop_duplicates("Ensembl Gene ID").set_index("Ensembl Gene ID")
            tpm287_nk = _counts_to_tpm(df287_nk[nk_cols], exonic_series)
            olfr287_nk = _olfr_share_ge5(tpm287_nk, olfr_gids)
            qc287_nk = [c for c in nk_cols if float(olfr287_nk.get(c, 0.0)) <= 0.025]
            if qc287_nk:
                peg287_nk = _project_peg_matrix(tpm287_nk[qc287_nk], groups_df)
                for pfx_nk in ("naiveNK", "IL15NK"):
                    cs_nk = [c for c in qc287_nk if c.startswith(pfx_nk)]
                    if len(cs_nk) >= 2:
                        pop_summaries[f"GSE287306:{pfx_nk}"] = (
                            "protected",
                            _replicated_p95(peg287_nk[cs_nk]),
                        )

        # 5. GSE336796 WT splenic & vehicle/ionomycin-stimulated NK cells (ENSMUSG linear TPM)
        p336_sp = extra_dir / "GSE336796_tpm_table_with_symbol.csv.gz"
        if p336_sp.is_file():
            source_rows.append(
                {
                    "source_id": "gse336796-splenic-nk",
                    "pillar": "normal_bulk_and_immune_safety",
                    "filename": p336_sp.name,
                    "unit": "TPM (native ENSMUSG)",
                    "strain": "C57BL/6 WT splenic NK subsets",
                    "bytes": p336_sp.stat().st_size,
                    "sha256": _sha256(p336_sp),
                }
            )
            df336_sp = pd.read_csv(p336_sp)
            df336_sp["X"] = df336_sp["X"].astype(str).str.split(".").str[0]
            wt_sp_cols = ["WT1_SP", "WT2_SP", "WT3_SP"]
            df336_sp = df336_sp.drop_duplicates("X").set_index("X")
            tpm336_sp = df336_sp[wt_sp_cols].astype(float)
            olfr336_sp = _olfr_share_ge5(tpm336_sp, olfr_gids)
            qc336_sp = [c for c in wt_sp_cols if float(olfr336_sp.get(c, 0.0)) <= 0.025]
            if len(qc336_sp) >= 2:
                peg336_sp = _project_peg_matrix(tpm336_sp[qc336_sp], groups_df)
                pop_summaries["GSE336796:WT_SP_NK"] = (
                    "protected",
                    _replicated_p95(peg336_sp[qc336_sp]),
                )

        p336_io = extra_dir / "GSE336796_vehicle_vs_iono_TPM.csv.gz"
        if p336_io.is_file():
            source_rows.append(
                {
                    "source_id": "gse336796-iono-nk",
                    "pillar": "normal_bulk_and_immune_safety",
                    "filename": p336_io.name,
                    "unit": "TPM (native ENSMUSG)",
                    "strain": "C57BL/6 WT vehicle & ionomycin-activated NK cells",
                    "bytes": p336_io.stat().st_size,
                    "sha256": _sha256(p336_io),
                }
            )
            df336_io = pd.read_csv(p336_io)
            df336_io["gene_id"] = df336_io["gene_id"].astype(str).str.split(".").str[0]
            io_cols = ["V1", "V2", "V3", "I1", "I2", "I3"]
            df336_io = df336_io.drop_duplicates("gene_id").set_index("gene_id")
            tpm336_io = df336_io[io_cols].astype(float)
            olfr336_io = _olfr_share_ge5(tpm336_io, olfr_gids)
            qc336_io = [c for c in io_cols if float(olfr336_io.get(c, 0.0)) <= 0.025]
            peg336_io = _project_peg_matrix(tpm336_io[qc336_io], groups_df)
            for lbl_io, cs_io in (
                ("WT_Vehicle_NK", [c for c in qc336_io if c.startswith("V")]),
                ("WT_Iono_NK", [c for c in qc336_io if c.startswith("I")]),
            ):
                if len(cs_io) >= 2:
                    pop_summaries[f"GSE336796:{lbl_io}"] = (
                        "protected",
                        _replicated_p95(peg336_io[cs_io]),
                    )

        # 6. GSE302070 WT splenic PBS & IL-2-IC-stimulated Tregs, NK cells & B cells
        p302 = extra_dir / "GSE302070_Genes_count_table.tsv.gz"
        if p302.is_file():
            source_rows.append(
                {
                    "source_id": "gse302070-il2-splenic-subsets",
                    "pillar": "normal_bulk_and_immune_safety",
                    "filename": p302.name,
                    "unit": "TPM (from sorted splenic raw counts)",
                    "strain": "Foxp3-GFP C57BL/6 splenic Tregs, NK & B cells (PBS & IL-2-IC)",
                    "bytes": p302.stat().st_size,
                    "sha256": _sha256(p302),
                }
            )
            df302 = pd.read_csv(p302, sep="\t")
            gid302 = df302["GeneSymbol"].astype(str).map(sym_to_ens)
            cols302 = [
                c
                for c in df302.columns
                if c != "GeneSymbol"
                and c.startswith(
                    (
                        "24hPBSTreg",
                        "24hPBSNK",
                        "24hrPBSBcell",
                        "24hrHF10Treg",
                        "24hrHF10NK",
                        "24hrHF10Bcell",
                    )
                )
            ]
            sub302 = df302.loc[gid302.notna(), cols302].copy()
            sub302.index = gid302[gid302.notna()].to_numpy()
            sub302 = sub302[~sub302.index.duplicated(keep="first")].reindex(
                sym_universe, fill_value=0.0
            )
            tpm302 = _counts_to_tpm(sub302, exonic_series)
            olfr302 = _olfr_share_ge5(tpm302, olfr_gids)
            qc302 = [c for c in cols302 if float(olfr302.get(c, 0.0)) <= 0.025]
            if qc302:
                peg302 = _project_peg_matrix(tpm302[qc302], groups_df)
                by_p302: dict[str, list[str]] = defaultdict(list)
                for c in qc302:
                    by_p302[c.rsplit("_", 1)[0]].append(c)
                for pop_k, cs_k in sorted(by_p302.items()):
                    if len(cs_k) >= 2:
                        tier_k = "manageable" if "Bcell" in pop_k else "protected"
                        pop_summaries[f"GSE302070:{pop_k}"] = (
                            tier_k,
                            _replicated_p95(peg302[cs_k]),
                        )

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

    # Exclude IL-2/IL-12/IL-18-inducible murine LAK & uterine/decidual NK cytolytic effector granzymes
    # (Jenne et al. 1988/1989; Revell et al. 2005; verified 0-2 TPM in pure cancer lines & 1000x induced with Prf1/Gzmb in GSE334005)
    lak_granzymes = {"Gzmd", "Gzme", "Gzmf", "Gzmg"}
    for r_g in groups_df.itertuples(index=False):
        lbl_parts = {p.strip() for p in str(r_g.protein_group_label).split("/") if p.strip()}
        hit_gzm = sorted(lbl_parts & lak_granzymes)
        if hit_gzm:
            gate_rows.append(
                {
                    "entity_id": str(r_g.entity_id),
                    "layer": "Curated murine activated NK/LAK cytolytic effector safety",
                    "unsafe": True,
                    "hold": False,
                    "reason": f"Essential IL-2/IL-12/IL-18 LAK & NK cytolytic granzyme ({', '.join(hit_gzm)})",
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
    extra_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    symbol_by_gene: dict[str, str],
    entrez_to_ens: dict[str, str],
    sym_to_ens: dict[str, str],
    risk_gids: set[str],
    unsafe_eids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    """Load 14 strain-explicit in-vivo murine tumor cohorts + 9 pure in-vitro cancer cell lines in linear TPM."""
    exonic_series = pd.Series(exonic_bp, dtype=float)
    sym_universe = pd.Index(sorted(set(sym_to_ens.values())), name="_gid")
    olfr_gids = _mouse_olfr_gene_ids(symbol_by_gene)
    cohort_registry: list[dict] = []
    cell_line_registry: list[dict] = []
    source_rows: list[dict] = []
    cohort_matrices: list[tuple[str, str, list[str], list[str], pd.DataFrame]] = []
    cell_line_matrices: list[tuple[str, str, str, list[str], pd.DataFrame]] = []

    def _qc_cols(tpm_frame: pd.DataFrame, cols: list[str]) -> list[str]:
        sh = _olfr_share_ge5(tpm_frame[cols], olfr_gids)
        return [c for c in cols if float(sh.get(c, 0.0)) <= 0.025]

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
    tpm_mc38 = _counts_to_tpm(
        df_mc38[["A1", "A2", "A3", "A4"]], df_mc38["gene_length"].astype(float)
    )
    mc38_cols = _qc_cols(tpm_mc38, ["A1", "A2", "A3", "A4"])
    peg_mc38 = _project_peg_matrix(tpm_mc38[mc38_cols], groups_df)
    mc38_sids = [f"MC38_{c}" for c in mc38_cols]
    peg_mc38.columns = mc38_sids
    cohort_matrices.append(
        ("gse355313", "MC38", mc38_sids, [f"GSE355313:{c}" for c in mc38_cols], peg_mc38)
    )
    cohort_registry.append(
        {
            "cancer_code": "MC38",
            "label": "MC38 subcutaneous colon adenocarcinoma",
            "model_kind": "syngeneic_graft",
            "strain": "C57BL/6",
            "mhc_haplotype": "H-2b",
            "source_accession": "GSE355313",
            "n_samples": len(mc38_cols),
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
    pbs_raw_cols = [c for c in df_334.columns if "_PBS_" in c]
    tpm_334 = df_334[pbs_raw_cols].astype(float)
    pbs_cols = _qc_cols(tpm_334, pbs_raw_cols)
    peg_334 = _project_peg_matrix(tpm_334[pbs_cols], groups_df)
    for code, prefix, strain, hap, label in [
        ("B16F10", "B16_PBS_", "C57BL/6", "H-2b", "B16-F10 subcutaneous melanoma"),
        ("CT26", "CT26_PBS_", "BALB/c", "H-2d", "CT26 subcutaneous colon carcinoma"),
        ("LLC", "LLC_PBS_", "C57BL/6", "H-2b", "Lewis lung carcinoma (LLC) subcutaneous graft"),
    ]:
        c_sub = [c for c in pbs_cols if c.startswith(prefix)]
        cohort_matrices.append(
            ("gse334005", code, c_sub, [f"GSE334005:{c}" for c in c_sub], peg_334[c_sub])
        )
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

    # 5. KP_LUAD (C57BL/6 primary Kras;Trp53 lung adenocarcinoma GEMM, GSE277076, n=10 sgControl primary tumors)
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
    kp_raw_cols = [
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
    tpm_kp = _counts_to_tpm(df_kp[kp_raw_cols].reindex(sym_universe, fill_value=0.0), exonic_series)
    kp_cols = _qc_cols(tpm_kp, kp_raw_cols)
    peg_kp = _project_peg_matrix(tpm_kp[kp_cols], groups_df)
    kp_sids = [f"KP_{c}" for c in kp_cols]
    peg_kp.columns = kp_sids
    cohort_matrices.append(
        ("gse277076", "KP_LUAD", kp_sids, [f"GSE277076:{c}" for c in kp_cols], peg_kp)
    )
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

    # 6. BMPC_PRAD (FVB primary & metastatic MYC;Pten prostate adenocarcinoma GEMM, GSE287166, n=6)
    df_bmpc = pd.read_csv(input_dir / "GSE287166_NormalizedCounts.csv.gz", index_col=0)
    df_bmpc["_gid"] = df_bmpc.index.astype(str).map(sym_to_ens)
    df_bmpc = df_bmpc.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    bmpc_raw_cols = ["BMPC_Pri1", "BMPC_Pri2", "BMPC_Pri3", "BMPC_Met1", "BMPC_Met2", "BMPC_Met3"]
    tpm_bmpc = _counts_to_tpm(
        df_bmpc[bmpc_raw_cols].reindex(sym_universe, fill_value=0.0), exonic_series
    )
    bmpc_cols = _qc_cols(tpm_bmpc, bmpc_raw_cols)
    peg_bmpc = _project_peg_matrix(tpm_bmpc[bmpc_cols], groups_df)
    cohort_matrices.append(
        ("gse287166", "BMPC_PRAD", bmpc_cols, [f"GSE287166:{c}" for c in bmpc_cols], peg_bmpc)
    )
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

    # Pure in-vitro cell lines moved out of in-vivo CA admission: 4T1 & P815
    p_4t1 = input_dir / "GSE332912_TPM_values.csv.gz"
    source_rows.append(
        {
            "source_id": "gse332912-4t1",
            "pillar": "tumor_rna",
            "filename": p_4t1.name,
            "unit": "TPM",
            "strain": "BALB/c (H-2d, 4T1 mammary carcinoma cell line)",
            "bytes": p_4t1.stat().st_size,
            "sha256": _sha256(p_4t1),
        }
    )
    df_4t1 = pd.read_csv(p_4t1)
    df_4t1["ID"] = df_4t1["ID"].astype(str).str.split(".").str[0]
    df_4t1 = df_4t1.drop_duplicates("ID").set_index("ID")
    tpm_4t1 = df_4t1[["XW1", "XW2", "XW3"]].astype(float)
    cols_4t1 = _qc_cols(tpm_4t1, ["XW1", "XW2", "XW3"])
    peg_4t1 = _project_peg_matrix(tpm_4t1[cols_4t1], groups_df)
    cell_line_matrices.append(("gse332912", "4T1", "MAMMARY", cols_4t1, peg_4t1))
    cell_line_registry.append(
        {
            "model_id": "4T1",
            "cancer_code": "MAMMARY",
            "label": "4T1 triple-negative mammary carcinoma cell line",
            "strain": "BALB/c",
            "mhc_haplotype": "H-2d",
            "source_accession": "GSE332912",
            "n_samples": len(cols_4t1),
            "native_feature_space": "Ensembl gene (ENSMUSG)",
            "quantification": "author linear TPM",
        }
    )

    p_p815 = input_dir / "GSE181773_gene_fpkm.txt.gz"
    source_rows.append(
        {
            "source_id": "gse181773-p815",
            "pillar": "tumor_rna",
            "filename": p_p815.name,
            "unit": "TPM (from native ENSMUSG FPKM)",
            "strain": "DBA/2 (H-2d, P815 mastocytoma cell line)",
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
    tpm_p815 = _fpkm_to_tpm(df_p815)
    p815_cols = _qc_cols(tpm_p815, ["FPKM.NC_1", "FPKM.NC_2", "FPKM.NC_3"])
    peg_p815 = _project_peg_matrix(tpm_p815[p815_cols], groups_df)
    cell_line_matrices.append(("gse181773", "P815", "MAST", p815_cols, peg_p815))
    cell_line_registry.append(
        {
            "model_id": "P815",
            "cancer_code": "MAST",
            "label": "P815 mastocytoma cell line (native DBA/2 source of P1A/Trap1a)",
            "strain": "DBA/2",
            "mhc_haplotype": "H-2d",
            "source_accession": "GSE181773",
            "n_samples": len(p815_cols),
            "native_feature_space": "Ensembl gene (ENSMUSG)",
            "quantification": "within-sample linear TPM from FPKM",
        }
    )

    if extra_dir.is_dir():
        # 7-8. EMT6 (BALB/c mammary graft) & RM1 (C57BL/6 prostate graft) in-vivo PBS tumors (GSE334006)
        p334_6 = extra_dir / "GSE334006_RNAseq_RM1_EMT6_tpm.tsv.gz"
        if p334_6.is_file():
            source_rows.append(
                {
                    "source_id": "gse334006-rm1-emt6-tpm",
                    "pillar": "tumor_rna",
                    "filename": p334_6.name,
                    "unit": "TPM",
                    "strain": "C57BL/6 (RM1) & BALB/c (EMT6) in-vivo grafts",
                    "bytes": p334_6.stat().st_size,
                    "sha256": _sha256(p334_6),
                }
            )
            df334_6 = pd.read_csv(p334_6, sep="\t")
            df334_6["_gid"] = df334_6["gene_id"].astype(str).map(entrez_to_ens)
            df334_6 = df334_6.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
            pbs334_6 = [c for c in df334_6.columns if "_PBS_" in c]
            tpm334_6 = df334_6[pbs334_6].astype(float)
            qc334_6 = _qc_cols(tpm334_6, pbs334_6)
            peg334_6 = _project_peg_matrix(tpm334_6[qc334_6], groups_df)
            for code, pfx, strain, hap, lbl in (
                ("EMT6", "EMT6_PBS_", "BALB/c", "H-2d", "EMT6 mammary carcinoma in-vivo graft"),
                ("RM1", "RM1_PBS_", "C57BL/6", "H-2b", "RM1 prostate carcinoma in-vivo graft"),
            ):
                cs_sub = [c for c in qc334_6 if c.startswith(pfx)]
                if cs_sub:
                    cohort_matrices.append(
                        (
                            "gse334006",
                            code,
                            cs_sub,
                            [f"GSE334006:{c}" for c in cs_sub],
                            peg334_6[cs_sub],
                        )
                    )
                    cohort_registry.append(
                        {
                            "cancer_code": code,
                            "label": lbl,
                            "model_kind": "syngeneic_graft",
                            "strain": strain,
                            "mhc_haplotype": hap,
                            "source_accession": "GSE334006",
                            "n_samples": len(cs_sub),
                            "native_feature_space": "NCBI Entrez GeneID -> MGI -> ENSMUSG",
                            "quantification": "author linear TPM",
                        }
                    )

        # 9. B78_MEL (GSE304809, C57BL/6 B78 melanoma in-vivo PBS grafts Sam3, Sam4, Sam5)
        p304 = extra_dir / "GSE304809_raw_counts_all_samples.tsv.gz"
        if p304.is_file():
            source_rows.append(
                {
                    "source_id": "gse304809-b78-melanoma",
                    "pillar": "tumor_rna",
                    "filename": p304.name,
                    "unit": "TPM (from RSEM raw counts)",
                    "strain": "C57BL/6 (H-2b, B78 melanoma in-vivo PBS grafts)",
                    "bytes": p304.stat().st_size,
                    "sha256": _sha256(p304),
                }
            )
            df304 = pd.read_csv(p304, sep="\t")
            df304["_gid"] = df304["Gene"].astype(str).map(sym_to_ens)
            df304 = df304.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
            b78_raw = ["Sam3", "Sam4", "Sam5"]
            tpm304 = _counts_to_tpm(
                df304[b78_raw].reindex(sym_universe, fill_value=0.0), exonic_series
            )
            b78_cols = _qc_cols(tpm304, b78_raw)
            if b78_cols:
                peg304 = _project_peg_matrix(tpm304[b78_cols], groups_df)
                b78_sids = [f"B78_{c}" for c in b78_cols]
                peg304.columns = b78_sids
                cohort_matrices.append(
                    (
                        "gse304809",
                        "B78_MEL",
                        b78_sids,
                        [f"GSE304809:{c}" for c in b78_cols],
                        peg304,
                    )
                )
                cohort_registry.append(
                    {
                        "cancer_code": "B78_MEL",
                        "label": "B78 amelanotic melanoma subcutaneous graft",
                        "model_kind": "syngeneic_graft",
                        "strain": "C57BL/6",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE304809",
                        "n_samples": len(b78_cols),
                        "native_feature_space": "MGI/Ensembl official symbol -> ENSMUSG",
                        "quantification": "linear TPM from raw counts and Ensembl 115 exonic_bp",
                    }
                )

        # 10-13. Orthotopic & autochthonous murine gliomas (GSE341165: GL261, CT-2A, SB28, MADR PBS tumors)
        p341 = extra_dir / "GSE341165_MurineGliomaBulkRNAseq_STAR_GeneCounts.GEO.txt.gz"
        if p341.is_file():
            source_rows.append(
                {
                    "source_id": "gse341165-murine-glioma-tumors",
                    "pillar": "tumor_rna",
                    "filename": p341.name,
                    "unit": "TPM (from STAR gene counts)",
                    "strain": "C57BL/6 orthotopic & autochthonous gliomas (GL261, CT2A, SB28, MADR)",
                    "bytes": p341.stat().st_size,
                    "sha256": _sha256(p341),
                }
            )
            df341 = pd.read_csv(p341, sep="\t")
            df341["_gid"] = df341["Sample"].astype(str).map(sym_to_ens)
            df341 = df341.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
            gli_pbs_cols = [c for c in df341.columns if c != "Sample" and "-PBS-" in c]
            tpm341 = _counts_to_tpm(
                df341[gli_pbs_cols].reindex(sym_universe, fill_value=0.0), exonic_series
            )
            qc341 = _qc_cols(tpm341, gli_pbs_cols)
            peg341 = _project_peg_matrix(tpm341[qc341], groups_df)
            for code, tok, mkind, lbl in (
                (
                    "GL261_TUMOR",
                    "-GL261-PBS-",
                    "syngeneic_graft",
                    "GL261 orthotopic intracranial glioma",
                ),
                (
                    "CT2A_TUMOR",
                    "-CT2A-PBS-",
                    "syngeneic_graft",
                    "CT-2A orthotopic high-grade astrocytoma",
                ),
                (
                    "SB28_TUMOR",
                    "-SB28-PBS-",
                    "syngeneic_graft",
                    "SB28 orthotopic poorly immunogenic glioblastoma",
                ),
                (
                    "MADR_GLIOMA",
                    "-MADR-GBO-PBS-",
                    "primary_gemm",
                    "MADR-mEGFRvIII autochthonous glioma",
                ),
            ):
                cs_g = [c for c in qc341 if tok in c]
                if cs_g:
                    cohort_matrices.append(
                        (
                            "gse341165",
                            code,
                            cs_g,
                            [f"GSE341165:{c}" for c in cs_g],
                            peg341[cs_g],
                        )
                    )
                    cohort_registry.append(
                        {
                            "cancer_code": code,
                            "label": lbl,
                            "model_kind": mkind,
                            "strain": "C57BL/6",
                            "mhc_haplotype": "H-2b",
                            "source_accession": "GSE341165",
                            "n_samples": len(cs_g),
                            "native_feature_space": "MGI/Ensembl official symbol -> ENSMUSG",
                            "quantification": "linear TPM from STAR raw counts and Ensembl 115 exonic_bp",
                        }
                    )

        # 14. DEN_HCC (GSE322844, carcinogen-induced primary murine hepatocellular carcinoma DEN1..6)
        p322_den = extra_dir / "GSE322844_DEN_FPKM.csv.gz"
        if p322_den.is_file():
            source_rows.append(
                {
                    "source_id": "gse322844-den-hcc",
                    "pillar": "tumor_rna",
                    "filename": p322_den.name,
                    "unit": "TPM (from native ENSMUSG FPKM)",
                    "strain": "C57BL/6 diethylnitrosamine-induced primary HCC",
                    "bytes": p322_den.stat().st_size,
                    "sha256": _sha256(p322_den),
                }
            )
            df322_d = pd.read_csv(p322_den, sep=";")
            df322_d["gene_id"] = df322_d["gene_id"].astype(str).str.split(".").str[0]
            den_raw = [c for c in df322_d.columns if c.startswith("DEN")]
            df322_d = df322_d.groupby("gene_id", sort=False)[den_raw].sum()
            tpm322_d = _fpkm_to_tpm(df322_d)
            den_cols = _qc_cols(tpm322_d, den_raw)
            if den_cols:
                tpm322_t = tpm322_d[den_cols].copy()
                tpm322_t.loc[tpm322_t.index.isin(risk_gids)] = np.nan
                peg322_d = _project_peg_matrix(tpm322_t, groups_df)
                cohort_matrices.append(
                    (
                        "gse322844",
                        "DEN_HCC",
                        den_cols,
                        [f"GSE322844:{c}" for c in den_cols],
                        peg322_d,
                    )
                )
                cohort_registry.append(
                    {
                        "cancer_code": "DEN_HCC",
                        "label": "Diethylnitrosamine (DEN) carcinogen-induced primary hepatocellular carcinoma",
                        "model_kind": "carcinogen_primary",
                        "strain": "C57BL/6",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE322844",
                        "n_samples": len(den_cols),
                        "native_feature_space": "Ensembl gene (ENSMUSG)",
                        "quantification": "within-sample linear TPM from FPKM",
                    }
                )

        # Additional pure in-vitro murine cancer cell lines: EO771, Hepa1-6, ID8, GL261_LINE, MB49, MOC1, MOC22
        p309 = extra_dir / "GSE309307_raw_counts-mouse.txt.gz"
        if p309.is_file():
            source_rows.append(
                {
                    "source_id": "gse309307-eo771-line",
                    "pillar": "tumor_rna",
                    "filename": p309.name,
                    "unit": "TPM (from native ENSMUSG raw counts)",
                    "strain": "C57BL/6 (H-2b, EO771 mammary carcinoma line WT untreated)",
                    "bytes": p309.stat().st_size,
                    "sha256": _sha256(p309),
                }
            )
            df309 = pd.read_csv(p309, sep="\t")
            df309.index = df309.index.astype(str).str.split(".").str[0]
            df309 = df309[~df309.index.duplicated(keep="first")]
            eo_raw = [c for c in df309.columns if str(c).startswith("WT_") and "_T" not in str(c)]
            tpm309 = _counts_to_tpm(df309[eo_raw], exonic_series)
            eo_qc = _qc_cols(tpm309, eo_raw)
            if eo_qc:
                peg309 = _project_peg_matrix(tpm309[eo_qc], groups_df)
                cell_line_matrices.append(("gse309307", "EO771", "MAMMARY", eo_qc, peg309))
                cell_line_registry.append(
                    {
                        "model_id": "EO771",
                        "cancer_code": "MAMMARY",
                        "label": "EO771 mammary adenocarcinoma cell line (WT untreated)",
                        "strain": "C57BL/6",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE309307",
                        "n_samples": len(eo_qc),
                        "native_feature_space": "Ensembl gene (ENSMUSG)",
                        "quantification": "linear TPM from raw counts and Ensembl 115 exonic_bp",
                    }
                )

        p329 = extra_dir / "GSE329164_gene.count.matrix.annotation.csv.gz"
        if p329.is_file():
            source_rows.append(
                {
                    "source_id": "gse329164-hepa16-line",
                    "pillar": "tumor_rna",
                    "filename": p329.name,
                    "unit": "TPM (from native ENSMUSG raw counts)",
                    "strain": "C57BL/6L (H-2b, Hepa1-6 hepatoma cell line DMSO vehicle)",
                    "bytes": p329.stat().st_size,
                    "sha256": _sha256(p329),
                }
            )
            df329 = pd.read_csv(p329, usecols=["Gene id", "Control1", "Control2", "Control3"])
            df329["Gene id"] = df329["Gene id"].astype(str).str.split(".").str[0]
            df329 = df329.drop_duplicates("Gene id").set_index("Gene id")
            hepa_raw = ["Control1", "Control2", "Control3"]
            tpm329 = _counts_to_tpm(df329[hepa_raw], exonic_series)
            hepa_qc = _qc_cols(tpm329, hepa_raw)
            if hepa_qc:
                peg329 = _project_peg_matrix(tpm329[hepa_qc], groups_df)
                cell_line_matrices.append(("gse329164", "Hepa1-6", "HCC", hepa_qc, peg329))
                cell_line_registry.append(
                    {
                        "model_id": "Hepa1-6",
                        "cancer_code": "HCC",
                        "label": "Hepa1-6 hepatocellular carcinoma cell line (DMSO controls)",
                        "strain": "C57BL/6L",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE329164",
                        "n_samples": len(hepa_qc),
                        "native_feature_space": "Ensembl gene (ENSMUSG)",
                        "quantification": "linear TPM from raw counts and Ensembl 115 exonic_bp",
                    }
                )

        p349_3 = extra_dir / "GSE349387_OVCAR5_OVCAR8_ID8_Vehicle_NDNB25_expression_matrix.tsv.gz"
        if p349_3.is_file():
            source_rows.append(
                {
                    "source_id": "gse349387-id8-line",
                    "pillar": "tumor_rna",
                    "filename": p349_3.name,
                    "unit": "TPM (from native ENSMUSG counts)",
                    "strain": "C57BL/6 (H-2b, ID8 ovarian carcinoma cell line vehicle)",
                    "bytes": p349_3.stat().st_size,
                    "sha256": _sha256(p349_3),
                }
            )
            id8_raw = [
                "ID8_Vehicle_Rep1_count",
                "ID8_Vehicle_Rep2_count",
                "ID8_Vehicle_Rep3_count",
            ]
            df349_3 = pd.read_csv(p349_3, sep="\t", usecols=["gene_id", *id8_raw])
            df349_3 = df349_3[df349_3["gene_id"].astype(str).str.startswith("ENSMUSG")].copy()
            df349_3["gene_id"] = df349_3["gene_id"].astype(str).str.split(".").str[0]
            df349_3 = df349_3.drop_duplicates("gene_id").set_index("gene_id")
            tpm349_3 = _counts_to_tpm(df349_3[id8_raw], exonic_series)
            id8_qc = _qc_cols(tpm349_3, id8_raw)
            if id8_qc:
                peg349_3 = _project_peg_matrix(tpm349_3[id8_qc], groups_df)
                cell_line_matrices.append(("gse349387", "ID8", "OV", id8_qc, peg349_3))
                cell_line_registry.append(
                    {
                        "model_id": "ID8",
                        "cancer_code": "OV",
                        "label": "ID8 ovarian surface epithelial carcinoma cell line (vehicle)",
                        "strain": "C57BL/6",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE349387",
                        "n_samples": len(id8_qc),
                        "native_feature_space": "Ensembl gene (ENSMUSG)",
                        "quantification": "linear TPM from expected counts and Ensembl 115 exonic_bp",
                    }
                )

        p349_6 = extra_dir / "GSE349688_Processed-Data-File.txt.gz"
        if p349_6.is_file():
            source_rows.append(
                {
                    "source_id": "gse349688-gl261-line",
                    "pillar": "tumor_rna",
                    "filename": p349_6.name,
                    "unit": "TPM (native ENSMUSG)",
                    "strain": "C57BL/6 (H-2b, GL261 glioma cell line untreated controls)",
                    "bytes": p349_6.stat().st_size,
                    "sha256": _sha256(p349_6),
                }
            )
            gl_cl_raw = ["TPM (Control-1)", "TPM (Control-2)", "TPM (Control-3)"]
            df349_6 = pd.read_csv(
                p349_6,
                sep="\t",
                encoding="utf-16",
                comment="#",
                usecols=["Gene ID", *gl_cl_raw],
            )
            df349_6["Gene ID"] = df349_6["Gene ID"].astype(str).str.split(".").str[0]
            df349_6 = df349_6.drop_duplicates("Gene ID").set_index("Gene ID")
            tpm349_6 = df349_6[gl_cl_raw].astype(float)
            gl_cl_qc = _qc_cols(tpm349_6, gl_cl_raw)
            if gl_cl_qc:
                peg349_6 = _project_peg_matrix(tpm349_6[gl_cl_qc], groups_df)
                cell_line_matrices.append(("gse349688", "GL261_LINE", "GLIOMA", gl_cl_qc, peg349_6))
                cell_line_registry.append(
                    {
                        "model_id": "GL261_LINE",
                        "cancer_code": "GLIOMA",
                        "label": "GL261 glioma cell line (in-vitro untreated controls)",
                        "strain": "C57BL/6",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE349688",
                        "n_samples": len(gl_cl_qc),
                        "native_feature_space": "Ensembl gene (ENSMUSG)",
                        "quantification": "author linear TPM",
                    }
                )

        p350 = extra_dir / "GSE350232_tpm.txt.gz"
        if p350.is_file():
            source_rows.append(
                {
                    "source_id": "gse350232-mb49-line",
                    "pillar": "tumor_rna",
                    "filename": p350.name,
                    "unit": "TPM (native ENSMUSG)",
                    "strain": "C57BL/6 (H-2b, MB49 bladder carcinoma cell line controls)",
                    "bytes": p350.stat().st_size,
                    "sha256": _sha256(p350),
                }
            )
            mb_raw = ["C_2", "C_3", "C_4"]
            df350 = pd.read_csv(p350, sep="\t", usecols=["GeneID", *mb_raw])
            df350["GeneID"] = df350["GeneID"].astype(str).str.split(".").str[0]
            df350 = df350.drop_duplicates("GeneID").set_index("GeneID")
            tpm350 = df350[mb_raw].astype(float)
            mb_qc = _qc_cols(tpm350, mb_raw)
            if mb_qc:
                peg350 = _project_peg_matrix(tpm350[mb_qc], groups_df)
                cell_line_matrices.append(("gse350232", "MB49", "BLCA", mb_qc, peg350))
                cell_line_registry.append(
                    {
                        "model_id": "MB49",
                        "cancer_code": "BLCA",
                        "label": "MB49 urothelial bladder carcinoma cell line (control medium)",
                        "strain": "C57BL/6",
                        "mhc_haplotype": "H-2b",
                        "source_accession": "GSE350232",
                        "n_samples": len(mb_qc),
                        "native_feature_space": "Ensembl gene (ENSMUSG)",
                        "quantification": "author linear TPM",
                    }
                )

        for moc_fn, moc_id, moc_cols_raw, moc_lbl in (
            (
                "GSE326892_MOC1expression.txt.gz",
                "MOC1",
                ["A1", "A2", "A3"],
                "MOC1 indolent immunogenic oral squamous carcinoma cell line",
            ),
            (
                "GSE326892_MOC22expression.txt.gz",
                "MOC22",
                ["A7", "A8", "A9"],
                "MOC22 aggressive oral squamous carcinoma cell line",
            ),
        ):
            p_moc = extra_dir / moc_fn
            if p_moc.is_file():
                source_rows.append(
                    {
                        "source_id": f"gse326892-{moc_id.lower()}-line",
                        "pillar": "tumor_rna",
                        "filename": p_moc.name,
                        "unit": "TPM (from log2(count+1) raw counts)",
                        "strain": f"C57BL/6 (H-2b, {moc_id} oral cancer cell line controls)",
                        "bytes": p_moc.stat().st_size,
                        "sha256": _sha256(p_moc),
                    }
                )
                df_moc = pd.read_csv(p_moc, sep="\t")
                df_moc["_gid"] = df_moc["id"].astype(str).map(sym_to_ens)
                df_moc = df_moc.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
                cnt_moc = np.round(np.power(2.0, df_moc[moc_cols_raw].astype(float)) - 1.0).clip(
                    lower=0.0
                )
                tpm_moc = _counts_to_tpm(
                    cnt_moc.reindex(sym_universe, fill_value=0.0), exonic_series
                )
                moc_qc = _qc_cols(tpm_moc, moc_cols_raw)
                if moc_qc:
                    peg_moc = _project_peg_matrix(tpm_moc[moc_qc], groups_df)
                    cell_line_matrices.append(("gse326892", moc_id, "HNSCC", moc_qc, peg_moc))
                    cell_line_registry.append(
                        {
                            "model_id": moc_id,
                            "cancer_code": "HNSCC",
                            "label": moc_lbl,
                            "strain": "C57BL/6",
                            "mhc_haplotype": "H-2b",
                            "source_accession": "GSE326892",
                            "n_samples": len(moc_qc),
                            "native_feature_space": "MGI/Ensembl official symbol -> ENSMUSG",
                            "quantification": "linear TPM from HTSeq raw counts and Ensembl 115 exonic_bp",
                        }
                    )

    eids_list = groups_df["entity_id"].tolist()
    non_unsafe_mask = np.array([eid not in unsafe_eids for eid in eids_list], dtype=bool)
    best_below10: dict[str, tuple[float, str, str, str, str]] = {}
    tumor_rows: list[dict] = []

    for src_id, code, s_cols, p_ids, peg_mat in cohort_matrices:
        if not s_cols:
            continue
        arr = peg_mat[s_cols].to_numpy(dtype=np.float64)
        r_idxs, c_idxs = np.where(non_unsafe_mask[:, None] & np.isfinite(arr) & (arr >= 10.0))
        for ri, ci in zip(r_idxs.tolist(), c_idxs.tolist()):
            tumor_rows.append(
                {
                    "entity_id": eids_list[ri],
                    "source_id": src_id,
                    "sample_id": s_cols[ci],
                    "participant_id": p_ids[ci],
                    "cancer_code": code,
                    "value": float(arr[ri, ci]),
                    "comparison_id": COMPARISON_ID,
                    "unit": "TPM",
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )
        safe_arr = np.where(np.isfinite(arr), arr, -np.inf)
        row_max = np.max(safe_arr, axis=1)
        row_arg = np.argmax(safe_arr, axis=1)
        for ri, mx in enumerate(row_max.tolist()):
            if not non_unsafe_mask[ri] or not np.isfinite(mx) or mx >= 10.0:
                continue
            eid = eids_list[ri]
            prev = best_below10.get(eid)
            if prev is None or mx > prev[0]:
                ci = int(row_arg[ri])
                best_below10[eid] = (mx, src_id, s_cols[ci], p_ids[ci], code)

    has_ge10 = {r["entity_id"] for r in tumor_rows}
    for eid, (mx, src_id, smp_id, part_id, code) in best_below10.items():
        if eid not in has_ge10:
            tumor_rows.append(
                {
                    "entity_id": eid,
                    "source_id": src_id,
                    "sample_id": smp_id,
                    "participant_id": part_id,
                    "cancer_code": code,
                    "value": float(mx),
                    "comparison_id": COMPARISON_ID,
                    "unit": "TPM",
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )

    cl_best_below10: dict[str, tuple[float, str, str, str, str]] = {}
    cl_rows: list[dict] = []
    for src_id, m_id, code, s_cols, peg_mat in cell_line_matrices:
        if not s_cols:
            continue
        arr_cl = peg_mat[s_cols].to_numpy(dtype=np.float64)
        r_idxs, c_idxs = np.where(non_unsafe_mask[:, None] & np.isfinite(arr_cl) & (arr_cl >= 10.0))
        for ri, ci in zip(r_idxs.tolist(), c_idxs.tolist()):
            cl_rows.append(
                {
                    "entity_id": eids_list[ri],
                    "source_id": src_id,
                    "sample_id": f"{m_id}:{s_cols[ci]}",
                    "model_id": m_id,
                    "cancer_code": code,
                    "value": float(arr_cl[ri, ci]),
                    "comparison_id": COMPARISON_ID,
                    "unit": "TPM",
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )
        safe_cl = np.where(np.isfinite(arr_cl), arr_cl, -np.inf)
        cl_max = np.max(safe_cl, axis=1)
        cl_arg = np.argmax(safe_cl, axis=1)
        for ri, mx in enumerate(cl_max.tolist()):
            if not non_unsafe_mask[ri] or not np.isfinite(mx) or mx >= 10.0:
                continue
            eid = eids_list[ri]
            prev = cl_best_below10.get(eid)
            if prev is None or mx > prev[0]:
                ci = int(cl_arg[ri])
                cl_best_below10[eid] = (mx, src_id, f"{m_id}:{s_cols[ci]}", m_id, code)

    cl_has_ge10 = {r["entity_id"] for r in cl_rows}
    for eid, (mx, src_id, smp_id, m_id, code) in cl_best_below10.items():
        if eid not in cl_has_ge10:
            cl_rows.append(
                {
                    "entity_id": eid,
                    "source_id": src_id,
                    "sample_id": smp_id,
                    "model_id": m_id,
                    "cancer_code": code,
                    "value": float(mx),
                    "comparison_id": COMPARISON_ID,
                    "unit": "TPM",
                    "qc_pass": True,
                    "identity_resolved": True,
                }
            )

    return (
        pd.DataFrame(tumor_rows),
        pd.DataFrame(cl_rows),
        pd.DataFrame(cohort_registry),
        pd.DataFrame(cell_line_registry),
        source_rows,
    )


def build_all(
    ens_dir: Path, input_dir: Path, out_root: Path, extra_dir: Path | None = None
) -> dict:
    eff_extra = extra_dir if extra_dir is not None else DEFAULT_EXTRA_DIR
    print("1/6 Qualifying Ensembl 115 GRCm39 proteome & 31-mer cDNA mapping risk...", flush=True)
    ref, groups_df, seq_df, exonic_bp, symbol_by_gene, mapping_risk_df, ens_sources = (
        load_ensembl115_mouse_reference(ens_dir)
    )
    coding_genes = set(symbol_by_gene)
    risk_gids = set(
        mapping_risk_df.loc[mapping_risk_df["mapping_risk"].astype(bool), "gene_id"].astype(str)
    )
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
        "3/6 Aggregating Bgee 15.2, GSE219045, WT liver & 16-compartment reproductive/embryo TPM panels...",
        flush=True,
    )
    comp_normals, comp_repro, rna_gates, confirmed_organs, norm_sources = (
        load_bgee_and_normal_panels(
            input_dir, eff_extra, groups_df, exonic_bp, symbol_by_gene, entrez_to_ens, sym_to_ens
        )
    )
    print(
        f"    normal_obs={len(comp_normals):,}, repro_obs={len(comp_repro):,}, rna_gates={len(rna_gates):,}",
        flush=True,
    )

    print(
        "4/6 Evaluating ImmGen/splenic/IL-15/IL-2 immune RNA, MGI GXD TS28 protein & Tabula Muris 10x UMI...",
        flush=True,
    )
    imm_gates, imm_sources = load_mouse_immune_panels(
        input_dir, eff_extra, groups_df, exonic_bp, symbol_by_gene, sym_to_ens
    )
    gxd_gates, gxd_sources = load_mgi_gxd_protein_gates(
        input_dir, groups_df, mgi_to_ens, confirmed_organs
    )
    umi_df, tm_sources = load_tabula_muris_rare_cell_umi(input_dir, sym_to_ens)
    print(
        f"    imm_gates={len(imm_gates):,}, gxd_gates={len(gxd_gates):,}, umi_rows={len(umi_df):,}",
        flush=True,
    )

    all_gates = pd.concat([rna_gates, imm_gates, gxd_gates], ignore_index=True)
    unsafe_eids = set(all_gates.loc[all_gates["unsafe"].astype(bool), "entity_id"].astype(str))

    print("5/6 Loading 14 in-vivo murine tumor cohorts & 9 pure cancer cell lines...", flush=True)
    comp_tumors, comp_cell_lines, cohort_df, cell_line_df, tumor_sources = load_mouse_tumor_cohorts(
        input_dir,
        eff_extra,
        groups_df,
        exonic_bp,
        symbol_by_gene,
        entrez_to_ens,
        sym_to_ens,
        risk_gids,
        unsafe_eids,
    )
    print(
        f"    tumor_obs={len(comp_tumors):,}, cohorts={len(cohort_df)} ({int(cohort_df['n_samples'].sum())} in-vivo samples), cell_lines={len(cell_line_df)}",
        flush=True,
    )

    print(
        "6/6 Running SpeciesScreeningSpec & 16-compartment HostContext sex-restriction...",
        flush=True,
    )
    spec = SpeciesScreeningSpec(
        species="mouse",
        reference=ref,
        policy_id=POLICY_SCHEMA_ID,
        comparison_id=COMPARISON_ID,
        comparability_reference="oncoref.mouse-grcm39-panel.2026-11-26.2",
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
        comparison_cell_lines=comp_cell_lines if not comp_cell_lines.empty else None,
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
    cell_line_df.to_csv(mouse_dir / "cell-lines.csv", index=False, lineterminator="\n")

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
        "cell_line_models": cell_line_df["model_id"].tolist(),
        "cell_line_samples": int(cell_line_df["n_samples"].sum()),
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
            "cta_activation_modes": {
                k: int(v)
                for k, v in cat.loc[cat["is_cta"].astype(bool), "tumor_activation_mode"]
                .value_counts()
                .sort_index()
                .items()
            },
        },
    }
    (mouse_dir / "target-selection.json").write_text(json.dumps(summary_payload, indent=2) + "\n")
    print("DONE! Summary:", json.dumps(summary_payload["summary"], indent=2), flush=True)
    return summary_payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ens-dir", type=Path, default=DEFAULT_ENS_DIR)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--extra-dir", type=Path, default=DEFAULT_EXTRA_DIR)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    build_all(args.ens_dir, args.input_dir, args.root, args.extra_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
