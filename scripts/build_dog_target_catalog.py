#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0.
"""Build the quantitative Ensembl 98 CanFam3.1 canine antigen catalog.

Evaluates all protein-coding loci in Canis lupus familiaris (taxon 9615,
Ensembl 98 CanFam3.1 / GCA_000002285.2, cross-referenced to NCBI RefSeq 106
UU_Cfam_GSD_1.0) through the unified 3-stage OncoRef target-selection architecture:

1. Reference proteome & PEG grouping: 43,727 Ensembl 98 CanFam3.1 peptide
   translations across 20,199 protein-coding genes collapsed into identical-sequence
   PEG:<sha256[:20]> groups, plus CanFam3.1 exon geometry, 31-mer cDNA paralog
   mapping-risk flags, and complete symbol/GeneID reconciliation for all 14 curated
   canine CTA literature nominations (hindriks2025, nemec2019).
2. Multi-modal normal safety gates:
   - Adult somatic bulk & leukocyte RNA (25 healthy Bgee 15.2 experiments with
     SRA run deduplication + GSE219045 300-sample sex-stratified multi-organ atlas
     in linear TPM; Vital <= 4, Protected <= 8, Manageable <= 16 hold / > 128
     Unsafe, canine Disposable <= 32 TPM across gallbladder, tonsil, and hair
     follicle; uniform <= 2.5% ectopic OR gDNA QC; zero negative-IHC multipliers).
   - 38-tissue quantitative mass-spectrometry proteomics (Pfizer MSV000099744,
     9,929 protein groups mapped to CanFam3.1 PEGs, converted to within-channel
     molar ppm with 10% same-plex reproductive TMT co-isolation bleed subtraction;
     Vital > 25 ppm Unsafe / 10-25 ppm hold; Protected > 50 ppm Unsafe / 20-50 ppm
     hold).
   - Healthy companion-dog circulating leukocyte 10x scRNA-seq (GSE225599, 7
     healthy donors H_1..H_7 across T/NK, monocyte, neutrophil, and platelet
     lineages), alveolar macrophages (GSE283373), activated-TIL checkpoint/cytokine
     gates, and somatic RT-PCR counterevidence holds (nemec2019 / hindriks2025).
3. Symmetric 10 x max(N, D, 1 TPM) Cancer Antigen (CA), Reproductive Antigen
   (RA), and CTA = CA & RA enrichment across 19 spontaneous canine patient tumor
   cohorts (935 primary tumors; 12 pure canine cancer cell lines segregated in a
   non-promoting purity/model-selection sidecar) and comparative canine reproductive,
   bicornuate-tract, mammary, placental, and d33-d44 embryonic compartments.
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
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oncoref import canine, target_selection  # noqa: E402
from oncoref.reference_identity import ReferenceIdentity  # noqa: E402
from oncoref.species_policies import POLICY_SCHEMA_ID  # noqa: E402
from oncoref.species_target_selection import (  # noqa: E402
    SpeciesScreeningSpec,
    anatomy_for_species,
    classify_sex_restriction,
    evaluate_species_catalog,
    group_reference_proteins,
)

DEFAULT_SCRATCH_DIR = Path(
    os.environ.get("ONCOREF_DOG_INPUTS_DIR", str(ROOT / "tmp" / "dog_catalog_inputs"))
)
DEFAULT_EXTRA_DIR = Path(
    os.environ.get("ONCOREF_DOG_EXTRA_DIR", str(DEFAULT_SCRATCH_DIR.parent / "dog_153_extra"))
)
COMPARISON_ID = "oncoref.dog-canfam31-tpm.2026-11-26.2"

# 22 adult somatic contexts measured across both Bgee 15.2 & GSE219045 on CanFam3.1
REQUIRED_DOG_NORMAL_CONTEXTS = frozenset(
    {
        "adipose tissue",
        "adrenal gland",
        "blood",
        "bone",
        "bone marrow",
        "brain",
        "colon",
        "eye",
        "heart",
        "kidney",
        "liver",
        "lung",
        "lymph node",
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
    "eye": "vital",
    "heart": "vital",
    "kidney": "vital",
    "liver": "vital",
    "lung": "vital",
    "pancreas": "vital",
    "peripheral nerve": "vital",
    "pituitary gland": "vital",
    "skeletal muscle": "vital",
    "small intestine": "vital",
    "smooth muscle": "vital",
    "spinal cord": "vital",
    # Protected (ceiling <= 8 TPM)
    "buccal mucosa": "protected",
    "cartilage": "protected",
    "colon": "protected",
    "joint synovium": "protected",
    "lacrimal gland": "protected",
    "lymph node": "protected",
    "parathyroid gland": "protected",
    "salivary gland": "protected",
    "skin": "protected",
    "spleen": "protected",
    "stomach": "protected",
    "tendon": "protected",
    "tongue": "protected",
    "ureter": "protected",
    "urinary bladder": "protected",
    # Manageable (hold > 16 TPM, Unsafe > 128 TPM)
    "adipose tissue": "manageable",
    "olfactory epithelium": "manageable",
    "thyroid gland": "manageable",
}

# Disease-model Bgee canine experiments excluded from healthy adult normal N
EXCLUDED_DISEASE_BGEE_EXPS = frozenset(
    {
        "SRP051588",  # degenerative mitral valve disease
        "DRP005124",  # tachypacing heart failure
        "ERP114741",  # hereditary nasal parakeratosis
        "SRP081080",  # juvenile dermatomyositis
        "SRP110851",  # atopic dermatitis skin
        "ERP015721",  # X-linked ectodermal dysplasia
    }
)

BGEE_DOG_ANATOMY_MAP: dict[str, tuple[str, str]] = {
    "brain": ("normal", "brain"),
    "cerebellum": ("normal", "brain"),
    "cerebral cortex": ("normal", "brain"),
    "frontal cortex": ("normal", "brain"),
    "occipital cortex": ("normal", "brain"),
    "occipital lobe": ("normal", "brain"),
    "prefrontal cortex": ("normal", "brain"),
    "temporal lobe": ("normal", "brain"),
    "hypothalamus": ("normal", "brain"),
    "olfactory bulb": ("normal", "brain"),
    "spinal cord": ("normal", "spinal cord"),
    "retina": ("normal", "eye"),
    "heart": ("normal", "heart"),
    "heart left ventricle": ("normal", "heart"),
    "heart right ventricle": ("normal", "heart"),
    "cardiac muscle of left ventricle": ("normal", "heart"),
    "left cardiac atrium": ("normal", "heart"),
    "right cardiac atrium": ("normal", "heart"),
    "lung": ("normal", "lung"),
    "liver": ("normal", "liver"),
    "adult mammalian kidney": ("normal", "kidney"),
    "kidney": ("normal", "kidney"),
    "metanephros cortex": ("normal", "kidney"),
    "renal medulla": ("normal", "kidney"),
    "pancreas": ("normal", "pancreas"),
    "stomach": ("normal", "stomach"),
    "small intestine": ("normal", "small intestine"),
    "jejunum": ("normal", "small intestine"),
    "colon": ("normal", "colon"),
    "large intestine": ("normal", "colon"),
    "skeletal muscle tissue": ("normal", "skeletal muscle"),
    "muscle tissue": ("normal", "skeletal muscle"),
    "hindlimb muscle": ("normal", "skeletal muscle"),
    "smooth muscle tissue": ("normal", "smooth muscle"),
    "bone marrow": ("normal", "bone marrow"),
    "cartilage tissue": ("normal", "cartilage"),
    "blood": ("normal", "blood"),
    "granulocyte": ("normal", "blood"),
    "spleen": ("normal", "spleen"),
    "lymph node": ("normal", "lymph node"),
    "urinary bladder": ("normal", "urinary bladder"),
    "mucosa of urinary bladder": ("normal", "urinary bladder"),
    "saliva-secreting gland": ("normal", "salivary gland"),
    "tongue": ("normal", "tongue"),
    "adrenal gland": ("normal", "adrenal gland"),
    "adrenal cortex": ("normal", "adrenal gland"),
    "pituitary gland": ("normal", "pituitary gland"),
    "thyroid gland": ("normal", "thyroid gland"),
    "adipose tissue": ("normal", "adipose tissue"),
    "olfactory segment of nasal mucosa": ("normal", "olfactory epithelium"),
    "nose": ("normal", "olfactory epithelium"),
    "zone of skin": ("normal", "skin"),
    "keratinocyte": ("normal", "skin"),
    # Disposable ("disposable")
    "hair follicle": ("disposable", "hair follicle"),
    # Thymic ("thymic")
    "thymus": ("thymic", "thymus"),
    # Reproductive ("reproductive")
    "testis": ("reproductive", "testis"),
    "ovary": ("reproductive", "ovary"),
    "placenta": ("reproductive", "placenta"),
    "mammary gland": ("reproductive", "resting_mammary"),
}

GSE219045_DOG_ORGAN_MAP: dict[str, tuple[str, str]] = {
    "AdrenalGlands": ("normal", "adrenal gland"),
    "AdrenalGland": ("normal", "adrenal gland"),
    "Aorta": ("normal", "blood vessel"),
    "FemoralVein": ("normal", "blood vessel"),
    "BoneMarrow": ("normal", "bone marrow"),
    "PhysisMetaphysis": ("normal", "bone"),
    "BuccalMucosa": ("normal", "buccal mucosa"),
    "CaudateNucleus": ("normal", "brain"),
    "Cerebellum": ("normal", "brain"),
    "FrontalLobe": ("normal", "brain"),
    "Hippocampus": ("normal", "brain"),
    "SubcorticalWhiteMatter": ("normal", "brain"),
    "Thalamus": ("normal", "brain"),
    "CervicalSpinalCord": ("normal", "spinal cord"),
    "DorsalRootGanglia": ("normal", "peripheral nerve"),
    "SciaticNerve": ("normal", "peripheral nerve"),
    "Cornea": ("normal", "eye"),
    "Lens": ("normal", "eye"),
    "Retina": ("normal", "eye"),
    "Uvea": ("normal", "eye"),
    "LacrimalGland": ("normal", "lacrimal gland"),
    "Heart": ("normal", "heart"),
    "Lung": ("normal", "lung"),
    "Liver": ("normal", "liver"),
    "Kidney": ("normal", "kidney"),
    "Pancreas": ("normal", "pancreas"),
    "FundicStomach": ("normal", "stomach"),
    "PyloricStomach": ("normal", "stomach"),
    "Duodenum": ("normal", "small intestine"),
    "Jejunum": ("normal", "small intestine"),
    "Ileum": ("normal", "small intestine"),
    "Colon": ("normal", "colon"),
    "MesentericAdiposeTissue": ("normal", "adipose tissue"),
    "MesentericAdipose": ("normal", "adipose tissue"),
    "MesentericLymphNode": ("normal", "lymph node"),
    "PeyersPatch": ("normal", "lymph node"),
    "Spleen": ("normal", "spleen"),
    "WholeBlood": ("normal", "blood"),
    "OlfactoryEpithelium": ("normal", "olfactory epithelium"),
    "Parathyroid": ("normal", "parathyroid gland"),
    "Thyroid": ("normal", "thyroid gland"),
    "ThyroidParathyroid": ("normal", "parathyroid gland"),
    "PituitaryGland": ("normal", "pituitary gland"),
    "Pituitary": ("normal", "pituitary gland"),
    "ParotidGland": ("normal", "salivary gland"),
    "ParotidSalivaryGland": ("normal", "salivary gland"),
    "SubmandibularSalivaryGland": ("normal", "salivary gland"),
    "Quadriceps": ("normal", "skeletal muscle"),
    "Soleus": ("normal", "skeletal muscle"),
    "GastrocnemiusTendon": ("normal", "tendon"),
    "Synovium": ("normal", "joint synovium"),
    "Tongue": ("normal", "tongue"),
    "UrinaryBladder": ("normal", "urinary bladder"),
    # Disposable
    "Gallbladder": ("disposable", "gallbladder"),
    # Male & female reproductive
    "Testes": ("reproductive", "testis"),
    "EpididymisHead": ("reproductive", "epididymis"),
    "EpididymisTail": ("reproductive", "epididymis"),
    "Prostate": ("reproductive", "prostate"),
    "CorpusSpongiosum": ("reproductive", "prostate"),
    "Ovaries": ("reproductive", "ovary"),
    "Uterus": ("reproductive", "uterus"),
    "Cervix": ("reproductive", "cervix"),
    "Vagina": ("reproductive", "vagina"),
    "MammaryTissue": ("reproductive", "resting_mammary"),
}

PFIZER_CONDITION_MAP: dict[str, tuple[str, str]] = {
    # Somatic Vital
    "Adrenal": ("normal", "adrenal gland"),
    "Bone_Marrow": ("normal", "bone marrow"),
    "Cerebral_Cortex": ("normal", "brain"),
    "Cerebellum": ("normal", "brain"),
    "White_Blood_Cells": ("normal", "blood"),
    "Pancreas": ("normal", "pancreas"),
    "GI_Duodenum": ("normal", "small intestine"),
    "GI_Jejunum": ("normal", "small intestine"),
    "GI_Ileum": ("normal", "small intestine"),
    "Heart": ("normal", "heart"),
    "Kidney": ("normal", "kidney"),
    "Skeletal_Muscle": ("normal", "skeletal muscle"),
    "Eye": ("normal", "eye"),
    "Endothelium": ("normal", "blood vessel"),
    "Lung": ("normal", "lung"),
    "Pituitary": ("normal", "pituitary gland"),
    "Spinal_Cord": ("normal", "spinal cord"),
    "Liver": ("normal", "liver"),
    "Peripheral_Nerve": ("normal", "peripheral nerve"),
    # Somatic Protected
    "GI_Stomach": ("normal", "stomach"),
    "Parotid_Gland": ("normal", "salivary gland"),
    "Bladder": ("normal", "urinary bladder"),
    "Ureter": ("normal", "ureter"),
    "Lymph_Node": ("normal", "lymph node"),
    "Parathyroid": ("normal", "parathyroid gland"),
    "Spleen": ("normal", "spleen"),
    "GI_Colon": ("normal", "colon"),
    "Skin": ("normal", "skin"),
    # Somatic Manageable
    "Thyroid": ("normal", "thyroid gland"),
    # Thymic & Disposable
    "Thymus": ("thymic", "thymus"),
    "Tonsil": ("disposable", "tonsil"),
    # Reproductive
    "Testis": ("reproductive", "testis"),
    "Prostate": ("reproductive", "prostate"),
    "Ovary": ("reproductive", "ovary"),
    "Fallopian_Tube": ("reproductive", "fallopian tube"),
    "Uterus_Endometrium": ("reproductive", "endometrium"),
    "Uterus_Cervix": ("reproductive", "cervix"),
    "Breast/Mammary": ("reproductive", "resting_mammary"),
}

DENTAL_HOME_SYMBOLS = frozenset(
    {"AMELX", "AMBN", "ENAM", "DSPP", "MMP20", "KLK4", "ODAPH", "ODAM", "AMTN"}
)
TASTE_HOME_SYMBOLS = frozenset({"OTOP1"})
HAIR_FOLLICLE_KERATINS = frozenset(
    {*(f"KRT{i}" for i in range(25, 41)), *(f"KRT{i}" for i in range(71, 87))}
)

# Verified coordinate/sequence crosswalk for the 14 curated canine nominations onto CanFam3.1
CURATED_NOMINATION_TO_ENSCAFG: dict[str, tuple[str, ...]] = {
    "AKR1E2": ("ENSCAFG00000005329",),
    "CAGE1": ("ENSCAFG00000009556",),
    "CTCFL": ("ENSCAFG00000012058",),
    "DDX53": ("ENSCAFG00000045979",),
    "IGF2BP3": ("ENSCAFG00000002766",),
    "MAGE": ("ENSCAFG00000013604", "ENSCAFG00000018988"),
    "MAGEA13P": ("ENSCAFG00000018988",),
    "MAGEB10": ("ENSCAFG00000013604",),
    "NT5C1B": ("ENSCAFG00000003760",),
    "PAGE3-like": ("ENSCAFG00000031355",),
    "PRAME": ("ENSCAFG00000013613",),
    "SPECC1": ("ENSCAFG00000018081",),
    "SSX6": ("ENSCAFG00000015307",),
    "TPX2": ("ENSCAFG00000007075",),
}
CURATED_SOMATIC_COUNTEREVIDENCE_NOMS = frozenset(
    {"AKR1E2", "CTCFL", "DDX53", "IGF2BP3", "NT5C1B", "SPECC1", "TPX2"}
)
# Essential canine activated-TIL checkpoints & effector cytokines (100% Unsafe in Human)
CURATED_CANINE_ACTIVATED_IMMUNE_SYMBOLS = frozenset(
    {
        "PDCD1",
        "TNFRSF4",
        "TNFSF4",
        "IFNG",
        "IFNB1",
        "IL17A",
        "IL17F",
        "CD70",
        "TNFRSF8",
        "SIGLEC15",
    }
)


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


def _counts_to_tpm(counts: pd.DataFrame, lengths_bp: pd.Series) -> pd.DataFrame:
    valid_len = lengths_bp.reindex(counts.index)
    keep = valid_len.notna() & valid_len.gt(0)
    sub = counts.loc[keep].astype(float).clip(lower=0.0)
    kb = valid_len.loc[keep].astype(float) / 1000.0
    rpk = sub.div(kb, axis=0)
    denom = rpk.sum(axis=0).replace(0.0, np.nan)
    return rpk.div(denom, axis=1) * 1e6


def _fpkm_to_tpm(fpkm: pd.DataFrame) -> pd.DataFrame:
    sub = fpkm.astype(float).clip(lower=0.0)
    denom = sub.sum(axis=0).replace(0.0, np.nan)
    return sub.div(denom, axis=1) * 1e6


def _project_peg_matrix(gene_tpm: pd.DataFrame, groups_df: pd.DataFrame) -> pd.DataFrame:
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


def _or_share_ge5(gene_tpm: pd.DataFrame, or_gids: set[str]) -> pd.Series:
    or_sub = gene_tpm.loc[gene_tpm.index.isin(or_gids)]
    if len(or_sub) < 150:
        return pd.Series(0.0, index=gene_tpm.columns, dtype=float)
    return (or_sub >= 5.0).mean(axis=0)


def load_ensembl98_dog_reference(
    dog_in: Path,
) -> tuple[
    ReferenceIdentity,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, int],
    dict[str, str],
    set[str],
    pd.DataFrame,
    pd.DataFrame,
    list[dict],
]:
    """Verify Ensembl 98 CanFam3.1 files, build PEGs, symbol index, OR set, and 31-mer risk."""
    reg = next(r for r in canine.references.registered_references() if r["id"] == "ensembl98")
    ref = ReferenceIdentity.from_dict(reg["identity"])
    gtf_path = dog_in / "Canis_familiaris.CanFam3.1.98.gtf.gz"
    pep_path = dog_in / "Canis_familiaris.CanFam3.1.pep.all.fa.gz"
    cdna_path = dog_in / "Canis_familiaris.CanFam3.1.cdna.all.fa.gz"
    ncrna_path = dog_in / "Canis_familiaris.CanFam3.1.ncrna.fa.gz"
    ft_path = dog_in / "GCF_011100685.1_UU_Cfam_GSD_1.0_feature_table.txt.gz"

    if _sha256(gtf_path) != ref.content_provenance["gtf"]:
        raise ValueError("CanFam3.1 GTF SHA-256 mismatch")
    if _sha256(pep_path) != ref.content_provenance["protein"]:
        raise ValueError("CanFam3.1 peptide FASTA SHA-256 mismatch")

    source_rows = [
        {
            "source_id": "ensembl98-gtf",
            "pillar": "reference_proteome",
            "filename": gtf_path.name,
            "unit": "annotation",
            "strain": "Boxer Tasha (CanFam3.1)",
            "bytes": gtf_path.stat().st_size,
            "sha256": _sha256(gtf_path),
        },
        {
            "source_id": "ensembl98-protein",
            "pillar": "reference_proteome",
            "filename": pep_path.name,
            "unit": "amino-acid sequence",
            "strain": "Boxer Tasha (CanFam3.1)",
            "bytes": pep_path.stat().st_size,
            "sha256": _sha256(pep_path),
        },
        {
            "source_id": "ensembl98-cdna",
            "pillar": "reference_proteome",
            "filename": cdna_path.name,
            "unit": "transcript sequence",
            "strain": "Boxer Tasha (CanFam3.1)",
            "bytes": cdna_path.stat().st_size,
            "sha256": _sha256(cdna_path),
        },
        {
            "source_id": "ensembl98-ncrna",
            "pillar": "reference_proteome",
            "filename": ncrna_path.name,
            "unit": "ncRNA sequence",
            "strain": "Boxer Tasha (CanFam3.1)",
            "bytes": ncrna_path.stat().st_size,
            "sha256": _sha256(ncrna_path),
        },
        {
            "source_id": "refseq106-feature-table",
            "pillar": "reference_proteome",
            "filename": ft_path.name,
            "unit": "NCBI GeneID crosswalk",
            "strain": "German Shepherd (UU_Cfam_GSD_1.0)",
            "bytes": ft_path.stat().st_size,
            "sha256": _sha256(ft_path),
        },
    ]

    # Parse RefSeq 106 feature_table GeneID -> symbol for unnamed Ensembl 98 genes carrying NCBI Acc
    refseq_gid_to_sym: dict[str, str] = {}
    with gzip.open(ft_path, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            ps = line.rstrip("\n").split("\t")
            if len(ps) > 15 and ps[0] == "gene":
                rs_sym = ps[14].strip()
                rs_gid = ps[15].strip()
                if rs_gid and rs_sym and not rs_sym.startswith("LOC"):
                    refseq_gid_to_sym[rs_gid] = rs_sym

    explicit_unnamed_labels = {
        "ENSCAFG00000045979": "DDX53",
        "ENSCAFG00000015307": "SSX6",
        "ENSCAFG00000031355": "PAGE3-like",
        "ENSCAFG00000018988": "MAGEA13P",
    }

    attr_re = re.compile(r'(\w+)\s+"([^"]+)"')
    gene_meta: dict[str, tuple[str, str, str]] = {}  # gid -> (chrom, strand, raw_sym)
    exons_by_gene: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with gzip.open(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feat, start_s, end_s, _score, strand, _frame, attrs_s = parts
            if feat not in ("gene", "exon"):
                continue
            attrs = dict(attr_re.findall(attrs_s))
            gid = attrs.get("gene_id", "").split(".")[0]
            if not gid:
                continue
            if feat == "gene":
                if attrs.get("gene_biotype") == "protein_coding":
                    sym = attrs.get("gene_name") or explicit_unnamed_labels.get(gid, "")
                    gene_meta[gid] = (chrom, strand, sym)
            elif feat == "exon":
                exons_by_gene[gid].append((int(start_s), int(end_s)))

    exonic_bp: dict[str, int] = {}
    merged_exons: dict[str, list[tuple[int, int]]] = {}
    for gid, ivs in exons_by_gene.items():
        m = _merge_intervals(ivs)
        merged_exons[gid] = m
        exonic_bp[gid] = sum(e - s + 1 for s, e in m)

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

    pep_records: list[dict] = []
    seen_coding_genes: set[str] = set()
    or_gids: set[str] = set()
    ncbi_acc_re = re.compile(r"\[Source:NCBI gene;Acc:(\d+)\]")
    or_sym_re = re.compile(r"^(?:COR\d+|OR\d+[A-Z]+\d+)", re.IGNORECASE)

    def _flush_pep(hdr: str, seq_parts: list[str]) -> None:
        if not hdr:
            return
        seq = "".join(seq_parts).strip()
        tokens = hdr[1:].split()
        pid = tokens[0].split(".")[0]
        kv = {}
        for tok in tokens[1:]:
            if ":" in tok and not tok.startswith("description:"):
                k, v = tok.split(":", 1)
                kv[k] = v
        loc_parts = tokens[2].split(":") if len(tokens) > 2 else []
        chrom = loc_parts[2] if len(loc_parts) >= 3 else ""
        gid = kv.get("gene", "").split(".")[0]
        tid = kv.get("transcript", "").split(".")[0]
        g_biotype = kv.get("gene_biotype", "")
        t_biotype = kv.get("transcript_biotype", "protein_coding")
        if g_biotype != "protein_coding" and gid not in gene_meta:
            return
        desc = hdr.split("description:", 1)[1] if "description:" in hdr else ""
        sym = kv.get("gene_symbol") or gene_meta.get(gid, ("", "", ""))[2]
        if not sym:
            m_acc = ncbi_acc_re.search(desc)
            if m_acc and m_acc.group(1) in refseq_gid_to_sym:
                sym = refseq_gid_to_sym[m_acc.group(1)]
            elif gid in explicit_unnamed_labels:
                sym = explicit_unnamed_labels[gid]
            else:
                sym = gid
        if or_sym_re.match(sym) or "olfactory receptor" in desc.lower():
            or_gids.add(gid)
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
    with gzip.open(pep_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush_pep(cur_hdr, cur_parts)
                cur_hdr = line.rstrip("\n")
                cur_parts = []
            else:
                cur_parts.append(line.strip())
    _flush_pep(cur_hdr, cur_parts)

    for gid, (_chrom, _strand, raw_sym) in sorted(gene_meta.items()):
        if gid not in seen_coding_genes:
            sym = raw_sym or explicit_unnamed_labels.get(gid, gid)
            if or_sym_re.match(sym):
                or_gids.add(gid)
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

    groups_df, seq_df = group_reference_proteins(pep_records, species="dog", reference=ref)
    symbol_by_gene: dict[str, str] = {}
    for r in pep_records:
        if r["gene_id"] not in symbol_by_gene or symbol_by_gene[r["gene_id"]].startswith("ENSCAFG"):
            symbol_by_gene[r["gene_id"]] = r["symbol"]

    # Assign canine home tissues; only dental pulp & taste bud lack quantitative normal TPM panels!
    home_by_eid: dict[str, str] = {}
    unmeasured_by_eid: dict[str, bool] = {}
    for row in groups_df.itertuples(index=False):
        gids = json.loads(row.member_gene_ids)
        syms_up = [symbol_by_gene.get(g, "").upper() for g in gids]
        homes: list[str] = []
        unmeas = False
        for g, su in zip(gids, syms_up):
            if su in DENTAL_HOME_SYMBOLS:
                homes.append("dental pulp")
                unmeas = True
            elif su.startswith(("TAS1R", "TAS2R")) or su in TASTE_HOME_SYMBOLS:
                homes.append("taste bud")
                unmeas = True
            elif su.startswith("VN1R") or re.match(r"^TAAR\d+$", su):
                homes.append("olfactory epithelium")
                unmeas = True
            elif su.startswith("KRTAP") or su in HAIR_FOLLICLE_KERATINS:
                homes.append("hair follicle")
            elif g in or_gids:
                homes.append("olfactory epithelium")
                unmeas = True
        if homes:
            home_by_eid[row.entity_id] = "; ".join(dict.fromkeys(homes))
        unmeasured_by_eid[row.entity_id] = unmeas
    groups_df["home_tissue"] = groups_df["entity_id"].map(home_by_eid).fillna("")
    groups_df["home_tissue_unmeasured"] = (
        groups_df["entity_id"].map(unmeasured_by_eid).fillna(False).astype(bool)
    )

    all_nom_gids = {g for gs in CURATED_NOMINATION_TO_ENSCAFG.values() for g in gs}
    groups_df["published_nomination"] = groups_df["member_gene_ids"].map(
        lambda s: bool(set(json.loads(s)) & all_nom_gids)
    )
    # Flag chimeric readthrough gene models in Ensembl 98 CanFam3.1 (e.g. ENSCAFG00000003760 merging NT5C1B + RDH14)
    chimeric_readthrough_gids = {"ENSCAFG00000003760"}
    groups_df["readthrough_locus"] = groups_df["member_gene_ids"].map(
        lambda s: bool(set(json.loads(s)) & chimeric_readthrough_gids)
    )
    single_exon_gids = {gid for gid in gene_meta if len(merged_exons.get(gid, ())) <= 1}

    # Build canine symbol-index.csv covering all symbols + 14 curated nominations + RefSeq GeneIDs
    sym_idx_rows: list[tuple[str, str]] = []
    for gid, sym in sorted(symbol_by_gene.items()):
        if sym and not sym.startswith("ENSCAFG"):
            sym_idx_rows.append((sym, gid))
    snap = canine.cta_snapshot()
    for cand in snap.candidates():
        cid = cand["id"]
        lbl = cand["label"]
        mapped_ens = CURATED_NOMINATION_TO_ENSCAFG.get(cid, ())
        for eg in mapped_ens:
            sym_idx_rows.append((cid, eg))
            sym_idx_rows.append((lbl, eg))
            for rgid in cand["gene_ids"]:
                sym_idx_rows.append((str(rgid), eg))
    sym_index_df = (
        pd.DataFrame(sym_idx_rows, columns=["symbol", "ensembl_gene_id"])
        .drop_duplicates()
        .sort_values(["symbol", "ensembl_gene_id"])
        .reset_index(drop=True)
    )

    # Two-pass phase-invariant 31-mer paralog check across all Ensembl 98 cDNA + ncRNA transcripts:
    # 1) Sample probe 31-mers (step=4) across protein-coding transcripts per PEG.
    # 2) Scan EVERY transcript in cdna.all.fa.gz + ncrna.fa.gz at step=1 (every base offset)
    #    so out-of-phase 5'-UTR shifts and processed pseudogenes/ncRNAs are caught identically.
    gene_to_entity = {
        g: r.entity_id
        for r in groups_df.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }
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

    with gzip.open(cdna_path, "rt", encoding="utf-8") as fh:
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
    kmer_owners: dict[bytes, str] = {}
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

    for tx_path in (cdna_path, ncrna_path):
        if not tx_path.exists():
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

    risk_rows = []
    for gid in sorted(gene_to_entity):
        seen_for_gene = probe_by_gene.get(gid, set())
        tot = len(seen_for_gene)
        amb = sum(1 for km_b in seen_for_gene if kmer_owners.get(km_b) == MULTI)
        shared_frac = (amb / tot) if tot > 0 else 0.0
        has_ov = gid in same_strand_overlap_genes
        risk_rows.append(
            {
                "gene_id": gid,
                "shared_kmer_fraction": round(shared_frac, 4),
                "same_strand_host_overlap": has_ov,
                "mapping_risk": bool(has_ov or shared_frac >= 0.05),
            }
        )
    return (
        ref,
        groups_df,
        seq_df,
        exonic_bp,
        symbol_by_gene,
        or_gids,
        single_exon_gids,
        pd.DataFrame(risk_rows),
        sym_index_df,
        source_rows,
    )


def load_dog_normal_and_repro_panels(
    scratch_dir: Path,
    dog_in: Path,
    extra_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    symbol_by_gene: dict[str, str],
    or_gids: set[str],
    risk_gids: set[str],
) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, set[str]], dict[str, float], list[dict]
]:
    """Aggregate canine normal N, Disposable D, thymic, and reproductive/embryonic TPM panels."""
    source_rows: list[dict] = []
    lib_tar = scratch_dir / "Canis_lupus_familiaris_RNA-Seq_experiments_libraries.tar.gz"
    source_rows.append(
        {
            "source_id": "bgee-15.2-dog-metadata",
            "pillar": "healthy_baseline_rna",
            "filename": lib_tar.name,
            "unit": "metadata",
            "strain": "multi-breed (explicit per library)",
            "bytes": lib_tar.stat().st_size,
            "sha256": _sha256(lib_tar),
        }
    )
    with tarfile.open(lib_tar) as tf:
        libs_meta = pd.read_csv(
            tf.extractfile("Canis_lupus_familiaris_RNA-Seq_libraries.tsv"), sep="\t"
        )
    for col in libs_meta.columns:
        if libs_meta[col].dtype == object:
            libs_meta[col] = libs_meta[col].astype(str).str.strip('"')

    # Order Bgee experiments with SRP181218 & SRP114662 first so SRA Run ID deduplication is deterministic
    all_exps = [
        e
        for e in libs_meta["Experiment ID"].drop_duplicates().tolist()
        if e not in EXCLUDED_DISEASE_BGEE_EXPS
    ]
    priority = {"SRP181218": 0, "SRP114662": 1, "SRP009687": 2}
    all_exps.sort(key=lambda x: (priority.get(x, 10), x))

    seen_run_keys: set[str] = set()
    somatic_peg_cols: dict[str, dict[str, pd.Series]] = defaultdict(dict)
    disposable_peg_cols: dict[str, dict[str, pd.Series]] = defaultdict(dict)
    thymic_peg_cols: dict[str, pd.Series] = {}
    repro_records: list[tuple[str, str, str, str, str, pd.Series]] = []
    male_skin_cols: list[pd.Series] = []

    for exp in all_exps:
        if exp == "SRP181218":
            fpath = dog_in / "bgee-srp181218.tsv.gz"
        elif exp == "SRP009687":
            fpath = dog_in / "bgee-srp009687-curl.tsv.gz"
        else:
            fpath = dog_in / f"Canis_lupus_familiaris_RNA-Seq_read_counts_TPM_{exp}.tsv.gz"
        if not fpath.is_file():
            continue

        exp_meta = libs_meta[libs_meta["Experiment ID"] == exp]
        keep_libs: dict[str, tuple[str, str, str, str, str]] = {}
        for _, mrow in exp_meta.iterrows():
            lid = str(mrow["Library ID"])
            run_key = str(mrow.get("Run IDs") or lid).strip()
            if run_key in seen_run_keys:
                continue
            anat = str(mrow["Anatomical entity name"])
            st_name = str(mrow["Stage name"])
            st_auth = str(mrow.get("Stage author annotation") or "").strip()
            strain = str(mrow.get("Strain") or "dog").strip()
            sex = str(mrow.get("Sex") or "na").strip()

            # BarkBase gestational d33/d36/d39/d44 libraries in SRP114662 are canine embryonic/fetal!
            if re.match(r"^d\d+$", st_auth, re.IGNORECASE) or any(
                k in st_name.lower() for k in ("embryo", "fetus", "neurula", "organogenesis")
            ):
                seen_run_keys.add(run_key)
                keep_libs[lid] = ("developmental", "embryo", st_auth or st_name, strain, sex)
                continue

            # Exclude neonatal/puppy/juvenile non-adult stages from adult somatic & reproductive panels
            if any(
                k in st_name.lower()
                for k in ("neonate", "neonatal", "puppy", "juvenile", "immature", "newborn")
            ):
                continue

            mapped = BGEE_DOG_ANATOMY_MAP.get(anat)
            if mapped is None:
                continue
            scope, target_ctx = mapped
            seen_run_keys.add(run_key)
            keep_libs[lid] = (scope, target_ctx, st_name, strain, sex)

        if not keep_libs:
            continue

        source_rows.append(
            {
                "source_id": f"bgee-dog-{exp.lower()}",
                "pillar": "reproductive_and_developmental_rna"
                if exp in {"SRP071663", "ERP019830", "SRP024250", "SRP061328"}
                else "healthy_baseline_rna",
                "filename": fpath.name,
                "unit": "TPM",
                "strain": "recorded per library",
                "bytes": fpath.stat().st_size,
                "sha256": _sha256(fpath),
            }
        )

        chunks = []
        for chunk in pd.read_csv(
            fpath,
            sep="\t",
            usecols=["Library ID", "Gene ID", "TPM"],
            dtype={"Library ID": str, "Gene ID": str, "TPM": float},
            chunksize=2_000_000,
        ):
            chunk["Library ID"] = chunk["Library ID"].str.strip('"')
            sub_c = chunk[chunk["Library ID"].isin(keep_libs)].copy()
            if not sub_c.empty:
                sub_c["Gene ID"] = sub_c["Gene ID"].str.strip('"').str.split(".").str[0]
                chunks.append(sub_c)
        if not chunks:
            continue
        piv = pd.concat(chunks, ignore_index=True).pivot(
            index="Gene ID", columns="Library ID", values="TPM"
        )
        or_sh = _or_share_ge5(piv, or_gids)
        qc_cols = [
            lid
            for lid, (_sc, t_ctx, *_rest) in keep_libs.items()
            if lid in piv.columns
            and (t_ctx in {"olfactory epithelium", "testis"} or float(or_sh.get(lid, 0.0)) <= 0.025)
        ]
        if not qc_cols:
            continue
        peg_mat = _project_peg_matrix(piv[qc_cols], groups_df)
        for lid in qc_cols:
            scope, target_ctx, _st_lbl, strain, sex = keep_libs[lid]
            col = peg_mat[lid]
            donor_key = f"bgee:{exp}:{strain}:{sex}:{lid}"
            if scope == "normal":
                somatic_peg_cols[target_ctx][donor_key] = col
                if target_ctx == "skin" and sex.lower() == "male":
                    male_skin_cols.append(col)
            elif scope == "disposable":
                disposable_peg_cols[target_ctx][donor_key] = col
            elif scope == "thymic":
                thymic_peg_cols[donor_key] = col
            elif scope in ("reproductive", "developmental"):
                repro_records.append(
                    (f"bgee-dog-{exp.lower()}", lid, donor_key, scope, target_ctx, col)
                )

    # 2. GSE219045 300-sample sex-stratified canine organ atlas (FPKM -> exact within-sample linear TPM)
    g219_path = dog_in / "GSE219045_Dog_FPKM_table_Jun22.xlsx"
    source_rows.append(
        {
            "source_id": "gse219045-dog-atlas",
            "pillar": "healthy_baseline_rna",
            "filename": g219_path.name,
            "unit": "TPM (normalized from within-sample FPKM)",
            "strain": "Beagle sex-stratified reference",
            "bytes": g219_path.stat().st_size,
            "sha256": _sha256(g219_path),
        }
    )
    g219_df = pd.read_excel(g219_path, sheet_name="FPKM_table", engine="openpyxl")
    g219_df["GeneID"] = g219_df["GeneID"].astype(str).str.split(".").str[0]
    g219_df = g219_df.groupby("GeneID", sort=False).sum(numeric_only=True)
    g219_tpm = _fpkm_to_tpm(g219_df)
    g219_or = _or_share_ge5(g219_tpm, or_gids)
    g219_peg = _project_peg_matrix(g219_tpm, groups_df)

    # Group multiple sub-regions from the same Beagle donor within each somatic organ by donor median
    g219_somatic_by_dog: dict[tuple[str, str], list[pd.Series]] = defaultdict(list)
    g219_disp_by_dog: dict[tuple[str, str], list[pd.Series]] = defaultdict(list)
    fem_pat = re.compile(r"^Dog_Female_([A-Za-z]+)_(\d+)$")
    male_pat = re.compile(r"^Dog_([A-Za-z]+)_Male_(\d+)$")

    for col_name in g219_peg.columns:
        mf = fem_pat.match(col_name)
        mm = male_pat.match(col_name)
        if mf:
            organ_tag, dog_id, sex_tag = mf.group(1), mf.group(2), "F"
        elif mm:
            organ_tag, dog_id, sex_tag = mm.group(1), mm.group(2), "M"
        else:
            continue
        mapped = GSE219045_DOG_ORGAN_MAP.get(organ_tag)
        if mapped is None:
            continue
        scope, target_ctx = mapped
        if (
            target_ctx not in {"olfactory epithelium", "testis"}
            and float(g219_or.get(col_name, 0.0)) > 0.025
        ):
            continue
        col_ser = g219_peg[col_name]
        donor_id = f"gse219045:{sex_tag}:{dog_id}"
        if scope == "normal":
            g219_somatic_by_dog[(target_ctx, donor_id)].append(col_ser)
        elif scope == "disposable":
            g219_disp_by_dog[(target_ctx, donor_id)].append(col_ser)
        elif scope == "reproductive":
            if target_ctx == "uterus":
                # Full-thickness canine bicornuate uterine necropsy sections span uterus, uterine_horn, & endometrium
                for u_comp in ("uterus", "uterine_horn", "endometrium"):
                    repro_records.append(
                        ("gse219045-dog-atlas", col_name, donor_id, "reproductive", u_comp, col_ser)
                    )
            else:
                repro_records.append(
                    ("gse219045-dog-atlas", col_name, donor_id, "reproductive", target_ctx, col_ser)
                )

    for (ctx, donor_id), sers in g219_somatic_by_dog.items():
        somatic_peg_cols[ctx][donor_id] = (
            sers[0] if len(sers) == 1 else pd.DataFrame(sers).median(axis=0)
        )
    for (ctx, donor_id), sers in g219_disp_by_dog.items():
        disposable_peg_cols[ctx][donor_id] = (
            sers[0] if len(sers) == 1 else pd.DataFrame(sers).median(axis=0)
        )

    # 3. GSE119810 normal female canine mammary biopsies (64 samples: resting vs milk-operon secretory/lactating)
    p119 = dog_in / "GSE119810_CMT_222S_FPKM.csv.gz"
    source_rows.append(
        {
            "source_id": "gse119810-normal-mammary",
            "pillar": "reproductive_and_developmental_rna",
            "filename": p119.name,
            "unit": "TPM (normalized from native ENSCAFG FPKM)",
            "strain": "companion dogs (64 normal female mammary glands)",
            "bytes": p119.stat().st_size,
            "sha256": _sha256(p119),
        }
    )
    df119 = pd.read_csv(p119, index_col=0)
    df119.index = df119.index.astype(str).str.split(".").str[0]
    norm_mam_cols = [c for c in df119.columns if str(c).endswith("-normal")]
    tpm119 = _fpkm_to_tpm(df119[norm_mam_cols])
    or119 = _or_share_ge5(tpm119, or_gids)
    qc_mam_cols = [c for c in norm_mam_cols if float(or119.get(c, 0.0)) <= 0.025]
    milk_gids = [
        g
        for g in (
            "ENSCAFG00000008950",  # LALBA
            "ENSCAFG00000002881",  # CSN2
            "ENSCAFG00000030389",  # CSN1S1
            "ENSCAFG00000002889",  # CSN3
        )
        if g in tpm119.index
    ]
    milk_sub = tpm119.loc[milk_gids, qc_mam_cols]
    # Mask paralog-ambiguous / host-overlapped genes on Cufflinks multi-mapping FPKM (GSE119810),
    # matching GSE107746 Cufflinks FPKM qualification while Bgee/GSE219045 supply unique resting_mammary
    tpm119_clean = tpm119[qc_mam_cols].copy()
    tpm119_clean.loc[tpm119_clean.index.isin(risk_gids)] = np.nan
    peg119 = _project_peg_matrix(tpm119_clean, groups_df)
    for c in qc_mam_cols:
        is_lact = int((milk_sub[c] >= 5.0).sum()) >= 3 and float(milk_sub[c].max()) >= 15.0
        m_comp = "lactating_mammary" if is_lact else "resting_mammary"
        dog_id = c.removesuffix("-normal")
        repro_records.append(
            (
                "gse119810-normal-mammary",
                c,
                f"gse119810:{dog_id}",
                "reproductive",
                m_comp,
                peg119[c],
            )
        )

    male_skin_p95: dict[str, float] = {}
    if male_skin_cols:
        ms_mat = pd.DataFrame(male_skin_cols)
        ms_q = target_selection.summarize_normal_tissue_tpm(ms_mat, method="p95")
        male_skin_p95 = {str(k): float(v) for k, v in ms_q.items() if pd.notna(v)}

    # 4. Supplementary companion-dog normal organ & primary macrophage RNA controls (<=2.5% OR gDNA QC)
    exonic_series_norm = pd.Series(exonic_bp, dtype=float)
    extra_immune_gate_rows: list[dict] = []
    if extra_dir.is_dir():
        # 4a. GSE196108 healthy canine adrenal cortex (NAC1..NAC10)
        p196_n = extra_dir / "GSE196108_Sanders_raw_counts.txt.gz"
        if p196_n.is_file():
            df196_n = pd.read_csv(p196_n, sep="\t")
            df196_n["gene"] = df196_n["gene"].astype(str).str.split(".").str[0]
            df196_n = df196_n.drop_duplicates("gene").set_index("gene")
            nac_cols = [c for c in df196_n.columns if c.startswith("NAC")]
            tpm196_n = _counts_to_tpm(df196_n[nac_cols], exonic_series_norm)
            or196_n = _or_share_ge5(tpm196_n, or_gids)
            qc_nac = [c for c in nac_cols if float(or196_n.get(c, 0.0)) <= 0.025]
            if qc_nac:
                peg196_n = _project_peg_matrix(tpm196_n[qc_nac], groups_df)
                for c in qc_nac:
                    somatic_peg_cols["adrenal gland"][f"gse196108:{c}"] = peg196_n[c]

        # 4b. GSE345838 normal canine peripheral blood (BL1..BL10) & bone marrow (BM1..BM5)
        p345_n = extra_dir / "GSE345838_Canine_Acute_Leukemia_featurecounts_raw_counts.txt.gz"
        if p345_n.is_file():
            df345_n = pd.read_csv(p345_n, sep="\t", comment="#")
            df345_n["Geneid"] = df345_n["Geneid"].astype(str).str.split(".").str[0]
            df345_n = df345_n.drop_duplicates("Geneid").set_index("Geneid")
            bl_bm_cols = [c for c in df345_n.columns if c.startswith(("BL", "BM"))]
            tpm345_n = _counts_to_tpm(df345_n[bl_bm_cols], df345_n["Length"].astype(float))
            or345_n = _or_share_ge5(tpm345_n, or_gids)
            qc_bl_bm = [c for c in bl_bm_cols if float(or345_n.get(c, 0.0)) <= 0.025]
            if qc_bl_bm:
                peg345_n = _project_peg_matrix(tpm345_n[qc_bl_bm], groups_df)
                for c in qc_bl_bm:
                    ctx_345 = "blood" if c.startswith("BL") else "bone marrow"
                    somatic_peg_cols[ctx_345][f"gse345838:{c}"] = peg345_n[c]

        # 4c. GSE317244 normal canine gingiva (Gingiva-*, buccal mucosa)
        p317_n = extra_dir / "GSE317244_DS1_rawCounts.txt.gz"
        if p317_n.is_file():
            df317_n = pd.read_csv(p317_n, sep="\t", index_col=0)
            df317_n.index = df317_n.index.astype(str).str.split(".").str[0]
            df317_n = df317_n[~df317_n.index.duplicated(keep="first")]
            ging_cols = [c for c in df317_n.columns if str(c).startswith("Gingiva-")]
            tpm317_n = _counts_to_tpm(df317_n[ging_cols], exonic_series_norm)
            or317_n = _or_share_ge5(tpm317_n, or_gids)
            qc_ging = [c for c in ging_cols if float(or317_n.get(c, 0.0)) <= 0.025]
            if qc_ging:
                peg317_n = _project_peg_matrix(tpm317_n[qc_ging], groups_df)
                for c in qc_ging:
                    somatic_peg_cols["buccal mucosa"][f"gse317244:{c}"] = peg317_n[c]

        # 4d. GSE288068 unaffected canine blood, liver, lung, spleen controls
        p288_n = extra_dir / "GSE288068_HS_MSU_RNAseq_Analysis_rawCounts.txt.gz"
        if p288_n.is_file():
            df288_n = pd.read_csv(p288_n, sep="\t")
            df288_n["Symbol"] = df288_n["Symbol"].astype(str).str.split(".").str[0]
            df288_n = df288_n.drop_duplicates("Symbol").set_index("Symbol")
            norm_288_map = {
                "HT5GWDSXX_18927780_S54": "blood",
                "HT5GWDSXX_18927782_S66": "blood",
                "HT5GWDSXX_18927808_S53": "blood",
                "HVTWLDSXX_18937521_S79": "liver",
                "HVTWLDSXX_18937523_S80": "lung",
                "HVTWLDSXX_18936778_S54": "spleen",
            }
            n288_cols = [c for c in norm_288_map if c in df288_n.columns]
            tpm288_n = _counts_to_tpm(df288_n[n288_cols], exonic_series_norm)
            or288_n = _or_share_ge5(tpm288_n, or_gids)
            qc_288 = [c for c in n288_cols if float(or288_n.get(c, 0.0)) <= 0.025]
            if qc_288:
                peg288_n = _project_peg_matrix(tpm288_n[qc_288], groups_df)
                for c in qc_288:
                    somatic_peg_cols[norm_288_map[c]][f"gse288068:{c}"] = peg288_n[c]

        # 4e. GSE160124 fresh-frozen healthy canine brain controls (*_FF_control.bam)
        p160_n = extra_dir / "GSE160124_FCresult.txt.gz"
        if p160_n.is_file():
            df160_n = pd.read_csv(p160_n, sep="\t", comment="#")
            df160_n["Geneid"] = df160_n["Geneid"].astype(str).str.split(".").str[0]
            df160_n = df160_n.drop_duplicates("Geneid").set_index("Geneid")
            ff_ctrl = [c for c in df160_n.columns if str(c).endswith("_FF_control.bam")]
            tpm160_n = _counts_to_tpm(df160_n[ff_ctrl], df160_n["Length"].astype(float))
            or160_n = _or_share_ge5(tpm160_n, or_gids)
            qc_160 = [c for c in ff_ctrl if float(or160_n.get(c, 0.0)) <= 0.025]
            if qc_160:
                peg160_n = _project_peg_matrix(tpm160_n[qc_160], groups_df)
                for c in qc_160:
                    somatic_peg_cols["brain"][f"gse160124:{c.split('_')[0]}"] = peg160_n[c]

        # 4f. GSE283373 mock primary canine alveolar macrophages (Protected <= 8 TPM replicated P95)
        p283 = extra_dir / "GSE283373_RawCountsAnnotated.xlsx"
        if p283.is_file():
            source_rows.append(
                {
                    "source_id": "gse283373-alveolar-macrophage",
                    "pillar": "normal_bulk_and_immune_safety",
                    "filename": p283.name,
                    "unit": "TPM (from primary canine alveolar macrophage raw counts)",
                    "strain": "healthy canine mock primary alveolar macrophages (n=4)",
                    "bytes": p283.stat().st_size,
                    "sha256": _sha256(p283),
                }
            )
            sym_cands_n: dict[str, list[str]] = defaultdict(list)
            for gid_s, sym_s in symbol_by_gene.items():
                if sym_s and not sym_s.startswith("ENSCAFG"):
                    sym_cands_n[sym_s.upper()].append(gid_s)
            sym_to_ens_n = {k: v[0] for k, v in sym_cands_n.items() if len(v) == 1}
            df283 = pd.read_excel(p283, sheet_name="Sheet 1", engine="openpyxl")
            df283["_gid"] = df283["SYMBOL"].astype(str).str.upper().map(sym_to_ens_n)
            df283 = df283.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
            am_cols = ["140_non", "43_non", "45_non", "48_non"]
            tpm283 = _counts_to_tpm(df283[am_cols].astype(float), exonic_series_norm)
            or283 = _or_share_ge5(tpm283, or_gids)
            qc_am = [c for c in am_cols if float(or283.get(c, 0.0)) <= 0.025]
            if len(qc_am) >= 2:
                peg283 = _project_peg_matrix(tpm283[qc_am], groups_df)
                arr283 = np.sort(
                    np.where(
                        np.isnan(peg283.to_numpy(dtype=float)),
                        -np.inf,
                        peg283.to_numpy(dtype=float),
                    ),
                    axis=1,
                )
                sec_hi = arr283[:, -2]
                p95_283 = peg283.quantile(0.95, axis=1)
                for ri, eid in enumerate(peg283.index):
                    if np.isfinite(sec_hi[ri]) and pd.notna(p95_283.iat[ri]):
                        bound_v = min(float(p95_283.iat[ri]), float(sec_hi[ri]))
                        if bound_v > 8.0:
                            extra_immune_gate_rows.append(
                                {
                                    "entity_id": str(eid),
                                    "layer": "Canine sorted primary macrophage RNA",
                                    "unsafe": True,
                                    "hold": False,
                                    "reason": f"alveolar macrophage (protected immune): {bound_v:.4g} TPM > 8",
                                }
                            )

    # Summarize somatic normal contexts by donor P95 -> comparison_normals + normal_gates
    normal_rows = []
    rna_gate_rows = []
    confirmed_organs: dict[str, set[str]] = defaultdict(set)

    for ctx, donor_series in sorted(somatic_peg_cols.items()):
        mat = pd.DataFrame(donor_series).T
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
                    "source_id": "dog_somatic_panel",
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
                        "layer": "Canine somatic bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{ctx} (vital): {v:.4g} TPM > 4",
                    }
                )
            elif tier == "protected" and v > 8.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Canine somatic bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{ctx} (protected): {v:.4g} TPM > 8",
                    }
                )
            elif tier == "manageable" and v > 128.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Canine somatic bulk RNA",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{ctx} (manageable): {v:.4g} TPM > 128",
                    }
                )
            elif tier == "manageable" and v > 16.0:
                rna_gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Canine somatic bulk RNA",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"{ctx} (manageable hold): {v:.4g} TPM > 16",
                    }
                )

    for disp_ctx, donor_series in sorted(disposable_peg_cols.items()):
        dmat = pd.DataFrame(donor_series).T
        dsum = target_selection.summarize_normal_tissue_tpm(dmat, method="p95")
        for eid, val in dsum.items():
            if pd.notna(val):
                normal_rows.append(
                    {
                        "entity_id": eid,
                        "source_id": "dog_somatic_panel",
                        "context_id": disp_ctx,
                        "tissue": disp_ctx,
                        "scope": "disposable",
                        "value": float(val),
                        "comparison_id": COMPARISON_ID,
                        "unit": "TPM",
                        "qc_pass": True,
                        "identity_resolved": True,
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
                        "source_id": "dog_somatic_panel",
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

    repro_by_comp: dict[tuple[str, str], list[tuple[str, str, str, pd.Series]]] = defaultdict(list)
    for src_id, smp_id, donor_id, scope, comp, ser in repro_records:
        repro_by_comp[(scope, comp)].append((src_id, smp_id, donor_id, ser))

    repro_rows = []
    for (scope, comp), items in sorted(repro_by_comp.items()):
        mat = pd.DataFrame({f"{i}": ser for i, (_src, _smp, _dnr, ser) in enumerate(items)})
        valid_any = mat.notna().any(axis=1)
        sub_mat = mat.loc[valid_any]
        max_vals = sub_mat.max(axis=1)
        best_col = sub_mat.idxmax(axis=1)
        # For large mammary/testis panels (>15 samples), use P99 for HostContext while keeping sample max for RA!
        he_vals = sub_mat.quantile(0.99, axis=1) if sub_mat.shape[1] >= 15 else max_vals
        for eid, val in max_vals.items():
            idx_i = int(best_col.at[eid])
            src_id, smp_id, donor_id, _ = items[idx_i]
            h_v = float(he_vals.at[eid])
            rec = {
                "entity_id": eid,
                "source_id": src_id,
                "sample_id": smp_id,
                "participant_id": donor_id,
                "scope": scope,
                "compartment": comp,
                "tissue": comp.replace("_", " "),
                "value": float(val),
                "host_max_tpm": h_v,
                "comparison_id": COMPARISON_ID,
                "unit": "TPM",
                "qc_pass": True,
                "identity_resolved": True,
            }
            if comp == "resting_mammary":
                # Female mammary biopsies only; male canine mammary RNA is unmeasured (fail-closed None)
                rec["female_max_tpm"] = h_v
            if scope == "developmental":
                rec["cells"] = 100
                rec["positive_cells"] = 100
                rec["developmental_stage"] = "embryonic"
                rec["source_context"] = "embryo"
            repro_rows.append(rec)

    return (
        pd.DataFrame(normal_rows),
        pd.DataFrame(repro_rows),
        pd.DataFrame(
            rna_gate_rows + extra_immune_gate_rows,
            columns=["entity_id", "layer", "unsafe", "hold", "reason"],
        ),
        dict(confirmed_organs),
        male_skin_p95,
        source_rows,
    )


def load_dog_healthy_pbmc_scrna_umi(
    extra_dir: Path,
    coding_genes: set[str],
) -> tuple[pd.DataFrame, list[dict]]:
    """Parse GSE225599 10x single-cell UMI matrices across 7 healthy companion dogs (H_1..H_7).

    Evaluates donor-resolved circulating T/NK lymphocytes, monocytes, neutrophils (Protected)
    and platelets (Vital) on native CanFam3.1 ENSCAFG IDs through assess_rare_cell_rna,
    excluding the 10 osteosarcoma dogs (OS_1..OS_10) and Manageable B lymphocytes.
    """
    tar_path = extra_dir / "GSE225599_RAW.tar"
    if not tar_path.is_file():
        return pd.DataFrame(), []
    source_rows = [
        {
            "source_id": "gse225599-healthy-dog-pbmc-10x",
            "pillar": "normal_cell_and_rare_umi",
            "filename": tar_path.name,
            "unit": "UMI counts",
            "strain": "healthy companion dogs (7 donors H_1..H_7 circulating leukocytes)",
            "bytes": tar_path.stat().st_size,
            "sha256": _sha256(tar_path),
        }
    ]

    tnk_gids = {
        "ENSCAFG00000012802",  # CD3E
        "ENSCAFG00000031719",  # CD3D
        "ENSCAFG00000041713",  # CD4
        "ENSCAFG00000007464",  # CD8A
        "ENSCAFG00000023787",  # NCR1
        "ENSCAFG00000014071",  # PRF1
        "ENSCAFG00000025287",  # GZMB
    }
    b_gids = {
        "ENSCAFG00000010381",  # MS4A1
        "ENSCAFG00000017303",  # CD19
        "ENSCAFG00000004980",  # CD79A
    }
    mono_gids = {
        "ENSCAFG00000000426",  # LYZ
        "ENSCAFG00000018219",  # CSF1R
    }
    neut_gids = {
        "ENSCAFG00000017557",  # S100A8
        "ENSCAFG00000023324",  # S100A12
    }
    plt_gids = {
        "ENSCAFG00000031174",  # PPBP
        "ENSCAFG00000014145",  # ITGA2B
    }

    pop_n_cells: dict[tuple[str, str, str, str], int] = defaultdict(int)
    pos_counts: dict[tuple[str, str, str, str, str], int] = defaultdict(int)
    multi_counts: dict[tuple[str, str, str, str, str], int] = defaultdict(int)
    assayed_genes: set[str] = set()

    from scipy.io import mmread

    with tarfile.open(tar_path) as tf:
        for i in range(7):
            donor_id = f"H_{i + 1}"
            prefix = f"GSM70509{8 + i:02d}_H_{i + 1}"
            with gzip.open(
                tf.extractfile(f"{prefix}_features.tsv.gz"), "rt", encoding="utf-8"
            ) as ffh:
                feat_gids = [ln.split("\t", 1)[0].split(".")[0] for ln in ffh]
            with gzip.open(tf.extractfile(f"{prefix}_matrix.mtx.gz"), "rb") as mfh:
                mat = mmread(mfh).tocsc()

            n_genes_per_cell = np.diff(mat.indptr)
            tot_umi_per_cell = np.asarray(mat.sum(axis=0)).ravel()
            gid_to_row = {g: r_i for r_i, g in enumerate(feat_gids)}
            for g in feat_gids:
                if g in coding_genes:
                    assayed_genes.add(g)

            def _score(
                marker_set: set[str], _g2r: dict[str, int] = gid_to_row, _m=mat
            ) -> np.ndarray:
                r_idxs = [_g2r[g] for g in marker_set if g in _g2r]
                if not r_idxs:
                    return np.zeros(_m.shape[1], dtype=np.float64)
                return np.asarray(_m[r_idxs, :].sum(axis=0)).ravel()

            s_tnk = _score(tnk_gids)
            s_b = _score(b_gids)
            s_mono = _score(mono_gids)
            s_neut = _score(neut_gids)
            s_plt = _score(plt_gids)

            # Standard droplet QC (>=200 genes & >=500 UMIs for nucleated leukocytes; >=50 genes & >=3 marker UMIs for anucleate platelets)
            qc_nuc = (n_genes_per_cell >= 200) & (tot_umi_per_cell >= 500)
            is_b = qc_nuc & (s_b >= 2) & (s_b >= s_tnk) & (s_b >= s_mono)
            is_plt = (
                (s_plt >= 3)
                & (n_genes_per_cell >= 50)
                & (n_genes_per_cell < 1200)
                & (s_tnk == 0)
                & (s_b == 0)
            )
            is_neut = (
                qc_nuc
                & (~is_b)
                & (~is_plt)
                & (s_neut >= 5)
                & (s_neut > 2.0 * s_mono)
                & (s_tnk <= 1)
            )
            is_mono = (
                qc_nuc
                & (~is_b)
                & (~is_plt)
                & (~is_neut)
                & ((s_mono >= 2) | (s_neut >= 2))
                & (s_mono + s_neut > s_tnk)
            )
            is_tnk = qc_nuc & (~is_b) & (~is_plt) & (~is_neut) & (~is_mono) & (s_tnk >= 1)

            lineage_masks = [
                ("T/NK lymphocyte", "protected", is_tnk),
                ("monocyte", "protected", is_mono),
                ("neutrophil", "protected", is_neut),
                ("platelet", "vital", is_plt),
            ]
            mat_csr = mat.tocsr()
            coding_rows = [(r_i, g) for r_i, g in enumerate(feat_gids) if g in coding_genes]
            for cell_type, tier, mask in lineage_masks:
                n_c = int(mask.sum())
                if n_c < 10:
                    continue
                pop = (donor_id, "blood", tier, cell_type)
                pop_n_cells[pop] = n_c
                sub_csr = mat_csr[:, mask]
                indptr = sub_csr.indptr
                data = sub_csr.data
                for r_i, gid in coding_rows:
                    s_idx, e_idx = indptr[r_i], indptr[r_i + 1]
                    if e_idx == s_idx:
                        continue
                    row_vals = data[s_idx:e_idx]
                    n_pos = int((row_vals >= 1).sum())
                    if n_pos > 0:
                        n_mul = int((row_vals >= 2).sum())
                        pos_counts[(gid, *pop)] = n_pos
                        if n_mul > 0:
                            multi_counts[(gid, *pop)] = n_mul

    valid_pops = {pop: n for pop, n in pop_n_cells.items() if n >= 10}
    if not valid_pops:
        return pd.DataFrame(), source_rows
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
                    "source_id": "gse225599_healthy_dog_pbmc",
                    "donor_id": donor_id,
                    "tissue": organ,
                    "cell_type": cell_type,
                    "tier": tier,
                    "n_cells": n_cells,
                    "n_positive": n_pos,
                    "n_multi": n_mul,
                }
            )
    s_donor, s_organ, s_tier, s_ct = sample_pop
    for gid in sorted(assayed_genes - genes_with_multi):
        umi_rows.append(
            {
                "gene_id": gid,
                "source_id": "gse225599_healthy_dog_pbmc",
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


def load_pfizer_dog_proteomics_gates(
    dog_in: Path,
    groups_df: pd.DataFrame,
    seq_df: pd.DataFrame,
    confirmed_organs: dict[str, set[str]],
) -> tuple[pd.DataFrame, dict[tuple[str, str], str], list[dict]]:
    """Evaluate Pfizer 38-tissue canine TMT proteomics (MSV000099744) with within-plex bleed QC."""
    abund_path = dog_in / "dog-protein-abundance.tsv"
    fasta_path = dog_in / "dog-protein-search.fasta"

    dog_assets = {a.id: a for a in canine.expression.assets() + canine.references.assets()}
    if _sha256(abund_path) != dog_assets["pfizer-dog-protein-matrix"].sha256:
        raise ValueError("Pfizer dog-protein-abundance.tsv SHA-256 mismatch")
    if _sha256(fasta_path) != dog_assets["pfizer-dog-search-fasta"].sha256:
        raise ValueError("Pfizer dog-protein-search.fasta SHA-256 mismatch")

    source_rows = [
        {
            "source_id": "pfizer-dog-protein",
            "pillar": "normal_protein",
            "filename": abund_path.name,
            "unit": "within-channel molar ppm (from rIBAQ)",
            "strain": "Beagle (3 male + 3 female pooled per tissue)",
            "bytes": abund_path.stat().st_size,
            "sha256": _sha256(abund_path),
        },
        {
            "source_id": "pfizer-dog-search-fasta",
            "pillar": "normal_protein",
            "filename": fasta_path.name,
            "unit": "searched UniProt canine FASTA",
            "strain": "Canis lupus familiaris (UP000002254)",
            "bytes": fasta_path.stat().st_size,
            "sha256": _sha256(fasta_path),
        },
    ]

    # Index CanFam3.1 sequences and unique single-gene symbols
    seq_to_eids: dict[str, set[str]] = defaultdict(set)
    ninemers_by_eid: dict[str, set[str]] = defaultdict(set)
    for r in seq_df.itertuples(index=False):
        s = str(r.sequence)
        if s:
            seq_to_eids[s].add(str(r.entity_id))
            if len(s) >= 9:
                ninemers_by_eid[str(r.entity_id)].update(s[i : i + 9] for i in range(len(s) - 8))

    sym_to_single_eid: dict[str, list[str]] = defaultdict(list)
    for r in groups_df.itertuples(index=False):
        if int(r.member_count) == 1 and "/" not in str(r.protein_group_label):
            sym_to_single_eid[str(r.protein_group_label).upper()].append(str(r.entity_id))
    unique_sym_eid = {k: v[0] for k, v in sym_to_single_eid.items() if len(v) == 1}

    # Parse dog-protein-search.fasta -> map UniProt accession to unique CanFam3.1 PEG
    gn_re = re.compile(r"\bGN=([^\s]+)")
    acc_to_eid: dict[str, str] = {}

    def _flush_pf(hdr: str, parts: list[str]) -> None:
        if not hdr:
            return
        seq = "".join(parts).strip()
        tok0 = hdr[1:].split()[0]
        acc = tok0.split("|")[1] if "|" in tok0 else tok0
        exact = seq_to_eids.get(seq, set())
        if len(exact) == 1:
            acc_to_eid[acc] = next(iter(exact))
            return
        m_gn = gn_re.search(hdr)
        if m_gn and len(seq) >= 30:
            gn = m_gn.group(1).upper()
            cand_eid = unique_sym_eid.get(gn)
            if cand_eid and cand_eid in ninemers_by_eid:
                q9 = {seq[i : i + 9] for i in range(len(seq) - 8)}
                denom_9 = min(len(q9), len(ninemers_by_eid[cand_eid]))
                if denom_9 > 0 and (len(q9 & ninemers_by_eid[cand_eid]) / denom_9) >= 0.80:
                    acc_to_eid[acc] = cand_eid

    f_hdr = ""
    f_parts: list[str] = []
    with fasta_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush_pf(f_hdr, f_parts)
                f_hdr = line.rstrip("\n")
                f_parts = []
            else:
                f_parts.append(line.strip())
    _flush_pf(f_hdr, f_parts)

    df = pd.read_csv(
        abund_path,
        sep="\t",
        usecols=["Protein", "Condition", "Mixture", "rIBAQ"],
        low_memory=False,
    )
    df["rIBAQ"] = pd.to_numeric(df["rIBAQ"], errors="coerce").fillna(0.0).clip(lower=0.0)
    # Channel totals use ALL quantified proteins per (Mixture, Condition) for true molar ppm
    chan_totals = (
        df.groupby(["Mixture", "Condition"])["rIBAQ"].transform("sum").replace(0.0, np.nan)
    )
    df["ppm"] = (df["rIBAQ"] / chan_totals) * 1e6

    # Map unambiguous single-PEG protein groups
    def _resolve_pg(pg_str: str) -> str | None:
        accs = [x.split("|")[1] if "|" in x else x for x in str(pg_str).split(";") if x.strip()]
        eids = {acc_to_eid[a] for a in accs if a in acc_to_eid}
        return next(iter(eids)) if len(eids) == 1 else None

    df["entity_id"] = [_resolve_pg(p) for p in df["Protein"]]
    mapped_df = df.dropna(subset=["entity_id"]).copy()
    agg = mapped_df.groupby(["entity_id", "Mixture", "Condition"], as_index=False)["ppm"].sum()

    repro_conds = {c for c, (sc, _) in PFIZER_CONDITION_MAP.items() if sc == "reproductive"}
    repro_by_eid_mix: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for rr in agg[agg["Condition"].isin(repro_conds)].itertuples(index=False):
        repro_by_eid_mix[(str(rr.entity_id), str(rr.Mixture))][str(rr.Condition)] = float(rr.ppm)
    repro_plex_max = {k: max(v.values()) for k, v in repro_by_eid_mix.items() if v}

    gate_rows: list[dict] = []
    repro_prot_status: dict[tuple[str, str], str] = {}
    somatic_max_by_eid_organ: dict[tuple[str, str], float] = defaultdict(float)

    for r in agg.itertuples(index=False):
        eid = str(r.entity_id)
        cond = str(r.Condition)
        mapped = PFIZER_CONDITION_MAP.get(cond)
        if mapped is None:
            continue
        scope, organ = mapped
        raw_ppm = float(r.ppm)
        if scope == "reproductive":
            # Subtract cross-organ reproductive TMT bleed within the same plex (e.g. Testis -> Ovary in Pool1)
            oth_vals = [
                v for c2, v in repro_by_eid_mix.get((eid, str(r.Mixture)), {}).items() if c2 != cond
            ]
            oth_max = max(oth_vals) if oth_vals else 0.0
            clean_repro_ppm = max(0.0, raw_ppm - 0.10 * oth_max)
            if clean_repro_ppm >= 10.0:
                repro_prot_status[(eid, organ)] = "positive"
            continue
        if scope != "normal":
            continue
        plex_repro = float(repro_plex_max.get((eid, r.Mixture), 0.0))
        clean_ppm = max(0.0, raw_ppm - 0.10 * plex_repro)
        if clean_ppm > somatic_max_by_eid_organ[(eid, organ)]:
            somatic_max_by_eid_organ[(eid, organ)] = clean_ppm

    for (eid, organ), ppm in sorted(somatic_max_by_eid_organ.items()):
        if ppm >= 2.0:
            confirmed_organs.setdefault(eid, set()).add(organ)
        tier = SOMATIC_CONTEXT_TIER[organ]
        if tier == "vital":
            if ppm > 25.0:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Pfizer canine normal proteomics",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{organ} (vital pooled protein): {ppm:.1f} ppm > 25",
                    }
                )
            elif ppm > 10.0:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Pfizer canine normal proteomics review",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"{organ} (vital pooled protein hold): {ppm:.1f} ppm > 10",
                    }
                )
        elif tier == "protected":
            if ppm > 50.0:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Pfizer canine normal proteomics",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{organ} (protected pooled protein): {ppm:.1f} ppm > 50",
                    }
                )
            elif ppm > 20.0:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Pfizer canine normal proteomics review",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"{organ} (protected pooled protein hold): {ppm:.1f} ppm > 20",
                    }
                )
        elif tier == "manageable":
            if ppm > 250.0:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Pfizer canine normal proteomics",
                        "unsafe": True,
                        "hold": False,
                        "reason": f"{organ} (manageable pooled protein): {ppm:.1f} ppm > 250",
                    }
                )
            elif ppm > 50.0:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Pfizer canine normal proteomics review",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"{organ} (manageable pooled protein hold): {ppm:.1f} ppm > 50",
                    }
                )

    # Preserve curated somatic RT-PCR counterevidence holds from nemec2019 / hindriks2025
    gene_to_eid = {
        g: r.entity_id
        for r in groups_df.itertuples(index=False)
        for g in json.loads(r.member_gene_ids)
    }
    for nom_id in sorted(CURATED_SOMATIC_COUNTEREVIDENCE_NOMS):
        for g in CURATED_NOMINATION_TO_ENSCAFG.get(nom_id, ()):
            eid = gene_to_eid.get(g)
            if eid:
                gate_rows.append(
                    {
                        "entity_id": eid,
                        "layer": "Curated canine somatic RT-PCR counterevidence",
                        "unsafe": False,
                        "hold": True,
                        "reason": f"Published normal canine somatic RT-PCR detection ({nom_id})",
                    }
                )

    # Exclude essential activated canine T/NK/Th17/myeloid checkpoint & cytokine effectors
    for r_g in groups_df.itertuples(index=False):
        lbl_parts = {
            p.strip().upper() for p in str(r_g.protein_group_label).split("/") if p.strip()
        }
        hit_syms = sorted(lbl_parts & CURATED_CANINE_ACTIVATED_IMMUNE_SYMBOLS)
        if hit_syms:
            gate_rows.append(
                {
                    "entity_id": str(r_g.entity_id),
                    "layer": "Curated canine activated T/NK/myeloid effector safety",
                    "unsafe": True,
                    "hold": False,
                    "reason": f"Essential canine activated immune checkpoint/cytokine ({', '.join(hit_syms)})",
                }
            )

    return (
        pd.DataFrame(gate_rows, columns=["entity_id", "layer", "unsafe", "hold", "reason"]),
        repro_prot_status,
        source_rows,
    )


def load_dog_tumor_cohorts(
    dog_in: Path,
    extra_dir: Path,
    groups_df: pd.DataFrame,
    exonic_bp: dict[str, int],
    symbol_by_gene: dict[str, str],
    or_gids: set[str],
    risk_gids: set[str],
    unsafe_eids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    """Load 19 spontaneous canine patient tumor cohorts + 12 pure canine cancer cell lines on CanFam3.1."""
    exonic_series = pd.Series(exonic_bp, dtype=float)
    sym_candidates: dict[str, list[str]] = defaultdict(list)
    for gid, sym in symbol_by_gene.items():
        if sym and not sym.startswith("ENSCAFG"):
            sym_candidates[sym.upper()].append(gid)
    sym_to_ens = {k: v[0] for k, v in sym_candidates.items() if len(v) == 1}

    cohort_registry: list[dict] = []
    cell_line_registry: list[dict] = []
    source_rows: list[dict] = []
    # Store compactly: (source_id, cancer_code, sample_ids, participant_ids, peg_matrix)
    cohort_matrices: list[tuple[str, str, list[str], list[str], pd.DataFrame]] = []
    cell_line_matrices: list[tuple[str, str, str, list[str], pd.DataFrame]] = []

    def _qc_cols(tpm_frame: pd.DataFrame, cols: list[str]) -> list[str]:
        sh = _or_share_ge5(tpm_frame[cols], or_gids)
        return [c for c in cols if float(sh.get(c, 0.0)) <= 0.025]

    # 1. OSA_NCI (GSE238110, primary & metastatic canine osteosarcoma, SYM_3 raw counts -> TPM)
    p238 = dog_in / "GSE238110_RawCountFile_combined.csv.gz"
    source_rows.append(
        {
            "source_id": "gse238110-osa",
            "pillar": "tumor_rna",
            "filename": p238.name,
            "unit": "TPM (from CanFam3 raw counts)",
            "strain": "spontaneous companion dogs (NCI COTC021/022)",
            "bytes": p238.stat().st_size,
            "sha256": _sha256(p238),
        }
    )
    df238 = pd.read_csv(p238)
    sym238 = df238["gene_id"].astype(str).str.replace(r"_3$", "", regex=True).str.upper()
    gid238 = sym238.map(sym_to_ens)
    df238 = df238.loc[gid238.notna()].copy()
    df238["_gid"] = gid238.loc[gid238.notna()]
    df238 = df238.drop_duplicates("_gid").set_index("_gid")
    osa238_raw = [
        c
        for c in df238.columns
        if c != "gene_id" and re.search(r"(?:_tumor$|__Tumor$|_Met_)", str(c), re.IGNORECASE)
    ]
    tpm238 = _counts_to_tpm(df238[osa238_raw], exonic_series)
    osa238_cols = _qc_cols(tpm238, osa238_raw)
    peg238 = _project_peg_matrix(tpm238[osa238_cols], groups_df)
    osa238_pids = [
        "GSE238110:"
        + re.sub(r"(?:__|_)(?:tumor|Met.*|Primary.*)$", "", str(c), flags=re.IGNORECASE)
        for c in osa238_cols
    ]
    cohort_matrices.append(
        (
            "gse238110",
            "OSA_NCI",
            osa238_cols,
            osa238_pids,
            peg238,
        )
    )
    cohort_registry.append(
        {
            "cancer_code": "OSA_NCI",
            "label": "Spontaneous primary & metastatic appendicular osteosarcoma (NCI COTC021/022)",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "multi-breed companion dogs",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE238110",
            "n_samples": len(osa238_cols),
            "native_feature_space": "CanFam3 symbol -> ENSCAFG",
            "quantification": "linear TPM from raw integer counts and CanFam3.1 exonic_bp",
        }
    )

    # 2. OSA_COTC030 (GSE304066, 26 primary amputation canine osteosarcomas, verified raw integer counts -> TPM)
    p304 = dog_in / "GSE304066_all_sample_raw_count.txt.gz"
    source_rows.append(
        {
            "source_id": "gse304066-osa-cotc030",
            "pillar": "tumor_rna",
            "filename": p304.name,
            "unit": "TPM (from verified raw integer counts)",
            "strain": "spontaneous companion dogs (NCI COTC030)",
            "bytes": p304.stat().st_size,
            "sha256": _sha256(p304),
        }
    )
    df304 = pd.read_csv(p304, sep="\t")
    df304["_sym"] = df304["gene_id"].astype(str).str.replace(r"_4$", "", regex=True).str.upper()
    df304["_gid"] = df304["_sym"].map(sym_to_ens)
    df304 = df304.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    osa304_raw = [c for c in df304.columns if c.startswith("COTC030_")]
    tpm304 = _counts_to_tpm(df304[osa304_raw], exonic_series)
    osa304_cols = _qc_cols(tpm304, osa304_raw)
    peg304 = _project_peg_matrix(tpm304[osa304_cols], groups_df)
    cohort_matrices.append(
        (
            "gse304066",
            "OSA_COTC030",
            osa304_cols,
            [f"GSE304066:{c}" for c in osa304_cols],
            peg304,
        )
    )
    cohort_registry.append(
        {
            "cancer_code": "OSA_COTC030",
            "label": "Spontaneous primary appendicular osteosarcoma (NCI COTC030 trial)",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "multi-breed companion dogs",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE304066",
            "n_samples": len(osa304_cols),
            "native_feature_space": "Canine symbol -> ENSCAFG",
            "quantification": "linear TPM from verified raw integer counts and CanFam3.1 exonic_bp",
        }
    )

    # 3. OSA_CCOGC (GSE239948, 43 CCOGC multi-breed osteosarcomas, RSEM expected counts -> TPM)
    p239 = dog_in / "GSE239948_CCOGC.txt.gz"
    source_rows.append(
        {
            "source_id": "gse239948-osa-ccogc",
            "pillar": "tumor_rna",
            "filename": p239.name,
            "unit": "TPM (from CanFam3.1 expected counts)",
            "strain": "multi-breed companion dogs (CCOGC)",
            "bytes": p239.stat().st_size,
            "sha256": _sha256(p239),
        }
    )
    df239 = pd.read_csv(p239, sep="\t")
    df239["_gid"] = df239["Name"].astype(str).str.upper().map(sym_to_ens)
    df239 = df239.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
    osa239_raw = [c for c in df239.columns if c != "Name"]
    tpm239 = _counts_to_tpm(df239[osa239_raw], exonic_series)
    osa239_cols = _qc_cols(tpm239, osa239_raw)
    peg239 = _project_peg_matrix(tpm239[osa239_cols], groups_df)
    cohort_matrices.append(
        ("gse239948", "OSA_CCOGC", osa239_cols, [f"GSE239948:{c}" for c in osa239_cols], peg239)
    )
    cohort_registry.append(
        {
            "cancer_code": "OSA_CCOGC",
            "label": "Canine Comparative Oncology & Genomics Consortium (CCOGC) osteosarcoma",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "Labrador, Greyhound, Rottweiler, Golden Retriever & mixed",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE239948",
            "n_samples": len(osa239_cols),
            "native_feature_space": "CanFam3.1 symbol -> ENSCAFG",
            "quantification": "linear TPM from expected counts and CanFam3.1 exonic_bp",
        }
    )

    # 4. CMT (GSE119810, 158 spontaneous canine mammary carcinomas on native ENSCAFG, FPKM -> TPM)
    p119 = dog_in / "GSE119810_CMT_222S_FPKM.csv.gz"
    df119 = pd.read_csv(p119, index_col=0)
    df119.index = df119.index.astype(str).str.split(".").str[0]
    cmt_raw = [c for c in df119.columns if str(c).endswith("-tumor")]
    tpm119 = _fpkm_to_tpm(df119[cmt_raw])
    cmt_cols = _qc_cols(tpm119, cmt_raw)
    tpm119_t = tpm119[cmt_cols].copy()
    tpm119_t.loc[tpm119_t.index.isin(risk_gids)] = np.nan
    peg119 = _project_peg_matrix(tpm119_t, groups_df)
    cohort_matrices.append(
        (
            "gse119810",
            "CMT",
            cmt_cols,
            [f"GSE119810:{c.removesuffix('-tumor')}" for c in cmt_cols],
            peg119,
        )
    )
    cohort_registry.append(
        {
            "cancer_code": "CMT",
            "label": "Spontaneous canine mammary gland carcinoma",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "multi-breed female companion dogs",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE119810",
            "n_samples": len(cmt_cols),
            "native_feature_space": "Ensembl gene (ENSCAFG)",
            "quantification": "within-sample linear TPM from FPKM",
        }
    )

    # 5. HSA (GSE95183, spontaneous canine visceral hemangiosarcoma DHSA-* on native ENSCAFG, FPKM -> TPM)
    p953 = dog_in / "GSE95183_fpkms_hemangiosarcoma.txt.gz"
    source_rows.append(
        {
            "source_id": "gse95183-hsa",
            "pillar": "tumor_rna",
            "filename": p953.name,
            "unit": "TPM (from native ENSCAFG FPKM)",
            "strain": "Golden Retriever & multi-breed companion dogs",
            "bytes": p953.stat().st_size,
            "sha256": _sha256(p953),
        }
    )
    df953 = pd.read_csv(p953, sep="\t")
    df953["tracking_id"] = df953["tracking_id"].astype(str).str.split(".").str[0]
    hsa_raw = [c for c in df953.columns if c.startswith("DHSA-")]
    df953 = df953.groupby("tracking_id", sort=False)[hsa_raw].sum()
    tpm953 = _fpkm_to_tpm(df953)
    hsa_cols = _qc_cols(tpm953, hsa_raw)
    tpm953_t = tpm953[hsa_cols].copy()
    tpm953_t.loc[tpm953_t.index.isin(risk_gids)] = np.nan
    peg953 = _project_peg_matrix(tpm953_t, groups_df)
    cohort_matrices.append(
        ("gse95183", "HSA", hsa_cols, [f"GSE95183:{c.split('.')[0]}" for c in hsa_cols], peg953)
    )
    cohort_registry.append(
        {
            "cancer_code": "HSA",
            "label": "Spontaneous canine visceral hemangiosarcoma",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "multi-breed companion dogs",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE95183",
            "n_samples": len(hsa_cols),
            "native_feature_space": "Ensembl gene (ENSCAFG)",
            "quantification": "within-sample linear TPM from FPKM",
        }
    )

    # 6. LYM (GSE95184, spontaneous canine multicentric lymphoma DLYM-* on native ENSCAFG, FPKM -> TPM)
    p954 = dog_in / "GSE95184_fpkms_lymphoma.txt.gz"
    source_rows.append(
        {
            "source_id": "gse95184-lym",
            "pillar": "tumor_rna",
            "filename": p954.name,
            "unit": "TPM (from native ENSCAFG FPKM)",
            "strain": "multi-breed companion dogs",
            "bytes": p954.stat().st_size,
            "sha256": _sha256(p954),
        }
    )
    df954 = pd.read_csv(p954, sep="\t")
    df954["tracking_id"] = df954["tracking_id"].astype(str).str.split(".").str[0]
    lym_raw = [c for c in df954.columns if c.startswith("DLYM-")]
    df954 = df954.groupby("tracking_id", sort=False)[lym_raw].sum()
    tpm954 = _fpkm_to_tpm(df954)
    lym_cols = _qc_cols(tpm954, lym_raw)
    tpm954_t = tpm954[lym_cols].copy()
    tpm954_t.loc[tpm954_t.index.isin(risk_gids)] = np.nan
    peg954 = _project_peg_matrix(tpm954_t, groups_df)
    cohort_matrices.append(
        ("gse95184", "LYM", lym_cols, [f"GSE95184:{c.split('.')[0]}" for c in lym_cols], peg954)
    )
    cohort_registry.append(
        {
            "cancer_code": "LYM",
            "label": "Spontaneous canine multicentric B- and T-cell lymphoma",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "multi-breed companion dogs",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE95184",
            "n_samples": len(lym_cols),
            "native_feature_space": "Ensembl gene (ENSCAFG)",
            "quantification": "within-sample linear TPM from FPKM",
        }
    )

    # 7. CLBL-1 DLBCL cell line (GSE287971, moved to pure in-vitro cell-line track)
    p287 = dog_in / "GSE287971_CountTable.tsv.gz"
    source_rows.append(
        {
            "source_id": "gse287971-dlbcl",
            "pillar": "tumor_rna",
            "filename": p287.name,
            "unit": "TPM (from native ENSCAFG raw counts)",
            "strain": "German Shepherd (CLBL-1 DLBCL cell line)",
            "bytes": p287.stat().st_size,
            "sha256": _sha256(p287),
        }
    )
    raw287 = pd.read_csv(p287, sep="\t", dtype=str)
    sub_hdr = raw287.iloc[0]
    dmso_raw = [c for c in raw287.columns[1:] if str(sub_hdr[c]).endswith("-D")]
    cnt287 = raw287.iloc[1:].copy()
    cnt287["#ID"] = cnt287["#ID"].astype(str).str.split(".").str[0]
    cnt287 = cnt287[cnt287["#ID"].str.startswith("ENSCAFG")].drop_duplicates("#ID").set_index("#ID")
    tpm287 = _counts_to_tpm(cnt287[dmso_raw].astype(float), exonic_series)
    dmso_cols = _qc_cols(tpm287, dmso_raw)
    peg287 = _project_peg_matrix(tpm287[dmso_cols], groups_df)
    cell_line_matrices.append(("gse287971", "CLBL-1", "DLBCL", dmso_cols, peg287))
    cell_line_registry.append(
        {
            "model_id": "CLBL-1",
            "cancer_code": "DLBCL",
            "label": "Canine diffuse large B-cell lymphoma line CLBL-1 (vehicle controls)",
            "strain": "German Shepherd",
            "mhc_haplotype": "DLA-88",
            "source_accession": "GSE287971",
            "n_samples": len(dmso_cols),
            "native_feature_space": "Ensembl gene (ENSCAFG)",
            "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
        }
    )

    # 8. GLIOMA (GSE284559, 18 spontaneous brachycephalic canine high-grade oligodendrogliomas, ENSCAFG counts -> TPM)
    p284 = dog_in / "GSE284559_RawCountsArnold.txt.gz"
    source_rows.append(
        {
            "source_id": "gse284559-glioma",
            "pillar": "tumor_rna",
            "filename": p284.name,
            "unit": "TPM (from native ENSCAFG raw counts)",
            "strain": "French Bulldog, Boxer & Boston Terrier",
            "bytes": p284.stat().st_size,
            "sha256": _sha256(p284),
        }
    )
    df284 = pd.read_csv(p284, sep="\t")
    df284["_gid"] = df284["NewName"].astype(str).str.split("_", n=1).str[0].str.split(".").str[0]
    gli_raw = [c for c in df284.columns if c not in {"NewName", "_gid"}]
    df284 = df284[df284["_gid"].str.startswith("ENSCAFG")].drop_duplicates("_gid").set_index("_gid")
    tpm284 = _counts_to_tpm(df284[gli_raw].astype(float), exonic_series)
    gli_cols = _qc_cols(tpm284, gli_raw)
    peg284 = _project_peg_matrix(tpm284[gli_cols], groups_df)
    cohort_matrices.append(
        ("gse284559", "GLIOMA", gli_cols, [f"GSE284559:{c}" for c in gli_cols], peg284)
    )
    cohort_registry.append(
        {
            "cancer_code": "GLIOMA",
            "label": "Spontaneous prosencephalic high-grade oligodendroglioma",
            "model_kind": "spontaneous_patient_tumor",
            "strain": "French Bulldog, Boxer & Boston Terrier",
            "mhc_haplotype": "diverse DLA",
            "source_accession": "GSE284559",
            "n_samples": len(gli_cols),
            "native_feature_space": "Ensembl gene (ENSCAFG)",
            "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
        }
    )

    if extra_dir.is_dir():
        # 9. LYM_CHOP (GSE130874, 25 de novo multicentric canine B/T-cell lymphomas, ENSCAFG raw counts)
        p130 = extra_dir / "GSE130874_read_count_matrix.tsv.gz"
        if p130.is_file():
            source_rows.append(
                {
                    "source_id": "gse130874-lym-chop",
                    "pillar": "tumor_rna",
                    "filename": p130.name,
                    "unit": "TPM (from native ENSCAFG raw counts)",
                    "strain": "multi-breed companion dogs (de novo multicentric lymphoma)",
                    "bytes": p130.stat().st_size,
                    "sha256": _sha256(p130),
                }
            )
            df130 = pd.read_csv(p130, sep="\t", index_col=0)
            df130.index = df130.index.astype(str).str.split(".").str[0]
            df130 = df130[~df130.index.duplicated(keep="first")]
            lym130_raw = list(df130.columns)
            tpm130 = _counts_to_tpm(df130[lym130_raw], exonic_series)
            lym130_cols = _qc_cols(tpm130, lym130_raw)
            peg130 = _project_peg_matrix(tpm130[lym130_cols], groups_df)
            cohort_matrices.append(
                (
                    "gse130874",
                    "LYM_CHOP",
                    lym130_cols,
                    [f"GSE130874:{c.split('_')[0]}" for c in lym130_cols],
                    peg130,
                )
            )
            cohort_registry.append(
                {
                    "cancer_code": "LYM_CHOP",
                    "label": "Spontaneous de novo multicentric B- and T-cell lymphoma (CHOP cohort)",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE130874",
                    "n_samples": len(lym130_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
                }
            )

        # 10. GLIOMA_G23 (GSE160124, spontaneous Grade II/III canine oligodendroglioma; deduplicate FF vs FFPE per dog)
        p160 = extra_dir / "GSE160124_FCresult.txt.gz"
        if p160.is_file():
            source_rows.append(
                {
                    "source_id": "gse160124-glioma-g23",
                    "pillar": "tumor_rna",
                    "filename": p160.name,
                    "unit": "TPM (from native ENSCAFG featureCounts)",
                    "strain": "Boxer, French Bulldog & multi-breed companion dogs",
                    "bytes": p160.stat().st_size,
                    "sha256": _sha256(p160),
                }
            )
            df160 = pd.read_csv(p160, sep="\t", comment="#")
            df160["Geneid"] = df160["Geneid"].astype(str).str.split(".").str[0]
            df160 = df160.drop_duplicates("Geneid").set_index("Geneid")
            od_all = [c for c in df160.columns[5:] if str(c).endswith("_od.bam")]
            ff_dogs = {c.split("_")[0] for c in od_all if "_FF_od.bam" in c}
            od_dedup = [
                c for c in od_all if ("_FF_od.bam" in c) or (c.split("_")[0] not in ff_dogs)
            ]
            tpm160 = _counts_to_tpm(df160[od_dedup], df160["Length"].astype(float))
            od_cols = _qc_cols(tpm160, od_dedup)
            peg160 = _project_peg_matrix(tpm160[od_cols], groups_df)
            cohort_matrices.append(
                (
                    "gse160124",
                    "GLIOMA_G23",
                    od_cols,
                    [f"GSE160124:{c.split('_')[0]}" for c in od_cols],
                    peg160,
                )
            )
            cohort_registry.append(
                {
                    "cancer_code": "GLIOMA_G23",
                    "label": "Spontaneous Grade II/III oligodendroglioma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "Boxer, French Bulldog & companion breeds",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE160124",
                    "n_samples": len(od_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from featureCounts and Length",
                }
            )

        # 11. URO (GSE199292, spontaneous invasive bladder urothelial carcinoma, ENSCAFG expected counts)
        p199 = extra_dir / "GSE199292_All_Samples_Raw_Counts.txt.gz"
        if p199.is_file():
            source_rows.append(
                {
                    "source_id": "gse199292-uro",
                    "pillar": "tumor_rna",
                    "filename": p199.name,
                    "unit": "TPM (from native ENSCAFG counts)",
                    "strain": "spontaneous companion dogs (invasive bladder urothelial carcinoma)",
                    "bytes": p199.stat().st_size,
                    "sha256": _sha256(p199),
                }
            )
            df199 = pd.read_csv(p199, sep="\t")
            df199["Gene_ID"] = df199["Gene_ID"].astype(str).str.split(".").str[0]
            df199 = df199.drop_duplicates("Gene_ID").set_index("Gene_ID")
            uro_raw = list(df199.columns)
            tpm199 = _counts_to_tpm(df199[uro_raw], exonic_series)
            uro_cols = _qc_cols(tpm199, uro_raw)
            peg199 = _project_peg_matrix(tpm199[uro_cols], groups_df)
            cohort_matrices.append(
                ("gse199292", "URO", uro_cols, [f"GSE199292:{c}" for c in uro_cols], peg199)
            )
            cohort_registry.append(
                {
                    "cancer_code": "URO",
                    "label": "Spontaneous invasive urinary bladder urothelial carcinoma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE199292",
                    "n_samples": len(uro_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from raw/expected counts and CanFam3.1 exonic_bp",
                }
            )

        # 12. STS (GSE214154, 29 spontaneous subcutaneous/oral soft-tissue sarcomas, ENSCAFG raw counts)
        p214 = extra_dir / "GSE214154_Canine_SoftTissueSarcoma_RNAseq_Raw_CountData.csv.gz"
        if p214.is_file():
            source_rows.append(
                {
                    "source_id": "gse214154-sts",
                    "pillar": "tumor_rna",
                    "filename": p214.name,
                    "unit": "TPM (from native ENSCAFG raw counts)",
                    "strain": "multi-breed companion dogs",
                    "bytes": p214.stat().st_size,
                    "sha256": _sha256(p214),
                }
            )
            df214 = pd.read_csv(p214)
            df214["ensembl_gene_id"] = df214["ensembl_gene_id"].astype(str).str.split(".").str[0]
            df214 = df214.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
            sts_raw = [c for c in df214.columns if c != "Hugo_Symbol"]
            tpm214 = _counts_to_tpm(df214[sts_raw], exonic_series)
            sts_cols = _qc_cols(tpm214, sts_raw)
            peg214 = _project_peg_matrix(tpm214[sts_cols], groups_df)
            cohort_matrices.append(
                ("gse214154", "STS", sts_cols, [f"GSE214154:{c}" for c in sts_cols], peg214)
            )
            cohort_registry.append(
                {
                    "cancer_code": "STS",
                    "label": "Spontaneous subcutaneous & oral soft-tissue sarcoma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE214154",
                    "n_samples": len(sts_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
                }
            )

        # 13. STS_SUB (GSE208664, spontaneous fibrosarcoma, PNST & perivascular wall tumor; deduplicate T1..T3 per dog)
        p208 = extra_dir / "GSE208664_cSTS_HTSeq_counts.csv.gz"
        if p208.is_file():
            source_rows.append(
                {
                    "source_id": "gse208664-sts-sub",
                    "pillar": "tumor_rna",
                    "filename": p208.name,
                    "unit": "TPM (from HTSeq raw counts)",
                    "strain": "multi-breed companion dogs (fibrosarcoma, PNST, PWT)",
                    "bytes": p208.stat().st_size,
                    "sha256": _sha256(p208),
                }
            )
            df208 = pd.read_csv(p208, index_col=0)
            df208["_gid"] = (
                df208.index.astype(str).str.removeprefix("gene-").str.upper().map(sym_to_ens)
            )
            df208 = df208.dropna(subset=["_gid"]).drop_duplicates("_gid").set_index("_gid")
            sts208_raw = list(df208.columns)
            tpm208 = _counts_to_tpm(df208[sts208_raw], exonic_series)
            sts208_cols = _qc_cols(tpm208, sts208_raw)
            peg208 = _project_peg_matrix(tpm208[sts208_cols], groups_df)
            sts208_pids = ["GSE208664:" + re.sub(r"-T\d+$", "", str(c)) for c in sts208_cols]
            cohort_matrices.append(("gse208664", "STS_SUB", sts208_cols, sts208_pids, peg208))
            cohort_registry.append(
                {
                    "cancer_code": "STS_SUB",
                    "label": "Spontaneous fibrosarcoma, peripheral nerve sheath & perivascular wall sarcoma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs (16 patients)",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE208664",
                    "n_samples": len(sts208_cols),
                    "native_feature_space": "Canine gene-<SYMBOL> -> ENSCAFG",
                    "quantification": "linear TPM from HTSeq raw counts and CanFam3.1 exonic_bp",
                }
            )

        # 14-15. COSCC & CAA (GSE317244, spontaneous oral squamous cell carcinoma & acanthomatous ameloblastoma)
        p317 = extra_dir / "GSE317244_DS1_rawCounts.txt.gz"
        if p317.is_file():
            source_rows.append(
                {
                    "source_id": "gse317244-oral-tumors",
                    "pillar": "tumor_rna",
                    "filename": p317.name,
                    "unit": "TPM (from native ENSCAFG raw counts)",
                    "strain": "multi-breed companion dogs (oral tumors & gingiva)",
                    "bytes": p317.stat().st_size,
                    "sha256": _sha256(p317),
                }
            )
            df317 = pd.read_csv(p317, sep="\t", index_col=0)
            df317.index = df317.index.astype(str).str.split(".").str[0]
            df317 = df317[~df317.index.duplicated(keep="first")]
            for o_code, o_pfx, o_lbl in (
                ("COSCC", "COSCC-", "Spontaneous canine oral squamous cell carcinoma"),
                ("CAA", "CAA-", "Spontaneous canine oral acanthomatous ameloblastoma"),
            ):
                o_raw = [c for c in df317.columns if str(c).startswith(o_pfx)]
                tpm317 = _counts_to_tpm(df317[o_raw], exonic_series)
                o_cols = _qc_cols(tpm317, o_raw)
                peg317 = _project_peg_matrix(tpm317[o_cols], groups_df)
                cohort_matrices.append(
                    ("gse317244", o_code, o_cols, [f"GSE317244:{c}" for c in o_cols], peg317)
                )
                cohort_registry.append(
                    {
                        "cancer_code": o_code,
                        "label": o_lbl,
                        "model_kind": "spontaneous_patient_tumor",
                        "strain": "multi-breed companion dogs",
                        "mhc_haplotype": "diverse DLA",
                        "source_accession": "GSE317244",
                        "n_samples": len(o_cols),
                        "native_feature_space": "Ensembl gene (ENSCAFG)",
                        "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
                    }
                )

        # 16. HS (GSE288068, 18 spontaneous histiocytic sarcoma primary tumors + 3 HS cell lines BD/OD/PJ)
        p288 = extra_dir / "GSE288068_HS_MSU_RNAseq_Analysis_rawCounts.txt.gz"
        if p288.is_file():
            source_rows.append(
                {
                    "source_id": "gse288068-hs-tumors",
                    "pillar": "tumor_rna",
                    "filename": p288.name,
                    "unit": "TPM (from native ENSCAFG HTSeq raw counts)",
                    "strain": "Bernese Mountain Dog (spontaneous histiocytic sarcoma)",
                    "bytes": p288.stat().st_size,
                    "sha256": _sha256(p288),
                }
            )
            df288 = pd.read_csv(p288, sep="\t")
            df288["Symbol"] = df288["Symbol"].astype(str).str.split(".").str[0]
            df288 = df288.drop_duplicates("Symbol").set_index("Symbol")
            norm_288_set = {
                "HT5GWDSXX_18927780_S54",
                "HT5GWDSXX_18927782_S66",
                "HT5GWDSXX_18927808_S53",
                "HVTWLDSXX_18937521_S79",
                "HVTWLDSXX_18937523_S80",
                "HVTWLDSXX_18936778_S54",
            }
            hs_raw = [c for c in df288.columns if c not in norm_288_set]
            tpm288 = _counts_to_tpm(df288[hs_raw], exonic_series)
            hs_cols = _qc_cols(tpm288, hs_raw)
            peg288 = _project_peg_matrix(tpm288[hs_cols], groups_df)
            cohort_matrices.append(
                ("gse288068", "HS", hs_cols, [f"GSE288068:{c}" for c in hs_cols], peg288)
            )
            cohort_registry.append(
                {
                    "cancer_code": "HS",
                    "label": "Spontaneous disseminated & hemophagocytic histiocytic sarcoma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "Bernese Mountain Dog",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE288068",
                    "n_samples": len(hs_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from HTSeq raw counts and CanFam3.1 exonic_bp",
                }
            )

        p288_cl = extra_dir / "GSE288068_HS_BD_OD_PJ_RNAseq_rawCounts.txt.gz"
        if p288_cl.is_file():
            source_rows.append(
                {
                    "source_id": "gse288068-hs-lines",
                    "pillar": "tumor_rna",
                    "filename": p288_cl.name,
                    "unit": "TPM (from native ENSCAFG HTSeq raw counts)",
                    "strain": "Bernese Mountain Dog & Rottweiler HS cell lines (BD, OD, PJ)",
                    "bytes": p288_cl.stat().st_size,
                    "sha256": _sha256(p288_cl),
                }
            )
            df288_cl = pd.read_csv(p288_cl, sep="\t")
            df288_cl["Symbol"] = df288_cl["Symbol"].astype(str).str.split(".").str[0]
            df288_cl = df288_cl.drop_duplicates("Symbol").set_index("Symbol")
            hs_cl_raw = ["BD cells", "OD cells", "PJ cells"]
            tpm288_cl = _counts_to_tpm(df288_cl[hs_cl_raw], exonic_series)
            hs_cl_qc = _qc_cols(tpm288_cl, hs_cl_raw)
            peg288_cl = _project_peg_matrix(tpm288_cl[hs_cl_qc], groups_df)
            for c_cl, m_id, breed_cl in (
                ("BD cells", "HS_BD", "Bernese Mountain Dog"),
                ("OD cells", "HS_OD", "Bernese Mountain Dog"),
                ("PJ cells", "HS_PJ", "Rottweiler"),
            ):
                if c_cl in hs_cl_qc:
                    cell_line_matrices.append(("gse288068", m_id, "HS", [c_cl], peg288_cl[[c_cl]]))
                    cell_line_registry.append(
                        {
                            "model_id": m_id,
                            "cancer_code": "HS",
                            "label": f"Canine histiocytic sarcoma cell line {m_id}",
                            "strain": breed_cl,
                            "mhc_haplotype": "DLA",
                            "source_accession": "GSE288068",
                            "n_samples": 1,
                            "native_feature_space": "Ensembl gene (ENSCAFG)",
                            "quantification": "linear TPM from HTSeq raw counts and CanFam3.1 exonic_bp",
                        }
                    )

        # 17. MNG (GSE250310, 62 spontaneous canine meningiomas, ENSCAFG featureCounts + Length)
        p250 = extra_dir / "GSE250310_counts.txt.gz"
        if p250.is_file():
            source_rows.append(
                {
                    "source_id": "gse250310-meningioma",
                    "pillar": "tumor_rna",
                    "filename": p250.name,
                    "unit": "TPM (from native ENSCAFG featureCounts)",
                    "strain": "Golden Retriever, Boxer & multi-breed companion dogs",
                    "bytes": p250.stat().st_size,
                    "sha256": _sha256(p250),
                }
            )
            df250 = pd.read_csv(p250, sep="\t", comment="#")
            df250["Geneid"] = df250["Geneid"].astype(str).str.split(".").str[0]
            df250 = df250.drop_duplicates("Geneid").set_index("Geneid")
            mng_raw = list(df250.columns[5:])
            tpm250 = _counts_to_tpm(df250[mng_raw], df250["Length"].astype(float))
            mng_cols = _qc_cols(tpm250, mng_raw)
            peg250 = _project_peg_matrix(tpm250[mng_cols], groups_df)
            mng_sids = [c.split("_R1_")[0] for c in mng_cols]
            peg250.columns = mng_sids
            cohort_matrices.append(
                ("gse250310", "MNG", mng_sids, [f"GSE250310:{s}" for s in mng_sids], peg250)
            )
            cohort_registry.append(
                {
                    "cancer_code": "MNG",
                    "label": "Spontaneous intracranial & spinal meningioma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE250310",
                    "n_samples": len(mng_sids),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from featureCounts and Length",
                }
            )

        # 18. ACC (GSE196108, 49 spontaneous cortisol-secreting adrenocortical tumors ACT1..49)
        p196 = extra_dir / "GSE196108_Sanders_raw_counts.txt.gz"
        if p196.is_file():
            source_rows.append(
                {
                    "source_id": "gse196108-acc",
                    "pillar": "tumor_rna",
                    "filename": p196.name,
                    "unit": "TPM (from native ENSCAFG raw counts)",
                    "strain": "spontaneous companion dogs (cortisol-secreting adrenocortical tumors)",
                    "bytes": p196.stat().st_size,
                    "sha256": _sha256(p196),
                }
            )
            df196 = pd.read_csv(p196, sep="\t")
            df196["gene"] = df196["gene"].astype(str).str.split(".").str[0]
            df196 = df196.drop_duplicates("gene").set_index("gene")
            acc_raw = [c for c in df196.columns if c.startswith("ACT")]
            tpm196 = _counts_to_tpm(df196[acc_raw], exonic_series)
            acc_cols = _qc_cols(tpm196, acc_raw)
            peg196 = _project_peg_matrix(tpm196[acc_cols], groups_df)
            cohort_matrices.append(
                ("gse196108", "ACC", acc_cols, [f"GSE196108:{c}" for c in acc_cols], peg196)
            )
            cohort_registry.append(
                {
                    "cancer_code": "ACC",
                    "label": "Spontaneous functional adrenocortical adenoma & carcinoma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE196108",
                    "n_samples": len(acc_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
                }
            )

        # 19. THCA (GSE289443, 30 spontaneous canine thyroid carcinomas TC*)
        p289 = extra_dir / "GSE289443_Canine_Thyroid_Carcinoma_Raw_CountData.csv.gz"
        if p289.is_file():
            source_rows.append(
                {
                    "source_id": "gse289443-thca",
                    "pillar": "tumor_rna",
                    "filename": p289.name,
                    "unit": "TPM (from native ENSCAFG raw counts)",
                    "strain": "multi-breed companion dogs (thyroid carcinoma)",
                    "bytes": p289.stat().st_size,
                    "sha256": _sha256(p289),
                }
            )
            df289 = pd.read_csv(p289)
            df289["ensembl_gene_id"] = df289["ensembl_gene_id"].astype(str).str.split(".").str[0]
            df289 = df289.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
            thca_raw = [c for c in df289.columns if c != "Hugo_Symbol"]
            tpm289 = _counts_to_tpm(df289[thca_raw], exonic_series)
            thca_cols = _qc_cols(tpm289, thca_raw)
            peg289 = _project_peg_matrix(tpm289[thca_cols], groups_df)
            cohort_matrices.append(
                ("gse289443", "THCA", thca_cols, [f"GSE289443:{c}" for c in thca_cols], peg289)
            )
            cohort_registry.append(
                {
                    "cancer_code": "THCA",
                    "label": "Spontaneous follicular & compact thyroid gland carcinoma",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE289443",
                    "n_samples": len(thca_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from raw counts and CanFam3.1 exonic_bp",
                }
            )

        # 20. LEUK (GSE345838, ~110 spontaneous canine acute leukemias AL*, ENSCAFG featureCounts + Length)
        p345 = extra_dir / "GSE345838_Canine_Acute_Leukemia_featurecounts_raw_counts.txt.gz"
        if p345.is_file():
            source_rows.append(
                {
                    "source_id": "gse345838-acute-leukemia",
                    "pillar": "tumor_rna",
                    "filename": p345.name,
                    "unit": "TPM (from native ENSCAFG featureCounts)",
                    "strain": "multi-breed companion dogs (acute leukemia + BL/BM controls)",
                    "bytes": p345.stat().st_size,
                    "sha256": _sha256(p345),
                }
            )
            df345 = pd.read_csv(p345, sep="\t", comment="#")
            df345["Geneid"] = df345["Geneid"].astype(str).str.split(".").str[0]
            df345 = df345.drop_duplicates("Geneid").set_index("Geneid")
            al_raw = [c for c in df345.columns[5:] if str(c).startswith("AL")]
            tpm345 = _counts_to_tpm(df345[al_raw], df345["Length"].astype(float))
            al_cols = _qc_cols(tpm345, al_raw)
            peg345 = _project_peg_matrix(tpm345[al_cols], groups_df)
            cohort_matrices.append(
                ("gse345838", "LEUK", al_cols, [f"GSE345838:{c}" for c in al_cols], peg345)
            )
            cohort_registry.append(
                {
                    "cancer_code": "LEUK",
                    "label": "Spontaneous canine acute myeloid & lymphoid leukemia",
                    "model_kind": "spontaneous_patient_tumor",
                    "strain": "multi-breed companion dogs",
                    "mhc_haplotype": "diverse DLA",
                    "source_accession": "GSE345838",
                    "n_samples": len(al_cols),
                    "native_feature_space": "Ensembl gene (ENSCAFG)",
                    "quantification": "linear TPM from featureCounts and Length",
                }
            )

        # Pure in-vitro mammary tumor cell lines (GSE230271, 8 MGT lines; exclude MDCK & TP53 K/O fibroblast)
        p230 = extra_dir / "GSE230271_MGT_10S_FPKM.csv.gz"
        if p230.is_file():
            source_rows.append(
                {
                    "source_id": "gse230271-mgt-lines",
                    "pillar": "tumor_rna",
                    "filename": p230.name,
                    "unit": "TPM (from native ENSCAFG FPKM)",
                    "strain": "canine mammary carcinoma cell lines (MGT 122..703#2)",
                    "bytes": p230.stat().st_size,
                    "sha256": _sha256(p230),
                }
            )
            df230 = pd.read_csv(p230, index_col=0)
            df230.index = df230.index.astype(str).str.split(".").str[0]
            mgt_raw = [c for c in df230.columns if str(c).startswith("MGT")]
            tpm230 = _fpkm_to_tpm(df230[mgt_raw])
            mgt_qc = _qc_cols(tpm230, mgt_raw)
            tpm230_t = tpm230[mgt_qc].copy()
            tpm230_t.loc[tpm230_t.index.isin(risk_gids)] = np.nan
            peg230 = _project_peg_matrix(tpm230_t, groups_df)
            for c_mgt in mgt_qc:
                m_id = c_mgt.replace(" ", "_").split("#")[0]
                cell_line_matrices.append(("gse230271", m_id, "CMT", [c_mgt], peg230[[c_mgt]]))
                cell_line_registry.append(
                    {
                        "model_id": m_id,
                        "cancer_code": "CMT",
                        "label": f"Canine mammary gland tumor cell line {c_mgt}",
                        "strain": "companion dog mammary line",
                        "mhc_haplotype": "DLA",
                        "source_accession": "GSE230271",
                        "n_samples": 1,
                        "native_feature_space": "Ensembl gene (ENSCAFG)",
                        "quantification": "within-sample linear TPM from FPKM",
                    }
                )

    # Emit exact comparison_tumors concisely for non-Unsafe PEGs: every sample observation >= 10.0 TPM
    # PLUS the global peak sample per PEG when max < 10.0 TPM (Unsafe PEGs exit in Stage 1 before tumor comparison).
    eids_list = groups_df["entity_id"].tolist()
    non_unsafe_mask = np.array([eid not in unsafe_eids for eid in eids_list], dtype=bool)
    best_below10: dict[str, tuple[float, str, str, str, str]] = {}
    tumor_rows: list[dict] = []

    for src_id, code, s_cols, p_ids, peg_mat in cohort_matrices:
        if not s_cols:
            continue
        arr = peg_mat[s_cols].to_numpy(dtype=np.float64)
        # 1. All non-Unsafe entries >= 10.0 TPM
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
        # 2. Track global peak per non-Unsafe PEG across all cohorts for rows whose max < 10.0
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


def build_all(scratch_dir: Path, out_root: Path, extra_dir: Path | None = None) -> dict:
    dog_in = scratch_dir / "dog_inputs"
    eff_extra = extra_dir if extra_dir is not None else DEFAULT_EXTRA_DIR
    print("1/6 Qualifying Ensembl 98 CanFam3.1 proteome & 31-mer cDNA mapping risk...", flush=True)
    (
        ref,
        groups_df,
        seq_df,
        exonic_bp,
        symbol_by_gene,
        or_gids,
        single_exon_gids,
        mapping_risk_df,
        sym_index_df,
        ens_sources,
    ) = load_ensembl98_dog_reference(dog_in)
    print(
        f"    PEGs={len(groups_df):,}, sequences={len(seq_df):,}, coding_genes={len(symbol_by_gene):,}, OR_genes={len(or_gids)}",
        flush=True,
    )

    risk_gids = set(
        mapping_risk_df.loc[mapping_risk_df["mapping_risk"].astype(bool), "gene_id"].astype(str)
    )
    print(
        "2/6 Aggregating Bgee 15.2, GSE219045 & companion-dog normal/reproductive TPM panels...",
        flush=True,
    )
    comp_normals, comp_repro, rna_gates, confirmed_organs, _male_skin_p95, norm_sources = (
        load_dog_normal_and_repro_panels(
            scratch_dir, dog_in, eff_extra, groups_df, exonic_bp, symbol_by_gene, or_gids, risk_gids
        )
    )
    print(
        f"    normal_obs={len(comp_normals):,}, repro_obs={len(comp_repro):,}, rna_gates={len(rna_gates):,}",
        flush=True,
    )

    print(
        "3/6 Evaluating Pfizer 38-tissue proteomics & GSE225599 healthy canine PBMC 10x scRNA-seq...",
        flush=True,
    )
    prot_gates, repro_prot_status, prot_sources = load_pfizer_dog_proteomics_gates(
        dog_in, groups_df, seq_df, confirmed_organs
    )
    umi_df, umi_sources = load_dog_healthy_pbmc_scrna_umi(eff_extra, set(symbol_by_gene))
    print(
        f"    prot_gates={len(prot_gates):,}, repro_prot_positives={len(repro_prot_status):,}, umi_rows={len(umi_df):,}",
        flush=True,
    )

    all_gates = pd.concat([rna_gates, prot_gates], ignore_index=True)
    unsafe_eids = set(all_gates.loc[all_gates["unsafe"].astype(bool), "entity_id"].astype(str))

    # Also mask single-exon intronless retrogenes sharing >=50% 9-mers with Unsafe proteins
    # on Cufflinks mismatch-tolerant FPKM cohorts (GSE119810, GSE95183, GSE95184, GSE230271)
    unsafe_9mers: set[str] = set()
    for useq in seq_df.loc[seq_df["entity_id"].isin(unsafe_eids), "sequence"].astype(str):
        if len(useq) >= 9:
            unsafe_9mers.update(useq[i : i + 9] for i in range(len(useq) - 8))
    single_exon_unsafe50_gids: set[str] = set()
    single_exon_unsafe50_eids: set[str] = set()
    seqs_by_eid: dict[str, list[str]] = defaultdict(list)
    for r_s in seq_df.itertuples(index=False):
        if r_s.sequence:
            seqs_by_eid[str(r_s.entity_id)].append(str(r_s.sequence))
    for r_g in groups_df.itertuples(index=False):
        eid_s = str(r_g.entity_id)
        if eid_s in unsafe_eids:
            continue
        m_gids = json.loads(r_g.member_gene_ids)
        if not all(g in single_exon_gids for g in m_gids):
            continue
        cand_9mers = {
            s[i : i + 9]
            for s in seqs_by_eid.get(eid_s, ())
            if len(s) >= 9
            for i in range(len(s) - 8)
        }
        if cand_9mers and (len(cand_9mers & unsafe_9mers) / len(cand_9mers)) >= 0.50:
            single_exon_unsafe50_gids.update(m_gids)
            single_exon_unsafe50_eids.add(eid_s)
    if single_exon_unsafe50_eids:
        retro_holds = pd.DataFrame(
            [
                {
                    "entity_id": eid_r,
                    "layer": "CanFam3.1 intronless Unsafe-paralog retrocopy review",
                    "unsafe": False,
                    "hold": True,
                    "reason": "Single-exon locus sharing >=50% 9-mers with Unsafe somatic paralog (processed retrocopy risk)",
                }
                for eid_r in sorted(single_exon_unsafe50_eids)
            ]
        )
        all_gates = pd.concat([all_gates, retro_holds], ignore_index=True)
        comp_repro = comp_repro[
            ~(
                comp_repro["entity_id"].isin(single_exon_unsafe50_eids)
                & comp_repro["source_id"].eq("gse119810-normal-mammary")
            )
        ].reset_index(drop=True)
    cufflinks_mask_gids = risk_gids | single_exon_unsafe50_gids

    print("4/6 Loading spontaneous canine patient tumor cohorts & cancer cell lines...", flush=True)
    comp_tumors, comp_cell_lines, cohort_df, cell_line_df, tumor_sources = load_dog_tumor_cohorts(
        dog_in,
        eff_extra,
        groups_df,
        exonic_bp,
        symbol_by_gene,
        or_gids,
        cufflinks_mask_gids,
        unsafe_eids,
    )
    print(
        f"    tumor_obs={len(comp_tumors):,}, cohorts={len(cohort_df)} ({int(cohort_df['n_samples'].sum())} patient samples), cell_lines={len(cell_line_df)}",
        flush=True,
    )

    print("5/6 Running SpeciesScreeningSpec & canine HostContext sex-restriction...", flush=True)
    spec = SpeciesScreeningSpec(
        species="dog",
        reference=ref,
        policy_id=POLICY_SCHEMA_ID,
        comparison_id=COMPARISON_ID,
        comparability_reference="oncoref.dog-canfam31-panel.2026-11-26.2",
        comparison_unit="TPM",
        required_normal_contexts=REQUIRED_DOG_NORMAL_CONTEXTS,
        anatomy=anatomy_for_species("dog"),
    )
    res = evaluate_species_catalog(
        spec,
        groups_df,
        sequences=seq_df,
        normal_gates=all_gates,
        rare_cell_counts=umi_df if not umi_df.empty else None,
        mapping_risk=mapping_risk_df,
        confirmed_rare_cell_organs=confirmed_organs,
        comparison_normals=comp_normals,
        comparison_reproductive=comp_repro,
        comparison_tumors=comp_tumors,
        comparison_cell_lines=comp_cell_lines if not comp_cell_lines.empty else None,
    )

    cat = res.catalog.copy()
    full_host_exposure: dict[str, dict[str, dict]] = {eid: {} for eid in cat["entity_id"]}
    for r in comp_repro.itertuples(index=False):
        eid = str(r.entity_id)
        comp = str(r.compartment)
        is_dev = r.scope == "developmental"
        h_val = float(getattr(r, "host_max_tpm", r.value))
        p_stat = repro_prot_status.get((eid, comp), "not_assessed")
        prev = full_host_exposure[eid].get(comp)
        if prev is None or h_val > prev["max_tpm"]:
            entry: dict[str, object] = {
                "coverage": "measured",
                "max_tpm": round(h_val, 6),
                "max_ntpm": None,
                "protein_status": p_stat,
                "scope_complete": not is_dev,
            }
            if comp == "resting_mammary":
                # Female mammary biopsies only; male canine mammary RNA is unmeasured (fail-closed None)
                entry["female_max_tpm"] = round(h_val, 6)
            full_host_exposure[eid][comp] = entry

    # Record honest proteome_only coverage on canine fallopian tube (assayed in Pfizer Pool3 126)
    for eid in cat["entity_id"]:
        eid_s = str(eid)
        if "fallopian tube" not in full_host_exposure[eid_s]:
            full_host_exposure[eid_s]["fallopian tube"] = {
                "coverage": "proteome_only",
                "max_tpm": None,
                "max_ntpm": None,
                "protein_status": repro_prot_status.get((eid_s, "fallopian tube"), "not_assessed"),
                "scope_complete": False,
            }

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

    sex_df = classify_sex_restriction(cat, species="dog", host_evidence=full_host_exposure)
    sex_class_map = (
        sex_df.set_index("entity_id")["sex_restriction_class"].to_dict() if not sex_df.empty else {}
    )
    sex_code_map = (
        sex_df.set_index("entity_id")["sex_restriction_code"].to_dict() if not sex_df.empty else {}
    )
    cat["sex_restriction_class"] = cat["entity_id"].map(sex_class_map).fillna("")
    cat["sex_restriction_code"] = cat["entity_id"].map(sex_code_map).fillna("")

    data_dir = out_root / "oncoref" / "data"
    dog_dir = data_dir / "species" / "dog"
    dog_dir.mkdir(parents=True, exist_ok=True)

    sym_index_df.to_csv(
        dog_dir / "symbol-index.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )

    cat_path = data_dir / "dog-target-selection.csv.gz"
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
    exp_path = data_dir / "dog-host-exposure.csv.gz"
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
    seq_path = data_dir / "dog-proteome-sequences.csv.gz"
    seq_out.to_csv(
        seq_path,
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )

    all_sources_df = pd.DataFrame(
        ens_sources + norm_sources + prot_sources + umi_sources + tumor_sources
    ).drop_duplicates("source_id")
    all_sources_df.to_csv(dog_dir / "screening-sources.csv", index=False, lineterminator="\n")
    cohort_df.to_csv(dog_dir / "tumor-cohorts.csv", index=False, lineterminator="\n")
    cell_line_df.to_csv(dog_dir / "cell-lines.csv", index=False, lineterminator="\n")

    cat_sha = _sha256(cat_path)
    exp_sha = _sha256(exp_path)
    seq_sha = _sha256(seq_path)
    if "representative_transcripts" in cat.columns:
        import gzip

        rep_doc = {
            "schema": "oncoref.species-representative-transcripts.v1",
            "species": "dog",
            "taxon_id": 9615,
            "catalog_data_sha256": cat_sha,
            "sequences_sha256": seq_sha,
            "rows": [json.loads(raw) for raw in cat["representative_transcripts"]],
        }
        (dog_dir / "representative-transcripts.json.gz").write_bytes(
            gzip.compress(
                json.dumps(rep_doc, separators=(",", ":")).encode("utf-8"),
                compresslevel=9,
                mtime=0,
            )
        )
    loci_count = int(cat["member_gene_ids"].map(json.loads).map(len).sum())

    summary_payload = {
        "schema": "oncoref.species-target-selection-metadata.v1",
        "taxon_id": 9615,
        "species": "dog",
        "scientific_name": "Canis lupus familiaris",
        "selection_status": "evaluated",
        "policy_id": POLICY_SCHEMA_ID,
        "comparison_id": COMPARISON_ID,
        "comparison_unit": "TPM",
        "symbol_index": "species/dog/symbol-index.csv.gz",
        "reference": ref.as_dict(),
        "data_sha256": cat_sha,
        "data_bytes": cat_path.stat().st_size,
        "host_exposure_sha256": exp_sha,
        "host_exposure_bytes": exp_path.stat().st_size,
        "sequences_sha256": seq_sha,
        "sequences_bytes": seq_path.stat().st_size,
        "required_normal_contexts": sorted(REQUIRED_DOG_NORMAL_CONTEXTS),
        "measured_reproductive_compartments": sorted(
            {r.compartment for r in comp_repro.itertuples(index=False)}
        ),
        "proteome_only_reproductive_compartments": ["fallopian tube"],
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
    (dog_dir / "target-selection.json").write_text(json.dumps(summary_payload, indent=2) + "\n")
    print("DONE! Summary:", json.dumps(summary_payload["summary"], indent=2), flush=True)
    return summary_payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH_DIR)
    parser.add_argument("--extra-dir", type=Path, default=DEFAULT_EXTRA_DIR)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    build_all(args.scratch_dir, args.root, args.extra_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
