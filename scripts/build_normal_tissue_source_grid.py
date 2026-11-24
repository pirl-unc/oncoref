#!/usr/bin/env python3
"""Build comprehensive human normal-tissue x source coverage grid & GTEx donor-sex table (#778)."""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import (  # noqa: E402
    gtex_deconvolution,
    normal_tissues,
    reference_data,
    unmeasured_home_panels,
)

# Map all 68 GTEx v11 median columns (52 primary shards + 14 histological subcompartments + 2 cell lines)
# onto canonical normal_tissues concept_ids.
GTEX_V11_MEDIAN_TO_CONCEPT: dict[str, str] = {
    "Adipose_Subcutaneous": "adipose_tissue",
    "Adipose_Visceral_Omentum": "adipose_tissue",
    "Adrenal_Gland": "adrenal_gland",
    "Artery_Aorta": "blood_vessel",
    "Artery_Coronary": "blood_vessel",
    "Artery_Tibial": "blood_vessel",
    "Bladder": "urinary_bladder",
    "Brain_Amygdala": "amygdala",
    "Brain_Anterior_cingulate_cortex_BA24": "cerebral_cortex",
    "Brain_Caudate_basal_ganglia": "caudate",
    "Brain_Cerebellar_Hemisphere": "cerebellum",
    "Brain_Cerebellum": "cerebellum",
    "Brain_Cortex": "cerebral_cortex",
    "Brain_Frontal_Cortex_BA9": "cerebral_cortex",
    "Brain_Hippocampus": "hippocampal_formation",
    "Brain_Hypothalamus": "hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia": "basal_ganglia",
    "Brain_Putamen_basal_ganglia": "basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1": "spinal_cord",
    "Brain_Substantia_nigra": "substantia nigra",
    "Breast_Mammary_Tissue": "mammary_tissue",
    "Cervix_Ectocervix": "cervix",
    "Cervix_Endocervix": "cervix",
    "Colon_Sigmoid": "colon",
    "Colon_Transverse": "colon",
    "Colon_Transverse_Mixed_Cell": "colon",
    "Colon_Transverse_Mucosa": "colon",
    "Colon_Transverse_Muscularis": "colon",
    "Esophagus_Gastroesophageal_Junction": "esophagus",
    "Esophagus_Mucosa": "esophagus",
    "Esophagus_Muscularis": "esophagus",
    "Fallopian_Tube": "fallopian_tube",
    "Heart_Atrial_Appendage": "heart muscle",
    "Heart_Left_Ventricle": "heart muscle",
    "Kidney_Cortex": "renal_cortex",
    "Kidney_Medulla": "kidney",
    "Liver": "liver",
    "Liver_Hepatocyte": "liver",
    "Liver_Mixed_Cell": "liver",
    "Liver_Portal_Tract": "liver",
    "Lung": "lung",
    "Minor_Salivary_Gland": "salivary_gland",
    "Muscle_Skeletal": "skeletal_muscle",
    "Nerve_Tibial": "peripheral_nerve",
    "Ovary": "ovary",
    "Pancreas": "pancreas",
    "Pancreas_Acini": "pancreas",
    "Pancreas_Islets": "pancreas",
    "Pancreas_Mixed_Cell": "pancreas",
    "Pituitary": "pituitary_gland",
    "Prostate": "prostate",
    "Skin_Not_Sun_Exposed_Suprapubic": "skin",
    "Skin_Sun_Exposed_Lower_leg": "skin",
    "Small_Intestine_Terminal_Ileum": "small_intestine",
    "Small_Intestine_Terminal_Ileum_Lymphoid_Aggregate": "small_intestine",
    "Small_Intestine_Terminal_Ileum_Mixed_Cell": "small_intestine",
    "Spleen": "spleen",
    "Stomach": "stomach",
    "Stomach_Mixed_Cell": "stomach",
    "Stomach_Mucosa": "stomach",
    "Stomach_Muscularis": "stomach",
    "Testis": "testis",
    "Thyroid": "thyroid_gland",
    "Uterus": "uterus",
    "Vagina": "vagina",
    "Whole_Blood": "blood",
}

SUBCOMPARTMENT_COLS = {
    "Colon_Transverse_Mixed_Cell",
    "Colon_Transverse_Mucosa",
    "Colon_Transverse_Muscularis",
    "Liver_Hepatocyte",
    "Liver_Mixed_Cell",
    "Liver_Portal_Tract",
    "Pancreas_Acini",
    "Pancreas_Islets",
    "Pancreas_Mixed_Cell",
    "Small_Intestine_Terminal_Ileum_Lymphoid_Aggregate",
    "Small_Intestine_Terminal_Ileum_Mixed_Cell",
    "Stomach_Mixed_Cell",
    "Stomach_Mucosa",
    "Stomach_Muscularis",
}


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_donor_sex() -> None:
    """Extract splash-robust chromosomal sex (XX vs XY) across all 946 GTEx v11 donors."""
    pins = sorted(reference_data.REFERENCE_SOURCES["gtex_individual_tissue_tpm"]["pins"])
    y_bytes = {g.encode("ascii"): g for g in gtex_deconvolution.Y_CHROMOSOME_SENTINEL_GENES}
    donor_shard_y: dict[str, list[float]] = {}

    for pk in pins:
        p = reference_data.pinned_path("gtex_individual_tissue_tpm", pk)
        with gzip.open(p, "rb") as fh:
            fh.readline()
            fh.readline()
            samps = fh.readline().decode("utf-8").rstrip("\n").split("\t")[2:]
            y_tot = np.zeros(len(samps), dtype=np.float64)
            for line in fh:
                rn, _, rest = line.partition(b"\t")
                if rn.split(b".", 1)[0] in y_bytes:
                    y_tot += np.asarray(rest.partition(b"\t")[2].split(), dtype=np.float64)
        for sid, val in zip(samps, y_tot):
            did = gtex_deconvolution.parse_gtex_sample_barcode(sid).donor_id
            donor_shard_y.setdefault(did, []).append(float(val))

    rows = []
    for did in sorted(donor_shard_y):
        vals = donor_shard_y[did]
        med_y = float(np.median(vals))
        sex = "XY" if med_y >= 20.0 else ("XX" if med_y <= 2.0 else "ambiguous")
        rows.append(
            {
                "donor_id": did,
                "chromosomal_sex": sex,
                "median_y_sentinel_tpm": round(med_y, 4),
                "sequenced_shards": len(vals),
            }
        )
    out_path = ROOT / "oncoref/data/gtex-donor-sex.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False, lineterminator="\n")

    receipt_path = ROOT / "oncoref/data/gtex-deconvolution-receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["files_sha256"]["gtex-donor-sex.csv"] = _sha256_file(out_path)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} donors to {out_path.name}")


def build_source_grid() -> None:
    """Build oncoref/data/normal-tissue-source-grid.csv across all human normal sources."""
    st_path = ROOT / "docs/audits/target-selection-20261110/rebuild/normal-rna/source-tissues.csv"
    st = pd.read_csv(st_path).fillna("")
    t_summary = gtex_deconvolution.gtex_deconvolution_tissue_summary().set_index("site")

    # Map concept label / tissue name to canonical concept_id
    tax_concepts = normal_tissues.normal_tissue_taxonomy()["concepts"]
    label_to_cid: dict[str, str] = {}
    cid_to_info: dict[str, dict] = {}

    def _reg_concept(label: str) -> str:
        ann = normal_tissues.resolve_normal_tissue(label, species="human")
        if ann.concept_id:
            cid = ann.concept_id
            if cid not in cid_to_info:
                det = normal_tissues.normal_tissue_screening_details(ann)
                cid_to_info[cid] = {
                    "concept_id": cid,
                    "tissue_label": tax_concepts[cid]["label"],
                    "tissue_category": "thymic" if cid == "thymus" else det["tissue_category"],
                }
            label_to_cid[label] = cid
            return cid
        return label

    for t in sorted(set(st["tissue"]) - {""}):
        if t == "cornea/conjunctiva":
            _reg_concept("cornea")
            label_to_cid[t] = "cornea"
        else:
            _reg_concept(t)

    for c_name in set(GTEX_V11_MEDIAN_TO_CONCEPT.values()):
        _reg_concept(c_name)
    for extra in (
        "olfactory epithelium",
        "nasal mucosa",
        "hair follicle",
        "taste bud",
        "dental pulp",
        "bone",
        "cartilage",
        "bronchus",
        "embryo",
    ):
        _reg_concept(extra)

    # Bucket collector per concept_id
    buckets: dict[str, dict[str, set[str]]] = {
        cid: {
            "gtex_v11_shards": set(),
            "gtex_v11_medians": set(),
            "gtex_v11_subcompartments": set(),
            "gtex_via_hpa_v8": set(),
            "hpa_surgical_biopsies": set(),
            "hpa_tissue_rna": set(),
            "hpa_rna_consensus": set(),
            "hpa_single_cell": set(),
            "cellxgene_census": set(),
            "sorted_immune": set(),
            "specialized_and_developmental": set(),
            "hpa_ihc_protein": set(),
        }
        for cid in cid_to_info
    }

    # 1. GTEx v11 medians & shards (hard-assert every column maps to a canonical concept_id)
    for col, raw_cname in GTEX_V11_MEDIAN_TO_CONCEPT.items():
        ann = normal_tissues.resolve_normal_tissue(raw_cname, species="human")
        if not ann.concept_id:
            raise AssertionError(
                f"GTEX_V11_MEDIAN_TO_CONCEPT[{col!r}] = {raw_cname!r} failed to resolve to a concept_id"
            )
        cid = _reg_concept(raw_cname)
        assert cid in buckets, f"Missing bucket for {cid!r} ({col!r})"
        buckets[cid]["gtex_v11_medians"].add(col)
        if col in SUBCOMPARTMENT_COLS:
            buckets[cid]["gtex_v11_subcompartments"].add(col)
        else:
            site = col.lower()
            if site in t_summary.index:
                n_s = int(t_summary.at[site, "total_samples"])
                buckets[cid]["gtex_v11_shards"].add(f"{site} ({n_s})")
                if site.startswith("brain_") and "brain" in buckets:
                    buckets["brain"]["gtex_v11_shards"].add(f"{site} ({n_s})")
                    buckets["brain"]["gtex_v11_medians"].add(col)
                if site == "uterus" and "endometrium" in buckets:
                    buckets["endometrium"]["gtex_v11_shards"].add(f"{site} ({n_s})")
                    buckets["endometrium"]["gtex_v11_medians"].add(col)

    # 2. source-tissues.csv rows
    for r in st.itertuples(index=False):
        cid = label_to_cid.get(r.tissue)
        if not cid or cid not in buckets:
            continue
        fam = str(r.family)
        ctx = str(r.source_context or r.tissue)
        if fam in ("rna_tissue_gtex", "rna_brain_gtex"):
            buckets[cid]["gtex_via_hpa_v8"].add(ctx)
        elif fam == "hpa_tissue_samples":
            buckets[cid]["hpa_surgical_biopsies"].add(ctx)
        elif fam in ("rna_tissue_hpa", "rna_brain_hpa", "hpa_brain"):
            buckets[cid]["hpa_tissue_rna"].add(ctx)
        elif fam == "hpa_rna_consensus":
            buckets[cid]["hpa_rna_consensus"].add(ctx)
        elif fam == "hpa_single_cell_tissue":
            buckets[cid]["hpa_single_cell"].add(str(r.tissue))
        elif fam == "census":
            buckets[cid]["cellxgene_census"].add(ctx)
        elif fam in ("dice", "hpa_immune_monaco", "hpa_immune_cell", "hpa_immune_cell_samples"):
            buckets[cid]["sorted_immune"].add(fam)
        elif fam in ("cornea_collin_adult_counts", "petropoulos_E-MTAB-3929", "tyser_E-MTAB-9388"):
            buckets[cid]["specialized_and_developmental"].add(fam)

    # Also record Mole peri-implantation embryo
    if "embryo" in buckets:
        buckets["embryo"]["specialized_and_developmental"].add("mole_E-MTAB-8060")

    # 3. #735 unmeasured home tissue panels
    hp = unmeasured_home_panels.unmeasured_home_tissue_panels()
    for r in hp.itertuples(index=False):
        cid = label_to_cid.get(str(r.home_tissue))
        if cid and cid in buckets:
            buckets[cid]["specialized_and_developmental"].add(f"{r.panel_id} ({r.modality})")

    # 4. HPA IHC normal_tissue mappings
    for row in normal_tissues.normal_tissue_mappings(
        species="human", source_id="hpa_normal_tissue", source_version="v23"
    ):
        cid = row["concept_id"]
        if cid in buckets and row["source_tissue"]:
            buckets[cid]["hpa_ihc_protein"].add(row["source_tissue"])

    out_rows = []
    for cid in sorted(buckets):
        b = buckets[cid]
        info = cid_to_info[cid]
        has_gtex = bool(b["gtex_v11_shards"] or b["gtex_v11_medians"] or b["gtex_via_hpa_v8"])
        has_hpa_bulk = bool(b["hpa_surgical_biopsies"] or b["hpa_tissue_rna"])
        has_sc = bool(b["hpa_single_cell"] or b["cellxgene_census"] or b["sorted_immune"])
        has_spec = bool(b["specialized_and_developmental"])
        has_ihc = bool(b["hpa_ihc_protein"])
        n_indep = sum((has_gtex, has_hpa_bulk, has_sc, has_spec, has_ihc))
        if n_indep == 0:
            continue

        notes = []
        if cid == "retina":
            notes.append(
                "Absent from 52 GTEx v11 shards; EyeGEx distributed via rna_tissue_gtex/HPA"
            )
        elif not b["gtex_v11_shards"] and b["gtex_via_hpa_v8"]:
            notes.append("Absent as a separate GTEx v11 shard; mapped via rna_tissue_gtex/HPA")
        elif not b["gtex_v11_shards"]:
            notes.append("Absent from GTEx v11; covered by orthogonal non-GTEx cohorts")
        if cid == "brain":
            notes.append("Rolls up all 13 regional GTEx v11 CNS shards")
        elif cid == "endometrium":
            notes.append(
                "Sampled within whole-organ GTEx v11 uterus shard; resolved separately in HPA"
            )
        if b["gtex_v11_subcompartments"]:
            notes.append(
                f"GTEx v11 adds {len(b['gtex_v11_subcompartments'])} histological subcompartment median column(s)"
            )

        out_rows.append(
            {
                "concept_id": cid,
                "tissue_label": info["tissue_label"],
                "tissue_category": info["tissue_category"],
                "gtex_v11_shards": ";".join(sorted(b["gtex_v11_shards"])),
                "gtex_v11_medians": ";".join(sorted(b["gtex_v11_medians"])),
                "gtex_v11_subcompartments": ";".join(sorted(b["gtex_v11_subcompartments"])),
                "gtex_via_hpa_v8": ";".join(sorted(b["gtex_via_hpa_v8"])),
                "hpa_surgical_biopsies": ";".join(sorted(b["hpa_surgical_biopsies"])),
                "hpa_tissue_rna": ";".join(sorted(b["hpa_tissue_rna"])),
                "hpa_rna_consensus": ";".join(sorted(b["hpa_rna_consensus"])),
                "hpa_single_cell": ";".join(sorted(b["hpa_single_cell"])),
                "cellxgene_census": ";".join(sorted(b["cellxgene_census"])),
                "sorted_immune": ";".join(sorted(b["sorted_immune"])),
                "specialized_and_developmental": ";".join(
                    sorted(b["specialized_and_developmental"])
                ),
                "hpa_ihc_protein": ";".join(sorted(b["hpa_ihc_protein"])),
                "independent_modality_families": n_indep,
                "notes": "; ".join(notes),
            }
        )

    grid_df = pd.DataFrame(out_rows)
    grid_path = ROOT / "oncoref/data/normal-tissue-source-grid.csv"
    grid_df.to_csv(grid_path, index=False, lineterminator="\n")
    print(f"Wrote {len(grid_df)} canonical tissues to {grid_path.name}")


if __name__ == "__main__":
    if not (ROOT / "oncoref/data/gtex-donor-sex.csv").is_file() or "--rebuild-sex" in sys.argv:
        build_donor_sex()
    build_source_grid()
