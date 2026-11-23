#!/usr/bin/env python3
"""Audit the exact adult-reproductive-bulk-zero tumor screen, from pinned sources."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_cta_candidate_evidence import replay_updates  # noqa: E402
from cta_report_common import input_hashes, sha256, write_csv  # noqa: E402

from oncoref import gene_knowledge  # noqa: E402

OUT = ROOT / "docs/audits/cancer-only-20261109"
BASE = ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz"
CQ = "pone.0096333.s002.xlsx"
COPIES = "pone.0096333.s004.xlsx"
NEURONS = "41593_2020_587_MOESM3_ESM.xlsx"
VERBEURGT = "https://doi.org/10.1371/journal.pone.0096333"
DURANTE = "https://doi.org/10.1038/s41593-020-0587-9"
FETAL = "41467_2026_71595_MOESM3_ESM.xlsx"
FETAL_URL = "https://doi.org/10.1038/s41467-026-71595-6"
PAIRS = {"OR52E6": "OR52E8", "OR52E8": "OR52E6", "OR5D18": "OR5D16"}


def native_assays(inputs):
    """The 26 donor columns are fixed by sample metadata, never inferred by dtype."""
    tables = {}
    for name in [CQ, COPIES]:
        raw = pd.read_excel(inputs / name, header=None)
        if raw.iloc[1, 2:28].tolist() != list(range(1, 27)):
            raise ValueError("Changed qPCR donor-column contract")
        rows = raw.iloc[5:]
        lookup = {}
        for _, row in rows.iterrows():
            label = str(row.iloc[1])
            if "-Hs" not in label:
                continue
            symbol, assay = label.split("-", 1)
            if not symbol.startswith("OR"):
                continue
            probes = lookup.setdefault(symbol, {})
            if assay in probes:
                raise ValueError("Duplicate exact-symbol qPCR assay")
            probes[assay] = pd.to_numeric(row.iloc[2:28], errors="raise")
        tables[name] = lookup
    neurons = pd.read_excel(inputs / NEURONS, sheet_name="SI Table 7", header=1)
    cells = pd.read_excel(inputs / NEURONS, sheet_name="SI Table 3", header=1)
    return tables, neurons.set_index("OR gene").iloc[:, 0].to_dict(), cells


def hpa_evidence(path):
    entry = ET.fromstring(gzip.decompress(path.read_bytes())).find("entry")
    antibodies, observations, rna = [], [], []
    overall = entry.findtext("tissueExpression/verification", "Not reported")
    for antibody in entry.findall("antibody"):
        aid = antibody.attrib["id"]
        for assay in antibody.findall("tissueExpression"):
            if assay.attrib.get("assayType") != "tissue":
                continue
            reliability = overall
            antibodies.append(
                {
                    "antibody": aid,
                    "reliability": reliability,
                    "antibody_validation": assay.findtext("verification", "Not reported"),
                }
            )
            for data in assay.findall("data"):
                tissue = data.findtext("tissue")
                for cell in data.findall("tissueCell"):
                    level = cell.find("level")
                    if level is not None:
                        observations.append(
                            {
                                "antibody": aid,
                                "tissue": tissue,
                                "cell": cell.findtext("cellType"),
                                "level": level.text,
                                "reliability": reliability,
                            }
                        )
    for assay in entry.findall("rnaExpression"):
        if "pig" in assay.attrib.get("assayType", "").lower():
            continue
        for data in assay.findall("data"):
            for level in data.findall("level"):
                if level.attrib.get("type") == "normalizedRNAExpression":
                    rna.append(
                        {
                            "assay": assay.attrib.get("assayType"),
                            "tissue": data.findtext("tissue"),
                            "value": float(level.attrib["expRNA"]),
                            "unit": level.attrib["unitRNA"],
                        }
                    )
    ms = entry.find("proteinEvidence/evidence[@source='MS']")
    return (
        antibodies,
        observations,
        rna,
        ms.attrib["evidence"] if ms is not None else "Not reported",
    )


NOTES = {
    "OTOR": "Human fetal cochlear cloning and mixed human/avian Northern evidence identify an inner-ear question; they do not establish adult human protein. Retain the stage-specific normal hold.",
    "TAAR6": "Human normal-brain RT-PCR supplies the existing hold. Separate exact human TAAR6 from rodent receptor studies and schizophrenia association.",
    "TAAR8": "Human and rodent expression experiments must be separated. Recombinant human TAAR8 surface expression in COS-7 cells verifies a construct, not normal human tissue production.",
    "TGM6": "The neural/ataxia annotation and uncertain HPA neuronal staining require direct adult human localization; neither genetic disease nor rodent antibody experiments alone establish human protein.",
    "IFNW1": "Activated normal immune-cell protein supplies the existing Unsafe finding. Resting bulk RNA misses induced production. Secretion does not identify every producing cell.",
    "IFNA4": "Activated normal immune-cell protein supplies the existing Unsafe finding. Highly homologous interferons and a shared HPA reagent require assay-specific attribution.",
    "MRGPRX4": "Human dorsal-root-ganglion neuronal RNA and protein supply the existing Unsafe finding. The rare sensory-neuron compartment is absent from the ordinary bulk comparison.",
    "OC90": "Inner-ear otoconial biology is strongly established in mice. The principal source cloned human PLA2L/OC90 but obtained native protein from mouse tissue; retain the human localization hold.",
    "TAS2R42": "Gene-specific adult human skin RNA supplies a normal-expression hold. The skin paper does not establish TAS2R42-specific protein, and other TAS2R staining cannot be borrowed.",
    "KRTAP21-1": "Hair-family membership identifies the compartment to examine but is not exact-protein proof. HPA lists MS protein-existence evidence without tissue or peptide attribution in its XML, whereas the retrieved UniProt entry remains homology-based. Zhang 2014 reports gene-specific qPCR in LoVo cancer cells; the paper's HCEC experiments are not a KRTAP21-1 normal-cell result.",
    "SCYGR2": "HPA assigns a hair-follicle-development RNA cluster and imports protein-existence evidence. This does not identify a tissue-specific endogenous peptide for the selected proteoform; exact normal protein localization remains unverified.",
    "SCYGR8": "The former KRTAP28-8 name and predicted structural role identify a skin/hair question. Homology and recombinant products do not establish endogenous human hair expression; inspect the exact selected sequence.",
    "USP17L13": "The current HPA antibody is uncertain and shared. A 2026 influenza study uses A549 cancer cells, CRISPR activation and overexpression; cytokine induction there cannot establish normal lung-cell production.",
    "USP17L8": "The HPA reagent is shared with related copies and gives uncertain normal staining. USP17L2/DUB3 functional literature cannot be transferred to USP17L8; retain an exact-copy normal-expression hold.",
    "CLDN25": "Durante reports enriched adult human olfactory-cell RNA, supplying the existing hold. Claudin-family function does not establish the localization or indispensability of this exact protein.",
    "OR11H1": "Zhao 2013 reports RNA in healthy-control PBMC and neurologically normal postmortem brain, with qPCR confirmation. This is normal non-olfactory RNA evidence; OR4M1 functional experiments do not establish OR11H1 protein or function.",
    "OR51B6": "The human leukocyte study includes OR51B6 in an OR56B4/OR52A4/OR51B6 antibody mixture. The positive mixed stain leaves exact-receptor attribution unresolved and supplies a review hold rather than an exact-protein exclusion.",
    "OR1D5": "The Parkinson nasal-marker paper validates OR1D5 RNA in treated RPMI-2650 carcinoma cells. Its patient/control nasal-fluid validation concerns other genes; it does not establish normal OR1D5 production.",
    "OR7G3": "A schizophrenia leukocyte RNA-seq study is an exact-symbol literature lead. A normal-control gene-specific result was not extracted, so this study supplies no normal-positive or normal-negative credit here.",
    "OR2B2": "The linked-receptor testicular-expression paper found in the search studies pigs. It cannot establish reproductive expression of human OR2B2. No exact-symbol assay occurs in the human olfactory qPCR table.",
}
EXTRA_SOURCES = {
    "TAS2R42": ["https://pubmed.ncbi.nlm.nih.gov/30332676/"],
    "OR5D18": ["https://doi.org/10.1053/j.gastro.2007.02.036"],
    "KRTAP21-1": [
        "https://doi.org/10.3892/ol.2014.2151",
        "https://doi.org/10.1186/s12862-023-02107-z",
    ],
    "USP17L13": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC13211412/"],
    "TAAR8": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC4264187/"],
    "OC90": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC28045/"],
    "OR11H1": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC3968322/"],
    "OR51B6": ["https://pubmed.ncbi.nlm.nih.gov/25624459/"],
    "OR1D5": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC7408990/"],
    "OR7G3": ["https://pubmed.ncbi.nlm.nih.gov/32038049/"],
    "OR2B2": ["https://pubmed.ncbi.nlm.nih.gov/39915713/"],
    "CLDN25": ["https://doi.org/10.1093/hmg/ddae038", "https://doi.org/10.1093/hmg/ddaf192"],
}
SOURCE_TITLES = {
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3968322/": "Zhao 2013: PBMC and normal brain olfactory-receptor RNA",
    "https://pubmed.ncbi.nlm.nih.gov/25624459/": "Class I odorant and taste receptors in circulating human leukocytes (2015)",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC7408990/": "Nasal epithelial Parkinson markers: carcinoma model and clinical validation (2020)",
    "https://pubmed.ncbi.nlm.nih.gov/32038049/": "Schizophrenia leukocyte RNA-seq (2020; discovery lead)",
    "https://pubmed.ncbi.nlm.nih.gov/39915713/": "MHC-linked receptors and pig testis (2025; species mismatch)",
    "https://doi.org/10.1093/hmg/ddae038": "CLDN25 variant and barrier-cell models (2024)",
    "https://doi.org/10.1093/hmg/ddaf192": "Correction to the CLDN25 variant study (2026)",
    VERBEURGT: "Verbeurgt 2014: normal human olfactory-mucosa qPCR",
    DURANTE: "Durante 2020: adult human olfactory single-cell atlas",
    FETAL_URL: "Mbouamboua 2026: fetal human nasal atlas",
    "https://pubmed.ncbi.nlm.nih.gov/30332676/": "Personalized expression of bitter receptors in human skin (2018)",
    "https://doi.org/10.1053/j.gastro.2007.02.036": "Braun 2007: human gut enterochromaffin-cell receptors",
    "https://doi.org/10.3892/ol.2014.2151": "Zhang 2014: parvovirus B19 and colorectal cancer cells",
    "https://doi.org/10.1186/s12862-023-02107-z": "Litman and Stein 2023: KRTAP evolution and skin transcriptomics",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC13211412/": "USP17L13 and influenza replication in A549 cells (2026)",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4264187/": "Human TAAR8 and murine Taar8b signaling and expression (2014)",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC28045/": "Otoconin-90 cloning and mammalian otoconial matrix protein (1998)",
}


def review_finding(row, suffix, finding, assay, url, location, limitations, value=None, unit=None):
    return {
        "id": f"GK-20261109-{row.protein_group_label}-{suffix}",
        "gene_ids": json.loads(row.member_gene_ids),
        "gene_symbols": [row.protein_group_label],
        "taxonomy_id": 9606,
        "stage": "adult",
        "tissue": "olfactory mucosa" if suffix == "olfactory" else "source-specified normal tissue",
        "cell_type": "human olfactory mucosa; producing cell unresolved"
        if suffix == "olfactory"
        else "source-specified normal cells",
        "state": "normal tissue",
        "assay": assay,
        "support": "unresolved_normal_expression",
        "tissue_tier": "unresolved",
        "finding": finding,
        "value": value,
        "unit": unit,
        "effect": "review",
        "interpretation_status": "unresolved",
        "decision_reason": "Normal expression or its exact-copy attribution requires review before CTA or tumor-restriction clearance.",
        "source_title": "Human normal-expression evidence: exact-gene cancer-only audit",
        "source_url": url,
        "source_location": location,
        "limitations": limitations,
        "review_status": "accepted",
        "review_date": "2026-11-09",
        "historical_finding_ids": [],
        "review_depth": "Individual source/assay review; preserved primary tables and HPA responses",
    }


def build(apply_curation=False):
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = OUT / "inputs"
    facts = pd.read_csv(BASE, low_memory=False)
    pool = facts[facts.tumor_restricted_profile]
    selected = pool[pool.reproductive_ntpm.eq(0)]
    tables, neurons, _cells = native_assays(inputs)
    fetal = pd.read_excel(inputs / FETAL, sheet_name="Supplementary Data 8", header=3)
    fetal = fetal.set_index("OR genes")["Nb cells expressing OR"].to_dict()
    knowledge = gene_knowledge.findings(include_updates=True)
    records, proposals, assay_rows, current_rna = [], [], [], []
    for row in selected.sort_values("protein_group_label").itertuples():
        name, gid = row.protein_group_label, row.representative_gene_id
        search = json.loads(
            gzip.decompress((inputs / f"{name}.literature-search.json.gz").read_bytes())
        )
        antibodies, staining, rna, ms = hpa_evidence(inputs / f"{name}.xml.gz")
        current_rna.extend({"protein": name, **v} for v in rna)
        cq, copies = tables[CQ].get(name), tables[COPIES].get(name)
        if cq and set(cq) != set(copies):
            raise ValueError("Cq/copy assay identities differ")
        donors = int(pd.DataFrame(copies).ge(5).any(axis=1).sum()) if copies else None
        # Cq <35 is an assay reproducibility statistic. The paper's expression
        # count uses >=5 estimated copies/20 ng RNA, which is a different test.
        for aid in cq or {None: None}:
            assay_rows.append(
                {
                    "protein": name,
                    "assay_id": aid,
                    "donors_cq_lt35": int(cq[aid].lt(35).sum()) if cq else None,
                    "donors_ge5_copies": int(copies[aid].ge(5).sum()) if copies else None,
                    "donors_measured": 26 if copies else None,
                    "max_copies_per_20ng": float(copies[aid].max()) if copies else None,
                    "cross_reactive_partner": PAIRS.get(name),
                    "positive_neurons": neurons.get(name),
                    "reported_fetal_cells": fetal.get(name),
                }
            )
        notes = [
            NOTES.get(
                name,
                "Exact-copy olfactory expression is the principal normal-tissue question. Odorant responses in transfected cells establish receptor activity, not endogenous normal or tumor protein.",
            )
        ]
        sources = [f"https://www.proteinatlas.org/{gid}/tissue", *EXTRA_SOURCES.get(name, [])]
        if not antibodies:
            notes.append(
                "The current HPA response contains no tissue-IHC antibody assay for this protein."
            )
        if ms == "Not available":
            notes.append(
                "HPA's MS field is unavailable; it supplies no protein-negative measurement."
            )
        if name in fetal:
            notes.append(
                f"The 2026 fetal nasal atlas reports {fetal[name]} cells expressing this exact gene label at 7-12 post-conceptional weeks. These counts are not TPM, protein or copy-specific peptide evidence."
            )
            sources.append(FETAL_URL)
        if donors is not None:
            sources.append(VERBEURGT)
            notes.append(
                f"Verbeurgt Table S3 reports {donors}/26 donors at >=5 estimated copies per 20 ng RNA; the raw Cq <35 count is recorded separately."
            )
            if len(copies) > 1:
                notes.append(
                    "Multiple assays are retained individually. The gene-level donor count is their union, with each donor counted once."
                )
            if name in PAIRS:
                notes.append(
                    f"The authors cannot discriminate {name} from {PAIRS[name]} with this assay. This is shared-assay evidence."
                )
            if donors:
                proposals.append(
                    review_finding(
                        row,
                        "olfactory",
                        f"Normal human olfactory-mucosa assay labeled {name}: {donors}/26 donors at >=5 estimated RNA copies/20 ng.",
                        "TaqMan RT-qPCR and calibrated estimated RNA copies",
                        VERBEURGT,
                        "Tables S1/S3; 26 donor columns; assay-specificity controls",
                        "Native copy estimates are not TPM or protein. Exact-symbol matching only; incomplete receptor coverage and donor demographics limit negative conclusions."
                        + (f" The assay cross-reacts with {PAIRS[name]}." if name in PAIRS else ""),
                        donors,
                        "donors with >=5 estimated RNA copies/20 ng",
                    )
                )
        else:
            notes.append(
                "No exact-symbol assay occurs in the qPCR table; this supplies no negative measurement. Ambiguous legacy names are not reassigned."
            )
        if name in neurons:
            notes.append(
                f"Durante Table 7 reports RNA in {neurons[name]} adult human olfactory sensory neurons."
            )
            sources.append(DURANTE)
        positive = [v for v in staining if str(v["level"]).lower() in {"low", "medium", "high"}]
        if positive:
            if any(v["reliability"].lower() != "uncertain" for v in positive):
                raise ValueError(
                    "Review changed antibody reliability before assigning a consequence"
                )
            notes.append(
                "Current HPA has positive normal staining with uncertain reagent reliability; it cannot establish exact protein production or absence."
            )
            proposals.append(
                review_finding(
                    row,
                    "HPA",
                    "Positive normal-tissue HPA staining remains unresolved under uncertain antibody validation.",
                    "HPA immunohistochemistry; " + ";".join(v["antibody"] for v in antibodies),
                    sources[0],
                    "Saved XML: antibody/tissueExpression/data; validation labels",
                    "Shared or uncertain antibodies cannot prove exact-copy production. Categorical staining is not TPM; the current response is separate from pinned HPA v23.",
                )
            )
            proposals[-1]["tissue"] = "; ".join(sorted({v["tissue"] for v in positive}))
            proposals[-1]["cell_type"] = "; ".join(sorted({v["cell"] for v in positive}))
            proposals[-1]["stage"] = (
                "placental"
                if all(v["tissue"].lower() == "placenta" for v in positive)
                else "adult and/or placental; see tissue"
            )
            if proposals[-1]["stage"] == "placental":
                proposals[-1].update(
                    effect="context",
                    tissue_tier="reproductive",
                    support="unresolved_reproductive_protein",
                    decision_reason="Uncertain placental staining leaves reproductive absence unresolved; it supplies no adult normal-cell exclusion or hold.",
                )
        if name == "TAS2R42":
            proposals.append(
                review_finding(
                    row,
                    "skin",
                    "Gene-specific TAS2R42 RNA is reported in normal adult human skin biopsies and donor skin RNA data.",
                    "Gene-specific RNA assays in adult human skin",
                    EXTRA_SOURCES[name][0],
                    "Abstract; skin expression and sun-exposure analysis",
                    "RNA does not establish exact-protein abundance, producing skin cell or therapeutic injury.",
                )
            )
            proposals[-1].update(tissue="skin", cell_type="producing skin cell unresolved")
        if name == "OR5D18":
            notes.append(
                "Braun 2007 reports OR73/OR5D18 in laser-microdissected normal human gut enterochromaffin cells; cultured-cell functional experiments are separate."
            )
        if name in {"OR11H1", "OR51B6"}:
            text = NOTES[name]
            proposals.append(
                review_finding(
                    row,
                    "blood",
                    text,
                    "RT-PCR/qPCR and normal brain RNA"
                    if name == "OR11H1"
                    else "Mixed-receptor immunocytochemistry in human blood PMN",
                    EXTRA_SOURCES[name][0],
                    "Figures 2-5 and normal-control cohorts"
                    if name == "OR11H1"
                    else "Figure 4: 1:1:1 receptor-antibody cocktail",
                    "RNA is not protein; a mixed stain cannot establish exact-copy protein. No native TPM or clinical injury is inferred.",
                )
            )
            proposals[-1].update(
                tissue="blood; brain" if name == "OR11H1" else "blood",
                cell_type="PBMC; postmortem brain producing cells unresolved"
                if name == "OR11H1"
                else "polymorphonuclear leukocytes",
            )
        if name == "CLDN25":
            notes.append(
                "The 2024 variant study (with a 2026 correction) examines a patient and transfected HEK293/canine MDCKII/mouse endothelial models. It adds a barrier-function question, without establishing adult human neural protein or cytotoxic-injury causality."
            )
        if name == "OR5D18":
            proposals.append(
                review_finding(
                    row,
                    "gut",
                    "OR73/OR5D18 RNA is reported in laser-microdissected human mucosal enterochromaffin cells.",
                    "Localized RT-PCR; human enterochromaffin cells",
                    EXTRA_SOURCES[name][0],
                    "Methods and receptor-expression results; OR73 historical label",
                    "Normal-cell RNA and cancer-cell functional assays are distinct. This review does not promote cell-line function to normal human protein evidence.",
                )
            )
            proposals[-1].update(tissue="intestine", cell_type="enterochromaffin cells")
        ids = set(json.loads(row.member_gene_ids))
        old = knowledge[knowledge.gene_ids.map(lambda v, ids=ids: bool(ids & set(v)))]
        sources.extend(old.source_url.tolist())
        notes.extend(old.finding.tolist())
        records.append(
            {
                "protein": name,
                "gene_id": gid,
                "entity_id": row.protein_equivalence_group_id,
                "released_classification": row.classification,
                "reproductive_bulk_ntpm": row.reproductive_ntpm,
                "reproductive_cell_ntpm": None
                if pd.isna(row.reproductive_cell_ntpm)
                else row.reproductive_cell_ntpm,
                "developmental_max_tpm": row.developmental_all_lineage_max_tpm,
                "tumor_tpm_max": row.oncoref_maximum_tpm,
                "tumor_support_samples_minimum": row.tumor_samples_ge10,
                "normal_ntpm_max": row.healthy_maximum_ntpm,
                "normal_tpm_max": row.healthy_maximum_tpm,
                "translation_qualified": bool(row.translation_qualified),
                "translation_qualification_reasons": row.translation_qualification_reasons
                if pd.notna(row.translation_qualification_reasons)
                else "",
                "representative_selection_approved": bool(row.representative_selection_approved),
                "hpa_antibodies": antibodies,
                "hpa_positive_staining": positive,
                "hpa_ms": ms,
                "reported_fetal_cells": fetal.get(name),
                "literature_search": {
                    "query": search["audit_query"],
                    "hits": search["hitCount"],
                    "retrieved": len(search["resultList"]["result"]),
                    "limit": 100,
                },
                "notes": notes,
                "sources": list(dict.fromkeys(sources)),
                "confirmed_cancer_only": False,
            }
        )
    gene_knowledge.validate(proposals)
    (OUT / "curation-proposals.json").write_text(json.dumps(proposals, indent=2) + "\n")
    if apply_curation:
        path = ROOT / "oncoref/data/gene-knowledge-updates.json"
        data = json.loads(path.read_text())
        data["findings"] = [
            r for r in data["findings"] if not r["id"].startswith("GK-20261109-")
        ] + proposals
        gene_knowledge.validate(data["findings"])
        path.write_text(json.dumps(data, indent=2) + "\n")
    revised = replay_updates(facts).classifications.set_index("protein_equivalence_group_id")
    for record in records:
        r = revised.loc[record["entity_id"]]
        record.update(
            classification=r.classification,
            assessment_status=r.assessment_status,
            admission_reasons=r.exclusion_reasons if pd.notna(r.exclusion_reasons) else "",
        )
        if not record["translation_qualified"]:
            record["notes"].append(
                "Translation review: the reference cannot supply a complete approved protein. "
                "Recorded failures: " + record["translation_qualification_reasons"] + "."
            )
    write_csv(pd.DataFrame(assay_rows), OUT / "olfactory-assays.csv")
    write_csv(pd.DataFrame(current_rna), OUT / "current-hpa-rna.csv.gz")
    write_csv(pool, OUT / "all-66-screen-matches.csv.gz")
    peak = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/normal-peak-evidence.csv.gz", low_memory=False
    )
    peak = peak[peak.gene_id.isin(selected.representative_gene_id)]
    top = (
        peak.sort_values("value", ascending=False)
        .groupby(["gene_id", "source", "unit"], dropna=False)
        .head(5)
    )
    write_csv(top, OUT / "normal-top-contexts.csv.gz")
    tumor_path = ROOT / "docs/audits/cta-selection-20261028/tumor-restricted-observations.csv.gz"
    tumor = pd.read_csv(tumor_path)
    tumor = tumor[tumor.protein_equivalence_group_id.isin(selected.protein_equivalence_group_id)]
    write_csv(tumor, OUT / "supporting-tumor-observations.csv")
    for record in records:
        record["normal_contexts"] = json.loads(
            top[top.gene_id.eq(record["gene_id"])][
                ["source", "tissue", "cell_type", "value", "unit", "Read count"]
            ].to_json(orient="records")
        )
        record["tumor_observations"] = json.loads(
            tumor[tumor.protein_equivalence_group_id.eq(record["entity_id"])].to_json(
                orient="records"
            )
        )
    write_csv(
        pd.DataFrame([{"protein": r["protein"], **r["literature_search"]} for r in records]),
        OUT / "literature-search-ledger.csv",
    )
    artifact = {
        "schema": "oncoref.cancer-only-audit",
        "version": 1,
        "reviewed_on": "2026-11-09",
        "baseline_release": "1.35.0",
        "source_titles": {
            **dict(zip(knowledge.source_url, knowledge.source_title, strict=True)),
            **SOURCE_TITLES,
        },
        "definition": {
            "numerical_screen": "Existing normal-low tumor profile with >=2 distinct native-TPM samples >=10 or reviewed tumor protein; no qualifying reproductive/developmental association.",
            "audited_subset": "Exactly zero adult reproductive bulk nTPM. Missing reproductive-cell measurements remain missing.",
            "cancer_only": "No confirmed cancer-only member. A normal-screen cutoff does not establish absence in all normal tissues or states.",
            "category_decision": "Retain CTA / Other / Unsafe. Cancer support and reproductive absence remain independent evidence fields; this screen does not establish a fourth category.",
        },
        "counts": {
            "numerical_screen": len(pool),
            "adult_reproductive_bulk_zero": len(selected),
            "adult_reproductive_bulk_nonzero": int(pool.reproductive_ntpm.gt(0).sum()),
            "zero_developmental_rna": int(selected.developmental_all_lineage_max_tpm.eq(0).sum()),
            "confirmed_cancer_only": 0,
            "reviewed_calls": pd.Series([r["classification"] for r in records])
            .value_counts()
            .to_dict(),
        },
        "records": records,
    }
    raw = json.dumps(artifact, indent=2, allow_nan=False).encode()
    (OUT / "audit.json").write_bytes(raw + b"\n")
    (ROOT / "oncoref/data/cancer-only-audit.json.gz").write_bytes(gzip.compress(raw, mtime=0))
    flat = pd.DataFrame(records)
    for c in [
        "notes",
        "sources",
        "hpa_antibodies",
        "hpa_positive_staining",
        "normal_contexts",
        "tumor_observations",
        "literature_search",
    ]:
        flat[c] = flat[c].map(json.dumps)
    write_csv(flat, OUT / "protein-reviews.csv")
    render(artifact)
    paths = [
        Path(__file__),
        BASE,
        tumor_path,
        ROOT / "docs/audits/cta-selection-20261028/normal-peak-evidence.csv.gz",
        ROOT / "scripts/build_cta_candidate_evidence.py",
        ROOT / "oncoref/cta_decision_flow.py",
        ROOT / "oncoref/gene_knowledge.py",
        ROOT / "oncoref/data/gene-knowledge.json",
        ROOT / "oncoref/data/gene-knowledge-updates.json",
        *sorted(inputs.iterdir()),
    ]
    manifest = {
        "inputs": input_hashes(paths),
        "outputs": {
            p.name: sha256(p)
            for p in sorted(OUT.iterdir())
            if p.is_file() and p.name != "manifest.json"
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(artifact["counts"], indent=2))


def render(artifact):
    import matplotlib

    matplotlib.use("Agg")
    import markdown
    import matplotlib.pyplot as plt

    plt.rcParams["svg.hashsalt"] = "cancer-only-20261109"
    fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
    ax.barh(
        [
            "Low normal RNA / high tumor RNA",
            "Zero adult reproductive bulk RNA",
            "Confirmed cancer-only",
        ],
        [66, 47, 0],
        color=["#356e81", "#d29b45", "#4b7868"],
    )
    ax.invert_yaxis()
    ax.set(
        xlabel="Proteins",
        xlim=(0, 72),
        title="Screening matches and confirmed cancer-only proteins",
    )
    for i, n in enumerate([66, 47, 0]):
        ax.text(n + 1, i, str(n), va="center")
    for suffix in ["png", "svg"]:
        fig.savefig(
            OUT / f"selection.{suffix}", dpi=150, metadata={"Date": None} if suffix == "svg" else {}
        )
    plt.close(fig)
    lines = [
        "# Cancer-only expression: individual normal-tissue audit\n",
        "**No protein in this screen qualifies as Cancer-Only.** The existing numerical screen has 66 matches. Forty-seven have zero adult reproductive bulk RNA and receive individual reviews below; the other 19 have measurable reproductive bulk RNA. All 47 have some recorded developmental RNA, and reproductive-cell coverage is incomplete. These small developmental values are observations, not proof of translated fetal protein.\n",
        "**Keep CTA, Other and Unsafe.** Use cancer expression and reproductive/developmental expression as separate table fields. The numerical screen mixes normal sensory receptors, inducible immune proteins, uncertain assays and poorly characterized copies. It does not define a coherent fourth category.\n",
        "Reviewed calls: "
        + ", ".join(f"**{n} {label}**" for label, n in artifact["counts"]["reviewed_calls"].items())
        + ".\n",
        "![Selection counts](selection.png)\n",
        "[Filterable individual reviews](table.html) · [All 47 reviews](protein-reviews.csv) · [All 66 numerical matches](all-66-screen-matches.csv.gz) · [Top five normal contexts per source](normal-top-contexts.csv.gz) · [Exact olfactory assay results](olfactory-assays.csv) · [Current HPA RNA](current-hpa-rna.csv.gz) · [Provenance](manifest.json)\n",
        "## How the set was selected\n",
        "The starting set requires the existing normal-low profile and at least two distinct tumors at >=10 native TPM, or reviewed tumor protein. The saved tumor count stops after two supporting observations: it is a **lower bound**, not prevalence. Zero adult reproductive bulk RNA defines the 47-protein audit, rather than calling a value below 5 nTPM absent. Placenta is included in the pinned reproductive bulk definition. Missing cell assays are never filled with zero.\n",
        "## Assay review\n",
        "Every protein has a retained exact-symbol [literature-search query and result count](literature-search-ledger.csv), capped at 100 retrieved records per query. Search hits identify leads; accepted findings below require review of the human tissue, assay and exact gene. Full abstracts are not redistributed in the query metadata.\n",
        f"[Verbeurgt 2014]({VERBEURGT}) supplies measured RNA-copy estimates in 26 human olfactory-mucosa donors. Expression counts below use the authors' >=5 copies per 20 ng RNA rule. Raw Cq <35 is recorded separately. OR52E6/OR52E8 and OR5D18/OR5D16 cross-react. OR5D18 has zero donors above the copy threshold in this study despite seven Cq values below 35; its gut-cell evidence is independent. [Durante 2020]({DURANTE}) supplies observed olfactory-neuron counts, not a complete negative receptor panel. Neither experiment proves protein production.\n",
        "HPA responses were retrieved separately from pinned v23 selection data. Uncertain or shared IHC supplies a review question; unavailable IHC/MS supplies no negative result. An ectopic-expression experiment, recombinant protein, animal ortholog or family-wide annotation cannot establish exact human normal protein. New reviewed normal findings enter the same knowledge gate used for every protein.\n",
        "| Source | Role |\n|---|---|\n"
        f"| [Verbeurgt 2014]({VERBEURGT}) | Donor-level olfactory qPCR, calibrated copy estimates and cross-reactivity controls. |\n"
        f"| [Durante 2020]({DURANTE}) | Positive adult olfactory-neuron counts and cell-enrichment tables. |\n"
        f"| [Mbouamboua 2026]({FETAL_URL}) | Exact gene labels and positive fetal nasal-cell counts; six matches in this set. |\n"
        "| [Human Protein Atlas](https://www.proteinatlas.org/) | Forty-seven saved XML responses: current RNA, normal-tissue IHC validation and MS evidence labels. Cancer staining and pig RNA are excluded from normal-human observations. |\n"
        "| [UniProt](https://www.uniprot.org/) | Five saved gene-specific responses for SCYGR2/8, KRTAP21-1 and USP17L8/13: protein-existence annotations and primary-reference leads. |\n"
        "| [Song 2026 adult olfactory atlas](https://doi.org/10.1038/s41467-026-75722-1) | Bibliographic discovery retained; processed gene-level tables were not retrieved. It contributes no positive or negative gene-level measurement here. |\n",
        "## Individual decisions\n",
    ]
    for r in artifact["records"]:
        name = r["protein"]
        lines.extend(
            [
                f"### {name}\n",
                f"**{r['classification']}** · [Ensembl](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r['gene_id']}) · [HPA](https://www.proteinatlas.org/{r['gene_id']}/tissue)\n",
                f"Pinned normal maxima: {r['normal_ntpm_max']:.3g} nTPM; {r['normal_tpm_max']:.3g} TPM. Adult reproductive bulk: 0 nTPM. Developmental maximum: {r['developmental_max_tpm']:.4g} TPM. Cancer maximum: {r['tumor_tpm_max']:.4g} TPM.\n",
                f"Admission reasons: {r['admission_reasons']}. Assessment status: {r['assessment_status']}.\n",
                *[f"- {note}" for note in r["notes"]],
                "\nSupporting tumor observations: "
                + "; ".join(
                    f"{o['cancer_code']} / {o['sample_id']}: {o['value']:.3g} TPM"
                    for o in r["tumor_observations"]
                )
                + ". Full source hashes and sample identities are in [the observation table](supporting-tumor-observations.csv).\n",
                "\n"
                + " · ".join(f"[{source_title(url, artifact)}]({url})" for url in r["sources"])
                + "\n",
            ]
        )
    report = "\n".join(lines)
    (OUT / "index.md").write_text(report)
    style = "body{font:16px/1.55 system-ui;max-width:1150px;margin:30px auto;padding:0 24px;color:#22313a}a{color:#126d80}img{max-width:100%}td,th{padding:9px;border-bottom:1px solid #dde3e5;text-align:left}table{border-collapse:collapse}h3{margin-top:2.5em}input{font:inherit;padding:8px;width:350px}"
    (OUT / "report.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>Cancer-only audit</title><style>'
        + style
        + "</style>"
        + markdown.markdown(report, extensions=["tables"])
    )
    rows = []
    for r in artifact["records"]:
        rows.append(
            "<tr>"
            + "".join(
                "<td>" + v + "</td>"
                for v in [
                    f'<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={r["gene_id"]}">{r["protein"]}</a>',
                    r["classification"],
                    f"{r['developmental_max_tpm']:.4g}",
                    f"{r['tumor_tpm_max']:.4g}",
                    html.escape("Admission reasons: " + r["admission_reasons"])
                    + "<br>"
                    + "<br>".join(html.escape(v) for v in r["notes"]),
                    " · ".join(
                        f'<a href="{html.escape(url)}">{html.escape(source_title(url, artifact))}</a>'
                        for url in r["sources"]
                    ),
                ]
            )
            + "</tr>"
        )
    (OUT / "table.html").write_text(
        '<!doctype html><meta charset="utf-8"><title>47 protein reviews</title><style>'
        + style
        + '</style><h1>47 proteins with cancer RNA and zero adult reproductive bulk RNA</h1><p><a href="report.html">Methods and full report</a> · <a href="protein-reviews.csv">Download table</a></p><label>Filter <input type="search" id="filter"></label><table><thead><tr>'
        + "".join(
            f"<th>{s}</th>"
            for s in [
                "Protein",
                "Decision",
                "Developmental TPM",
                "Tumor maximum TPM",
                "Normal evidence / review",
                "Sources",
            ]
        )
        + "</tr></thead><tbody>"
        + "".join(rows)
        + '</tbody></table><script>document.getElementById("filter").addEventListener("input",e=>{const q=e.target.value.toLowerCase();document.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script>'
    )


def source_title(url, artifact):
    if "proteinatlas.org" in url:
        return "HPA normal-tissue RNA, IHC and protein-evidence record"
    return artifact["source_titles"].get(url, url)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply-curation", action="store_true")
    build(parser.parse_args().apply_curation)
