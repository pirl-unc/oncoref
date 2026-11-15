"""Fold the 2026-11-12 candidate literature sweep into the central literature review.

Reads the sweep's findings and coverage (docs/audits/target-filters-20261112) and
rewrites oncoref/data/normal-expression-literature.json, the one literature store the
catalog reads. Each swept protein gets one entry under the review's existing rubric.

An exclusion needs more than the sweep recorded: a verbatim quote retrieved from the
source. EXCLUSIONS below holds the fourteen proteins for which one was retrieved on
2026-11-13, with the quote. Every other finding is recorded with its tissue and tier
and has no category effect; RUBRIC_GAPS says why for those in a Vital or Protected
tissue. Re-running replaces the sweep's entries and leaves the earlier review as is.

Usage: PYTHONPATH=. python scripts/fold_normal_function_sweep.py
"""

import json
import re
from pathlib import Path

import pandas as pd

from oncoref import target_selection, tissue_tiers

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/audits/target-filters-20261112"
LITERATURE = ROOT / "oncoref/data/normal-expression-literature.json"
SWEEP, EARLIER = "candidate_sweep", "rare_cell_dependents"
OLENDER = (
    "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2960-3/"
    "MediaObjects/12864_2016_2960_MOESM2_ESM.xlsx"
)
OLENDER_TITLE = (
    "Olender et al. 2016, BMC Genomics 17:619: The human olfactory transcriptome "
    "(Additional file 2, Table S2, FPKM per gene)"
)


def pubmed(pmid):
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def pmc(pmcid):
    return f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"


def item(title, url, assay, tissue, finding, quote, limitations):
    return {
        "source_title": title,
        "url": url,
        "assay": assay,
        "tissue": tissue,
        "finding": finding,
        "quote": quote,
        "paralog_specific": True,
        "limitations": limitations,
        "review": SWEEP,
        "retrieved": "2026-11-13",
    }


#: Registry tissue (and cell population) of each compartment a finding was assigned to.
COMPARTMENTS = {
    "skin": ("skin", None),
    "sweat_gland": ("skin", "sweat gland cells"),
    "melanocytes": ("skin", "melanocytes"),
    "hair_follicle": ("hair follicle", None),
    "cns": ("brain", None),
    "eye": ("eye", None),
    "cartilage": ("cartilage", None),
    "stomach": ("stomach", None),
    "lung": ("lung", None),
    "pancreas": ("pancreas", None),
    "adrenal": ("adrenal gland", None),
    "kidney": ("kidney", None),
    "thyroid_c": ("thyroid gland", "parafollicular cells"),
    "adipose": ("adipose tissue", None),
    "breast": ("breast", None),
    "skeletal_muscle": ("skeletal muscle", None),
    "myeloid": ("blood", "macrophages"),
    "megakaryocyte": ("blood", "platelets"),
    "lymph_node": ("lymph node", None),
    "salivary": ("salivary gland", None),
    "oral": ("tongue", None),
    "olfactory_epithelium": ("olfactory epithelium", None),
    "taste_buds": ("taste bud", None),
    "dental_pulp": ("dental pulp", None),
    "sensory_ganglion": ("dorsal root ganglion", None),
    "nasal_mucosa": ("nasal mucosa", None),
    "lacrimal": ("lacrimal gland", None),
    "bone": ("bone", None),
}

#: R2 exclusions: (decision tissue, reason, evidence items with retrieved quotes).
EXCLUSIONS = {
    "ADAM30": (
        "brain",
        "RNA measured in human frontal cortex of non-demented controls by array and probe "
        "assay, with neuronal immunostaining. Absolute level is not given; controls are "
        "aged post-mortem brains.",
        [
            item(
                "Letronne et al. 2016, EBioMedicine 9:278: ADAM30 Downregulates APP-Linked "
                "Defects Through Cathepsin D Activation in Alzheimer's Disease",
                pmc("PMC4972530"),
                "custom oligonucleotide array and probe-based mRNA quantification; IHC",
                "human frontal cortex, 42 controls and 51 Alzheimer's disease cases",
                "ADAM30 mRNA is measured in control brains and is about half as high in cases.",
                "By transcriptomic analyses, low levels of ADAM30 expression was observed in "
                "the brain of AD cases (compared with controls).",
                "Values are relative. The antibody's validation was not checked.",
            ),
            item(
                "Same paper, immunohistochemistry",
                pmc("PMC4972530"),
                "IHC",
                "human brain",
                "ADAM30 immunostaining in neurons.",
                "immunohistochemistry experiments were consistent with the expression of human "
                "ADAM30 in neurons in human brain tissue",
                "One polyclonal antibody.",
            ),
        ],
    ),
    "ASCL4": (
        "skin",
        "Gene-specific RT-PCR restricted to skin, quantified in adult skin. Read at abstract "
        "level; the adult level is one seventh of fetal and its absolute value is not given.",
        [
            item(
                "Jonsson et al. 2004, Genomics 84:859: Hash4, a novel human achaete-scute "
                "homologue found in fetal skin (PubMed abstract)",
                pubmed(15475265),
                "RT-PCR tissue panel and real-time PCR",
                "human fetal and adult skin",
                "Expression restricted to skin; adult skin is positive at a lower level.",
                "the expression of Hash4 was restricted to skin. Quantification of the Hash4 "
                "expression by real-time PCR revealed a 7-fold higher expression in fetal skin "
                "compared to adult skin.",
                "Abstract only; source and normality of the adult skin RNA not checked.",
            )
        ],
    ),
    "BECN2": (
        "brain",
        "Protein detected in adult human brain lysate by western blot with a peptide "
        "antibody. One laboratory; no cell-type resolution.",
        [
            item(
                "He et al. 2013, Cell 154:1085: Beclin 2 functions in autophagy, degradation "
                "of G protein-coupled receptors, and metabolism",
                pmc("PMC4231430"),
                "western blot",
                "human fetal and adult brain lysates",
                "Beclin 2 protein in adult human brain.",
                "We also detected human Beclin 2 protein expression (predicted molecular weight "
                "49 kD; apparent molecular weight 53 kD) in samples from human fetal and adult "
                "brains using an antibody generated against amino acids 49-67 of human Beclin 2",
                "Single blot; the tissue RNA panel in the paper is mouse.",
            )
        ],
    ),
    "BPIFB3": (
        "olfactory epithelium",
        "Gene-level RNA-seq: 242 to 1,234 FPKM in all four human olfactory epithelium "
        "samples, 0.06 in respiratory epithelium. Donor ages are not stated.",
        [
            item(
                OLENDER_TITLE,
                OLENDER,
                "RNA-seq",
                "human olfactory epithelium, four samples, against 16 control tissues",
                "High in every olfactory sample; at the floor in the control tissues.",
                "BPIFB3: OE7 648.385, OE12 332.77, OE15 242.483, OE17 1234.29, respiratory 0.061",
                "Table row, not prose. Biopsy and autopsy samples of unstated age.",
            )
        ],
    ),
    "BPIFB6": (
        "olfactory epithelium",
        "Gene-level RNA-seq: 5 to 224 FPKM in all four human olfactory epithelium samples, "
        "0.07 in respiratory epithelium. Three of the four samples are below 10 FPKM.",
        [
            item(
                OLENDER_TITLE,
                OLENDER,
                "RNA-seq",
                "human olfactory epithelium, four samples, against 16 control tissues",
                "Detected in every olfactory sample; at the floor in the control tissues.",
                "BPIFB6: OE7 9.902, OE12 5.003, OE15 6.841, OE17 224.132, respiratory 0.066",
                "Table row, not prose. Biopsy and autopsy samples of unstated age.",
            )
        ],
    ),
    "CCN6": (
        "cartilage",
        "Protein in chondrocytes of normal adult articular cartilage by "
        "immunohistochemistry. Read at abstract level; one 2004 study.",
        [
            item(
                "Sen et al. 2004, Arthritis Rheum 50:488: WISP3-dependent regulation of type II "
                "collagen and aggrecan production in chondrocytes (PubMed abstract)",
                pubmed(14872491),
                "IHC",
                "normal adult human articular cartilage",
                "WISP3 (CCN6) protein in midzone chondrocytes.",
                "WISP3 protein is expressed in the midzone chondrocytes of normal adult "
                "articular cartilage",
                "Abstract only; antibody validation against other CCN proteins not checked.",
            )
        ],
    ),
    "DEFB106A/B": (
        "lung",
        "PCR on human tissue cDNA panels reports expression in lung, with in situ "
        "hybridisation. Read at abstract level; the lung cell type is not named.",
        [
            item(
                "Kao et al. 2003, Am J Respir Cell Mol Biol 29:71: ORFeome-based search of "
                "airway epithelial cell-specific novel human beta-defensin genes (PubMed abstract)",
                pubmed(12600824),
                "PCR on tissue cDNA panels; in situ hybridisation",
                "human lung",
                "DEFB106 expressed in lung.",
                "Two of them, DEFB106 and DEFB109, were expressed in the lung.",
                "Abstract only. DEFB106A and DEFB106B are identical copies.",
            )
        ],
    ),
    "GALR3": (
        "adrenal gland",
        "The human cloning paper reports mRNA detectable in adrenal gland and pancreas. "
        "Read at abstract level; a 1998 assay whose method is not given there.",
        [
            item(
                "Kolakowski et al. 1998, J Neurochem 71:2239: Molecular characterization and "
                "expression of cloned human galanin receptors GALR2 and GALR3 (PubMed abstract)",
                pubmed(9832121),
                "mRNA detection, method not stated in the abstract",
                "human adrenal gland and pancreas",
                "GALR3 mRNA detectable in adrenal gland and pancreas.",
                "GALR3 mRNA was highly expressed in human testis and was detectable in adrenal "
                "gland and pancreas.",
                "Detectable only; no cell type and no protein.",
            )
        ],
    ),
    "MRGPRX1": (
        "dorsal root ganglion",
        "qPCR of whole human dorsal root ganglia from four donors and single-cell RT-PCR of "
        "85 human sensory neurons. Donor ages are not stated; the neurons picked were mostly "
        "small.",
        [
            item(
                "Castro et al. 2019, JCI Insight 4:e131712: Activation of pruritogenic TGR5, "
                "MrgprA3, and MrgprC11 on colon-innervating afferents induces visceral "
                "hypersensitivity",
                pmc("PMC6824308"),
                "single-cell RT-PCR with TaqMan primers; qPCR",
                "human thoracolumbar dorsal root ganglia",
                "MRGPRX1 in 79% of picked human sensory neurons.",
                "Single-cell RT-PCR from 85 individual human DRG neurons, of predominately "
                "smaller diameter, demonstrated that 38% expressed TGR5 , 79% expressed MrgprX1",
                "Primer specificity against MRGPRX2-4 was not checked.",
            )
        ],
    ),
    "OBP2A": (
        "lung",
        "The cloning paper reports strong expression in lung, nasal structures, salivary and "
        "lachrymal glands, and separates it from the 95%-identical OBP2B. Read at abstract "
        "level. RNA-seq of four human olfactory samples and a lung sample shows it below "
        "0.1 FPKM, which the rubric does not count against a gene-specific report.",
        [
            item(
                "Lacazette et al. 2000, Hum Mol Genet 9:289: A novel human odorant-binding "
                "protein gene family resulting from genomic duplicons at 9q34 (PubMed abstract)",
                pubmed(10607840),
                "expression analysis, method not stated in the abstract",
                "human nasal structures, salivary and lachrymal glands, lung",
                "OBP2A (hOBP IIa) strongly expressed in lung and the oral sphere.",
                "hOBP (IIa) was strongly expressed in the nasal structures, salivary and "
                "lachrymal glands, and lung",
                "Abstract only; donor age and tissue normality not checked.",
            )
        ],
    ),
    "SCGB1C1": (
        "nasal mucosa",
        "Sequence-confirmed RT-PCR finds it in normal inferior turbinate mucosa of people "
        "without sinus disease. Gene-level RNA-seq also gives 322 to 3,126 FPKM in all four "
        "human olfactory epithelium samples.",
        [
            item(
                OLENDER_TITLE,
                OLENDER,
                "RNA-seq",
                "human olfactory epithelium, four samples, against 16 control tissues",
                "High in every olfactory sample; at the floor in respiratory epithelium.",
                "SCGB1C1: OE7 322.312, OE12 907.047, OE15 742.842, OE17 3126.19, respiratory 0.003",
                "Table row, not prose. Biopsy and autopsy samples of unstated age.",
            ),
            item(
                "Lu et al. 2011, Respir Res 12:28: Expression and regulation of secretoglobins "
                "in human upper airway",
                pmc("PMC3063214"),
                "quantitative RT-PCR, product confirmed by sequencing",
                "inferior turbinate mucosa of 16 patients without sinus disease or rhinitis",
                "SCGB1C1 mRNA in normal nasal mucosa.",
                "Their mRNA expression could be found in normal nasal mucosa except SCGB1D2 (LIPB)",
                "RNA only, relative quantification; the eight genes are named in the methods.",
            ),
        ],
    ),
    "SRY": (
        "brain",
        "Immunostaining in dopamine neurons of the adult male substantia nigra and absent in "
        "female. One laboratory and one antibody; no RNA.",
        [
            item(
                "Czech et al. 2012, J Neurochem 122:260: The human testis-determining factor "
                "SRY localizes in midbrain dopamine neurons (PubMed abstract)",
                pubmed(22568433),
                "IHC",
                "human male and female substantia nigra, post-mortem controls",
                "SRY protein in a subset of male nigral dopamine neurons.",
                "SRY immunoreactivity was detected in the human male, but not female substantia "
                "nigra pars compacta, within a sub-population of tyrosine hydroxylase (TH) "
                "positive neurons.",
                "One sheep polyclonal antibody; two female negative controls.",
            )
        ],
    ),
    "SULT6B1": (
        "kidney",
        "RT-PCR reports expression in human kidney and testis. Read at abstract level; one "
        "non-quantitative report, RNA only.",
        [
            item(
                "Takahashi et al. 2009, Chem Biol Interact 181:253: Identification of a novel "
                "thyroid hormone-sulfating cytosolic sulfotransferase, SULT6B1 (PubMed abstract)",
                pubmed(19505954),
                "RT-PCR",
                "human kidney",
                "SULT6B1 expressed in kidney.",
                "human SULT6B1 was specifically expressed in kidney and testis",
                "Abstract only; tissue source and donors not stated.",
            )
        ],
    ),
    "TCHHL1": (
        "skin",
        "Protein in the basal layer of normal human epidermis by immunohistochemistry, with "
        "mRNA in normal skin. Read at abstract level; one in-house antibody.",
        [
            item(
                "Yamakoshi et al. 2013, Biochem Biophys Res Commun 432:66: Trichohyalin-like 1 "
                "protein, a member of fused S100 proteins, is expressed in normal and pathologic "
                "human skin (PubMed abstract)",
                pubmed(23376073),
                "IHC; mRNA detection",
                "normal human epidermis",
                "TCHHL1 protein in basal keratinocytes.",
                "An immunohistochemical study showed that TCHHL1 proteins were expressed in the "
                "basal layer of the normal epidermis.",
                "Abstract only; donor age not stated.",
            )
        ],
    ),
}

#: Why a confirmed finding in a Vital or Protected tissue does not meet R1 to R3.
RUBRIC_GAPS = {
    "AADACL3": "Atlas annotation alone; no gene-specific assay.",
    "DEFB104A/B": "Database annotation alone; the primary paper could not be opened.",
    "DEFB118": "Database annotation alone; the cited abstract does not mention pancreas.",
    "DEFB130A/B": "Macrophages differentiated and stimulated in culture, not tissue.",
    "GALP": "The assay detects alarin, a splice product sharing only the peptide's first "
    "residues with GALP; antibody specificity is not established.",
    "LRRC52": "The abstract names testis only; the skeletal-muscle statement is in a full "
    "text that could not be retrieved for a quote.",
    "NLRP11": "Cultured, primed macrophages; not an in situ tissue measurement.",
    "OR51A7": "One healthy donor; the probe was not tested against OR51 paralogs.",
    "PRKACG": "R3 needs expression in the affected tissue; megakaryocyte expression could "
    "not be checked. One family with one missense variant.",
    "RTP2": "The abstract reports the RTP and REEP families in taste papillae; the "
    "gene-specific statement is a database annotation.",
    "SCGB1D4": "Induced expression; the resting pattern conflicts between sources.",
    "TAS2R1": "Family-level source for all 25 bitter receptors; not opened per gene.",
    "TAS2R7": "Family-level source for all 25 bitter receptors; not opened per gene.",
    "TAS2R13": "Family-level source for all 25 bitter receptors; not opened per gene.",
    "TAS2R41": "Family-level source for all 25 bitter receptors; not opened per gene.",
    "TAS2R50": "Family-level source for all 25 bitter receptors; not opened per gene.",
    "TAS2R60": "Family-level source for all 25 bitter receptors; not opened per gene.",
    "ZNRF4": "Cultured human monocytes and macrophages; expression is induced.",
}

#: Name patterns whose members share a published tissue. Recorded; no category effect.
FAMILY_RULES = [
    {
        "id": "olfactory_receptor",
        "label": "Olfactory receptors",
        "pattern": r"^OR\d+[A-Z]+\d+",
        "tissue": "olfactory epithelium",
        "filter_acts": True,
        "basis": "356 olfactory receptor loci were profiled by TaqMan array in whole "
        "olfactory mucosa from 26 adults; 273 were detected on average. Per-gene values "
        "were not reviewed for each candidate.",
        "sources": ["PMID 24800820"],
    },
    {
        "id": "bitter_taste_receptor",
        "label": "Bitter taste receptors",
        "pattern": r"^TAS2R",
        "tissue": "taste bud",
        "filter_acts": True,
        "basis": "All 25 human TAS2R genes are reported in taste receptor cells of "
        "circumvallate papillae by in situ hybridisation. Not opened per gene.",
        "sources": ["PMID 18003842"],
    },
    {
        "id": "hair_keratin",
        "label": "Hair keratins and keratin-associated proteins",
        "pattern": r"^(KRTAP|KRT3[1-9]|KRT40|KRT8[1-6])",
        "tissue": "hair follicle",
        "filter_acts": True,
        "basis": "Gene-specific in situ hybridisation and RT-PCR place KAP5, KAP10, KAP12 "
        "and KAP24 members and the hair keratins K37 and K39 in the hair follicle. Extended "
        "to the other members of these families.",
        "sources": [
            "PMID 15028290",
            "PMID 14962103",
            "PMID 15144888",
            "PMID 17235325",
            "PMID 15086535",
            "PMID 17301834",
        ],
    },
    {
        "id": "trace_amine_receptor",
        "label": "Trace-amine receptors",
        "pattern": r"^TAAR",
        "tissue": "olfactory epithelium",
        "filter_acts": False,
        "basis": "UniProt annotates an olfactory function by similarity to mouse only. No "
        "human tissue evidence was found.",
        "sources": ["UniProt Q9P1P5", "UniProt Q96RI9"],
    },
    {
        "id": "sweet_umami_taste_receptor",
        "label": "Sweet and umami taste receptors",
        "pattern": r"^TAS1R",
        "tissue": "taste bud",
        "filter_acts": False,
        "basis": "No human gene-specific source was opened.",
        "sources": [],
    },
    {
        "id": "late_cornified_envelope",
        "label": "Late cornified envelope proteins",
        "pattern": r"^LCE",
        "tissue": "skin",
        "filter_acts": False,
        "basis": "Induction in normal epidermis after barrier disruption is reported for the "
        "LCE group with a pan-LCE3 antibody. Not gene-specific.",
        "sources": ["PMID 19169253", "PMID 15854049", "PMID 26556599"],
    },
]

RUBRIC_RECORDED = (
    "A finding that a second reader confirmed as human, adult, normal-tissue and "
    "gene-specific but that does not meet R1 to R3: cultured or induced cells, database "
    "annotation alone, a family-level source, no retrievable quote, or a Manageable or "
    "Disposable tissue. Recorded with its tissue and tier. No category effect; the "
    "candidate filters can act on it."
)
SWEEP_METHOD = (
    " Candidate sweep: on 2026-11-12, 37 Claude research agents searched the literature "
    "for normal-tissue function or expression of every Target candidate of policy "
    "oncoref.target-selection.2026-11-10.1, 25 proteins each, and 36 second readers each "
    "tried to refute one batch's findings. On 2026-11-13 the abstract, open full text or "
    "supplementary table behind each proposed exclusion was retrieved from PubMed, PMC or "
    "the publisher and the quote recorded; an exclusion without a retrieved quote was not "
    "made."
)
SWEEP_LIMITS = [
    "Candidate-sweep findings other than the exclusions carry the second reader's "
    "paraphrase, not a quote",
    "A candidate-sweep null is one agent's search, not proof of absence; FSCB was not searched",
    "Olfactory, taste and hair-keratin family rules rest on family-level sources and were "
    "not reviewed gene by gene",
    "Tiers for olfactory epithelium, taste buds, nasal mucosa, lacrimal gland, sensory "
    "ganglia, bone and tooth come from the 2026-11-13 registry updates, written from "
    "abstracts",
]


def tier_of(compartments):
    """Most protective registry tier among a finding's compartments, with its tissue."""
    best = None
    for key in compartments:
        tissue, cell = COMPARTMENTS[key]
        found = (
            tissue_tiers.assess(tissue, cell, resolution="cell")
            if cell
            else tissue_tiers.assess(tissue)
        )
        tier = found.get("fallback_tier", found["tier"])
        if tier not in tissue_tiers.TIERS:
            raise ValueError(f"Unresolved tissue for a recorded finding: {tissue}")
        rank = tissue_tiers.TIERS.index(tier)
        if best is None or rank < best[0]:
            best = (rank, tissue if not cell else f"{tissue}: {cell}", tier)
    return (best[1], best[2]) if best else (None, None)


def url_of(row):
    for text in (row.best_source, row.source):
        if isinstance(text, str):
            if found := re.search(r"https?://[^\s;,)]+", text):
                return found.group(0).rstrip(".")
            if found := re.search(r"PMID:?\s*(\d+)", text):
                return pubmed(found.group(1))
    return None


def sweep_items(rows):
    return [
        {
            "source_title": r.best_source if isinstance(r.best_source, str) else r.source,
            "url": url_of(r),
            "assay": r.evidence_type,
            "tissue": r.tissue
            if not isinstance(r.cell_type, str)
            else f"{r.tissue}; {r.cell_type}",
            "finding": r.claim,
            "quote": None,
            "paralog_specific": bool(r.gene_specific),
            "limitations": r.verify_reason,
            "review": SWEEP,
            "second_reader_verdict": r.verdict,
            "human_adult_normal": bool(r.human_adult_normal),
        }
        for r in rows.itertuples()
    ]


def build():
    document = json.loads(LITERATURE.read_text())
    sweep = pd.read_csv(AUDIT / "normal-function-sweep.csv")
    coverage = pd.read_csv(AUDIT / "normal-function-sweep-coverage.csv")
    assigned = pd.read_csv(AUDIT / "normal-function-sweep-tissues.csv")
    assigned = assigned.groupby("protein_group_label").compartment.agg(list).to_dict()
    catalog = target_selection.results().set_index("protein_group_label")
    confirmed = (
        sweep.verdict.eq("confirmed")
        & sweep.human_adult_normal.eq(True)
        & sweep.gene_specific.eq(True)
    )
    confirmed_genes = set(sweep.gene[confirmed])
    if missing := set(EXCLUSIONS) - confirmed_genes:
        raise ValueError(f"An exclusion has no confirmed sweep finding: {sorted(missing)}")

    # Start from the earlier review alone, so a re-run replaces the sweep.
    earlier = [e for e in document["entries"] if e.get("review", EARLIER) == EARLIER]
    for entry in earlier:
        entry["review"] = EARLIER
        entry["evidence"] = [x for x in entry["evidence"] if x.get("review") != SWEEP]
        entry["recorded_tissue"] = entry["recorded_tier"] = None
    by_protein = {e["protein"]: e for e in earlier}

    statuses = coverage.set_index("protein_group_label").status.to_dict()
    # One sweep row describes a whole family, not a protein; the family rules hold it.
    sweep = sweep[sweep.gene.isin(catalog.index)]
    statuses.update(dict.fromkeys(sweep.gene, "finding"))
    added = []
    for protein in sorted(statuses):
        rows = sweep[sweep.gene.eq(protein)]
        recorded_tissue, recorded_tier = (
            tier_of(assigned.get(protein, [])) if protein in confirmed_genes else (None, None)
        )
        if protein in by_protein:
            # Reviewed before: keep its decision, add what the sweep saw.
            entry = by_protein[protein]
            entry["evidence"] += sweep_items(rows)
            entry["recorded_tissue"], entry["recorded_tier"] = recorded_tissue, recorded_tier
            continue
        record = catalog.loc[protein]
        entry = {
            "entity_id": record.entity_id,
            "protein": protein,
            "gene_ids": json.loads(record.member_gene_ids),
            "verdict": None,
            "reported_tissue": "; ".join(dict.fromkeys(rows.tissue.dropna())) or None,
            "summary": None,
            "effect": "none",
            "rule": None,
            "decision_tissue": None,
            "tissue_tier": None,
            "decision_reason": None,
            "evidence": sweep_items(rows),
            "toxicity_reports": None,
            "review": SWEEP,
            "recorded_tissue": recorded_tissue,
            "recorded_tier": recorded_tier,
        }
        verdicts = rows.verdict.value_counts().to_dict()
        counted = ", ".join(f"{n} {v}" for v, n in verdicts.items())
        if statuses[protein] == "no_literature":
            entry["verdict"] = "no_literature"
            entry["summary"] = "The candidate sweep found no publication to assess."
            entry["decision_reason"] = "No publication on normal adult expression was found."
        elif rows.empty:
            entry["verdict"] = "no_normal_expression_evidence"
            entry["summary"] = (
                "The candidate sweep found publications but no report of function or "
                "expression in normal adult non-reproductive tissue."
            )
            entry["decision_reason"] = (
                "No report of expression or function in normal adult non-reproductive tissue."
            )
        elif protein in EXCLUSIONS:
            tissue, reason, evidence = EXCLUSIONS[protein]
            tier = tissue_tiers.assess(tissue)["tier"]
            if tier not in {"vital", "protected"}:
                # The registry re-tiered the tissue: the finding stays, recorded, and acts no more.
                reason += (
                    f" The registry tiers {tissue} {tier}; only Vital or Protected tissues exclude."
                )
                entry.update(
                    verdict="specialized_tissue_expression",
                    decision_reason=reason,
                    recorded_tissue=tissue,
                    recorded_tier=tier,
                    summary=f"Candidate sweep: {len(rows)} finding(s) ({counted}). {reason}",
                )
                entry["evidence"] = evidence + entry["evidence"]
                added.append(entry)
                continue
            entry.update(
                verdict="vital_expression_supported"
                if tier == "vital"
                else "protected_expression_supported",
                effect="exclude",
                rule="R2",
                decision_tissue=tissue,
                tissue_tier=tier,
                decision_reason=reason,
                recorded_tissue=tissue,
                recorded_tier=tier,
                summary=f"Candidate sweep: {len(rows)} finding(s) ({counted}). {reason}",
            )
            entry["evidence"] = evidence + entry["evidence"]
        elif recorded_tier:
            gap = RUBRIC_GAPS.get(protein)
            if recorded_tier in {"vital", "protected"} and not gap:
                raise ValueError(f"{protein}: a {recorded_tier} finding needs a decision")
            entry["verdict"] = "specialized_tissue_expression"
            entry["decision_reason"] = gap or (
                f"Reported in {recorded_tissue}, which the registry tiers {recorded_tier}; "
                "only Vital or Protected tissues exclude."
            )
            entry["summary"] = (
                f"Candidate sweep: {len(rows)} finding(s) ({counted}). Confirmed in "
                f"{recorded_tissue} ({recorded_tier}). {entry['decision_reason']}"
            )
        else:
            entry["verdict"] = "conflicting_or_nonspecific"
            entry["decision_reason"] = (
                "No finding survived a second reading as human, adult, normal-tissue and "
                "gene-specific."
            )
            entry["summary"] = (
                f"Candidate sweep: {len(rows)} finding(s) ({counted}). " + entry["decision_reason"]
            )
        added.append(entry)

    families = []
    for rule in FAMILY_RULES:
        tier = tissue_tiers.assess(rule["tissue"])["tier"]
        members = sorted(p for p in catalog.index if re.search(rule["pattern"], p))
        families.append({**rule, "tier": tier, "effect": "none", "catalog_proteins": len(members)})

    entries = sorted(earlier + added, key=lambda e: e["protein"])
    if len({e["entity_id"] for e in entries}) != len(entries):
        raise ValueError("A protein has two literature entries")
    method = document["method"].split(" Candidate sweep:")[0] + SWEEP_METHOD
    limits = [x for x in document["limitations"] if x not in SWEEP_LIMITS] + SWEEP_LIMITS
    rubric = {k: v for k, v in document["rubric"].items() if k != "recorded"}
    rubric["R2"] = rubric["R2"].replace(
        "the tissue-tier registry tiers Vital or Protected,",
        "the tissue-tier registry or its updates tier Vital or Protected,",
    )
    rubric["recorded"] = RUBRIC_RECORDED
    out = {
        "schema": document["schema"],
        "version": 2,
        "review_date": "2026-11-15",
        "scope": (
            "Every protein that would be a Target candidate while its rare-cell signal is "
            "dismissed, or after release from the pooled 1% detection flag, and every Target "
            "candidate of the 2026-11-10 catalog."
        ),
        "reviews": {
            EARLIER: {"date": "2026-11-11", "entries": len(earlier)},
            SWEEP: {"date": "2026-11-12", "entries": len(added), "quotes_retrieved": "2026-11-13"},
        },
        "rubric": rubric,
        "method": method,
        "limitations": limits,
        "effects": pd.Series([e["effect"] for e in entries]).value_counts().to_dict(),
        "family_rules": families,
        "entries": entries,
    }
    LITERATURE.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    return out


if __name__ == "__main__":
    result = build()
    acting = [e["protein"] for e in result["entries"] if e["effect"] == "exclude"]
    print(len(result["entries"]), "entries;", result["effects"], "; excluded:", ", ".join(acting))
