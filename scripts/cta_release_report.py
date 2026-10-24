#!/usr/bin/env python3
"""Build the 1.19.1 release review without changing any screening decision.

The HTML documentation and PDF share the same report blocks. Immune findings
are joined by exact gene IDs; absence of a trigger is never absence of expression.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from pathlib import Path

import pandas as pd
from render_cta_atlas import WIDTH, Report, cover, link, paragraph, section, table
from reportlab.platypus import PageBreak, Spacer

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "docs/audits/cta-atlas-20261024"
CELLS = ROOT / "docs/audits/cta-normal-cells-20261023"
OUT = ROOT / "docs/audits/cta-release-1.19.1"
PDF = ROOT / "output/pdf/cta-release-review.pdf"
ATLAS_URL = "https://pirl-unc.github.io/oncoref/audits/cta-atlas-20261024/"
POLICY = "oncoref.target-sets.2026-10-24.1"

REFERENCES = {
    "atlas": ("Complete CTA evidence atlas and gene catalogs", ATLAS_URL),
    "rules": (
        "Current high-normal-RNA policy and candidate follow-up",
        "https://pirl-unc.github.io/oncoref/audits/cta-neither-review-20261024/",
    ),
    "cells": (
        "Normal-cell screen: source-specific rules and measurements",
        "https://pirl-unc.github.io/oncoref/audits/cta-normal-cells-20261023/",
    ),
    "dice": ("Schmiedel et al., DICE, Cell 2018", "https://pubmed.ncbi.nlm.nih.gov/30449622/"),
    "eosinophils": (
        "Korn et al., benralizumab treatment for up to five years, 2021",
        "https://pubmed.ncbi.nlm.nih.gov/34487870/",
    ),
    "basophils": (
        "Benralizumab randomized allergic-asthma study, 2024",
        "https://publications.ersnet.org/content/erj/64/3/2400512",
    ),
    "b_cells": (
        "CD19 CAR T therapy, hypogammaglobulinemia and infection, clinical cohort",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC11805655/",
    ),
    "treg": (
        "Kim et al., regulatory T-cell depletion and autoimmunity in mice, 2007",
        "https://www.nature.com/articles/ni1428",
    ),
    "il5": (
        "Normal human IL-5-producing T-helper-cell biology",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC3445433/",
    ),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def immune_tables(genes, dice, immune):
    """Summarize measured review triggers, never infer cell abundance or tolerance."""
    ids = set(genes["Gene ID"])
    dice = dice.loc[dice.gene_id.isin(ids)].copy()
    immune = immune.loc[immune.gene_id.isin(ids)].copy()
    assert set(dice.state) <= {"activated", "unstimulated"}
    assert dice.unit.eq("TPM").all()
    assert dice.review.eq(dice.high_measurements.ge(3)).all()
    assert not immune.loc[immune.source.eq("hpa_immune_schmiedel"), "review"].any()
    active = set(dice.loc[dice.review & dice.state.eq("activated"), "gene_id"])
    resting = set(dice.loc[dice.review & dice.state.ne("activated"), "gene_id"])
    other = set(immune.loc[immune.review, "gene_id"])
    assessed = set(dice.gene_id) | set(immune.gene_id)
    result = genes[["Gene", "Gene ID", "Category", "Detail", "Immune RNA findings"]].copy()
    result["DICE activated trigger"] = result["Gene ID"].isin(active)
    result["DICE unstimulated trigger"] = result["Gene ID"].isin(resting)
    result["Other immune trigger"] = result["Gene ID"].isin(other)
    result["Lower-threshold immune evidence assessed"] = result["Gene ID"].isin(assessed)
    result["Activated DICE triggers only"] = result["Gene ID"].isin(active - resting - other)
    result["Interpretation"] = "No lower-threshold immune trigger in assessed tables; not absence"
    result.loc[~result["Lower-threshold immune evidence assessed"], "Interpretation"] = (
        "Not assessed in the imported lower-threshold tables"
    )
    result.loc[result["Gene ID"].isin(active | resting | other), "Interpretation"] = (
        "Normal immune RNA review; not demonstrated toxicity"
    )
    result.loc[result["Activated DICE triggers only"], "Interpretation"] = (
        "Activated DICE triggers only among assessed lower-threshold tables; other gates still apply"
    )
    details = []
    for row in dice.loc[dice.review].itertuples():
        details.append(
            {
                "Gene ID": row.gene_id,
                "Source": row.source,
                "Cell type": row.cell_type,
                "State": row.state,
                "Unit": row.unit,
                "Maximum": row.maximum,
                "Values >=5": row.high_measurements,
                "Reported values": row.measurements,
                "Basis": "At least three reported values >=5 TPM; not a fraction of cells",
            }
        )
    for row in immune.loc[immune.review].itertuples():
        details.append(
            {
                "Gene ID": row.gene_id,
                "Source": row.source,
                "Cell type": row.cell_type,
                "State": row.state,
                "Unit": row.unit,
                "Maximum": row.value,
                "Values >=5": "",
                "Reported values": "",
                "Basis": "Sorted-cell summary >=5 in native units; not a donor/cell frequency",
            }
        )
    details = genes[["Gene", "Gene ID", "Category", "Detail"]].merge(
        pd.DataFrame(details), on="Gene ID", validate="one_to_many"
    )
    return result, details


def report_blocks(genes, immune, summary):
    blocks = []

    def heading(title, new_page=True):
        blocks.append({"kind": "heading", "text": title, "new_page": new_page})

    def prose(text, *refs):
        blocks.append({"kind": "paragraph", "text": text, "refs": list(refs)})

    def grid(headers, rows, proportions):
        blocks.append({"kind": "table", "headers": headers, "rows": rows, "widths": proportions})

    def named(names):
        return genes.set_index("Gene").loc[names].reset_index().to_dict("records")

    heading("Release snapshot", False)
    prose(
        "oncoref 1.19.1 documents the current CTA list and its caveats. It adds this detailed "
        "review and an immune-state audit; it changes no gene, protein, tissue category or "
        "threshold. The 399 retained genes remain the frozen screening result under policy "
        "2026-10-24.1, not a list with established clinical safety."
    )
    grid(
        ["Scope", "Retained CTA", "Neither", "Excluded"],
        [
            ["Mapped CTA nominations", "399", "561", "1,580"],
            ["Canonical human genes", "399", "1,999", "17,699"],
            ["Complete unique proteins", "723", "2,878", "126,384"],
        ],
        [0.43, 0.19, 0.19, 0.19],
    )
    prose(
        "CTA means selected by this research screen. Unsafe means membership in the forbidden "
        "background, not necessarily observed clinical toxicity. Neither supplies neither targets "
        "nor default forbidden proteins. It includes uncertainty, weak target evidence and "
        "reconsiderable findings; it is not a safe reserve list.",
        "rules",
    )
    prose(
        "The most consequential new literature finding is that 37 retained genes have normal-tissue "
        "or isoform conflicts. ABCC11 and MAGEC2 need priority neural protein review; YBX2 has "
        "vascular/adrenal staining to adjudicate. Their flags are visible beside the frozen "
        "decisions. This report does not silently clear those conflicts.",
        "atlas",
    )
    prose(
        "Read this review for rules and interpretation; use the companion atlas for every gene's "
        "full evidence. Its 138-page retained, 189-page Neither and 588-page excluded catalogs "
        "cover all 2,540 mapped nominations, with bookmarks and PDF page links.",
        "atlas",
    )
    prose(
        "Evidence labels matter: a trial registration is a nomination, human tumor recognition "
        "is evidence of recognition, and a clinical response applies to the studied product and "
        "patient selection. None automatically proves reproductive restriction or general safety."
    )

    heading("How a gene enters the retained set")
    prose(
        "Start from a mapped curated CTA nomination with usable normal-tissue coverage. Require "
        "measured reproductive expression and either reproductive bulk RNA >=5 nTPM, reliable "
        "positive reproductive protein, or cited reproductive expression plus tumor/antigen-recognition "
        "support. Good reproductive protein evidence can satisfy admission below 5 nTPM. "
        "CT55, with reproductive bulk RNA 3.7 nTPM, qualifies through published support. "
        "Tumor expression alone cannot override an exclusion.",
        "cells",
    )
    prose(
        "All normal-tissue, donor, antibody, immune/cellular and safety gates still apply. "
        "Insufficient reproductive evidence or selectivity alone leads to Neither, not automatically "
        "Unsafe. Unresolved identity and missing observations remain explicit gaps. BAGE and GAGE2B "
        "remain recorded in intake but outside the 2,540 mapped-gene denominator."
    )
    grid(
        ["Tissue category", "Examples", "Treatment in this screen"],
        [
            [
                "Reproductive",
                "Testis and male tract; ovarian/female tract; mammary; placenta, extraembryonic and explicitly fetal/embryonic compartments",
                "Required positive expression; permitted normal context",
            ],
            [
                "Disposable",
                "Mapped hair follicles, cutaneous melanocytes, tonsil, appendix, gallbladder",
                "Permitted under conditional project assumptions",
            ],
            [
                "Manageable",
                "Other skin, spleen, thymus, thyroid, adipose",
                "Still screened; positive normal IHC excludes",
            ],
            [
                "Protected",
                "Salivary glands, lymph nodes, cartilage, several mucosal and urinary compartments",
                "Stricter screened normal contexts",
            ],
            [
                "Vital",
                "Heart, lung, CNS, liver, kidney, pancreas, eye, muscle, blood/marrow and other critical compartments",
                "Strictest normal RNA tier",
            ],
        ],
        [0.17, 0.52, 0.31],
    )
    prose(
        "These examples are not the full anatomy registry. A fetal liver sample is reproductive "
        "by developmental scope; adult maternal liver is not. Immune cells inside a permitted "
        "organ retain their immune identity. The broad set is not a pregnancy- or fertility-sparing "
        "selection. Surgical removability does not establish tolerability of systemic immune injury."
    )

    heading("RNA ceilings, IHC confidence and selectivity")
    prose(
        "Accepted positive normal IHC, including Low, vetoes selection in Vital, Protected or "
        "Manageable compartments. Uncertain/cross-reactive primary staining is handled as a "
        "separate antibody review, not accepted as a definitive protein finding."
    )
    prose(
        "RNA relaxation requires complete negative IHC across the tissue's routine measured "
        "compartments plus a reliable positive reproductive control in the curated gene profile. "
        "Missing/uncertain or incomplete local staining, or no reproductive control, uses baseline "
        "limits. The export does not prove matched antibody/batch/donor across tissues or calibrate "
        "RNA-to-protein conversion."
    )
    grid(
        ["Category", "Baseline", "Approved negative", "Supported negative", "Enhanced negative"],
        [
            ["Vital", "2 / 80%", "4 / 75%", "6 / 70%", "12 / 65%"],
            ["Protected", "6 / 75%", "8 / 70%", "12 / 65%", "24 / 60%"],
            ["Manageable", "10 / 70%", "20 / 65%", "30 / 60%", "30* / 55%"],
        ],
        [0.2, 0.17, 0.21, 0.21, 0.21],
    )
    prose(
        "Cells show effective maximum nTPM / selectivity setting. *The underlying Enhanced "
        "Manageable ceiling is 60; the independent >30 rule makes the effective ceiling 30. "
        "The thresholds are project research settings, not clinically validated safety limits."
    )
    prose(
        "For each screened bulk tissue: N <= C; N <= R(1-s)/s + 0.1; and N < R. "
        "R is the highest reproductive nTPM, N the normal tissue's nTPM, C its absolute ceiling, "
        "and s its selectivity setting. All three comparisons must pass. For example, with "
        "R=5 and Vital baseline s=0.80, the selectivity ceiling is 1.35 nTPM. A normal value "
        "of 1.5 passes the absolute ceiling of 2 but fails selectivity."
    )
    prose(
        "This is a maximum-based comparison, not a global expression-breadth statistic. One "
        "high reproductive compartment can dominate R. The fixed 0.1 nTPM allowance softens "
        "near-zero ratios; it is not a fitted uncertainty interval. Broad low-level normal "
        "expression still merits scrutiny even when each individual comparison passes."
    )

    heading("Cellular evidence and the independent hard ceiling")
    prose(
        "Any observed >30 TPM or >30 nTPM outside explicitly mapped Reproductive/Disposable "
        "contexts yields Unsafe / High RNA. Exactly 30 does not trigger this gate. A single high "
        "observation suffices; no average, strong reproductive signal or negative IHC can waive it. "
        "Unresolved anatomy does not establish permission. This also overrides a prior Variable hold.",
        "rules",
    )
    grid(
        ["Source", "Lower review trigger", "Usual result without a stronger exclusion"],
        [
            [
                "HPA tissue/cell clusters",
                ">=5 nTPM and >=20 reads; normal cluster at least as high as the reproductive-cell maximum within that source. Substantial unresolved context or absent reproductive-cell coverage also prompts review.",
                "Neither / Cellular",
            ],
            [
                "Original DICE",
                ">=5 TPM in at least three reported values in any immune cell/state distribution",
                "Neither / Cellular",
            ],
            [
                "HPA/Monaco sorted immune summaries",
                ">=5 in the source's native nTPM or TPM; HPA's DICE summary is descriptive, not independent replication",
                "Neither / Cellular",
            ],
            [
                "Adult corneal atlas",
                "Mean per-cell CP10K >=1, >=5 positive cells and >=20 counts within a group; support in at least two groups",
                "Neither / Cellular",
            ],
            [
                "Primary antibodies",
                "Unresolved positive normal staining / curated-profile conflicts",
                "Neither / Antibody",
            ],
            [
                "Curated normal protein literature",
                "Accepted production by nonmalignant cells outside permitted contexts, including activation-dependent production",
                "Unsafe / Normal",
            ],
        ],
        [0.22, 0.56, 0.22],
    )
    prose(
        "Bulk RNA, sorted-cell RNA and cluster RNA are separate evidence sources. Even when "
        "labelled nTPM, their aggregation and normalization differ. TPM, nTPM and corneal CP10K "
        "are never pooled or numerically converted. The hard ceiling deliberately applies the "
        "same number separately to native TPM/nTPM, not to CP10K. Low read support, contamination "
        "and multimapping require explicit adjudication rather than silent exemptions.",
        "cells",
    )

    heading("What remains: named targets and cancer occurrence")
    prose(
        "The following are examples, not a ranking of clinical benefit. Cancer RNA shows the "
        "highest reported cohort fraction at >=5 native pTPM among cohorts with at least 20 "
        "measured samples. TCGA and validation cohorts remain separate. Counts are samples, "
        "not the fraction of tumor cells expressing a target."
    )
    panel = [
        "MAGEA4",
        "PRAME",
        "CTAG1A",
        "CTAG1B",
        "XAGE1A",
        "XAGE1B",
        "CT83",
        "CT55",
        "SAGE1",
        "CTCFL",
        "GAGE1",
        "SSX4",
        "SSX5",
    ]
    grid(
        ["Gene", "Screen result", "Example cancer RNA occurrence"],
        [
            [
                r["Gene"],
                r["Category"] + " / " + r["Detail"],
                r["Cancer RNA at >=5 pTPM"].split("; ")[0],
            ]
            for r in named(panel)
        ],
        [0.15, 0.2, 0.65],
    )
    prose(
        "RNA occurrence is not protein detection, peptide presentation, HLA eligibility or a "
        "response rate. Bulk tumor signals can originate from tumor-associated normal cells. "
        "The full tables include smaller cohorts, thresholds 0.1/1/2/5 pTPM, and separate "
        "Medium/High cancer-IHC counts. RNA and IHC subjects are unpaired; the aggregate "
        "cancer-IHC source does not supply antibody reliability.",
        "atlas",
    )

    heading("Why candidates are excluded or left in Neither")
    detail_names = {
        "Protein": "Accepted positive normal IHC",
        "RNA": "Bulk RNA above applicable ceiling",
        "High RNA": "Additional >30 normal RNA rule",
        "Normal": "Published normal-cell protein",
        "Safety": "Independent safety hold",
        "Antibody": "Unresolved antibody evidence",
        "Cellular": "Cell / immune RNA concerns",
        "Low": "Weak affirmative reproductive evidence",
        "Screening": "Other admission / curation requirements",
        "Selectivity": "Insufficient reproductive enrichment",
        "Unassessed": "Coverage / identity gap",
        "Variable": "Vital donor zero/high pattern",
    }
    grid(
        ["Result", "Principal recorded reason", "Genes"],
        [
            [r["Category"], detail_names[r["Detail"]], str(r["Genes"])]
            for r in summary["details"]
            if r["Category"] != "CTA"
        ],
        [0.18, 0.68, 0.14],
    )
    prose(
        "These are mutually exclusive final labels, not counts of every failing observation. "
        "A gene may fail several gates. Existing exclusions were preserved: a gene already "
        "excluded for normal protein need not be relabelled High RNA when another source is high. "
        "Earlier Unsafe assignments were not all re-audited across every expanded source."
    )
    grid(
        ["Examples", "Why not retained"],
        [
            ["MAGEA3", "Colon undifferentiated cluster 53.1 nTPM; >30 rule"],
            ["ROPN1", "Salivary serous cells 32.9 nTPM; near the hard ceiling but still excluded"],
            [
                "ACRBP / TEX101",
                "High normal-cell/immune RNA; source-native maxima are in the follow-up table",
            ],
            [
                "IL5 / IL13 / PAPPA",
                "Published normal activated immune or vascular protein production",
            ],
            ["MAGEA12", "Independent reviewed safety hold"],
            ["SPAG9", "Mapped candidate with positive normal protein evidence"],
            [
                "FOXR2",
                "Strong cancer biology supports interest, but uncertain primary antibodies and normal staining remain unresolved",
            ],
            ["CTAG2", "Variable: donor-level vital RNA fails despite negative curated local IHC"],
        ],
        [0.24, 0.76],
    )

    heading("Retained genes with higher normal RNA")
    retained = genes.loc[genes.Category.eq("CTA")].copy()
    numeric = pd.to_numeric(retained["Normal RNA maximum (nTPM)"], errors="coerce")
    top = retained.loc[numeric.nlargest(10).index]
    grid(
        ["Gene", "Highest non-permitted nTPM observation", "Interpretation"],
        [
            [
                r["Gene"],
                r["Highest normal nTPM context"],
                "Below the hard ceiling; not evidence of absent normal protein",
            ]
            for r in top.to_dict("records")
        ],
        [0.14, 0.56, 0.3],
    )
    prose(
        "These are cell-context observations, not comparable bulk-tissue concentrations. The "
        "lower tissue-cluster review also requires a comparison with reproductive-cell RNA, so "
        "a substantial normal signal can remain when the reproductive-cell maximum is larger. "
        "This is an explicit limitation of the current rule, not a conclusion that those cells "
        "could safely be targeted."
    )
    prose(
        f"{summary['retained_reproductive_control']} retained genes have a reliable curated reproductive "
        f"IHC positive control; {summary['retained_without_reproductive_protein']} have neither that "
        "control nor additional reliable reproductive protein recorded. The latter can qualify "
        "through RNA or cited support, but have a weaker direct protein foundation. "
        "The retained-review table exposes these flags per gene."
    )

    heading("Retained literature conflicts: priority cases")
    conflicts = retained.loc[retained["Literature conflict"].ne("")]
    prose(
        "All 37 flagged retained genes are listed in the companion table. These flags are "
        "independent of the frozen pass/fail result. A conflicting annotation requires checking "
        "the original tissue, antibody/probe, isoform, developmental stage and patient context. "
        "A shared family record is not independent evidence for every locus."
    )
    for r in named(["ABCC11", "MAGEC2", "YBX2", "PATE2", "TPTE2", "SUN5"]):
        prose(r["Gene"] + ": " + r["Literature conflict"] + " Action: " + r["Literature follow-up"])
        key = "conflict_" + r["Gene"]
        REFERENCES[key] = (r["Gene"] + " normal-tissue evidence", r["Literature conflict source"])
        prose("Source evidence.", key)
    prose(
        "Other flags include ACTL8/pancreas, FAM170A/brain, FKBP6/broad normal expression, "
        "PRAME/low normal annotations, MAGEA11/adrenal, CT47 shared-family brain attribution, "
        "SSX1/2 thyroid annotations and developmental-versus-adult compartment ambiguity. "
        "These findings are not all the same evidentiary strength and should not be collapsed "
        "into a single toxicity claim."
    )
    assert len(conflicts) == 37

    heading("Activated immune expression: what the audit actually shows")
    prose(
        "The current mapping sends sorted immune populations to blood/Vital. The additional "
        "immune RNA triggers are also shared across lineages. They do not estimate cell rarity, "
        "depletion tolerance or reversibility. This is conservative and biologically coarse. "
        "DICE's activated CD4/CD8 measurements are stimulated normal T cells; activated is a "
        "state, not a fixed disposable lineage.",
        "dice",
    )
    grid(
        ["Lower-threshold immune finding", "Nominated genes", "Meaning"],
        [
            [
                "Activated DICE trigger",
                str(summary["dice_activated"]),
                "Includes genes with other immune flags",
            ],
            [
                "Unstimulated DICE trigger",
                str(summary["dice_unstimulated"]),
                "Can overlap activated and summary flags",
            ],
            [
                "Other HPA/Monaco immune trigger",
                str(summary["other_immune"]),
                "Source-native summaries; not independent replication",
            ],
            [
                "Activated DICE triggers only",
                str(summary["activated_only"]),
                "No other lower-threshold immune trigger in assessed tables; not absence elsewhere",
            ],
        ],
        [0.36, 0.17, 0.47],
    )
    only = immune.loc[immune["Activated DICE triggers only"]]
    grid(
        ["Gene", "Current result", "Observed activated distribution(s)"],
        [
            [r["Gene"], r["Category"] + " / " + r["Detail"], r["Immune RNA findings"]]
            for r in only.to_dict("records")
        ],
        [0.16, 0.23, 0.61],
    )
    prose(
        "This seven-gene flag is not a rescue list. Three remain Unsafe / High RNA, two have "
        "antibody holds, and two are Neither / Cellular. Removing activation-related review "
        "would not automatically admit them. All other gates need a complete recomputation."
    )
    prose(
        "A fraction such as 3/88 is a fraction of reported expression values meeting a threshold, "
        "not a fraction of cells. DICE columns are anonymous per-gene distributions here and "
        "cannot be paired across genes as stable donor IDs. Rarity across people, rarity of "
        "a cell lineage and rarity of expressing cells within that lineage are different questions."
    )

    heading("Immune depletion: a caveat, not a new permission rule")
    prose(
        "Clinical experience supports distinguishing immune lineages rather than equating "
        "all blood expression with the same consequence. It does not establish that a "
        "gene expressed in a depleted lineage is reproductive-restricted. The following "
        "interpretations are recommendations for future review, not implemented exemptions."
    )
    prose(
        "Eosinophils: benralizumab follow-up in selected asthma patients found no new safety "
        "signals during treatment extending up to five years. This supports considering "
        "conditional manageability of eosinophil depletion; it does not validate permanent "
        "elimination in every tissue, population or treatment context.",
        "eosinophils",
    )
    prose(
        "Basophils: a small population with observed reduction in benralizumab studies, including "
        "blood and marrow findings. These observations support separate review rather than "
        "an unconditional Disposable label. Depletion extent, tissue compartment and duration matter.",
        "basophils",
    )
    prose(
        "B cells: depletion is an accepted therapeutic trade-off in some settings, but "
        "hypogammaglobulinemia and infection can require monitoring and immunoglobulin replacement. "
        "A clinical cohort links post-CD19-CAR-T antibody deficiency to infection; associated "
        "conditioning, prior treatment and other cytopenias complicate attribution. "
        "Conditionally manageable is more defensible than harmless or disposable.",
        "b_cells",
    )
    prose(
        "Activated CD4/CD8 cells: ordinary T cells can enter these states during immune responses. "
        "A low resting signal cannot establish that depletion during activation is tolerable. "
        "The relevant lineage and state distribution need review; the RNA assay alone does not "
        "show presentation of a target peptide or actual cell killing.",
        "dice",
    )
    prose(
        "Regulatory T cells: small numbers do not imply dispensability. Experimental depletion "
        "in neonatal and adult mice showed a critical role in preventing autoimmune pathology. "
        "That is mechanistic mouse evidence, not a numerical human safety limit.",
        "treg",
    )
    prose(
        "IL5 illustrates a modality distinction: neutralizing a secreted cytokine, depleting "
        "IL-5-receptor-bearing eosinophils, and killing IL5-producing cells are different actions. "
        "The latter includes normal activated T-helper biology; tolerating eosinophil depletion "
        "does not by itself clear an IL5-producing-cell target.",
        "il5",
    )
    prose(
        "A future revision should record lineage, activation state, evidence of normal protein "
        "or peptide, depletion extent/duration, clinical manageability and nonimmune expression "
        "separately. This release leaves blood/Vital mapping and all immune thresholds unchanged."
    )

    heading("Neither follow-up and numerical boundaries")
    for r in named(["CTAG2", "FOXR2", "HORMAD1", "TDRD1", "BRDT", "DAZL"]):
        prose(
            r["Gene"]
            + " - "
            + r["Category"]
            + " / "
            + r["Detail"]
            + ". "
            + r["Review question"]
            + " "
            + r["Reconsideration evidence"]
        )
    prose(
        "CTAG2's pinned heart RNA values are 13.5, 0, 0 and 7.1 nTPM. The applicable "
        "absolute ceiling is 12 and the selectivity-derived ceiling is 5.538; both nonzero "
        "donors fail at least one rule. Negative IHC specimens are different people. Variable "
        "requires reported zero and failing RNA in the same tissue/cohort/version/units and "
        "appropriate negative local IHC, not a claim that RNA-zero patients are treatment-eligible.",
        "rules",
    )
    prose(
        "The atlas flags 55 genes near a numerical boundary: retained maxima 24-30, high-RNA "
        "exclusions with all audited native-unit maxima <=33, or failed bulk numeric comparisons "
        "within 20% of their limits. These descriptive bands are arbitrary review aids, not new "
        "selection thresholds. ROPN1 at 32.9 is an example; it remains excluded."
    )
    prose(
        "Reconsideration should resolve the particular uncertainty: assay/paralog mapping, "
        "ambient contamination, cell identity, isoform, activation state or donor dependence. "
        "Familiarity, a trial registration or tumor abundance alone cannot waive a normal finding."
    )

    heading("Clinical efficacy, safety and biological evidence")
    prose(
        "Nineteen mapped genes received targeted primary-source reviews in the atlas; 30 have "
        "curated trial/therapy links. This is not a systematic clinical review of all 2,540 genes. "
        "The source-linked biological annotations cover 2,505 genes, with function annotations "
        "for 1,873. Missing review is not evidence of no publications, no benefit or no risk."
    )
    reviews = pd.read_csv(ATLAS / "literature-reviews.csv", keep_default_na=False)
    for identifier in [
        "tecelra-2026",
        "ima203-2025",
        "nyeso-2015",
        "lete-mrcls-2025",
        "magea3-titin-2013",
        "magea3-neuro-2013",
        "magrit-2016",
        "ct55-2017",
        "foxr2-2022",
    ]:
        if identifier == "magea3-titin-2013":
            heading("Clinical safety and biological interpretation")
        r = reviews.loc[reviews.id.eq(identifier)].iloc[0]
        REFERENCES[identifier] = (r.title, r.source)
        prose(r.title + ": " + r.finding, identifier)
    prose(
        "Mechanism matters: an engineered receptor's recognition of an unrelated normal protein "
        "does not establish native CTA expression in that organ. Conversely, a tolerated product "
        "does not establish safety for every receptor, epitope, modality, dose or population. "
        "Cytokine-release syndrome and conditioning toxicity must remain separate from established "
        "target-mediated destruction of a normal compartment."
    )

    heading("Limits of this release and downstream use")
    for text in [
        "Frozen selection versus new concern: the 37 retained literature conflicts do not automatically modify the target API or protein FASTAs. Users must inspect the separate review flags; the raw 399-gene label alone does not resolve them.",
        "Expression breadth: the current selectivity rule uses a reproductive maximum. It has no continuous whole-profile breadth penalty, no fitted RNA uncertainty model and no empirical false-negative IHC calibration.",
        "Measurement coverage: no missing observation is zero. Adult, fetal, stimulated, healthy-volunteer and tumor-adjacent contexts are not interchangeable. Earlier Unsafe genes have incomplete expanded-source audit coverage because the review focused on potential targets and Neither.",
        "Donor scope: original DICE/immune distributions are available, but the explicit reviewed vital-organ zero/high donor rule is currently populated only for CTAG2 heart. Other organs and genes are not thereby donor-negative.",
        "Gene versus protein: gene-level RNA does not allocate expression to isoforms. Exact complete sequences shared by multiple genes use Unsafe > Neither > CTA. Full protein sets are disjoint but can share k-mers; the published exclusion recipe still requires sequence-level subtraction.",
        "Species: anatomy definitions include human, dog and mouse mappings, but this release's HPA/DICE/corneal measurements and numerical limits are human-specific. Human absence cannot establish absence in another species.",
        "Clinical scope: reproductive and Disposable permissions are research assumptions. The broad list does not protect pregnancy or fertility; separate sex/developmental selections are needed. Manageable still means screened, not permitted.",
    ]:
        prose(text)
    prose(
        "Release contents: this PDF, retained-review.csv, immune-summary.csv, immune-trigger-details.csv, "
        "activated-only-review.csv, summary.json and a provenance manifest. The all-gene atlas and "
        "v1.18.0 exact-protein assets remain the authoritative frozen evidence/sequence exports. "
        "The release manifest hashes the inputs and outputs; no raw source download or threshold "
        "change is needed to reproduce this supplement."
    )
    prose(
        "Reproduce: python scripts/cta_release_report.py. Use oncoref.target_sets.results() for current gene decisions; earlier cta_screen exports preserve their dated policy."
    )
    heading("Sources and further reading")
    for key, (title, url) in REFERENCES.items():
        blocks.append({"kind": "source", "key": key, "text": title, "url": url})
    return blocks


def render(blocks):
    markdown = [
        "# CTA release 1.19.1: the list, filtering rules and caveats\n",
        "[Download this review PDF](cta-release-review.pdf) · [All supporting tables](cta-release-tables.zip) · "
        "[Search every gene](../cta-atlas-20261024/explore.html) · [Complete gene catalogs](../cta-atlas-20261024/index.md)\n",
        "[Retained-gene review flags](retained-review.csv) · [Immune audit](immune-summary.csv) · "
        "[Activated-only review](activated-only-review.csv) · [Input provenance](provenance.json) · "
        "[Output checksums](artifacts.json) · [Review record](review.md)\n",
    ]
    story = cover(
        "The CTA list and its caveats",
        "Release 1.19.1 | 24 October 2026 | 399 retained candidates; unchanged policy 2026-10-24.1",
    )
    for block in blocks:
        kind = block["kind"]
        if kind == "heading":
            markdown.append("## " + block["text"] + "\n")
            if block["new_page"]:
                story.append(PageBreak())
            story.append(section(block["text"]))
        elif kind == "paragraph":
            refs = [REFERENCES[key] for key in block["refs"]]
            markdown.append(
                block["text"]
                + (" " + " · ".join(f"[{title}]({url})" for title, url in refs) if refs else "")
                + "\n"
            )
            story.append(paragraph(block["text"]))
            if refs:
                story.append(
                    paragraph(
                        " | ".join(link(url, title) for title, url in refs), "Small", raw=True
                    )
                )
        elif kind == "table":

            def cell(value):
                return str(value).replace("|", " / ").replace("\n", " ")

            markdown += [
                "| " + " | ".join(map(cell, block["headers"])) + " |",
                "|" + "---|" * len(block["headers"]),
            ]
            markdown += ["| " + " | ".join(map(cell, row)) + " |" for row in block["rows"]]
            markdown.append("")
            story.append(
                table(block["headers"], block["rows"], [WIDTH * p for p in block["widths"]])
            )
        elif kind == "source":
            markdown.append(f"- [{block['text']}]({block['url']})")
            story.append(paragraph(link(block["url"], block["text"]), "Small", raw=True))
            story.append(Spacer(1, 4))
    (OUT / "index.md").write_text("\n".join(markdown) + "\n")
    Report(PDF, "CTA release 1.19.1: the list and caveats").build(story)
    shutil.copy2(PDF, OUT / PDF.name)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PDF.parent.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ATLAS / "all-genes.csv.gz", keep_default_na=False)
    assert genes.Policy.eq(POLICY).all()
    assert genes.Category.value_counts().to_dict() == {"Unsafe": 1580, "Neither": 561, "CTA": 399}
    assert not genes["Gene ID"].duplicated().any()
    dice = pd.read_csv(CELLS / "dice-distributions.csv.gz", keep_default_na=False)
    immune = pd.read_csv(CELLS / "immune-rna.csv.gz", keep_default_na=False)
    immune_summary, details = immune_tables(genes, dice, immune)
    immune_summary.to_csv(OUT / "immune-summary.csv", index=False)
    details.to_csv(OUT / "immune-trigger-details.csv", index=False)
    active_only = immune_summary.loc[immune_summary["Activated DICE triggers only"]]
    active_only.to_csv(OUT / "activated-only-review.csv", index=False)
    retained = genes.loc[genes.Category.eq("CTA")].copy()
    retained["Reliable reproductive protein recorded"] = (
        retained["Reproductive positive control"] | retained["Additional reproductive protein"]
    )
    retained["Normal RNA within 20 percent below hard ceiling"] = pd.to_numeric(
        retained["Normal RNA maximum (TPM)"], errors="coerce"
    ).between(24, 30) | pd.to_numeric(
        retained["Normal RNA maximum (nTPM)"], errors="coerce"
    ).between(24, 30)
    retained["Review flag is a classification change"] = False
    retained.to_csv(OUT / "retained-review.csv", index=False)
    summary = json.loads((ATLAS / "summary.json").read_text())
    summary.update(
        {
            "release": "1.19.1",
            "classification_changes": 0,
            "dice_activated": int(immune_summary["DICE activated trigger"].sum()),
            "dice_unstimulated": int(immune_summary["DICE unstimulated trigger"].sum()),
            "other_immune": int(immune_summary["Other immune trigger"].sum()),
            "activated_only": len(active_only),
            "retained_reproductive_control": int(retained["Reproductive positive control"].sum()),
            "retained_without_reproductive_protein": int(
                (~retained["Reliable reproductive protein recorded"]).sum()
            ),
        }
    )
    assert len(active_only) == 7
    assert active_only.groupby(["Category", "Detail"]).size().to_dict() == {
        ("Unsafe", "High RNA"): 3,
        ("Neither", "Antibody"): 2,
        ("Neither", "Cellular"): 2,
    }
    summary_path = OUT / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    render(report_blocks(genes, immune_summary, summary))
    inputs = [
        ATLAS / "all-genes.csv.gz",
        ATLAS / "summary.json",
        ATLAS / "literature-reviews.csv",
        CELLS / "dice-distributions.csv.gz",
        CELLS / "immune-rna.csv.gz",
        Path(__file__).resolve(),
        ROOT / "scripts/render_cta_atlas.py",
    ]
    provenance = {
        "release": "1.19.1",
        "policy": POLICY,
        "classification_changes": 0,
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "references": REFERENCES,
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with zipfile.ZipFile(OUT / "cta-release-tables.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.iterdir()):
            if path.suffix in {".csv", ".json"} and path.name != "artifacts.json":
                archive.write(path, path.name)
    artifacts = {
        p.name: {"bytes": p.stat().st_size, "sha256": digest(p)}
        for p in sorted(OUT.iterdir())
        if p.is_file() and p.name != "artifacts.json"
    }
    (OUT / "artifacts.json").write_text(json.dumps(artifacts, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in summary.items() if k not in {"details", "coverage"}}, indent=2)
    )


if __name__ == "__main__":
    main()
