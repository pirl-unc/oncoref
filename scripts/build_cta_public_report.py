#!/usr/bin/env python3
"""Generate one textual, visual and machine-readable account of the current CTA policy."""

from __future__ import annotations

import html
import json
import re
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_cta_contract_report import digest, write_table  # noqa: E402

from oncoref import (  # noqa: E402
    cta_selection,
    cta_selection_audit,
    gene_knowledge,
    sequence_bundle,
    tissue_hazards,
    translation_qualification,
)
from oncoref.version import __version__  # noqa: E402

OUT = ROOT / "docs/audits/cta-definition-20261104"
PRESENTATION_REVISION = "2026-11-05.3"
INK = "#173B55"
TEAL = "#087E8B"
AMBER = "#A76708"
RUST = "#AC493B"
MUTED = "#536675"
CASES = ["PRAME", "MAGEC2", "LGALS13", "PLCZ1", "MAGEA12", "SERPINA2", "MROH5", "SELENOV"]


def figure(name, fig):
    for suffix in ["svg", "png"]:
        fig.savefig(
            OUT / "figures" / f"{name}.{suffix}",
            dpi=170,
            bbox_inches="tight",
            metadata={"Date": None} if suffix == "svg" else {},
        )
    plt.close(fig)


def figures(facts, roles, coverage):
    (OUT / "figures").mkdir(exist_ok=True)
    # The walkthrough renderer may run in the same process with its own style.
    plt.rcdefaults()
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "svg.fonttype": "none",
            "svg.hashsalt": "oncoref-cta-definition",
        }
    )
    counts = facts.category.value_counts()
    u = coverage["universe_counts"]
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis("off")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    boxes = [
        (0.02, 0.59, 0.25, 0.26, f"{u['input_groups']:,}\nCandidate\nprotein records", INK),
        (0.37, 0.59, 0.25, 0.26, f"{len(facts):,}\nClassified proteins", TEAL),
    ]
    for i, (category, color) in enumerate([("CTA", TEAL), ("Other", AMBER), ("Unsafe", RUST)]):
        boxes.append((0.75, 0.73 - i * 0.28, 0.23, 0.2, f"{counts[category]:,}\n{category}", color))
    for x, y, w, h, t, c in boxes:
        ax.add_patch(
            FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01", facecolor=c, edgecolor="none")
        )
        ax.text(
            x + w / 2,
            y + h / 2,
            t,
            color="white",
            ha="center",
            va="center",
            fontsize=17,
            fontweight="bold",
        )
    for a, b in [
        ((0.28, 0.72), (0.35, 0.72)),
        ((0.63, 0.72), (0.73, 0.83)),
        ((0.63, 0.69), (0.73, 0.55)),
        ((0.63, 0.66), (0.73, 0.27)),
    ]:
        ax.annotate("", xy=b, xytext=a, arrowprops={"arrowstyle": "->", "color": MUTED, "lw": 2})
    ax.set_title(
        "Every candidate record receives a classification",
        loc="left",
        fontsize=22,
        fontweight="bold",
        color=INK,
        pad=20,
    )
    restored = facts.loc[~facts.translation_support_eligible, "category"].value_counts()
    ax.text(
        0.02,
        0.32,
        f"Includes all {u['translation_support_insufficient_groups']:,} records with insufficient translation support:\n"
        f"{restored.get('Other', 0):,} Other:Translation and {restored.get('Unsafe', 0):,} Unsafe.\n"
        "Translation controls CTA admission; normal exclusions take precedence.",
        color=INK,
        fontsize=12,
        linespacing=1.7,
    )
    ax.text(
        0.02,
        0.02,
        "Identical selected reference sequences count together. Member loci remain separately traceable.",
        color=MUTED,
        fontsize=11,
    )
    figure("universe-accounting", fig)
    funnel = pd.read_csv(
        ROOT / "docs/audits/cta-source-walkthrough-20261028/figures/category-funnel.csv"
    )
    fig, ax = plt.subplots(figsize=(12, 6))
    labels = funnel.stage.str.replace(
        "Reproductive/developmental presence", "Reproductive/developmental\npresence", regex=False
    )
    ax.barh(labels[::-1], funnel.eligible[::-1], color=TEAL)
    for i, value in enumerate(funnel.eligible[::-1]):
        ax.text(value + 10, i, f"{value:,}", va="center")
    ax.set_xlim(0, max(funnel.eligible) * 1.16)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel("Proteins still eligible after this presentation stage")
    ax.set_title(
        "Admission gates: every applicable requirement must pass",
        loc="left",
        fontsize=18,
        color=INK,
        pad=18,
    )
    fig.tight_layout()
    figure("admission-flow", fig)
    for suffix in ["png", "svg"]:
        shutil.copyfile(
            ROOT
            / "docs/audits/cta-source-walkthrough-20261028/figures"
            / f"publication-funnel.{suffix}",
            OUT / "figures" / f"publication-funnel.{suffix}",
        )
    holds = facts[[c for c in facts if c.startswith("hold:")]].astype(bool)
    names = holds.sum().nlargest(9).index.tolist()
    combinations = holds[names].apply(lambda r: tuple(r.astype(int)), axis=1).value_counts().head(9)
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(12, 8), gridspec_kw={"height_ratios": [2, 3]}, sharex=True
    )
    top.bar(range(len(combinations)), combinations.values, color=INK)
    for i, n in enumerate(combinations.values):
        top.text(i, n, f"{n:,}", ha="center", va="bottom", fontsize=10)
    top.set_ylabel("Proteins")
    top.set_title(
        "Overlapping reasons survive the funnel", loc="left", fontsize=21, color=INK, pad=18
    )
    matrix = np.array(combinations.index.tolist()).T
    bottom.imshow(
        matrix,
        cmap=matplotlib.colors.ListedColormap(["#E8EDF1", RUST]),
        vmin=0,
        vmax=1,
        aspect="auto",
    )
    bottom.set_yticks(range(len(names)), [n.removeprefix("hold:") for n in names])
    bottom.set_xticks(range(len(combinations)), [str(i + 1) for i in range(len(combinations))])
    bottom.set_xlabel("Most common combinations across the nine largest exclusion/hold layers")
    bottom.set_title(
        "Red: recorded hold or exclusion; gray: no recorded trigger in this layer",
        fontsize=10,
        pad=10,
    )
    fig.tight_layout()
    figure("independent-evidence", fig)
    table = (
        pd.crosstab(roles.category, roles.shared_sequence)
        .reindex(["CTA", "Other", "Unsafe"], fill_value=0)
        .reindex(columns=["allow", "review", "exclude"], fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.imshow(np.log1p(table), cmap="Blues")
    ax.set_xticks(range(3), ["Allow under\npolicy", "Review\nrequired", "Exclude"])
    ax.set_yticks(range(3), table.index)
    for i in range(3):
        for j in range(3):
            ax.text(
                j,
                i,
                f"{table.iloc[i, j]:,}",
                ha="center",
                va="center",
                fontsize=19,
                color="white" if table.iloc[i, j] > 1000 else INK,
            )
    ax.set_title(
        "Classification and sequence role answer different questions",
        fontsize=17,
        color=INK,
        pad=18,
    )
    fig.tight_layout()
    figure("category-role-matrix", fig)
    write_table(funnel, OUT / "admission-flow.csv")
    write_table(table.reset_index(), OUT / "category-role-matrix.csv")
    write_table(
        pd.DataFrame(
            [
                {
                    "combination": i + 1,
                    "groups": int(n),
                    **{name.removeprefix("hold:"): bool(flag) for name, flag in zip(names, key)},
                }
                for i, (key, n) in enumerate(combinations.items())
            ]
        ),
        OUT / "independent-evidence.csv",
    )


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    inventory = json.loads((ROOT / "docs/audits/cta-inventory-20261105/summary.json").read_text())
    zero = inventory["stored_zero_rna_screen"]
    preliminary = inventory["preliminary_expression_profile"]
    bundle_manifest = json.loads(
        (ROOT / "docs/audits/sequence-bundle-20261104/manifest.json").read_text()
    )
    if (
        bundle_manifest["bundle_id"] != sequence_bundle.descriptor()["bundle_id"]
        or sequence_bundle.content_id(bundle_manifest) != bundle_manifest["bundle_id"]
    ):
        raise ValueError("Public report requires the current verified sequence bundle")
    for name, expected in bundle_manifest["package_dependencies"].items():
        if digest(ROOT / "oncoref/data" / name) != expected:
            raise ValueError("Rebuild the sequence bundle before the public report: " + name)
    facts = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/protein-group-classifications.csv.gz",
        keep_default_na=False,
        low_memory=False,
    )
    roles = cta_selection.sequence_dispositions()
    coverage = cta_selection.coverage()
    outside_path = ROOT / "docs/audits/cta-selection-20261028/out-of-universe.csv.gz"
    outside = pd.read_csv(outside_path)
    outside_counts = outside.reason.value_counts().to_dict()
    support_counts = (
        facts.loc[~facts.translation_support_eligible, "translation_support_reasons"]
        .value_counts()
        .to_dict()
    )
    if len(outside) != coverage["universe_counts"]["out_of_universe_groups"]:
        raise ValueError("Out-of-universe records disagree with current coverage")
    contract = cta_selection.contract()
    counts = facts.category.value_counts().to_dict()
    publication_funnel = pd.read_csv(
        ROOT / "docs/audits/cta-source-walkthrough-20261028/figures/publication-funnel.csv"
    ).fillna("")
    publication_union = publication_funnel.set_index("source_id").loc["registered_union"]
    unannotated = publication_funnel.set_index("source_id").loc["unannotated"]
    admission = pd.read_csv(
        ROOT / "docs/audits/cta-source-walkthrough-20261028/figures/category-funnel.csv"
    ).set_index("stage")
    levels = tissue_hazards.assessments().level.value_counts().to_dict()
    migrated_count = len(gene_knowledge.historical_dispositions())
    figures(facts, roles, coverage)
    criteria = cta_selection_audit.audit().criteria
    write_table(criteria, OUT / "rule-catalog.csv")
    selected = roles.set_index("label")
    profile = facts.set_index("protein_group_label")
    knowledge = gene_knowledge.findings()
    worked = []
    for gene in CASES:
        row = selected.loc[gene]
        fact = profile.loc[gene]
        members = set(row.gene_ids)
        findings = knowledge[
            knowledge.gene_ids.map(lambda ids, members=members: bool(set(ids) & members))
        ]
        worked.append(
            {
                "label": gene,
                "entity_id": row.entity_id,
                "classification": row.classification,
                "assessment_status": row.assessment_status,
                "sequence_role": row.shared_sequence,
                "reasons": list(row.reasons),
                "translation_qualified": bool(fact.translation_qualified),
                "translation_reasons": fact.translation_qualification_reasons,
                "findings": findings[
                    [
                        "id",
                        "effect",
                        "interpretation_status",
                        "finding",
                        "decision_reason",
                        "source_title",
                        "source_url",
                        "limitations",
                        "review_depth",
                    ]
                ].to_dict("records"),
            }
        )
    write_table(pd.DataFrame(worked), OUT / "worked-decisions.csv.gz")
    sections = [
        (
            "What OncoRef calls a CTA",
            [
                contract["definition"],
                f"The current human policy produces {counts['CTA']:,} CTA, {counts['Unsafe']:,} Unsafe and {counts['Other']:,} Other proteins. It evaluates {len(facts):,} candidate proteins in the complete inventory, counting identical selected reference sequences together and retaining every member locus. This is not a complete inventory of splice isoforms or modified proteoforms.",
                "CTA is a research candidate designation; it does not establish peptide-HLA presentation, receptor recognition, therapeutic abundance or clinical safety. Reproductive/developmental support includes adult reproductive tissues, placenta, extraembryonic tissues, embryo and fetus. These routes remain distinguishable; developmental association does not imply adult testis restriction.",
            ],
            "universe-accounting",
        ),
        (
            "From candidates to CTAs",
            [
                f"All {coverage['universe_counts']['input_groups']:,} candidate records enter the same normal-evidence and admission path. The {sum(support_counts.values()):,} records formerly removed at entry stay visible: {support_counts.get('nmd_or_cds_not_defined_only', 0):,} have only NMD/CDS-not-defined reference transcripts and {support_counts.get('single_study_ribo_seq', 0):,} have only single-study Ribo-seq support. Insufficient translation support does not establish absence or normal-tissue safety.",
                "CTA admission requires qualifying translation support and the separately qualified exact reference record. Applicable normal exclusion takes precedence and establishes Unsafe; otherwise insufficient translation qualification yields Other:Translation and sequence review. An annotated coding sequence is not experimental proof of protein production.",
                "Apply all normal-tissue evidence independently. Accepted exclusion evidence establishes Unsafe. Unresolved findings or missing required coverage prevent admission and require review; they do not establish Unsafe. Every target must also pass reproductive/developmental association, cancer presence, concentration, translation qualification and sequence-containment rules.",
                f"The chart uses cumulative counts. Of {admission.loc['Normal safety + coverage', 'eligible']:,} proteins passing normal assessment, cancer presence retains {admission.loc['Cancer presence', 'eligible']:,}; reproductive/developmental presence then narrows these to {admission.loc['Reproductive/developmental presence', 'eligible']:,}. The previous order showed 645 at both gates because every reproductive/developmental survivor already met cancer presence. One compatible positive tumor RNA sample or reviewed tumor protein makes cancer presence a broad gate, not an abundance or prevalence requirement.",
                "Every applicable evidence layer is still evaluated, including after an early exclusion. Unsafe takes precedence. Display order changes neither final classifications nor the priority used to assign Other subtypes.",
            ],
            "admission-flow",
        ),
        (
            "Published nominations through the funnel",
            [
                f"All {len(publication_funnel[publication_funnel.source_kind.eq('paper')])} registered nomination papers appear in the companion funnel. Their union contains {publication_union['Inventory']:,} proteins: {publication_union['Unsafe sequence containment']:,} finish as CTA, {publication_union['Final Unsafe']:,} as Unsafe and {publication_union['Final Other']:,} as Other. A protein can occur in several papers, so individual paper rows must not be summed.",
                f"The {unannotated['Inventory']:,} candidates without a registered-paper annotation undergo the same assessment and contribute {unannotated['Unsafe sequence containment']:,} CTAs. Prior nomination is recorded provenance, not a required admission gate. Absence from this finite registry does not establish novelty or absence from the literature.",
                "A publication's measured evidence can support reproductive/developmental association, tumor expression, translation or normal-tissue counterevidence when it meets the corresponding evidence rule. A list label alone does not establish those facts. Germ-cell, testis-biased, placental and tumor-vascular lists retain their original scope alongside CTA nominations.",
            ],
            "publication-funnel",
        ),
        (
            "Normal-tissue consequences and expression evidence",
            [
                "Vital, Protected, Manageable and Disposable describe the consequence and tolerability of injury to a specified tissue/cell compartment. They are qualitative research priorities, not measured fatality probabilities. Immediate function, regenerative reserve, persistent depletion, anatomy and inflammatory amplification all matter.",
                f"The central registry distinguishes {levels['bulk']} bulk tissues and {levels['cell']} cell populations, with source-linked rationale. A cell-level exception requires resolved tissue and cell identity. An unresolved or mixed label cannot earn a narrower allowance; reproductive/developmental membership is a separate anatomical axis.",
                "RNA detection, accessible surface protein, peptide-HLA presentation, functional recognition and demonstrated injury are distinct evidence levels. CAR and TCR target availability are different questions. A weak initiating target can still produce inflammatory injury in other cells. Missing evidence is never measured zero.",
            ],
            None,
        ),
        (
            "Evidence has multiple independent meanings",
            [
                "Review completion, interpretation and policy effect are separate. A completed review can retain an unresolved interpretation. A rejected protein claim remains visible with its rationale; rejecting that claim does not prove normal-tissue absence.",
                f"All {migrated_count} historical findings have linked current dispositions. MAGEC2, LGALS13 and PLCZ1 retain material holds. Their evidence access depth, reagent identity, adult/fetal and isoform limits are stated explicitly. The broader literature review remains incomplete.",
                "Fatal neurological injury in the MAGE-A3-reactive TCR cohort remains product-specific evidence. The MAGEA12 and EPS8L2 explanations are both retained as unresolved. MAGEA12 remains Unsafe under an independent rare-cell finding. Neither receptor causality nor favorable clinical credit transfers across a gene family.",
            ],
            "independent-evidence",
        ),
        (
            "Translation, identity and downstream sequence use",
            [
                "Translation qualification asks whether the exact selected primary-locus record supports target sequence use. Complete coding annotation, start/stop codons, biotype and sequence integrity are assessed separately from normal-tissue hazard. An alternate-locus protein cannot silently rescue an incomplete primary record.",
                "The verified bundle retains original sequences, every source occurrence, qualification reasons and coverage gaps. It supplies contiguous AA20 segments with original coordinates and never joins across ambiguity, internal stops or different records. Selenocysteine-containing records require review under this AA20 contract; this is not a claim that the protein is unexpressed.",
                "Across occurrences, exclude takes precedence over review, then allow. Unsafe sequence remains exclusion evidence even in an incomplete record. Allow means no exclusion from the assessed record under this policy. Exact matching does not validate an approximate similarity cutoff or patient-HLA decision. The historical Unsafe plus additional-forbidden assets retain their older policy identity.",
            ],
            "category-role-matrix",
        ),
        (
            "Research parameters and coverage limits",
            [
                "The nominal RNA bases are 4, 8 and 16 TPM/nTPM for Vital, Protected and Manageable contexts. Matched negative IHC can raise a ceiling by 0-3 doubling steps; reviewed product evidence adds at most one step within that cap. Corneal CP10K keeps its separately declared scale. These are policy choices, not calibrated clinical toxicity limits.",
                "One compatible tumor sample above zero establishes the current cancer-presence criterion, not sufficient target abundance or cohort prevalence. Reproductive RNA uses a 5-unit criterion with assay-specific support requirements; reliable reproductive IHC and reviewed protein provide separate routes. Adult/placental diffuseness above 10% can block identity, never establish Unsafe.",
                "The rare-cell panel uses a 1% detection screen with cell-support rules. Detection frequency is not abundance. Coverage is limited by sampled organs, donor age, activation states, paralog mapping and assay resolution. Optional panels remain disclosed when missing. A 95% Unsafe 9-mer overlap blocks practical candidate admission, but is not a clinical cross-reactivity validation.",
                "The complete gate table below retains the actual executable rules and source IDs. Human thresholds are not transferred to dog or mouse; their native nominations remain unassessed until supported by a species-specific policy.",
            ],
            None,
        ),
        (
            "Does CTA mean expressed only in cancer?",
            [
                "No. CTA admission requires cancer expression plus qualifying reproductive or developmental association. The CTA total is not a count of proteins absent from all normal tissues, including testis, placenta, embryo and fetus. Passing normal-tissue thresholds is not a claim of zero normal expression.",
                f"The current data do support an observed-expression count: {zero['adult_and_reproductive_bulk_zero']} cancer-expressed records have complete required bulk coverage and zero stored adult/reproductive RNA summaries. Of these, {zero['also_developmental_measured']} have developmental measurements and {zero['also_developmental_zero']} also have zero developmental RNA; {zero['developmental_measurement_missing']} lack developmental measurements. Stored zeros are not proof of absent protein, and tiny positive RNA values are not validated protein production.",
                f"The separate low-expression screen finds {preliminary['matches']} preliminary matches. Independent normal-cell evidence excludes {preliminary['disposition_counts'].get('excluded_by_normal_evidence', 0)} of them. A numerical profile match is not a final cancer-only designation; remaining matches retain explicit review or coverage requirements. The inventory audit lists every record, counterevidence source and measurement gap.",
                "Cancer-specific altered sequences and peptides require a separate analysis from these reference-protein counts. A mutation, fusion or alternative transcript can create a target sequence even when the source gene is normally expressed; evidence for one such sequence cannot establish exclusivity of every product of the gene.",
            ],
            None,
        ),
        (
            "Reproduction and version changes",
            [
                f"Selection policy: {coverage['policy_id']}. Role policy: {contract['role_policy_id']}. Translation policy: {translation_qualification.POLICY_ID}.",
                "The 1.29.0 comparison remains pinned. The evidence and translation PRs also retain adjacent baselines, separating changes in category, subtype, role, assessment status and reasons. No target CTA count is an acceptance criterion.",
                "Use the public decision, finding and sequence APIs. The final release receipt binds the decision exports, ontology, curation, qualification, bundle and generated figures. Complete and cached scientific replays must be byte-identical. Historical reports keep their original inputs and policies.",
                f"Presentation revision {PRESENTATION_REVISION} places the broader cancer-presence gate first after normal assessment and follows every registered nomination paper through the same cumulative gates. Classifications, thresholds and sequence roles are unchanged. Earlier release assets remain frozen.",
            ],
            None,
        ),
    ]
    data = {
        "release": __version__,
        "presentation_revision": PRESENTATION_REVISION,
        "publication_funnel": publication_funnel.to_dict("records"),
        "out_of_universe_reasons": outside_counts,
        "translation_support_reasons": support_counts,
        "contract": contract,
        "coverage": coverage,
        "counts": counts,
        "inventory_audit": inventory,
        "sections": sections,
        "worked": worked,
        "bundle": sequence_bundle.descriptor(),
        "criteria": criteria.to_dict("records"),
        "registry_references": len(tissue_hazards.references()),
        "tier_definitions": tissue_hazards.registry()["tier_definitions"],
        "rows": roles[
            [
                "entity_id",
                "label",
                "classification",
                "assessment_status",
                "shared_sequence",
                "gene_ids",
                "reasons",
            ]
        ].to_dict("records"),
    }
    (OUT / "report-data.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    render_html(data)
    render_pdf(data)
    update_current_docs(data)
    inputs = [
        Path(__file__),
        ROOT / "oncoref/cta_contract.py",
        ROOT / "oncoref/cta_selection_policy.py",
        ROOT / "oncoref/data/cta-selection.json",
        ROOT / "oncoref/data/cta-selection-human.csv.gz",
        ROOT / "oncoref/data/sequence-bundle.json",
        ROOT / "oncoref/data/translation-qualification.json",
        ROOT / "oncoref/data/gene-knowledge.json",
        ROOT / "oncoref/data/tcr-trial-evidence.json",
        ROOT / "oncoref/data/normal-tissue-hazards.yaml",
        OUT / "manifest.json",
        ROOT / "docs/audits/cta-selection-20261028/manifest.json",
        ROOT / "docs/audits/cta-selection-20261028/reproducibility.json",
        ROOT / "docs/audits/cta-inventory-20261105/manifest.json",
        outside_path,
        ROOT / "docs/audits/cta-source-walkthrough-20261028/manifest.json",
        ROOT / "docs/audits/cta-evidence-migration-20261104/manifest.json",
        ROOT / "docs/audits/translation-qualification-20261104/manifest.json",
        ROOT / "docs/audits/normal-breadth-20261104/manifest.json",
        ROOT / "docs/audits/sequence-bundle-20261104/manifest.json",
        ROOT / "docs/audits/sequence-bundle-20261104/historical-transition.json",
    ]
    receipt = {
        "release": __version__,
        "presentation_revision": PRESENTATION_REVISION,
        "selection_policy_id": coverage["policy_id"],
        "bundle_id": data["bundle"]["bundle_id"],
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "outputs": {
            str(p.relative_to(OUT)): digest(p)
            for p in sorted(OUT.rglob("*"))
            if p.is_file() and p.name not in ["manifest.json", "release-receipt.json"]
        },
        "generated_documentation": {
            name: digest(ROOT / name)
            for name in ["README.md", "docs/cta-selection.md", "docs/index.md"]
        },
    }
    (OUT / "release-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps(
            {"counts": counts, "release": __version__, "bundle_id": data["bundle"]["bundle_id"]}
        )
    )


def update_current_docs(data):
    """Keep the two public entry points tied to the same current receipt."""
    counts = data["counts"]
    universe = data["coverage"]["universe_counts"]
    policy = data["coverage"]["policy_id"]
    summary = (
        f"**{counts['CTA']:,} CTA / {counts['Unsafe']:,} Unsafe / {counts['Other']:,} Other** "
        f"across **{universe['classified_groups']:,} proteins** "
        "(identical selected reference sequences counted together; member loci retained). "
        f"All {universe['input_groups']:,} candidate protein records enter classification; "
        "translation support is a CTA admission requirement. "
        "Every classified human protein has a category, with independent evidence status "
        "and sequence role. CTA denotes a research candidate under the declared policy."
    )
    path = ROOT / "README.md"
    text = path.read_text()
    start = text.index("The [current CTA/oncofetal selection]")
    end = text.index("\n\nThe [normal-tissue hazard ontology]", start)
    text = (
        text[:start]
        + (
            "The [current CTA/oncofetal selection](docs/audits/cta-definition-20261104/index.md) supplies "
            + summary
            + "\n[Illustrated report](https://pirl-unc.github.io/oncoref/audits/cta-definition-20261104/report.html) "
            + "· [PDF](docs/audits/cta-definition-20261104/cta-definition.pdf)."
        )
        + text[end:]
    )
    path.write_text(text)
    path = ROOT / "docs/cta-selection.md"
    text = path.read_text()
    start = text.index("The current policy supplies")
    end = text.index("\n\n##", start)
    text = (
        text[:start]
        + (
            "The current policy supplies "
            + summary
            + f" Selection policy: `{policy}`. "
            + "[Illustrated, searchable report](audits/cta-definition-20261104/report.html), "
            + "[PDF](audits/cta-definition-20261104/cta-definition.pdf), "
            + "[source walkthrough](audits/cta-source-walkthrough-20261028/index.md) and "
            + "[verified sequences](sequence-bundle.md) describe the same versioned decisions."
        )
        + text[end:]
    )
    text = re.sub(r"`oncoref\.cta-selection\.2026-\d\d-\d\d\.\d+`", f"`{policy}`", text)
    path.write_text(text)
    path = ROOT / "docs/index.md"
    text = path.read_text()
    start = text.index("[Normal-tissue hazard ontology]")
    end = text.index("The [343-group harmonized funnel]", start)
    text = (
        text[:start]
        + (
            "[Normal-tissue hazard ontology](tissue-hazards.md) centralizes cited bulk and cell-subset calls.\n\n"
            + "The [current CTA and oncofetal definition](audits/cta-definition-20261104/index.md) supplies "
            + summary
            + "\n\n[Illustrated report](audits/cta-definition-20261104/report.html) · "
            + "[PDF](audits/cta-definition-20261104/cta-definition.pdf) · "
            + "[Verified sequence records](sequence-bundle.md).\n\n"
            + "The normal screen, cancer presence, reproductive/developmental association, "
            + "coverage, translation and sequence-overlap gates remain distinct. Accepted "
            + "exclusion evidence establishes Unsafe; unresolved evidence requires review. "
            + f"[Current release notes](cta-release-{__version__}.md) and the "
            + "[source walkthrough](audits/cta-source-walkthrough-20261028/index.md) explain the decisions.\n\n"
        )
        + text[end:]
    )
    text = text.replace(
        "[Continuous normal-expression breadth](audits/normal-breadth-20261026/index.md)",
        "[Current normal-expression breadth](audits/normal-breadth-20261104/index.md)",
    )
    text = text.replace(
        "[Completed antibody adjudications](audits/cta-adjudication-20261025/index.md) cover all 178 remaining discrepancy cases. Older reports below retain their dated policies.",
        "[Historical antibody reviews](audits/cta-adjudication-20261025/index.md) retain their original findings. The [current migration ledger](audits/cta-evidence-migration-20261104/index.md) explains all revised dispositions without claiming experimental resolution of every reagent conflict. Older reports retain their dated policies.",
    )
    path.write_text(text)


def render_html(data):
    esc = html.escape
    chunks = [
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>OncoRef CTA definition</title><style>body{margin:0;background:#f4f6f7;color:#173b55;font:17px/1.6 system-ui,sans-serif}main{max-width:1120px;margin:auto;padding:50px 30px}header{border-top:8px solid #087e8b;padding-top:24px}h1{font-size:46px;line-height:1.1}h2{font-size:28px;margin-top:48px}p{max-width:83ch}figure{margin:30px 0;background:white;padding:16px}img{width:100%;height:auto}small,.muted{color:#536675}table{border-collapse:collapse;width:100%;background:white;font-size:14px}td,th{padding:12px;text-align:left;border-bottom:1px solid #d9e2e8;vertical-align:top}th{background:#173b55;color:white}a{color:#087e8b}input,select{font:inherit;padding:10px;border:1px solid #8995a4;border-radius:4px}details{background:white;padding:18px;margin:16px 0;border-left:4px solid #a76708}summary{font-weight:700;cursor:pointer}.counts{display:flex;gap:30px;flex-wrap:wrap;font-size:24px}.counts strong{font-size:36px;display:block}code{overflow-wrap:anywhere;font-size:13px}.scroller{overflow:auto}footer{margin:60px 0 10px;border-top:1px solid #8995a4;padding-top:20px}</style><main><header><small>ONCOREF / RESEARCH DEFINITION / '
        + esc(data["release"])
        + '</small><h1>From the candidate universe<br>to CTAs</h1><p>Traceable evidence. Explicit uncertainty. Qualified sequence records.</p><div class="counts">'
    ]
    for category in ["CTA", "Unsafe", "Other"]:
        chunks.append(
            f"<div><strong>{data['counts'][category]:,}</strong>{category} proteins</div>"
        )
    chunks.append(
        '</div><p><a href="cta-definition.pdf">Download PDF</a> · <a href="release-receipt.json">Report receipt</a> · <a href="report-data.json">Report data</a> · <a href="../cta-inventory-20261105/">Full inventory and expression-profile audit</a> · <a href="../../tissue-hazards/">Tissue ontology</a></p></header>'
    )
    for title, paragraphs, fig in data["sections"]:
        chunks.append(
            "<section><h2>"
            + esc(title)
            + "</h2>"
            + "".join("<p>" + esc(p) + "</p>" for p in paragraphs)
        )
        if title == "Normal-tissue consequences and expression evidence":
            chunks.append(
                "<dl>"
                + "".join(
                    "<dt><b>" + esc(tier.title()) + "</b></dt><dd>" + esc(definition) + "</dd>"
                    for tier, definition in data["tier_definitions"].items()
                )
                + "</dl>"
            )
        if fig:
            chunks.append(f'<figure><img src="figures/{fig}.svg" alt="{esc(title)}"></figure>')
        if title == "Published nominations through the funnel":
            chunks.append(
                '<p><a href="../cta-source-walkthrough-20261028/figures/publication-funnel.csv">All paper-by-gate counts</a> · <a href="../cta-source-walkthrough-20261028/human-memberships.csv.gz">Protein memberships</a></p><table><thead><tr><th>Paper</th><th>Extracted scope</th></tr></thead><tbody>'
            )
            for row in data["publication_funnel"]:
                if row["source_kind"] == "paper":
                    chunks.append(
                        '<tr><td><a href="'
                        + esc(row["source_url"], quote=True)
                        + '">'
                        + esc(row["citation"])
                        + "</a></td><td>"
                        + esc(row["evidence_scope"])
                        + "</td></tr>"
                    )
            chunks.append("</tbody></table>")
        chunks.append("</section>")
    chunks.append(
        '<h2>Complete gate rules</h2><div class="scroller"><table><thead><tr><th>Gate</th><th>Rule</th><th>Sources</th></tr></thead><tbody>'
    )
    for r in data["criteria"]:
        chunks.append(
            "<tr><td>"
            + esc(r["stage"])
            + "</td><td>"
            + esc(r["rule"])
            + "</td><td>"
            + esc(", ".join(r["source_ids"]))
            + "</td></tr>"
        )
    chunks.append("</tbody></table></div><h2>Worked decisions</h2>")
    for case in data["worked"]:
        chunks.append(
            "<details><summary>"
            + esc(
                case["label"]
                + " · "
                + case["classification"]
                + " · sequence "
                + case["sequence_role"]
            )
            + "</summary><p>"
            + esc(case["assessment_status"] + ": " + "; ".join(case["reasons"]))
            + "</p><p>Translation: "
            + esc(case["translation_reasons"] or "qualified selected record")
            + "</p>"
        )
        for f in case["findings"]:
            chunks.append(
                "<p><b>"
                + esc(f["interpretation_status"] + " / " + f["effect"])
                + "</b>: "
                + esc(f["decision_reason"])
                + ' <a href="'
                + esc(f["source_url"], quote=True)
                + '">'
                + esc(f["source_title"])
                + "</a><br><small>"
                + esc(f["review_depth"] + "; " + f["limitations"])
                + "</small></p>"
            )
        chunks.append("</details>")
    chunks.append(
        '<h2>Explore every classified protein</h2><p>Search labels, gene IDs or reasons. Results retain all overlapping reasons.</p><input id="query" aria-label="Search proteins" placeholder="Protein, gene, reason"><select id="category" aria-label="Classification"><option value="">All categories</option><option>CTA</option><option>Unsafe</option><option>Other</option></select><p id="matched" aria-live="polite"></p><div class="scroller"><table><thead><tr><th>Protein</th><th>Classification</th><th>Status / role</th><th>All recorded reasons</th></tr></thead><tbody id="results"></tbody></table></div>'
    )
    rows = json.dumps(data["rows"], separators=(",", ":")).replace("<", "\\u003c")
    chunks.append(
        '<script type="application/json" id="records">'
        + rows
        + '</script><script>const rows=JSON.parse(document.getElementById("records").textContent);function render(){const q=document.getElementById("query").value.toLowerCase();const c=document.getElementById("category").value;const found=rows.filter(r=>(!c||r.classification.split(":")[0]===c)&&JSON.stringify(r).toLowerCase().includes(q));document.getElementById("matched").textContent=found.length.toLocaleString()+" matching proteins; showing the first 100.";const target=document.getElementById("results");target.replaceChildren();for(const r of found.slice(0,100)){const tr=document.createElement("tr");for(const t of [r.label+" · "+r.gene_ids.join(", "),r.classification,r.assessment_status+" / "+r.shared_sequence,r.reasons.join("; ")||"All declared gates pass"]){const td=document.createElement("td");td.textContent=t;tr.append(td)}target.append(tr)}}document.getElementById("query").addEventListener("input",render);document.getElementById("category").addEventListener("change",render);render();</script>'
    )
    chunks.append(
        '<footer>Sources and audits: <a href="../cta-evidence-migration-20261104/">historical findings and clinical causality</a> · <a href="../translation-qualification-20261104/">translation qualification</a> · <a href="../sequence-bundle-20261104/">verified sequence bundle</a> · <a href="../normal-breadth-20261104/">normal breadth</a> · <a href="../../tissue-hazards/">tissue registry and primary references</a><p><code>'
        + esc(data["coverage"]["policy_id"])
        + "</code></p></footer></main></html>"
    )
    (OUT / "report.html").write_text("".join(chunks))


def render_pdf(data):
    font_dir = Path(matplotlib.get_data_path()) / "fonts/ttf"
    for name, filename in [
        ("DejaVu", "DejaVuSans.ttf"),
        ("DejaVu-Bold", "DejaVuSans-Bold.ttf"),
        ("DejaVu-Italic", "DejaVuSans-Oblique.ttf"),
        ("DejaVu-BoldItalic", "DejaVuSans-BoldOblique.ttf"),
    ]:
        pdfmetrics.registerFont(TTFont(name, str(font_dir / filename)))
    pdfmetrics.registerFontFamily(
        "DejaVu",
        normal="DejaVu",
        bold="DejaVu-Bold",
        italic="DejaVu-Italic",
        boldItalic="DejaVu-BoldItalic",
    )
    styles = getSampleStyleSheet()
    for style in styles.byName.values():
        if hasattr(style, "fontName"):
            style.fontName = "DejaVu-Bold" if "Bold" in style.fontName else "DejaVu"
    styles["Heading1"].keepWithNext = True
    styles["Heading2"].keepWithNext = True
    styles.add(
        ParagraphStyle(
            name="BodyReport",
            fontName="DejaVu",
            fontSize=10.5,
            leading=15,
            spaceAfter=10,
            textColor=colors.HexColor(INK),
        )
    )
    styles.add(
        ParagraphStyle(name="SmallReport", fontName="DejaVu", fontSize=8, leading=10, spaceAfter=6)
    )
    body = styles["BodyReport"]
    small = styles["SmallReport"]
    story = []

    def paragraph(text, style=body):
        return Paragraph(html.escape(text), style)

    story += [
        Paragraph("From the candidate universe to CTAs", styles["Title"]),
        paragraph("OncoRef " + data["release"] + " | 5 November 2026"),
        paragraph("Research candidate selection, evidence status and verified sequence use."),
        Spacer(1, 16),
    ]
    for title, paragraphs, fig in data["sections"]:
        story.append(Paragraph(html.escape(title), styles["Heading1"]))
        story.extend(paragraph(p) for p in paragraphs)
        if title == "Normal-tissue consequences and expression evidence":
            for tier, definition in data["tier_definitions"].items():
                story.append(
                    Paragraph(
                        "<b>" + html.escape(tier.title()) + ":</b> " + html.escape(definition), body
                    )
                )
            story.append(
                Paragraph(
                    '<a href="https://pirl-unc.github.io/oncoref/tissue-hazards/">Central tissue registry, compartment rationales and primary references</a>',
                    small,
                )
            )
        if title == "Evidence has multiple independent meanings":
            story.append(
                Paragraph(
                    '<a href="https://pubmed.ncbi.nlm.nih.gov/23377668/">Morgan et al. 2013: neurological injury report</a> · <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7946352/">Martin et al. 2021: MAGEA12/EPS8L2 re-examination</a>',
                    small,
                )
            )
        if fig:
            story.append(
                Image(
                    str(OUT / "figures" / f"{fig}.png"),
                    width=500,
                    height=350
                    if fig == "publication-funnel"
                    else 330
                    if fig == "independent-evidence"
                    else 250,
                    kind="proportional",
                )
            )
        if title == "Published nominations through the funnel":
            story.append(
                Paragraph(
                    '<a href="https://pirl-unc.github.io/oncoref/audits/cta-source-walkthrough-20261028/figures/publication-funnel.csv">Download every paper-by-gate count and its source scope</a>',
                    small,
                )
            )
            story.append(PageBreak())
            story.append(
                Paragraph("Registered nomination papers and extraction scope", styles["Heading1"])
            )
            for row in data["publication_funnel"]:
                if row["source_kind"] == "paper":
                    story.append(
                        Paragraph(
                            '<a href="'
                            + html.escape(row["source_url"], quote=True)
                            + '">'
                            + html.escape(row["citation"])
                            + "</a>: "
                            + html.escape(row["evidence_scope"]),
                            small,
                        )
                    )
        if title != "Research parameters and coverage limits":
            story.append(PageBreak())
    story.append(Paragraph("Complete executable gate rules", styles["Heading1"]))
    for r in data["criteria"]:
        source_text = (
            ", ".join(r["source_ids"])
            or "Declared policy; see linked source walkthrough and coverage audit."
        )
        story.append(
            KeepTogether(
                [
                    Paragraph(html.escape(r["stage"]), styles["Heading2"]),
                    paragraph(r["rule"]),
                    paragraph("Sources: " + source_text, small),
                ]
            )
        )
    story.append(Spacer(1, 14))
    story.append(Paragraph("Worked decisions and primary sources", styles["Heading1"]))
    for case in data["worked"]:
        story.append(
            Paragraph(
                html.escape(case["label"] + " - " + case["classification"]), styles["Heading2"]
            )
        )
        story.append(
            paragraph(
                case["assessment_status"]
                + "; sequence role: "
                + case["sequence_role"]
                + ". Reasons: "
                + ("; ".join(case["reasons"]) or "All declared gates pass.")
            )
        )
        story.append(
            paragraph(
                "Translation: "
                + (case["translation_reasons"] or "Qualified selected primary-locus record.")
            )
        )
        for f in case["findings"]:
            story.append(
                paragraph(
                    f["interpretation_status"] + " / " + f["effect"] + ": " + f["decision_reason"]
                )
            )
            story.append(
                Paragraph(
                    '<a href="'
                    + html.escape(f["source_url"], quote=True)
                    + '">'
                    + html.escape(f["source_title"])
                    + "</a>",
                    small,
                )
            )
            story.append(paragraph(f["review_depth"] + ". " + f["limitations"], small))

    def footer(canvas, doc):
        canvas.setFont("DejaVu", 7)
        canvas.setFillColor(colors.HexColor(MUTED))
        canvas.drawString(
            42,
            25,
            "OncoRef " + data["release"] + " | CTA definition | " + data["coverage"]["policy_id"],
        )
        canvas.drawRightString(570, 25, str(doc.page))

    path = ROOT / "output/pdf/cta-definition.pdf"
    path.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(
        str(path),
        pagesize=(612, 792),
        leftMargin=44,
        rightMargin=44,
        topMargin=42,
        bottomMargin=44,
        title="OncoRef CTA definition",
        author="OncoRef",
        invariant=1,
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    (OUT / "cta-definition.pdf").write_bytes(path.read_bytes())


if __name__ == "__main__":
    build()
