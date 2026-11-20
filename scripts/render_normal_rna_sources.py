"""Render the measured source-by-tissue inventory and the antigen membership calculation."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cta_report_common import write_csv
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

from oncoref import target_selection

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "docs/audits/target-selection-20261110/rebuild"
OUT = BUILD / "normal-rna"
SOURCES = {
    "gtex_donor_p95": (
        "GTEx P95",
        "GTEx v11. Linear 95th percentile within each anatomical site across donors.",
    ),
    "gtex_individual_tissue_tpm": (
        "GTEx samples",
        "GTEx v11. Individual reproductive samples can witness RA.",
    ),
    "gtex_tissue_rna": (
        "GTEx median",
        "GTEx v11. Published site-level donor medians, with exact member-locus mapping.",
    ),
    "hpa_tissue_samples": (
        "HPA biopsies",
        "HPA v23. Within-biopsy transcript TPM summed by locus. Within-tissue P95 sets the healthy biopsy reference and each reproductive biopsy can witness RA. Donor IDs are unavailable.",
    ),
    "rna_tissue_hpa": ("HPA tissue", "HPA v23 tissue means in TPM."),
    "rna_tissue_gtex": (
        "GTEx via HPA",
        "GTEx tissue summaries distributed in HPA v23. These overlap the GTEx resource.",
    ),
    "rna_brain_hpa": (
        "HPA brain",
        "HPA brain subregion means in TPM, mapped to anatomical regions.",
    ),
    "rna_brain_gtex": ("GTEx brain via HPA", "GTEx brain region summaries distributed by HPA v23."),
    "hpa_immune_cell_samples": (
        "HPA sorted immune",
        "P95 of donor-level native TPM within each sorted blood population.",
    ),
    "hpa_immune_monaco": (
        "Monaco",
        "Mean TPM in sorted blood immune populations. The immune screen also evaluates these measurements.",
    ),
    "dice": (
        "DICE",
        "Supported blood immune-cell TPM, using within-population P95 for single loci and the conservative bound for unpaired protein-group members. The immune screen also uses DICE.",
    ),
    "petropoulos_E-MTAB-3929": (
        "Petropoulos",
        "Preimplantation embryo-lineage means in TPM; stage, cell and independent embryo support apply.",
    ),
    "tyser_E-MTAB-9388": (
        "Tyser",
        "Carnegie stage 7 embryo-lineage means in TPM; cell and embryo support apply.",
    ),
    "hpa_rna_consensus": (
        "HPA consensus",
        "Normal tissue consensus nTPM defines anatomical coverage and contributes to the shared expression screen.",
    ),
    "hpa_single_cell_tissue": (
        "HPA cell clusters",
        "Tissue cell-type cluster nTPM contributes to shared normal-expression and reproductive-context evidence.",
    ),
    "hpa_immune_cell": (
        "HPA immune",
        "Blood immune-cell nTPM evaluated against the tissue-tier limits.",
    ),
    "hpa_brain": (
        "HPA brain nTPM",
        "Brain subregion nTPM contributes to the independent normal-expression screen.",
    ),
    "cornea_collin_adult_counts": (
        "Collin cornea",
        "Adult cornea/conjunctiva cell populations evaluated in counts per 10,000.",
    ),
    "census": (
        "Census donors",
        "CELLxGENE Census donor-resolved single-cell and single-nucleus counts. The rare-cell gate uses within-donor frequency and study replication.",
    ),
}
ROLES = ["Normal", "Disposable", "RA witness", "Reproductive context", "Shared screen", "Thymic"]
COLORS = ["#f1f3f4", "#336689", "#b07a35", "#277763", "#a9cdba", "#785780", "#a05a66"]


def measured_contexts():
    columns = [
        "source_id",
        "source_context",
        "tissue",
        "scope",
        "qc_pass",
        "identity_resolved",
        "observation_kind",
        "cell_type",
    ]
    parts = []
    for filename in ["normal-comparison.csv.gz", "reproductive-comparison.csv.gz"]:
        for chunk in pd.read_csv(
            BUILD / filename, usecols=columns, chunksize=250000, keep_default_na=False
        ):
            d = chunk.loc[chunk.qc_pass & chunk.identity_resolved].drop_duplicates().copy()
            d["family"] = d.source_id.str.split(":").str[0]
            d.loc[d.family.str.startswith("dice_"), "family"] = "dice"
            d["role"] = d.scope.map(
                {
                    "normal": "Normal",
                    "thymic": "Thymic",
                    "disposable": "Disposable",
                    "developmental": "RA witness",
                    "reproductive": "Reproductive context",
                }
            )
            sample = d.observation_kind.isin(["individual_tissue_sample", "individual_gtex_sample"])
            d.loc[d.scope.eq("reproductive") & sample, "role"] = "RA witness"
            parts.append(
                d[["family", "source_id", "tissue", "source_context", "cell_type", "role"]]
            )
    context = pd.concat(parts).drop_duplicates()
    if context.role.isna().any() or set(context.family) - SOURCES.keys():
        raise ValueError("Unrecognized RNA source or anatomical scope")
    return context


def independent_contexts():
    frames = []

    def add(key, frame):
        frames.append(frame.assign(family=key, source_id=key, role="Shared screen"))

    frames.append(pd.read_csv(BUILD / "hpa-normal-source-contexts.csv", keep_default_na=False))
    frozen = (
        ROOT
        / "docs/audits/cta-transcript-selection-20261108/inputs/frozen/immune-cornea-contexts.csv.gz"
    )
    immune = pd.read_csv(frozen, usecols=["source", "tissue", "cell_type"]).drop_duplicates()
    for key, group in immune.groupby("source"):
        label = "dice" if key.startswith("dice_") else key
        frames.append(
            group.rename(columns={"source": "source_id"}).assign(
                family=label, source_context=group.cell_type, role="Shared screen"
            )
        )
    brain = pd.read_csv(
        ROOT / "docs/audits/cta-selection-20261028/alternative-evidence/hpa-brain-anatomy.csv"
    )
    add(
        "hpa_brain",
        brain.rename(columns={"Brain region": "tissue", "Subregion": "source_context"})[
            ["tissue", "source_context"]
        ],
    )
    census = pd.read_csv(ROOT / "docs/audits/normal-cell-atlas-20261110/rare-cell/panel.csv")
    add(
        "census",
        census.assign(source_context=census.tissue)[["tissue", "source_context"]].drop_duplicates(),
    )
    return pd.concat(frames, ignore_index=True).fillna("")


TISSUE_ALIASES = {
    "hippocampus": "hippocampal formation",
}
TISSUE_GROUPS = [
    (
        "Somatic organs",
        [
            "adipose tissue",
            "adrenal gland",
            "blood vessel",
            "bronchus",
            "colon",
            "duodenum",
            "esophagus",
            "heart",
            "heart muscle",
            "kidney",
            "liver",
            "lung",
            "pancreas",
            "parathyroid gland",
            "rectum",
            "salivary gland",
            "skeletal muscle",
            "skin",
            "small intestine",
            "smooth muscle",
            "stomach",
            "thyroid gland",
            "tongue",
            "ureter",
            "urinary bladder",
            "vascular",
        ],
    ),
    (
        "Blood & immune",
        ["blood", "bone marrow", "lymph node", "pbmc", "spleen"],
    ),
    (
        "Disposable & thymic",
        ["appendix", "gallbladder", "thymus", "tonsil"],
    ),
    (
        "Brain, nerve & eye",
        [
            "amygdala",
            "basal ganglia",
            "brain",
            "caudate",
            "cerebellum",
            "cerebral cortex",
            "choroid plexus",
            "cornea/conjunctiva",
            "eye",
            "hippocampal formation",
            "hypothalamus",
            "medulla oblongata",
            "midbrain",
            "peripheral nerve",
            "pituitary gland",
            "pons",
            "retina",
            "spinal cord",
            "substantia nigra",
            "thalamus",
            "white matter",
        ],
    ),
    (
        "Reproductive & developmental",
        [
            "breast",
            "cervix",
            "embryo",
            "endometrium",
            "epididymis",
            "fallopian tube",
            "ovary",
            "placenta",
            "prostate",
            "seminal vesicle",
            "testis",
            "uterus",
            "vagina",
        ],
    ),
]


def _grouped_tissue_order(tissues):
    """Order tissues by anatomical scope so coverage blocks align across panels."""
    present = set(tissues)
    ordered = []
    assigned = set()
    for group_name, members in TISSUE_GROUPS:
        hits = [t for t in members if t in present]
        if hits:
            ordered.append((group_name, hits))
            assigned.update(hits)
    leftover = sorted(present - assigned)
    if leftover:
        ordered.append(("Other tissues", leftover))
    return ordered


def _draw_coverage_panel(ax, grid, grp_list, keys):
    row_order = [t for _, members in grp_list for t in members]
    panel = grid.reindex(index=row_order)
    ax.imshow(
        panel,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(COLORS),
        norm=BoundaryNorm(np.arange(-0.5, len(COLORS) + 0.5), len(COLORS)),
    )
    ax.set_xticks(
        np.arange(len(keys)),
        [SOURCES[k][0] for k in keys],
        rotation=48,
        ha="left",
        fontsize=9.0,
    )
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(panel)), panel.index, fontsize=9.0)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(keys)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(panel)), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.95)
    ax.tick_params(which="minor", length=0)
    cursor = 0
    for idx, (group_name, members) in enumerate(grp_list):
        if idx > 0:
            ax.axhline(cursor - 0.5, color="#334854", linewidth=1.3)
        mid = cursor + (len(members) - 1) / 2
        ax.text(
            len(keys) - 0.35,
            mid,
            f"  {group_name}",
            va="center",
            ha="left",
            fontsize=8.4,
            fontweight="bold",
            color="#3b505c",
            clip_on=False,
        )
        cursor += len(members)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_grid(data, keys, filename, *, title="", split_stems=None):
    cells = data.assign(
        tissue=data.tissue.replace(TISSUE_ALIASES),
        code=data.role.map({r: i + 1 for i, r in enumerate(ROLES)}),
    )
    grid = (
        cells.pivot_table(index="tissue", columns="family", values="code", aggfunc="min")
        .reindex(columns=keys)
        .fillna(0)
    )
    groups = _grouped_tissue_order(grid.index)
    # Partition groups cleanly into two balanced halves without splitting a group.
    total_rows = sum(len(members) for _, members in groups)
    left_groups, right_groups, left_count = [], [], 0
    for grp in groups:
        if left_groups and left_count + len(grp[1]) > total_rows / 2 + 4:
            right_groups.append(grp)
        else:
            left_groups.append(grp)
            left_count += len(grp[1])
    if not right_groups:
        half = max(1, len(groups) // 2)
        left_groups, right_groups = groups[:half], groups[half:]

    used = sorted(cells.code.unique())
    panel_specs = [left_groups, right_groups]

    if split_stems is not None:
        for (split_name, split_title), grp_list in zip(split_stems, panel_specs, strict=True):
            n_rows = sum(len(m) for _, m in grp_list)
            fig_s, ax_s = plt.subplots(
                figsize=(max(8.8, len(keys) * 0.68), 2.6 + 0.27 * n_rows),
                layout="constrained",
            )
            _draw_coverage_panel(ax_s, grid, grp_list, keys)
            ax_s.set_title(split_title, loc="left", fontsize=11.5, fontweight="bold", pad=16)
            fig_s.legend(
                handles=[Patch(color=COLORS[c], label=ROLES[c - 1]) for c in used],
                loc="outside lower center",
                ncol=len(used),
                frameon=False,
                fontsize=9.2,
            )
            fig_s.savefig(OUT / split_name, dpi=165, bbox_inches="tight")
            plt.close(fig_s)

    max_rows = max(sum(len(m) for _, m in g) for g in panel_specs)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(max(11.5, len(keys) * 1.18), 3.3 + 0.25 * max_rows),
        layout="constrained",
    )
    for ax, grp_list in zip(axes, panel_specs, strict=True):
        _draw_coverage_panel(ax, grid, grp_list, keys)
    if title:
        fig.suptitle(title, fontsize=12.5, fontweight="bold", color="#1d2f3a")
    fig.legend(
        handles=[Patch(color=COLORS[c], label=ROLES[c - 1]) for c in used],
        loc="outside lower center",
        ncol=len(used),
        frameon=False,
        fontsize=9,
    )
    fig.savefig(OUT / filename, dpi=165, bbox_inches="tight")
    plt.close(fig)


def plot_reference_sources(catalog):
    fields = ["entity_id", "source_id", "scope", "value", "qc_pass", "identity_resolved"]
    pieces = []
    for chunk in pd.read_csv(BUILD / "normal-comparison.csv.gz", usecols=fields, chunksize=250000):
        normal = chunk.loc[
            chunk.scope.isin(["normal", "disposable"]) & chunk.qc_pass & chunk.identity_resolved
        ].copy()
        normal["family"] = normal.source_id.str.split(":").str[0]
        normal.loc[normal.family.str.startswith("dice_"), "family"] = "dice"
        pieces.append(normal.groupby(["entity_id", "family"]).value.max())
    peaks = pd.concat(pieces).groupby(level=[0, 1]).max().unstack()
    peaks["floor"] = 1.0
    expected = (
        catalog.set_index("entity_id")
        .loc[peaks.index, ["normal_reference_tpm", "disposable_max_tpm"]]
        .max(axis=1)
        .clip(lower=1)
    )
    if not np.allclose(peaks.max(axis=1), expected):
        raise ValueError("Source peaks disagree with catalog denominator")
    wins = peaks.eq(peaks.max(axis=1), axis=0)
    counts = pd.DataFrame(
        {
            "unique": wins.loc[wins.sum(axis=1).eq(1)].sum(),
            "tied": wins.loc[wins.sum(axis=1).gt(1)].sum(),
        }
    )
    counts.index.name = "source"
    write_csv(counts.reset_index(), OUT / "denominator-sources.csv")
    counts = counts.loc[counts.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(8, 4.6), layout="constrained")
    labels = ["1 TPM floor" if k == "floor" else SOURCES[k][0] for k in counts.index]
    ax.barh(labels, counts["unique"], color=COLORS[1], label="Sole source")
    ax.barh(
        labels,
        counts["tied"],
        left=counts["unique"],
        color="#a3bfcf",
        label="Tied with another source",
    )
    for i, total in enumerate(counts.sum(axis=1)):
        ax.annotate(
            str(total),
            (total, i),
            xytext=(4, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )
    ax.set_xlim(0, counts.sum(axis=1).max() * 1.13)
    ax.set_xlabel(f"Protein groups attaining the combined reference ({len(wins):,} assessed)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", frameon=False)
    fig.savefig(OUT / "denominator-sources.png", dpi=165)
    plt.close(fig)


def build():
    OUT.mkdir(exist_ok=True, parents=True)
    fold = measured_contexts()
    shared = independent_contexts()
    medians = fold.loc[fold.family.eq("gtex_tissue_rna")].assign(role="Shared screen")
    shared = pd.concat([shared, medians], ignore_index=True)
    inventory = (
        pd.concat([fold, shared])
        .drop_duplicates()
        .sort_values(["family", "tissue", "source_id", "source_context", "cell_type", "role"])
    )
    write_csv(inventory, OUT / "source-tissues.csv")
    plot_grid(
        fold,
        [k for k in SOURCES if k in set(fold.family)],
        "tpm-coverage.png",
        title="Native-TPM tissue coverage and analytical role across sources",
        split_stems=[
            (
                "tpm-coverage-somatic.png",
                "Native-TPM coverage: somatic, immune, disposable, and thymic tissues",
            ),
            (
                "tpm-coverage-neural-repro.png",
                "Native-TPM coverage: brain, eye, nerve, and reproductive/developmental tissues",
            ),
        ],
    )
    plot_grid(
        shared,
        [k for k in SOURCES if k in set(shared.family)],
        "shared-coverage.png",
        title="Independent normal RNA screens: tissue coverage by source",
    )
    cat = target_selection.results()
    plot_reference_sources(cat)
    ca, ra = cat.is_cancer_antigen, cat.is_reproductive_antigen
    lines = [
        "# Normal RNA sources and antigen membership",
        "",
        "Every native-TPM reference is evaluated within its anatomical context. The highest supported within-site statistic across sources sets N. Disposable tissue supplies a separate D. Both antigen routes use **H₀ = max(N, D, 1 TPM)**. The source panels overlap; measurements are never added across studies.",
        "",
        "[Cell and isoform evidence](../../../../curation.md#cell-and-isoform-evidence) describes expression in olfactory neurons, developing teeth, bone, and peri-implantation embryos, with the assay and age for each observation.",
        "",
        "![Tissue coverage and role by native TPM source](tpm-coverage.png)",
        "",
        "A colored cell indicates measured coverage for at least one assessed group. Reproductive context denotes a tissue summary; qualifying RA evidence requires an individual reproductive sample or supported embryo-lineage mean. Gray denotes no qualifying coverage in that source and role. The inventory lists every source-native site and cell type used in these comparison records.",
        "",
        "[Download source, tissue, native context, cell type and role](source-tissues.csv)",
        "",
        "| Source | Measurement and use |",
        "|---|---|",
    ]
    lines.extend(f"| {label} | {detail} |" for label, detail in SOURCES.values())
    lines.extend(
        [
            "",
            "## Tissues by source",
            "",
            "Tissue names follow each panel's anatomical mapping. The downloadable inventory preserves finer site, subregion and cell-type labels.",
            "",
            "| Source | Tissues represented |",
            "|---|---|",
        ]
    )
    for key, (label, _) in SOURCES.items():
        tissues = sorted(inventory.loc[inventory.family.eq(key), "tissue"].unique())
        lines.append(f"| {label} | {', '.join(tissues)} |")
    qc = pd.read_csv(OUT / "carryover-summary.csv")
    lines.extend(
        [
            "",
            "## Cross-tissue RNA in individual libraries",
            "",
            "P denotes a pancreatic enzyme signature (PRSS1, CLPS, PNLIP, CELA3A). E denotes a squamous signature observed in esophageal mucosa (KRT4, KRT13). These markers come from [Nieuwenhuis et al.](https://www.nature.com/articles/s41467-020-15821-9).",
            "",
            "OncoRef compares each library with samples from the same anatomical site or sorted immune population. The panel requires at least three complete libraries. Every marker must exceed its within-site interpolated P95 and exceed the site median by at least 0.1% of its native GTEx median. The source organ is skipped. A blank flag means the panel could not be assessed. The thresholds here are OncoRef choices.",
            "",
            "These are advisory sample flags. Normal epithelial content can produce E, particularly in tissues with squamous epithelium. Shared marker expression alone cannot identify a transfer event. Tissue means, unpaired DICE ranks and single-cell UMI counts cannot support this within-library bulk TPM comparison.",
            "",
            "| Source | Libraries | P flagged / assessed | E flagged / assessed | Either |",
            "|---|---:|---:|---:|---:|",
            *(
                f"| {row.source} | {row.libraries:,} | {row.P:,} / {row.P_assessed:,} | {row.E:,} / {row.E_assessed:,} | {row.either:,} |"
                for row in qc.itertuples()
            ),
            "",
            "[GTEx libraries](../gtex-library-qc.csv.gz) · [HPA tissue biopsies](../normal-library-qc.csv) · [HPA sorted immune libraries](immune-library-qc.csv) · [GTEx batch matches and percentile sensitivity](../../../gtex-cross-tissue-signals/index.md)",
        ]
    )
    lines.extend(
        [
            "",
            "## Independent normal RNA screens",
            "",
            "![Independent RNA screen coverage](shared-coverage.png)",
            "",
            "The tissue-tier, immune, corneal and rare-cell screens can prevent shared eligibility before enrichment is evaluated. Some native TPM sources also inform these screens. Cluster nTPM, CP10K, donor UMIs and protein staining retain their own evidence rules. A bulk percentile cannot override an exclusion.",
            "",
            "## From measurements to memberships",
            "",
            "For GTEx v11 normal sites, sum the measured loci encoding each identical protein within a library. Give each donor one value per site (median of replicate libraries), then take the linearly interpolated 95th percentile of the measured donors. Interpolation produces a tissue statistic without a sample or donor ID. HPA biopsies contribute the linear P95 of measured biopsy TPM per healthy tissue. Missing required tissue coverage blocks qualification.",
            "",
            "![Sources setting the reference for each assessed protein](denominator-sources.png)",
            "",
            "| Route after shared eligibility | Definition | Protein groups |",
            "|---|---|---:|",
            f"| CA | At least one QC-passing tumor reaches 10 × H₀ | {ca.sum():,} |",
            f"| RA | At least one reproductive sample or supported embryo-lineage observation reaches 10 × H₀ | {ra.sum():,} |",
            f"| CTA | Both CA and RA | {(ca & ra).sum():,} |",
            f"| CA ∪ RA | Either route | {(ca | ra).sum():,} |",
            "",
            "![Cancer and reproductive antigen overlap](../summary/venn.png)",
            "",
            "[Set overlap and eligibility funnel](../summary.md) · [Exact comparison contract](../comparison.json) · [Census study inventory](../../../normal-cell-atlas-20261110/index.md)",
            "",
        ]
    )
    (OUT / "index.md").write_text("\n".join(lines))


if __name__ == "__main__":
    build()
