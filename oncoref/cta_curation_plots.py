# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""CTA curation documentation figures over oncoref's packaged CTA table.

Twelve figures describe complete published candidate sets, source intersections,
the full intake funnel, per-source attrition, and the protein/RNA evidence axes.
Source memberships include unmapped and noncoding entries; Venn diagrams use
unique mapped protein-coding loci before HPA filtering.

Two frames are in play and they are not interchangeable:

``_evidence()``
    The full coding union of the minimum-cover papers, including rows that fail
    curation. Historical-only nominations remain in a separate audit export.

``_curated()``
    The same table after :mod:`oncoref.cta` drops explicitly non-CTA genes and
    joins the machine-readable specificity adjudication. Its
    ``specificity_action == "include_default"`` rows are exactly what
    ``cta_gene_ids()`` returns, so the stage funnel lands on the shipped set
    rather than on the raw ``passes_filters`` column (which is one stage short).

The standard oncoref installation includes ``matplotlib_venn`` for the
source-overlap figure.
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd

from . import figure_style
from .cta_tissues import HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS
from .figure_style import ACCENT, DROP, KEPT, THRESHOLD, WEAK
from .load_dataset import get_data
from .version import DATA_VERSION, __version__

# Every complete imported list contributes upstream of the HPA/default gates.
# da Silva's broad nominations and final tumor-positive subset remain distinct.
LANDSCAPE_SOURCES = {
    "Wang 2016": lambda tags: "Wang2016_CT" in tags,
    "Bruggeman 2018": lambda tags: "Bruggeman2018_GC" in tags,
    "da Silva 2017 · CT": lambda tags: "daSilva2017_CT" in tags,
    "da Silva 2017 · testis-biased": lambda tags: "daSilva2017_testis_biased" in tags,
    "Jamin 2021": lambda tags: bool(tags & {"Jamin2021_CT", "Jamin2021_core"}),
    "Chang 2019 · TGCT": lambda tags: "Chang2019_TGCT" in tags,
    "Carter 2023": lambda tags: "Carter2023_CT" in tags,
    "Seager 2024 · panel": lambda tags: "Seager2024_CTA" in tags,
}
PRIMARY_SOURCES = {
    "Wang 2016": "10.1038/ncomms10499",
    "Bruggeman 2018": "10.1038/s41388-018-0357-2",
    "da Silva 2017": "10.18632/oncotarget.21715",
    "Jamin 2021": "10.1002/1878-0261.12900",
    "Chang 2019": "10.1002/cam4.2223",
    "Carter 2023": "10.1136/jitc-2023-007935",
    "Seager 2024": "10.1186/s12967-024-04918-0",
    "Gong 2021": "10.1038/s41467-021-22695-y",
    "Loriot 2025": "10.1371/journal.pgen.1011734",
    "Bai 2016": "10.1158/0008-5472.CAN-16-0225",
}


# Deflated-RNA-fraction threshold each protein-reliability tier must clear.
RELIABILITY_THRESHOLD = {
    **{k: v for k, v in HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS.items() if k != "Missing"},
    "no data": HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS["Missing"],
}
RELIABILITY_ORDER = ["no data", "Uncertain", "Approved", "Supported", "Enhanced"]

# Hyphenated filenames match pirlygenes/docs references.
FILENAMES = {
    "source_venn": "cta-source-venn.png",
    "stage_funnel": "cta-stage-funnel.png",
    "filter_funnel": "cta-filter-funnel.png",
    "filter_outcome": "cta-filter-outcome.png",
    "deflated_dist": "cta-deflated-frac-dist.png",
    "protein_vs_rna": "cta-protein-vs-rna.png",
    "legacy_source_venn": "cta-legacy-source-venn.png",
    "landscape_source_venn": "cta-landscape-source-venn.png",
    "source_overlap": "cta-source-overlap.png",
    "placental_source_overlap": "cta-placental-source-overlap.png",
    "publication_funnel": "cta-publication-funnel.png",
    "placental_evidence_coverage": "cta-placental-evidence-coverage.png",
}


def _evidence():
    """Complete coding union of the minimum-cover papers, before filtering."""
    from .cta import cta_evidence
    from .cta_provenance import candidate_evidence

    raw = candidate_evidence()
    reviewed = cta_evidence()
    columns = ["Ensembl_Gene_ID", *[c for c in reviewed if c.startswith("specificity_")]]
    return raw.merge(reviewed[columns], on="Ensembl_Gene_ID", how="left", validate="one_to_one")


def _curated():
    """CTA table with non-CTA exclusions dropped and specificity decisions joined."""
    from . import cta

    return cta._cta_with_specificity_frame().copy()


def _bool_series(series):
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def _tag_sets(df):
    """Exact paper-level membership, restricted to the supplied gene frame."""
    from .cta_provenance import selected_membership

    refs = selected_membership()
    unknown = set(refs.doi) - set(PRIMARY_SOURCES.values())
    if unknown:
        raise ValueError(f"Selected papers need plot labels: {sorted(unknown)}")
    ids = set(df.Ensembl_Gene_ID)
    return {
        name: set(refs.loc[refs.doi.eq(doi), "Ensembl_Gene_ID"]) & ids
        for name, doi in PRIMARY_SOURCES.items()
        if doi in set(refs.doi)
    }


def _historical_sets(df):
    tags = df.source_databases.fillna("").str.split(";").map(set)
    predicates = {
        "CTpedia": lambda x: "CTpedia" in x,
        "CTexploreR": lambda x: bool(x & {"CTexploreR_CT", "CTexploreR_CTP"}),
        "daSilva2017_protein": lambda x: "daSilva2017_protein" in x,
        "placental_antigen": lambda x: "placental_antigen" in x,
    }
    return {k: set(df.loc[tags.map(pred), "Ensembl_Gene_ID"]) for k, pred in predicates.items()}


def _per_source_counts(df):
    from . import cta

    sets = _tag_sets(df)
    sets["Other prior nominations"] = set(df.Ensembl_Gene_ID) - set().union(*sets.values())
    default = cta.cta_gene_ids()
    rows = []
    for name, members in sets.items():
        if not members:
            continue
        sub = df[df["Ensembl_Gene_ID"].astype(str).isin(members)]
        passes = sub.Ensembl_Gene_ID.isin(default)
        hpa = cta.passes_filters_mask(sub) & sub.Ensembl_Gene_ID.isin(cta.cta_unfiltered_gene_ids())
        weak = _bool_series(sub["never_expressed"])
        rows.append(
            {
                "source": name,
                "total": len(sub),
                "kept_confident": int((passes & ~weak).sum()),
                "kept_weak": int((passes & weak).sum()),
                "excluded": int((~passes).sum()),
                "default_panel": int(passes.sum()),
                "hpa_pass_outside_default": int((hpa & ~passes).sum()),
                "family_or_hpa_excluded": int((~hpa).sum()),
            }
        )
    return sorted(rows, key=lambda r: r["total"], reverse=True)


STAGES = {
    "nominated": "source union",
    "mapped": "canonical ID mapped",
    "protein_coding": "protein-coding candidates",
    "family_eligible": "non-CTA removed",
    "hpa_restriction": "HPA restriction",
    "default_panel": "specificity + expression",
}


def stage_membership(df=None):
    """Distinct source identities, including unmapped/noncoding nominations.

    Mapped aliases collapse to one canonical locus. Unknown identities stay in the
    first-stage audit; they never become measured zeros or pass downstream filters.
    """
    from . import cta
    from .cta_provenance import selected_membership

    raw = _evidence() if df is None else df.copy()
    if raw.Ensembl_Gene_ID.duplicated().any():
        raise ValueError("Duplicate candidate gene IDs")
    refs = selected_membership().fillna("")
    rows = []
    for r in refs.to_dict("records"):
        gid = r["Ensembl_Gene_ID"]
        native = r["source_gene_id"] or r["source_symbol"]
        kind = r.get("source_id_type") or (
            "ensembl" if str(native).startswith("ENSG") else "symbol"
        )
        rows.append(
            {
                "identity": gid or f"unmapped:{kind}:{native}",
                "Ensembl_Gene_ID": gid,
                "Symbol": r["Symbol"] or r["source_symbol"],
                "biotype": r["biotype"],
                "source_databases": r["source_tag"],
            }
        )
    rows.extend(
        {
            "identity": r.Ensembl_Gene_ID,
            "Ensembl_Gene_ID": r.Ensembl_Gene_ID,
            "Symbol": r.Symbol,
            "biotype": r.biotype,
            "source_databases": r.source_databases,
        }
        for r in raw.itertuples(index=False)
    )
    records = []
    for _, group in pd.DataFrame(rows).groupby("identity", sort=True):
        record = group.iloc[0].to_dict()
        record["source_databases"] = ";".join(
            sorted(
                {
                    tag
                    for value in group.source_databases
                    for tag in str(value).split(";")
                    if tag and tag != "nan"
                }
            )
        )
        records.append(record)
    result = pd.DataFrame(records)
    ids = result.Ensembl_Gene_ID
    result["nominated"] = True
    result["mapped"] = ids.ne("")
    result["protein_coding"] = result.mapped & result.biotype.eq("protein_coding")
    # Every mapped coding nomination must actually have been evaluated.
    missing = set(ids[result.protein_coding]) - set(raw.Ensembl_Gene_ID)
    if missing:
        raise ValueError(f"Unassessed publication candidates: {sorted(missing)[:8]}")
    result["family_eligible"] = result.protein_coding & ids.isin(cta.cta_unfiltered_gene_ids())
    result["hpa_restriction"] = result.family_eligible & ids.isin(
        raw.loc[cta.passes_filters_mask(raw), "Ensembl_Gene_ID"]
    )
    result["default_panel"] = ids.isin(cta.cta_gene_ids())
    if not cta.cta_gene_ids() <= set(ids) or (result.default_panel & ~result.hpa_restriction).any():
        raise ValueError("Default panel cannot be represented by a nested intake funnel")
    return result


def stage_counts():
    """Unique-identity funnel from complete source intake to the public default."""
    membership = stage_membership()
    rows, previous = [], len(membership)
    for key, label in STAGES.items():
        n = int(membership[key].sum())
        rows.append((label, n, previous - n))
        previous = n
    return rows


def _save(fig, path, plt):
    from matplotlib.text import Text

    for text in fig.findobj(Text):
        text.set_fontsize(text.get_fontsize() * _FONT_SCALE)
    figure_style.save(fig, path, keep_titles=True)
    figure_style.save(fig, Path(path).with_suffix(".pdf"), keep_titles=True)
    plt.close(fig)


def _fig_source_venn(df, path, plt):
    from matplotlib_venn import venn3

    sets = _tag_sets(df)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    keys = ("Wang 2016", "Bruggeman 2018", "da Silva 2017")
    v = venn3(
        [sets[k] for k in keys],
        set_labels=keys,
        set_colors=(KEPT, ACCENT, WEAK),
        alpha=0.55,
        ax=ax,
    )
    for text in list(v.set_labels or []) + list(v.subset_labels or []):
        if text is not None:
            text.set_fontsize(10)
    ax.set_title(
        "Published CT / germ-cell cancer sets\nMapped protein-coding nominations, before HPA filtering"
    )
    _save(fig, path, plt)


def _fig_stage_funnel(df, path, plt):
    """Sequential attrition from source intake through the curation stages."""
    stages = stage_counts()
    labels = [s[0] for s in stages]
    remaining = np.array([s[1] for s in stages])
    dropped = np.array([s[2] for s in stages])
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.barh(y, remaining, color=KEPT, height=0.62)
    ax.barh(y, dropped, left=remaining, color=DROP, height=0.62)
    for i, (n, d) in enumerate(zip(remaining, dropped)):
        ax.text(n - max(remaining) * 0.015, i, str(n), va="center", ha="right", color="white")
        if d:
            ax.text(n + d + max(remaining) * 0.012, i, f"−{d}", va="center", fontsize=9)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Distinct source identities / canonical genes")
    ax.set_title("Full union of the ten minimum-cover papers")
    ax.set_xlim(0, max(remaining) * 1.1)
    _save(fig, path, plt)


def _fig_filter_funnel(df, path, plt):
    rows = _per_source_counts(df)
    labels = [r["source"] for r in rows]
    kept = np.array([r["kept_confident"] + r["kept_weak"] for r in rows])
    dropped = np.array([r["excluded"] for r in rows])
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.62 * len(labels) + 1.4))
    ax.barh(y, kept, color=KEPT, label="default panel", height=0.62)
    ax.barh(y, dropped, left=kept, color=DROP, label="not retained", height=0.62)
    for i, r in enumerate(rows):
        ax.text(
            r["total"] + max(kept + dropped) * 0.012,
            i,
            f"{kept[i]}/{r['total']}",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("genes")
    ax.set_title("Default-panel retention by source")
    ax.set_xlim(0, max(kept + dropped) * 1.14)
    ax.legend(loc="lower right")
    _save(fig, path, plt)


def _fig_filter_outcome(df, path, plt):
    rows = _per_source_counts(df)
    labels = [r["source"] for r in rows]
    conf = np.array([r["default_panel"] for r in rows])
    weak = np.array([r["hpa_pass_outside_default"] for r in rows])
    excl = np.array([r["family_or_hpa_excluded"] for r in rows])
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.62 * len(labels) + 1.4))
    ax.barh(y, conf, color=KEPT, label="default panel", height=0.62)
    ax.barh(y, weak, left=conf, color=WEAK, label="HPA pass, outside default", height=0.62)
    ax.barh(y, excl, left=conf + weak, color=DROP, label="family / HPA exclusion", height=0.62)
    for i, row in enumerate(rows):
        ax.text(
            row["total"] + max(r["total"] for r in rows) * 0.012,
            i,
            f"{row['default_panel']}/{row['total']}",
            va="center",
            fontsize=9,
        )
    ax.set_xlim(0, max(r["total"] for r in rows) * 1.18)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("genes")
    ax.set_title("Default-panel outcomes by source")
    ax.legend(loc="lower right")
    _save(fig, path, plt)


def _fig_deflated_dist(df, path, plt):
    frac = df["rna_deflated_reproductive_frac"].astype(float).to_numpy()
    passes = _bool_series(df["passes_filters"]).to_numpy()
    bins = np.linspace(0, 1, 41)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.hist(
        [frac[passes], frac[~passes]],
        bins=bins,
        stacked=True,
        color=[KEPT, DROP],
        label=["HPA gate pass", "HPA gate fail"],
    )
    for thr in sorted(set(RELIABILITY_THRESHOLD.values())):
        ax.axvline(thr, color=THRESHOLD, ls="--", lw=0.7, alpha=0.7)
    ax.set_xlabel("deflated reproductive fraction")
    ax.set_ylabel("genes")
    ax.set_title("HPA restriction gate across all coding candidates")
    ax.legend(loc="upper left")
    fig.text(
        0.5,
        -0.025,
        f"{int(np.isnan(frac).sum())} candidates lack RNA fractions and are omitted here. HPA pass precedes final-panel curation.",
        ha="center",
        fontsize=8,
    )
    _save(fig, path, plt)


def _fig_protein_vs_rna(df, path, plt):
    from matplotlib.lines import Line2D

    frac = df["rna_deflated_reproductive_frac"].astype(float)
    rel = df["protein_reliability"].fillna("no data").astype(str)
    passes = _bool_series(df["passes_filters"])
    fig, ax = plt.subplots(figsize=(7, 4.2))
    rng = np.random.default_rng(0)
    for i, tier in enumerate(RELIABILITY_ORDER):
        m = rel == tier
        if not m.any():
            continue
        x = i + (rng.random(int(m.sum())) - 0.5) * 0.6
        ax.scatter(x, frac[m], s=12, alpha=0.6, linewidths=0, c=np.where(passes[m], KEPT, DROP))
        thr = RELIABILITY_THRESHOLD.get(tier)
        if thr is not None:
            ax.plot([i - 0.42, i + 0.42], [thr, thr], color=THRESHOLD, lw=1.8)
    ax.set_xticks(range(len(RELIABILITY_ORDER)), RELIABILITY_ORDER)
    ax.set_xlabel("HPA protein reliability")
    ax.set_ylabel("deflated reproductive fraction")
    ax.set_title("HPA protein reliability and reproductive RNA")
    ax.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="", color=KEPT, label="HPA gate pass"),
            Line2D([], [], marker="o", linestyle="", color=DROP, label="HPA gate fail"),
            Line2D([], [], color=THRESHOLD, label="adaptive RNA threshold"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=3,
        fontsize=8,
    )
    fig.text(
        0.5,
        -0.10,
        f"{int(frac.isna().sum())} candidates lack RNA fractions and are omitted here. HPA pass precedes final-panel curation.",
        ha="center",
        fontsize=8,
    )
    _save(fig, path, plt)


def source_overlap_counts(df=None):
    """Pairwise intersections count unique canonical loci, never memberships."""
    from . import cta

    sets = _tag_sets(_evidence() if df is None else df)
    default = cta.cta_gene_ids()
    return pd.DataFrame(
        [
            {
                "source_a": a,
                "source_b": b,
                "candidate_overlap": len(x & y),
                "default_overlap": len(x & y & default),
            }
            for a, x in sets.items()
            for b, y in sets.items()
        ]
    )


def _venn(ax, sets, title):
    from matplotlib_venn import venn3

    diagram = venn3(
        list(sets.values()),
        set_labels=tuple(sets),
        ax=ax,
        set_colors=(KEPT, ACCENT, WEAK),
        alpha=0.55,
    )
    for label in list(diagram.set_labels or []) + list(diagram.subset_labels or []):
        if label is not None:
            label.set_fontsize(9)
    if diagram.set_labels and diagram.set_labels[2] is not None:
        label = diagram.set_labels[2]
        x, y = label.get_position()
        label.set_position((x, y - 0.04))
    ax.set_title(title, fontsize=11, pad=20)


def _fig_legacy_source_venn(df, path, plt):
    sets = _historical_sets(get_data("cancer-testis-antigens"))
    keys = ("CTpedia", "CTexploreR", "daSilva2017_protein")
    fig, ax = plt.subplots(figsize=(7, 6))
    _venn(
        ax,
        {k: sets[k] for k in keys},
        "Historical tags (audit only)\nNot the paper-union intake or independent validation",
    )
    _save(fig, path, plt)


def _fig_landscape_source_venn(df, path, plt):
    sets = _tag_sets(df)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
    _venn(
        axes[0],
        {k: sets[k] for k in ("Jamin 2021", "Chang 2019", "Carter 2023")},
        "Additional complete candidate lists",
    )
    _venn(
        axes[1],
        {k: sets[k] for k in ("Gong 2021", "Loriot 2025", "Seager 2024")},
        "Reproductive enrichment and cancer-germline lists",
    )
    fig.text(
        0.5,
        0.015,
        "Unique mapped coding genes before HPA filtering. Bai's single EGFL6 nomination is shown in the all-paper matrix.\n"
        "Normal-only reproductive nominations remain candidates; shared probes do not establish locus-specific validation.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.09, 1, 1), w_pad=2)
    _save(fig, path, plt)


def placental_source_sets(df=None):
    from .cta_sources import BRADLEY, GONG_NC, GONG_PC, publication_membership

    refs = publication_membership()
    coding = refs[refs.biotype.eq("protein_coding")]
    raw = get_data("cancer-testis-antigens")
    return {
        "Gong 2021\nS5 + S6": set(
            coding.loc[coding.source_tag.isin([GONG_PC, GONG_NC]), "Ensembl_Gene_ID"]
        ),
        "Bradley 2020\nFigure 3": set(coding.loc[coding.source_tag.eq(BRADLEY), "Ensembl_Gene_ID"]),
        "Prior placental\nnominations": _historical_sets(raw)["placental_antigen"],
    }


def _fig_placental_source_overlap(df, path, plt):
    fig, ax = plt.subplots(figsize=(8, 6.2))
    _venn(
        ax,
        placental_source_sets(df),
        "Placental source audit (includes unselected Bradley)\nMapped coding candidates before HPA filtering",
    )
    _save(fig, path, plt)


def _fig_source_overlap(df, path, plt):
    counts = source_overlap_counts(df)
    keys = list(_tag_sets(df))
    fig, axes = plt.subplots(1, 2, figsize=(19, 10.5))
    for ax, field, title in zip(
        axes,
        ("candidate_overlap", "default_overlap"),
        ("Complete coding candidate pool", "Retained by current default rules"),
    ):
        values = (
            counts.pivot(index="source_a", columns="source_b", values=field)
            .loc[keys, keys]
            .to_numpy()
        )
        ax.pcolormesh(
            np.arange(len(keys) + 1) - 0.5,
            np.arange(len(keys) + 1) - 0.5,
            values,
            cmap="Blues",
            vmin=0,
            vmax=max(1, values.max()),
        )
        ax.invert_yaxis()
        ax.set_aspect("equal")
        for (i, j), n in np.ndenumerate(values):
            ax.text(
                j,
                i,
                str(n),
                ha="center",
                va="center",
                fontsize=7.5,
                color="white" if n > values.max() * 0.55 else "#182c3e",
            )
        ax.set_xticks(range(len(keys)), keys, rotation=65, ha="right", fontsize=9)
        ax.set_yticks(range(len(keys)), keys, fontsize=9)
        ax.set_title(title, pad=18)
        ax.spines[:].set_visible(False)
    fig.text(
        0.5,
        0.02,
        "Cells count unique shared genes; diagonals are source totals. Nested and overlapping sources must not be summed.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1), w_pad=3)
    _save(fig, path, plt)


def _fig_publication_funnel(df, path, plt):
    from .cta_provenance import paper_intake_counts

    counts = paper_intake_counts()
    tags = list(counts.source_tag.unique())
    stages = [
        "published",
        "mapped",
        "protein_coding",
        "family_eligible",
        "hpa_restriction",
        "default_panel",
    ]
    values = (
        counts.pivot(index="source_tag", columns="stage", values="remaining")
        .loc[tags, stages]
        .to_numpy()
    )
    fractions = values / values[:, :1]
    fig, ax = plt.subplots(figsize=(12, max(6, len(tags) * 0.52)))
    ax.pcolormesh(
        np.arange(len(stages) + 1) - 0.5,
        np.arange(len(tags) + 1) - 0.5,
        fractions,
        cmap="Greens",
        vmin=0,
        vmax=1,
    )
    ax.invert_yaxis()
    for (i, j), n in np.ndenumerate(values):
        ax.text(
            j,
            i,
            f"{n:,}",
            ha="center",
            va="center",
            fontsize=11,
            color="white" if fractions[i, j] > 0.6 else "#182c3e",
        )
    ax.set_yticks(
        range(len(tags)),
        [{doi: label for label, doi in PRIMARY_SOURCES.items()}[tag] for tag in tags],
    )
    ax.set_xticks(
        range(6),
        [
            "Unique source\nidentities",
            "Unique IDs\nmapped",
            "Protein\ncoding",
            "Family\nfilter",
            "HPA\nrestriction",
            "Default\nrules",
        ],
    )
    ax.set_title("Ten selected papers: full lists through curation", pad=18)
    ax.set_xlabel(
        "Counts after each stage; shading is the fraction retained within that source", labelpad=12
    )
    ax.spines[:].set_visible(False)
    fig.text(
        0.5,
        0.015,
        "Source identities collapse mapped aliases and repeated entries across tables within each paper; unresolved probes remain upstream.\nPaper lists overlap and must not be summed. RNA/protein nomination does not establish HLA presentation.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    _save(fig, path, plt)


def _fig_placental_evidence_coverage(df, path, plt):
    from matplotlib.colors import ListedColormap

    from .cta_sources import placental_source_coverage

    data = placental_source_coverage().sort_values("Symbol")
    papers = {
        "Gong 2021": {"Gong2021_placenta_PC", "Gong2021_placenta_ncRNA"},
        "Bradley 2020": {"Bradley2020_CPA"},
        "Rull 2005": {"Rull2005_CGB_placenta"},
        "Rull 2008*": {"Rull2008_CGB_RNA"},
        "Kubiczak 2013*": {"Kubiczak2013_CGB_ovarian"},
        "Białas 2020": {"Bialas2020_CGB_cancer"},
        "McKellar 2025†": {"McKellar2025_CGB7_cancer"},
    }
    tags = data.publication_sources.str.split(";").map(set)
    values = np.array([[bool(t & p) for p in papers.values()] for t in tags], dtype=int)
    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    ax.pcolormesh(
        np.arange(len(papers) + 1) - 0.5,
        np.arange(len(data) + 1) - 0.5,
        values,
        cmap=ListedColormap(["#f0f0f0", "#398c70"]),
        vmin=0,
        vmax=1,
    )
    ax.invert_yaxis()
    for (i, j), present in np.ndenumerate(values):
        if present:
            ax.text(j, i, "●", color="white", ha="center", va="center", fontsize=11)
    ax.set_xticks(range(len(papers)), papers, rotation=35, ha="left")
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
    ax.set_yticks(
        range(len(data)), [r.Symbol + ("  ✓" if r.default_panel else "") for r in data.itertuples()]
    )
    ax.spines[:].set_visible(False)
    ax.set_title(
        "All 19 earlier placental nominations now have publication provenance", pad=88, fontsize=12
    )
    fig.text(
        0.5,
        0.025,
        "✓ Retained in default CTA panel (9/19). Blank cells mean no record in these curated source rows.\n"
        "* Combined CGB1/CGB2 assays, not separate positives. † Preprint.\n"
        "Expression / nomination evidence is distinct from gene-specific antigen validation.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.095, 1, 1))
    _save(fig, path, plt)


_BUILDERS = {
    "source_venn": _fig_source_venn,
    "stage_funnel": _fig_stage_funnel,
    "filter_funnel": _fig_filter_funnel,
    "filter_outcome": _fig_filter_outcome,
    "deflated_dist": _fig_deflated_dist,
    "protein_vs_rna": _fig_protein_vs_rna,
    "legacy_source_venn": _fig_legacy_source_venn,
    "landscape_source_venn": _fig_landscape_source_venn,
    "source_overlap": _fig_source_overlap,
    "placental_source_overlap": _fig_placental_source_overlap,
    "publication_funnel": _fig_publication_funnel,
    "placental_evidence_coverage": _fig_placental_evidence_coverage,
}


_FONT_SCALE = 1.0


def render(out_dir="cta_curation_out", *, kinds=None, font_scale=1.0) -> dict:
    """Write the CTA-curation figures into ``out_dir``.

    Returns ``{"n_genes": int, "stages": [...], "paths": {key: Path}}``.
    """
    global _FONT_SCALE
    if not np.isfinite(font_scale) or font_scale <= 0:
        raise ValueError("font_scale must be positive and finite")
    _FONT_SCALE = font_scale
    figure_style.apply()
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = _evidence()
    paths = {}
    import matplotlib

    from . import cta
    from .cta_sources import (
        intake_counts,
        intake_membership,
        publication_membership,
        publication_sources,
    )

    with matplotlib.rc_context({"font.size": 10}):
        for key in kinds if kinds is not None else _BUILDERS:
            path = out / FILENAMES[key]
            _BUILDERS[key](df, path, plt)
            paths[key] = path
    from .cta_provenance import (
        candidate_provenance,
        gene_citation_evidence_report,
        historical_tag_audit,
        legacy_only_candidates,
        paper_intake_counts,
        selected_membership,
        source_cover,
    )

    (out / "cta-minimum-source-cover.json").write_text(json.dumps(source_cover(), indent=2) + "\n")
    candidate_provenance().to_csv(out / "cta-candidate-provenance.csv", index=False)
    gene_citation_evidence_report().to_csv(out / "cta-gene-citation-evidence.csv", index=False)
    pd.DataFrame(_per_source_counts(df)).to_csv(out / "cta-source-outcome-counts.csv", index=False)
    historical_tag_audit().to_csv(out / "cta-historical-tag-audit.csv", index=False)
    paper_intake_counts().to_csv(out / "cta-selected-paper-intake-counts.csv", index=False)
    from .cta_sources import gene_publication_evidence

    gene_publication_evidence().to_csv(out / "cta-gene-publication-evidence.csv", index=False)
    from .cta_sources import placental_source_coverage

    placental_source_coverage().to_csv(out / "placental-nomination-provenance.csv", index=False)
    selected_membership().to_csv(out / "cta-selected-paper-membership.csv", index=False)
    legacy_only_candidates().to_csv(out / "cta-legacy-only-candidates.csv", index=False)
    stages = stage_counts()
    pd.DataFrame(stages, columns=["stage", "remaining", "dropped"]).to_csv(
        out / "cta-stage-counts.csv", index=False
    )
    stage_membership(df).to_csv(out / "cta-stage-membership.csv", index=False)
    intake_counts().to_csv(out / "cta-publication-intake-counts.csv", index=False)
    intake_membership().to_csv(out / "cta-publication-intake-membership.csv", index=False)
    publication_sources().to_csv(out / "cta-publication-sources.csv", index=False)
    source_overlap_counts(df).to_csv(out / "cta-source-overlap-counts.csv", index=False)
    df.to_csv(out / "cta-candidate-evidence.csv", index=False)
    default = cta.cta_gene_ids()
    df[df.Ensembl_Gene_ID.isin(default)].to_csv(out / "cta-default-panel.csv", index=False)
    sets = _tag_sets(df)
    pd.DataFrame(
        [
            {"source": source, "Ensembl_Gene_ID": gid, "default_panel": gid in default}
            for source, ids in sets.items()
            for gid in sorted(ids)
        ]
    ).to_csv(out / "cta-source-membership.csv", index=False)
    counts = {}
    for gid in set().union(*sets.values()):
        signature = tuple(name for name, ids in sets.items() if gid in ids)
        counts[signature] = counts.get(signature, 0) + 1
    pd.DataFrame(
        [
            {"sources": ";".join(key), "genes": n}
            for key, n in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        ]
    ).to_csv(out / "cta-exact-intersections.csv", index=False)
    inputs = {
        name: hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()
        for name, frame in [
            ("candidate_table", df),
            ("publication_membership", publication_membership()),
            ("publication_sources", publication_sources()),
        ]
    }
    manifest = {
        "oncoref_version": __version__,
        "oncoref_data_version": DATA_VERSION,
        "input_sha256": {
            name: hashlib.sha256(
                files("oncoref").joinpath("data").joinpath(name).read_bytes()
            ).hexdigest()
            for name in (
                "cancer-testis-antigens.csv",
                "cta-specificity-audit.csv",
                "cta-publication-sources.csv",
                "cta-publication-membership.csv",
                "cta-gene-publication-evidence.csv",
            )
        },
        "default_gene_ids": sorted(default),
        "default_panel_sha256": hashlib.sha256(
            ("\n".join(sorted(default)) + "\n").encode()
        ).hexdigest(),
        "figure_kinds": list(paths),
        "font_scale": font_scale,
        "inputs": inputs,
        "candidate_genes": len(df),
        "default_genes": len(default),
        "stages": stages,
        "minimum_source_cover": source_cover(),
        "outputs": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(out.iterdir())
            if p.is_file()
            and p.suffix in {".png", ".pdf", ".csv", ".json"}
            and p.name != "run-manifest.json"
        },
    }
    (out / "run-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = [
        "# CTA source intake and curation",
        "",
        f"{len(df):,} coding candidates; {len(default):,} pass the default rules, including the normal-reproductive-only nomination holdback.",
        "",
        "Sources are upstream nominations, not antigen-validation claims. Complete memberships and exact intersections accompany every figure.",
        "",
    ]
    for key, path in paths.items():
        lines.extend(
            [
                f"## {key.replace('_', ' ').title()}",
                "",
                f"[Vector PDF]({path.with_suffix('.pdf').name})",
                "",
                f"![{key}]({path.name})",
                "",
            ]
        )
    (out / "index.md").write_text("\n".join(lines))
    return {"n_genes": len(df), "stages": stages, "paths": paths}
