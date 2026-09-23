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

Six figures describing how the CTA panel is defined: where the candidate genes
come from (source venn), how many survive each curation stage (stage funnel),
how that attrition splits by source (per-source funnel and outcome), and the two
evidence axes the filter actually thresholds on (deflated reproductive-fraction
distribution, protein-reliability-vs-RNA).

Two frames are in play and they are not interchangeable:

``_evidence()``
    The raw packaged table — every candidate ever considered, including rows that
    fail curation. This is the denominator the source and filter figures describe.

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

from pathlib import Path

import numpy as np

from . import figure_style
from .cta_tissues import HPA_ADAPTIVE_PROTEIN_RNA_THRESHOLDS
from .figure_style import ACCENT, DROP, KEPT, THRESHOLD, WEAK
from .load_dataset import get_data

# Primary gene-contributing sources. The broad ``daSilva2017`` cross-reference tag
# (the full 1,103-gene 0.9-threshold set) and the tiny paralog/candidate tags are
# deliberately excluded: only ``daSilva2017_protein``, the mass-spec-validated
# subset, is a primary source. Both tags coexist in the packaged table, so the
# predicate has to test for the protein suffix specifically.
PRIMARY_SOURCES = {
    "Gong 2021": lambda tags: bool(tags & {"Gong2021_placenta_PC", "Gong2021_placenta_ncRNA"}),
    "Bradley 2020": lambda tags: "Bradley2020_CPA" in tags,
    "CTpedia": lambda tags: "CTpedia" in tags,
    "CTexploreR": lambda tags: "CTexploreR_CT" in tags or "CTexploreR_CTP" in tags,
    "daSilva2017_protein": lambda tags: "daSilva2017_protein" in tags,
    "placental_antigen": lambda tags: "placental_antigen" in tags,
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
}


def _evidence():
    """Raw packaged CTA evidence table, including rows that fail curation filters."""
    return get_data("cancer-testis-antigens").copy()


def _curated():
    """CTA table with non-CTA exclusions dropped and specificity decisions joined."""
    from . import cta

    return cta._cta_with_specificity_frame().copy()


def _bool_series(series):
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def _tag_sets(df):
    """Return ``{source_label: set(Ensembl_Gene_ID)}`` for primary sources."""
    out = {name: set() for name in PRIMARY_SOURCES}
    for ensg, raw in zip(df["Ensembl_Gene_ID"], df["source_databases"].fillna("")):
        tags = {t.strip() for t in str(raw).split(";") if t.strip()}
        for name, pred in PRIMARY_SOURCES.items():
            if pred(tags):
                out[name].add(str(ensg))
    return out


def _per_source_counts(df):
    sets = _tag_sets(df)
    rows = []
    for name, members in sets.items():
        sub = df[df["Ensembl_Gene_ID"].astype(str).isin(members)]
        passes = _bool_series(sub["passes_filters"])
        weak = _bool_series(sub["never_expressed"])
        rows.append(
            {
                "source": name,
                "total": len(sub),
                "kept_confident": int((passes & ~weak).sum()),
                "kept_weak": int((passes & weak).sum()),
                "excluded": int((~passes).sum()),
            }
        )
    return sorted(rows, key=lambda r: r["total"], reverse=True)


def stage_counts():
    """Sequential curation stages, as ``[(label, n_remaining, n_dropped), ...]``.

    The four stages a candidate has to survive, in the order the pipeline applies
    them: the raw source union, removal of genes reclassified as non-CTA, the
    tiered HPA protein/RNA restriction filter, and the specificity adjudication
    that produces the shipped default set.
    """
    raw = _evidence()
    cur = _curated()
    passes = _bool_series(cur["passes_filters"])
    action = cur["specificity_action"].astype(str) if "specificity_action" in cur else None
    n_source = len(raw)
    n_kept_cta = len(cur)
    n_passes = int(passes.sum())
    n_default = int((action == "include_default").sum()) if action is not None else n_passes
    return [
        ("source union", n_source, 0),
        ("non-CTA removed", n_kept_cta, n_source - n_kept_cta),
        ("HPA restriction", n_passes, n_kept_cta - n_passes),
        ("specificity", n_default, n_passes - n_default),
    ]


def _save(fig, path, plt):
    figure_style.save(fig, path)
    plt.close(fig)


def _fig_source_venn(df, path, plt):
    from matplotlib_venn import venn3

    sets = _tag_sets(df)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    keys = ("CTpedia", "CTexploreR", "daSilva2017_protein")
    v = venn3(
        [sets[k] for k in keys],
        set_labels=("CTpedia", "CTexploreR", "da Silva 2017"),
        set_colors=(KEPT, ACCENT, WEAK),
        alpha=0.55,
        ax=ax,
    )
    for text in list(v.set_labels or []) + list(v.subset_labels or []):
        if text is not None:
            text.set_fontsize(10)
    _save(fig, path, plt)


def _fig_stage_funnel(df, path, plt):
    """Sequential attrition through the four curation stages."""
    stages = stage_counts()
    labels = [s[0] for s in stages]
    remaining = np.array([s[1] for s in stages])
    dropped = np.array([s[2] for s in stages])
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.barh(y, remaining, color=KEPT, height=0.62)
    ax.barh(y, dropped, left=remaining, color=DROP, height=0.62)
    for i, (n, d) in enumerate(zip(remaining, dropped)):
        ax.text(n - max(remaining) * 0.015, i, str(n), va="center", ha="right", color="white")
        if d:
            ax.text(n + d + max(remaining) * 0.012, i, f"−{d}", va="center", fontsize=9)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("genes")
    ax.set_xlim(0, max(remaining) * 1.1)
    _save(fig, path, plt)


def _fig_filter_funnel(df, path, plt):
    rows = _per_source_counts(df)
    labels = [r["source"] for r in rows]
    kept = np.array([r["kept_confident"] + r["kept_weak"] for r in rows])
    dropped = np.array([r["excluded"] for r in rows])
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.62 * len(labels) + 1.4))
    ax.barh(y, kept, color=KEPT, label="kept", height=0.62)
    ax.barh(y, dropped, left=kept, color=DROP, label="excluded", height=0.62)
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
    ax.set_xlim(0, max(kept + dropped) * 1.14)
    ax.legend(loc="lower right")
    _save(fig, path, plt)


def _fig_filter_outcome(df, path, plt):
    rows = _per_source_counts(df)
    labels = [r["source"] for r in rows]
    conf = np.array([r["kept_confident"] for r in rows])
    weak = np.array([r["kept_weak"] for r in rows])
    excl = np.array([r["excluded"] for r in rows])
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7, 0.62 * len(labels) + 1.4))
    ax.barh(y, conf, color=KEPT, label="kept", height=0.62)
    ax.barh(y, weak, left=conf, color=WEAK, label="kept, weak evidence", height=0.62)
    ax.barh(y, excl, left=conf + weak, color=DROP, label="excluded", height=0.62)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("genes")
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
        label=["kept", "excluded"],
    )
    for thr in sorted(set(RELIABILITY_THRESHOLD.values())):
        ax.axvline(thr, color=THRESHOLD, ls="--", lw=0.7, alpha=0.7)
    ax.set_xlabel("deflated reproductive fraction")
    ax.set_ylabel("genes")
    ax.legend(loc="upper left")
    _save(fig, path, plt)


def _fig_protein_vs_rna(df, path, plt):
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
    _save(fig, path, plt)


_BUILDERS = {
    "source_venn": _fig_source_venn,
    "stage_funnel": _fig_stage_funnel,
    "filter_funnel": _fig_filter_funnel,
    "filter_outcome": _fig_filter_outcome,
    "deflated_dist": _fig_deflated_dist,
    "protein_vs_rna": _fig_protein_vs_rna,
}


def render(out_dir="cta_curation_out") -> dict:
    """Write the CTA-curation figures into ``out_dir``.

    Returns ``{"n_genes": int, "stages": [...], "paths": {key: Path}}``.
    """
    figure_style.apply()
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = _evidence()
    paths = {}
    for key, builder in _BUILDERS.items():
        path = out / FILENAMES[key]
        builder(df, path, plt)
        paths[key] = path
    return {"n_genes": len(df), "stages": stage_counts(), "paths": paths}
