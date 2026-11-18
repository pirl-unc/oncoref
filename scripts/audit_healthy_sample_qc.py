"""Cross-tissue HPA biopsy QC and GTEx donor sensitivity for the healthy TPM ceiling."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from antigen_sample_evidence import HPA_SOURCE, tissue_sample_genes
from audit_gtex_cross_tissue_signals import INVENTORY, OUT, ROOT, source_panels
from cta_report_common import sha256, write_csv
from gtex_sample_evidence import SOURCE, SOURCE_ID, sample_matrix

from oncoref import reference_data

REBUILD = INVENTORY.parent
BASE = ROOT / "docs/audits/cta-selection-20261028"
HPA_NATIVE = {
    "Adrenal": ["adrenal gland"],
    "Brain": ["cerebral cortex", "choroid plexus"],
    "Heart": ["heart muscle"],
    "Kidney": ["kidney"],
    "Liver": ["liver"],
    "Muscle": ["skeletal muscle"],
    "Pancreas": ["pancreas"],
    "Salivary gland": ["salivary gland"],
    "Skin": ["skin"],
    "Small intestine": ["small intestine", "duodenum"],
    "Stomach": ["stomach"],
    "Testis": ["testis"],
    "Thyroid": ["thyroid gland"],
}


def hpa_screen(abundance, panels, q):
    genes = tissue_sample_genes()
    genes.columns = genes.columns.str.removeprefix("TPM.")
    if set(genes.columns) != set(q.index):
        raise ValueError("Incomplete HPA biopsy inventory")
    decisions = pd.read_csv(BASE / "alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz")
    med = pd.read_csv(
        ROOT / "tmp/cta-source-consolidation/gtex-v11-median-tpm.gct.gz",
        sep="\t",
        skiprows=2,
        usecols=["Name"],
    )
    chosen = decisions.dropna(subset=["row_position"])
    lookup = pd.Series(
        chosen.gene_id.to_numpy(), index=med.Name.iloc[chosen.row_position.astype(int)]
    )
    hpa_ids = pd.read_csv(BASE / "hpa_rna_consensus-accepted-identities.csv.gz")
    genes.index = genes.index.map(hpa_ids.set_index("source_gene_id").canonical_gene_id)
    genes = genes.loc[genes.index.notna()]
    if not genes.index.is_unique:
        raise ValueError("Duplicate mapped HPA genes")
    hits, coverage = [], []
    for panel, ids in panels.items():
        mapped = lookup.reindex(ids).dropna()
        mapped = mapped[mapped.isin(genes.index)]
        values = genes.loc[mapped].copy()
        values.index = mapped.index
        values = values.dropna()
        coverage.append((panel, len(ids), len(values)))
        if len(values) < 10 or len(values) < 0.5 * len(ids):
            continue
        # External scale preserves a common definition even for organs absent in HPA.
        scale = abundance.loc[values.index, panel]
        for tissue, group in q.groupby("tissue"):
            if tissue in HPA_NATIVE.get(panel, []) or len(group) < 3:
                continue
            for sample in group.index:
                peers = values.loc[:, group.index.drop(sample)]
                signal = values[sample]
                excess = signal.sub(peers.median(axis=1)).clip(lower=0).div(scale)
                score = excess.median()
                breadth = signal.gt(peers.max(axis=1)).mean()
                if score >= 0.01 and breadth >= 0.75:
                    hits.append(
                        (sample, tissue, panel, len(values), len(peers.columns), score, breadth)
                    )
    flags = pd.DataFrame(
        hits, columns=["sample_id", "tissue", "panel", "markers", "peers", "score", "breadth"]
    )
    return (
        flags.sort_values("score", ascending=False),
        pd.DataFrame(coverage, columns=["panel", "gtex_markers", "hpa_markers"]),
    )


def donor_statistics(raw, q, normal, catalog, probable, flagged):
    ids = pd.read_csv(BASE / "alternative-evidence/gtex_tissue_rna-row-decisions.csv.gz")
    ids = ids.dropna(subset=["row_position"])
    med = pd.read_csv(
        ROOT / "tmp/cta-source-consolidation/gtex-v11-median-tpm.gct.gz",
        sep="\t",
        skiprows=2,
        usecols=["Name"],
    )
    groups = {g: r.entity_id for r in catalog.itertuples() for g in json.loads(r.member_gene_ids)}
    selected = ids.loc[ids.gene_id.isin(groups)]
    canonical = dict(
        zip(
            med.Name.iloc[selected.row_position.astype(int)].str.encode("utf-8"),
            selected.gene_id,
            strict=True,
        )
    )
    sizes = catalog.member_gene_ids.map(lambda s: len(json.loads(s)))
    parts, maxima = [], []
    for site, table in normal[normal.source_id.eq(SOURCE_ID)].groupby("source_context"):
        file = raw / (site.lower() + ".gct.gz")
        if (
            sha256(file)
            != reference_data.REFERENCE_SOURCES[SOURCE]["pins"]["v11-" + site.lower()]["sha256"]
        ):
            raise ValueError(f"Changed GTEx shard {site}")
        genes, donors = sample_matrix(file, canonical)
        assigned = q[q.source_tissue.eq(site)]
        if (
            set(assigned.index) != set(genes.columns)
            or not donors.is_unique
            or not np.array_equal(donors, assigned.loc[genes.columns, "participant_id"])
        ):
            raise ValueError(f"Unexpected donor assignments in {site}")
        grouped = genes.groupby(genes.index.map(groups))
        vals = grouped.sum(min_count=1).where(grouped.count().eq(sizes, axis=0)).dropna(how="all")
        ref = table.set_index("entity_id").value.sort_index()
        if not ref.index.equals(vals.index.sort_values()) or not np.allclose(
            vals.max(axis=1).loc[ref.index], ref
        ):
            raise ValueError(f"Cannot reproduce original sample maxima in {site}")
        clean = vals.loc[:, ~vals.columns.isin(probable)]
        v = vals.to_numpy()
        if not np.isfinite(v).all() or clean.empty:
            raise ValueError(f"Incomplete readings in {site}")
        stats = pd.DataFrame(
            {
                "max": np.max(v, axis=1),
                "without_M": clean.max(axis=1),
                "second_donor": np.partition(v, -2, axis=1)[:, -2],
                "p95": np.quantile(v, 0.95, axis=1),
                "p75": np.quantile(v, 0.75, axis=1),
            },
            index=vals.index,
        )
        parts.append(stats)
        peaked = table.loc[table.sample_id.isin(flagged)].copy()
        if not peaked.empty:
            idx = peaked.entity_id
            peaked["site_donors"] = len(donors)
            peaked["second_donor_tpm"] = stats.second_donor.reindex(idx).to_numpy()
            peaked["site_p95_tpm"] = stats.p95.reindex(idx).to_numpy()
            peaked["site_without_M_tpm"] = stats.without_M.reindex(idx).to_numpy()
            maxima.append(peaked)
    return pd.concat(parts).groupby(level=0).max(), pd.concat(maxima)


def run(raw, out):
    q = pd.read_csv(INVENTORY).set_index("sample_id")
    h = pd.read_csv(REBUILD / "normal-library-qc.csv").set_index("sample_id")
    abundance, panels = source_panels(sorted(q.source_tissue.unique()))
    hflags, coverage = hpa_screen(abundance, panels, h)
    write_csv(hflags, out / "hpa-signals.csv")
    write_csv(coverage, out / "hpa-panels.csv")
    gcounts = q.groupby("source_tissue").agg(
        samples=("participant_id", "size"), known_donors=("participant_id", "nunique")
    )
    hcounts = h.groupby("tissue").size().to_frame("samples")
    counts = pd.concat({"GTEx": gcounts, "HPA": hcounts}, names=["source", "site"])
    write_csv(counts.astype({"known_donors": "Int64"}).reset_index(), out / "sample-coverage.csv")
    normal = pd.read_csv(
        REBUILD / "normal-comparison.csv.gz",
        usecols=[
            "entity_id",
            "source_id",
            "scope",
            "sample_id",
            "source_context",
            "qc_pass",
            "identity_resolved",
            "value",
        ],
    )
    normal = normal[
        normal.qc_pass & normal.identity_resolved & normal.scope.isin(["normal", "disposable"])
    ]
    cat = pd.read_csv(ROOT / "oncoref/data/target-selection.csv.gz").set_index(
        "entity_id", drop=False
    )
    cat = cat.loc[sorted(normal.entity_id.unique())]
    signals = pd.read_csv(out / "signals.csv.gz")
    matched = set(signals.loc[signals.probable, "sample_id"])
    flags = set(hflags.sample_id) | set(signals.loc[signals.score.ge(0.01), "sample_id"])
    stats, maxima = donor_statistics(raw, q, normal, cat, matched, flags)
    fixed = (
        normal.loc[normal.source_id.ne(SOURCE_ID)]
        .groupby("entity_id")
        .value.max()
        .reindex(cat.index)
    )
    ceiling = stats.reindex(cat.index).apply(lambda s: pd.concat([s, fixed], axis=1).max(axis=1))
    actual = cat[["normal_reference_tpm", "disposable_max_tpm"]].max(axis=1)
    if not np.allclose(ceiling["max"], actual):
        raise ValueError("Current catalog ceiling did not reproduce")
    # Also vary HPA independently within each biopsy tissue, retaining maxima in
    # tissues with fewer than three observations. Anonymous biopsies cannot be deduplicated.
    individual_sources = {SOURCE_ID, HPA_SOURCE["source_id"]}
    other = normal.loc[~normal.source_id.isin(individual_sources)].groupby("entity_id").value.max()
    biopsies = normal.loc[normal.source_id.eq(HPA_SOURCE["source_id"])].copy()
    biopsies["site"] = biopsies.sample_id.map(h.tissue)
    if biopsies.site.isna().any():
        raise ValueError("Unknown HPA biopsy")
    tissue_values = biopsies.groupby(["entity_id", "site"]).value
    for label, quantile in [("p95", 0.95), ("p75", 0.75)]:
        hp = tissue_values.quantile(quantile).where(tissue_values.count() >= 3, tissue_values.max())
        hp = hp.groupby("entity_id").max()
        ceiling["both_" + label] = pd.concat([stats[label], hp, other], axis=1).max(axis=1)
    thresholds = ceiling.clip(lower=1) * 10
    passes = thresholds.le(cat.tumor_max_tpm, axis=0)
    rows = []
    for label in thresholds:
        rows.append(
            (
                label,
                len(cat),
                int(passes[label].sum()),
                int((passes[label] & ~passes["max"]).sum()),
                int((thresholds[label] < thresholds["max"] - 1e-9).sum()),
            )
        )
    write_csv(
        pd.DataFrame(
            rows,
            columns=["reference_statistic", "proteins", "tumor_10x", "gained", "lower_threshold"],
        ),
        out / "normal-sensitivity.csv",
    )
    detail = cat[["protein_group_label", "classification", "tumor_max_tpm"]].join(
        ceiling.add_prefix("ceiling_")
    )
    write_csv(detail.reset_index(), out / "protein-sensitivity.csv.gz")
    top = normal[
        normal.sample_id.isin(flags)
        & normal.value.eq(normal.entity_id.map(actual))
        & normal.value.gt(0)
    ].copy()
    top["flag"] = np.where(top.sample_id.isin(matched), "M", "O")
    top["protein"] = top.entity_id.map(cat.protein_group_label)
    top["ceiling_without_M"] = top.entity_id.map(ceiling.without_M)
    top = top.merge(
        maxima[
            [
                "entity_id",
                "sample_id",
                "site_donors",
                "second_donor_tpm",
                "site_p95_tpm",
                "site_without_M_tpm",
            ]
        ],
        on=["entity_id", "sample_id"],
        how="left",
        validate="many_to_one",
    )
    write_csv(
        top[
            [
                "protein",
                "source_id",
                "sample_id",
                "source_context",
                "flag",
                "value",
                "ceiling_without_M",
                "site_donors",
                "second_donor_tpm",
                "site_p95_tpm",
                "site_without_M_tpm",
            ]
        ].sort_values(["protein", "sample_id"]),
        out / "flagged-maxima.csv",
    )
    print(pd.DataFrame(rows).to_string(index=False, header=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    run(args.raw, args.out)
