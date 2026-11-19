"""Survey cross-tissue RNA outliers across the 52 primary GTEx v11 tissue matrices."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from cta_report_common import sha256, write_csv
from gtex_sample_evidence import SOURCE, sample_matrix
from scipy.stats import rankdata

from oncoref import reference_data

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "docs/audits/target-selection-20261110/rebuild/gtex-library-qc.csv.gz"
OUT = ROOT / "docs/audits/gtex-cross-tissue-signals"
THRESHOLDS = (0.001, 0.01, 0.03, 0.05, 0.1)
LABELS = {"Small": "Small intestine", "Whole": "Whole blood", "Minor": "Salivary gland"}


def organ(site):
    prefix = site.split("_")[0]
    return LABELS.get(prefix, prefix)


def normalized_ranks(a):
    ranks = rankdata(a, axis=0, method="average")
    ranks -= ranks.mean(axis=0)
    norms = np.linalg.norm(ranks, axis=0)
    return ranks / np.where(norms > 0, norms, np.nan)


def source_panels(sites):
    median_path = reference_data.pinned_path("gtex_tissue_rna", "v11")
    med = pd.read_csv(median_path, sep="\t", skiprows=2).set_index("Name")[sites]
    abundance = med.T.groupby(organ).max().T
    panels = {}
    for src in abundance:
        highest_other = abundance.drop(columns=src).max(axis=1).clip(lower=1)
        ids = abundance.index[abundance[src].ge(20 * highest_other)].sort_values()
        if len(ids) >= 10:
            panels[src] = ids
    return abundance, panels


def survey(raw, annotations, out):
    q = pd.read_csv(INVENTORY, keep_default_na=False).set_index("sample_id")
    if sha256(annotations) != "98c246b3855edf62cb050982a7ca50f3069d010624baea47ef84a9465900ca4b":
        raise ValueError("Unexpected v11 annotation file")
    meta = pd.read_csv(annotations, sep="\t", low_memory=False).set_index("SAMPID").loc[q.index]
    sites = sorted(q.source_tissue.unique())
    abundance, panels = source_panels(sites)
    wanted = {g.encode(): g for ids in panels.values() for g in ids}
    pieces = []
    for site in sites:
        f = raw / (site.lower() + ".gct.gz")
        if (
            sha256(f)
            != reference_data.REFERENCE_SOURCES[SOURCE]["pins"]["v11-" + site.lower()]["sha256"]
        ):
            raise ValueError(f"Changed GTEx shard {site}")
        m = sample_matrix(f, wanted)[0]
        if set(m.columns) != set(q.index[q.source_tissue.eq(site)]) or len(m) != len(wanted):
            raise ValueError(f"Incomplete tissue matrix {site}")
        pieces.append(m)
    data = pd.concat(pieces, axis=1)
    if data.isna().any().any():
        raise ValueError("Missing GTEx panel measurement")
    rows = []
    for src, ids in panels.items():
        scale = abundance.loc[ids, src]
        native = q.index[q.source_tissue.map(organ).eq(src)]
        donors = data.loc[ids, native].div(scale, axis=0).to_numpy()
        split = [np.arange(len(ids)), np.arange(0, len(ids), 2), np.arange(1, len(ids), 2)]
        ranked = [normalized_ranks(donors[idx, :]) for idx in split] if len(ids) >= 20 else []
        for site in sites:
            if organ(site) == src:
                continue
            obs = data.loc[ids, q.index[q.source_tissue.eq(site)]]
            excess = obs.sub(obs.median(axis=1), axis=0).clip(lower=0).div(scale, axis=0)
            score = excess.median(axis=0)
            breadth = obs.gt(obs.quantile(0.95, axis=1), axis=0).mean(axis=0)
            for sample in score.index[(score >= min(THRESHOLDS)) & (breadth >= 0.75)]:
                rec = {
                    "sample_id": sample,
                    "tissue": site,
                    "source_panel": src,
                    "markers": len(ids),
                    "score": score[sample],
                    "breadth": breadth[sample],
                    "best_source": "",
                    "rho": np.nan,
                    "split_match": False,
                    "shared_extraction": False,
                    "shared_expression": False,
                    "different_donor": False,
                    "probable": False,
                }
                if score[sample] >= 0.01 and ranked:
                    matches = []
                    for idx, ref in zip(split, ranked, strict=True):
                        corr = normalized_ranks(excess[sample].to_numpy()[idx]) @ ref
                        matches.append(corr)
                    if all(np.isfinite(r).any() for r in matches):
                        i = np.nanargmax(matches[0])
                        best = native[i]
                        rec.update(
                            best_source=best,
                            rho=matches[0][i],
                            split_match=all(
                                np.flatnonzero(r == np.nanmax(r)).tolist() == [i] for r in matches
                            ),
                            different_donor=q.at[best, "participant_id"]
                            != q.at[sample, "participant_id"],
                        )
                        for field, col in [
                            ("shared_extraction", "SMNABTCH"),
                            ("shared_expression", "SMGEBTCH"),
                        ]:
                            a, b = meta.at[best, col], meta.at[sample, col]
                            rec[field] = bool(
                                pd.notna(a)
                                and pd.notna(b)
                                and set(a.split(",")) & set(b.split(","))
                            )
                        rec["probable"] = (
                            rec["different_donor"]
                            and rec["split_match"]
                            and rec["rho"] >= 0.5
                            and (rec["shared_extraction"] or rec["shared_expression"])
                        )
                rows.append(rec)
    flags = pd.DataFrame(rows).sort_values(["score", "sample_id"], ascending=[False, True])
    out.mkdir(parents=True, exist_ok=True)
    write_csv(flags, out / "signals.csv.gz")
    matched = flags.loc[flags.probable].copy()
    matched["recipient_extraction_batch"] = matched.sample_id.map(meta.SMNABTCH)
    matched["source_extraction_batch"] = matched.best_source.map(meta.SMNABTCH)
    matched["recipient_expression_batch"] = matched.sample_id.map(meta.SMGEBTCH)
    matched["source_expression_batch"] = matched.best_source.map(meta.SMGEBTCH)
    write_csv(matched, out / "batch-matched.csv")
    raw_notes = q.loc[q.source_tissue.ne("Testis"), "testis_tissue_content"]
    counts = []
    for th in THRESHOLDS:
        part = flags[flags.score >= th]
        counts.append(
            {
                "threshold": th,
                "screened_libraries": len(q),
                "any_panel": part.sample_id.nunique(),
                "broad_testis_panel": part[part.source_panel.eq("Testis")].sample_id.nunique(),
                "six_testis_markers": int(raw_notes.ge(th).sum()),
            }
        )
    write_csv(pd.DataFrame(counts), out / "thresholds.csv")
    write_csv(
        pd.DataFrame(
            [(src, len(ids)) for src, ids in panels.items()], columns=["panel", "markers"]
        ),
        out / "panels.csv",
    )
    print(pd.DataFrame(counts).to_string(index=False))
    print("Probable libraries", flags.loc[flags.probable, "sample_id"].nunique())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("annotations", type=Path)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    survey(args.raw, args.annotations, args.out)
