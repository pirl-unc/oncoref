"""Native TPM references computed within a normal anatomical site or sorted population."""

import json

import pandas as pd

from oncoref import reference_data, target_selection


def tissue_percentiles(normal, *, method="p95"):
    mask = normal.source_id.eq("hpa_tissue_samples:v23") & normal.scope.isin(
        ["normal", "disposable", "thymic"]
    )
    raw = normal.loc[mask].copy()
    if raw.empty:
        return normal
    keys = ["entity_id", "tissue"]
    raw["value"] = raw.value.where(raw.identity_resolved & raw.qc_pass)
    if raw.duplicated([*keys, "sample_id"]).any():
        raise ValueError("Duplicate HPA tissue biopsy for a protein")
    measured = raw.groupby(keys).value
    summary = raw.drop_duplicates(keys).set_index(keys)
    if method == "mean_plus_2std":
        summary["value"] = measured.apply(
            lambda s: target_selection.summarize_normal_tissue_tpm(s, method=method)
        )
        summary["statistic"] = "within-tissue biopsy mean + 2 SD"
        summary["observation_kind"] = "biopsy_mean_plus_2std"
    else:
        summary["value"] = measured.quantile(0.95)
        summary["statistic"] = "within-tissue biopsy P95, linear interpolation"
        summary["observation_kind"] = "biopsy_percentile"
    summary["n_observations"] = measured.count()
    summary["source_context"] = summary.index.get_level_values("tissue")
    for field in ["sample_id", "participant_id", "library_qc_review", "mammary_state_inference"]:
        summary[field] = ""
    summary["qc_pass"] = summary.value.notna()
    summary["identity_resolved"] = summary.value.notna()
    return pd.concat([normal.loc[~mask], summary.reset_index()], ignore_index=True)


def _members(facts):
    return {g: r.entity_id for r in facts.itertuples() for g in json.loads(r.member_gene_ids)}


def _population(frame):
    return frame.assign(
        tissue="blood",
        context_id="blood",
        cell_type=frame.source_context,
        scope="normal",
        qc_pass=frame.value.notna(),
        identity_resolved=frame.value.notna(),
        unit="TPM",
    )


def _input(name, receipts):
    version = reference_data.resolve_version(name)
    receipts[f"reference:{name}:{version}"] = reference_data.REFERENCE_SOURCES[name]["pins"][
        version
    ]["sha256"]
    return reference_data.ensure(name, version)


def dice_percentiles(facts, receipts, *, method="p95"):
    members = _members(facts)
    sizes = facts.set_index("entity_id").member_count
    rows = []
    for name, spec in reference_data.REFERENCE_SOURCES.items():
        if not name.startswith("dice_"):
            continue
        frame = pd.read_csv(_input(name, receipts))
        ids = frame.Feature_name.str.split(".").str[0]
        valid = ids.isin(members) & ~ids.duplicated(keep=False)
        values = frame.loc[valid, [c for c in frame if c.isdigit()]]
        values.index = ids[valid].map(members)
        single = values.loc[sizes.reindex(values.index).eq(1).to_numpy()]
        # Rank columns are not specimen identifiers across genes. The sum of
        # measured locus maxima bounds any possible within-specimen group sum.
        multi = values.max(axis=1).groupby(level=0).agg(["sum", "count"])
        multi = multi.loc[sizes.reindex(multi.index).gt(1)]
        multi["value"] = multi["sum"].where(multi["count"].eq(sizes.reindex(multi.index)))
        stat_values = target_selection.summarize_normal_tissue_tpm(single.T, method=method)
        stat_label = (
            "within-population mean + 2 SD"
            if method == "mean_plus_2std"
            else "within-population P95"
        )
        estimates = pd.DataFrame({"value": stat_values, "n_observations": single.count(axis=1)})
        result = pd.concat(
            [
                estimates.assign(statistic=stat_label),
                multi[["value"]].assign(
                    statistic="sum of member maxima; unpaired rank distributions"
                ),
            ]
        )
        rows.append(
            result.rename_axis("entity_id")
            .reset_index()
            .assign(source_id=name + ":" + spec["version"], source_context=spec["cell_type"])
        )
    return _population(pd.concat(rows, ignore_index=True))


def hpa_immune_percentiles(facts, accepted_ids, receipts, *, method="p95"):
    name = "hpa_immune_cell_samples"
    members = _members(facts)
    lookup = accepted_ids.set_index("source_gene_id").canonical_gene_id
    sizes = facts.set_index("entity_id").member_count
    parts = []
    for chunk in pd.read_csv(
        _input(name, receipts),
        sep="\t",
        usecols=["ENSG ID", "Donor", "Sample ID", "Immune cell", "TPM"],
        chunksize=200000,
    ):
        chunk["entity_id"] = chunk["ENSG ID"].map(lookup).map(members)
        parts.append(chunk.dropna(subset=["entity_id"]))
    raw = pd.concat(parts)
    keys = ["entity_id", "Immune cell", "Donor", "Sample ID"]
    if raw.duplicated([*keys, "ENSG ID"]).any() or raw[keys].isna().any().any():
        raise ValueError("Duplicate locus or missing HPA sorted sample identity")
    grouped = raw.groupby(keys).TPM.agg(["sum", "count"])
    v = grouped["sum"].where(
        grouped["count"].eq(grouped.index.get_level_values("entity_id").map(sizes))
    )
    donors = v.groupby(level=keys[:3]).median().groupby(level=keys[:2])
    if method == "mean_plus_2std":
        agg_val = donors.apply(
            lambda s: target_selection.summarize_normal_tissue_tpm(s, method=method)
        )
        stat_label = "within-population donor mean + 2 SD"
    else:
        agg_val = donors.quantile(0.95)
        stat_label = "within-population donor P95"
    frame = (
        pd.DataFrame({"value": agg_val, "n_observations": donors.count()})
        .reset_index()
        .rename(columns={"Immune cell": "source_context"})
    )
    return _population(frame.assign(source_id=name + ":v23", statistic=stat_label))
