"""Individual native GTEx gene TPM, grouped within a sample using reviewed locus identities."""

import gzip
import json

import numpy as np
import pandas as pd
from antigen_sample_evidence import LACTATION_MARKERS, mammary_states_from_panel

from oncoref import curation, library_qc, reference_data, target_selection
from oncoref.gene_matching import gene_identity_catalog

PANEL = curation.manifest()["reference_panels"]["gtex_individual_samples"]
SOURCE = PANEL["source"]
RELEASE = PANEL["release"]
SOURCE_ID = f"{SOURCE}:{RELEASE}"


def sample_matrix(path, canonical):
    """Read exact versioned rows, leaving absent or unmeasured loci missing."""
    values = {}
    with gzip.open(path, "rb") as handle:
        version = handle.readline().strip()
        nrows, ncols = map(int, handle.readline().split())
        columns = handle.readline().decode().rstrip().split("\t")
        samples = pd.Index(columns[2:], name="sample_id")
        donors = samples.str.extract(r"^(GTEX-[^-]+)-.+$", expand=False)
        if (
            version != b"#1.2"
            or columns[:2] != ["Name", "Description"]
            or len(samples) != ncols
            or samples.has_duplicates
            or donors.hasnans
        ):
            raise ValueError(f"Unexpected GTEx sample header in {path}")
        count = 0
        for line in handle:
            count += 1
            name, _, rest = line.partition(b"\t")
            if name not in canonical:
                continue
            gene = canonical[name]
            if gene in values:
                raise ValueError(f"Duplicate GTEx locus: {gene}")
            vector = np.asarray(rest.partition(b"\t")[2].split(), dtype=float)
            if len(vector) != ncols or np.isinf(vector).any() or (vector < 0).any():
                raise ValueError(f"Invalid TPM values for {gene}")
            values[gene] = vector
        if count != nrows:
            raise ValueError(f"Unexpected GTEx row count in {path}")
    return pd.DataFrame.from_dict(values, orient="index", columns=samples), donors


def donor_p99(values, donors):
    """Linear P99 across measured donors after combining member loci within each sample."""
    if len(donors) != len(values.columns) or pd.isna(donors).any():
        raise ValueError("GTEx P99 requires a donor for every sample")
    # Give each donor one value if a later release includes replicate libraries.
    donor_values = values.T.groupby(donors, sort=False).median()
    return donor_values.quantile(0.99, interpolation="linear").dropna()


def individual_samples(facts, mapped, tissue_columns):
    labels = gene_identity_catalog().set_index("canonical_gene_id").canonical_symbol.fillna("")
    group = {g: r.entity_id for r in facts.itertuples() for g in json.loads(r.member_gene_ids)}
    counts = facts.set_index("entity_id").member_count
    markers = labels[
        labels.str.match(library_qc.OLFACTORY_PATTERN) | labels.isin(LACTATION_MARKERS)
    ]
    wanted = set(group) | set(markers.index) | set(library_qc.marker_panels().gene_id)
    loci = mapped[mapped.gene_id.isin(wanted)]
    if loci.Name.duplicated().any() or loci.gene_id.duplicated().any():
        raise ValueError("Ambiguous GTEx source-to-canonical mapping")
    canonical = {r.Name.encode(): r.gene_id for r in loci.itertuples()}
    rows, audits, seen = [], [], set()
    prefix = RELEASE + "-"
    pins = reference_data.REFERENCE_SOURCES[SOURCE]["pins"]
    if {prefix + name.lower() for name in tissue_columns} != {
        v for v in pins if v.startswith(prefix)
    }:
        raise ValueError("GTEx primary-tissue registry and anatomy crosswalk disagree")
    for native, tissue in sorted(tissue_columns.items()):
        version = prefix + native.lower()
        path = reference_data.pinned_path(SOURCE, version)
        genes, donors = sample_matrix(path, canonical)
        samples = genes.columns
        if seen.intersection(samples):
            raise ValueError(f"Repeated GTEx sample across tissue shards: {native}")
        seen.update(samples)
        audit = library_qc.signatures(
            genes,
            list(samples),
            gene_ids=pd.Series(genes.index),
            symbols=labels.reindex(genes.index),
        )
        audit["source_id"] = SOURCE_ID
        audit["tissue"] = tissue
        audit["source_tissue"] = native
        audit["participant_id"] = donors
        audit["library_qc_review"] = audit.apply(lambda r: ";".join(library_qc.reasons(r)), axis=1)
        audit["library_qc_notes"] = audit.apply(lambda r: ";".join(library_qc.notes(r)), axis=1)
        compartments = pd.Series(tissue, index=samples)
        if tissue == "breast":
            panel = genes.reindex(labels[labels.isin(LACTATION_MARKERS)].index)
            if len(panel) != len(LACTATION_MARKERS):
                raise ValueError("Incomplete canonical mammary marker panel")
            states = mammary_states_from_panel(panel)
            audit["mammary_state_inference"] = [
                json.dumps(states[s], allow_nan=False) for s in samples
            ]
            compartments = pd.Series(
                {
                    s: {
                        "secretory_signature": "lactating_breast",
                        "no_high_secretory_signature": "resting_breast",
                        "unknown": "breast",
                    }[states[s]["state"]]
                    for s in samples
                }
            )
        audits.append(audit.reset_index())
        members = genes.loc[genes.index.isin(group)]
        grouped = members.groupby(members.index.map(group))
        values = (
            grouped.sum(min_count=1).where(grouped.count().eq(counts, axis=0)).reindex(counts.index)
        )
        values.index.name = "entity_id"
        if tissue in {*target_selection.REPRODUCTIVE_COMPARTMENTS, "breast"}:
            frame = values.melt(
                ignore_index=False, var_name="sample_id", value_name="value"
            ).reset_index()
            frame["participant_id"] = frame.sample_id.map(audit.participant_id)
            frame["library_qc_review"] = frame.sample_id.map(audit.library_qc_review)
            frame["compartment"] = frame.sample_id.map(compartments)
            statistic = "individual sample gene TPM"
            kind = "individual_gtex_sample"
            source = SOURCE_ID
        else:
            frame = donor_p99(values, donors).rename("value").reset_index()
            frame = frame.assign(
                sample_id="", participant_id="", library_qc_review="", compartment=tissue
            )
            statistic = "within-site donor P99, linear interpolation"
            kind = "donor_percentile"
            source = PANEL["normal_source_id"]
        frame["member_count"] = frame.entity_id.map(counts)
        frame["measured_members"] = frame.member_count.where(frame.value.notna())
        rows.append(
            frame.assign(
                tissue=tissue,
                source_context=native,
                source_id=source,
                cell_type="bulk",
                statistic=statistic,
                observation_kind=kind,
                unit="TPM",
            )
        )
    return pd.concat(rows, ignore_index=True), pd.concat(audits, ignore_index=True)
