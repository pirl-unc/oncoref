"""Human olfactory measurements in provider units and source feature names."""

import numpy as np
import pandas as pd

BULK = "olfactory_olender_transcriptome"
RECEPTORS = "olfactory_olender_receptors"
NUCLEI = "olfactory_mssm_nuclei"
CQ = "olfactory_verbeurgt_cq"
COPIES = "olfactory_verbeurgt_copies"
DURANTE = "olfactory_durante_supplement"
SOURCES = {BULK, RECEPTORS, NUCLEI, CQ, COPIES, DURANTE}
BIOPSIES = ("OE7", "OE12", "OE15", "OE17")


def _select(frame, genes):
    return frame if genes is None else frame[frame.gene_name.isin(genes)]


def _biopsies(path, source, genes):
    receptor = source == RECEPTORS
    feature = "geneName" if receptor else "symbol"
    samples = [s + "_Unique" if receptor else s for s in BIOPSIES]
    fields = [feature, "locus", *samples]
    data = pd.read_excel(path, sheet_name="Table S7" if receptor else "Table S2", usecols=fields)
    data = data.dropna(subset=[feature, "locus"]).rename(columns={feature: "gene_name"})
    data = _select(data, genes).melt(
        id_vars=["gene_name", "locus"], value_vars=samples, var_name="sample", value_name="value"
    )
    data["sample"] = data["sample"].str.removesuffix("_Unique")
    data["unit"] = "FPKM"
    data["provider_zero_floor"] = data.value.eq(0.003)
    return data


def _qpcr(path, source, genes):
    sheet = pd.read_excel(path, header=None)
    samples = list(range(1, 27))
    if list(sheet.iloc[1, 2:28]) != samples:
        raise ValueError("Unexpected olfactory qPCR sample columns")
    native = sheet.iloc[5:, 1:28].copy()
    native.columns = ["detector", *samples]
    labels = native.detector.str.extract(r"^(OR[^-]+)-(Hs.+)$")
    native[["gene_name", "assay_id"]] = labels
    data = _select(native.dropna(subset=["assay_id"]), genes).melt(
        id_vars=["gene_name", "assay_id"], value_vars=samples, var_name="donor", value_name="value"
    )
    data["value"] = pd.to_numeric(data.value, errors="raise")
    data["unit"] = "Cq" if source == CQ else "estimated RNA copies per 20 ng"
    return data


def _neurons(path, genes):
    data = pd.read_excel(path, sheet_name="SI Table 7", header=1, usecols=[0, 1])
    if data.columns[0] != "OR gene":
        raise ValueError("Unexpected Durante receptor table")
    data.columns = ["gene_name", "value"]
    return _select(data.dropna(subset=["gene_name", "value"]), genes).assign(
        unit="reported positive neurons", cell_type="olfactory sensory neuron"
    )


def _nuclei(path, genes):
    try:
        import h5py
    except ImportError as error:
        raise ImportError("Reading the olfactory H5AD requires oncoref[olfactory]") from error

    with h5py.File(path) as f:
        matrix = f["X"]
        labels = f["var/gene_name"].asstr()[:]
        if matrix.attrs["encoding-type"] != "csr_matrix" or not pd.Index(labels).is_unique:
            raise ValueError("Expected a CSR count matrix with unique feature names")
        selected = (
            np.arange(len(labels))
            if genes is None
            else np.flatnonzero(pd.Index(labels).isin(genes))
        )
        cell, donor = f["obs/Cell identity"], f["obs/donor_name"]
        osns = cell["codes"][:] == list(cell["categories"].asstr()[:]).index("OSNs")
        donors = donor["codes"][:]
        if (donors[osns] < 0).any():
            raise ValueError("Unassigned olfactory neuron donor")
        pointers = matrix["indptr"][:]
        rows = []
        for code, name in enumerate(donor["categories"].asstr()[:]):
            cells = np.flatnonzero(osns & (donors == code))
            if not len(cells):
                continue
            total = np.zeros(len(labels), dtype=np.int64)
            positive = np.zeros_like(total)
            for i in cells:
                span = slice(pointers[i], pointers[i + 1])
                indices = matrix["indices"][span]
                values = matrix["data"][span]
                if (
                    not np.isfinite(values).all()
                    or (values < 0).any()
                    or (values != np.floor(values)).any()
                ):
                    raise ValueError("Expected finite nonnegative integer UMI counts")
                np.add.at(total, indices, values.astype(np.int64))
                positive[np.unique(indices[values > 0])] += 1
            rows.append(
                pd.DataFrame(
                    {
                        "gene_name": labels[selected],
                        "donor": name,
                        "nuclei": len(cells),
                        "positive_nuclei": positive[selected],
                        "value": total[selected],
                    }
                )
            )
        return pd.concat(rows, ignore_index=True).assign(
            unit="UMI sum", cell_type="olfactory sensory neuron"
        )


def measurements(path, source, genes):
    if source in (BULK, RECEPTORS):
        return _biopsies(path, source, genes)
    if source in (CQ, COPIES):
        return _qpcr(path, source, genes)
    if source == DURANTE:
        return _neurons(path, genes)
    if source == NUCLEI:
        return _nuclei(path, genes)
    raise ValueError(f"Unknown olfactory source: {source}")
