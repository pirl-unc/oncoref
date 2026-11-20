"""Advisory cross-tissue signatures for individual bulk RNA libraries."""

import pandas as pd

# Nieuwenhuis et al. (2020), doi:10.1038/s41467-020-15821-9.
# Single-site TPM excess is an advisory measure of a coordinated foreign program.
CARRYOVER_PANELS = {
    "P": ("Pancreas", ("PRSS1", "CLPS", "PNLIP", "CELA3A")),
    "E": ("Esophagus_Mucosa", ("KRT4", "KRT13")),
}


def carryover_reference(mapped):
    """Native site medians from the pinned GTEx gene crosswalk, keyed by locus."""
    panels = []
    for code, (native, symbols) in CARRYOVER_PANELS.items():
        rows = mapped.loc[mapped.Description.isin(symbols), ["gene_id", "Description", native]]
        if sorted(rows.Description) != sorted(symbols) or not rows[native].gt(0).all():
            raise ValueError(f"Missing or ambiguous carryover marker for {native}")
        panels.append(
            rows.set_index("gene_id")[[native]].rename(columns={native: "scale"}).assign(code=code)
        )
    return pd.concat(panels, verify_integrity=True)


def published_carryover(values, site, reference):
    """Screen a single native-TPM site against its peers. NA means unassessed."""
    result = pd.DataFrame(index=values.columns)
    for code, (native, _) in CARRYOVER_PANELS.items():
        field = "carryover_" + code
        result[field] = pd.Series(pd.NA, index=result.index, dtype="boolean")
        if site.lower().startswith(native.split("_")[0].lower()):
            continue
        scale = reference.loc[reference.code.eq(code), "scale"]
        panel = values.reindex(scale.index).dropna(axis=1)
        if panel.shape[1] < 3:
            continue
        excess = panel.sub(panel.median(axis=1), axis=0).div(scale, axis=0)
        high = panel.gt(panel.quantile(0.95, axis=1), axis=0)
        result.loc[panel.columns, field] = (high & excess.ge(0.001)).all(axis=0)
    return result
