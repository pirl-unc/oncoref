"""Published antigen evidence. Source lists and controls are independent of catalog membership."""

import pandas as pd

from . import catalog, cta_literature, cta_sources, reference_data, tissue_tiers
from .gene_matching import gene_identity_catalog

COLUMNS = ("list_id", "gene_id", "gene_name", "reported_name", "evidence", "source_url", "scope")


def tissues(*, tier: str | None = None) -> pd.DataFrame:
    """Human anatomical concepts with parent IDs, selection scope and adult injury tier."""
    data = tissue_tiers.anatomy_inventory()
    if tier is None:
        return data
    if tier not in {*tissue_tiers.TIERS, "unresolved"}:
        raise ValueError(f"Unknown injury tier {tier!r}")
    return data.loc[data.tier.eq(tier)].reset_index(drop=True)


def tissue(name: str, cell_type: str | None = None, *, stage: str = "adult") -> dict:
    """Loss consequences and linked evidence for an organ or a named cell population."""
    return tissue_tiers.consequences(
        name, cell_type, resolution="cell" if cell_type else "bulk", stage=stage
    )


def tiers() -> pd.DataFrame:
    """Ordered injury tier IDs and their definitions."""
    definitions = tissue_tiers.registry()["tier_definitions"]
    return pd.DataFrame(
        [(key, definitions[key]) for key in tissue_tiers.TIERS], columns=["tier", "definition"]
    )


def gene_lists() -> pd.DataFrame:
    """List identifiers for curated controls, source tables, and reviewed paper claims."""
    controls = catalog.read_table("cta-literature-benchmarks")
    tables = cta_sources.publication_sources()
    papers = cta_literature.publication_registry()["sources"]
    rows = [("benchmark/" + key, "classification control", "") for key in controls.list_id.unique()]
    rows.extend(
        ("table/" + r.source_tag, r.evidence_scope, r.source_url) for r in tables.itertuples()
    )
    rows.extend(("paper/" + r["source_id"], r["title"], r["url"]) for r in papers)
    rows.append(
        (
            "markers/human_primary_osteoblast_markers",
            "Primary human osteoblast subtype markers",
            reference_data.REFERENCE_SOURCES["human_primary_osteoblast_markers"]["publication_url"],
        )
    )
    return pd.DataFrame(rows, columns=["list_id", "description", "source_url"])


def gene_list(list_id: str) -> pd.DataFrame:
    """Read source-reported symbols, identity, finding and scope without dropping unmapped rows."""
    family, _, key = list_id.partition("/")
    if family == "benchmark":
        data = catalog.read_table("cta-literature-benchmarks")
        rows = data[data.list_id.eq(key)].rename(
            columns={"basis": "evidence", "qualification": "scope", "list_id": "control_basis"}
        )
        rows = rows.assign(reported_name=rows.gene_name)
    elif family == "table":
        data = cta_sources.publication_membership()
        source = cta_sources.publication_sources().set_index("source_tag")
        if key not in source.index:
            raise KeyError(list_id)
        rows = (
            data[data.source_tag.eq(key)]
            .rename(
                columns={
                    "Ensembl_Gene_ID": "gene_id",
                    "Symbol": "gene_name",
                    "source_symbol": "reported_name",
                    "evidence_type": "evidence",
                    "validation_scope": "scope",
                }
            )
            .assign(source_url=source.at[key, "source_url"])
        )
    elif list_id == "markers/human_primary_osteoblast_markers":
        spec = reference_data.REFERENCE_SOURCES[key]
        native = pd.read_excel(reference_data.pinned_path(key, spec["version"]), header=1)
        identity = gene_identity_catalog()[
            ["canonical_symbol", "canonical_gene_id"]
        ].drop_duplicates()
        lookup = (
            identity[~identity.canonical_symbol.duplicated(keep=False)]
            .set_index("canonical_symbol")
            .canonical_gene_id
        )
        rows = native.assign(
            gene_id=native.gene.map(lookup),
            gene_name=native.gene,
            reported_name=native.gene,
            evidence=native.cluster,
            source_url=spec["publication_url"],
            scope="Author-reported subtype enrichment in one adult donor",
        )
    elif family == "paper":
        papers = cta_literature.publication_registry()["sources"]
        source = next((r for r in papers if r["source_id"] == key), None)
        if source is None:
            raise KeyError(list_id)
        rows = (
            pd.DataFrame(source["gene_claims"])
            .explode("gene_ids")
            .rename(
                columns={
                    "gene_ids": "gene_id",
                    "protein": "reported_protein",
                    "reported_as": "evidence",
                }
            )
            .assign(source_url=source["url"])
        )
        symbols = gene_identity_catalog().set_index("canonical_gene_id").canonical_symbol
        rows["gene_name"] = rows.gene_id.map(symbols)
    else:
        raise KeyError(list_id)
    if rows.empty:
        raise KeyError(list_id)
    rows = rows.assign(list_id=list_id)
    return rows.reindex(columns=[*COLUMNS, *rows.columns.difference(COLUMNS)]).reset_index(
        drop=True
    )


def expression(panel: str, *, entity_ids=None, version=None) -> pd.DataFrame:
    """Catalog comparison observations with original TPM, source and sample identifiers.

    Normal, disposable and thymic return their respective scopes. Healthy
    returns the full tissue comparison, including reproductive tissues.
    Tumor rows are qualifying or above-normal witnesses. Use per-sample cohort
    matrices for all tumor expression values, including zeroes.
    """
    names = {
        "normal": "normal-comparison",
        "disposable": "normal-comparison",
        "healthy": "normal-comparison",
        "thymic": "normal-comparison",
        "reproductive": "reproductive-comparison",
        "tumor": "tumor-witnesses",
    }
    if panel not in names:
        raise ValueError(f"Unknown comparison panel {panel!r}. Choose from {tuple(names)}")
    dataset = "antigen-" + names[panel]
    wanted = (
        None
        if entity_ids is None
        else {entity_ids}
        if isinstance(entity_ids, str)
        else set(entity_ids)
    )
    parts = []
    for chunk in pd.read_csv(
        catalog.ensure(dataset, version=version), chunksize=100000, low_memory=False
    ):
        if panel in {"normal", "disposable", "thymic"}:
            chunk = chunk.loc[chunk.scope.eq(panel)]
        parts.append(chunk if wanted is None else chunk.loc[chunk.entity_id.isin(wanted)])
    return pd.concat(parts, ignore_index=True)


def sources() -> pd.DataFrame:
    """Input, comparison, assay and gene-list files with catalog selectors and SHA-256."""
    return catalog.read_table("antigen-source-inventory", keep_default_na=False)


def cell_atlases(*, version=None) -> pd.DataFrame:
    """Census datasets with provider-versioned native H5AD downloads and study identifiers."""
    data = catalog.read_table("antigen-dataset-versions", version=version)
    urls = data.citation.str.extract(
        r"(https://datasets\.cellxgene\.cziscience\.com/[^ ]+\.h5ad)", expand=False
    )
    if (
        urls.isna().any()
        or not urls.str.rsplit("/", n=1).str[-1].eq(data.dataset_version_id + ".h5ad").all()
    ):
        raise ValueError("Census citation and dataset version disagree")
    return data.assign(download_url=urls)


def download(dataset: str, *, version=None):
    """Download the exact file used in curation and verify its inventory digest."""
    from .curation import digest

    rows = sources()
    rows = rows[rows.dataset.eq(dataset)]
    if version is not None:
        rows = rows[rows.version.eq(version)]
    if rows.empty:
        raise KeyError((dataset, version))
    if len(rows) != 1:
        raise ValueError(f"Choose a version from antigen_data.sources() for {dataset}")
    record = rows.iloc[0]
    path = catalog.ensure(dataset, version=record.version or None)
    if digest(path) != record.sha256:
        raise ValueError(f"Content differs from the curated source {dataset} {record.version}")
    return path
