# Licensed under the Apache License, Version 2.0.
"""Additional, source-aware normal-cell evidence for the CTA research screen.

Normal protein production excludes; strong RNA or antibody conflicts prompt
review (Neither), not a claim of clinical toxicity. Native RNA units stay separate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import normal_tissues

POLICY_ID = "oncoref.target-sets.2026-10-23.2"
BASE_POLICY_ID = "oncoref.target-sets.2026-10-23.1"
RELIABLE = {"Approved", "Supported", "Enhanced"}
POSITIVE = {"Low", "Medium", "High", "Ascending", "Descending"}
# Evidence-review triggers, not transplanted bulk ceilings or validated safety limits.
CLUSTER_NTPM = 5.0
CLUSTER_READS = 20
IMMUNE_TPM = 5.0
IMMUNE_SPECIMENS = 3
CORNEA_CP10K = 1.0
CORNEA_CELLS = 5
CORNEA_COUNTS = 20
CORNEA_GROUPS = 2


def contexts(frame):
    """Resolve each source/tissue/cell context explicitly, retaining missing maps."""
    keys = ["source", "version", "tissue", "cell_type"]
    if not frame.taxonomy_id.eq(9606).all():
        raise ValueError("This evidence screen is human-specific; do not transfer assay thresholds")
    rows = []
    for source, version, tissue, cell in (
        frame[keys].drop_duplicates().itertuples(index=False, name=None)
    ):
        a = normal_tissues.resolve_normal_tissue(
            tissue,
            cell_type=cell,
            species="human",
            source_id=source,
            source_version=version,
            stage="adult",
        )
        rows.append(
            dict(zip(keys, [source, version, tissue, cell], strict=True))
            | {
                "mapping_status": a.mapping_status,
                "concept_id": a.concept_id or "",
                "category": normal_tissues.normal_tissue_category(a),
            }
        )
    return frame.merge(pd.DataFrame(rows), on=keys, how="left", validate="many_to_one")


def cluster_findings(frame):
    """Review substantial normal-cell RNA at least as high as reproductive-cell RNA.

    Comparisons use ONLY this single-cell source. The read-count guard limits
    tiny-cluster/isolated-read artifacts; this remains a review flag, not proof
    of somatic protein. Missing reproductive cell coverage is explicitly retained.
    """
    d = contexts(frame)
    if not d.source.eq("hpa_single_cell_tissue").all() or not d.unit.eq("nTPM").all():
        raise ValueError("Cluster comparison requires the pinned HPA cluster nTPM source")
    r = d[d.category.eq("reproductive")].groupby("gene_id").value.max()
    d["reproductive_cell_max"] = d.gene_id.map(r)
    d["review"] = (
        d.value.ge(CLUSTER_NTPM)
        & d.read_count.ge(CLUSTER_READS)
        & (
            ~d.mapping_status.eq("mapped")
            | (
                d.category.isin(["vital", "protected", "manageable"])
                & (d.value.ge(d.reproductive_cell_max) | d.reproductive_cell_max.isna())
            )
        )
    )
    return d


def dice_findings(frame):
    """Review recurrent high TPM in an original DICE specimen distribution."""
    if not frame.taxonomy_id.eq(9606).all():
        raise ValueError("Human immune data required")
    if not frame.source.str.startswith("dice_").all() or not frame.unit.eq("TPM").all():
        raise ValueError("Use original DICE TPM distributions, not the HPA DICE summary")
    d = frame.assign(high=frame.value.ge(IMMUNE_TPM))
    result = d.groupby(
        ["gene_id", "gene_name", "source", "cell_type", "state"], as_index=False
    ).agg(
        measurements=("value", "size"),
        minimum=("value", "min"),
        median=("value", "median"),
        maximum=("value", "max"),
        high_measurements=("high", "sum"),
    )
    result["review"] = result.high_measurements.ge(IMMUNE_SPECIMENS)
    result["unit"] = "TPM"
    return result


def immune_findings(frame):
    """Flag substantial sorted-cell RNA; keep each study and its units separate.

    HPA DICE summaries are descriptive only: original DICE distributions supply
    that study's decision, so importing both never counts as replication.
    """
    if not frame.taxonomy_id.eq(9606).all():
        raise ValueError("Human immune data required")
    if not frame.source.isin(
        ["hpa_immune_cell", "hpa_immune_monaco", "hpa_immune_schmiedel"]
    ).all():
        raise ValueError("Expected HPA immune summary tables")
    if not frame.unit.isin(["nTPM", "TPM"]).all():
        raise ValueError("Native immune RNA units required")
    d = frame.copy()
    d["review"] = d.value.ge(5) & ~d.source.eq("hpa_immune_schmiedel")
    return d


def cornea_findings(frame):
    """Review supported adult ocular RNA in multiple annotated sample groups."""
    if not frame.taxonomy_id.eq(9606).all():
        raise ValueError("Human corneal data required")
    if not frame.unit.eq("CP10K").all() or not frame.stage.eq("adult").all():
        raise ValueError("Adult corneal CP10K evidence required")
    d = frame.assign(
        supported=frame.value.ge(CORNEA_CP10K)
        & frame.positive_cells.ge(CORNEA_CELLS)
        & frame.read_count.ge(CORNEA_COUNTS)
    )
    result = d.groupby(["gene_name", "cell_type"], as_index=False).agg(
        maximum=("value", "max"),
        positive_cells=("positive_cells", "sum"),
        supported_groups=("supported", "sum"),
    )
    result["review"] = result.supported_groups.ge(CORNEA_GROUPS)
    result["unit"] = "CP10K"
    return result


def antibody_findings(primary, *, ignored_genes=()):
    """Primary staining conflicts remain distinct from curated protein expression.

    An uncertain antibody within a reliable curated profile does not invalidate
    that whole profile. Existing explicitly reviewed unreliable gene assays are
    retained as observations but cannot supply exclusion or admission evidence.
    """
    d = primary.rename(
        columns={
            "Gene": "gene_id",
            "Gene name": "gene_name",
            "Tissue": "tissue",
            "Cell type": "cell_type",
            "Level": "level",
        }
    ).copy()
    d["source"], d["version"], d["taxonomy_id"] = "hpa_antibody_xml", "v23", 9606
    d = contexts(d)
    positive = d.level.isin(POSITIVE)
    normal = d.category.isin(["vital", "protected", "manageable"])
    ignored = d.gene_id.isin(ignored_genes)
    d["normal_signal"] = positive & normal & ~ignored
    # Conflicting/uncertain profiles need explicit review; never call them proven toxicity.
    d["uncertain_normal"] = d.normal_signal & ~d["Profile reliability"].isin(RELIABLE)
    d["reproductive_protein"] = (
        positive
        & d.category.eq("reproductive")
        & ~ignored
        & d["Profile reliability"].isin(RELIABLE)
        & d["Antibody validation"].isin(RELIABLE)
    )
    return d


def classify(base, summary, publications):
    """Apply evidence facts to a frozen baseline; no gene-specific decision overrides.

    Publications supply affirmative findings, not an allowlist. Reproductive +
    tumor evidence bypasses the admission floor only. Normal production and
    unresolved conflicts still win. Variable and prior absolute/protein/safety
    exclusions cannot be rescued by tumor evidence.
    """
    if base["Gene ID"].duplicated().any() or summary["Gene ID"].duplicated().any():
        raise ValueError("Unique gene identities required")
    if not base.Policy.eq(BASE_POLICY_ID).all():
        raise ValueError("Supply the frozen v1.16 target-set baseline")
    allowed = {"reproductive", "tumor", "recognition", "normal_protein"}
    if not set(publications.Finding) <= allowed or publications.Source.eq("").any():
        raise ValueError("Cited affirmative publication findings required")
    d = base.merge(summary, on="Gene ID", how="left", validate="one_to_one")
    for col in [
        "Cell RNA review",
        "Immune RNA review",
        "Corneal RNA review",
        "Antibody review",
        "Additional reproductive protein",
        "Cell evidence assessed",
    ]:
        d[col] = d[col].fillna(False).astype(bool)
    for col in [
        "Cell RNA findings",
        "Immune RNA findings",
        "Corneal RNA findings",
        "Antibody findings",
    ]:
        d[col] = d[col].fillna("")
    findings = publications.groupby("Gene ID").Finding.agg(set)
    support = {
        gid
        for gid, f in findings.items()
        if "reproductive" in f and ("tumor" in f or "recognition" in f)
    }
    normal = {gid for gid, f in findings.items() if "normal_protein" in f}
    citation = publications.groupby("Gene ID").Source.agg(
        lambda values: "; ".join(dict.fromkeys(values))
    )
    d["Previous category"] = d.Category
    d["Previous detail"] = d.Detail
    d["Baseline policy"] = d.Policy
    d["Publication support"] = d["Gene ID"].isin(support)
    d["Published normal protein"] = d["Gene ID"].isin(normal)
    d["Publication sources"] = d["Gene ID"].map(citation).fillna("")
    eligible = (
        d["Previous CTA screen"].eq("Retained")
        & ~d.Category.eq("Unsafe")
        & ~d.Detail.isin(["Variable", "Unassessed"])
    )
    r = pd.to_numeric(d["Reproductive RNA (nTPM)"], errors="coerce")
    weak = (
        d.Category.eq("CTA")
        & r.lt(5)
        & ~d["Reproductive positive control"]
        & ~d["Additional reproductive protein"]
        & ~d["Publication support"]
    )
    d.loc[weak, ["Category", "Detail", "Reason"]] = [
        "Neither",
        "Admission",
        "Reproductive bulk RNA below 5 nTPM without reliable reproductive protein or cited reproductive/tumor support",
    ]
    admission = (
        eligible & r.gt(0) & (d["Publication support"] | d["Additional reproductive protein"])
    )
    d.loc[admission, ["Category", "Detail", "Reason"]] = [
        "CTA",
        "Supported",
        "Affirmative reproductive and tumor/protein evidence supports admission below the bulk-RNA floor; normal gates still apply",
    ]
    review = d[
        ["Cell RNA review", "Immune RNA review", "Corneal RNA review", "Antibody review"]
    ].any(axis=1)
    held = review & ~d.Category.eq("Unsafe") & ~d.Detail.eq("Variable")
    d.loc[held, ["Category", "Detail", "Reason"]] = [
        "Neither",
        "Cellular",
        "Normal-cell RNA or uncertain antibody evidence requires review; not labelled proven clinical toxicity or added to forbidden proteins",
    ]
    d.loc[held & d["Antibody review"], "Detail"] = "Antibody"
    d.loc[d["Published normal protein"], ["Category", "Detail", "Reason"]] = [
        "Unsafe",
        "Normal",
        "Published protein production by nonmalignant cells outside permitted reproductive/disposable compartments",
    ]
    d["Target"], d["Forbidden"] = d.Category.eq("CTA"), d.Category.eq("Unsafe")
    d["Policy"] = POLICY_ID
    d["Change"] = np.where(
        d.Category.eq(d["Previous category"]),
        "Unchanged",
        d["Previous category"] + " → " + d.Category,
    )
    return d.sort_values(["Category", "Gene ID"]).reset_index(drop=True)
