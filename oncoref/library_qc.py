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
"""Library exclusion screens and advisory tissue-program notes.

Three signatures that can flag a technical problem or admixed normal tissue:

* **Broad background.** Olfactory receptors are hundreds of genes that tumors do not
  express together. Genomic DNA or unspliced background gives many of them a few TPM
  at once; real expression involves a few receptors, usually the same ones across a
  tumor type.
* **Coordinated testis program (note only).** Reproductive-specific normal RNA
  markers can also be expressed by cancer. Their coexpression cannot establish
  normal-cell admixture and never excludes a library or changes its QC status.
* **Several unrelated tissues.** One other tissue is common (adjacent normal tissue,
  the metastasis site); markers of four or more unrelated tissues in one library are
  not.

A tissue-program score is the *lowest* ratio, across its panel, of library TPM to
the tissue's GTEx median. It is not a calibrated cell fraction or purity estimate.
Panels come from :file:`tissue-marker-panels.csv` (``scripts/build_tissue_marker_panels.py``).
Exclusion thresholds apply only to linear TPM sources. Reproductive content is
reported separately by :func:`notes`, never by :func:`reasons`.
"""

from __future__ import annotations

from functools import cache

import numpy as np
import pandas as pd

from .load_dataset import get_data

OLFACTORY_PATTERN = r"^OR\d+[A-Z]+\d+$"
#: A receptor at this TPM counts toward broad background.
BACKGROUND_TPM = 5.0
#: Fewer measured receptors than this cannot be judged.
BACKGROUND_MIN_RECEPTORS = 100
#: Fail when more than this share of measured receptors is at ``BACKGROUND_TPM`` (>9 of ~390 ORs).
BACKGROUND_MAX_SHARE = 0.025
#: Annotate a coordinated program; this ratio is not a normal-cell fraction.
TESTIS_PROGRAM_NOTE_THRESHOLD = 0.001
#: A non-testis tissue counts as present at this content...
UNRELATED_TISSUE_CONTENT = 0.003
#: ...when the cohort's median content is below this (the tissue is not the cohort's own).
NATIVE_TISSUE_MEDIAN = 0.0003
#: Fail with this many unrelated tissues or more.
UNRELATED_TISSUE_MAX = 4

REASONS = {
    "broad_background": "More than 2.5% of measured olfactory receptors at 5 TPM or more",
    "several_unrelated_tissues": "Four or more unrelated tissues at 0.3% of their level or more",
}
NOTES = {
    "R": (
        "Coordinated reproductive RNA: all six testis markers ≥0.1% of their reference; "
        "cell of origin unresolved. Advisory only."
    ),
}
COLUMNS = [
    "olfactory_receptors_measured",
    "olfactory_share_ge5",
    "testis_tissue_content",
    "unrelated_tissues",
]


@cache
def _panels() -> pd.DataFrame:
    return get_data("tissue-marker-panels.csv")


def marker_panels() -> pd.DataFrame:
    """Tissue marker panels: panel, gene_id, symbol, reference_tpm and source."""
    return _panels().copy()


def signatures(
    raw: pd.DataFrame,
    sample_cols: list[str],
    *,
    gene_ids: pd.Series | None = None,
    symbols: pd.Series | None = None,
) -> pd.DataFrame:
    """Per-sample background and tissue-content signatures for one cohort matrix.

    ``raw`` holds linear TPM with ``Ensembl_Gene_ID`` and ``Symbol`` columns unless
    ``gene_ids``/``symbols`` (aligned to ``raw``) are given. Versions are ignored.
    Cohort medians decide which tissues are the cohort's own, so pass a whole cohort.
    """
    gene_ids = (raw["Ensembl_Gene_ID"] if gene_ids is None else gene_ids).astype(str)
    gene_ids = gene_ids.str.split(".").str[0].to_numpy()
    symbols = (raw["Symbol"] if symbols is None else symbols).astype(str)
    values = raw[sample_cols].apply(pd.to_numeric, errors="coerce")
    values.index = gene_ids

    receptors = values[symbols.str.match(OLFACTORY_PATTERN).to_numpy()]
    measured = receptors.notna().sum()
    share = receptors.ge(BACKGROUND_TPM).sum() / measured.replace(0, np.nan)

    content = {}
    panels = _panels()
    for panel, rows in panels.groupby("panel", sort=True):
        reference = rows.set_index("gene_id").reference_tpm
        present = reference.index[reference.index.isin(values.index)]
        if len(present) < len(reference):
            content[panel] = pd.Series(np.nan, index=sample_cols)
            continue
        level = values.loc[present].groupby(level=0).sum(min_count=1)
        content[panel] = level.div(reference.loc[level.index], axis=0).min(skipna=False)
    content = pd.DataFrame(content, index=sample_cols)

    others = [p for p in content if p != "testis"]
    native = content[others].median().ge(NATIVE_TISSUE_MEDIAN)
    present = content[others].ge(UNRELATED_TISSUE_CONTENT) & ~native
    out = pd.DataFrame(index=pd.Index(sample_cols, name="sample_id"))
    out["olfactory_receptors_measured"] = measured.astype(int)
    out["olfactory_share_ge5"] = share
    out["testis_tissue_content"] = content.get("testis", pd.Series(np.nan, index=sample_cols))
    out["unrelated_tissues"] = present.apply(
        lambda row: ";".join(p for p in others if row[p]), axis=1
    )
    for panel in content:
        out[f"tissue_content:{panel}"] = content[panel]
    return out


def reasons(signature: pd.Series) -> list[str]:
    """QC failure reasons for one row of :func:`signatures`."""
    out = []
    if (
        signature.olfactory_receptors_measured >= BACKGROUND_MIN_RECEPTORS
        and signature.olfactory_share_ge5 > BACKGROUND_MAX_SHARE
    ):
        out.append("broad_background")
    value = signature.unrelated_tissues
    tissues = [] if pd.isna(value) else [t for t in str(value).split(";") if t]
    if len(tissues) >= UNRELATED_TISSUE_MAX:
        out.append("several_unrelated_tissues")
    return out


def notes(signature: pd.Series) -> list[str]:
    """Advisory observations that never affect QC status or sample eligibility.

    ``testis_tissue_content`` is a legacy column name for the minimum marker/reference
    expression ratio. It does not measure the fraction of normal testis cells.
    """
    if pd.notna(signature.testis_tissue_content) and (
        signature.testis_tissue_content >= TESTIS_PROGRAM_NOTE_THRESHOLD
    ):
        return ["R"]
    return []


#: Atomic library-rejection operations for bulk tumor cohorts (#780).
ALL_TUMOR_QC: tuple[str, ...] = (
    "broad_background",
    "unrelated_tissues",
    "normal_swaps",
)

TUMOR_QC_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_TUMOR_QC),
    "none": frozenset(),
    "raw": frozenset(),
}

#: Atomic continuous background-adjustment operations for bulk tumor cohorts (#780).
#: Reproductive / CTA / oncofetal programs are strictly protected and NEVER subtracted.
ALL_TUMOR_CLEANUP: tuple[str, ...] = (
    "index_hopping",
    "adjacent_normal_stroma",
)

TUMOR_CLEANUP_ALIASES: dict[str, frozenset[str]] = {
    "all": frozenset(ALL_TUMOR_CLEANUP),
    "none": frozenset(),
    "raw": frozenset(),
}

TUMOR_PROTECTED_REPRODUCTIVE_PANELS: frozenset[str] = frozenset(
    {
        "testis",
        "ovary",
        "placenta",
        "embryo",
        "epididymis",
        "fallopian_tube",
        "uterus",
        "cervix",
        "vagina",
        "testis_germline",
        "germ_cell",
        "trophoblast",
    }
)

NORMAL_SWAP_CONTENT_THRESHOLD = 0.25
NORMAL_SWAP_NATIVE_CEILING = 0.01

#: Explicit anatomical organ-of-origin panels for tumor codes so dedifferentiated cohorts
#: (e.g. PDAC/PAAD with low median islet INS) always recognize their own organ as native.
CANCER_CODE_NATIVE_PANELS: dict[str, frozenset[str]] = {
    "ACC": frozenset({"adrenal_gland"}),
    "PCPG": frozenset({"adrenal_gland", "brain"}),
    "PAAD": frozenset({"pancreas"}),
    "ACINIC": frozenset({"pancreas", "salivary_gland"}),
    "PANET": frozenset({"pancreas"}),
    "LIHC": frozenset({"liver"}),
    "CHOL": frozenset({"liver"}),
    "HB": frozenset({"liver"}),
    "FLC": frozenset({"liver"}),
    "LUAD": frozenset({"lung"}),
    "LUSC": frozenset({"lung", "squamous_epithelium"}),
    "SCLC": frozenset({"lung"}),
    "LCNEC": frozenset({"lung"}),
    "MESO": frozenset({"lung"}),
    "KIRC": frozenset({"kidney"}),
    "KIRP": frozenset({"kidney"}),
    "KICH": frozenset({"kidney"}),
    "WT": frozenset({"kidney"}),
    "THCA": frozenset({"thyroid"}),
    "MTC": frozenset({"thyroid"}),
    "GBM": frozenset({"brain"}),
    "LGG": frozenset({"brain"}),
    "MB": frozenset({"brain"}),
    "ATRT": frozenset({"brain"}),
    "EPN": frozenset({"brain"}),
    "DIPG": frozenset({"brain"}),
    "NB": frozenset({"brain", "adrenal_gland"}),
    "STAD": frozenset({"stomach"}),
    "ESCA": frozenset({"stomach", "squamous_epithelium"}),
    "ESCC": frozenset({"squamous_epithelium"}),
    "EAC": frozenset({"stomach"}),
    "HNSC": frozenset({"squamous_epithelium", "salivary_gland"}),
    "CESC": frozenset({"squamous_epithelium"}),
    "cSCC": frozenset({"squamous_epithelium"}),
    "BCC": frozenset({"squamous_epithelium"}),
    "ADCC": frozenset({"salivary_gland"}),
    "MEC": frozenset({"salivary_gland", "squamous_epithelium"}),
    "SARC": frozenset({"striated_muscle"}),
    "RMS": frozenset({"striated_muscle"}),
    "TGCT": frozenset({"testis"}),
}

_ONCOFETAL_AND_TROPHOBLAST_SYMBOL_PREFIXES: tuple[str, ...] = (
    "MAGE",
    "CTAG",
    "GAGE",
    "PAGE",
    "XAGE",
    "BAGE",
    "SSX",
    "SPANX",
    "PRAME",
    "CT45",
    "CT47",
    "CT83",
    "VCX",
    "DPPA",
    "PIWIL",
    "TDRD",
    "SPATA",
    "TEX",
    "HORMAD",
    "SYCP",
    "PSG",
    "CGB",
    "ERV",
    "IGF2BP",
    "LIN28",
)

_ONCOFETAL_AND_TROPHOBLAST_EXACT_SYMBOLS: frozenset[str] = frozenset(
    {
        "CGA",
        "LHB",
        "FSHB",
        "POMC",
        "INSM1",
        "NCAM1",
        "PLAC1",
        "TBXT",
        "AFP",
        "GPC3",
        "SALL4",
        "CLDN6",
        "SOX2",
        "POU5F1",
        "NANOG",
        "DLL3",
        "ASCL1",
        "SYP",
        "CHGA",
        "FOXE1",
        "DUX4",
        "ZAR1",
        "ZAR1L",
        "KHDC3L",
        "NLRP13",
        "BTN1A1",
        "CSN2",
        "LALBA",
    }
)


@cache
def _protected_reproductive_and_cta_gene_ids() -> frozenset[str]:
    """All Ensembl gene IDs belonging to CTA, reproductive, placental, or oncofetal programs.

    Enforces the hard biological invariant that ectopic tumor expression of CTAs
    (including published nominations that failed somatic safety screens such as ``MAGEA3``),
    placental/gonadotropin subunits (``CGA``, ``CGB*``, ``PSG*``, ``ERV*``), or oncofetal
    factors (``AFP``, ``GPC3``, ``IGF2BP*``, ``LIN28*``, ``SALL4``, ``CLDN6``) is NEVER
    subtracted as 'normal background' during tumor ``cleanup``.
    """
    import json

    from .gene_ids import canonical_gene_space

    protected: set[str] = set()
    panels_df = _panels()
    protected.update(
        panels_df.loc[
            panels_df["panel"].isin(TUMOR_PROTECTED_REPRODUCTIVE_PANELS), "gene_id"
        ].astype(str)
    )
    anchors_df = get_data("gtex-deconvolution-anchors")
    protected.update(
        anchors_df.loc[
            anchors_df["source_family"].isin(TUMOR_PROTECTED_REPRODUCTIVE_PANELS)
            | anchors_df["source_site"].isin(TUMOR_PROTECTED_REPRODUCTIVE_PANELS),
            "gene_id",
        ].astype(str)
    )
    ts_df = get_data("target-selection")
    repro_cta_ca_or_pub = (
        ts_df["is_cta"].astype(bool)
        | ts_df["is_reproductive_antigen"].astype(bool)
        | ts_df["is_cancer_antigen"].astype(bool)
        | ts_df["published_nomination"].fillna(False).astype(bool)
    )
    for members_json in ts_df.loc[repro_cta_ca_or_pub, "member_gene_ids"].dropna():
        for gid in json.loads(str(members_json)):
            protected.add(str(gid).split(".")[0])

    cspace = canonical_gene_space()
    sym_col = cspace["symbol"].fillna("").astype(str).str.upper()
    prefix_hit = sym_col.map(lambda s: s.startswith(_ONCOFETAL_AND_TROPHOBLAST_SYMBOL_PREFIXES))
    exact_hit = sym_col.isin(_ONCOFETAL_AND_TROPHOBLAST_EXACT_SYMBOLS)
    protected.update(
        cspace.loc[prefix_hit | exact_hit, "ensembl_gene_id"].astype(str).str.split(".").str[0]
    )
    return frozenset(protected)


def _somatic_nonreproductive_panels() -> tuple[dict[str, pd.Series], pd.DataFrame]:
    from . import bulk_qc

    protected_genes = _protected_reproductive_and_cta_gene_ids()
    panels_df = _panels()
    somatic_panels: dict[str, pd.Series] = {}
    prof_rows: dict[str, dict[str, float]] = {}
    for pname, grp in panels_df.groupby("panel", sort=True):
        pname_str = str(pname)
        if pname_str in TUMOR_PROTECTED_REPRODUCTIVE_PANELS:
            continue
        sub = grp[~grp["gene_id"].astype(str).isin(protected_genes)]
        if sub.empty:
            continue
        ser = sub.set_index("gene_id")["reference_tpm"].astype(float)
        ser = ser[~ser.index.duplicated(keep="first")]
        somatic_panels[pname_str] = ser
        for gid, val in ser.items():
            prof_rows.setdefault(str(gid), {})[pname_str] = float(val)

    anchors_df = get_data("gtex-deconvolution-anchors")
    for extra_site in ("adrenal_gland", "pituitary"):
        if extra_site not in somatic_panels:
            sub_a = anchors_df[
                anchors_df["source_site"].eq(extra_site)
                & ~anchors_df["gene_id"].astype(str).isin(protected_genes)
            ]
            if not sub_a.empty:
                ser_a = sub_a.set_index("gene_id")["source_median_tpm"].astype(float)
                ser_a = ser_a[~ser_a.index.duplicated(keep="first")]
                somatic_panels[extra_site] = ser_a
                for gid, val in ser_a.items():
                    prof_rows.setdefault(str(gid), {})[extra_site] = float(val)

    for row in anchors_df.itertuples(index=False):
        gid = str(row.gene_id)
        if gid in protected_genes:
            continue
        mapped = bulk_qc.GTEX_ANCHOR_SITE_TO_PANEL.get(str(row.source_site))
        if (
            mapped
            and mapped not in TUMOR_PROTECTED_REPRODUCTIVE_PANELS
            and mapped in somatic_panels
        ):
            prof_rows.setdefault(gid, {})[mapped] = max(
                prof_rows.get(gid, {}).get(mapped, 0.0), float(row.source_median_tpm)
            )
    assert not (set(somatic_panels) & TUMOR_PROTECTED_REPRODUCTIVE_PANELS)
    prof_df = pd.DataFrame.from_dict(prof_rows, orient="index").fillna(0.0)
    return somatic_panels, prof_df


def apply_tumor_bulk_spec(
    raw: pd.DataFrame,
    *,
    cancer_code: str | None = None,
    qc=None,
    cleanup=None,
    normalize: str | None = None,
) -> pd.DataFrame:
    """Apply ``(qc, cleanup, normalize)`` to a linear-TPM tumor matrix (#780).

    Hard biological invariants:
      - ``qc`` drops libraries with gDNA ``broad_background``, ``>= 4`` ``unrelated_tissues``,
        or gross non-reproductive ``normal_swaps``; coordinated reproductive/CTA expression
        (``testis_tissue_content``, Note ``R``) **never** rejects a tumor library.
      - ``cleanup`` (``index_hopping``, ``adjacent_normal_stroma``) subtracts strictly
        non-reproductive somatic spillover; reproductive / CTA / oncofetal genes are
        **never** subtracted.
    """
    from . import bulk_qc
    from .expression_engine import id_columns, sample_columns
    from .normalization import clean_tpm

    qc_set = bulk_qc.expand_ops(
        qc,
        allowed=ALL_TUMOR_QC,
        aliases=TUMOR_QC_ALIASES,
        param_name="qc",
        domain_label="tumor bulk",
    )
    cleanup_set = bulk_qc.expand_ops(
        cleanup,
        allowed=ALL_TUMOR_CLEANUP,
        aliases=TUMOR_CLEANUP_ALIASES,
        param_name="cleanup",
        domain_label="tumor bulk",
    )
    norm = bulk_qc.resolve_normalize_option(normalize, domain_label="tumor bulk")

    out = raw.copy()
    id_cols = id_columns(out)
    samp_cols = sample_columns(out)
    dropped_samples: list[str] = []
    code_key = str(cancer_code or "").strip().split("_")[0]
    native_panels_set: set[str] = set(
        CANCER_CODE_NATIVE_PANELS.get(str(cancer_code or ""), frozenset())
    ) | set(CANCER_CODE_NATIVE_PANELS.get(code_key, frozenset()))

    if samp_cols and (qc_set or cleanup_set):
        sig = signatures(out, samp_cols)
        somatic_content_cols = [
            c
            for c in sig.columns
            if c.startswith("tissue_content:")
            and c.removeprefix("tissue_content:") not in TUMOR_PROTECTED_REPRODUCTIVE_PANELS
        ]
        # Only infer native somatic panels empirically when cancer_code is unknown,
        # requiring >= 5% cohort median (never the 0.03% lane-bleed floor!).
        if not native_panels_set:
            for col_name in somatic_content_cols:
                pname = col_name.removeprefix("tissue_content:")
                if float(sig[col_name].median()) >= 0.05:
                    native_panels_set.add(pname)

        if qc_set:
            keep: list[str] = []
            for sid in samp_cols:
                srow = sig.loc[sid]
                drop = False
                if (
                    "broad_background" in qc_set
                    and srow.olfactory_receptors_measured >= BACKGROUND_MIN_RECEPTORS
                    and srow.olfactory_share_ge5 > BACKGROUND_MAX_SHARE
                ):
                    drop = True
                if not drop and "unrelated_tissues" in qc_set:
                    ut = srow.unrelated_tissues
                    t_list = [] if pd.isna(ut) else [t for t in str(ut).split(";") if t]
                    if len(t_list) >= UNRELATED_TISSUE_MAX:
                        drop = True
                if (
                    not drop
                    and "normal_swaps" in qc_set
                    and somatic_content_cols
                    and native_panels_set
                    and "R" not in notes(srow)
                ):
                    foreign_vals = [
                        float(srow[c])
                        for c in somatic_content_cols
                        if c.removeprefix("tissue_content:") not in native_panels_set
                        and pd.notna(srow[c])
                    ]
                    native_vals = [
                        float(srow[c])
                        for c in somatic_content_cols
                        if c.removeprefix("tissue_content:") in native_panels_set
                        and pd.notna(srow[c])
                    ]
                    max_foreign = max(foreign_vals, default=0.0)
                    max_native = max(native_vals, default=0.0)
                    if (
                        max_foreign >= NORMAL_SWAP_CONTENT_THRESHOLD
                        and max_native < NORMAL_SWAP_NATIVE_CEILING
                    ):
                        drop = True
                if drop:
                    dropped_samples.append(sid)
                else:
                    keep.append(sid)
            samp_cols = keep
            out = out.loc[:, [*id_cols, *samp_cols]].copy()

    if samp_cols and cleanup_set:
        somatic_panels, prof_df = _somatic_nonreproductive_panels()
        unver_gids = (
            out["Ensembl_Gene_ID"].astype(str).str.replace(r"\.\d+(?=_PAR_Y$|$)", "", regex=True)
        )
        protected_set = set(_protected_reproductive_and_cta_gene_ids())
        if "Symbol" in out.columns:
            sym_upper = out["Symbol"].fillna("").astype(str).str.strip().str.upper()
            extra_prot = sym_upper.map(
                lambda s: (
                    s.startswith(_ONCOFETAL_AND_TROPHOBLAST_SYMBOL_PREFIXES)
                    or s in _ONCOFETAL_AND_TROPHOBLAST_EXACT_SYMBOLS
                )
            )
            protected_set.update(unver_gids.loc[extra_prot].astype(str))
        protected_genes = frozenset(protected_set)
        vals = out[samp_cols].apply(pd.to_numeric, errors="coerce")
        orig_index = vals.index
        vals.index = pd.Index(unver_gids.to_numpy(), name="gene_id")
        sample_native_map = {sid: frozenset(native_panels_set) for sid in samp_cols}
        if "index_hopping" in cleanup_set:
            vals = bulk_qc.subtract_panel_index_hopping(
                vals,
                somatic_panels,
                sample_native_panels=sample_native_map,
                protected_panels=TUMOR_PROTECTED_REPRODUCTIVE_PANELS,
                protected_features=protected_genes,
                full_profiles=prof_df,
            )
        if "adjacent_normal_stroma" in cleanup_set:
            vals = bulk_qc.deconvolve_panel_admixture(
                vals,
                somatic_panels,
                sample_native_panels=sample_native_map,
                protected_panels=TUMOR_PROTECTED_REPRODUCTIVE_PANELS,
                protected_features=protected_genes,
                full_profiles=prof_df,
                max_subtract_fraction=0.50,
            )
        vals.index = orig_index
        out.loc[:, samp_cols] = vals

    if samp_cols and norm == "clean_tpm":
        cleaned = clean_tpm(out[samp_cols], gene_table=out[id_cols], on_empty_compartment="record")
        out = pd.concat(
            [out[id_cols].reset_index(drop=True), cleaned.reset_index(drop=True)], axis=1
        )
        if "unfilled_compartments" in cleaned.attrs:
            out.attrs["unfilled_compartments"] = dict(cleaned.attrs["unfilled_compartments"])

    out.attrs["bulk_spec"] = {
        "cancer_code": cancer_code,
        "qc": sorted(qc_set),
        "cleanup": sorted(cleanup_set),
        "normalize": norm,
        "dropped_samples": dropped_samples,
        "reproductive_cta_protected": True,
    }
    return out


def tumor_gdna_background_deltas(*, only_applied: bool = False) -> pd.DataFrame:
    """Return SHA-256-verified per-protein tumor witness adjustments under ``BACKGROUND_MAX_SHARE`` (2.5% ORs >= 5 TPM)."""
    import hashlib
    import json
    from importlib.resources import files

    data_root = files("oncoref").joinpath("data")
    csv_res = data_root.joinpath("tumor-gdna-background-deltas.csv")
    meta_res = data_root.joinpath("tumor-gdna-background-deltas.json")
    raw_bytes = csv_res.read_bytes()
    receipt = json.loads(meta_res.read_text(encoding="utf-8"))
    if hashlib.sha256(raw_bytes).hexdigest() != receipt["sha256"]:
        raise ValueError("tumor-gdna-background-deltas.csv does not match its SHA-256 receipt")
    df = get_data("tumor-gdna-background-deltas")
    if only_applied:
        df = df[df["applied"].astype(bool)].reset_index(drop=True)
    return df
