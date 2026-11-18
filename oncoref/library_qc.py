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
#: Fail when more than this share of measured receptors is at ``BACKGROUND_TPM``.
BACKGROUND_MAX_SHARE = 0.20
#: Annotate a coordinated program; this ratio is not a normal-cell fraction.
TESTIS_PROGRAM_NOTE_THRESHOLD = 0.001
#: A non-testis tissue counts as present at this content...
UNRELATED_TISSUE_CONTENT = 0.003
#: ...when the cohort's median content is below this (the tissue is not the cohort's own).
NATIVE_TISSUE_MEDIAN = 0.0003
#: Fail with this many unrelated tissues or more.
UNRELATED_TISSUE_MAX = 4

REASONS = {
    "broad_background": "More than 20% of measured olfactory receptors at 5 TPM or more",
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
