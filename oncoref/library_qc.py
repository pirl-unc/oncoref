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
"""Conservative tumor-library screens for background and tissue programs.

Three signatures that can flag a technical problem or admixed normal tissue:

* **Broad background.** Olfactory receptors are hundreds of genes that tumors do not
  express together. Genomic DNA or unspliced background gives many of them a few TPM
  at once; real expression involves a few receptors, usually the same ones across a
  tumor type.
* **Normal testis tissue.** Seminiferous tubules bring the whole late-spermatid
  program (protamines, transition proteins) at the same fraction of its testis level.
  A coordinated program is treated conservatively as possible admixture. RNA
  alone cannot establish its cell of origin or exclude tumor differentiation.
* **Several unrelated tissues.** One other tissue is common (adjacent normal tissue,
  the metastasis site); markers of four or more unrelated tissues in one library are
  not.

A tissue-program score is the *lowest* ratio, across its panel, of library TPM to
the tissue's GTEx median. It is not a calibrated cell fraction or purity estimate.
Panels come from :file:`tissue-marker-panels.csv` (``scripts/build_tissue_marker_panels.py``).
The thresholds remove whole libraries; they apply only to linear TPM sources.
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
#: Fail at this testis-tissue content: testis transcripts at ~10,000 TPM then reach
#: 10 TPM, the cancer rule's floor, on their own.
TESTIS_TISSUE_MAX_CONTENT = 0.001
#: A non-testis tissue counts as present at this content...
UNRELATED_TISSUE_CONTENT = 0.003
#: ...when the cohort's median content is below this (the tissue is not the cohort's own).
NATIVE_TISSUE_MEDIAN = 0.0003
#: Fail with this many unrelated tissues or more.
UNRELATED_TISSUE_MAX = 4

REASONS = {
    "broad_background": "More than 20% of measured olfactory receptors at 5 TPM or more",
    "normal_testis_tissue": "All six late-spermatid transcripts at 0.1% of testis or more",
    "several_unrelated_tissues": "Four or more unrelated tissues at 0.3% of their level or more",
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
    if pd.notna(signature.testis_tissue_content) and (
        signature.testis_tissue_content >= TESTIS_TISSUE_MAX_CONTENT
    ):
        out.append("normal_testis_tissue")
    value = signature.unrelated_tissues
    tissues = [] if pd.isna(value) else [t for t in str(value).split(";") if t]
    if len(tissues) >= UNRELATED_TISSUE_MAX:
        out.append("several_unrelated_tissues")
    return out
