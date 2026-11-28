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

"""Functional biology & adult normal micro-niche review across all 66 baseline ``CA only`` proteins.

Documents what every ``Cancer Antigen only`` (``CA \\ RA``) protein in ``oncoref.target-selection.2026-11-26.2``
does biologically and adjudicates 13 adult somatic, neuroendocrine, sensory, osteoclast,
or near-identical ``Unsafe``-paralog proteins out of ``CA only``:

- **7 confirmed adult human ``Vital`` / ``Protected`` normal-tissue exclusions (``Unsafe``)**:
  ``GALP`` (adult human hypothalamic infundibular nucleus neuroendocrine neurons, ``HYPOMAP``),
  ``BHLHE23`` (adult human retinal rod and ON cone bipolar neurons),
  ``IL36B`` (adult human epidermal keratinocytes and lingual/esophageal squamous mucosa),
  ``ZACN`` (adult human pancreas, brain, heart, and olfactory/gustatory mucosa channel),
  ``HES3`` (adult human pancreatic islet alpha and beta cells),
  ``OCSTAMP`` (adult human RANKL-differentiated multinucleated bone-resorbing osteoclasts), and
  ``PROKR1`` (adult human adrenal medulla, endocrine microvascular, and GI smooth-muscle GPCR).
- **6 adult sensory / manageable-cell / near-identical ``Unsafe``-paralog holds (``Other:Review``)**:
  ``SRMS`` (HPA v25.1 gastric chief cells ``19.7 nCPM`` and goblet cells ``15.7 nCPM``),
  ``OCM2`` (``100.0%`` 1-mismatch 9-mer overlap / ``0`` safe ``>=2-mm`` 9-mers vs ``Unsafe`` ``OCM``,
  plus retinal ganglion cells ``7.7 nCPM`` and inner-ear outer hair cells),
  ``CIB3`` (unscreened ``Vital`` inner-ear vestibular hair-cell ``MET`` auxiliary subunit),
  ``OTOP1`` (diluted Type III sour-taste receptor & inner-ear vestibular otoconia proton channel),
  ``GFRA4`` (adult human thyroid calcitonin-producing parafollicular C-cells, ``Manageable`` cell
  type diluted ``~100x`` in whole-thyroid bulk RNA), and
  ``MRGPRG`` (Mas-related sensory GPCR with unscreened adult human dorsal root ganglia).
- **53 retained post-adjudication ``CA only`` proteins**, partitioned across five concrete
  somatic-silent biological classes (sub-``10x`` male germline/spermatogenic proteins,
  sub-``10x`` ``CT-X`` ``CSAG``/``MAGE``/``SSX`` paralogs, post-gastrulation embryonic homeobox
  and fetal ameloblast morphogens, ``Disposable`` hair-cortex ``KRTAP`` proteins, and
  hominoid pericentromeric/``RS447`` macrosatellite duplications and retrogenes).
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files

import pandas as pd

from .load_dataset import get_data

CA_ONLY_REVIEW_POLICY_ID = "oncoref.ca-only-functional-review.2026-11-27.1"

CA_ONLY_DISPOSITIONS: tuple[str, ...] = (
    "adjudicated_unsafe_adult_human",
    "adjudicated_review_niche_or_paralog",
    "retained_ca_only",
)

CA_ONLY_FUNCTIONAL_CLASSES: tuple[str, ...] = (
    "sub_10x_male_germline_and_spermatogenic",
    "sub_10x_ctx_cancer_germline_paralog",
    "post_gastrulation_fetal_morphogen_or_homeobox",
    "disposable_hair_cortex_keratin",
    "pericentromeric_macrosatellite_duplication_or_retrogene",
    "adult_neuroendocrine_sensory_or_somatic_niche",
)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _verify_receipt() -> dict:
    data_root = files("oncoref").joinpath("data")
    csv_bytes = data_root.joinpath("ca-only-functional-review.csv").read_bytes()
    receipt = json.loads(
        data_root.joinpath("ca-only-functional-review.json").read_text(encoding="utf-8")
    )
    if _sha256_bytes(csv_bytes) != receipt["sha256"]:
        raise ValueError("ca-only-functional-review.csv does not match its SHA-256 receipt")
    return receipt


def ca_only_review_metadata() -> dict:
    """Return SHA-256 provenance and summary counts for the 66-gene ``CA only`` review."""
    return _verify_receipt()


def ca_only_functional_review(
    protein: str | None = None,
    *,
    disposition: str | None = None,
    functional_class: str | None = None,
) -> pd.DataFrame:
    """Return the SHA-256-verified 66-protein ``CA only`` functional & micro-niche review table."""
    _verify_receipt()
    df = get_data("ca-only-functional-review").copy()
    if disposition is not None:
        if disposition not in CA_ONLY_DISPOSITIONS:
            raise ValueError(
                f"Unknown disposition {disposition!r}; expected one of {CA_ONLY_DISPOSITIONS}"
            )
        df = df[df["adjudication_disposition"].eq(disposition)].reset_index(drop=True)
    if functional_class is not None:
        if functional_class not in CA_ONLY_FUNCTIONAL_CLASSES:
            raise ValueError(
                f"Unknown functional_class {functional_class!r}; expected one of {CA_ONLY_FUNCTIONAL_CLASSES}"
            )
        df = df[df["functional_class"].eq(functional_class)].reset_index(drop=True)
    if protein is not None:
        q = str(protein).strip()
        if not q:
            raise ValueError("protein query must be nonempty")
        q_fold = q.casefold()
        mask = (
            df["entity_id"].astype(str).eq(q)
            | df["protein_group_label"].astype(str).str.casefold().eq(q_fold)
            | df["member_gene_ids"]
            .astype(str)
            .map(lambda s: q in json.loads(s) if s.startswith("[") else q in s.split("|"))
        )
        df = df[mask].reset_index(drop=True)
    return df


def ca_only_niche_adjudications() -> pd.DataFrame:
    """Return the 13 ``CA only`` proteins adjudicated to ``Unsafe`` (7) or ``Other:Review`` (6)."""
    df = ca_only_functional_review()
    return df[~df["adjudication_disposition"].eq("retained_ca_only")].reset_index(drop=True)


def ca_only_class_summary() -> pd.DataFrame:
    """Summarize all 66 baseline ``CA only`` proteins by disposition and biological mechanism."""
    df = ca_only_functional_review()
    grouped = (
        df.groupby(
            ["adjudication_disposition", "post_adjudication_classification", "functional_class"],
            as_index=False,
        )
        .agg(
            n_proteins=("entity_id", "count"),
            median_tumor_fold=("tumor_fold", "median"),
            max_patients_at_tenfold=("independent_patients_at_tenfold", "max"),
            examples=(
                "protein_group_label",
                lambda s: ", ".join(list(s.astype(str))[:6]),
            ),
        )
        .sort_values(["adjudication_disposition", "n_proteins"], ascending=[True, False])
        .reset_index(drop=True)
    )
    grouped["median_tumor_fold"] = grouped["median_tumor_fold"].round(1)
    return grouped


def annotate_catalog(frame: pd.DataFrame, *, enabled: bool = True) -> pd.DataFrame:
    """Attach ``ca_only_review_disposition`` and ``ca_only_functional_class`` columns."""
    out = frame.copy()
    csv_res = files("oncoref").joinpath("data/ca-only-functional-review.csv")
    if not enabled or not csv_res.is_file():
        out["ca_only_review_disposition"] = ""
        out["ca_only_functional_class"] = ""
        return out
    rev = ca_only_functional_review().set_index("entity_id")
    disp_map = rev["adjudication_disposition"].astype(str).to_dict()
    cls_map = rev["functional_class"].astype(str).to_dict()
    eids = out["entity_id"].astype(str)
    out["ca_only_review_disposition"] = [disp_map.get(e, "") for e in eids]
    out["ca_only_functional_class"] = [cls_map.get(e, "") for e in eids]
    return out


__all__ = [
    "CA_ONLY_DISPOSITIONS",
    "CA_ONLY_FUNCTIONAL_CLASSES",
    "CA_ONLY_REVIEW_POLICY_ID",
    "annotate_catalog",
    "ca_only_class_summary",
    "ca_only_functional_review",
    "ca_only_niche_adjudications",
    "ca_only_review_metadata",
]
