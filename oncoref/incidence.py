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

"""Cancer incidence / mortality burden reference data and the registry-driven
resolution of a cancer type to its anatomic burden category (and explicit
multi-bucket / subsite scope metadata for mixed cohorts such as ``CHOL``)."""

from __future__ import annotations

from functools import lru_cache

from .cancer_types import cancer_type_registry, resolve_cancer_type
from .load_dataset import _register_derived_cache, get_data

_BURDEN_METRICS = (
    "us_incidence_pct",
    "us_mortality_pct",
    "world_incidence_pct",
    "world_mortality_pct",
)


def cancer_burden_df():
    """Return the curated ``cancer-incidence-mortality.csv`` reference: each
    cancer **burden category**'s share (%) and extracted counts of annual cancer
    incidence and mortality for the US (ACS Cancer Facts & Figures 2024) and
    worldwide (GLOBOCAN 2024, Sung et al. 2026 ``PMID:42417444`` /
    ``DOI:10.3322/caac.70090``; minor sites not in the 34-site Table 1 retain
    their documented GLOBOCAN 2022 or literature vintage).

    The percentage columns are the primary lookup surface. Companion ``*_count``
    and ``*_total`` columns, ``*_source_locator`` fields, ICD-10 site fields
    (``source_site_codes``, ``included_source_sites``, ``excluded_source_sites``),
    ``derivation_basis``, and ``rounding_rule`` expose the structured audit
    contract. ``aggregation`` and ``source_anchor`` preserve the pirlygenes
    compatibility surface. Rows flagged as ``residual`` or overlapping minor
    sites (such as ``anus`` ``C21``, folded into ``colorectal`` ``C18-C21`` in
    GLOBOCAN 2024 Table 1) must not be summed blindly."""
    return get_data("cancer-incidence-mortality")


def cancer_code_burden_map():
    """Return ``{cancer_code: burden_category}`` from
    ``cancer-code-burden-map.csv``. This holds only the small set of
    **overrides** the registry ontology can't express on its own (e.g.
    ``SARC_KS`` -> Kaposi rather than soft-tissue; ``LAML`` -> AML rather than
    other-leukemia; ``HL`` -> Hodgkin; ``CTCL`` -> non-Hodgkin). Everything else
    is resolved by :func:`burden_category` from the registry's family +
    primary_tissue."""
    df = get_data("cancer-code-burden-map")
    return dict(zip(df["cancer_code"].astype(str), df["burden_category"].astype(str)))


def cancer_burden(category=None, *, metric="us_incidence_pct"):
    """Burden share (%) for one category and metric, or the whole
    ``{category: pct}`` map. ``metric`` is one of ``us_incidence_pct``,
    ``us_mortality_pct``, ``world_incidence_pct``, ``world_mortality_pct``."""
    if metric not in _BURDEN_METRICS:
        raise ValueError(f"metric must be one of {_BURDEN_METRICS}")
    df = cancer_burden_df()
    mapping = dict(zip(df["burden_category"].astype(str), df[metric].astype(float)))
    if category is None:
        return mapping
    return mapping.get(category)


_BONE_SARCOMA_TISSUES = {"bone", "cartilage", "notochord"}


@lru_cache(maxsize=1)
def _tissue_burden_map() -> dict[str, str]:
    """``{primary_tissue: burden_category}`` from ``tissue-burden-map.csv``."""
    df = get_data("tissue-burden-map")
    return dict(zip(df["primary_tissue"].astype(str), df["burden_category"].astype(str)))


@lru_cache(maxsize=1)
def _family_burden_map() -> dict[str, str]:
    """``{family: burden_category}`` fallback from ``family-burden-map.csv``."""
    df = get_data("family-burden-map")
    return dict(zip(df["family"].astype(str), df["burden_category"].astype(str)))


@lru_cache(maxsize=1)
def _burden_scope_overrides() -> dict[str, dict]:
    """Explicit scope records from ``cancer-burden-scope.csv`` keyed by cancer code."""
    df = get_data("cancer-burden-scope")
    out: dict[str, dict] = {}
    for row in df.itertuples(index=False):
        code = str(row.cancer_code).strip()
        sec_raw = (
            ""
            if row.secondary_burden_categories != row.secondary_burden_categories
            else str(row.secondary_burden_categories or "")
        )
        secondary = tuple(s.strip() for s in sec_raw.split(";") if s.strip())
        counts_raw = (
            ""
            if row.expression_source_site_counts != row.expression_source_site_counts
            else str(row.expression_source_site_counts or "")
        )
        site_counts: dict[str, int] = {}
        for item in counts_raw.split(";"):
            item = item.strip()
            if ":" in item:
                k, v = item.split(":", 1)
                site_counts[k.strip()] = int(v.strip())
        notes = "" if row.scope_notes != row.scope_notes else str(row.scope_notes or "")
        out[code] = {
            "scope_status": str(row.scope_status).strip(),
            "secondary_burden_categories": secondary,
            "expression_source_site_counts": site_counts,
            "notes": notes,
        }
    return out


_register_derived_cache(_tissue_burden_map.cache_clear)
_register_derived_cache(_family_burden_map.cache_clear)
_register_derived_cache(_burden_scope_overrides.cache_clear)


def burden_category(cancer_type):
    """Map a cancer type to its primary **incidence/mortality (cancer-burden)
    category** — the coarse anatomic bucket under which ACS / GLOBOCAN report
    population burden (e.g. ``SARC_OS`` → ``"bone_and_joint"``).

    For cohorts spanning multiple mortality buckets (such as ``CHOL``, which
    contains both intrahepatic ``liver`` ``C22.1`` and extrahepatic
    ``gallbladder_biliary`` ``C24`` cases), :func:`burden_category` returns the
    single primary bucket without duplicating population totals; use
    :func:`burden_categories`, :func:`burden_scope`, or
    :func:`cohorts_for_burden_category` when testing whether an expression
    cohort covers a mortality category."""
    try:
        code = resolve_cancer_type(cancer_type)
    except ValueError:
        return None
    if code is None:
        return None
    override = cancer_code_burden_map()
    tissue_map = _tissue_burden_map()
    family_map = _family_burden_map()
    reg = cancer_type_registry().set_index("code")
    # 1. explicit override (true exceptions only), walking up the parent chain
    cur, seen = code, set()
    while cur and cur not in seen:
        seen.add(cur)
        if cur in override:
            return override[cur]
        if cur not in reg.index:
            break
        cur = str(reg.loc[cur].get("parent_code", "") or "").strip() or None
    # 2. registry-driven, walking up the parent chain for blank tissue/family
    cur, seen = code, set()
    while cur and cur not in seen:
        seen.add(cur)
        if cur not in reg.index:
            break
        row = reg.loc[cur]
        family = str(row.get("family", "") or "")
        tissue = str(row.get("primary_tissue", "") or "")
        if family == "sarcoma":
            return "bone_and_joint" if tissue in _BONE_SARCOMA_TISSUES else "soft_tissue_sarcoma"
        if family == "heme-plasma":
            return "multiple_myeloma"
        if tissue in tissue_map:
            return tissue_map[tissue]
        if family in family_map:
            return family_map[family]
        cur = str(row.get("parent_code", "") or "").strip() or None
    return None


def burden_scope(cancer_type) -> dict | None:
    """Return structured anatomic-burden scope metadata for ``cancer_type``,
    resolving multi-bucket cohorts (``CHOL``, ``BTC``) and explicit subsite
    inclusions (salivary gland cohorts ``SGC`` / ``ADCC`` inside
    ``head_and_neck`` ``C07-C08``).

    Walks ``parent_code`` only while the child retains the same primary
    ``burden_category`` (so ``ADCC`` inherits ``SGC``'s salivary subsite note,
    whereas ``GBC`` keeps ``gallbladder`` ``exact`` scope instead of inheriting
    ``BTC``'s pan-biliary ``mixed_sites`` scope)."""
    try:
        code = resolve_cancer_type(cancer_type)
    except ValueError:
        return None
    if code is None:
        return None
    primary = burden_category(code)
    if primary is None:
        return None
    overrides = _burden_scope_overrides()
    reg = cancer_type_registry().set_index("code")
    matched = None
    cur, seen = code, set()
    while cur and cur not in seen:
        seen.add(cur)
        if cur in overrides:
            # Only inherit a parent's scope override if the child shares the
            # same primary burden category and has not narrowed to a single site
            # (e.g. GBC has primary_tissue='gallbladder' under parent BTC).
            if cur == code:
                matched = overrides[cur]
                break
            child_tissue = (
                str(reg.loc[code].get("primary_tissue", "") or "") if code in reg.index else ""
            )
            parent_tissue = (
                str(reg.loc[cur].get("primary_tissue", "") or "") if cur in reg.index else ""
            )
            if child_tissue == parent_tissue and burden_category(cur) == primary:
                matched = overrides[cur]
                break
        if cur not in reg.index:
            break
        cur = str(reg.loc[cur].get("parent_code", "") or "").strip() or None

    if matched is None:
        scope_status = "unassigned_umbrella" if primary == "other_and_unknown_primary" else "exact"
        return {
            "code": code,
            "primary_burden_category": primary,
            "burden_categories": (primary,),
            "scope_status": scope_status,
            "expression_source_site_counts": {},
            "notes": "",
        }
    categories = tuple(dict.fromkeys((primary, *matched["secondary_burden_categories"])))
    return {
        "code": code,
        "primary_burden_category": primary,
        "burden_categories": categories,
        "scope_status": matched["scope_status"],
        "expression_source_site_counts": dict(matched["expression_source_site_counts"])
        if cur == code
        else {},
        "notes": matched["notes"],
    }


def burden_categories(cancer_type) -> tuple[str, ...]:
    """Return all ``burden_category`` buckets spanned by ``cancer_type``
    (primary bucket first, followed by any secondary anatomic buckets for
    mixed-site cohorts such as ``CHOL`` → ``('gallbladder_biliary', 'liver')``).
    Returns ``()`` if ``cancer_type`` does not resolve."""
    scope = burden_scope(cancer_type)
    return () if scope is None else scope["burden_categories"]


def cohorts_for_burden_category(category: str, *, include_mixed: bool = True) -> tuple[str, ...]:
    """Return canonical registry cancer codes that map to ``category``.

    When ``include_mixed=True`` (default for expression-coverage queries),
    mixed-site codes whose secondary buckets include ``category`` (such as
    ``CHOL`` and ``BTC`` for ``"liver"``) are included alongside primary-mapped
    codes. Pass ``include_mixed=False`` when partitioning codes into disjoint
    buckets so population burden shares are never counted twice."""
    valid = set(cancer_burden_df()["burden_category"].astype(str))
    if category not in valid:
        raise ValueError(f"Unknown burden_category {category!r}")
    reg = cancer_type_registry()
    out = []
    for code in reg["code"].astype(str):
        cats = burden_categories(code) if include_mixed else (burden_category(code),)
        if category in cats:
            out.append(code)
    return tuple(out)
