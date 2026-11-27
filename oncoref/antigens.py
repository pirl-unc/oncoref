# Licensed under the Apache License, Version 2.0.
"""Unified, plain-language antigen discovery and filtering API across human, dog, and mouse.

``oncoref.antigens`` is the single front-door entry point for querying biological
antigen sets (``CTA``, ``CA``, ``RA``, ``Unsafe``, ``Other``), physiological sex
restrictions (``CTA-M``, ``CTA-F``, ``CTA-P``, ``CTA-MF``), patient-specific organ
protection (``sex="female" | "male"``, ``pregnant``, ``lactating``, ``protect``,
``accept_loss``, ``absent``), ectopic tumor activation, and self-proteome peptide
subtraction—without needing to instantiate internal filter or host classes.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import replace

import pandas as pd

from . import antigen_sets, host_context, target_filters
from .species import available_species, resolve_species
from .species_target_selection import (
    anatomy_for_species,
    packaged_species_host_exposure,
    packaged_species_results,
)

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY*")


def _build_host_context(
    *,
    species: str = "human",
    sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect: Iterable[str] = (),
    accept_loss: Iterable[str] = (),
    absent: Iterable[str] = (),
) -> host_context.HostContext | None:
    """Construct a species-aware HostContext from plain-language sex and organ arguments."""
    resolved = resolve_species(species)
    anat = anatomy_for_species(resolved)
    protect_t = anat.normalize_compartments(tuple(protect))
    allow_t = anat.normalize_compartments(tuple(accept_loss))
    absent_t = anat.normalize_compartments(tuple(absent))

    prot_s, allow_s, absent_s = set(protect_t), set(allow_t), set(absent_t)
    if (prot_s & allow_s) or (prot_s & absent_s) or (allow_s & absent_s):
        raise ValueError(
            "Compartment choices conflict: a tissue cannot appear in more than one of protect, accept_loss, absent"
        )

    assumptions: tuple[str, ...] = ()
    if sex is not None:
        sex_norm = str(sex).strip().lower()
        if sex_norm in {"f", "female", "xx"}:
            sex_protected = set(anat.female_adult_protected_compartments)
            sex_absent = set(anat.male_sex_compartments)
            assumptions = ("female-patient: explicit sex='female' adult organ protection",)
        elif sex_norm in {"m", "male", "xy"}:
            if pregnant or lactating:
                raise ValueError(
                    "pregnant=True and lactating=True are incompatible with sex='male'"
                )
            sex_protected = set(anat.male_adult_protected_compartments)
            sex_absent = set(anat.female_sex_compartments)
            if resolved.taxon_id == 10090:
                sex_absent.add(anat.resting_mammary_compartment)
            assumptions = ("male-patient: explicit sex='male' adult organ protection",)
        else:
            raise ValueError(f"sex must be 'female', 'male', or None; got {sex!r}")
    else:
        sex_protected = set()
        sex_absent = set()

    eff_protect = (sex_protected | prot_s) - allow_s - absent_s
    eff_absent = (sex_absent | absent_s) - allow_s - prot_s
    eff_allow = set(allow_s)

    lact_comp = anat.lactating_mammary_compartment
    if pregnant:
        preg_comps = {"embryo", "placenta", lact_comp} - allow_s - absent_s
        eff_absent.difference_update(preg_comps)
        eff_protect.update(preg_comps)
    elif sex is not None:
        eff_absent.update({"embryo", "placenta"} - eff_protect - eff_allow)

    if lactating:
        if lact_comp not in absent_s and lact_comp not in allow_s:
            eff_absent.discard(lact_comp)
            eff_protect.add(lact_comp)
    elif (
        sex is not None
        and not pregnant
        and lact_comp not in eff_protect
        and lact_comp not in eff_allow
    ):
        eff_absent.add(lact_comp)

    if not (eff_protect or eff_absent or eff_allow or assumptions):
        return None
    return host_context.HostContext(
        protect=tuple(sorted(eff_protect)),
        accept_loss=tuple(sorted(eff_allow)),
        absent=tuple(sorted(eff_absent)),
        assumptions=assumptions,
        species=resolved.name,
    )


def table(
    group: str = "CTA",
    *,
    species: str = "human",
    sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect: Iterable[str] = (),
    accept_loss: Iterable[str] = (),
    absent: Iterable[str] = (),
    ectopic_only: bool = False,
    replicated_only: bool = False,
    min_fold: float = 10.0,
    min_patients: int = 0,
    mode: str = "all",
) -> pd.DataFrame:
    """Return the requested antigen or proteome table for ``species``.

    Parameters
    ----------
    group
        ``"CTA"``, ``"CA"``, ``"RA"``, ``"CA|RA"``, ``"CTA-M"``, ``"CTA-F"``,
        ``"CTA-P"``, ``"CTA-MF"``, ``"RA-M"``, ``"RA-F"``, ``"RA-P"``, ``"RA-MF"``,
        ``"Unsafe"``, ``"Other"``, ``"Unsafe|Other"``, or ``"all"`` (complete proteome).
    species
        ``"human"`` (9606), ``"dog"`` / ``"canine"`` (9615), or ``"mouse"`` (10090).
    sex
        ``"female"`` protects all non-pregnant, non-lactating adult female reproductive
        and resting mammary organs of ``species``; ``"male"`` protects all adult male
        reproductive and resting mammary organs (minus any ``absent`` or ``accept_loss``).
    ectopic_only
        Keep only rows with confirmed ectopic tumor activation outside homologous
        source organs (``tumor_activation_mode != "lineage_retained"``).
    replicated_only
        Keep only cancer antigens supported by ``>= 2`` independent patients or tumor samples.
    """
    resolved = resolve_species(species)
    host = _build_host_context(
        species=resolved.name,
        sex=sex,
        pregnant=pregnant,
        lactating=lactating,
        protect=protect,
        accept_loss=accept_loss,
        absent=absent,
    )
    if group == "all":
        if (
            host is not None
            or ectopic_only
            or replicated_only
            or min_fold > 10.0
            or min_patients > 0
            or mode != "all"
        ):
            raise ValueError("Candidate filters do not apply when group='all'")
        return packaged_species_results(resolved)

    if group in antigen_sets.EXCLUSION_CATEGORIES:
        if (
            host is not None
            or ectopic_only
            or replicated_only
            or min_fold > 10.0
            or min_patients > 0
            or mode != "all"
        ):
            raise ValueError("Candidate filters do not apply to Unsafe/Other backgrounds")
        return antigen_sets.proteins(group, species=resolved.name)

    overrides: dict[str, object] = {}
    if min_fold > 10.0:
        overrides["min_tumor_fold"] = float(min_fold)
    if min_patients > 0:
        overrides["min_tumor_patients"] = int(min_patients)
    sel = antigen_sets.select(group, mode=mode, host=host, species=resolved.name, **overrides)
    df = sel.proteins()
    if ectopic_only and "tumor_activation_mode" in df.columns:
        df = df[
            df["tumor_activation_mode"].isin(["ectopic", "ectopic_and_lineage_retained"])
        ].reset_index(drop=True)
    if replicated_only:
        if "tumor_support" in df.columns:
            df = df[df["tumor_support"].eq("replicated")].reset_index(drop=True)
        else:
            p_cnt = df.get("independent_patients_at_tenfold", pd.Series(0, index=df.index)).fillna(
                0
            )
            s_cnt = df.get("samples_at_tenfold", pd.Series(0, index=df.index)).fillna(0)
            df = df[(p_cnt >= 2) | (s_cnt >= 2)].reset_index(drop=True)
    return df


def gene_ids(group: str = "CTA", **kwargs) -> frozenset[str]:
    """Return member Ensembl gene IDs for ``table(group, **kwargs)``."""
    df = table(group, **kwargs)
    return frozenset(g for raw in df["member_gene_ids"] for g in json.loads(raw))


def explain(
    gene: str,
    *,
    species: str = "human",
    sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect: Iterable[str] = (),
    accept_loss: Iterable[str] = (),
    absent: Iterable[str] = (),
) -> dict[str, object]:
    """Return a plain-language, structured explanation for a gene symbol or Ensembl ID."""
    resolved = resolve_species(species)
    rows = packaged_species_results(resolved, gene)
    if rows.empty:
        raise ValueError(f"Unknown gene or protein {gene!r} in species {resolved.name!r}")
    r = rows.iloc[0]
    host = _build_host_context(
        species=resolved.name,
        sex=sex,
        pregnant=pregnant,
        lactating=lactating,
        protect=protect,
        accept_loss=accept_loss,
        absent=absent,
    )
    host_info: dict[str, object] = {"status": "not_requested", "protected_organs": ()}
    if host is not None:
        ev = packaged_species_host_exposure(resolved)
        raw_host = host_context.assess(rows, host, evidence=ev).iloc[0]
        assessed = target_filters.apply(
            replace(target_filters.TargetFilter(), host=host), rows, host_evidence=ev
        ).iloc[0]
        organ_ok = bool(raw_host["pass_host"])
        keep_ok = bool(assessed["keep"])
        eff_fold = (
            round(float(assessed["effective_tumor_fold"]), 2)
            if pd.notna(assessed.get("effective_tumor_fold"))
            else None
        )
        reasons_str = str(assessed.get("host_reasons") or "").strip(" |")
        host_info = {
            "status": str(assessed["host_status"]),
            "organ_ceiling_pass": organ_ok,
            "pass_host": bool(assessed["pass_host"]),
            "keep": keep_ok,
            "effective_tumor_fold": eff_fold,
            "protected_organs": list(host.protect),
            "limiting_compartment": str(assessed.get("host_limiting_compartment") or ""),
            "reasons": reasons_str,
        }

    reason_val = str(
        r.get("decision_reason")
        or r.get("Unsafe_reasons")
        or (
            "Qualifies both in-vivo cancer and reproductive/developmental >=10x routes"
            if bool(r.get("is_cta"))
            else str(r.get("classification", ""))
        )
    )
    out: dict[str, object] = {
        "gene": str(r["protein_group_label"]),
        "entity_id": str(r["entity_id"]),
        "member_gene_ids": json.loads(str(r["member_gene_ids"])),
        "species": resolved.name,
        "taxon_id": resolved.taxon_id,
        "category": str(r["category"]),
        "classification": str(r["classification"]),
        "sex_subset": str(r.get("sex_restriction_code") or ""),
        "activation_mode": str(r.get("tumor_activation_mode") or ""),
        "tumor_support": str(r.get("tumor_support") or ""),
        "reason": reason_val,
        "normal_baseline_tpm": (
            round(float(r["normal_reference_tpm"]), 4)
            if pd.notna(r.get("normal_reference_tpm"))
            else None
        ),
        "disposable_max_tpm": (
            round(float(r["disposable_max_tpm"]), 4)
            if pd.notna(r.get("disposable_max_tpm"))
            else None
        ),
        "tumor_peak_tpm": (
            round(float(r["tumor_max_tpm"]), 4) if pd.notna(r.get("tumor_max_tpm")) else None
        ),
        "tumor_top_cohort": str(r.get("tumor_top_cancer_code") or ""),
        "tumor_fold": round(float(r["tumor_fold"]), 2) if pd.notna(r.get("tumor_fold")) else None,
        "samples_at_10x": int(r.get("samples_at_tenfold") or 0)
        if pd.notna(r.get("samples_at_tenfold"))
        else 0,
        "patients_at_10x": int(r.get("independent_patients_at_tenfold") or 0)
        if pd.notna(r.get("independent_patients_at_tenfold"))
        else 0,
        "reproductive_peak_tpm": (
            round(float(r["reproductive_supported_tpm"]), 4)
            if pd.notna(r.get("reproductive_supported_tpm"))
            else None
        ),
        "reproductive_organs": str(r.get("reproductive_qualifying_compartments") or ""),
        "reproductive_fold": (
            round(float(r["reproductive_fold"]), 2)
            if pd.notna(r.get("reproductive_fold"))
            else None
        ),
        "cell_line_peak_tpm": (
            round(float(r["cell_line_max_tpm"]), 4)
            if pd.notna(r.get("cell_line_max_tpm"))
            else None
        ),
        "qualifying_cell_lines": str(r.get("qualifying_cell_lines") or ""),
        "host": host_info,
    }
    if resolved.taxon_id == 9606:
        from . import tcr_trials

        gs = tcr_trials.gold_standard()
        mem_set = set(out["member_gene_ids"])  # type: ignore[arg-type]
        cross_match = (
            gs["cross_recognized_gene_ids"]
            .fillna("")
            .astype(str)
            .map(lambda s: bool({g.strip() for g in s.split(";") if g.strip()} & mem_set))
        )
        gs = gs[gs["entity_id"].eq(str(r["entity_id"])) | cross_match]
        out["clinical_tcr_benchmarks"] = (
            json.loads(gs.to_json(orient="records")) if not gs.empty else []
        )
    return out


def summary(species: str = "all") -> pd.DataFrame:
    """Return a concise overview table across ``human``, ``dog``, and ``mouse`` (or one species)."""
    from .species import for_species

    targets = (
        [s.name for s in available_species()]
        if str(species).strip().lower() == "all"
        else [resolve_species(species).name]
    )
    rows = []
    for sp_name in targets:
        sp_obj = resolve_species(sp_name)
        df = packaged_species_results(sp_obj)
        ctas = df[df["is_cta"].astype(bool)]
        ect_n = (
            int(
                ctas["tumor_activation_mode"]
                .isin(["ectopic", "ectopic_and_lineage_retained"])
                .sum()
            )
            if "tumor_activation_mode" in ctas.columns
            else len(ctas)
        )
        rep_n = (
            int(ctas["tumor_support"].eq("replicated").sum())
            if "tumor_support" in ctas.columns
            else int(
                (
                    ctas.get("independent_patients_at_tenfold", pd.Series(0, index=ctas.index))
                    .fillna(0)
                    .ge(2)
                    | ctas.get("samples_at_tenfold", pd.Series(0, index=ctas.index)).fillna(0).ge(2)
                ).sum()
            )
        )
        sp_data = for_species(sp_obj)
        try:
            n_cohorts = len(sp_data.antigens.tumor_cohorts())
        except Exception:
            n_cohorts = 0
        try:
            n_lines = len(sp_data.antigens.cell_lines())
        except Exception:
            n_lines = 0
        rows.append(
            {
                "species": sp_obj.name,
                "taxon_id": sp_obj.taxon_id,
                "proteins": len(df),
                "CTA": len(ctas),
                "CTA_ectopic": ect_n,
                "CTA_replicated": rep_n,
                "CTA_M": int(ctas["sex_restriction_code"].eq("CTA-M").sum()),
                "CTA_P": int(ctas["sex_restriction_code"].eq("CTA-P").sum()),
                "CTA_F": int(ctas["sex_restriction_code"].eq("CTA-F").sum()),
                "CTA_MF": int(ctas["sex_restriction_code"].eq("CTA-MF").sum()),
                "CA_only": int(df["classification"].eq("Cancer Antigen only").sum()),
                "RA_only": int(df["classification"].eq("Reproductive Antigen only").sum()),
                "Other": int(df["category"].eq("Other").sum()),
                "Unsafe": int(df["category"].eq("Unsafe").sum()),
                "tumor_cohorts": n_cohorts,
                "cell_lines": n_lines,
            }
        )
    return pd.DataFrame(rows)


def funnel(
    species: str = "human",
    *,
    mode: str = "all",
    sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect: Iterable[str] = (),
    accept_loss: Iterable[str] = (),
    absent: Iterable[str] = (),
    ectopic_only: bool = False,
    replicated_only: bool = False,
) -> pd.DataFrame:
    """Return the sequential screening and filter attrition table for ``species``."""
    resolved = resolve_species(species)
    cat = packaged_species_results(resolved)
    host = _build_host_context(
        species=resolved.name,
        sex=sex,
        pregnant=pregnant,
        lactating=lactating,
        protect=protect,
        accept_loss=accept_loss,
        absent=absent,
    )
    filt = target_filters.FILTER_MODES[mode] if isinstance(mode, str) else mode
    if host is not None:
        from dataclasses import replace

        filt = replace(filt, host=host)
    ev = packaged_species_host_exposure(resolved) if host is not None else None
    steps = [
        {"stage": "All coding protein groups", "remaining": len(cat), "removed": 0},
        {
            "stage": "Pass somatic & immune safety (not Unsafe)",
            "remaining": int((~cat["category"].eq("Unsafe")).sum()),
            "removed": int(cat["category"].eq("Unsafe").sum()),
        },
        {
            "stage": "Pass coverage, sequence & home-niche eligibility",
            "remaining": int(
                (
                    ~cat["category"].eq("Unsafe")
                    & ~cat["classification"].isin(
                        ["Other:Review", "Other:Sequence", "Other:Unscreened tissue"]
                    )
                ).sum()
            ),
            "removed": int(
                cat["classification"]
                .isin(["Other:Review", "Other:Sequence", "Other:Unscreened tissue"])
                .sum()
            ),
        },
        {
            "stage": "Candidate antigens (CA ∪ RA >= 10x H0)",
            "remaining": int(cat["is_target_candidate"].sum()),
            "removed": int(cat["classification"].eq("Other:Expression").sum()),
        },
        {
            "stage": "Cancer-Testis Antigens (CTA = CA ∩ RA)",
            "remaining": int(cat["is_cta"].sum()),
            "removed": int((cat["is_target_candidate"] & ~cat["is_cta"]).sum()),
        },
    ]
    cur = cat[cat["is_cta"].astype(bool)].copy()
    if mode != "all" or host is not None:
        applied = target_filters.apply(filt, cur, host_evidence=ev)
        kept = applied[applied["keep"]].copy()
        steps.append(
            {
                "stage": f"Host / filter mode ({mode})",
                "remaining": len(kept),
                "removed": len(cur) - len(kept),
            }
        )
        cur = kept
    if ectopic_only and "tumor_activation_mode" in cur.columns:
        nxt = cur[cur["tumor_activation_mode"].isin(["ectopic", "ectopic_and_lineage_retained"])]
        steps.append(
            {
                "stage": "Ectopic tumor activation",
                "remaining": len(nxt),
                "removed": len(cur) - len(nxt),
            }
        )
        cur = nxt
    if replicated_only and "tumor_support" in cur.columns:
        nxt = cur[cur["tumor_support"].eq("replicated")]
        steps.append(
            {
                "stage": "Replicated across >=2 patients/samples",
                "remaining": len(nxt),
                "removed": len(cur) - len(nxt),
            }
        )
    return pd.DataFrame(steps)


def sequences(
    group: str = "CTA",
    *,
    species: str = "human",
    sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect: Iterable[str] = (),
    accept_loss: Iterable[str] = (),
    absent: Iterable[str] = (),
    sequence_source: str = "auto",
) -> pd.DataFrame:
    """Return annotated protein sequences for ``group`` under the requested host constraints."""
    resolved = resolve_species(species)
    host = _build_host_context(
        species=resolved.name,
        sex=sex,
        pregnant=pregnant,
        lactating=lactating,
        protect=protect,
        accept_loss=accept_loss,
        absent=absent,
    )
    if group in antigen_sets.EXCLUSION_CATEGORIES:
        if host is not None:
            return antigen_sets.select(
                "CA|RA", species=resolved.name, host=host, sequence_source=sequence_source
            ).background_sequences(group)
        return antigen_sets.sequences(group, species=resolved.name, sequence_source=sequence_source)
    return antigen_sets.select(
        group, species=resolved.name, host=host, sequence_source=sequence_source
    ).sequences()


def subtract(
    query: str,
    *,
    species: str = "human",
    sex: str | None = None,
    pregnant: bool = False,
    lactating: bool = False,
    protect: Iterable[str] = (),
    accept_loss: Iterable[str] = (),
    absent: Iterable[str] = (),
    background: str = "Unsafe|Other",
    max_mismatches: float = 0.0,
    lengths: Iterable[int] = range(8, 12),
    sequence_source: str = "auto",
) -> pd.DataFrame:
    """Subtract host self-proteome k-mers from a protein sequence (or gene symbol)."""
    from .peptide_subtraction import MatchPolicy, PeptideIndex

    resolved = resolve_species(species)
    lengths_list = list(lengths)
    if not lengths_list or min(lengths_list) < 1:
        raise ValueError("lengths must contain positive integers")
    q_str = str(query).strip()
    if not q_str:
        raise ValueError("query must be a non-empty gene identifier or protein sequence")

    query_records: list[tuple[str, str]] = []
    if not (set(q_str.upper()) <= _AMINO_ACIDS and len(q_str) >= 64):
        rows = packaged_species_results(resolved, q_str)
        if not rows.empty:
            seq_df = antigen_sets._sequences(
                rows, "CA|RA", species=resolved.name, sequence_source=sequence_source
            )
            seen_seqs: set[tuple[str, str]] = set()
            for srow in seq_df.itertuples(index=False):
                tid = str(getattr(srow, "transcript_id", "") or q_str)
                sval = str(srow.sequence).upper()
                if (tid, sval) not in seen_seqs:
                    seen_seqs.add((tid, sval))
                    query_records.append((tid, sval))

    if not query_records:
        q_up = q_str.upper()
        if not set(q_up) <= _AMINO_ACIDS or len(q_up) < min(lengths_list):
            raise ValueError(
                f"Query {query!r} is neither a known gene in {resolved.name!r} nor a valid "
                f"amino-acid sequence of length >= {min(lengths_list)}"
            )
        query_records.append(("query", q_up))

    host = _build_host_context(
        species=resolved.name,
        sex=sex,
        pregnant=pregnant,
        lactating=lactating,
        protect=protect,
        accept_loss=accept_loss,
        absent=absent,
    )
    policy = MatchPolicy(max_distance=float(max_mismatches))
    if host is not None:
        bg_records = antigen_sets.select(
            "CA|RA", species=resolved.name, host=host, sequence_source=sequence_source
        ).background_sequences(background)
    else:
        bg_records = antigen_sets.sequences(
            background, species=resolved.name, sequence_source=sequence_source
        )
    if bg_records.attrs.get("missing_gene_ids"):
        raise ValueError(
            "Background has missing gene translations; subtraction coverage incomplete"
        )
    idx = PeptideIndex(bg_records, policy=policy, lengths=lengths_list)
    parts = [
        idx.subtract(seq, lengths=lengths_list).assign(query_id=qid) for qid, seq in query_records
    ]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


__all__ = [
    "explain",
    "funnel",
    "gene_ids",
    "sequences",
    "subtract",
    "summary",
    "table",
]
