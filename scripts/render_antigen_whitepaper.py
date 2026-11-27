"""Render the antigen whitepaper and its evidence tables from a source checkout."""

import argparse
import json
import re
import sys
from pathlib import Path
from string import Template

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_target_filter_audit import FAMILIAR, lookup  # noqa: E402

from oncoref import (  # noqa: E402
    antigen_data,
    canine,
    catalog,
    cta_selection_policy,
    human,
    mouse,
    species_policies,
    species_target_selection,
    target_filters,
    target_selection,
    tcr_trials,
)

DOCS = ROOT / "docs"
BUILD = DOCS / "audits/target-selection-20261110/rebuild"
ASSETS = "antigen-whitepaper"
NORMAL_PEAKS = DOCS / "audits/cta-selection-20261028/normal-peak-evidence.csv.gz"
FEMALE_TUMOR_PEAKS = ROOT / "scripts/cta_evidence_assets/female-lineage-tumor-peaks.csv"
TEMPLATE = ROOT / "scripts/cta_evidence_assets/antigen-whitepaper.md"
SPECIES_TEMPLATE = ROOT / "scripts/cta_evidence_assets/species-antigen-whitepaper.md"
PEAK = ["tumor_max_tpm", "protein_group_label", "entity_id"]


def table(headers, rows):
    def cell(x):
        if pd.isna(x) or x == "":
            return "—"
        return str(x).replace("|", " / ").replace("\n", " ")

    return "\n".join(
        "| " + " | ".join(map(cell, line)) + " |"
        for line in [headers, ["---"] * len(headers), *rows]
    )


def number(x):
    return "—" if pd.isna(x) else f"{x:,.2f}".rstrip("0").rstrip(".")


def contexts(raw, field):
    observations = json.loads(raw) if pd.notna(raw) else []
    return " | ".join(dict.fromkeys(str(o[field]) for o in observations if o.get(field)))


def disposition(r):
    if r.category == "Unsafe":
        result = r.Unsafe_reasons
    elif r.subtype == "Unscreened tissue":
        result = f"Coverage gap in {r.home_tissue}"
    else:
        result = r.decision_reason
    if r.rare_cell_status == "established" or (
        r.subtype == "Review" and r.rare_cell_status == "unresolved"
    ):
        result += f" — {r.finding_tissue}, {r.finding_top_cell_type} ({r.rare_cell_basis})"
    return result


def study_table(proteins):
    counts = target_filters.paper_counts()
    data = proteins.join(counts, on="entity_id")
    data[counts.columns] = data[counts.columns].fillna(0).astype(int)
    familiar = {}
    for members in FAMILIAR.values():
        for name in members.split():
            row = lookup(data, name)
            if row is None:
                raise ValueError(f"Familiar protein missing from catalog: {name}")
            familiar.setdefault(row.entity_id, []).append(name)
    data["familiar_names"] = data.entity_id.map(
        {key: " | ".join(names) for key, names in familiar.items()}
    )
    trials = pd.read_csv(ROOT / "oncoref/data/cta-trial-evidence.csv")
    by_gene = data.set_index("entity_id").member_gene_ids.map(json.loads).explode()
    trials = trials.merge(
        by_gene.rename("Ensembl_Gene_ID").reset_index(),
        on="Ensembl_Gene_ID",
        how="left",
        validate="many_to_one",
    )
    if trials.entity_id.isna().any():
        raise ValueError("Clinical evidence with no catalog locus")
    for source, name in [("source_url", "clinical_sources"), ("evidence_kind", "clinical_scope")]:
        data[name] = data.entity_id.map(
            trials.groupby("entity_id")[source].agg(lambda x: " | ".join(sorted(set(x))))
        )
    data = data[data.papers_cta.gt(0) | data.familiar_names.notna() | data.clinical_sources.notna()]
    data = data.copy()
    data["disposition"] = [disposition(r) for r in data.itertuples()]
    # Quote the strongest excluding normal-expression measurement whose specific
    # gate and current injury tier still exclude the group.
    measurements = pd.read_csv(
        NORMAL_PEAKS,
        usecols=[
            "gene_id",
            "tissue",
            "cell_type",
            "source",
            "unit",
            "value",
            "candidate_ceiling",
            "selection_tier",
            "category",
            "taxonomy_id",
            "stage",
            "mapping_status",
            "unsafe_evidence",
            "assay",
        ],
        low_memory=False,
    )
    flags = data.loc[data.category.eq("Unsafe"), ["entity_id", "Unsafe_reasons"]]
    measured = (
        measurements.loc[
            measurements.unsafe_evidence
            & measurements.selection_tier.eq(cta_selection_policy.compartment_tiers(measurements))
        ]
        .merge(by_gene.rename("gene_id").reset_index(), on="gene_id", validate="many_to_one")
        .merge(flags, on="entity_id", validate="many_to_one")
    )
    matched = measured.loc[
        [r.assay in r.Unsafe_reasons.split(" | ") for r in measured.itertuples()]
    ].copy()
    matched["ceiling_ratio"] = matched.value / matched.candidate_ceiling
    strongest = matched.sort_values(
        ["ceiling_ratio", "source", "tissue", "cell_type"], ascending=[False, True, True, True]
    ).drop_duplicates("entity_id")
    locations = {
        r.entity_id: f"{r.tissue}, {r.cell_type} · {number(r.value)} {r.unit} "
        f"(ceiling {number(r.candidate_ceiling)})"
        for r in strongest.itertuples()
    }
    data["matched_normal_assay"] = data.entity_id.map(locations)
    return data.sort_values(
        ["papers_cta", "protein_group_label", "entity_id"], ascending=[False, True, True]
    )


def ctag2_evidence(data):
    row = data.loc[data.protein_group_label.eq("CTAG2")].iloc[0]
    path = DOCS / "audits/normal-cell-atlas-20261110/rare-cell"
    findings = pd.read_csv(path / "findings.csv")
    findings = findings[findings.entity_id.eq(row.entity_id)]
    risk = catalog.read_table("rare-cell-evidence").set_index("entity_id")
    return {
        "ctag2_cells": table(
            [
                "Organ",
                "Strongest population",
                "Independent studies",
                "Donor IDs with signal (linked where known)",
                "Cells with ≥2 molecules in top donor",
                "Grading",
            ],
            [
                [
                    r.tissue,
                    r.top_cell_type,
                    r.studies_multi_molecule,
                    r.donors_multi_molecule,
                    f"{r.top_n_multi:,}/{r.top_n_cells:,} ({r.top_multi_molecule_pct:.2f}%)",
                    r.status,
                ]
                for r in findings.itertuples()
            ],
        ),
        "ctag2_shared": f"{100 * risk.at[row.entity_id, 'shared_kmer_fraction']:.1f}%",
        "ctag2_classification": row.classification,
        "ctag2_site": contexts(row.normal_reference_contexts, "tissue"),
        "ctag2_reference": number(row.normal_reference_tpm),
        "ctag2_bar": number(row.single_tumor_threshold_tpm),
        "ctag2_repro": number(row.reproductive_supported_tpm),
        "ctag2_fold": number(row.reproductive_fold),
        "ctag2_tumor": number(row.tumor_max_tpm),
        "ctag2_patients": int(row.independent_patients_at_tenfold),
    }


def tissue_hierarchy():
    from oncoref import tissue_tiers

    sections = [
        "# Human tissue injury tiers",
        "The adult human rules describe loss of an organ or a named population. Expression, antigen presentation and exposure are assessed separately. Each cell rule applies to the organs listed in its assessment.",
    ]
    rules = tissue_tiers.assessments()
    references = tissue_tiers.references().set_index("reference_id")
    for tier in antigen_data.tiers().itertuples():
        sections.extend([f"## {tier.tier.title()}", tier.definition])
        records = []
        for rule in rules[rules.tier.eq(tier.tier)].sort_values(["level", "label"]).itertuples():
            links = ", ".join(f"[{key}]({references.at[key, 'url']})" for key in rule.references)
            records.append(
                [rule.label, rule.level, ", ".join(rule.tissues), rule.persistent_loss, links]
            )
        sections.append(
            table(
                [
                    "Structure or population",
                    "Resolution",
                    "Organ IDs",
                    "Persistent loss",
                    "Evidence",
                ],
                records,
            )
        )
    return "\n\n".join(sections) + "\n"


def build_species_sections() -> tuple[dict[str, str], dict[str, str]]:
    """Build cross-species tables and CSV exports from live species_policies and snapshots."""
    policies = species_policies.all_species_policies()
    all_bindings = pd.concat([p.dataset_table() for p in policies], ignore_index=True)

    summary_rows = []
    anat_rows = []
    for mod, pol in ((human, policies[0]), (canine, policies[1]), (mouse, policies[2])):
        cov = mod.target_selection_coverage()
        a = pol.anatomy
        active_roles = sum(
            any(b.role == r and b.status == "active_in_catalog" for b in pol.datasets)
            for r in species_policies.ROLE_IDS
        )
        counts = cov["counts"]
        if cov["selection_status"] == "evaluated":
            cat_text = (
                f"{counts['CTA']:,} CTA · {counts['Cancer Antigen only']:,} CA only · "
                f"{counts['Reproductive Antigen only']:,} RA only · "
                f"{counts['Other']:,} Other · {counts['Unsafe']:,} Unsafe"
            )
        else:
            cat_text = (
                f"{counts['Other']:,} native nominations in Other:Review "
                f"({cov['selection_status']})"
            )
        summary_rows.append(
            [
                f"{pol.species.scientific_name} (`{pol.species.name}`, {pol.species.taxon_id})",
                ", ".join(f"`{r}`" for r in pol.registered_reference_aliases),
                ", ".join(f"`{u}`" for u in pol.thresholds.comparison_units),
                f"{active_roles}/8 active",
                cat_text,
            ]
        )
        anat_rows.append(
            [
                f"{pol.species.scientific_name} (`{pol.species.taxon_id}`)",
                ", ".join(a.core_reproductive_compartments),
                ", ".join(a.male_accessory_compartments),
                ", ".join(a.female_tract_compartments),
                ", ".join(a.mammary_compartments),
                ", ".join(a.disposable_compartments),
            ]
        )

    role_comp = species_policies.role_comparison_table()
    eight_pillar_rows = [
        [
            f"**{r.role_label}**",
            f"{r.human_sources} — *{r.human_status}*",
            f"{r.dog_sources} — *{r.dog_status}*",
            f"{r.mouse_sources} — *{r.mouse_status}*",
        ]
        for r in role_comp.itertuples(index=False)
    ]

    def per_species_dataset_md(pol: species_policies.SpeciesPolicy) -> str:
        df = pol.dataset_table()
        return table(
            ["Canonical role", "Dataset(s)", "Native unit", "Status", "Threshold / screening rule"],
            [
                [
                    r.role_label,
                    f"`{r.source_id}` — {r.label}",
                    r.native_unit,
                    f"`{r.status}`",
                    r.threshold_or_use,
                ]
                for r in df.itertuples(index=False)
            ],
        )

    nom_csv_rows = []
    nom_tables: dict[str, str] = {}
    for mod, key in ((canine, "dog"), (mouse, "mouse")):
        snap = mod.cta_snapshot()
        ev_by_cand: dict[str, int] = {}
        for ev in snap.evidence():
            for cid in ev["candidate_ids"]:
                ev_by_cand[cid] = ev_by_cand.get(cid, 0) + 1
        ts = mod.target_selection().set_index("entity_id")
        md_rows = []
        for cand in snap.candidates():
            eid = f"nomination:{snap.taxon_id}:{cand['id']}"
            if eid in ts.index:
                call_str = str(ts.loc[eid].classification)
            else:
                hits = mod.target_selection(cand["id"])
                call_str = " / ".join(dict.fromkeys(hits["classification"].astype(str)))
            gids = ", ".join(cand["gene_ids"]) or "unresolved"
            n_obs = ev_by_cand.get(cand["id"], 0)
            sources = ", ".join(dict.fromkeys(d["source_id"] for d in cand["derivation"]))
            missing = "; ".join(cand["evidence_coverage"]["missing"])
            md_rows.append(
                [
                    f"`{cand['label']}`",
                    f"`{gids}`",
                    cand["mapping"]["status"],
                    sources,
                    n_obs,
                    f"`{call_str}`",
                    missing,
                ]
            )
            nom_csv_rows.append(
                {
                    "taxon_id": snap.taxon_id,
                    "species": key,
                    "candidate_id": cand["id"],
                    "label": cand["label"],
                    "gene_ids": "|".join(cand["gene_ids"]),
                    "mapping_status": cand["mapping"]["status"],
                    "mapping_complete": bool(cand["mapping"]["complete"]),
                    "nomination_sources": "|".join(
                        dict.fromkeys(d["source_id"] for d in cand["derivation"])
                    ),
                    "native_observations": n_obs,
                    "classification": call_str,
                    "missing_requirements": missing,
                }
            )
        nom_tables[key] = table(
            [
                "Candidate",
                "Reference gene ID(s)",
                "Mapping",
                "Source(s)",
                "Native obs.",
                "Catalog call",
                "Missing for quantitative CA/RA",
            ],
            md_rows,
        )

    template_vars = {
        "species_summary_table": table(
            [
                "Species (taxon)",
                "Pinned references",
                "Comparison unit(s)",
                "Evidence pillars",
                "Packaged catalog status",
            ],
            summary_rows,
        ),
        "comparative_anatomy_table": table(
            [
                "Species",
                "Core reproductive / embryonic",
                "Male accessory",
                "Female tract",
                "Mammary",
                "Non-thymic Disposable (D)",
            ],
            anat_rows,
        ),
        "eight_pillar_matrix": table(
            [
                "Screening pillar",
                "Human (9606)",
                "Domestic dog (9615)",
                "Laboratory mouse (10090)",
            ],
            eight_pillar_rows,
        ),
        "human_datasets_table": per_species_dataset_md(policies[0]),
        "dog_datasets_table": per_species_dataset_md(policies[1]),
        "mouse_datasets_table": per_species_dataset_md(policies[2]),
        "dog_nominations_table": nom_tables["dog"],
        "mouse_nominations_table": nom_tables["mouse"],
    }
    h_cat = human.target_selection()
    d_cat = canine.target_selection()
    d_sex = canine.sex_restricted_antigens()
    d_sex_cta = d_sex[d_sex["is_cta"].astype(bool)] if not d_sex.empty else d_sex
    d_loci = int(d_cat.member_gene_ids.map(json.loads).map(len).sum())
    d_ca = d_cat.is_cancer_antigen.astype(bool)
    d_ra = d_cat.is_reproductive_antigen.astype(bool)
    d_cohorts = canine.tumor_cohorts()
    d_cell_lines = canine.cell_lines()
    d_cta_df = d_cat[d_ca & d_ra]
    d_ca_only_df = d_cat[d_ca & ~d_ra]
    d_top_cta = (
        d_cta_df.dropna(subset=["tumor_max_tpm"])
        .sort_values(PEAK, ascending=[False, True, True])
        .head(15)
    )
    prame_dog = canine.target_selection("PRAME").iloc[0]

    m_cat = mouse.target_selection()
    m_sex = mouse.sex_restricted_antigens()
    m_sex_cta = m_sex[m_sex["is_cta"].astype(bool)] if not m_sex.empty else m_sex
    m_loci = int(m_cat.member_gene_ids.map(json.loads).map(len).sum())
    m_ca = m_cat.is_cancer_antigen.astype(bool)
    m_ra = m_cat.is_reproductive_antigen.astype(bool)
    m_cohorts = mouse.tumor_cohorts()
    m_cell_lines = mouse.cell_lines()
    m_cta_df = m_cat[m_ca & m_ra]
    m_ca_only_df = m_cat[m_ca & ~m_ra]
    m_top_cta = (
        m_cta_df.dropna(subset=["tumor_max_tpm"])
        .sort_values(PEAK, ascending=[False, True, True])
        .head(15)
    )
    p1a_row = mouse.target_selection("P1A").iloc[0]

    def _fmt_qcl(raw_qcl: object) -> str:
        if pd.isna(raw_qcl) or not str(raw_qcl).strip():
            return "—"
        items = [x.strip() for x in str(raw_qcl).split("|") if x.strip()]
        if not items:
            return "—"
        shown = ", ".join(f"`{x}`" for x in items[:3])
        return shown + (f" (+{len(items) - 3})" if len(items) > 3 else "")

    def _mode_vars(prefix: str, cta_sub: pd.DataFrame, ca_only_sub: pd.DataFrame) -> dict[str, str]:
        ect = int(cta_sub["tumor_activation_mode"].eq("ectopic").sum())
        both = int(cta_sub["tumor_activation_mode"].eq("ectopic_and_lineage_retained").sum())
        lin = int(cta_sub["tumor_activation_mode"].eq("lineage_retained").sum())
        any_ect = ect + both
        tot = max(len(cta_sub), 1)
        ca_ect = int(
            ca_only_sub["tumor_activation_mode"]
            .isin(["ectopic", "ectopic_and_lineage_retained"])
            .sum()
        )
        ca_lin = int(ca_only_sub["tumor_activation_mode"].eq("lineage_retained").sum())
        return {
            f"{prefix}cta_ectopic_n": f"{ect:,}",
            f"{prefix}cta_both_mode_n": f"{both:,}",
            f"{prefix}cta_ectopic_any_n": f"{any_ect:,}",
            f"{prefix}cta_ectopic_any_pct": f"{100.0 * any_ect / tot:.1f}",
            f"{prefix}cta_lineage_only_n": f"{lin:,}",
            f"{prefix}ca_only_ectopic_n": f"{ca_ect:,}",
            f"{prefix}ca_only_lineage_n": f"{ca_lin:,}",
        }

    # Build complete cross-species ortholog concordance ledger across 100% of human, dog, and mouse CTAs
    import re
    from collections import defaultdict

    from oncoref.load_dataset import get_data

    syn_df = get_data("ncbi-symbol-synonyms")
    human_official_up = set(syn_df["official_symbol"].astype(str).str.strip().str.upper()) | {
        str(lbl).split("/")[0].strip().upper()
        for lbl in h_cat["protein_group_label"]
        if not str(lbl).startswith("ENSG0")
    }
    alias_to_officials_up: dict[str, set[str]] = defaultdict(set)
    for r_syn in syn_df.itertuples(index=False):
        a_up = str(r_syn.alias).strip().upper()
        o_up = str(r_syn.official_symbol).strip().upper()
        if a_up and o_up:
            alias_to_officials_up[a_up].add(o_up)
    uniq_syn_to_human = {
        a_up: next(iter(offs))
        for a_up, offs in alias_to_officials_up.items()
        if len(offs) == 1 and a_up not in human_official_up
    }
    canine_synteny_re = re.compile(r"^C(?:\d+|X)H(\d+|X)ORF(\d+)$")
    mouse_defb_re = re.compile(r"^DEFB(\d{1,2})$")

    def _norm_orth_key(lbl: str, sp_tag: str) -> str:
        raw = str(lbl).split("/")[0].strip()
        up = raw.upper()
        if up.startswith(("ENSG0", "ENSCAFG0", "ENSMUSG0")) or re.match(r"^GM\d+$", up):
            return f"{sp_tag.upper()}:{raw}"
        # Strip murine .1/.2 tandem paralog suffix when base exists
        if sp_tag == "mouse" and re.search(r"\.\d+$", up):
            up = re.sub(r"\.\d+$", "", up)
        # Normalize canine syntenic ORF symbols C{dog}H{human}orf{n} -> C{human}ORF{n}
        m_syn = canine_synteny_re.match(up)
        if m_syn:
            up = f"C{m_syn.group(1)}ORF{m_syn.group(2)}"
        # Normalize murine beta-defensins Defb{n} -> DEFB1{n:02d} when human ortholog exists
        if sp_tag == "mouse":
            m_db = mouse_defb_re.match(up)
            if m_db:
                cand_db = f"DEFB1{int(m_db.group(1)):02d}"
                if cand_db in human_official_up:
                    up = cand_db
        explicit_map = {
            "NY-ESO-1": "CTAG1",
            "CTAG1A": "CTAG1",
            "CTAG1B": "CTAG1",
            "TRAP1A": "TRAP1A",
        }
        if up in explicit_map:
            return explicit_map[up]
        return uniq_syn_to_human.get(up, up)

    def _index_by_orth(df: pd.DataFrame, sp_tag: str) -> dict[str, pd.Series]:
        out_idx: dict[str, pd.Series] = {}
        rank_cat = {
            "CTA": 0,
            "Reproductive Antigen only": 1,
            "Cancer Antigen only": 2,
            "Other": 3,
            "Unsafe": 4,
        }
        for _, r in df.iterrows():
            lbl = str(r["protein_group_label"])
            k = _norm_orth_key(lbl, sp_tag)
            prev = out_idx.get(k)
            if prev is None or rank_cat.get(str(r["category"]), 9) < rank_cat.get(
                str(prev["category"]), 9
            ):
                out_idx[k] = r
        return out_idx

    h_orth = _index_by_orth(h_cat, "human")
    d_orth = _index_by_orth(d_cat, "dog")
    m_orth = _index_by_orth(m_cat, "mouse")
    cta_keys = sorted(
        {
            *(k for k, r in h_orth.items() if bool(r["is_cta"])),
            *(k for k, r in d_orth.items() if bool(r["is_cta"])),
            *(k for k, r in m_orth.items() if bool(r["is_cta"])),
        }
    )
    orth_rows: list[dict] = []
    for k in cta_keys:
        hr, dr, mr = h_orth.get(k), d_orth.get(k), m_orth.get(k)
        cta_sp = [
            sp
            for sp, rr in (("human", hr), ("dog", dr), ("mouse", mr))
            if rr is not None and bool(rr["is_cta"])
        ]
        ra_or_cta_sp = [
            sp
            for sp, rr in (("human", hr), ("dog", dr), ("mouse", mr))
            if rr is not None and bool(rr["is_reproductive_antigen"])
        ]
        present_sp = [
            sp for sp, rr in (("human", hr), ("dog", dr), ("mouse", mr)) if rr is not None
        ]
        if len(cta_sp) >= 2:
            regime = "Conserved multi-species CTA"
        elif len(cta_sp) == 1 and len(ra_or_cta_sp) >= 2:
            regime = "CTA + conserved RA-only ortholog"
        elif len(present_sp) >= 2:
            regime = "Divergent somatic / review gate"
        else:
            regime = "Lineage-specific / no 1:1 symbol ortholog"
        orth_rows.append(
            {
                "ortholog_key": k,
                "cta_species": "|".join(cta_sp),
                "human_symbol": str(hr["protein_group_label"]) if hr is not None else "",
                "human_classification": str(hr["classification"]) if hr is not None else "absent",
                "human_sex_code": str(hr.get("sex_restriction_code", "") or "")
                if hr is not None
                else "",
                "dog_symbol": str(dr["protein_group_label"]) if dr is not None else "",
                "dog_classification": str(dr["classification"]) if dr is not None else "absent",
                "dog_sex_code": str(dr.get("sex_restriction_code", "") or "")
                if dr is not None
                else "",
                "mouse_symbol": str(mr["protein_group_label"]) if mr is not None else "",
                "mouse_classification": str(mr["classification"]) if mr is not None else "absent",
                "mouse_sex_code": str(mr.get("sex_restriction_code", "") or "")
                if mr is not None
                else "",
                "concordance_regime": regime,
            }
        )
    orth_df = pd.DataFrame(orth_rows)
    h_cta_orth = orth_df[orth_df["human_classification"].eq("CTA")]
    conserved_cta_rows = orth_df[orth_df["concordance_regime"].eq("Conserved multi-species CTA")]
    conserved_cta_examples = ", ".join(
        f"`{r.human_symbol or r.dog_symbol}`"
        + (
            f"/`{r.mouse_symbol}`"
            if r.mouse_classification == "CTA"
            and r.mouse_symbol
            and r.mouse_symbol != (r.human_symbol or r.dog_symbol)
            else ""
        )
        for r in conserved_cta_rows.itertuples(index=False)
    )
    conserved_sub = orth_df[
        orth_df["concordance_regime"].isin(
            ["Conserved multi-species CTA", "CTA + conserved RA-only ortholog"]
        )
    ]

    template_vars.update(
        {
            "dog_n": f"{len(d_cat):,}",
            "dog_loci": f"{d_loci:,}",
            "dog_cta": f"{int((d_ca & d_ra).sum()):,}",
            "dog_ca_only": f"{int((d_ca & ~d_ra).sum()):,}",
            "dog_ra": f"{int(d_ra.sum()):,}",
            "dog_ra_only": f"{int((d_ra & ~d_ca).sum()):,}",
            "dog_union": f"{int((d_ca | d_ra).sum()):,}",
            "dog_other_n": f"{int(d_cat.category.eq('Other').sum()):,}",
            "dog_unsafe_n": f"{int(d_cat.category.eq('Unsafe').sum()):,}",
            "dog_mra_n": f"{int(d_sex.sex_restriction_class.eq('Male-restricted').sum()):,}",
            "dog_fra_n": f"{int(d_sex.sex_restriction_class.eq('Female-restricted').sum()):,}",
            "dog_era_n": f"{int(d_sex.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "dog_sra_n": f"{int(d_sex.sex_restriction_class.eq('Shared male & female').sum()):,}",
            "dog_cta_m_n": f"{int(d_sex_cta.sex_restriction_class.eq('Male-restricted').sum()):,}",
            "dog_cta_f_n": f"{int(d_sex_cta.sex_restriction_class.eq('Female-restricted').sum()):,}",
            "dog_cta_e_n": f"{int(d_sex_cta.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "dog_cta_mf_n": f"{int(d_sex_cta.sex_restriction_class.eq('Shared male & female').sum()):,}",
            "dog_cta_p_n": f"{int(d_sex_cta.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "dog_ra_p_n": f"{int(d_sex.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "dog_prame_call": prame_dog.classification,
            "dog_cohorts_n": f"{len(d_cohorts):,}",
            "dog_tumor_samples_n": f"{int(d_cohorts['n_samples'].sum()):,}",
            "dog_cell_lines_n": f"{len(d_cell_lines):,}",
            "dog_cta_in_cell_lines_n": f"{int(d_cta_df['qualifying_cell_lines'].fillna('').astype(str).str.strip().ne('').sum()):,}",
            **_mode_vars("dog_", d_cta_df, d_ca_only_df),
            "dog_cohorts_table": table(
                ["Cohort", "Model", "Breeds / DLA", "Accession", "Samples", "Quantification"],
                [
                    [
                        f"`{r.cancer_code}` — {r.label}",
                        r.model_kind,
                        f"{r.strain} (`{r.mhc_haplotype}`)",
                        f"`{r.source_accession}`",
                        r.n_samples,
                        r.quantification,
                    ]
                    for r in d_cohorts.itertuples(index=False)
                ],
            ),
            "dog_cell_lines_table": table(
                [
                    "Cell line",
                    "Histology",
                    "Origin / breed",
                    "Accession",
                    "Samples",
                    "Quantification",
                ],
                [
                    [
                        f"`{r.model_id}` — {r.label}",
                        f"`{r.cancer_code}`",
                        f"{r.strain} (`{r.mhc_haplotype}`)",
                        f"`{r.source_accession}`",
                        r.n_samples,
                        r.quantification,
                    ]
                    for r in d_cell_lines.itertuples(index=False)
                ],
            ),
            "dog_top_cta_table": table(
                [
                    "Rank",
                    "Protein",
                    "Sex subset",
                    "Mode",
                    "Peak cohort",
                    "Tumor TPM",
                    "N TPM",
                    "RA peak TPM",
                    "Cell lines ≥10×",
                    "Qualifying reproductive compartments",
                ],
                [
                    [
                        i,
                        f"`{r.protein_group_label}`",
                        f"`{r.sex_restriction_code}`",
                        f"`{r.tumor_activation_mode}`",
                        f"`{r.tumor_top_cancer_code}`",
                        number(r.tumor_max_tpm),
                        number(r.normal_reference_tpm),
                        number(r.reproductive_supported_tpm),
                        _fmt_qcl(r.qualifying_cell_lines),
                        r.reproductive_qualifying_compartments,
                    ]
                    for i, r in enumerate(d_top_cta.itertuples(index=False), 1)
                ],
            ),
            "mouse_n": f"{len(m_cat):,}",
            "mouse_loci": f"{m_loci:,}",
            "mouse_cta": f"{int((m_ca & m_ra).sum()):,}",
            "mouse_ca_only": f"{int((m_ca & ~m_ra).sum()):,}",
            "mouse_ra": f"{int(m_ra.sum()):,}",
            "mouse_ra_only": f"{int((m_ra & ~m_ca).sum()):,}",
            "mouse_union": f"{int((m_ca | m_ra).sum()):,}",
            "mouse_other_n": f"{int(m_cat.category.eq('Other').sum()):,}",
            "mouse_unsafe_n": f"{int(m_cat.category.eq('Unsafe').sum()):,}",
            "mouse_mra_n": f"{int(m_sex.sex_restriction_class.eq('Male-restricted').sum()):,}",
            "mouse_fra_n": f"{int(m_sex.sex_restriction_class.eq('Female-restricted').sum()):,}",
            "mouse_era_n": f"{int(m_sex.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "mouse_sra_n": f"{int(m_sex.sex_restriction_class.eq('Shared male & female').sum()):,}",
            "mouse_cta_m_n": f"{int(m_sex_cta.sex_restriction_class.eq('Male-restricted').sum()):,}",
            "mouse_cta_f_n": f"{int(m_sex_cta.sex_restriction_class.eq('Female-restricted').sum()):,}",
            "mouse_cta_e_n": f"{int(m_sex_cta.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "mouse_cta_mf_n": f"{int(m_sex_cta.sex_restriction_class.eq('Shared male & female').sum()):,}",
            "mouse_cta_p_n": f"{int(m_sex_cta.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "mouse_ra_p_n": f"{int(m_sex.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
            "mouse_p1a_call": p1a_row.classification,
            "mouse_cohorts_n": f"{len(m_cohorts):,}",
            "mouse_tumor_samples_n": f"{int(m_cohorts['n_samples'].sum()):,}",
            "mouse_cell_lines_n": f"{len(m_cell_lines):,}",
            "mouse_cta_in_cell_lines_n": f"{int(m_cta_df['qualifying_cell_lines'].fillna('').astype(str).str.strip().ne('').sum()):,}",
            **_mode_vars("mouse_", m_cta_df, m_ca_only_df),
            "mouse_cohorts_table": table(
                ["Cohort", "Model", "Host strain (MHC)", "Accession", "Samples", "Quantification"],
                [
                    [
                        f"`{r.cancer_code}` — {r.label}",
                        r.model_kind,
                        f"{r.strain} (`{r.mhc_haplotype}`)",
                        f"`{r.source_accession}`",
                        r.n_samples,
                        r.quantification,
                    ]
                    for r in m_cohorts.itertuples(index=False)
                ],
            ),
            "mouse_cell_lines_table": table(
                [
                    "Cell line",
                    "Histology",
                    "Strain (MHC)",
                    "Accession",
                    "Samples",
                    "Quantification",
                ],
                [
                    [
                        f"`{r.model_id}` — {r.label}",
                        f"`{r.cancer_code}`",
                        f"{r.strain} (`{r.mhc_haplotype}`)",
                        f"`{r.source_accession}`",
                        r.n_samples,
                        r.quantification,
                    ]
                    for r in m_cell_lines.itertuples(index=False)
                ],
            ),
            "mouse_top_cta_table": table(
                [
                    "Rank",
                    "Protein",
                    "Sex subset",
                    "Mode",
                    "Peak cohort",
                    "Tumor TPM",
                    "N TPM",
                    "RA peak TPM",
                    "Cell lines ≥10×",
                    "Qualifying reproductive compartments",
                ],
                [
                    [
                        i,
                        f"`{r.protein_group_label}`",
                        f"`{r.sex_restriction_code}`",
                        f"`{r.tumor_activation_mode}`",
                        f"`{r.tumor_top_cancer_code}`",
                        number(r.tumor_max_tpm),
                        number(r.normal_reference_tpm),
                        number(r.reproductive_supported_tpm),
                        _fmt_qcl(r.qualifying_cell_lines),
                        r.reproductive_qualifying_compartments,
                    ]
                    for i, r in enumerate(m_top_cta.itertuples(index=False), 1)
                ],
            ),
            "ortholog_total_keys_n": f"{len(orth_df):,}",
            "human_cta_with_dog_orth_n": f"{int(h_cta_orth['dog_classification'].ne('absent').sum()):,}",
            "human_cta_with_mouse_orth_n": f"{int(h_cta_orth['mouse_classification'].ne('absent').sum()):,}",
            "conserved_cta_pairs_n": f"{len(conserved_cta_rows):,}",
            "conserved_cta_examples": conserved_cta_examples,
            "conserved_cta_or_ra_n": f"{len(conserved_sub):,}",
            "ortholog_summary_table": table(
                [
                    "Cross-species orthology regime (union of human, dog & mouse CTAs)",
                    "Ortholog families",
                    "Biological mechanism",
                ],
                [
                    [
                        "**Conserved multi-species CTA**",
                        len(conserved_cta_rows),
                        f"Independently passes Stage 1 somatic safety AND 10× cancer + reproductive enrichment in ≥2 species ({conserved_cta_examples})",
                    ],
                    [
                        "**CTA + conserved RA-only ortholog**",
                        int(
                            orth_df["concordance_regime"]
                            .eq("CTA + conserved RA-only ortholog")
                            .sum()
                        ),
                        "Conserved germline/reproductive restriction (`RA`) in ≥2 species, with ≥10× tumor reactivation observed in the broader cohort panel of one species (`DPPA5`, `MAGEB18`, `CYLC1`, `ACTL8`, `ADAM2`, `AKAP4`)",
                    ],
                    [
                        "**Divergent somatic / review gate**",
                        int(
                            orth_df["concordance_regime"]
                            .eq("Divergent somatic / review gate")
                            .sum()
                        ),
                        "1:1 ortholog exists across species, but somatic expression diverges evolutionary or is caught by deeper human single-cell/IHC screens (`HORMAD1`, `DDX43`, `CLDN6`, `IGF2BP1`)",
                    ],
                    [
                        "**Lineage-specific / paralog-diverged CTA**",
                        int(
                            orth_df["concordance_regime"]
                            .eq("Lineage-specific / no 1:1 symbol ortholog")
                            .sum()
                        ),
                        "Primate-amplified CT-X/subtelomeric families (`NY-ESO-1`, `SPANX*`, `VCX/Y`, `GAGE*`, `XAGE*`, `CT45A*`, `CT47A*`, `PRAMEF*`, `DUX4`), canine-specific loci, or murine `Gm*`/`Pramel*`/`Mageb*`/`Pbsn`/`Sbp` expansions",
                    ],
                ],
            ),
            "conserved_ortholog_ctas_table": table(
                [
                    "Ortholog family",
                    "Human call (subset)",
                    "Dog call (subset)",
                    "Mouse call (subset)",
                    "Regime",
                ],
                [
                    [
                        f"`{r.ortholog_key}`",
                        f"`{r.human_classification}`"
                        + (f" (`{r.human_sex_code}`)" if r.human_sex_code else ""),
                        f"`{r.dog_classification}`"
                        + (f" (`{r.dog_sex_code}`)" if r.dog_sex_code else ""),
                        f"`{r.mouse_classification}`"
                        + (f" (`{r.mouse_sex_code}`)" if r.mouse_sex_code else ""),
                        r.concordance_regime,
                    ]
                    for r in conserved_sub.itertuples(index=False)
                ],
            ),
        }
    )
    csv_outputs = {
        f"{ASSETS}/species-datasets.csv": all_bindings.to_csv(index=False, lineterminator="\n"),
        f"{ASSETS}/species-nominations.csv": pd.DataFrame(nom_csv_rows).to_csv(
            index=False, lineterminator="\n"
        ),
        f"{ASSETS}/cross-species-orthologs.csv": orth_df.to_csv(index=False, lineterminator="\n"),
    }
    return template_vars, csv_outputs


FEMALE_LINEAGE_PROGRAMS: list[tuple[str, list[str]]] = [
    (
        "Mammary & lactation",
        [
            "LALBA",
            "CSN2",
            "BTN1A1",
            "CSN1S1",
            "CSN3",
            "SCGB2A2",
            "SCGB1D2",
            "SCGB2A1",
            "SCGB1D4",
            "PIP",
            "MUCL1",
            "ANKRD30A",
            "ANKRD30B",
            "PIGR",
            "LTF",
            "LPO",
            "CEL",
            "XDH",
            "MFGE8",
            "PRLR",
            "OXTR",
            "ELF5",
            "TRPS1",
            "GATA3",
            "ESR1",
            "PGR",
            "VTCN1",
            "FOLR1",
            "SLC34A2",
            "ROPN1",
            "ROPN1B",
            "CT83",
        ],
    ),
    (
        "Oocyte, ZP, SCMC & follicle",
        [
            "ZP1",
            "ZP2",
            "ZP3",
            "ZP4",
            "GDF9",
            "BMP15",
            "MOS",
            "FIGLA",
            "NOBOX",
            "LHX8",
            "SOHLH1",
            "SOHLH2",
            "FOXL2",
            "AMH",
            "AMHR2",
            "FSHR",
            "LHCGR",
            "CYP19A1",
            "ASTL",
            "FETUB",
            "H1-8",
            "TUBB8",
            "WEE2",
            "PATL2",
            "BTG4",
            "PABPN1L",
            "ZAR1",
            "ZAR1L",
            "NLRP5",
            "OOEP",
            "TLE6",
            "KHDC3L",
            "PADI6",
            "NLRP2",
            "NLRP4",
            "NLRP7",
            "NLRP8",
            "NLRP9",
            "NLRP11",
            "NLRP13",
            "NLRP14",
            "OOSP1",
            "OOSP2",
            "OOSP4A",
            "DPPA3",
            "KPNA7",
            "NPM2",
            "ZBED3",
            "TCL1B",
        ],
    ),
    (
        "Shared meiosis & germ cell",
        [
            "HORMAD1",
            "HORMAD2",
            "SPO11",
            "SYCP1",
            "SYCP2",
            "SYCP3",
            "SMC1B",
            "STAG3",
            "REC8",
            "MEIOC",
            "STRA8",
            "DAZL",
            "DDX4",
            "MAEL",
            "PIWIL1",
            "PIWIL2",
            "PIWIL3",
            "PIWIL4",
            "TDRD1",
            "TDRD9",
            "BOLL",
        ],
    ),
    (
        "Endometrium, oviduct & Müllerian",
        [
            "OVGP1",
            "PAEP",
            "MMP26",
            "MMP7",
            "MUC16",
            "WFDC2",
            "MSLN",
            "PAX8",
            "LEFTY2",
            "PRL",
            "IGFBP1",
            "LIF",
            "HOXA10",
            "HOXA11",
            "WNT7A",
            "CLDN6",
            "PRAME",
            "TACSTD2",
        ],
    ),
    (
        "Placenta & trophoblast",
        [
            "ERVV-1",
            "ERVV-2",
            "ERVW-1",
            "ERVFRD-1",
            "PLAC1",
            "PLAC8",
            "PSG1",
            "PSG2",
            "PSG3",
            "PSG4",
            "PSG5",
            "PSG6",
            "PSG7",
            "PSG8",
            "PSG9",
            "PSG11",
            "CGB3/5/8",
            "CGB1",
            "CGB2",
            "CGB7",
            "CGA",
            "CSHL1",
            "GH2",
            "GCM1",
            "PAPPA",
            "PAPPA2",
            "KISS1",
            "LGALS13",
            "LGALS14",
            "LGALS16",
            "INSL4",
            "ALPP",
            "HLA-G",
            "VGLL1",
            "XAGE2",
            "XAGE3",
            "PAGE4",
            "MAGEA10",
            "HSD3B1",
        ],
    ),
]

FEMALE_HIGHLIGHT_SYMBOLS: tuple[str, ...] = (
    "CSN2",
    "BTN1A1",
    "LALBA",
    "CSN1S1",
    "CSN3",
    "SCGB2A2",
    "PIP",
    "MUCL1",
    "ANKRD30A",
    "ROPN1",
    "ELF5",
    "ZP4",
    "ZP1",
    "ZP2",
    "ZP3",
    "BMP15",
    "WEE2",
    "PABPN1L",
    "KHDC3L",
    "ZAR1L",
    "FIGLA",
    "PADI6",
    "NLRP5",
    "TUBB8",
    "DPPA3",
    "OOSP2",
    "NLRP4",
    "HORMAD1",
    "SYCP1",
    "SCGB1D4",
    "PRAME",
    "OVGP1",
    "PAEP",
    "MUC16",
    "WFDC2",
    "PAX8",
    "FOLR1",
    "SLC34A2",
    "CLDN6",
    "ERVV-1",
    "ERVV-2",
    "ERVW-1",
    "PLAC1",
    "CGA",
    "CGB3/5/8",
    "XAGE3",
    "CT83",
)


def build_sex_restriction_interrogation(
    data: pd.DataFrame, sex_table: pd.DataFrame
) -> tuple[dict[str, str], pd.DataFrame]:
    """Build symmetric RA-*/CTA-* tables and the 159-gene female-lineage audit."""
    he = species_target_selection.packaged_species_host_exposure("human")
    by_lbl = data.set_index("protein_group_label", drop=False)
    sex_by_eid = sex_table.set_index("entity_id")

    gyn_peaks_df = (
        pd.read_csv(FEMALE_TUMOR_PEAKS).drop_duplicates("entity_id").set_index("entity_id")
    )
    gyn_col_map = {
        "BRCA": "brca_max_tpm",
        "BRCA_TNBC": "brca_tnbc_max_tpm",
        "BRCA_Basal": "brca_basal_max_tpm",
        "OV": "ov_max_tpm",
        "UCEC": "ucec_max_tpm",
        "UCS": "ucs_max_tpm",
        "CESC": "cesc_max_tpm",
    }

    # Strongest Unsafe normal-tissue assay signal across all candidate member loci
    by_gene = data.set_index("entity_id").member_gene_ids.map(json.loads).explode()
    measurements = pd.read_csv(
        NORMAL_PEAKS,
        usecols=[
            "gene_id",
            "tissue",
            "cell_type",
            "source",
            "unit",
            "value",
            "candidate_ceiling",
            "selection_tier",
            "category",
            "taxonomy_id",
            "stage",
            "mapping_status",
            "unsafe_evidence",
            "assay",
        ],
        low_memory=False,
    )
    cand_eids = {
        by_lbl.loc[g].iloc[0].entity_id
        if isinstance(by_lbl.loc[g], pd.DataFrame)
        else by_lbl.loc[g].entity_id
        for _, genes in FEMALE_LINEAGE_PROGRAMS
        for g in genes
        if g in by_lbl.index
    }
    flags = data.loc[
        data.entity_id.isin(cand_eids) & data.category.eq("Unsafe"), ["entity_id", "Unsafe_reasons"]
    ]
    measured = (
        measurements.loc[
            measurements.unsafe_evidence
            & measurements.selection_tier.eq(cta_selection_policy.compartment_tiers(measurements))
        ]
        .merge(by_gene.rename("gene_id").reset_index(), on="gene_id", validate="many_to_one")
        .merge(flags, on="entity_id", validate="many_to_one")
    )
    matched = measured.loc[
        [r.assay in r.Unsafe_reasons.split(" | ") for r in measured.itertuples()]
    ].copy()
    matched["ceiling_ratio"] = matched.value / matched.candidate_ceiling
    strongest = matched.sort_values(
        ["ceiling_ratio", "source", "tissue", "cell_type"], ascending=[False, True, True, True]
    ).drop_duplicates("entity_id")
    strongest_hit_by_eid = {
        r.entity_id: f"{r.tissue}/{r.cell_type}: {number(r.value)} {r.unit} (ceil {number(r.candidate_ceiling)})"
        for r in strongest.itertuples()
    }

    # Build 159-gene audit rows + 5-program summary counts
    audit_rows: list[dict] = []
    prog_summary_rows: list[list] = []
    for prog_name, genes in FEMALE_LINEAGE_PROGRAMS:
        c_unsafe = c_other = c_m = c_mf = c_e = c_f = 0
        f_names: list[str] = []
        for sym in genes:
            if sym not in by_lbl.index:
                continue
            r = by_lbl.loc[sym]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            eid = str(r.entity_id)
            sr = sex_by_eid.loc[eid] if eid in sex_by_eid.index else None
            sex_code = str(sr["sex_restriction_code"]) if sr is not None else ""
            sex_class = str(sr["sex_restriction_class"]) if sr is not None else ""
            e_he = he.get(eid, {})
            rb = e_he.get("resting_breast", {})
            rb_xx = rb.get("female_max_tpm")
            rb_xy = rb.get("male_max_tpm")
            testis_ntpm = e_he.get("testis", {}).get("max_ntpm")
            ovary_ntpm = e_he.get("ovary", {}).get("max_ntpm")

            if r.category == "Unsafe":
                c_unsafe += 1
            elif sr is None:
                c_other += 1
            elif sex_class == "Male-restricted":
                c_m += 1
            elif sex_class == "Shared male & female":
                c_mf += 1
            elif sex_class == "Pregnancy / peripartum-only":
                c_e += 1
            else:
                c_f += 1
                f_names.append(f"`{sym}` ({sex_code})")

            gyn_vals: dict[str, float] = {}
            p_row = gyn_peaks_df.loc[eid] if eid in gyn_peaks_df.index else None
            if p_row is not None:
                for c_code, col_name in gyn_col_map.items():
                    val = pd.to_numeric(p_row.get(col_name), errors="coerce")
                    if pd.notna(val) and float(val) > 0:
                        gyn_vals[c_code] = float(val)
            pan_peak_val = (
                float(r.tumor_max_tpm)
                if pd.notna(r.tumor_max_tpm)
                else (
                    float(p_row["pan_max_tpm"])
                    if p_row is not None and pd.notna(p_row["pan_max_tpm"])
                    else float("nan")
                )
            )
            pan_peak_cohort = (
                str(r.tumor_top_cancer_code)
                if pd.notna(r.tumor_max_tpm) and str(r.tumor_top_cancer_code)
                else (
                    str(p_row["pan_top_cohort"])
                    if p_row is not None and pd.notna(p_row["pan_top_cohort"])
                    else ""
                )
            )
            top_gyn_str = ", ".join(
                f"{k}:{number(v)}"
                for k, v in sorted(gyn_vals.items(), key=lambda kv: kv[1], reverse=True)[:3]
            )

            if r.category == "Unsafe":
                hit = strongest_hit_by_eid.get(eid, "")
                if not hit and getattr(r, "rare_cell_status", "") == "established":
                    hit = f"{r.finding_tissue}/{r.finding_top_cell_type}: {r.finding_top_multi_molecule_pct:.1f}% ≥2 UMI"
                why = f"Unsafe ({r.Unsafe_reasons})" + (f" — {hit}" if hit else "")
            elif sr is None:
                if r.subtype == "Review" and getattr(r, "rare_cell_status", "") == "unresolved":
                    why = f"{r.classification} — {r.finding_tissue}/{r.finding_top_cell_type} paralog/snRNA review"
                else:
                    why = f"{r.classification} — {r.decision_reason}"
            elif sex_class == "Female-restricted":
                if bool(r.is_cta):
                    why = (
                        f"Qualifies CTA-F via {r.reproductive_qualifying_compartments}; "
                        f"male repro max {number(sr['male_sex_organ_max_tpm'])} TPM, "
                        f"XY resting breast {number(rb_xy)} TPM (<8 / <16 IHC−)"
                    )
                else:
                    why = (
                        f"RA-F-only via {r.reproductive_qualifying_compartments}; "
                        f"no somatic tumor reaches 10× H₀ (peak {number(pan_peak_val)} TPM in {pan_peak_cohort or '—'})"
                    )
            elif sex_class == "Pregnancy / peripartum-only":
                if bool(r.is_cta):
                    why = f"Adult-absent in non-pregnant/non-lactating men & women (CTA-P via {r.reproductive_qualifying_compartments})"
                else:
                    why = (
                        f"RA-P-only via {r.reproductive_qualifying_compartments}; "
                        f"adult-absent, tumor peak {number(pan_peak_val)} TPM (<10× H₀)"
                    )
            elif sex_class == "Male-restricted":
                why = (
                    f"Active in male repro ({number(sr['male_sex_organ_max_tpm'])} TPM / "
                    f"{number(testis_ntpm)} nTPM testis) while non-pregnant adult female <8 / <16 IHC− → {sex_code}"
                )
            else:
                why = (
                    f"Shared male ({number(sr['male_sex_organ_max_tpm'])} TPM / {number(testis_ntpm)} nTPM testis; "
                    f"XY breast {number(rb_xy)} TPM) & female ({number(sr['female_sex_organ_max_tpm'])} TPM / "
                    f"{number(ovary_ntpm)} nTPM ovary) → {sex_code}"
                )

            audit_rows.append(
                {
                    "program": prog_name,
                    "gene_symbol": sym,
                    "entity_id": eid,
                    "classification": r.classification,
                    "sex_restriction_code": sex_code,
                    "normal_reference_tpm": r.normal_reference_tpm,
                    "reproductive_supported_tpm": r.reproductive_supported_tpm,
                    "reproductive_qualifying_compartments": r.reproductive_qualifying_compartments,
                    "male_sex_organ_max_tpm": sr["male_sex_organ_max_tpm"]
                    if sr is not None
                    else float("nan"),
                    "testis_sc_max_ntpm": testis_ntpm,
                    "female_sex_organ_max_tpm": sr["female_sex_organ_max_tpm"]
                    if sr is not None
                    else float("nan"),
                    "ovary_sc_max_ntpm": ovary_ntpm,
                    "resting_breast_female_tpm": rb_xx,
                    "resting_breast_male_tpm": rb_xy,
                    "pan_tumor_max_tpm": pan_peak_val,
                    "pan_tumor_top_cohort": pan_peak_cohort,
                    "brca_max_tpm": gyn_vals.get("BRCA"),
                    "brca_tnbc_max_tpm": gyn_vals.get("BRCA_TNBC"),
                    "ov_max_tpm": gyn_vals.get("OV"),
                    "ucec_max_tpm": gyn_vals.get("UCEC"),
                    "ucs_max_tpm": gyn_vals.get("UCS"),
                    "cesc_max_tpm": gyn_vals.get("CESC"),
                    "top_brca_gyn_cohorts": top_gyn_str,
                    "disposition_mechanism": why,
                }
            )

        tot_p = c_unsafe + c_other + c_m + c_mf + c_e + c_f
        prog_summary_rows.append(
            [
                f"**{prog_name}**",
                tot_p,
                f"{c_unsafe} ({100 * c_unsafe / max(tot_p, 1):.0f}%)",
                c_other,
                c_m,
                c_mf,
                c_e,
                f"**{c_f}** ({', '.join(f_names) if f_names else 'none'})",
            ]
        )

    audit_df = pd.DataFrame(audit_rows)

    # 1. Symmetric 4-way mutually exclusive RA-{M,P,F,MF} / CTA-{M,P,F,MF} physiological partition
    sex_cta = sex_table[sex_table["is_cta"].astype(bool)]
    sex_ra_only = sex_table[~sex_table["is_cta"].astype(bool)]
    peri_cta_ids = species_target_selection.peripartum_antigen_ids(
        data[data.is_cta.astype(bool)], species="human"
    )
    peri_ra_ids = species_target_selection.peripartum_antigen_ids(data, species="human")
    cta_p_eids = set(
        sex_cta.loc[sex_cta.sex_restriction_class.eq("Pregnancy / peripartum-only"), "entity_id"]
    )
    ra_p_eids = set(
        sex_table.loc[
            sex_table.sex_restriction_class.eq("Pregnancy / peripartum-only"), "entity_id"
        ]
    )
    peri_cta_n = len(peri_cta_ids)
    peri_ra_n = len(peri_ra_ids)
    peri_dual_cta_n = len(peri_cta_ids - cta_p_eids)
    peri_dual_ra_n = len(peri_ra_ids - ra_p_eids)
    partition_rows = [
        [
            "**Male-restricted**",
            f"`RA-M` (`MRA`): **{int(sex_table.sex_restriction_class.eq('Male-restricted').sum()):,}**",
            f"`CTA-M`: **{int(sex_cta.sex_restriction_class.eq('Male-restricted').sum()):,}**",
            f"`RA-M-only`: **{int(sex_ra_only.sex_restriction_class.eq('Male-restricted').sum()):,}**",
            "Clears all non-pregnant, non-lactating female adult organs (`ovary`, female tract, `XX resting_breast` < 8; < 16 with `IHC−`)",
            "`NY-ESO-1`, `MAGEA4`, `XAGE1A/B`, `VCY/1B`, `CT83`, `SEMG2`, `PLAC1`, `XAGE3`, `DKKL1`, `SYCP1`, `KRTAP20-4`",
        ],
        [
            "**Pregnancy / peripartum-only**",
            f"`RA-P` (`PRA`/`RA-E`): **{int(sex_table.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}**",
            f"`CTA-P` (`CTA-E`): **{int(sex_cta.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}**",
            f"`RA-P-only`: **{int(sex_ra_only.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}**",
            "Clears all non-pregnant, non-lactating adult reproductive and resting-mammary organs in **both** sexes (< 8; < 16 with `IHC−`)",
            "`DUX4`, `ERVV-1`, `ERVV-2`, `LGALS16`, `BTN1A1`, `CER1`, `LEUTX`, `PSG7` (`RA-P-only`), `NLRP8` (`RA-P-only`)",
        ],
        [
            "**Female-restricted**",
            f"`RA-F` (`FRA`): **{int(sex_table.sex_restriction_class.eq('Female-restricted').sum()):,}**",
            f"`CTA-F`: **{int(sex_cta.sex_restriction_class.eq('Female-restricted').sum()):,}**",
            f"`RA-F-only`: **{int(sex_ra_only.sex_restriction_class.eq('Female-restricted').sum()):,}**",
            "Clears all male adult organs (`testis`, `epididymis`, `prostate`, `seminal vesicle`, non-secretory `XY resting_breast` < 8; < 16 with `IHC−`)",
            "`CSN2`, `LALBA`, `ZP4`, `DPPA3`, `CPSF4L`, `IZUMO1R`, `ACCSL`, `TBXT`, `KHDC3L`, `ZAR1L`, `VRTN`, `PRR23E`, `BMP15` (`RA-F-only`), `WEE2` (`RA-F-only`), `PABPN1L` (`RA-F-only`)",
        ],
        [
            "**Shared male & female**",
            f"`RA-MF` (`SRA`): **{int(sex_table.sex_restriction_class.eq('Shared male & female').sum()):,}**",
            f"`CTA-MF`: **{int(sex_cta.sex_restriction_class.eq('Shared male & female').sum()):,}**",
            f"`RA-MF-only`: **{int(sex_ra_only.sex_restriction_class.eq('Shared male & female').sum()):,}**",
            "≥ 8 TPM/nTPM (or ≥ 16 with `IHC−`, or positive IHC) in adult organs of both sexes",
            "`PRAME`, `SCGB1D4`, `CTAG2`, `MOS`, `ACTL8`, `NLRP4`, `OOSP2`, `SPO11`, `DMRTC2`, `PIWIL3`, `FOXR1`, `NLRP13` (`RA-MF-only`)",
        ],
    ]

    # 2. Complete table of all RA-F (CTA-F + RA-F-only) antigens
    raf_rows = []
    for r in sex_table[sex_table.sex_restriction_class.eq("Female-restricted")].itertuples():
        cat_r = data.loc[data.entity_id.eq(r.entity_id)].iloc[0]
        e_he = he.get(r.entity_id, {})
        rb = e_he.get("resting_breast", {})
        t_ntpm = e_he.get("testis", {}).get("max_ntpm")
        raf_rows.append(
            [
                f"`{r.protein_group_label}`",
                f"`{r.sex_restriction_code}`",
                number(cat_r.normal_reference_tpm),
                f"{number(cat_r.reproductive_supported_tpm)} · {r.reproductive_qualifying_compartments}",
                f"{number(r.male_sex_organ_max_tpm)} TPM ({number(t_ntpm)} nTPM)",
                number(rb.get("male_max_tpm")),
                f"{number(r.tumor_max_tpm)} · `{r.tumor_top_cancer_code}`"
                if pd.notna(r.tumor_max_tpm)
                else "—",
            ]
        )

    # 3. Highlighted gene-by-gene audit table across the 5 female biological programs
    audit_by_sym = audit_df.set_index("gene_symbol")
    highlight_rows = []
    for sym in FEMALE_HIGHLIGHT_SYMBOLS:
        if sym not in audit_by_sym.index:
            continue
        ar = audit_by_sym.loc[sym]
        call_tag = (
            f"`{ar['sex_restriction_code']}`"
            if ar["sex_restriction_code"]
            else f"`{ar['classification']}`"
        )
        pan_str = (
            f"{number(ar['pan_tumor_max_tpm'])} (`{ar['pan_tumor_top_cohort']}`)"
            if pd.notna(ar["pan_tumor_max_tpm"]) and ar["pan_tumor_top_cohort"]
            else "—"
        )
        highlight_rows.append(
            [
                ar["program"],
                f"`{sym}`",
                call_tag,
                ar["top_brca_gyn_cohorts"] or pan_str,
                pan_str,
                ar["disposition_mechanism"],
            ]
        )

    ra_m_cnt = int(sex_table.sex_restriction_class.eq("Male-restricted").sum())
    ra_f_cnt = int(sex_table.sex_restriction_class.eq("Female-restricted").sum())
    cta_m_cnt = int(sex_cta.sex_restriction_class.eq("Male-restricted").sum())
    cta_f_cnt = int(sex_cta.sex_restriction_class.eq("Female-restricted").sum())
    cta_p_cnt = int(sex_cta.sex_restriction_class.eq("Pregnancy / peripartum-only").sum())
    extra_vars = {
        "ra_m_n": f"{ra_m_cnt:,}",
        "ra_f_n": f"{ra_f_cnt:,}",
        "ra_e_n": f"{int(sex_table.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
        "ra_mf_n": f"{int(sex_table.sex_restriction_class.eq('Shared male & female').sum()):,}",
        "ra_p_n": f"{int(sex_table.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
        "cta_p_n": f"{cta_p_cnt:,}",
        "peri_ra_n": f"{peri_ra_n:,}",
        "peri_cta_n": f"{peri_cta_n:,}",
        "peri_dual_ra_n": f"{peri_dual_ra_n:,}",
        "peri_dual_cta_n": f"{peri_dual_cta_n:,}",
        "ra_m_only_n": f"{int(sex_ra_only.sex_restriction_class.eq('Male-restricted').sum()):,}",
        "ra_f_only_n": f"{int(sex_ra_only.sex_restriction_class.eq('Female-restricted').sum()):,}",
        "ra_e_only_n": f"{int(sex_ra_only.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
        "ra_p_only_n": f"{int(sex_ra_only.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
        "ra_mf_only_n": f"{int(sex_ra_only.sex_restriction_class.eq('Shared male & female').sum()):,}",
        "cta_female_safe_pct": f"{100.0 * (cta_m_cnt + cta_p_cnt) / max(len(sex_cta), 1):.1f}",
        "cta_male_safe_pct": f"{100.0 * (cta_f_cnt + cta_p_cnt) / max(len(sex_cta), 1):.1f}",
        "cta_m_f_ratio": f"{cta_m_cnt / max(cta_f_cnt, 1):.1f}",
        "ra_m_f_ratio": f"{ra_m_cnt / max(ra_f_cnt, 1):.1f}",
        "female_audit_n": f"{len(audit_df):,}",
        "sex_partition_table": table(
            [
                "Sex-restriction class",
                "All RAs (`RA-*`)",
                "CTAs (`CTA-*`)",
                "RA-only (`RA-*-only`)",
                "Adult host-protection criterion",
                "Representative proteins",
            ],
            partition_rows,
        ),
        "ra_f_all_table": table(
            [
                "Protein",
                "Subset",
                "Somatic N (TPM)",
                "Qualifying RA peak (TPM) · compartments",
                "Male repro max (testis scRNA)",
                "XY resting breast (TPM)",
                "Pan-cancer tumor peak",
            ],
            raf_rows,
        ),
        "female_audit_summary_table": table(
            [
                "Canonical female program",
                "Genes",
                "Stage 1 Unsafe",
                "Other / CA-only",
                "RA-M / CTA-M",
                "RA-MF / CTA-MF",
                "RA-P / CTA-P",
                "RA-F (`CTA-F` / `RA-F-only`)",
            ],
            prog_summary_rows,
        ),
        "female_audit_highlights_table": table(
            [
                "Biological program",
                "Gene",
                "Call",
                "Top BRCA / TNBC / Gyn cohorts (TPM)",
                "Pan-cancer peak (TPM)",
                "Quantitative disposition & gating measurement",
            ],
            highlight_rows,
        ),
    }
    return extra_vars, audit_df


def build():
    data = target_selection.results()
    meta = target_selection.metadata()["summary"]
    data["loci"] = data.member_gene_ids.map(json.loads).map(len)
    loci = data.member_gene_ids.map(json.loads).explode()
    if not loci.is_unique or not data.entity_id.is_unique:
        raise ValueError("Catalog is not a partition of source loci")
    ca, ra = data.is_cancer_antigen, data.is_reproductive_antigen
    if not ((ca & ra).equals(data.is_cta) and (ca | ra).equals(data.is_target_candidate)):
        raise ValueError("Catalog membership flags disagree")
    for source, prefix in [
        ("normal_reference_contexts", "normal"),
        ("reproductive_supported_observations", "reproductive"),
    ]:
        for field in ["source_id", "tissue", "source_context"]:
            data[f"{prefix}_{field}"] = data[source].map(lambda x, field=field: contexts(x, field))
    funnel_steps = [
        ("Normal-tissue safety exclusions (Unsafe)", int(data.category.eq("Unsafe").sum())),
        (
            "Coverage, identity, or rare-cell review (Other:Review)",
            int(data.classification.eq("Other:Review").sum()),
        ),
        (
            "Nested, readthrough, or Unsafe 9-mer sequence (Other:Sequence)",
            int(data.classification.eq("Other:Sequence").sum()),
        ),
        (
            "Unscreened adult home sub-niche (Other:Unscreened tissue)",
            int(data.classification.eq("Other:Unscreened tissue").sum()),
        ),
        (
            "Below 10x enrichment or Disposable cap (Other:Expression)",
            int(data.classification.eq("Other:Expression").sum()),
        ),
    ]
    rem = len(data)
    funnel_rows = []
    for step_label, removed_n in funnel_steps:
        rem -= removed_n
        funnel_rows.append({"step": step_label, "removed": removed_n, "remaining": rem})
    funnel = pd.DataFrame(funnel_rows)
    if not funnel.remaining.eq(len(data) - funnel.removed.cumsum()).all() or (
        funnel.remaining.iloc[-1] != data.is_target_candidate.sum()
    ):
        raise ValueError("Funnel does not reconcile with protein membership")

    sets = {"CA": ca, "RA": ra, "CTA": ca & ra, "CA only": ca & ~ra, "RA only": ra & ~ca}
    candidates = data[ca | ra].sort_values(PEAK, ascending=[False, True, True]).copy()
    for name, mask in sets.items():
        ranked = candidates.index[
            candidates.index.isin(data.index[mask]) & candidates.tumor_max_tpm.notna()
        ]
        candidates.loc[ranked, "rank_" + name.lower().replace(" ", "_")] = pd.array(
            range(1, len(ranked) + 1), dtype="Int64"
        )
    witnesses = pd.read_csv(
        BUILD / "tumor-witnesses.csv.gz",
        usecols=["entity_id", "source_id", "sample_id", "value", "qc_pass", "identity_resolved"],
    )
    in_range = witnesses.loc[
        witnesses.entity_id.isin(candidates.entity_id)
        & witnesses.qc_pass
        & witnesses.identity_resolved
        & witnesses.value.between(136, 1300, inclusive="neither")
    ]
    peak_sections = []
    for label, mask in [
        ("CTA", ca & ra),
        ("CA only", ca & ~ra),
        ("RA only", ra & ~ca),
        ("Other", data.category.eq("Other")),
    ]:
        group = data[mask]
        measured = group.dropna(subset=["tumor_max_tpm"]).sort_values(
            PEAK, ascending=[False, True, True]
        )
        if len(measured) == len(group):
            note = f"Top 10 of {len(group):,} groups by peak tumor sample TPM."
        else:
            note = (
                f"Top 10 of {len(measured):,} groups with complete baseline TPM "
                f"(out of {len(group):,} total `{label}` groups)."
            )
        peak_sections.append(f"### {label}\n\n{note}")
        if measured.empty:
            continue
        peak_sections.append(
            table(
                [
                    "Rank",
                    "Protein",
                    "Peak cohort",
                    "Tumor TPM",
                    "N TPM",
                    "Reproductive peak TPM",
                    "Tumor fold",
                    "Known patients at 10×",
                ],
                [
                    [
                        i,
                        r.protein_group_label,
                        r.tumor_top_cancer_code,
                        number(r.tumor_max_tpm),
                        number(r.normal_reference_tpm),
                        number(r.reproductive_supported_tpm),
                        number(r.tumor_fold),
                        number(r.independent_patients_at_tenfold),
                    ]
                    for i, r in enumerate(measured.head(10).itertuples(), 1)
                ],
            )
        )
    study = study_table(data)
    inventory = antigen_data.sources().copy()
    inventory["contents"] = inventory.dataset.map(lambda name: catalog.dataset(name).description)

    # Expose the actual schema for shipped tables, including evidence fields
    # that a source title alone cannot describe.
    def columns(source):
        path = ROOT / source
        if source.startswith("oncoref/data/") and source.endswith((".csv", ".csv.gz")):
            return " | ".join(pd.read_csv(path, nrows=0).columns)
        return ""

    inventory["table_columns"] = inventory.source_file.map(columns)
    coverage = pd.read_csv(BUILD / "tumor-source-coverage.csv")
    fields = [
        "entity_id",
        "protein_group_label",
        "member_gene_ids",
        "classification",
        "normal_reference_tpm",
        "normal_source_id",
        "normal_tissue",
        "normal_source_context",
        "disposable_max_tpm",
        "single_tumor_threshold_tpm",
        "tumor_max_tpm",
        "tumor_top_cancer_code",
        "tumor_top_sample",
        "tumor_fold",
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "anonymous_samples_at_tenfold",
        "cancer_types_at_tenfold",
        "reproductive_supported_tpm",
        "reproductive_source_id",
        "reproductive_tissue",
        "reproductive_source_context",
        "reproductive_fold",
    ]
    studied_fields = [
        "familiar_names",
        "papers_cta",
        "papers_counted",
        "nomination_source_ids",
        "clinical_sources",
        "clinical_scope",
        "disposition",
        "rare_cell_status",
        "finding_tissue",
        "finding_top_cell_type",
        "healthy_maximum_ntpm",
        "matched_normal_assay",
    ]
    outputs = {
        "tissue-hierarchy.md": tissue_hierarchy(),
        f"{ASSETS}/sources.csv": inventory.to_csv(index=False, lineterminator="\n"),
        f"{ASSETS}/sex-restricted-antigens.csv": human.sex_restricted_antigens().to_csv(
            index=False, lineterminator="\n"
        ),
        f"{ASSETS}/tumor-ranks.csv": candidates[
            fields + [c for c in candidates if c.startswith("rank_")]
        ].to_csv(index=False, lineterminator="\n"),
        f"{ASSETS}/studied-proteins.csv": study[fields + studied_fields].to_csv(
            index=False, lineterminator="\n"
        ),
    }
    examples = study[study.familiar_names.notna() | study.clinical_scope.notna()]
    held = examples[examples.category.isin(["Unsafe", "Other"])].copy()
    other_examples = study[
        study.protein_group_label.isin(["CAGE1", "DAZL", "POTEC"])
        & study.category.isin(["Other", "Unsafe"])
    ].copy()
    highlighted_holds = pd.concat([held, other_examples]).drop_duplicates("entity_id").copy()

    gate_labels = {
        "Normal bulk": "Bulk tissue RNA",
        "HPA tissue clusters": "Single-cell clusters (HPA)",
        "HPA normal protein": "Normal IHC",
        "Accepted normal-protein evidence": "Curated protein staining",
        "HPA immune": "Sorted immune RNA",
        "Monaco immune": "Sorted immune RNA",
        "DICE immune": "Sorted immune RNA",
        "Adult cornea": "Adult cornea scRNA",
        "Established rare normal-cell RNA": "Donor rare-cell RNA",
        "Reviewed normal-cell knowledge": "Curated normal literature",
        "Published normal-tissue expression": "Curated normal literature",
    }
    rare_basis_short = {
        "detection confirmed by a published gene-specific assay": "qPCR-confirmed",
        "replicated signal at a locus without mapping risk": "replicated, unique mapping",
        "replicated normal-cell signal with unresolved locus attribution": "replicated, paralog-ambiguous",
        "frequent single-study normal-cell signal with unresolved locus attribution": "single-study, paralog-ambiguous",
    }
    context_labels = {
        "B_CELL_NAIVE": "naive B cells",
        "CD4_STIM": "stimulated CD4 T cells",
        "CD8_STIM": "stimulated CD8 T cells",
        "NK_CELL": "NK cells",
        "MONOCYTES": "monocytes",
        "Bladder": "urinary bladder",
    }

    def clean_context_label(raw):
        if pd.isna(raw) or not raw:
            return "—"
        parts = []
        for token in str(raw).split(" | "):
            t = token.strip()
            if t in context_labels:
                parts.append(context_labels[t])
                continue
            base, dot, suffix = t.rpartition(".")
            if dot and suffix.isdigit():
                t = base
            if t.startswith("Brain_"):
                t = t.removeprefix("Brain_")
            if " day / unspecified lineage" in t:
                t = t.replace(" day / unspecified lineage", "-day embryo")
            elif " day / " in t:
                t = t.replace(" day / ", "-day ")
            parts.append(t.replace("_", " ").lower())
        return " | ".join(dict.fromkeys(parts))

    def concise_gate(r):
        if r.category == "Unsafe":
            raw_parts = [
                re.sub(r"^#\d+\s+", "", p.strip())
                for p in str(r.Unsafe_reasons).split(" | ")
                if p.strip()
            ]
            mapped = list(
                dict.fromkeys(
                    gate_labels.get(p, p[:1].upper() + p[1:] if p else p) for p in raw_parts
                )
            )
            if raw_parts == ["Established rare normal-cell RNA"] and pd.notna(r.rare_cell_basis):
                note = rare_basis_short.get(r.rare_cell_basis, r.rare_cell_basis)
                return f"Donor rare-cell RNA ({note})"
            return " + ".join(mapped)
        if r.subtype == "Review" and r.rare_cell_status == "unresolved":
            note = rare_basis_short.get(r.rare_cell_basis, r.rare_cell_basis)
            return f"Donor rare-cell RNA ({note})"
        return r.disposition

    def strongest_signal(r):
        if pd.notna(r.matched_normal_assay) and r.matched_normal_assay:
            return r.matched_normal_assay
        if pd.notna(r.finding_tissue) and pd.notna(r.finding_top_multi_molecule_pct):
            return (
                f"{r.finding_tissue}, {r.finding_top_cell_type} · "
                f"{r.finding_top_multi_molecule_pct:.2f}% of cells (≥2 molecules)"
            )
        return "—"

    highlighted_holds["concise_gate"] = [concise_gate(r) for r in highlighted_holds.itertuples()]
    highlighted_holds["strongest_signal"] = [
        strongest_signal(r) for r in highlighted_holds.itertuples()
    ]
    most = study.head(10)
    # Include familiar reference points even when they have fewer curated CTA papers.
    anchors = study[study.protein_group_label.isin(["NY-ESO-1", "PRAME", "CT83", "CTAG2", "CAGE1"])]
    highlighted = pd.concat([most, anchors]).drop_duplicates("entity_id")
    comparison_fields = ["normal_reference_tpm", "reproductive_supported_tpm", "tumor_max_tpm"]
    skipped = (
        highlighted.category.eq("Unsafe")
        & highlighted.cancer_evidence.eq("unassessed")
        & highlighted[comparison_fields].isna().all(axis=1)
    )
    skipped_ids = set(highlighted.loc[skipped, "entity_id"])

    numeric_cols = {"normal": "normal_reference_tpm", "reproductive": "reproductive_supported_tpm"}

    def measure(r, prefix):
        if r.entity_id in skipped_ids:
            return "NR"
        v = getattr(r, numeric_cols[prefix])
        where = clean_context_label(getattr(r, prefix + "_source_context"))
        return f"{number(v)} · {where}" if pd.notna(v) else "—"

    partition = data.groupby("classification").loci.agg(["size", "sum"])
    context_path = BUILD / "normal-rna/source-tissues.csv"
    source_coverage = pd.read_csv(context_path)
    lit_sources = pd.read_csv(BUILD / "summary/literature-sources.csv")
    published = data.published_nomination.fillna(False).astype(bool)
    shared_ok_mask = ~data.category.eq("Unsafe") & data[
        list(target_selection.SHARED_CHECKS)
    ].astype(bool).all(axis=1)
    other_compared_n = int(
        data[data.category.eq("Other") & ~shared_ok_mask]["tumor_max_tpm"].notna().sum()
    )
    profiles = target_filters.profile_summary(data).set_index("profile")
    testis_profile = profiles.loc["male-patient (protect testis)"]
    pregnant_profile = profiles.loc["female-patient (pregnant)"]
    sex_table = human.sex_restricted_antigens()
    sex_cta = sex_table[sex_table["is_cta"].astype(bool)]
    unsafe_reasons = data.Unsafe_reasons.fillna("").map(
        lambda v: set(v.split(" | ")) if v else set()
    )
    bulk_unsafe = unsafe_reasons.map(lambda s: "Normal bulk" in s)
    any_unsafe = data.category.eq("Unsafe")
    values = {
        "policy": target_selection.POLICY_ID,
        "other_compared": f"{other_compared_n:,}",
        "n": f"{len(data):,}",
        "loci": f"{len(loci):,}",
        "missing_seq": f"{data.representative_sequence_sha256s.map(json.loads).map(len).eq(0).sum():,}",
        "ca": f"{int(ca.sum()):,}",
        "ca_only": f"{int((ca & ~ra).sum()):,}",
        "ra": f"{int(ra.sum()):,}",
        "ra_only": f"{int((ra & ~ca).sum()):,}",
        "cta": f"{int((ca & ra).sum()):,}",
        "cta_published": f"{int((ca & ra & published).sum()):,}",
        "cta_new": f"{int((ca & ra & ~published).sum()):,}",
        "published_noms": f"{int(published.sum()):,}",
        "published_surviving": f"{int(((ca | ra) & published).sum()):,}",
        "published_outside": f"{int((~ca & ~ra & published).sum()):,}",
        "union": f"{int((ca | ra).sum()):,}",
        "unsafe_n": f"{int(data.category.eq('Unsafe').sum()):,}",
        "bulk_unsafe_n": f"{int(bulk_unsafe.sum()):,}",
        "non_bulk_added_n": f"{int((any_unsafe & ~bulk_unsafe).sum()):,}",
        "other_n": f"{int(data.category.eq('Other').sum()):,}",
        "testis_protected_cta": f"{int(testis_profile.CTA):,}",
        "pregnant_retained": f"{int(pregnant_profile.CTA):,}",
        "pregnant_unresolved": f"{int(pregnant_profile.CTA_unresolved):,}",
        "pregnant_excluded": f"{int(pregnant_profile.CTA_excluded):,}",
        "mra_n": f"{int(sex_table.sex_restriction_class.eq('Male-restricted').sum()):,}",
        "fra_n": f"{int(sex_table.sex_restriction_class.eq('Female-restricted').sum()):,}",
        "era_n": f"{int(sex_table.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
        "sra_n": f"{int(sex_table.sex_restriction_class.eq('Shared male & female').sum()):,}",
        "cta_m_n": f"{int(sex_cta.sex_restriction_class.eq('Male-restricted').sum()):,}",
        "cta_f_n": f"{int(sex_cta.sex_restriction_class.eq('Female-restricted').sum()):,}",
        "cta_e_n": f"{int(sex_cta.sex_restriction_class.eq('Pregnancy / peripartum-only').sum()):,}",
        "cta_mf_n": f"{int(sex_cta.sex_restriction_class.eq('Shared male & female').sum()):,}",
        "cta_female_safe_n": f"{int(sex_cta.sex_restriction_code.isin(['CTA-M', 'CTA-P']).sum()):,}",
        "cta_male_safe_n": f"{int(sex_cta.sex_restriction_code.isin(['CTA-F', 'CTA-P']).sum()):,}",
        "cta_m_sex_only_n": f"{int((sex_cta.female_sex_organ_only_eligible & ~sex_cta.male_all_adult_eligible).sum()):,}",
        "fold": f"{target_selection.ENRICHMENT_FOLD:g}",
        "floor": f"{target_selection.DETECTION_FLOOR:g}",
        "cap": f"{target_selection.DISPOSABLE_CAP:g}",
        "cells": target_selection.RARE_CELL_MIN_CELLS,
        "pct": f"{target_selection.RARE_CELL_DETECTION_PCT:g}",
        "frequent": f"{target_selection.RARE_CELL_SINGLE_STUDY_PCT:g}",
        "rna_bases": ", ".join(
            f"{v:g} for {k.title()}" for k, v in cta_selection_policy.TIER_BASES.items()
        ),
        "manageable_ceiling": f"{cta_selection_policy.MANAGEABLE_UNSAFE_FLOOR:g}",
        "shared": f"{meta['shared_eligible']:,}",
        "assessed": f"{meta['assessed_for_cancer_evidence']:,}",
        "required": len(meta["comparison_required_contexts"]),
        "files": f"{len(inventory):,}",
        "cohorts": f"{int(coverage.status.eq('assessed').sum()):,}",
        "tissue_rows": f"{len(source_coverage):,}",
        "lit_sources": f"{len(lit_sources):,}",
        "studied": f"{len(study):,}",
        "study_extra": f"{len(study) - int(published.sum()):,}",
        "held": f"{len(held):,}",
        "anchors": f"{len(examples):,}",
        "nr_proteins": ", ".join(highlighted.loc[skipped, "protein_group_label"]),
        "sources": table(
            ["Inventory role", "Versioned files"],
            [[role, f"{count:,}"] for role, count in sorted(inventory.role.value_counts().items())],
        ),
        "partition": table(
            ["Exclusive category", "Protein records", "Member loci"],
            [[k, f"{r['size']:,}", f"{r['sum']:,}"] for k, r in partition.iterrows()],
        ),
        "funnel": table(
            ["First attributed exclusion", "Removed", "Remaining"],
            [[r.step, f"{r.removed:,}", f"{r.remaining:,}"] for r in funnel.itertuples()],
        ),
        "peaks": "\n\n".join(peak_sections),
        "union_ge_100": f"{int(candidates.tumor_max_tpm.ge(100).sum()):,}",
        "union_ge_1000": f"{int(candidates.tumor_max_tpm.ge(1000).sum()):,}",
        "union_patients_ge_10": f"{int(candidates.independent_patients_at_tenfold.ge(10).sum()):,}",
        "cta_median_tpm": number(data.loc[ca & ra, "tumor_max_tpm"].median()),
        "ca_only_median_tpm": number(data.loc[ca & ~ra, "tumor_max_tpm"].median()),
        "range_observations": f"{len(in_range):,}",
        "range_samples": f"{len(in_range[['source_id', 'sample_id']].drop_duplicates()):,}",
        "midrange": int(candidates.tumor_max_tpm.between(136, 1300, inclusive="neither").sum()),
        "publications": table(
            [
                "Protein",
                "CTA papers",
                "Catalog",
                "N TPM · site",
                "RA peak TPM · context",
                "Tumor TPM · cohort",
            ],
            [
                [
                    r.protein_group_label,
                    r.papers_cta,
                    r.classification,
                    measure(r, "normal"),
                    measure(r, "reproductive"),
                    f"{number(r.tumor_max_tpm)} · {r.tumor_top_cancer_code}"
                    if pd.notna(r.tumor_max_tpm)
                    else ("NR" if r.entity_id in skipped_ids else "—"),
                ]
                for r in highlighted.itertuples()
            ],
        ),
        "excluded": table(
            [
                "Studied protein",
                "Category",
                "Excluding or review gate",
                "Strongest measured normal signal",
            ],
            highlighted_holds.loc[
                highlighted_holds.protein_group_label.isin(
                    [
                        "MAGEA3",
                        "MAGEA12",
                        "CTCFL",
                        "CTAG2",
                        "PLAC1",
                        "SPA17",
                        "DDX43",
                        "CLDN6",
                        "CAGE1",
                        "DAZL",
                        "POTEC",
                    ]
                ),
                ["protein_group_label", "classification", "concise_gate", "strongest_signal"],
            ].itertuples(index=False, name=None),
        ),
    }
    h_cta_df = data[ca & ra]
    h_ca_only_df = data[ca & ~ra]
    h_ect = int(h_cta_df["tumor_activation_mode"].eq("ectopic").sum())
    h_both = int(h_cta_df["tumor_activation_mode"].eq("ectopic_and_lineage_retained").sum())
    h_lin = int(h_cta_df["tumor_activation_mode"].eq("lineage_retained").sum())
    h_any_ect = h_ect + h_both
    h_ca_ect = int(
        h_ca_only_df["tumor_activation_mode"]
        .isin(["ectopic", "ectopic_and_lineage_retained"])
        .sum()
    )
    h_ca_lin = int(h_ca_only_df["tumor_activation_mode"].eq("lineage_retained").sum())
    h_cta_lin_rows = h_cta_df[h_cta_df["tumor_activation_mode"].eq("lineage_retained")].sort_values(
        ["tumor_max_tpm", "protein_group_label"], ascending=[False, True]
    )
    h_gpcr_rx = re.compile(r"^(?:OR\d+[A-Z]+\d+P?|TAAR\d+P?|TAS1R\d+|TAS2R\d+|VN1R\d+)$")
    d_gpcr_rx = re.compile(
        r"^(?:COR\d+|OR\d+[A-Z]*\d*|TAAR\d+|TAS1R\d+|TAS2R\d+|VN1R\d+)$", re.IGNORECASE
    )
    m_gpcr_rx = re.compile(
        r"^(?:Olfr\d+|Or\d+[a-z]+\d+|Taar\d+[a-z]*|Tas1r\d+|Tas2r\d+|Vmn1r\d+|Vmn2r\d+)$",
        re.IGNORECASE,
    )
    d_cat_full = canine.target_selection()
    m_cat_full = mouse.target_selection()
    values.update(
        {
            "total_cohorts": "135",
            "baseline_cohorts": "132",
            "cta_ectopic_n": f"{h_ect:,}",
            "cta_both_mode_n": f"{h_both:,}",
            "cta_ectopic_any_n": f"{h_any_ect:,}",
            "cta_ectopic_any_pct": f"{100.0 * h_any_ect / max(len(h_cta_df), 1):.1f}",
            "cta_lineage_only_n": f"{h_lin:,}",
            "cta_lineage_only_genes": ", ".join(
                f"`{lbl}`" for lbl in h_cta_lin_rows["protein_group_label"].astype(str)
            ),
            "ca_only_ectopic_n": f"{h_ca_ect:,}",
            "ca_only_lineage_n": f"{h_ca_lin:,}",
            "unscreened_n": f"{int(data.classification.eq('Other:Unscreened tissue').sum()):,}",
            "human_gpcr_n": f"{int(data['protein_group_label'].astype(str).map(lambda x: bool(h_gpcr_rx.match(x.split('/')[0]))).sum()):,}",
            "dog_gpcr_n": f"{int(d_cat_full['protein_group_label'].astype(str).map(lambda x: bool(d_gpcr_rx.match(x.split('/')[0]))).sum()):,}",
            "mouse_gpcr_n": f"{int(m_cat_full['protein_group_label'].astype(str).map(lambda x: bool(m_gpcr_rx.match(x.split('/')[0]))).sum()):,}",
        }
    )
    values.update(ctag2_evidence(data))
    sex_interrog_vars, female_audit_df = build_sex_restriction_interrogation(data, sex_table)
    values.update(sex_interrog_vars)
    outputs[f"{ASSETS}/female-lineage-audit.csv"] = female_audit_df.to_csv(
        index=False, lineterminator="\n"
    )
    ra_all_rows = data[data.is_reproductive_antigen.astype(bool)].copy()
    ra_all_rows["has_ep"] = (
        ra_all_rows["reproductive_qualifying_compartments"]
        .fillna("")
        .map(
            lambda s: bool(
                {"embryo", "placenta", "lactating_breast"} & {c.strip() for c in str(s).split("|")}
            )
        )
    )
    ep_rows = ra_all_rows[ra_all_rows["has_ep"]]
    ep_cta_mask = ep_rows["is_cancer_antigen"].astype(bool)
    values.update(
        {
            "ep_ra_n": f"{len(ep_rows):,}",
            "ep_cta_n": f"{int(ep_cta_mask.sum()):,}",
            "ep_m_ra_n": f"{int(ep_rows['sex_restriction_class'].eq('Male-restricted').sum()):,}",
            "ep_m_cta_n": f"{int((ep_rows['sex_restriction_class'].eq('Male-restricted') & ep_cta_mask).sum()):,}",
            "ep_mf_ra_n": f"{int(ep_rows['sex_restriction_class'].eq('Shared male & female').sum()):,}",
            "ep_mf_cta_n": f"{int((ep_rows['sex_restriction_class'].eq('Shared male & female') & ep_cta_mask).sum()):,}",
            "ep_f_ra_n": f"{int(ep_rows['sex_restriction_class'].eq('Female-restricted').sum()):,}",
            "ep_f_cta_n": f"{int((ep_rows['sex_restriction_class'].eq('Female-restricted') & ep_cta_mask).sum()):,}",
        }
    )
    tcr_eval_df = tcr_trials.evaluate_against_gold_standard()
    tcr_summary = tcr_trials.gold_standard_summary(data)
    stratum_labels = {
        "safe_cta_potency_proven": "Potent, safe CTA",
        "toxic_on_target": "On-target organ toxicity",
        "toxic_cross_reactivity": "Fatal cross-reactivity",
        "conservative_gene_level_unsafe": "Tolerated non-CTA epitope",
        "excluded_not_potency_proven": "Excluded (0% ORR / pending)",
    }
    values.update(
        {
            "tcr_gold_n": f"{tcr_summary['total_entries']:,}",
            "tcr_potent_n": f"{tcr_summary['potency_proven_entries']:,}",
            "tcr_nonpotent_n": f"{tcr_summary['excluded_non_potent_or_pending_entries']:,}",
            "tcr_toxic_n": f"{tcr_summary['severe_specificity_toxic_entries']:,}",
            "tcr_toxic_caught_n": f"{tcr_summary['severe_specificity_toxic_flagged_unsafe']:,}",
            "tcr_safe_cta_n": f"{tcr_summary['potency_proven_safe_cta_entries']:,}",
            "tcr_safe_cta_kept_n": f"{tcr_summary['potency_proven_safe_cta_retained']:,}",
            "tcr_conservative_n": f"{tcr_summary['conservative_non_cta_unsafe_entries']:,}",
            "tcr_gold_standard_table": table(
                [
                    "Target (cross-hit)",
                    "Product · HLA · peptide",
                    "Potency (n, response)",
                    "Observed human clinical toxicity",
                    "OncoRef call & strongest proxy signal",
                    "Stratum",
                ],
                [
                    [
                        f"**{r.target_symbol}**"
                        + (
                            f" (`→ {r.cross_recognized_symbols}`)"
                            if pd.notna(r.cross_recognized_symbols) and r.cross_recognized_symbols
                            else ""
                        ),
                        f"{r.product_name} · `{r.hla_restriction}` · `{r.peptide}` "
                        f"([{r.citation.split(';')[0]}]({r.source_url}))",
                        (
                            f"n={int(r.patients_treated)}, **100%** RFS (44 mo)"
                            if r.benchmark_id == "wt1_ttcr_c4"
                            else (
                                f"n={int(r.patients_treated)}, **{r.orr_pct:g}%** ORR"
                                if pd.notna(r.patients_treated)
                                and pd.notna(r.orr_pct)
                                and r.orr_pct > 0
                                else (
                                    f"n={int(r.patients_treated)}, lethal d4–5 lysis"
                                    if r.benchmark_id == "magea3_a1_titin"
                                    else (
                                        f"n={int(r.patients_treated)}, **0%** ORR"
                                        if pd.notna(r.patients_treated)
                                        else "Phase 1 pending"
                                    )
                                )
                            )
                        ),
                        r.observed_clinical_toxicities,
                        f"`{r.oncoref_classification}` · {r.strongest_proxy_context} "
                        f"(`{r.strongest_proxy_ratio:g}×`)",
                        stratum_labels[r.clinical_stratum],
                    ]
                    for r in tcr_eval_df.itertuples(index=False)
                ],
            ),
        }
    )
    outputs[f"{ASSETS}/clinical-tcr-gold-standard.csv"] = tcr_eval_df.to_csv(
        index=False, lineterminator="\n"
    )
    species_vars, species_csvs = build_species_sections()
    values.update(species_vars)
    outputs.update(species_csvs)
    outputs["antigen-whitepaper.md"] = Template(TEMPLATE.read_text()).substitute(values)
    outputs["species-antigen-whitepaper.md"] = Template(SPECIES_TEMPLATE.read_text()).substitute(
        values
    )
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Compare to the committed artifacts")
    args = parser.parse_args()
    for name, content in build().items():
        path = DOCS / name
        if args.check:
            if not path.exists() or path.read_bytes() != content.encode("utf-8"):
                parser.exit(1, f"Stale whitepaper artifact: {name}\n")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content.encode("utf-8"))


if __name__ == "__main__":
    main()
