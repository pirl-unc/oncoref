#!/usr/bin/env python3
"""Build schema-validated clinical organ-clearances ledger & audit report (#766).

Verifies exact Ensembl 115 coding-transcript peptide containment across all reviewed
TCR-T products in ``oncoref.tcr_trials``, records negative-clearance toxicity anchors
(``magea3_a3a``, ``magea3_a2_neuro``, ``tbi1301``), and records the Ensembl 115 Xq28
locus/paralog attribution for ``CTAG2`` (`#711`/`#766`).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pyensembl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oncoref import clinical_organ_clearance as coc  # noqa: E402
from oncoref import tcr_trials  # noqa: E402

DATA_DIR = ROOT / "oncoref/data"
AUDIT_DIR = ROOT / "docs/audits/clinical-organ-clearance-20261123"


def _partition_isoforms(
    ens: pyensembl.EnsemblRelease, gene_ids: list[str], peptide: str
) -> tuple[list[str], list[str]]:
    cleared, uncleared = [], []
    for gid in sorted(gene_ids):
        g = ens.gene_by_id(gid)
        for tx in sorted(g.transcripts, key=lambda t: t.transcript_id):
            if not (tx.is_protein_coding and tx.contains_start_codon and tx.contains_stop_codon):
                continue
            if peptide in tx.protein_sequence:
                cleared.append(tx.transcript_id)
            else:
                uncleared.append(tx.transcript_id)
    return cleared, uncleared


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    ens = pyensembl.EnsemblRelease(115)
    tcr = tcr_trials.evidence()
    sources = tcr["sources"]
    prods = {p["product_id"]: p for p in tcr["products"]}
    ts = pd.read_csv(DATA_DIR / "target-selection.csv.gz", low_memory=False)

    def _eid_and_label(gid: str) -> tuple[str, str, list[str]]:
        sub = ts[ts.member_gene_ids.map(lambda s: gid in json.loads(s))]
        r = sub.iloc[0]
        return str(r.entity_id), str(r.protein_group_label), json.loads(r.member_gene_ids)

    entries: list[dict] = []

    # 1. CTAG2 (LAGE-1a, ENSG00000126890) via lete_cel_c259
    ctag2_eid, ctag2_lbl, ctag2_gids = _eid_and_label("ENSG00000126890")
    ctag2_cl, ctag2_un = _partition_isoforms(ens, ctag2_gids, "SLLMWITQC")
    entries.append(
        {
            "clearance_id": "COC-20261123-CTAG2-LETE-CEL",
            "entity_id": ctag2_eid,
            "protein_group_label": ctag2_lbl,
            "gene_id": "ENSG00000126890",
            "product_id": "lete_cel_c259",
            "modality": "TCR-T",
            "status": "cleared",
            "clearance_scope": "epitope",
            "exposure_route": "systemic_iv",
            "barrier_compartment_exposure_verified": False,
            "objective_tumor_responses": True,
            "cohort_ids": ["lete_mrcls_20"],
            "hla_matched_patients_treated": 20,
            "cleared_organs": ["lung", "urinary bladder"],
            "cleared_epitopes": ["SLLMWITQC"],
            "hla_alleles": prods["lete_cel_c259"]["hla"],
            "cleared_transcripts": ctag2_cl,
            "uncleared_coding_transcripts": ctag2_un,
            "source_keys": ["lete_mrcls_2025", "lete_protocol"],
            "trial_source_ids": [
                sources["lete_mrcls_2025"]["source_id"],
                sources["lete_protocol"]["source_id"],
            ],
            "rationale": (
                "Intravenously infused letetresgene autoleucel (c259 TCR-T) recognizes HLA-A*02 / SLLMWITQC "
                "at positions 157-165 of LAGE-1a (ENST00000369585, 180 aa), achieves ~40% objective tumor "
                "responses across two sequential MRCLS cohorts (n=20 infused, JCO-24-01466 / NCT02992743), "
                "traverses the pulmonary capillary/alveolar bed on first pass, and shows zero attributed "
                "alveolar pneumonitis or cystitis. Explicitly excludes LAGE-1b (ENST00000247306, 210 aa, "
                "which lacks SLLMWITQC), heart (confounded day-161 cardiac arrest retained per #544), and "
                "brain (blood-brain barrier; brain snRNA-seq hold resolved separately by #711 Xq28 attribution)."
            ),
        }
    )

    # 2. NY-ESO-1 (CTAG1A/CTAG1B) via lete_cel_c259
    ny_eid, ny_lbl, ny_gids = _eid_and_label("ENSG00000184033")
    ny_cl, ny_un = _partition_isoforms(ens, ny_gids, "SLLMWITQC")
    entries.append(
        {
            "clearance_id": "COC-20261123-NYESO1-LETE-CEL",
            "entity_id": ny_eid,
            "protein_group_label": ny_lbl,
            "gene_id": "ENSG00000184033",
            "product_id": "lete_cel_c259",
            "modality": "TCR-T",
            "status": "cleared",
            "clearance_scope": "epitope",
            "exposure_route": "systemic_iv",
            "barrier_compartment_exposure_verified": False,
            "objective_tumor_responses": True,
            "cohort_ids": ["lete_mrcls_20"],
            "hla_matched_patients_treated": 20,
            "cleared_organs": ["lung", "urinary bladder"],
            "cleared_epitopes": ["SLLMWITQC"],
            "hla_alleles": prods["lete_cel_c259"]["hla"],
            "cleared_transcripts": ny_cl,
            "uncleared_coding_transcripts": ny_un,
            "source_keys": ["lete_mrcls_2025", "lete_protocol"],
            "trial_source_ids": [
                sources["lete_mrcls_2025"]["source_id"],
                sources["lete_protocol"]["source_id"],
            ],
            "rationale": (
                "Shared HLA-A*02 / SLLMWITQC target of letetresgene autoleucel (n=20 infused MRCLS patients); "
                "clears full-length 180 aa NY-ESO-1 transcripts while excluding truncated 168 aa ENST00000359887."
            ),
        }
    )

    # 3. PRAME via ima203 (n=40 infused, kidney + lung vascular perfusion)
    pr_eid, pr_lbl, pr_gids = _eid_and_label("ENSG00000185686")
    pr_cl, pr_un = _partition_isoforms(ens, pr_gids, "SLLQHLIGL")
    entries.append(
        {
            "clearance_id": "COC-20261123-PRAME-IMA203",
            "entity_id": pr_eid,
            "protein_group_label": pr_lbl,
            "gene_id": "ENSG00000185686",
            "product_id": "ima203",
            "modality": "TCR-T",
            "status": "cleared",
            "clearance_scope": "epitope",
            "exposure_route": "systemic_iv",
            "barrier_compartment_exposure_verified": False,
            "objective_tumor_responses": True,
            "cohort_ids": ["ima203_phase1_40"],
            "hla_matched_patients_treated": 40,
            "cleared_organs": ["kidney", "lung"],
            "cleared_epitopes": ["SLLQHLIGL"],
            "hla_alleles": prods["ima203"]["hla"],
            "cleared_transcripts": pr_cl,
            "uncleared_coding_transcripts": pr_un,
            "source_keys": ["ima203_2025"],
            "trial_source_ids": [sources["ima203_2025"]["source_id"]],
            "rationale": (
                "IMA203 autologous TCR-T targeting HLA-A*02:01 / SLLQHLIGL (n=40 infused solid-tumor patients, "
                "Nat Med 2025) achieved confirmed tumor responses with zero clinical renal or pulmonary on-target "
                "toxicity; liver excluded from cleared_organs due to CRS-associated grade 4 transaminitis."
            ),
        }
    )

    # 4. MAGEA4 via afami_cel_c1032 (n=90 infused across phase 1 + SPEARHEAD-1)
    m4_eid, m4_lbl, m4_gids = _eid_and_label("ENSG00000147381")
    m4_cl, m4_un = _partition_isoforms(ens, m4_gids, "GVYDGREHTV")
    entries.append(
        {
            "clearance_id": "COC-20261123-MAGEA4-AFAMI-CEL",
            "entity_id": m4_eid,
            "protein_group_label": m4_lbl,
            "gene_id": "ENSG00000147381",
            "product_id": "afami_cel_c1032",
            "modality": "TCR-T",
            "status": "cleared",
            "clearance_scope": "epitope",
            "exposure_route": "systemic_iv",
            "barrier_compartment_exposure_verified": False,
            "objective_tumor_responses": True,
            "cohort_ids": ["afami_phase1_38", "spearhead1_cohort1_52"],
            "hla_matched_patients_treated": 90,
            "cleared_organs": ["lung", "liver", "kidney"],
            "cleared_epitopes": ["GVYDGREHTV"],
            "hla_alleles": prods["afami_cel_c1032"]["hla"],
            "cleared_transcripts": m4_cl,
            "uncleared_coding_transcripts": m4_un,
            "source_keys": ["afami_phase1_2023", "spearhead1_2024", "afami_protocol"],
            "trial_source_ids": [
                sources["afami_phase1_2023"]["source_id"],
                sources["spearhead1_2024"]["source_id"],
                sources["afami_protocol"]["source_id"],
            ],
            "rationale": (
                "FDA-approved afamitresgene autoleucel (c1032 TCR-T, n=90 infused across NCT03132922 and "
                "SPEARHEAD-1 NCT04044768) clears GVYDGREHTV-containing isoforms in perfused visceral organs; "
                "bone marrow and brain excluded due to recorded conditioning aplasia and stroke fatalities."
            ),
        }
    )

    # 5-7. Negative-clearance clinical toxicity & unresolved-event anchors
    m3_eid, m3_lbl, _ = _eid_and_label("ENSG00000221867")
    m12_eid, m12_lbl, _ = _eid_and_label("ENSG00000213401")
    entries.extend(
        [
            {
                "clearance_id": "COC-20261123-NEG-MAGEA3-TITIN",
                "entity_id": m3_eid,
                "protein_group_label": m3_lbl,
                "gene_id": "ENSG00000221867",
                "product_id": "magea3_a3a",
                "modality": "TCR-T",
                "status": "blocked_clinical_toxicity",
                "clearance_scope": "",
                "hla_matched_patients_treated": 2,
                "cleared_organs": [],
                "cleared_epitopes": [],
                "hla_alleles": prods["magea3_a3a"]["hla"],
                "cleared_transcripts": [],
                "uncleared_coding_transcripts": [],
                "source_keys": ["magea3_titin_2013"],
                "trial_source_ids": [sources["magea3_titin_2013"]["source_id"]],
                "rationale": (
                    "Negative-clearance control: fatal cardiac toxicity in 2/2 treated HLA-A*01:01 patients "
                    "via titin (TTN ESDPIVAQY) cross-recognition (Linette et al. 2013)."
                ),
            },
            {
                "clearance_id": "COC-20261123-NEG-MAGEA12-NEURO",
                "entity_id": m12_eid,
                "protein_group_label": m12_lbl,
                "gene_id": "ENSG00000213401",
                "product_id": "magea3_a2_neuro",
                "modality": "TCR-T",
                "status": "blocked_clinical_toxicity",
                "clearance_scope": "",
                "hla_matched_patients_treated": 9,
                "cleared_organs": [],
                "cleared_epitopes": [],
                "hla_alleles": prods["magea3_a2_neuro"]["hla"],
                "cleared_transcripts": [],
                "uncleared_coding_transcripts": [],
                "source_keys": ["magea3_neuro_2013", "magea3_reexamination_2021"],
                "trial_source_ids": [
                    sources["magea3_neuro_2013"]["source_id"],
                    sources["magea3_reexamination_2021"]["source_id"],
                ],
                "rationale": (
                    "Negative-clearance control: fatal necrotizing leukoencephalopathy (2 deaths among 3 "
                    "neurological events in 9 patients, Morgan et al. 2013); MAGEA12 remains Unsafe."
                ),
            },
            {
                "clearance_id": "COC-20261123-NEG-TBI1301-LUNG",
                "entity_id": ctag2_eid,
                "protein_group_label": ctag2_lbl,
                "gene_id": "ENSG00000126890",
                "product_id": "tbi1301",
                "modality": "TCR-T",
                "status": "blocked_unresolved_event",
                "clearance_scope": "",
                "hla_matched_patients_treated": 9,
                "cleared_organs": [],
                "cleared_epitopes": [],
                "hla_alleles": prods["tbi1301"]["hla"],
                "cleared_transcripts": [],
                "uncleared_coding_transcripts": [],
                "source_keys": ["tbi1301_2022"],
                "trial_source_ids": [sources["tbi1301_2022"]["source_id"]],
                "rationale": (
                    "Product-isolation control: TBI-1301 (n=9) reported an unresolved Grade 3 steroid-responsive "
                    "lung injury and is never pooled or credited; CTAG2 lung clearance rests solely on lete_cel_c259."
                ),
            },
        ]
    )

    payload = {
        "schema": "oncoref.clinical-organ-clearances",
        "schema_version": 1,
        "policy_id": coc.CLINICAL_ORGAN_CLEARANCE_POLICY_ID,
        "ensembl_release": 115,
        "locus_attributions": {
            "CTAG2": {
                "id": "RCA-20261123-CTAG2-XQ28",
                "entity_id": ctag2_eid,
                "protein": "CTAG2",
                "gene_id": "ENSG00000126890",
                "locus": "X:154651972-154653579 (-)",
                "selected_transcript": "ENST00000369585",
                "uncleared_isoform": "ENST00000247306",
                "shared_kmer_fraction": 0.50856,
                "self_family_partners": ["CTAG1A (ENSG00000268651)", "CTAG1B (ENSG00000184033)"],
                "external_unsafe_overlap_pct": 0.40,
                "xq28_neighbors": [
                    "ENSG00000306286 (antisense + strand lncRNA X:154619409-154655372 spanning 100% of CTAG2)",
                    "GAB3 (same - strand protein_coding X:154675249-154751583, 21,670 bp upstream)",
                ],
                "brain_snrnaseq_resolution": (
                    "Brain snRNA-seq signal (3.18% L5 cortical neurons in 2 donors) is nucleus-only with "
                    "0.00 TPM bulk brain RNA, negative brain IHC (0/4), and paralog_ambiguous driven 100% "
                    "by closed self-family CTA NY-ESO-1 (CTAG1A/CTAG1B, itself 0.00 TPM brain bulk) plus "
                    "Xq28 GAB3 / antisense ENSG00000306286 transcription."
                ),
                "lung_scrnaseq_resolution": (
                    "Lung AT2 scRNA-seq signal (4.29% in 3 donors, 0.00 TPM bulk lung, 0.2 nTPM HPA lung) "
                    "is cleared at the epitope/isoform level by first-pass intravenous lete-cel c259 TCR-T "
                    "clinical safety (COC-20261123-CTAG2-LETE-CEL)."
                ),
            }
        },
        "entries": entries,
    }
    coc.validate(payload)
    out_path = DATA_DIR / coc.CLEARANCE_FILE
    out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    md = [
        "# Clinical pMHC-Directed Organ-Safety Evidence Tier (`#766`)",
        "",
        f"**Policy ID:** `{coc.CLINICAL_ORGAN_CLEARANCE_POLICY_ID}`  ",
        "**Ensembl Release:** `115` (`GRCh38`)",
        "",
        "## Guardrails & `tsarina` Isoform Contract",
        "",
        "- **Isoform- & epitope-scoped clearance (`CTAG2` `LAGE-1a` vs. `LAGE-1b`):**",
        "  `ENST00000369585` (`LAGE-1a`, `180 aa`) contains `SLLMWITQC` (`157–165`), whereas",
        "  `ENST00000247306` (`LAGE-1b`, `210 aa`) lacks `SLLMWITQC` and is explicitly listed in",
        "  `uncleared_coding_transcripts`.",
        "- **Vascular perfusion & BBB protection:** Systemic IV `TCR-T` clears `lung` and `urinary bladder`",
        "  (`CTAG2` also exceeds its `6.92 TPM` GTEx Urinary Bladder baseline by `87.3x` at `604.33 TPM` in `OV`),",
        "  while never clearing `brain` clinically (resolved separately via `#711` `Xq28` attribution) or `heart`",
        "  (`#544` confounded day-161 cardiac event preserved).",
        "",
        "## Reviewed Entries",
        "",
        "| Clearance ID | Target | Product | Status | Patients | Cleared Organs | Cleared Transcripts | Uncleared Transcripts |",
        "| --- | --- | --- | --- | ---: | --- | --- | --- |",
    ]
    for e in entries:
        md.append(
            f"| `{e['clearance_id']}` | `{e['protein_group_label']}` | `{e['product_id']}` | `{e['status']}` | "
            f"`{e['hla_matched_patients_treated']}` | `{','.join(e['cleared_organs']) or '-'}` | "
            f"`{len(e['cleared_transcripts'])}` | `{len(e['uncleared_coding_transcripts'])}` |"
        )
    (AUDIT_DIR / "index.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(
        f"Wrote {len(entries)} entries to {out_path}; CTAG2 cleared={ctag2_cl}, uncleared={ctag2_un}"
    )


if __name__ == "__main__":
    main()
