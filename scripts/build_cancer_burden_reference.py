#!/usr/bin/env python3
"""Rebuild ``oncoref/data/cancer-incidence-mortality.csv`` with GLOBOCAN 2024 world
counts, shares, and ICD-10 site provenance while preserving independent US ACS
Cancer Facts & Figures 2024 estimates and the frozen pirlygenes ``aggregation``
compatibility column.

Primary sources:
- World (GLOBOCAN 2024): Sung et al., *Global Cancer Statistics 2024: GLOBOCAN
  Estimates of Incidence and Mortality Worldwide for 34 Cancers in 186 Countries*,
  CA Cancer J Clin 2026; PMID:42417444; DOI:10.3322/caac.70090; PMC13343830
  (Table 1 and Data Sources and Methods; XML SHA-256
  ``60c96e5a6b8de9fea1fc225f22c4c440c77fa70038e14fc055987bf8b2fcc55c``).
- US (ACS 2024): Siegel et al., *Cancer statistics, 2024*, CA Cancer J Clin 2024;
  PMID:38230766 (Table 1: 2,001,140 total cases; 611,720 total deaths).
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "oncoref/data/cancer-incidence-mortality.csv"

WORLD_TOTAL_CASES = 20_636_778
WORLD_TOTAL_DEATHS = 9_762_507
WORLD_TOTAL_CASES_EXCL_NMSC = 19_498_464
WORLD_TOTAL_DEATHS_EXCL_NMSC = 9_703_786

US_TOTAL_CASES_2024 = 2_001_140
US_TOTAL_DEATHS_2024 = 611_720

# Exact GLOBOCAN 2024 Table 1 rows (PMC13343830, caac70090-tbl-0001):
# key -> (table_label, icd10_codes, cases, deaths)
GLOBOCAN_2024_TABLE1 = {
    "lung": ("Lung, bronchus, and trachea", "C33-C34", 2_637_005, 1_861_839),
    "breast": ("Female breast", "C50", 2_434_087, 693_660),
    "colorectum": ("Colorectum", "C18-C21", 2_041_007, 917_895),
    "prostate": ("Prostate", "C61", 1_546_112, 419_849),
    "stomach": ("Stomach", "C16", 980_286, 641_554),
    "thyroid": ("Thyroid", "C73", 959_281, 48_070),
    "liver": ("Liver and intrahepatic bile duct", "C22", 843_045, 732_489),
    "bladder": ("Bladder", "C67", 635_264, 227_626),
    "cervix": ("Cervix uteri", "C53", 604_196, 279_583),
    "non_hodgkin_lymphoma": ("Non-Hodgkin lymphoma", "C82-C86,C96", 562_793, 234_903),
    "pancreas": ("Pancreas", "C25", 531_318, 490_786),
    "leukemia": ("Leukemia", "C91-C95", 495_113, 292_234),
    "esophagus": ("Esophagus", "C15", 494_320, 441_994),
    "lip_oral_cavity": ("Lip, oral cavity", "C00-C06", 452_205, 194_108),
    "kidney": ("Kidney and renal pelvis", "C64-C65", 442_570, 144_871),
    "uterus_endometrium": ("Corpus uteri", "C54", 434_620, 100_680),
    "melanoma": ("Melanoma of skin", "C43", 337_969, 58_398),
    "ovary": ("Ovary", "C56", 330_731, 203_850),
    "brain_cns": ("Brain, central nervous system", "C70-C72", 324_095, 248_964),
    "multiple_myeloma": ("Multiple myeloma", "C88,C90", 196_157, 119_029),
    "larynx": ("Larynx", "C32", 193_937, 96_619),
    "gallbladder": ("Gallbladder", "C23", 126_384, 92_029),
    "nasopharynx": ("Nasopharynx", "C11", 119_750, 62_060),
    "oropharynx": ("Oropharynx", "C09-C10", 107_887, 51_823),
    "hypopharynx": ("Hypopharynx", "C12-C13", 91_464, 49_748),
    "hodgkin_lymphoma": ("Hodgkin lymphoma", "C81", 84_555, 19_515),
    "testicular_germ_cell": ("Testis", "C62", 70_863, 9_240),
    "salivary_glands": ("Salivary glands", "C07-C08", 55_229, 20_707),
    "vulva": ("Vulva", "C51", 48_894, 18_868),
    "penis": ("Penis", "C60", 34_765, 12_346),
    "kaposi_sarcoma": ("Kaposi sarcoma", "C46", 34_230, 10_727),
    "mesothelioma": ("Mesothelioma", "C45", 27_350, 22_715),
    "vagina": ("Vagina", "C52", 18_143, 7_484),
    "non_melanoma_skin": ("Nonmelanoma skin cancer", "C44", 1_138_314, 58_721),
}

HEAD_AND_NECK_KEYS = (
    "lip_oral_cavity",
    "salivary_glands",
    "oropharynx",
    "nasopharynx",
    "hypopharynx",
    "larynx",
)

SITE_METADATA = {
    "lung": (
        "Lung, bronchus, and trachea",
        "C33-C34",
        "trachea (C33); bronchus and lung (C34)",
        "mesothelioma (C45); thymus (C37)",
    ),
    "breast": (
        "Female breast",
        "C50",
        "female breast (C50)",
        "male breast (included in US ACS breast total; excluded from GLOBOCAN female breast row)",
    ),
    "prostate": ("Prostate", "C61", "prostate gland (C61)", ""),
    "colorectal": (
        "Colorectum",
        "C18-C21",
        "colon (C18); rectosigmoid junction and rectum (C19-C20); GLOBOCAN 2024 Table 1 also includes anus (C21)",
        "small intestine (C17); US ACS-CFF2024 tabulates anus (C21) separately",
    ),
    "pancreas": ("Pancreas", "C25", "pancreas (C25)", ""),
    "liver": (
        "Liver and intrahepatic bile duct",
        "C22",
        "hepatocellular carcinoma (C22.0); intrahepatic bile duct cholangiocarcinoma (C22.1); other primary liver (C22.2-C22.9)",
        "gallbladder (C23); extrahepatic bile duct and ampulla of Vater (C24)",
    ),
    "stomach": ("Stomach", "C16", "stomach (C16)", "esophagus (C15); small intestine (C17)"),
    "esophagus": ("Esophagus", "C15", "esophagus (C15)", "stomach cardia (C16.0)"),
    "bladder": (
        "Urinary bladder",
        "C67",
        "urinary bladder (C67)",
        "renal pelvis and ureter (C65-C66); urethra (C68)",
    ),
    "kidney": (
        "Kidney and renal pelvis",
        "C64-C65",
        "kidney parenchyma (C64); renal pelvis (C65)",
        "ureter (C66); urinary bladder (C67)",
    ),
    "melanoma": (
        "Melanoma of skin",
        "C43",
        "cutaneous malignant melanoma (C43)",
        "uveal/ocular melanoma (C69); mucosal melanoma; keratinocyte skin cancer (C44)",
    ),
    "ovary": (
        "Ovary",
        "C56",
        "ovary (C56)",
        "fallopian tube and adnexa (C57, untabulated in GLOBOCAN 34-site Table 1)",
    ),
    "uterus_endometrium": (
        "Corpus uteri",
        "C54",
        "corpus uteri / endometrium (C54); US ACS uterine corpus also includes uterus NOS (C55)",
        "cervix uteri (C53)",
    ),
    "cervix": ("Cervix uteri", "C53", "cervix uteri (C53)", "corpus uteri (C54)"),
    "head_and_neck": (
        "Lip/oral cavity, salivary glands, pharynx, and larynx",
        "C00-C13,C32",
        "lip and oral cavity (C00-C06); major salivary glands (C07-C08: 55,229 world cases, 20,707 world deaths; included in US SEER/ACS oral cavity & pharynx C00-C14); oropharynx (C09-C10); nasopharynx (C11); hypopharynx (C12-C13); larynx (C32)",
        "ill-defined lip/oral/pharynx (C14 untabulated in GLOBOCAN 2024 Table 1; included in US ACS C00-C14); nasal cavity and paranasal sinuses (C30-C31); thymus (C37)",
    ),
    "thyroid": ("Thyroid", "C73", "thyroid gland (C73)", "parathyroid gland (C75.0)"),
    "brain_cns": (
        "Brain and central nervous system",
        "C70-C72",
        "meninges (C70); brain (C71); spinal cord, cranial nerves, and other CNS (C72)",
        "pituitary and craniopharyngeal duct (C75.1-C75.2)",
    ),
    "non_hodgkin_lymphoma": (
        "Non-Hodgkin lymphoma",
        "C82-C86,C96",
        "follicular, diffuse, peripheral/cutaneous T-cell, and other NHL (C82-C86, C96)",
        "Hodgkin lymphoma (C81); multiple myeloma (C88, C90)",
    ),
    "hodgkin_lymphoma": ("Hodgkin lymphoma", "C81", "Hodgkin lymphoma (C81)", ""),
    "multiple_myeloma": (
        "Multiple myeloma and immunoproliferative diseases",
        "C88,C90",
        "immunoproliferative diseases (C88); multiple myeloma and plasma cell neoplasms (C90)",
        "",
    ),
    "leukemia_AML": (
        "Acute myeloid leukemia (subset of leukemia)",
        "C92.0",
        "acute myeloid leukemia (C92.0, subset of GLOBOCAN C91-C95 leukemia)",
        "lymphoid, chronic myeloid, and other leukemias (in leukemia_all_other)",
    ),
    "leukemia_all_other": (
        "Leukemia excluding AML",
        "C91-C95_excl_AML",
        "ALL, CLL, CML, and other non-AML leukemias (C91-C95 excluding AML)",
        "acute myeloid leukemia (in leukemia_AML)",
    ),
    "mesothelioma": (
        "Mesothelioma",
        "C45",
        "pleural, peritoneal, and other mesothelioma (C45)",
        "",
    ),
    "testicular_germ_cell": ("Testis", "C62", "testis (C62)", "extragonadal germ cell tumors"),
    "adrenal": (
        "Adrenal gland and paraganglia",
        "C74,C75.5",
        "adrenal cortex and medulla (C74); paraganglia (C75.5)",
        "",
    ),
    "bone_and_joint": (
        "Bones and articular cartilage",
        "C40-C41",
        "primary bone and joint neoplasms (C40-C41)",
        "soft tissue sarcoma (C47, C49); bone metastases",
    ),
    "soft_tissue_sarcoma": (
        "Connective and soft tissue (incl. heart)",
        "C38.0,C47,C49",
        "peripheral nerves (C47); connective and soft tissue (C49); heart (C38.0 in US ACS)",
        "primary bone sarcoma (C40-C41); Kaposi sarcoma (C46)",
    ),
    "gallbladder_biliary": (
        "Gallbladder (world Table 1 C23; US C23-C24)",
        "C23-C24",
        "gallbladder (C23); US ACS-CFF2024 also includes extrahepatic bile duct and ampulla of Vater (C24)",
        "intrahepatic bile duct (C22.1, tabulated under liver); GLOBOCAN 2024 34-site Table 1 leaves extrahepatic biliary (C24) in untabulated residual",
    ),
    "small_intestine": ("Small intestine", "C17", "small intestine (C17)", "colorectum (C18-C20)"),
    "anus": (
        "Anus and anal canal",
        "C21",
        "anus and anal canal (C21)",
        "overlapping in GLOBOCAN 2024 Table 1 with colorectal (C18-C21); separate in US ACS-CFF2024",
    ),
    "vulva": ("Vulva", "C51", "vulva (C51)", "vagina (C52)"),
    "vagina": ("Vagina", "C52", "vagina (C52)", "vulva (C51); cervix (C53)"),
    "penis": ("Penis", "C60", "penis (C60)", ""),
    "eye_ocular": ("Eye and orbit", "C69", "eye and adnexa (C69)", "cutaneous melanoma (C43)"),
    "kaposi_sarcoma": ("Kaposi sarcoma", "C46", "Kaposi sarcoma (C46)", ""),
    "other_and_unknown_primary": (
        "Untabulated sites and ill-defined/unknown primary (residual)",
        "C00-C97_residual",
        "all cancers (C00-C97) minus the 34 GLOBOCAN 2024 Table 1 rows (includes unknown primary C76-C80 plus untabulated minor sites such as C14, C17, C24, C30-C31, C37-C41, C47-C49, C57-C58, C63, C66, C68-C69, C74-C75)",
        "overlaps the separately curated minor-site rows (small_intestine, bone_and_joint, soft_tissue_sarcoma, eye_ocular, adrenal)",
    ),
    "non_melanoma_skin": (
        "Nonmelanoma skin cancer",
        "C44",
        "cutaneous squamous cell and other nonmelanoma skin cancer (C44)",
        "basal cell carcinoma excluded from GLOBOCAN C44 incidence; BCC and cSCC excluded from US reportable incidence",
    ),
}


def _pct(count: int, total: int) -> float:
    return round(100.0 * count / total, 1)


def build_frame() -> pd.DataFrame:
    # Verify exact GLOBOCAN 2024 Table 1 arithmetic before touching any row.
    nmsc_cases, nmsc_deaths = GLOBOCAN_2024_TABLE1["non_melanoma_skin"][2:]
    assert WORLD_TOTAL_CASES - nmsc_cases == WORLD_TOTAL_CASES_EXCL_NMSC
    assert WORLD_TOTAL_DEATHS - nmsc_deaths == WORLD_TOTAL_DEATHS_EXCL_NMSC

    sum_34_cases = sum(v[2] for v in GLOBOCAN_2024_TABLE1.values())
    sum_34_deaths = sum(v[3] for v in GLOBOCAN_2024_TABLE1.values())
    residual_cases = WORLD_TOTAL_CASES - sum_34_cases
    residual_deaths = WORLD_TOTAL_DEATHS - sum_34_deaths
    assert residual_cases == 1_202_839
    assert residual_deaths == 877_523

    df = pd.read_csv(CSV_PATH, dtype=str, keep_default_na=False)
    frozen_agg_sha = hashlib.sha256(
        "\n".join(
            f"{r.burden_category}\t{r.aggregation}"
            for r in df.sort_values("burden_category").itertuples()
        ).encode()
    ).hexdigest()
    assert frozen_agg_sha == "0a1e3af0cf10da3fb43e180ca0d0c6d0f7cd081f53af3a4c58f3bcab48cbd589"

    direct_map = {
        "lung": "lung",
        "breast": "breast",
        "prostate": "prostate",
        "colorectal": "colorectum",
        "pancreas": "pancreas",
        "liver": "liver",
        "stomach": "stomach",
        "esophagus": "esophagus",
        "bladder": "bladder",
        "kidney": "kidney",
        "melanoma": "melanoma",
        "ovary": "ovary",
        "uterus_endometrium": "uterus_endometrium",
        "cervix": "cervix",
        "thyroid": "thyroid",
        "brain_cns": "brain_cns",
        "non_hodgkin_lymphoma": "non_hodgkin_lymphoma",
        "hodgkin_lymphoma": "hodgkin_lymphoma",
        "multiple_myeloma": "multiple_myeloma",
        "mesothelioma": "mesothelioma",
        "testicular_germ_cell": "testicular_germ_cell",
        "vulva": "vulva",
        "vagina": "vagina",
        "penis": "penis",
        "kaposi_sarcoma": "kaposi_sarcoma",
        "non_melanoma_skin": "non_melanoma_skin",
    }

    df = df.set_index("burden_category", drop=False)

    for cat, (labels, codes, inc_sites, exc_sites) in SITE_METADATA.items():
        df.at[cat, "source_site_labels"] = labels
        df.at[cat, "source_site_codes"] = codes
        df.at[cat, "included_source_sites"] = inc_sites
        df.at[cat, "excluded_source_sites"] = exc_sites

    # 1. Direct single-row GLOBOCAN 2024 categories
    for cat, gkey in direct_map.items():
        tlabel, icd, cases, deaths = GLOBOCAN_2024_TABLE1[gkey]
        df.at[cat, "world_incidence_count"] = str(cases)
        df.at[cat, "world_incidence_total"] = str(WORLD_TOTAL_CASES)
        df.at[cat, "world_incidence_pct"] = f"{_pct(cases, WORLD_TOTAL_CASES):.1f}"
        df.at[cat, "world_mortality_count"] = str(deaths)
        df.at[cat, "world_mortality_total"] = str(WORLD_TOTAL_DEATHS)
        df.at[cat, "world_mortality_pct"] = f"{_pct(deaths, WORLD_TOTAL_DEATHS):.1f}"
        df.at[cat, "world_source_locator"] = f"PMC13343830 Table 1 ({tlabel}, {icd})"
        df.at[cat, "world_source_locator_status"] = "extracted"
        df.at[cat, "rounding_rule"] = "round_half_even_1dp_over_all_cancers_incl_nmsc"
        if cat == "kaposi_sarcoma":
            df.at[cat, "source"] = "GLOBOCAN2024;literature"
            df.at[cat, "source_anchor"] = (
                "PMID:42417444;DOI:10.3322/caac.70090;DOI:10.1016/S2214-109X(23)00349-2"
            )
            df.at[cat, "provenance_notes"] = (
                "World counts and shares from GLOBOCAN 2024 Table 1; US values remain curated rare-cancer approximations"
            )
        else:
            df.at[cat, "source"] = "ACS-CFF2024;GLOBOCAN2024"
            df.at[cat, "source_anchor"] = "PMID:38230766;PMID:42417444;DOI:10.3322/caac.70090"

    df.at["liver", "notes"] = (
        "liver+intrahepatic bile duct; world>>US (GLOBOCAN2024 843k=4.1% cases, 732k=7.5% deaths)"
    )
    df.at["thyroid", "notes"] = (
        "lowest mortality:incidence; overdiagnosis (GLOBOCAN2024 959k=4.6% cases, 48k=0.5% deaths)"
    )

    # 2. Six-site head_and_neck aggregate (resolves salivary C07-C08 inclusion on both US & World)
    hn_cases = sum(GLOBOCAN_2024_TABLE1[k][2] for k in HEAD_AND_NECK_KEYS)
    hn_deaths = sum(GLOBOCAN_2024_TABLE1[k][3] for k in HEAD_AND_NECK_KEYS)
    assert hn_cases == 1_020_472 and hn_deaths == 475_065
    df.at["head_and_neck", "us_incidence_count"] = str(58_450 + 12_650)
    df.at["head_and_neck", "us_incidence_total"] = str(US_TOTAL_CASES_2024)
    df.at["head_and_neck", "us_mortality_count"] = str(12_230 + 3_880)
    df.at["head_and_neck", "us_mortality_total"] = str(US_TOTAL_DEATHS_2024)
    df.at["head_and_neck", "us_source_locator"] = (
        "ACS-CFF2024 Table 1: Oral cavity & pharynx (58,450 cases; 12,230 deaths; SEER C00-C14 includes salivary glands C07-C08) + Larynx (12,650 cases; 3,880 deaths)"
    )
    df.at["head_and_neck", "us_source_locator_status"] = "extracted"
    df.at["head_and_neck", "world_incidence_count"] = str(hn_cases)
    df.at["head_and_neck", "world_incidence_total"] = str(WORLD_TOTAL_CASES)
    df.at["head_and_neck", "world_incidence_pct"] = f"{_pct(hn_cases, WORLD_TOTAL_CASES):.1f}"
    df.at["head_and_neck", "world_mortality_count"] = str(hn_deaths)
    df.at["head_and_neck", "world_mortality_total"] = str(WORLD_TOTAL_DEATHS)
    df.at["head_and_neck", "world_mortality_pct"] = f"{_pct(hn_deaths, WORLD_TOTAL_DEATHS):.1f}"
    df.at["head_and_neck", "world_source_locator"] = (
        "PMC13343830 Table 1 sum of 6 sites: Lip/oral cavity (C00-C06) + Salivary glands (C07-C08) + Oropharynx (C09-C10) + Nasopharynx (C11) + Hypopharynx (C12-C13) + Larynx (C32)"
    )
    df.at["head_and_neck", "world_source_locator_status"] = "extracted"
    df.at["head_and_neck", "source"] = "ACS-CFF2024;GLOBOCAN2024"
    df.at["head_and_neck", "source_anchor"] = "PMID:38230766;PMID:42417444;DOI:10.3322/caac.70090"
    df.at["head_and_neck", "rounding_rule"] = "round_half_even_1dp_over_all_cancers_incl_nmsc"
    df.at["head_and_neck", "notes"] = (
        "oral+salivary+pharynx+larynx aggregate (incl C07-C08 salivary & C32 larynx; world also includes C11 nasopharynx)"
    )

    # 3. Gallbladder / biliary (C23 in GLOBOCAN 2024 Table 1; C23-C24 in US ACS)
    gb_cases, gb_deaths = GLOBOCAN_2024_TABLE1["gallbladder"][2:]
    df.at["gallbladder_biliary", "world_incidence_count"] = str(gb_cases)
    df.at["gallbladder_biliary", "world_incidence_total"] = str(WORLD_TOTAL_CASES)
    df.at["gallbladder_biliary", "world_incidence_pct"] = f"{_pct(gb_cases, WORLD_TOTAL_CASES):.1f}"
    df.at["gallbladder_biliary", "world_mortality_count"] = str(gb_deaths)
    df.at["gallbladder_biliary", "world_mortality_total"] = str(WORLD_TOTAL_DEATHS)
    df.at["gallbladder_biliary", "world_mortality_pct"] = (
        f"{_pct(gb_deaths, WORLD_TOTAL_DEATHS):.1f}"
    )
    df.at["gallbladder_biliary", "world_source_locator"] = (
        "PMC13343830 Table 1 (Gallbladder, C23; extrahepatic biliary C24 untabulated in 34-site table)"
    )
    df.at["gallbladder_biliary", "world_source_locator_status"] = "extracted_partial_scope"
    df.at["gallbladder_biliary", "source"] = "ACS-CFF2024;GLOBOCAN2024"
    df.at["gallbladder_biliary", "source_anchor"] = (
        "PMID:38230766;PMID:42417444;DOI:10.3322/caac.70090"
    )
    df.at["gallbladder_biliary", "rounding_rule"] = "round_half_even_1dp_over_all_cancers_incl_nmsc"
    df.at["gallbladder_biliary", "provenance_notes"] = (
        "World counts reflect GLOBOCAN 2024 Gallbladder (C23); US ACS-CFF2024 covers Gallbladder & other biliary (C23-C24). Intrahepatic bile duct (C22.1) is in liver in both sources."
    )

    # 4. Leukemia split (parent C91-C95 extracted in GLOBOCAN 2024; subtype counts not split in Table 1)
    leuk_cases, leuk_deaths = GLOBOCAN_2024_TABLE1["leukemia"][2:]
    leuk_locator = f"PMC13343830 Table 1 parent Leukemia (C91-C95: {leuk_cases:,} cases = 2.4%, {leuk_deaths:,} deaths = 3.0%); subtype split approximated from US/hematologic proportions"
    for lcat, wi, wm in (("leukemia_AML", "0.9", "1.5"), ("leukemia_all_other", "1.5", "1.5")):
        df.at[lcat, "world_incidence_pct"] = wi
        df.at[lcat, "world_mortality_pct"] = wm
        df.at[lcat, "world_incidence_total"] = str(WORLD_TOTAL_CASES)
        df.at[lcat, "world_mortality_total"] = str(WORLD_TOTAL_DEATHS)
        df.at[lcat, "world_source_locator"] = leuk_locator
        df.at[lcat, "world_source_locator_status"] = "parent_extracted_subtype_split_approximated"
        df.at[lcat, "source"] = "ACS-CFF2024;GLOBOCAN2024"
        df.at[lcat, "source_anchor"] = "PMID:38230766;PMID:42417444;DOI:10.3322/caac.70090"
        df.at[lcat, "rounding_rule"] = "subtype_shares_sum_to_parent_1dp"

    # 5. Residual bucket (All cancers C00-C97 minus the 34 GLOBOCAN 2024 Table 1 rows)
    df.at["other_and_unknown_primary", "world_incidence_count"] = str(residual_cases)
    df.at["other_and_unknown_primary", "world_incidence_total"] = str(WORLD_TOTAL_CASES)
    df.at["other_and_unknown_primary", "world_incidence_pct"] = (
        f"{_pct(residual_cases, WORLD_TOTAL_CASES):.1f}"
    )
    df.at["other_and_unknown_primary", "world_mortality_count"] = str(residual_deaths)
    df.at["other_and_unknown_primary", "world_mortality_total"] = str(WORLD_TOTAL_DEATHS)
    df.at["other_and_unknown_primary", "world_mortality_pct"] = (
        f"{_pct(residual_deaths, WORLD_TOTAL_DEATHS):.1f}"
    )
    df.at["other_and_unknown_primary", "world_source_locator"] = (
        "PMC13343830 Table 1: All cancers (20,636,778 cases; 9,762,507 deaths) minus sum of 34 tabulated rows (19,433,939 cases; 8,884,984 deaths)"
    )
    df.at["other_and_unknown_primary", "world_source_locator_status"] = "extracted"
    df.at["other_and_unknown_primary", "source"] = "ACS-CFF2024;GLOBOCAN2024"
    df.at["other_and_unknown_primary", "source_anchor"] = (
        "PMID:38230766;PMID:42417444;DOI:10.3322/caac.70090"
    )
    df.at["other_and_unknown_primary", "rounding_rule"] = (
        "round_half_even_1dp_over_all_cancers_incl_nmsc"
    )
    df.at["other_and_unknown_primary", "provenance_notes"] = (
        "World residual equals All cancers (C00-C97) minus the 34 GLOBOCAN 2024 Table 1 rows and therefore overlaps separately curated minor-site rows (small_intestine, bone_and_joint, soft_tissue_sarcoma, eye_ocular, adrenal) and extrahepatic biliary (C24)."
    )

    # 6. Six minor categories not separately tabulated in GLOBOCAN 2024's 34-site Table 1:
    # preserve their existing vintage rather than mixing 2022/literature shares with 2024 counts.
    for mcat in (
        "adrenal",
        "bone_and_joint",
        "soft_tissue_sarcoma",
        "small_intestine",
        "anus",
        "eye_ocular",
    ):
        df.at[mcat, "world_source_locator_status"] = "not_in_globocan2024_34site_table"
        if mcat == "anus":
            df.at[mcat, "provenance_notes"] = (
                "Retained GLOBOCAN 2022 standalone C21 share; GLOBOCAN 2024 Table 1 folds anus (C21) into colorectum (C18-C21) — do not sum anus and colorectal world shares."
            )

    return df.reset_index(drop=True)


def main() -> None:
    out = build_frame()
    out.to_csv(CSV_PATH, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
