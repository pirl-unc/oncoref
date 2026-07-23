# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Canonical source identities for reference-expression summary rows."""

TREEHOUSE_TCGA_LEGACY_COHORT = "TREEHOUSE_POLYA_25_01_TCGA_SUBSET"
TREEHOUSE_TCGA_SAMPLES_COHORT = "TREEHOUSE_POLYA_25_01_TCGA_SAMPLES"
TREEHOUSE_TCGA_SARC_HISTOLOGY_COHORT = "TREEHOUSE_POLYA_25_01_TCGA_SARC_HISTOLOGY"

# These rows are derived by joining TCGA-SARC samples to GDC primary diagnoses.
# They are not members of the generic project-level TCGA sample cohort.
TREEHOUSE_TCGA_SARC_HISTOLOGY_CODES = frozenset({"SARC_DDLPS", "SARC_PLEOLPS", "SARC_WDLPS"})


def canonical_treehouse_tcga_summary_cohort(cancer_code: str) -> str:
    """Return the exact cohort for a legacy Treehouse TCGA summary row."""

    if cancer_code in TREEHOUSE_TCGA_SARC_HISTOLOGY_CODES:
        return TREEHOUSE_TCGA_SARC_HISTOLOGY_COHORT
    return TREEHOUSE_TCGA_SAMPLES_COHORT
