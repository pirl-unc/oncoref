# cancerdata

[![Tests](https://github.com/pirl-unc/cancerdata/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/cancerdata/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/cancerdata.svg)](https://pypi.org/project/cancerdata/)

Curated cancer reference data — cancer-type ontology, tumor mutational burden
(TMB), incidence/mortality, and per-cohort RNA-seq expression — behind one small
Python API and a data fetch/cache CLI.

`cancerdata` sits at the **bottom of the openvax/PIRL stack**: it depends only on
pandas/numpy/pyarrow and is consumed by `pirlygenes` (gene-set curation and
analysis) and `tsarina` (personalized target selection). The small curated
tables ship in the wheel; the heavy per-cohort expression bundle downloads on
first use from the matching GitHub Release.

> **Status:** early. Ships the cancer-type ontology, TMB, incidence/mortality
> tables, and the per-cohort expression **percentile** read-accessors over a
> lazy-downloaded data bundle, plus the full fetch/cache CLI. The within-sample
> percentile signal lands in a following milestone.

## Install

```bash
pip install cancerdata
```

## Python API

```python
import cancerdata as cd

cd.resolve_cancer_type("prostate")        # -> "PRAD"
cd.cancer_type_info("SARC_RMS_ARMS")      # full registry record + burden + tmb
cd.cancer_tmb("LUAD_EGFR")                # 6.9  (inherited from LUAD)
cd.cancer_burden("pancreas", metric="us_mortality_pct")
cd.burden_category("SARC_OS")             # -> "bone_and_joint"
cd.cohort_aggregate_members("SARC")       # pan-sarcoma grand union

# Per-cohort expression percentiles (downloads the data bundle on first use):
cd.available_percentile_cohorts()         # 118 cohorts with per-sample data
cd.cohort_gene_percentiles("PRAD")        # per-gene p0…p100 vector (within-cohort)
cd.within_sample_top_fraction("PRAD")     # per-gene frac of samples where it's
                                          # top-5% expressed (within-sample)
```

### Domains

- **Ontology** — `cancer_type_registry`, `resolve_cancer_type`,
  `cancer_type_info`, `cancer_types_in_family`, `cancer_type_subtypes_of`,
  `viral_status`, `fusion_status`, the cohort vocabulary (`cohort_registry`,
  `cohort_aggregates`).
- **TMB** — `cancer_tmb`, `cancer_tmb_df` (parent-chain inheritance).
- **Incidence / mortality** — `cancer_burden`, `cancer_burden_df`,
  `burden_category` (ACS Cancer Facts & Figures / GLOBOCAN).

## CLI

```bash
cancerdata cancer-type prostate     # registry info as JSON
cancerdata tmb LUAD_EGFR            # 6.9
cancerdata burden pancreas --metric us_mortality_pct

# data bundle (per-cohort expression):
cancerdata fetch                   # download the ~340 MB bundle
cancerdata status                  # which bundle paths are cached (no download)
cancerdata cache-dir               # where the data bundle is cached
cancerdata prune --yes             # delete stale version caches
cancerdata version
```

## Development

```bash
./develop.sh   # editable install with dev extras
./format.sh    # ruff format
./lint.sh      # ruff check + format --check
./test.sh      # lint + pytest with coverage
```

## License

Apache 2.0.
