# oncoref

`oncoref` provides curated reference data for cancer genomics and antigen discovery:
cancer-type ontology, RNA-seq cohorts, clean-TPM normalization, tumor mutational burden
(TMB), checkpoint-inhibitor (ICI) response, normal-tissue injury tiers, and proteome-wide
<span class="cat cat-ca">CA</span>, <span class="cat cat-ra">RA</span>, and
<span class="cat cat-cta">CTA</span> classification.

[Antigen whitepaper](antigen-whitepaper.md) ·
[Interactive catalog](audits/target-selection-20261110/rebuild/report.html) ·
[Classification rules](target-selection.md) ·
[Patient profiles & filters](target-filters.md) ·
[CLI & builds](curation.md) ·
[Sequences & peptide subtraction](antigen-sequences.md) ·
[Tissue injury tiers](tissue-tiers.md) ·
[API guide](api.md)

## Antigen Classification

Across **19,893** human protein groups (**20,099** Ensembl 115 coding loci), OncoRef
classifies **332 CTA, 117 Cancer Antigen only, 133 Reproductive Antigen only, 1,062 Other and
18,249 Unsafe**:

1. **Normal-tissue safety screen:** Bulk RNA (GTEx v11, HPA v23), sorted immune subsets (HPA, DICE, Monaco),
   single-cell clusters, donor-resolved rare-cell UMI atlases (CELLxGENE + Siletti), adult cornea scRNA,
   normal protein staining (HPA IHC), and reviewed normal-tissue literature exclude proteins expressed in healthy
   adult organs as <span class="cat cat-unsafe">Unsafe</span>. Proteins with missing coverage, unresolved paralog signals,
   or `Unsafe` sequence overlap land in <span class="cat cat-other">Other</span> rather than being treated as safe.
2. **Symmetric 10× enrichment:** Against one shared non-thymic healthy baseline
   **H₀ = max(Normal N, Disposable D, 1 TPM)**,
   <span class="cat cat-ca">Cancer Antigens (CA)</span> require **≥ 10 × H₀** in a QC-passing tumor sample and
   <span class="cat cat-ra">Reproductive Antigens (RA)</span> require **≥ 10 × H₀** in a reproductive or embryonic sample.
   Proteins passing both routes form the <span class="cat cat-cta">CTA</span> set.
3. **Patient-specific & cross-sex filtering:** Biological <span class="cat cat-cta">CTA</span>, <span class="cat cat-ca">CA</span>,
   and <span class="cat cat-ra">RA</span> labels record where a protein is expressed, not what a particular patient can
   safely lose. Downstream patient profiles (`--profile` / `target_filters`) mark reproductive compartments as
   <span class="cat cat-protect">Protect</span>, <span class="cat cat-loss">Accept loss</span>, or
   <span class="cat cat-absent">Absent</span>—identifying **276** female-safe CTAs (`CTA-M + CTA-E`: `232 + 44`) and **50** male-safe
   CTAs (`CTA-F + CTA-E`: `6 + 44`) under strict protection of all opposite-sex adult organs
   ([cross-sex CTAs](antigen-whitepaper.md#sex-restricted-reproductive-antigens-and-cross-sex-ctas),
   [rules and limits](host-selection.md)).

![Antigen classification overview](audits/target-selection-20261110/rebuild/summary/source-to-tiers-architecture.png)

[Antigen whitepaper](antigen-whitepaper.md) ·
[Interactive protein catalog](audits/target-selection-20261110/rebuild/report.html) ·
[Summary figures](audits/target-selection-20261110/rebuild/summary.md) ·
[Data & downloads](antigen-data.md) ·
[Donor-level rare-cell RNA](audits/normal-cell-atlas-20261110/rare-cell/index.md)

## Install & Quick Start

```bash
pip install oncoref
```

```python
from oncoref import antigen_sets, cancer_ontology, ici_response

# Biological selection and consistent decision / sequence views
selection = antigen_sets.select("CTA", mode="standard")
ctas = selection.proteins()
decisions = selection.decisions()
forbidden_bg = antigen_sets.proteins("Unsafe|Other")

# Cancer ontology and ICI response
crc_msi = cancer_ontology.cancer_type_records(subtype_group="MSI", under="CRC")
ici_response.best_available_ici_response("COAD_MSI")
```

## Stack Boundary

- **oncoref**: empirical facts, source-linked curation, and canonical identifiers
  (cancer codes, Ensembl gene IDs, reference expression and normalization, epidemiology, TMB,
  ICI response, tissue-injury tiers, and `CA` / `RA` / `CTA` / `Unsafe` classifications).
- [**pirlygenes**](https://github.com/pirl-unc/pirlygenes): curated gene
  panels (lineage, TME, surfaceome, therapy signatures) keyed to `oncoref` IDs.
- [**trufflepig**](https://github.com/pirl-unc/trufflepig): per-sample QC,
  deconvolution, scoring, and sample-level interpretation.

## Data Domains

| Domain | Preferred modules | Use for |
| --- | --- | --- |
| Antigen classification | `oncoref.target_selection`, `oncoref.target_filters`, `oncoref.antigen_sets`, `oncoref.peptide_subtraction` | Proteome-wide `CA`, `RA`, `CTA`, `Other`, `Unsafe` calls, patient-specific organ filters, Ensembl 115 sequences, exact/edit-distance peptide subtraction |
| Tissue injury tiers | `oncoref.tissue_tiers`, `oncoref.target_curation` | Cited bulk and cell-level injury tiers (`Vital`, `Protected`, `Manageable`, `Disposable`), acute/persistent loss consequences, curated normal-tissue homes |
| Cancer vocabulary | `oncoref.cancer_ontology`, `oncoref.cohorts` | Registry records, aliases, hierarchy, subtype axes, MMR/MSI status, cohort IDs, matched normal tissues |
| Expression & QC | `oncoref.expression`, `oncoref.normalization`, `oncoref.library_qc`, `oncoref.samples` | Source-matrix ingestion, percentiles, clean TPM, housekeeping normalization, library-contamination screens, patient/donor deduplication |
| Burden & ICI response | `oncoref.ici_response`, `oncoref.tmb`, `oncoref.incidence`, `oncoref.fusions` | Anti-PD-1 and broader ICI response estimates, TMB, incidence/mortality, defining fusions |
| Genes, species & cache | `oncoref.gene_matching`, `oncoref.gene_ids`, `oncoref.human`, `oncoref.canine`, `oncoref.catalog`, `oncoref.reference_data` | Canonical Ensembl identities, species-scoped datasets, SHA-256-verified reference caches |
