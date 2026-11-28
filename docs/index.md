# oncoref

`oncoref` provides curated reference data for cancer genomics and antigen discovery:
cancer-type ontology, RNA-seq cohorts, clean-TPM normalization, tumor mutational burden
(TMB), checkpoint-inhibitor (ICI) response, normal-tissue injury tiers, and proteome-wide
<span class="cat cat-ca">CA</span>, <span class="cat cat-ra">RA</span>, and
<span class="cat cat-cta">CTA</span> classification.

[Interactive antigen explorer](antigen-explorer.html) ·
[Antigen whitepaper](antigen-whitepaper.md) ·
[Multi-species whitepaper](species-antigen-whitepaper.md) ·
[Classification rules](target-selection.md) ·
[Patient profiles & filters](target-filters.md) ·
[CLI & builds](curation.md) ·
[Sequences & peptide subtraction](antigen-sequences.md) ·
[Tissue injury tiers](tissue-tiers.md) ·
[API guide](api.md)

## Antigen Classification

Across **19,893** human protein groups (**20,099** Ensembl 115 coding loci), OncoRef
classifies **322 CTA, 53 Cancer Antigen only, 146 Reproductive Antigen only, 1,116 Other and
18,256 Unsafe**—using the same mathematical screening architecture across **human**, **domestic dog** (**118** CTAs), and **laboratory mouse** (**8** *in vivo* CTAs):

1. **Normal-tissue safety screen:** Bulk RNA (GTEx v11, HPA v23), sorted and activated immune subsets,
   single-cell clusters, donor-resolved rare-cell UMI atlases (CELLxGENE + Siletti), adult cornea scRNA,
   normal protein staining (HPA IHC), and reviewed normal-tissue literature exclude proteins expressed in healthy
   adult organs as <span class="cat cat-unsafe">Unsafe</span>. Proteins with missing coverage, unresolved paralog signals,
   unscreened sensory micro-niches, or `Unsafe` sequence overlap land in <span class="cat cat-other">Other</span> rather than being treated as safe.
2. **Symmetric 10× enrichment:** Against one shared non-thymic healthy baseline
   **H₀ = max(Normal N, Disposable D, 1 TPM)**,
   <span class="cat cat-ca">Cancer Antigens (CA)</span> require **≥ 10 × H₀** in a QC-passing *in vivo* tumor sample and
   <span class="cat cat-ra">Reproductive Antigens (RA)</span> require **≥ 10 × H₀** in a reproductive or embryonic sample.
   Proteins passing both routes form the <span class="cat cat-cta">CTA</span> set (**311/322 = 96.6%** with ectopic activation outside homologous source organs).
3. **Patient-specific & cross-sex filtering:** Biological <span class="cat cat-cta">CTA</span>, <span class="cat cat-ca">CA</span>,
   and <span class="cat cat-ra">RA</span> labels record where a protein is expressed, not what a particular patient can
   safely lose. Downstream patient protection (`antigens.table("CTA", sex="female" | "male")` or `--sex female|male`) marks reproductive compartments as
   <span class="cat cat-protect">Protect</span>, <span class="cat cat-loss">Accept loss</span>, or
   <span class="cat cat-absent">Absent</span>—partitioning all **322** human CTAs into four mutually exclusive host-physiology classes:
   **222** `CTA-M` (`270` non-pregnant female-safe `CTA-M + CTA-P`), **12** `CTA-F` (`60` male-safe `CTA-F + CTA-P`),
   **48** adult-absent pregnancy/peripartum `CTA-P`, and **40** shared `CTA-MF`
   ([cross-sex CTAs](antigen-whitepaper.md#sex-restricted-reproductive-antigens-and-cross-sex-ctas),
   [multi-species whitepaper](species-antigen-whitepaper.md),
   [rules and limits](host-selection.md)).

![Antigen classification overview](audits/target-selection-20261110/rebuild/summary/source-to-tiers-architecture.png)

[Interactive antigen explorer](antigen-explorer.html) ·
[Antigen whitepaper](antigen-whitepaper.md) ·
[Summary figures](audits/target-selection-20261110/rebuild/summary.md) ·
[Data & downloads](antigen-data.md) ·
[Donor-level rare-cell RNA](audits/normal-cell-atlas-20261110/rare-cell/index.md)

## Install & Quick Start

```bash
pip install oncoref
```

```python
from oncoref import antigens, cancer_ontology, ici_response

# Cross-species summary, patient-safe antigen filtering, and peptide subtraction
antigens.summary()
female_safe_ctas = antigens.table("CTA", sex="female", ectopic_only=True)
dog_ctas = antigens.table("CTA", species="dog", replicated_only=True)
antigens.explain("PRAME", sex="female")
unique_peptides = antigens.subtract("MAGEA4", sex="female", max_mismatches=1)

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
| Antigen classification | `oncoref.antigens`, `oncoref.antigen_sets`, `oncoref.target_selection`, `oncoref.peptide_subtraction` | Cross-species `CA`, `RA`, `CTA`, `Other`, `Unsafe` discovery, patient sex/organ protection, Ensembl sequences, exact/edit-distance peptide subtraction |
| Tissue injury tiers | `oncoref.tissue_tiers`, `oncoref.target_curation` | Cited bulk and cell-level injury tiers (`Vital`, `Protected`, `Manageable`, `Disposable`), acute/persistent loss consequences, curated normal-tissue homes |
| Cancer vocabulary | `oncoref.cancer_ontology`, `oncoref.cohorts` | Registry records, aliases, hierarchy, subtype axes, MMR/MSI status, cohort IDs, matched normal tissues |
| Expression & QC | `oncoref.expression`, `oncoref.normalization`, `oncoref.library_qc`, `oncoref.samples` | Source-matrix ingestion, percentiles, clean TPM, housekeeping normalization, library-contamination screens, patient/donor deduplication |
| Burden & ICI response | `oncoref.ici_response`, `oncoref.tmb`, `oncoref.incidence`, `oncoref.fusions` | Anti-PD-1 and broader ICI response estimates, TMB, incidence/mortality, defining fusions |
| Genes, species & cache | `oncoref.gene_matching`, `oncoref.gene_ids`, `oncoref.human`, `oncoref.canine`, `oncoref.mouse`, `oncoref.catalog`, `oncoref.reference_data` | Canonical Ensembl identities, evaluated human/canine/murine antigen catalogs, SHA-256-verified reference caches |
