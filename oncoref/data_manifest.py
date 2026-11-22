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

"""The complete, classified inventory of oncoref-domain data.

This is the single declarative source of truth for *what data oncoref owns*
and how each dataset is held. ``catalog.inventory()`` exposes that current
self-owned inventory. A separately named, version-pinned pirlygenes migration
snapshot preserves the June 2026 takeover audit without pretending to follow
later additions or removals in another repository.

Ownership rule:
  oncoref owns empirical base facts and reference mechanics: canonical gene IDs,
  cancer ontology, reference expression/normalization, epidemiology, TMB, ICI/aPD1,
  source-anchored CTA facts, and measurement/provenance tables. pirlygenes owns
  purpose-specific gene sets/panels and target-to-therapy registries. trufflepig
  owns per-sample rules and QC interpretation.

Held buckets (oncoref domain):
  WHEEL     — small base-layer tables shipped in the wheel (no fetch). A
              ``legacy-compat`` category marks shipped historical surfaces that
              remain importable but are frozen. ``legacy-dataset-dispositions``
              records the replacement owner/surface for every such table.
  BUNDLE    — heavy expression artifacts in the version-pinned release tarball.
  HPA       — Human Protein Atlas reference tables fetched per-source on demand.
  SOURCE    — the raw per-sample TPM matrices (a separate large optional bundle).
  PLANNED   — oncoref-domain fact tables still to port.
  SUPERSEDED— a pirlygenes table oncoref replaced with its own regenerated one.

OUT_OF_SCOPE — pirlygenes or trufflepig data that is NOT oncoref's domain:
  target-selection panels, marker/signature gene sets, therapy modality catalogs,
  sample-rule tables, and artifact/QC interpretation rules.

The catalog (:mod:`oncoref.catalog`) manages the fetchable subset; this module
is the inventory + classification.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

#: {name: (category, description)} — small curated tables shipped in the wheel.
WHEEL: dict[str, tuple[str, str]] = {
    "cancer-type-registry": ("ontology", "cancer-type codes, hierarchy, family/tissue"),
    "cancer-subtype-groupings": ("ontology", "cross-cutting MSI/POLE/HPV/MYCN axes"),
    "cancer-cohort-aggregates": ("ontology", "computed histology rollups (SARC_RMS, …)"),
    "cancer-lineage-groups": ("ontology", "registry family -> coarse histogenesis lineage group"),
    "cancer-lineage-group-overrides": ("ontology", "per-code lineage-group overrides (NBL)"),
    "cohort-registry": ("ontology", "first-class cohort vocabulary"),
    "cancer-fusions": ("ontology", "characteristic driver fusions per cancer type"),
    "cancer-code-burden-map": ("ontology", "anatomic burden-category overrides"),
    "cancer-incidence-mortality": ("epidemiology", "incidence/mortality by burden category"),
    "cancer-tmb": ("genomics", "median tumor mutational burden per type"),
    "cancer-apd1-response": ("response", "anti-PD-1 monotherapy ORR per type"),
    "therapy-benefit-toxicity-evidence": (
        "therapy-evidence",
        "source-anchored clinical benefit/toxicity fact rows",
    ),
    "cancer-reference-expression-samples": ("expression", "per-sample curation manifest"),
    "expression_sources": ("expression", "cohort expression-source registry"),
    # normalization references (R-norm; consumed by the clean TPM engine)
    "housekeeping-genes": ("normalization", "housekeeping panel for normalization"),
    "censored-gene-reference-tpm": (
        "normalization",
        "legacy symbol projection of the global clean-TPM censored profile",
    ),
    "clean-tpm-censored-genes": (
        "normalization",
        "global censored-gene membership and clean-PolyA reference composition",
    ),
    "histone-genes": ("normalization", "histone-cluster genes"),
    "ribosomal-protein-genes": ("normalization", "ribosomal-protein genes"),
    "ribosomal-protein-pseudogenes": ("normalization", "ribosomal-protein pseudogenes"),
    "mitochondrial-genes": ("normalization", "mtDNA-encoded genes"),
    "rrna-and-pseudogenes": ("normalization", "rRNA + rRNA pseudogenes"),
    "numt-pseudogenes": ("normalization", "NUMT-like nuclear-mito pseudogenes"),
    "small-noncoding-rnas": ("normalization", "small non-coding RNA loci"),
    "nuclear-retained-lncrnas": ("normalization", "polyA-bias nuclear-retained lncRNAs"),
    "hemoglobin-genes": ("normalization", "hemoglobin genes"),
    # gene-id / symbol resolution (R-resolve)
    "ensembl-id-aliases": ("gene-id", "alt-haplotype → primary-contig Ensembl gene ids"),
    "ncbi-symbol-synonyms": ("gene-id", "NCBI gene symbol synonyms → official symbol"),
    "extra-tx-mappings": ("gene-id", "supplemental transcript→gene mappings"),
    "cdna-identical-gene-groups": ("gene-id", "cDNA-identical gene groups"),
    "proteoform-collapse-overrides": ("gene-id", "manual proteoform-collapse overrides"),
    # Legacy compatibility tables imported before the stack boundary was clarified.
    # Keep shipped/readable for old callers, but do not treat them as the future
    # oncoref ownership surface for marker panels, target selection, or sample rules.
    "cancer-key-genes": (
        "legacy-compat",
        "legacy per-type biomarker/target rows; future purpose-specific panels live downstream",
    ),
    "cancer-driver-genes": (
        "legacy-compat",
        "frozen Bailey Table S1 schema; rows migrated to driver-gene-evidence",
    ),
    "cancer-driver-variants": (
        "legacy-compat",
        "frozen Bailey Table S4 schema; rows migrated to driver-variant-evidence",
    ),
    "cancer-type-genes": (
        "legacy-compat",
        "legacy role-stratified per-type gene rows; not a general panel-ownership surface",
    ),
    "cancer-viral-antigens": (
        "legacy-compat",
        "legacy per-oncovirus targetable-antigen rows; target-selection curation lives downstream",
    ),
    "disease-state-rules": (
        "legacy-compat",
        "legacy narrative disease-state rules; per-sample rule ownership belongs downstream",
    ),
    "narrative-gene-sets": (
        "legacy-compat",
        "legacy named narrative gene sets; purpose-specific gene-set ownership belongs to pirlygenes",
    ),
    "degenerate-subtype-pairs": ("ontology", "expression-degenerate subtype pairs"),
    "rare-cancer-fusion-rules": (
        "legacy-compat",
        "legacy direct fusion rules for rare cancers; per-sample interpretation belongs downstream",
    ),
    "fusion-surrogate-expression": (
        "legacy-compat",
        "legacy fusion-expression surrogate rows; sample-rule use belongs downstream",
    ),
    "fusion-expression-effects": (
        "legacy-compat",
        "legacy downstream-expression rules per fusion; sample-rule use belongs downstream",
    ),
    # expression-source metadata + genomics (R-exprmeta)
    "cancer-expression-source-candidates": ("expression", "candidate expression sources per type"),
    "cancer-frameshift-burden": ("genomics", "per-type frameshift-indel burden"),
}

#: {name: (category, description)} — pirlygenes-originated heavy artifacts in
#: the release tarball.
BUNDLE_PIRLYGENES: dict[str, tuple[str, str]] = {
    "cancer-reference-expression": ("expression", "per-cohort RNA-seq summary shards"),
    "cancer-reference-expression-percentiles": ("expression", "per-gene percentile vectors"),
    "cancer-reference-expression-representatives": ("expression", "per-cohort medoid samples"),
    "pan-cancer-expression": ("expression", "pan-cancer HPA-tissue + TCGA matrix"),
    "hpa-cell-type-expression": ("hpa", "HPA cell-type nTPM matrix"),
}

#: {name: (category, description)} — oncoref-originated bundle artifacts produced
#: by the rebuilt expression-artifact pipeline. These are not in the frozen
#: pirlygenes data inventory, but they are part of oncoref's downstream bundle
#: contract.
BUNDLE_ORIGINATED: dict[str, tuple[str, str]] = {
    "cancer-reference-expression-within-sample-top5": (
        "expression",
        "within-sample high-expression prevalence by cohort",
    ),
    "cancer-reference-expression-percentiles-proteoform-cta": (
        "expression",
        "CTA-scope proteoform-summed percentile vectors",
    ),
    "cancer-reference-expression-within-sample-top5-proteoform-cta": (
        "expression",
        "CTA-scope proteoform-summed within-sample prevalence",
    ),
    "source-matrix-sample-qc": (
        "expression",
        "per-sample QC manifest for rebuilt expression artifacts",
    ),
    "expression-artifact-build-metadata": (
        "expression",
        "per-cohort build/QC provenance for rebuilt expression artifacts",
    ),
    "expression-artifact-build-metadata.json": (
        "expression",
        "bundle-level build/QC provenance for rebuilt expression artifacts",
    ),
}

#: {name: (category, description)} — artifacts migrated from Trufflepig. The
#: importer and release provenance live in oncoref, but the pinned TCGA
#: decomposition output and legacy subtype source table originated downstream.
BUNDLE_TRUFFLEPIG: dict[str, tuple[str, str]] = {
    "tcga-deconvolved-expression.csv.gz": (
        "expression",
        "TCGA tumor-attributed TPM summaries migrated from pinned Trufflepig output",
    ),
    "subtype-deconvolved-expression.csv.gz": (
        "expression",
        "source-separated tumor-reference subtype and cohort summaries",
    ),
    "tumor-reference-expression-provenance": (
        "expression",
        "source-level derivation metadata for tumor-reference summary artifacts",
    ),
}

#: {name: (category, description)} — every heavy artifact in the release tarball.
#: Names match ``data_bundle.DOWNLOADABLE_PATHS`` (sans extension for ``.csv``);
#: kept in sync by a test.
BUNDLE: dict[str, tuple[str, str]] = {
    **BUNDLE_PIRLYGENES,
    **BUNDLE_ORIGINATED,
    **BUNDLE_TRUFFLEPIG,
}

#: {name: (category, description)} — HPA tables fetched per-source (not pirlygenes files).
HPA: dict[str, tuple[str, str]] = {
    "hpa_rna_consensus": ("hpa", "HPA RNA consensus per-tissue nTPM"),
    "hpa_normal_tissue": ("hpa", "HPA IHC protein detection per tissue"),
    "hpa_single_cell": ("hpa", "HPA single-cell-type RNA nTPM"),
    "hpa_cancer_ihc": ("hpa", "HPA cancer IHC scored-patient counts"),
    "hpa_cancer_rna": ("hpa", "HPA cancer RNA pTPM counts and prevalence by cohort"),
    "hpa_cancer_rna_samples": ("hpa", "HPA cancer sample pTPM source (gzip; 1.36 GB)"),
}

#: The raw per-sample TPM matrices — every derived artifact is built from these.
SOURCE: dict[str, tuple[str, str]] = {
    "per-sample-tpm-matrices": ("expression", "raw per-sample cohort TPM (build inputs)"),
}

#: {name: (category, description)} — oncoref-domain fact tables still to port.
PLANNED: dict[str, tuple[str, str]] = {}

#: pirlygenes tables oncoref replaced with its own regenerated equivalent.
SUPERSEDED: dict[str, str] = {
    "protein-identical-gene-groups": "proteoform-groups-genome (byte-identical, regenerated)",
    "cta-protein-groups": "proteoform-groups (byte-identical; pirlygenes' is ≥90% identity)",
}

#: {name: (category, description)} — oncoref-ORIGINATED wheel tables: derived or
#: regenerated here rather than copied from pirlygenes, so they aren't in the
#: pirlygenes snapshot but DO ship in the wheel and belong in the inventory.
CANCERDATA_ORIGINATED: dict[str, tuple[str, str]] = {
    "antigen-source-inventory": (
        "antigens",
        "Versioned input, assay and comparison files used by the antigen catalog",
    ),
    "cta-literature-benchmarks": (
        "antigens",
        "Independent CTA classification controls with exact-gene evidence and coverage limitations",
    ),
    "marcinkowski2019-ct83-expression": (
        "expression",
        "Tumor-localized CT83 expression by histology and assay, Figure 4",
    ),
    "transcript-selection": (
        "protein",
        "Pinned MANE/canonical representative decisions, exact versions, evidence holds and support flags",
    ),
    "gene-knowledge": (
        "protein",
        "Reviewed gene findings with exact identities, tissue/state, evidence limits and screening consequences",
    ),
    "normal-induction-evidence": (
        "expression",
        "Cited normal induction contexts; qualitative reproductive-association counterevidence",
    ),
    "external-ihc-evidence": (
        "protein",
        "categorical external staining observations, source evidence, reagent identity and explicit validation gaps",
    ),
    "hpa-cancer-crosswalk": ("hpa", "IHC/RNA cancer correspondence and scope limitations"),
    "hpa-cancer-rna-cohorts": ("hpa", "HPA TCGA/validation RNA labels and sample counts"),
    "legacy-dataset-dispositions": (
        "governance",
        "reviewed owner/replacement and freeze policy for every legacy-compat dataset",
    ),
    "rna-protein-hpa-prior-sources": (
        "protein-calibration",
        "checksum-pinned HPA v23 RNA and normal-tissue IHC weak-prior sources",
    ),
    "rna-protein-hpa-priors": (
        "protein-calibration",
        "gene-specific HPA v23 normal-tissue RNA/IHC weak priors",
    ),
    "rna-protein-calibration-samples": (
        "protein-calibration",
        "matched CPTAC patient manifest used by the RNA/protein calibration models",
    ),
    "rna-protein-calibration-sources": (
        "protein-calibration",
        "checksum-pinned matched CPTAC RNA/protein source files and measurement scales",
    ),
    "hpa-safety-tissue-map": (
        "hpa",
        "source/version-aware conceptual-to-native safety-tissue mappings",
    ),
    "normal-tissue-taxonomy": (
        "cta",
        "cited species-aware reproductive and conditional disposable normal tissue scopes",
    ),
    "normal-tissue-mappings": (
        "cta",
        "versioned source tissue/cell-type mappings to the shared anatomical scopes",
    ),
    "cancer-entity-driver-spectrum": (
        "genomics",
        "structured observed driver distributions for heterogeneous cancer entities",
    ),
    "driver-gene-evidence": (
        "genomics",
        "source-anchored Bailey Table S1 cancer-driver gene evidence",
    ),
    "driver-variant-evidence": (
        "genomics",
        "source-anchored Bailey Table S4 pan-cancer recurrent driver variants",
    ),
    "cancer-reference-sample-molecular-provenance": (
        "genomics",
        "sample-level diagnosis and orthogonal molecular evidence provenance",
    ),
    "cancer-reference-expression-availability": (
        "expression",
        "one-row-per-source reference-expression availability manifest",
    ),
    "cancer-testis-antigens": ("cta", "CTA definition (HPA tissue-restriction over candidates)"),
    "cancer-testis-antigens-extended": (
        "cta",
        "extended reproductive CTA definition over the same candidate universe",
    ),
    "cta-gene-publication-evidence": (
        "cta",
        "targeted primary gene evidence with exact locations, assay resolution and validation limits",
    ),
    "cta-publication-membership": (
        "cta",
        "complete CTA landscape and placental nominations with source locations and canonical mapping",
    ),
    "cta-publication-sources": (
        "cta",
        "publication citations, exact source tables, checksums and evidence scope",
    ),
    "cta-intake-identity-reviews": (
        "cta-evidence",
        "authority-backed identities without canonical Ensembl mappings",
    ),
    "cta-trial-evidence": (
        "cta-evidence",
        "primary clinical-trial and approved-therapy target evidence with alias scope",
    ),
    "cta-external-safety-reviews": (
        "cta-evidence",
        "external normal-tissue injury and receptor cross-reactivity evidence reviews",
    ),
    "cta-screening-review-overrides": (
        "cta-evidence",
        "policy-scoped reviews of historical CTA holdbacks",
    ),
    "cta-clinical-targets": (
        "cta",
        "explicit clinical/canonical CTA target tier joined to HPA evidence",
    ),
    "cta-specificity-audit": (
        "cta",
        "machine-readable CTA specificity demotion and candidate-only audit decisions",
    ),
    "cta-candidate-references": (
        "cta",
        "referenced CTA candidate watchlist with historical citations for subsequently assessed rows",
    ),
    "cta-warning-reviews": (
        "cta",
        "reviewed CTA discovery exceptions with unresolved normal-tissue evidence",
    ),
    "cta-reviewed-evidence": (
        "cta",
        "source-anchored supplemental per-modality observations for CTA candidates",
    ),
    "cta-donor-rna": (
        "cta-evidence",
        "reviewed individual HPA v23 vital-tissue RNA donors; currently CTAG2 heart only",
    ),
    "cta-donor-ihc": (
        "cta-evidence",
        "supplementary CTAG2 heart IHC donor metadata recording RNA/specimen mismatch",
    ),
    "cta-policy-findings": (
        "cta-evidence",
        "reference-bound reviewed CTA findings, acceptance status, citations and limitations",
    ),
    "cta-research-gene-sets": (
        "cta",
        "historical 215-gene expression/clinical research shortlist; see cta-membership-gene-sets for biological membership",
    ),
    "cta-selection-human": (
        "cta",
        "current research selection across the full human protein-group reference, with decisions, member identities and retained flags",
    ),
    "cta-selection-sources": (
        "cta-evidence",
        "species-scoped evidence roles, native units, source integration status and coverage gaps",
    ),
    "cta-protein-universe": (
        "cta-evidence",
        "coding or established-ncORF status per gene, the stage-0 universe for CTA selection, with evidence and exclusion reason",
    ),
    "target-selection": (
        "cta",
        "Cancer and Reproductive Antigen catalog for every human protein record, with set memberships, evidence receipts and paper links",
    ),
    "host-exposure": (
        "cta-evidence",
        "per-compartment host RNA/protein observations and coverage, independent of antigen membership and explicit tissue-loss choices",
    ),
    "tissue-marker-panels": (
        "expression-qc",
        "GTEx v11 protein-coding tissue marker panels and native TPM reference levels for conservative tumor-library QC",
    ),
    "tgct-reproductive-localization-audit": (
        "expression-qc",
        "somatic Leydig/Sertoli niche and spermatid stoichiometry audit across all 154 TCGA-TGCT libraries (#713)",
    ),
    "tgct-aafl-protein-support": (
        "cta-evidence",
        "cross-cohort witness sensitivity audit for all 34 catalog proteins peaking in TCGA-2G-AAFL-01 (#713)",
    ),
    "nose-skin-reference": (
        "cta-evidence",
        "normal nose and skin expression of each Target candidate: olfactory mucosa bulk RNA (Olender 2016) and adult nose, skin, scalp and tongue single cells (CELLxGENE Census), with the tissue recorded when several people or datasets agree",
    ),
    "readthrough-loci": (
        "cta-evidence",
        "read-through loci whose every coding transcript joins two separate protein-coding genes, with parent genes; the catalog's sequence check moves them to Other",
    ),
    "rare-cell-evidence": (
        "cta-evidence",
        "per-donor rare normal-cell RNA grades with multi-molecule support, study replication, mapping-risk flags and bulk confirmation",
    ),
    "normal-expression-literature": (
        "cta-evidence",
        "reviewed publications on normal adult expression for every Target candidate and rare-cell dependent, with rubric outcome, recorded tissue tiers, family rules, URLs and quotes",
    ),
    "gene-mapping-risk": (
        "cta-evidence",
        "per-gene shared transcript sequence and same-strand host overlap that can mis-assign sparse single-cell reads",
    ),
    "rare-cell-attribution-reviews": (
        "cta-evidence",
        "central Ensembl 115 exon-geometry and paralog 31-mer attribution adjudications for rare normal-cell RNA review candidates",
    ),
    "cta-ihc-compartments": (
        "cta-evidence",
        "source-specific RNA/IHC compartment correspondences with citations and coverage limits",
    ),
    "tumor-protein-evidence": (
        "cta-evidence",
        "review ledger for reliable tumor protein support with assay, source and identity scope",
    ),
    "tumor-restriction-reviews": (
        "cta-evidence",
        "normal-context coverage reviews for tumor-supported reconsideration leads",
    ),
    "cta-publication-context": (
        "cta-evidence",
        "author-recognized prior CT antigens with protein-group identities and source locations",
    ),
    "cta-membership-evidence": (
        "cta-evidence",
        "reference-bound reviewed reproductive/developmental and tumor evidence with primary citations",
    ),
    "cta-membership-controls": (
        "cta-evidence",
        "positive biological controls and negative comparators; evaluation labels only",
    ),
    "cta-membership-source-roles": (
        "cta-evidence",
        "publication evidence roles distinguishing cancer/reproductive nominations from normal-only enrichment",
    ),
    "cta-research-clinical-evidence": (
        "cta-evidence",
        "reviewed clinical qualification and comparator evidence bound to normal-finding hashes",
    ),
    "cta-antibody-adjudications": (
        "cta-evidence",
        "completed primary-antibody adjudications with unresolved conflicts and reconsideration requirements",
    ),
    "normal-cell-publication-evidence": (
        "cta-evidence",
        "cited human reproductive/tumor support and normal-cell protein counterevidence",
    ),
    "cta-followup-notes": (
        "cta-evidence",
        "review questions and reconsideration requirements for watched CTA candidates; no decision overrides",
    ),
    "cta-followup-reviews": (
        "cta-evidence",
        "current decisions, native-unit RNA maxima and evidence for watched CTA candidates",
    ),
    "sequence-gene-rna-reviews": (
        "cta-evidence",
        "RNA exclusions and coverage for exact reference sequence loci outside the canonical gene table",
    ),
    "proteoform-groups": ("gene-id", "byte-identical CTA proteoform groups"),
    "proteoform-groups-genome": ("gene-id", "byte-identical proteoform groups (genome-wide)"),
    "source-matrices": (
        "expression",
        "per-cohort raw-matrix registry (code/source/n_samples/release)",
    ),
    "expression-artifact-gene-universe-deltas": (
        "expression",
        "known pirlygenes/oncoref expression artifact row-universe deltas",
    ),
    "representative-source-adjudications": (
        "expression",
        "reviewed source-sample roles and benchmark eligibility for representative artifacts",
    ),
    "ncbi-entrez-gene-mappings": (
        "gene-id",
        "filtered NCBI Entrez/GeneID -> canonical Ensembl gene mappings",
    ),
    "tissue-burden-map": ("ontology", "primary_tissue -> anatomic burden category"),
    "family-burden-map": ("ontology", "registry family -> burden category (fallback)"),
    "cancer-burden-scope": (
        "ontology",
        "multi-bucket and subsite anatomic-burden scope metadata for mixed cohorts",
    ),
    "cancer-mismatch-repair-statuses": (
        "ontology",
        "cross-cancer MMRd/MSI-H vs pMMR/MSS-like classifier state semantics",
    ),
    "cancer-who-soft-tissue-bone-audit": (
        "ontology",
        "WHO soft-tissue/bone category and behavior parity audit",
    ),
    "cta-ihc-unreliable": ("cta", "CTAs whose HPA IHC is cross-reactive (RNA-only fallback)"),
    "clean-tpm-biological-housekeeping-genes": (
        "normalization",
        "HPA-stable biological housekeeping denominator candidates for clean TPM",
    ),
    "cancer-response-signatures": (
        "legacy-compat",
        "transitional response-signature surface; new/extended therapy panels belong in pirlygenes",
    ),
    "cancer-fusion-citation-audit": (
        "ontology",
        "per-fusion primary-literature citation audit behind cancer-fusions",
    ),
    "cancer-ici-response": (
        "response",
        "curated ICI objective response rates by cancer code and regimen",
    ),
    "cancer-ici-response-estimates": (
        "response",
        "trial-level ICI response estimates, denominators and confidence intervals",
    ),
    "cancer-ici-source-locator-audit": (
        "response",
        "primary-source table/figure locator audit for curated ICI response rows",
    ),
    "cancer-key-gene-citation-audit": (
        "governance",
        "per-gene primary-literature citation audit behind legacy cancer-key-genes",
    ),
    "cancer-normal-tissue-map": (
        "ontology",
        "cancer type -> matched normal tissue of origin map",
    ),
    "cancer-tmb-source-audit": (
        "genomics",
        "per-cancer-type primary-source audit behind cancer-tmb",
    ),
    "canonical-gene-space": (
        "gene-id",
        "authoritative primary-assembly Ensembl canonical gene-ID universe",
    ),
    "gse294016-sample-diagnoses": (
        "expression",
        "curated per-sample diagnosis mapping for GSE294016 cohorts",
    ),
}

#: pirlygenes/trufflepig data that is NOT oncoref's domain — target selection /
#: therapy modality catalogs, marker/signature panels, surfaceome, and per-sample
#: interpretation/rule tables. These stay downstream and key to oncoref IDs.
OUT_OF_SCOPE: frozenset[str] = frozenset(
    {
        # therapeutic modality / drug data
        "ADC-approved",
        "ADC-trials",
        "ADC-withdrawn",
        "bispecific-antibodies-approved",
        "CAR-T-approved",
        "TCR-T-approved",
        "TCR-T-trials",
        "multispecific-tcell-engager-trials",
        "radioligand-targets",
        # surfaceome / targetability
        "surface-proteins",
        "cancer-surfaceome",
        # therapy signatures / evidence
        "therapy-response-signatures",
        "estimate-signatures",
        "immune-receptor-segments",
        # analysis / QC gene sets
        "tme-markers",
        "stem-cell-marker-panels",
        "culture-stress-genes",
        "ffpe-sensitive-markers",
        "degradation-gene-pairs",
        "mutation-expression-effects",
        "rare-cancer-rna-surrogates",
        # differential-expression signatures
        "tumor-up-vs-matched-normal",
        "heme-tumor-up-vs-matched-normal",
        # gene-set panels (analysis / target selection)
        "cancer-lineage-panels",
        "cancer-family-panels",
        "lineage-genes",
        "gene-sets",
        # build QC
        "artifact-expectations",
    }
)

#: Provenance for :data:`PIRLYGENES_MIGRATION_SNAPSHOT`. This exact commit's
#: ``pirlygenes/data`` tree normalizes to the snapshot below. It is historical
#: migration evidence, not a live cross-repository completeness promise.
PIRLYGENES_MIGRATION_SNAPSHOT_METADATA: dict[str, str] = {
    "repository": "pirl-unc/pirlygenes",
    "ref": "v5.22.80",
    "commit": "bb6679e39029218e93acb2e54317a8d4ff98f424",
    "data_path": "pirlygenes/data",
    "captured_at": "2026-06-12",
}

#: Exact normalized inventory of ``pirlygenes/data`` at the pinned migration
#: commit above. Use ``scripts/audit_pirlygenes_inventory.py`` to reproduce this
#: snapshot or compare it with a newer pirlygenes ref.
PIRLYGENES_MIGRATION_SNAPSHOT: frozenset[str] = frozenset(
    {
        "ADC-approved",
        "ADC-trials",
        "ADC-withdrawn",
        "artifact-expectations",
        "bispecific-antibodies-approved",
        "cancer-apd1-response",
        "cancer-code-burden-map",
        "cancer-cohort-aggregates",
        "cancer-driver-genes",
        "cancer-driver-variants",
        "cancer-expression-source-candidates",
        "cancer-family-panels",
        "cancer-frameshift-burden",
        "cancer-fusions",
        "cancer-incidence-mortality",
        "cancer-key-genes",
        "cancer-lineage-group-overrides",
        "cancer-lineage-groups",
        "cancer-lineage-panels",
        "cancer-reference-expression",
        "cancer-reference-expression-percentiles",
        "cancer-reference-expression-representatives",
        "cancer-reference-expression-samples",
        "cancer-subtype-groupings",
        "cancer-surfaceome",
        "cancer-tmb",
        "cancer-type-genes",
        "cancer-type-registry",
        "cancer-viral-antigens",
        "CAR-T-approved",
        "cdna-identical-gene-groups",
        "censored-gene-reference-tpm",
        "clean-tpm-censored-genes",
        "cohort-registry",
        "cta-protein-groups",
        "culture-stress-genes",
        "degenerate-subtype-pairs",
        "degradation-gene-pairs",
        "disease-state-rules",
        "ensembl-id-aliases",
        "estimate-signatures",
        "expression_sources",
        "extra-tx-mappings",
        "ffpe-sensitive-markers",
        "fusion-expression-effects",
        "fusion-surrogate-expression",
        "gene-sets",
        "heme-tumor-up-vs-matched-normal",
        "hemoglobin-genes",
        "histone-genes",
        "housekeeping-genes",
        "hpa-cell-type-expression",
        "immune-receptor-segments",
        "lineage-genes",
        "mitochondrial-genes",
        "multispecific-tcell-engager-trials",
        "mutation-expression-effects",
        "narrative-gene-sets",
        "ncbi-symbol-synonyms",
        "nuclear-retained-lncrnas",
        "numt-pseudogenes",
        "pan-cancer-expression",
        "protein-identical-gene-groups",
        "proteoform-collapse-overrides",
        "radioligand-targets",
        "rare-cancer-fusion-rules",
        "rare-cancer-rna-surrogates",
        "ribosomal-protein-genes",
        "ribosomal-protein-pseudogenes",
        "rrna-and-pseudogenes",
        "small-noncoding-rnas",
        "stem-cell-marker-panels",
        "surface-proteins",
        "TCR-T-approved",
        "TCR-T-trials",
        "therapy-benefit-toxicity-evidence",
        "therapy-response-signatures",
        "tme-markers",
        "tumor-up-vs-matched-normal",
    }
)


def captured() -> set[str]:
    """Datasets already held or shipped for compatibility."""
    return set(WHEEL) | set(BUNDLE) | set(SUPERSEDED)


def in_scope() -> set[str]:
    """Base-layer inventory: captured datasets plus still-planned fact tables."""
    return captured() | set(PLANNED)


RETENTION_TIERS: tuple[str, ...] = (
    "core_runtime",
    "active_calibration_or_audit",
    "archival_compatibility_snapshot",
)

#: Exact file stems under ``oncoref/data/`` belonging to frozen historical CTA
#: releases (v1.19–v1.45 ``oncoref.cta`` surfaces) or superseded legacy-compat tables.
_ARCHIVAL_EXACT_STEMS: frozenset[str] = frozenset(
    {
        "cancer-testis-antigens",
        "cancer-testis-antigens-extended",
        "cta-antibody-adjudications",
        "cta-candidate-references",
        "cta-clinical-targets",
        "cta-donor-ihc",
        "cta-donor-rna",
        "cta-external-safety-reviews",
        "cta-finding-migrations",
        "cta-followup-notes",
        "cta-followup-reviews",
        "cta-gene-publication-evidence",
        "cta-ihc-compartments",
        "cta-ihc-unreliable",
        "cta-intake-identity-reviews",
        "cta-literature",
        "cta-membership-controls",
        "cta-membership-evidence",
        "cta-membership-gene-sets",
        "cta-membership-source-roles",
        "cta-policy-findings",
        "cta-protein-universe",
        "cta-publication-context",
        "cta-publication-membership",
        "cta-publication-sources",
        "cta-quantitative-gene-sets",
        "cta-rare-cell-screen",
        "cta-research-clinical-evidence",
        "cta-research-gene-sets",
        "cta-reviewed-evidence",
        "cta-screen-results",
        "cta-screening-review-overrides",
        "cta-selection",
        "cta-selection-audit",
        "cta-selection-audit-human",
        "cta-selection-human",
        "cta-selection-profiles-human",
        "cta-selection-sources",
        "cta-specificity-audit",
        "cta-warning-reviews",
        "target-gene-sets",
    }
)

#: Heavy empirical calibration / qualification / audit tables used by build-time
#: verification or optional calibration APIs rather than default runtime queries.
_CALIBRATION_OR_AUDIT_STEMS: frozenset[str] = frozenset(
    {
        "translation-qualification",
        "rna-protein-calibrations",
        "rna-protein-calibration-samples",
        "rna-protein-calibration-sources",
        "rna-protein-hpa-priors",
        "rna-protein-hpa-prior-sources",
        "curation-build-receipt",
        "normal-expression-literature",
        "tgct-reproductive-localization-audit",
        "tgct-aafl-protein-support",
    }
)


#: Heavy archival compatibility and detailed calibration/qualification tables
#: excluded from the PyPI wheel (`pyproject.toml` `exclude-package-data`) and
#: managed via `oncoref.reference_data` (`reference-releases.json`, #738).
EXTERNALIZED_ORIGINATED: dict[str, tuple[str, str]] = {
    "translation-qualification": (
        "protein",
        "Pinned Ensembl peptide records with primary-locus CDS and sequence qualification, source identities and review reasons",
    ),
    "rna-protein-calibrations": (
        "protein-calibration",
        "cohort-specific matched CPTAC RNA/protein detection and quantitative models",
    ),
    "cta-selection-audit-human": (
        "cta-evidence",
        "Frozen 1.36.0 human CTA selection audit explanation, membership, overlap, funnel, and criteria tables",
    ),
    "target-gene-sets": (
        "cta",
        "disjoint CTA, Unsafe and Neither gene sets; only Unsafe is the forbidden background",
    ),
    "cta-membership-gene-sets": (
        "cta",
        "historical literature-inclusive biological CTA/oncofetal catalogue and separate expression-risk flags",
    ),
    "cta-quantitative-gene-sets": (
        "cta",
        "experimental uniform quantitative CTA/developmental candidates, measured gates, coverage and retained expression-risk flags",
    ),
    "cta-selection-profiles-human": (
        "protein",
        "Current protein evidence attributes: native-unit RNA, tumor support and adult Disposable pattern",
    ),
    "cta-rare-cell-screen": (
        "cta-evidence",
        "unreviewed single-cell screen of rare normal-cell expression for CTA and Held proteins, with source rules and limits",
    ),
    "cta-screen-results": (
        "cta",
        "ordered CTA research decisions; Variable genes are excluded pending reconsideration",
    ),
}


def _stem_for_data_file(rel_path: str) -> str:
    name = rel_path.rsplit("/", 1)[-1]
    for suffix in (".csv.gz", ".json.gz", ".tsv.gz", ".csv", ".json", ".yaml", ".yml"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def package_data_size_audit() -> pd.DataFrame:
    """Classify every ``oncoref/data/`` file by byte size, wheel status, and retention tier (#738).

    Partitions all data files into three explicit tiers:

    - ``core_runtime``: ontology, gene-ID canonicalization, clean-TPM normalization,
      source registries, and the active ``target-selection`` / multi-species catalogs.
    - ``active_calibration_or_audit``: heavy translation-qualification, CPTAC/HPA
      RNA-protein calibration models, and curation audit tables.
    - ``archival_compatibility_snapshot``: frozen ``legacy-compat`` tables and
      historical ``v1.19–v1.45`` ``oncoref.cta`` snapshots retained for reproducibility.
    """
    from pathlib import Path

    import pandas as pd

    data_root = Path(__file__).resolve().parent / "data"
    rows = []
    for path in sorted(p for p in data_root.rglob("*") if p.is_file()):
        rel = path.relative_to(data_root).as_posix()
        stem = _stem_for_data_file(rel)
        size_bytes = path.stat().st_size
        shipped_in_wheel = stem not in EXTERNALIZED_ORIGINATED

        if stem in WHEEL:
            bucket = "WHEEL"
            category = WHEEL[stem][0]
        elif stem in CANCERDATA_ORIGINATED:
            bucket = "CANCERDATA_ORIGINATED"
            category = CANCERDATA_ORIGINATED[stem][0]
        elif stem in EXTERNALIZED_ORIGINATED:
            bucket = "EXTERNALIZED_ORIGINATED"
            category = EXTERNALIZED_ORIGINATED[stem][0]
        elif rel.startswith("species/"):
            bucket = "SPECIES_DATA"
            category = "species"
        else:
            bucket = "METADATA_SIDECAR"
            category = "provenance"

        if (
            category == "legacy-compat"
            or stem in _ARCHIVAL_EXACT_STEMS
            or rel == "species/human/cta-reference-mappings.json.gz"
        ):
            tier = "archival_compatibility_snapshot"
            packaging = (
                "externalized_reference_release"
                if not shipped_in_wheel
                else "candidate_for_legacy_archive_bundle"
            )
        elif stem in _CALIBRATION_OR_AUDIT_STEMS or category == "protein-calibration":
            tier = "active_calibration_or_audit"
            packaging = (
                "externalized_reference_release"
                if not shipped_in_wheel
                else "candidate_for_optional_calibration_bundle"
            )
        else:
            tier = "core_runtime"
            packaging = "keep_in_wheel"

        rows.append(
            {
                "relative_path": rel,
                "dataset_stem": stem,
                "manifest_bucket": bucket,
                "category": category,
                "retention_tier": tier,
                "recommended_packaging": packaging,
                "shipped_in_wheel": bool(shipped_in_wheel),
                "size_bytes": int(size_bytes),
            }
        )
    df = pd.DataFrame(rows).sort_values(["size_bytes", "relative_path"], ascending=[False, True])
    return df.reset_index(drop=True)


def package_data_size_summary() -> pd.DataFrame:
    """Summarize packaged ``oncoref/data/`` file counts and byte totals by retention tier (#738)."""
    audit = package_data_size_audit()
    total_bytes = int(audit["size_bytes"].sum()) or 1
    grouped = (
        audit.groupby(["retention_tier", "recommended_packaging"], as_index=False)
        .agg(
            n_files=("relative_path", "count"),
            total_bytes=("size_bytes", "sum"),
        )
        .sort_values("total_bytes", ascending=False)
        .reset_index(drop=True)
    )
    grouped["total_mb"] = (grouped["total_bytes"] / 1_000_000.0).round(3)
    grouped["share_of_wheel_data"] = (grouped["total_bytes"] / total_bytes).round(4)
    return grouped
