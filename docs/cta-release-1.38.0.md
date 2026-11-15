# OncoRef 1.38.0: one cancer rule and one set of curation records

**A Cancer Antigen now needs one tumor sample at 10 times every non-reproductive normal
value, and nothing else. The curation that the 1.37.0 filters kept in their own files
is part of the catalog, proteins of tissues the panels do not resolve need
reproductive evidence, and the filters are down to seven knobs. Catalog policy
`oncoref.target-selection.2026-11-15.1` has 670 Target candidates, down from 902.**
The 1.36.0 classification snapshot, sequence bundle and expression data (5.23.27) are
unchanged. These are outputs of screening rules, not validated tumor antigens.

[Catalog rules](target-selection.md) · [filter rules](target-filters.md) ·
[tissue tiers](tissue-tiers.md) · [web summary](audits/target-filters-20261112/report.html)

## Counts

| Category | 1.37.0 | 1.38.0 |
|---|---:|---:|
| CTA | 626 | 453 |
| Cancer Antigen only | 245 | 37 |
| Reproductive Antigen only | 31 | 180 |
| Other | 750 | 970 |
| Unsafe | 18,241 | 18,253 |

| Filter mode | 1.37.0 | 1.38.0 |
|---|---:|---:|
| `all` | 902 | 670 |
| `loose` | 473 | 452 |
| `standard` | 297 | 295 |
| `strict` | 145 | 145 |

## One cancer rule

**Cancer RNA qualifies when one tumor sample has `T ≥ 10 × max(N, D, 1)`:** N is the
highest non-reproductive normal value (tissues, brain subregions, immune cells), D the
highest Disposable value. Any of the 129 comparable TPM cohort views counts. The
two-patient route is gone; patient counts stay on each row as description. Germ cell
tumors no longer count, because their germline program is lineage rather than
cancer-specific activation. Each row has a new `tumor_fold` column.

- **All 24 familiar CTAs that were CTAs stay CTAs.** The lowest is AKAP4 at 16x;
  PRAME is at 1,231x, NY-ESO-1 at 944x and MAGE-A4 at 3,202x.
- **153 CTAs become Reproductive Antigen only.** 138 have no tumor at 10x, mostly
  low-expressed testis genes, and 15 qualified only through germ cell tumors (DAZL,
  PRDM14, PGK2, PRM2 and others). 19 of them are on the CTpedia list, among them
  MAGEB3, MAGEB10, SPANXN5 and CABS1. They stay Target candidates.
- **81 Cancer Antigens only become Other:Expression,** and 32 proteins held for
  unresolved patient identity under the old two-patient route are now Other:Expression.
- One sample can admit a protein. One Hodgkin lymphoma sample with no documented
  patient (GSE120328) is the best tumor for 31 Cancer Antigens; `tumor_top_sample`
  shows it.

## Home tissue the panels do not resolve

**A protein whose recorded normal home is olfactory epithelium, nasal mucosa, taste
bud, hair follicle or dental pulp is a candidate only with reproductive evidence;
otherwise it is Other:Tissue, a new subtype.** Its tumor value was never compared with
that tissue. 119 proteins, all formerly Cancer Antigen only: 87 olfactory receptors,
20 hair keratins and keratin-associated proteins, 6 bitter taste receptors, TAS1R2,
TAAR9, BPIFB3, BPIFB6, DSPP and PANX3. Every family rule names a home, whether or not
the filters act on its tier. 41 proteins with reproductive evidence keep their category.
[Rule](target-selection.md#home-tissue-the-panels-do-not-resolve).

## Olfactory epithelium is Manageable

Tiers now use the oncology precedent the same way for every tissue: a loss graded at
most Grade 2 by CTCAE v5 that standard cancer treatment causes without stopping for it
supports Manageable. Anosmia is Grade 1 only and is accepted with chemotherapy and
head and neck radiotherapy, so olfactory epithelium moves from Protected to
Manageable. Taste buds stay Protected because taste loss stops vismodegib in 6% of
patients. The sources are catalogued in the [tissue tiers page](tissue-tiers.md#how-a-tier-is-chosen)
and in the registry update file.

- BPIFB3 and BPIFB6, Unsafe on olfactory expression in the draft, are recorded instead
  and are Other:Tissue.
- SCGB1C1 stays Unsafe on its sequence-confirmed expression in normal nasal respiratory
  mucosa, which is Protected.

## What moved where

| Record | 1.37.0 | 1.38.0 |
|---|---|---|
| Sweep findings | `target-normal-function-findings.csv`, read by the filters | Entries in `normal-expression-literature.json`, read by the catalog |
| Family rules | `target-curation.yaml` | `family_rules` in the same literature file |
| Tissue tiers | Proposals in `target-curation.yaml` | `normal-tissue-tier-updates.yaml`, replayed over the tissue-tier registry |
| Read-through loci | `target-readthrough-loci.csv`, a filter knob | `readthrough-loci.csv`, the catalog's sequence check |

## Other category changes

- **Twelve proteins are Unsafe on published normal expression:** ADAM30, ASCL4, BECN2,
  CCN6, DEFB106A/B, GALR3, MRGPRX1, OBP2A, SCGB1C1, SRY, SULT6B1 and TCHHL1. Each has a
  gene-specific report in a Vital or Protected tissue and a quote retrieved from the
  source.
- **Twenty read-through loci are Other:Sequence,** including ATF7-NPFF, ENSG00000273155
  and ENSG00000267618.

## Normal nose and skin

**Four candidates are recorded in normal nose or skin, from references the panels
lacked: GLYATL3 and PABPC3 (olfactory mucosa), FKBP6 (nasal ciliated cells) and EFCAB3
(skin cells).** Nasal mucosa and skin are Protected, so `strict` drops FKBP6 and
EFCAB3. RNA values are graded low, medium or
high against HPA staining with Enhanced-reliability antibodies.
[Method, table and limits](target-filters.md#normal-nose-and-skin-which-the-expression-panels-do-not-sample).

## API changes from 1.37.0

- **`TargetFilter` has seven knobs in five layers:** `max_off_target_tpm`,
  `normal_function`, `min_tumor_fold` (default 10, the catalog's rule), `membership`,
  `min_papers`, `exclude_prefixes`, `except_labels`. Removed: `max_off_target_ntpm`,
  `max_disposable`, `rare_cell`, `literature`, `mapping` (the `warnings` layer),
  `min_patients`, `min_reproductive`, `families`, `min_families`,
  `exact_members_only` and `readthrough`. `standard` no longer counts patients and
  `strict` no longer asks for five patients or exact nomination members.
- `target_selection.SINGLE_TUMOR_FOLD` and `REPEATED_TUMORS` are replaced by
  `TUMOR_FOLD`; `GERM_CELL_TUMORS` is new. Receipts name one route,
  `tumor_tenfold_above_nonreproductive_normal`, and carry `disposable_maximum`.
- `target_curation.registry()` and `target_curation.table(labels)` are replaced by
  `literature()`, `entries()`, `findings()` and `table(catalog)`.
- **`tissue_hazards` is now `tissue_tiers`.** The module, `normal-tissue-tiers.yaml`,
  `normal-tissue-tier-updates.yaml` and the fields `tier_registry_id`, `tier_policy_id`,
  `tier_assessment_id`, `requires_tier_review` and `functional_tier` replace the
  "hazard" names. There is no alias. `tissue_tiers.assess()` resolves the registry and
  its seven added tissues; `registry_assess()` answers from the base registry alone.
- Registry ids now start `oncoref:tissue-tiers:v1`. The registry's assessments are
  unchanged.
- `target_curation.nose_skin()`, `nose_skin_provenance()`, `rna_ihc_calibration()` and
  `rna_grade()` are new, with catalog columns `reference_recorded_tissue`,
  `reference_recorded_tier` and `reference_recorded_basis`.
- The catalog has new columns: `literature_recorded_tissue`, `literature_recorded_tier`,
  `literature_family`, `literature_family_tier`, `readthrough_locus`,
  `readthrough_parents`.

## Limits

- The sweep was automated. The abstract, open full text or supplementary table behind
  each of the fourteen exclusions was retrieved and quoted; the other findings carry a
  second automated reader's paraphrase.
- Several exclusions rest on one report or an abstract. The rubric counts that.
- **Olfactory receptors are not reviewed gene by gene.** The family rule names their
  home; 26 with reproductive evidence remain candidates. Whole-tissue RNA dilutes a
  receptor that one neuron in hundreds expresses, so the olfactory RNA-seq table cannot
  clear them.
- **The cancer rule rests on one sample.** It does not ask for replication, and cell
  lines or samples without a documented patient count.
- The oncology precedent for olfactory epithelium rests on abstracts, reviews and
  regulatory summaries; full texts were not read.
- The seven added tissue tiers were written from abstracts. Nasal mucosa is low
  confidence.
- **The registry and taxonomy files were renamed and reworded, not rebuilt.**
  `data-file-renames.json` lists the exact substitutions, and
  `data_file_renames.as_released()` reconstructs the released bytes and checks their
  hash. The pinned sequence bundle and dated receipts verify against that. Bound
  scripts and modules that changed are kept as released under
  `tests/fixtures/release-1.37-inputs`.
- Dated audits and earlier release notes keep the old wording, including the 1.29.0
  registry report under `docs/audits/tissue-hazards-1.29.0/`. The old page address
  `tissue-hazards/` is a pointer to the new page because published reports link to it.
- The seven added tissues are still dated updates beside the registry.
- [Nose and skin references](target-filters.md#normal-nose-and-skin-which-the-expression-panels-do-not-sample)
  record four candidates and remove two from `strict`. They hold no olfactory neurons
  at single-cell level and no hair-shaft cortex.
