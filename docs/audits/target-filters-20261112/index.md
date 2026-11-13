# Target filters: categories, knobs and where proteins land

**The catalog assigns five categories and they are not changed here. Filters only
remove proteins from the 902 Target candidates, through six layers that
follow the catalog's own shape.** Filter policy `oncoref.target-filters.2026-11-12.4`, applied to
catalog `oncoref.target-selection.2026-11-10.1`. Every number below is computed by
`scripts/build_target_filter_audit.py`. [Rules reference](../../target-filters.md).

These are outputs of screening rules. They are not validated tumor antigens.

## 1. The categories

| Category | Proteins | Rule |
|---|---|---|
| CTA | 626 | Passes the shared screen; cancer evidence and reproductive evidence both qualify |
| Cancer Antigen only | 245 | Passes the shared screen; only cancer evidence qualifies |
| Reproductive Antigen only | 31 | Passes the shared screen; only reproductive evidence qualifies |
| Other | 750 | Not excluded, but fails coverage, identity, expression bounds, sequence or both evidence routes |
| Unsafe | 18,241 | An accepted normal-tissue or clinical exclusion |

Target candidates are the first three rows: 902 proteins. A filter never moves
a protein into or out of Other or Unsafe.

## 2. The six filter layers

A protein is kept only if it passes all six. `filter_failed` names the first one it
fails. `max_*` knobs mean **below** the value, `min_*` knobs mean **at least** the
value, and a missing measurement never passes.

| Layer | Question | Knobs | Removed by `strict`, alone |
|---|---|---|---|
| 1 `off_target` | How low must expression be outside tumors and reproductive tissue? | `max_off_target_tpm`, `max_off_target_ntpm`, `max_disposable` | 63 |
| 2 `warnings` | Should flags the catalog records, but does not exclude on, remove a protein? | `rare_cell`, `literature`, `mapping` | 0 |
| 3 `curation` | Does cited literature place the protein in a normal tissue, or is it a read-through? | `normal_function`, `readthrough` | 218 |
| 4 `sets` | How strong must each set's evidence be, and which set do you want? | `min_patients`, `min_tumor_fold`, `min_reproductive`, `membership` | 452 |
| 5 `papers` | How much published support? | `min_papers`, `families`, `min_families`, `exact_members_only` | 691 |
| 6 `exclusions` | Which gene families do you rule out by name? | `exclude_prefixes`, `except_labels` | 17 |

**Layer 4 works on the two sets separately.** `min_patients` and `min_tumor_fold`
raise the bar for being a Cancer Antigen. `min_reproductive` raises it for being a
Reproductive Antigen. A protein that loses one membership keeps the other, so a CTA
with weak tumor evidence becomes Reproductive Antigen only and stays in the list.
`membership` then picks the set: `membership="cta"` requires both.

## 3. Each knob alone

Target candidates kept with only that knob changed. Categories are the ones after
layer 4 has tightened the evidence.

| Knob | Value | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|---|
| `max_off_target_tpm` | 5.0 | 866 | 603 | 235 | 28 |
| `max_off_target_tpm` | 2.0 | 692 | 482 | 194 | 16 |
| `max_off_target_tpm` | 1.0 | 501 | 339 | 150 | 12 |
| `max_off_target_tpm` | 0.5 | 349 | 229 | 110 | 10 |
| `max_off_target_tpm` | 0.1 | 135 | 78 | 50 | 7 |
| `max_off_target_ntpm` | 5.0 | 877 | 609 | 238 | 30 |
| `max_off_target_ntpm` | 2.0 | 696 | 475 | 202 | 19 |
| `max_off_target_ntpm` | 1.0 | 522 | 344 | 165 | 13 |
| `max_off_target_ntpm` | 0.5 | 377 | 239 | 127 | 11 |
| `max_off_target_ntpm` | 0.1 | 155 | 89 | 59 | 7 |
| `max_disposable` | 10.0 | 901 | 626 | 244 | 31 |
| `max_disposable` | 1.0 | 889 | 619 | 241 | 29 |
| `rare_cell` | no_dismissed | 834 | 568 | 237 | 29 |
| `literature` | no_conflict | 870 | 597 | 242 | 31 |
| `mapping` | no_risk | 452 | 322 | 118 | 12 |
| `normal_function` | no_vital | 890 | 614 | 245 | 31 |
| `normal_function` | no_protected | 705 | 580 | 96 | 29 |
| `normal_function` | no_manageable | 700 | 578 | 93 | 29 |
| `readthrough` | exclude | 881 | 615 | 239 | 27 |
| `min_patients` | 2 | 882 | 614 | 225 | 43 |
| `min_patients` | 5 | 813 | 560 | 156 | 97 |
| `min_patients` | 10 | 765 | 509 | 108 | 148 |
| `min_patients` | 50 | 692 | 305 | 35 | 352 |
| `min_tumor_fold` | 10.0 | 821 | 486 | 164 | 171 |
| `min_tumor_fold` | 30.0 | 720 | 286 | 63 | 371 |
| `min_tumor_fold` | 100.0 | 673 | 146 | 16 | 511 |
| `min_reproductive` | 5.0 | 895 | 565 | 306 | 24 |
| `min_reproductive` | 10.0 | 891 | 530 | 341 | 20 |
| `membership` | cancer_antigen | 871 | 626 | 245 | 0 |
| `membership` | reproductive_antigen | 657 | 626 | 0 | 31 |
| `membership` | cta | 626 | 626 | 0 | 0 |
| `min_papers` | 1 | 473 | 442 | 25 | 6 |
| `min_papers` | 2 | 309 | 302 | 3 | 4 |
| `min_papers` | 3 | 211 | 209 | 0 | 2 |
| `min_papers` | 4 | 74 | 74 | 0 | 0 |
| `min_families` | 2 | 19 | 19 | 0 | 0 |
| `exact_members_only` | True | 896 | 621 | 245 | 30 |

## 4. Filter modes

`all` → `loose` → `standard` → `strict`: each mode keeps a subset of the one before.
Every mode keeps the positive controls PRAME, NY-ESO-1, MAGEA4; a test fails
if a change to a mode drops one.

| Filter mode | Kept | CTA | Cancer Antigen only | Reproductive Antigen only |
|---|---|---|---|---|
| `all` | 902 | 626 | 245 | 31 |
| `loose` | 473 | 442 | 25 | 6 |
| `standard` | 297 | 290 | 3 | 4 |
| `strict` | 145 | 145 | 0 | 0 |

### `loose`

- Curation: not a read-through of other protein-coding genes.
- Sets: keep proteins that are then in either antigen set.
- Papers: at least 1 distinct paper(s) in cta/oncofetal/tumor.

Target candidates **902** → off_target **902** → warnings **902** → curation **881** → sets **881** → papers **473** → exclusions **473**

### `standard`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable tissue included, is below 5 TPM.
- Curation: no cited normal function or expression in a Vital tissue.
- Curation: not a read-through of other protein-coding genes.
- Sets: a Cancer Antigen needs at least 2 independent tumors above its normal reference.
- Sets: keep proteins that are then in either antigen set.
- Papers: at least 2 distinct paper(s) in cta/oncofetal/tumor.

Target candidates **902** → off_target **866** → warnings **866** → curation **834** → sets **820** → papers **297** → exclusions **297**

### `strict`

- Off target: the highest tissue value outside tumors and reproductive tissue, Disposable tissue included, is below 4 TPM.
- Curation: no cited normal function or expression in a Vital or Protected tissue.
- Curation: not a read-through of other protein-coding genes.
- Sets: a Cancer Antigen needs at least 5 independent tumors above its normal reference.
- Sets: a Cancer Antigen's highest tumor is at least 10 x max(normal reference, 1).
- Sets: keep proteins that are then in cta.
- Papers: at least 3 distinct paper(s) in cta/oncofetal/tumor.
- Papers: no nomination that names only some of the protein's member genes.
- Exclusions: drop labels starting with MAGE, except MAGEA4.

Target candidates **902** → off_target **839** → warnings **839** → curation **628** → sets **394** → papers **155** → exclusions **145**

## 5. Where the familiar CTAs land

46 well-known names. In the catalog, 27 are Target
candidates, 18 are Unsafe and
1 Other. A mode column shows ✓ or the first layer
that removes the protein. **16 pass every mode:** MAGEA4, NY-ESO-1, CTAG2, PRAME, SSX1, SSX2, GAGE1, GAGE2A, PAGE1, SYCP1, AKAP4, ADAM2, TPTE, CT45A1, CT83, NUTM1.

| Name | Catalog | Off-target TPM | in | Off-target nTPM | Tumors above normal | CTA papers | Flags | `loose` | `standard` | `strict` |
|---|---|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 0.0 |  | 0.0 | 1394 | 7 |  | ✓ | ✓ | exclusions |
| MAGEA2 | CTA | 0.6 | white matter | 0.8 | 1677 | 4 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA4 | CTA | 0.7 | lymph node | 0.6 | 1711 | 8 |  | ✓ | ✓ | ✓ |
| MAGEA6 | CTA | 0.2 | adipose tissue | 0.2 | 2150 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEA10 | CTA | 2.2 | white matter | 2.4 | 890 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB1 | CTA | 0.0 |  | 0.0 | 552 | 6 | mapping risk | ✓ | ✓ | exclusions |
| MAGEB4 | CTA | 0.0 |  | 0.0 | 48 | 1 |  | ✓ | papers | papers |
| MAGEC1 | CTA | 0.5 | white matter | 0.6 | 890 | 5 |  | ✓ | ✓ | exclusions |
| NY-ESO-1 | CTA | 0.3 | cerebellum | 0.5 | 976 | 5 | mapping risk | ✓ | ✓ | ✓ |
| CTAG2 | CTA | 3.9 | thymus | 24.9 | 947 | 3 | dismissed rare-cell, mapping risk, literature conflict | ✓ | ✓ | ✓ |
| PRAME | CTA | 3.0 | adrenal gland | 20.3 | 4851 | 6 | dismissed rare-cell, literature conflict | ✓ | ✓ | ✓ |
| SSX1 | CTA | 1.2 | cerebellum | 1.3 | 1164 | 5 | mapping risk | ✓ | ✓ | ✓ |
| SSX2 | CTA | 2.4 | hypothalamus | 2.8 | 677 | 6 | mapping risk | ✓ | ✓ | ✓ |
| SSX4 | CTA | 0.2 | cerebral cortex | 0.2 | 687 | 2 | mapping risk | ✓ | ✓ | papers |
| GAGE1 | CTA | 0.7 | retina | 0.5 | 859 | 6 | mapping risk | ✓ | ✓ | ✓ |
| GAGE2A | CTA | 1.1 | retina | 0.8 | 830 | 6 | mapping risk | ✓ | ✓ | ✓ |
| PAGE1 | CTA | 0.0 |  | 0.0 | 464 | 4 |  | ✓ | ✓ | ✓ |
| XAGE1A/B | CTA | 1.8 | lung | 1.6 | 1843 | 2 |  | ✓ | ✓ | papers |
| SYCP1 | CTA | 1.9 | blood | 1.8 | 103 | 4 |  | ✓ | ✓ | ✓ |
| AKAP4 | CTA | 0.1 | blood | 0.0 | 77 | 3 |  | ✓ | ✓ | ✓ |
| ODF1 | CTA | 0.6 | basal ganglia | 3.4 | 35 | 2 |  | ✓ | ✓ | papers |
| ADAM2 | CTA | 0.6 | midbrain | 0.5 | 421 | 4 |  | ✓ | ✓ | ✓ |
| TPTE | CTA | 2.0 | white matter | 2.5 | 357 | 3 | mapping risk | ✓ | ✓ | ✓ |
| CT45A1 | CTA | 0.9 | retina | 0.8 | 785 | 3 | mapping risk | ✓ | ✓ | ✓ |
| CT83 | CTA | 1.7 | salivary gland | 2.6 | 1187 | 5 |  | ✓ | ✓ | ✓ |
| CSAG1 | CTA | 1.4 | basal ganglia | 1.9 | 2164 | 1 | mapping risk, literature conflict | ✓ | papers | papers |
| NUTM1 | CTA | 0.9 | cerebellum | 1.1 | 109 | 3 |  | ✓ | ✓ | ✓ |

Not candidates in the catalog:

| Name | Catalog | Off-target TPM | Off-target nTPM | Unsafe reasons |
|---|---|---|---|---|
| MAGEA3 | Unsafe |  | 53.1 | HPA tissue clusters |
| MAGEA12 | Unsafe |  | 0.5 | Established rare normal-cell RNA |
| MAGEB2 | Unsafe |  | 118.4 | HPA tissue clusters |
| MAGEC2 | Other |  | 0.1 |  |
| LAGE3 | Unsafe |  | 136.4 | HPA normal protein | Normal bulk | HPA tissue clusters | HPA immune | Monaco immune | DICE immune | Adult cornea |
| PAGE4 | Unsafe |  | 23.3 | Normal bulk |
| TEX14 | Unsafe |  | 28.4 | HPA tissue clusters |
| TDRD1 | Unsafe |  | 10.8 | Established rare normal-cell RNA |
| HORMAD1 | Unsafe |  | 14.8 | HPA tissue clusters | DICE immune | Adult cornea | Established rare normal-cell RNA |
| MAEL | Unsafe |  | 57.0 | Normal bulk | HPA tissue clusters |
| PIWIL2 | Unsafe |  | 16.0 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |
| DDX43 | Unsafe |  | 6.6 | HPA normal protein | HPA tissue clusters | Monaco immune | DICE immune | Established rare normal-cell RNA |
| CTCFL | Unsafe |  | 15.0 | Normal bulk |
| SPA17 | Unsafe |  | 345.8 | HPA normal protein | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |
| SPAG9 | Unsafe |  | 480.1 | HPA normal protein | Normal bulk | HPA tissue clusters | HPA immune | Monaco immune | DICE immune | Adult cornea | Established rare normal-cell RNA |
| PBK | Unsafe |  | 157.7 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |
| CEP55 | Unsafe |  | 106.0 | HPA normal protein | Normal bulk | HPA tissue clusters | DICE immune | Adult cornea | Established rare normal-cell RNA |
| PLAC1 | Unsafe |  | 2.4 | Established rare normal-cell RNA |
| LY6K | Unsafe |  | 106.2 | Normal bulk | HPA tissue clusters | Adult cornea | Established rare normal-cell RNA |

## 6. The MAGE family

`strict` drops every label starting with `MAGE` except MAGEA4:
17 of 18 MAGE candidates. **This is a judgment by family name, not
a measurement.** The last column shows what `strict` would decide from each protein's
own data with the name rule switched off: 10 of the 17 it drops
would pass (MAGEA1, MAGEA10, MAGEA11, MAGEA2/B, MAGEA6, MAGEA8, MAGEA9/B, MAGEB1, MAGEB16, MAGEC1).

| MAGE | Catalog | Off-target TPM | Off-target nTPM | Tumors above normal | CTA papers | Mapping risk | `standard` | `strict` on its own data |
|---|---|---|---|---|---|---|---|---|
| MAGEA1 | CTA | 0.0 | 0.0 | 1394 | 7 |  | ✓ | ✓ |
| MAGEA10 | CTA | 2.2 | 2.4 | 890 | 6 | yes | ✓ | ✓ |
| MAGEA11 | CTA | 0.1 | 0.1 | 1253 | 6 |  | ✓ | ✓ |
| MAGEA2/B | CTA | 0.6 | 0.8 | 1677 | 4 | yes | ✓ | ✓ |
| MAGEA4 | CTA | 0.7 | 0.6 | 1711 | 8 |  | ✓ | ✓ |
| MAGEA6 | CTA | 0.2 | 0.2 | 2150 | 6 | yes | ✓ | ✓ |
| MAGEA8 | CTA | 1.0 | 1.0 | 308 | 4 |  | ✓ | ✓ |
| MAGEA9/B | CTA | 0.2 | 0.2 | 1441 | 4 |  | ✓ | ✓ |
| MAGEB1 | CTA | 0.0 | 0.0 | 552 | 6 | yes | ✓ | ✓ |
| MAGEB10 | CTA | 1.9 | 2.4 | 19 | 2 | yes | ✓ | sets |
| MAGEB16 | CTA | 0.0 | 0.0 | 145 | 4 |  | ✓ | ✓ |
| MAGEB17 | CTA | 1.9 | 2.8 | 474 | 0 |  | papers | papers |
| MAGEB18 | CTA | 0.0 | 0.0 | 16 | 2 |  | ✓ | papers |
| MAGEB3 | CTA | 0.2 | 0.3 | 47 | 3 |  | ✓ | sets |
| MAGEB4 | CTA | 0.0 | 0.0 | 48 | 1 |  | papers | papers |
| MAGEB5 | Reproductive Antigen only | 0.0 | 0.0 | 1 | 0 |  | papers | sets |
| MAGEB6B | Cancer Antigen only | 0.0 | 0.0 | 4 | 0 | yes | papers | sets |
| MAGEC1 | CTA | 0.5 | 0.6 | 890 | 5 |  | ✓ | ✓ |

The data-driven alternative is `mapping="no_risk"`, which drops proteins whose reads
could belong to a paralog or a host gene. It removes the MAGEs marked above and keeps
MAGEA4, but it also removes 450 of the 902
candidates, including most SSX, GAGE and CT45 proteins.

## 7. Disposable tissue and the "below 1 everywhere" list

- **Disposable.** 13 of 902 candidates reach 1 TPM or 1 nTPM in
  appendix, gallbladder or tonsil. The highest is
  LCE3C at 21.4 TPM.
  [List](disposable-at-least-1.csv).
- **Below 1 everywhere off target.** `max_off_target_tpm=1.0` with
  `max_off_target_ntpm=1.0` keeps 455 candidates.
  [List](off-target-below-1.csv). It is a filter, not a catalog category.

## 8. The curation layer

Layer 3 acts on cited literature. A finding acts only when a second reader confirmed
it as human, adult, normal-tissue and gene-specific. Its tissue is tiered by the
hazard registry, or by a cited proposal where the registry has no call.
[Findings and method](normal-function-sweep.md) ·
[every finding with its sources](normal-function-sweep.csv).

| Curated tier | Candidates | Tissues | Removed by `normal_function` |
|---|---|---|---|
| Vital | 12 | Adrenal gland, Brain neurons, Eye, Kidney, Lung, Megakaryocytes and platelets, Pancreas, Dorsal root and trigeminal ganglion sensory neurons, Skeletal muscle | `no_vital` and stricter |
| Protected | 185 | Cartilage, Lymph node, Monocytes and macrophages, Nasal respiratory mucosa and its glands, Olfactory epithelium and olfactory sensory neurons, Tongue and oral mucosa, Salivary gland, Skin, Stomach, Eccrine sweat gland, Taste buds and taste receptor cells | `no_protected` and stricter |
| Manageable | 5 | Adipose tissue, Lactating breast, Tooth, dental pulp and odontoblasts, Thyroid C cells | `no_manageable` |
| Disposable | 42 | Hair follicle, Cutaneous melanocytes | never |

Compartments the registry did not tier, with the tier proposed here:

| Compartment | Proposed tier | Confidence | References |
|---|---|---|---|
| Olfactory epithelium and olfactory sensory neurons | Protected | medium | PMID 15023839; PMID 25271633; PMID 29934351 |
| Taste buds and taste receptor cells | Protected | medium | PMID 33545253; PMID 29133390; PMID 21481196 |
| Tooth, dental pulp and odontoblasts | Manageable | medium | PMID 20158529; PMID 38945200; PMID 29573048 |
| Dorsal root and trigeminal ganglion sensory neurons | Vital | high | PMID 33491284; PMID 11353730; PMID 19506068 |
| Nasal respiratory mucosa and its glands | Protected | low | PMID 11777241; PMID 39126290; PMID 37449389 |
| Lacrimal gland | Protected | medium | PMID 24204119; PMID 34805617; PMID 35753665 |
| Bone, osteoblasts and osteocytes | Vital | medium | PMID 17550781; PMID 10946902 |

`readthrough="exclude"` removes 21 loci that share exons with two or more other
protein-coding genes, including ATF7-NPFF, ENSG00000273155 and ENSG00000267618:
ZNF816-ZNF321P, ENSG00000241690, ENSG00000248235, PRR5-ARHGAP8, EPPIN-WFDC6, ENSG00000249773, ENSG00000255432, ENSG00000255835, ENSG00000256349, ENSG00000256825, ENSG00000258417, ENSG00000258691, ENSG00000266202, ATF7-NPFF, ENSG00000267618, ENSG00000270299, URGCP-MRPS24, ENSG00000273155, ENSG00000274944, ENSG00000279983, ENSG00000284779.

## 9. What to keep in mind

- **The paper layer does most of the cutting.** 432 of 902 candidates
  are nominated by any CTA-list paper. `min_papers=3` keeps 211, and
  170 of those are named by all of the three largest lists
  (da Silva et al. (2017), Wang et al. (2016), Chang et al. (2019)).
- **Paper counts are not independent replication.** The nominating sources are
  published catalogs and screens, and their lists overlap. Several lists naming a gene
  shows that several screens picked it, not that separate experiments confirmed it.
  [Sources and recorded scopes](paper-families.csv).
- **No mode reads the nTPM knob or the warning knobs.** PRAME is
  3.0 TPM at tissue level (adrenal gland) but
  20.3 nTPM on one single-cell cluster, and carries a dismissed
  rare-cell signal and a conflicting literature verdict. CTAG2 is
  3.9 TPM (thymus) and
  24.9 nTPM. A mode using `max_off_target_ntpm`,
  `rare_cell="no_dismissed"` or `literature="no_conflict"` would drop PRAME, so those
  stay knobs you set yourself.
- **46 candidates are below 1 TPM at tissue level but not below 1 nTPM.**
  `max_off_target_tpm=1.0` with `max_off_target_ntpm=1.0` is the "below 1 everywhere"
  filter; it removes PRAME by construction.
