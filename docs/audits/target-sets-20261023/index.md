# CTA, Unsafe and Neither

**589 CTA candidates**, down from 740. **151 previous candidates move to Neither**, not the forbidden background. The rule is **at least 5 nTPM reproductive RNA when there is no positive IHC anywhere measured**. All 21 currently retained clinical/requested benchmark genes remain CTA candidates.

[Find any gene](explore.html) · [Removed candidates and reasons](removed-cta-candidates.csv) · [Gene sets](#the-three-sets) · [Protein downloads](#protein-sets-and-k-mer-exclusions) · [Threshold comparison](admission-sensitivity.csv)

## The three sets

| Set | Meaning | Target? | Forbidden background? |
|---|---|---|---|
| **CTA** | Curated candidate passing the existing RNA/IHC/selectivity rules and the new admission floor | Yes | No |
| **Unsafe** | Positive IHC in Vital/Protected/Manageable tissue, above-ceiling normal RNA, or an independent reviewed safety hold | No | Yes |
| **Neither** | Low-expression, Variable, unresolved, uncurated or selectivity-only failures without decisive normal-expression exclusion | No | **No** |

These are research labels. CTA means eligible under this screen, not proven safe treatment; Unsafe means operationally forbidden under the policy, not demonstrated clinical toxicity. Neither is not evidence of protein absence.

| Category | Genes | Exact protein sequences |
| --- | --- | --- |
| CTA | 589 | 1138 |
| Unsafe | 17405 | 125578 |
| Neither | 2103 | 3269 |

The gene universe is **20,097 canonical human protein-coding genes**, not just the 2,540 source-nominated CTA candidates. [All genes](gene-results.csv.gz) · [CTA](cta-genes.csv.gz) · [Unsafe](unsafe-genes.csv.gz) · [Neither](neither-genes.csv.gz) · [Candidate funnel](cta-candidate-results.csv.gz) · [Neither subgroups](subgroups.csv).

## What changed

- **150** previously retained genes have reproductive RNA <5 nTPM and no measured IHC at all. Their RNA is also <5 nTPM in every measured tissue.
- **IL5** adds one more removal: reproductive maximum 4.9 nTPM, measured IHC but no positive IHC anywhere. “No positive IHC” therefore removes **151**, rather than the 150 estimated for “no IHC.”
- **CTAG2 remains Neither / Variable**, with its reviewed donor RNA and negative local IHC preserved. Negative IHC and RNA donors are different people; only CTAG2 has imported individual-donor RNA evidence.
- **ZNF730 is Neither / Low**: it fails CTA selectivity and has low RNA without positive protein evidence, but does not exceed its normal RNA ceilings. It is not automatically a forbidden protein.
- Positive normal IHC, including Low, and independent safety holds still take precedence. **MAGEA12 stays Unsafe.**

| Gene | Category | Detail | Reproductive RNA (nTPM) | IHC coverage |
| --- | --- | --- | --- | --- |
| MAGEA4 | CTA | Pass | 21.5 | Positive somewhere |
| CTAG1B | CTA | Pass | 7.4 | Positive somewhere |
| PRAME | CTA | Pass | 102.2 | Positive somewhere |
| CT83 | CTA | Pass | 62.2 | Missing |
| XAGE1A | CTA | Pass | 45.7 | Missing |
| XAGE1B | CTA | Pass | 53.1 | Missing |
| MAGEA3 | CTA | Pass | 7.8 | Missing |
| CTAG1A | CTA | Pass | 8.1 | Positive somewhere |
| IL5 | Neither | Low | 4.9 | Negative where measured |
| CTAG2 | Neither | Variable | 10.1 | Positive somewhere |
| ZNF730 | Neither | Low | 4.9 | Missing |
| MAGEA12 | Unsafe | Safety | 2.7 | Missing |

The 5 nTPM value is an **admission floor**, not a replacement normal-tissue safety ceiling. A gene below 5 nTPM everywhere can still be Unsafe if, for example, brain RNA exceeds its applicable 2 nTPM baseline. Missing measurements never become zero. Low requires usable routine RNA coverage; missing/unreliable evidence remains explicitly Unassessed within Neither.

IHC counts describe curated tissue/cell-type annotations, not independent donor counts. Missing IHC and measured-negative IHC are reported separately.

## The existing tissue rules remain

Each cell is **maximum nTPM / nominal minimum selectivity**. Relaxation requires reliable negative IHC across routine local cell types plus reliable positive reproductive IHC in the curated profile.

| Tissue | Baseline | Approved negative | Supported negative | Enhanced negative |
|---|---:|---:|---:|---:|
| Vital | 2 / 80% | 4 / 75% | 6 / 70% | 12 / 65% |
| Protected | 6 / 75% | 8 / 70% | 12 / 65% | 24 / 60% |
| Manageable | 10 / 70% | 20 / 65% | 30 / 60% | 60 / 55% |

CTA selection still requires `N ≤ C`, `N ≤ R × (1−s)/s + 0.1`, and `N < R`. A selectivity failure alone removes CTA eligibility; it does **not** create an Unsafe assignment. Reviewed Variable evidence is a separate reconsideration route, excluded from both target and forbidden sets unless competing protein/safety evidence wins. [Prior ordered rules and donor evidence](../cta-screen-20261023/index.md).

## Protein sets and k-mer exclusions

Complete identical amino-acid sequences are assigned once, with every encoding gene/transcript/protein occurrence retained. Priority is **Unsafe → Neither → CTA**. Thus an Unsafe encoding gene makes the shared protein forbidden; an unresolved/Neither encoding gene blocks target selection without making that protein forbidden. XAGE1A/B remain separate loci while their identical proteins are grouped.

| Set | Protein inventory | FASTA |
|---|---|---|
| CTA | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/cta-proteins.csv.gz) | [Sequences](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/cta.fasta.gz) |
| Unsafe | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/unsafe-proteins.csv.gz) | [Sequences](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/unsafe.fasta.gz) |
| Neither | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/neither-proteins.csv.gz) | [Sequences](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/neither.fasta.gz) |

**The forbidden union is `unsafe.fasta.gz` plus `additional-forbidden.fasta.gz`. Neither, Variable, Low and Unassessed are not added.** The additional file contains incomplete/ambiguous translation records only from Unsafe genes. The historical Normal + Unassessed exclusion recipe is not the current default.

Full proteins are disjoint; k-mers can still overlap. Subtract all forbidden k-mers from CTA k-mers and never join across record boundaries, ambiguous residues or stops. **34,361 translation records** require review; **307 classified genes** have no complete usable protein. Their gaps are preserved, not interpreted as absence. Gene RNA does not establish isoform-specific translation or peptide-HLA presentation.

**Gene and protein counts differ.** For example, MAGEA3 passes at the gene level, but its complete reference protein is also encoded by unresolved alternate-contig gene `ENSG00000292186`. That shared protein stays Neither pending identity review, as it did in the previous Unassessed sequence export. Matching the gene name alone does not resolve the locus. XAGE1A/B and CTAG1A/B shared proteins remain in CTA.

[Machine-readable recipe](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/exclusion-basis.json) · [Additional forbidden records](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/additional-forbidden.fasta.gz) · [All proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/target-proteins.csv.gz) · [All occurrences](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/target-protein-occurrences.csv.gz) · [Translation review](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/translations-for-review.csv.gz) · [Missing complete proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/genes-without-complete-proteins.csv.gz) · [Unclassified sequence genes](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/unclassified-sequence-genes.csv.gz).

## Provenance and API

Pinned **human HPA v23** expression and **Ensembl 115 / GRCh38.p14** complete-sequence inventory. Existing source tissue/cell mappings and antibody confidence are preserved; human thresholds are not transferred to dogs or mice. Original dated reports remain available. [HPA sources and antibodies](../cta-curation-20261021-r2/index.md) · [Reference identity](https://github.com/pirl-unc/oncoref/releases/download/v1.16.0/reference.json) · [Input/code/output checksums](run-manifest.json).

```python
from oncoref import target_sets
genes = target_sets.results()
targets = target_sets.gene_ids("CTA")
forbidden = target_sets.gene_ids("Unsafe")
neither = target_sets.gene_ids("Neither")
```

Use `target_sets` for the current three-way sets. Older `cta_screen` results and normal-sequence exports retain their explicitly dated policies.
