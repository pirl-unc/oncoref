# Core and extended CTA curation audit

[Both panels: figures and audits](cta-core-extended-20260929/index.md) ·
[All 25 figures, vector PDF](cta-core-extended-20260929/oncoref-cta-landscape-figures.pdf)

The default remains **core CTAs** (624 genes). **Extended reproductive CTAs**
are an explicit second panel (688 genes), adding 64 genes and removing none.
Both use the same 2,537 active coding candidates from ten selected papers and
pinned HPA v23 RNA/IHC evidence. The two bundled snapshots also preserve nine
legacy-only candidates, for 2,546 rows each. Core membership is unchanged from
current main (oncoref 1.8.206); the earlier 298/315 counts described the smaller
placental-source snapshot and are superseded by this integrated audit.

## Tissue scopes and decisions

Core RNA restriction counts testis, ovary and placenta in the numerator.
Extended restriction also counts cervix, endometrium, epididymis, fallopian tube,
prostate, seminal vesicle and vagina. Neither numerator adds breast or thymus.
Both denominators include all observed normal tissues. Deflation subtracts
1 nTPM from each tissue, flooring at zero; fractions are stored at four decimals,
as in the existing gate. Missing RNA remains missing and fails the gate.
See [RNA expression units](../cta-expression-units.md) for nTPM versus TPM and
the distinction between expression measurements and restriction fractions.

Regeneration clears prior RNA evidence for genes absent from its HPA input,
including when converting between scopes. Their fractions, nTPM measurements,
somatic detection counts and thymus detection are missing; RNA restriction axes
are `NO_DATA`, and RNA qualification/filter flags are false. An empty tissue
name or safety-flag list does not establish absence of expression. Observed
zero expression remains numeric zero and is subject to the separate expression
floor. The tissue audit likewise recomputes somatic maxima from its RNA input.

Protein evidence, its reproductive-tissue convention, reliability thresholds,
family exclusions and low-expression policy are shared. The somatic RNA maximum
excludes the permissive reproductive set and thymus; it is a different measure
from the signal outside the selected numerator. These anatomical definitions
are not claims of equal immune privilege or target safety.

The review table explicitly scopes decisions to `all`, `core` or `extended`.
PAGE4's exclusion is a **core** fraction failure: 0.4539 < 0.80. Its extended
fraction is 0.9851, so it enters the extended panel. Epididymis at 850.2 nTPM is
its largest outside-core signal, prostate is 166.6, and the somatic maximum is
smooth muscle at 23.3. Historical prostate watchlist evidence is preserved.
Its clinical evidence helper continues to describe the core default.

MAGEA11 remains excluded by its shared specificity review. TRIM64 remains
candidate-only in both panels. The existing holdback for new nominations based
only on normal reproductive expression also applies to both panels.

PATE1, PATE4 and LIPI pass the extended RNA gate but carry only historical
CTpedia tags, without membership in the imported complete-paper registry. They
remain visible as candidates under explicit extended-scope provenance reviews;
they are not silently promoted into the paper-defined panel. This does not
assert that these genes lack cancer-associated expression.

The permissive (breast-inclusive) calculation remains a sensitivity audit only,
not a third public panel.

## Reconciled funnels

| Stage | Core | Extended |
|---|---:|---:|
| source union | 3,895 | 3,895 |
| canonical ID mapped | 3,654 | 3,654 |
| protein-coding candidates | 2,537 | 2,537 |
| non-CTA removed | 2,523 | 2,523 |
| HPA restriction | 880 | 970 |
| specificity + expression | 624 | 688 |

The source union includes unresolved and noncoding identities. The selected
coding union has 883 core / 973 extended raw HPA passes, before removing three
family-excluded genes from each. Across all 2,546 snapshot rows (including the
nine archived legacy candidates), raw HPA passes are 885 / 978. These are
intermediate evidence counts, not public-panel membership.

Source lists overlap: per-source counts cannot be added. Both scopes select the
same minimum set of ten papers; each scope exports its own source-cover
certificate, full paper intake and final membership. See
[complete landscape intake](../cta-landscape-sources.md) and
[placental-source provenance](../cta-publication-sources.md).

## Extended panel additions

| Gene | Core fraction | Extended fraction | Required fraction | Largest outside-core tissue | nTPM |
|---|---:|---:|---:|---|---:|
| ADAM12 | 0.7208 | 0.8046 | 0.80 | cervix | 12.7 |
| ADAM7 | 0.0003 | 1.0000 | 0.80 | epididymis | 3638.5 |
| APOBEC4 | 0.4579 | 0.9947 | 0.97 | fallopian tube | 11.2 |
| BSPH1 | 0.0000 | 0.9834 | 0.80 | epididymis | 373.6 |
| CABYR | 0.7758 | 0.8105 | 0.80 | choroid plexus | 12.7 |
| CATSPERE | 0.9171 | 1.0000 | 0.97 | epididymis | 2.6 |
| CFAP141 | 0.3724 | 0.9712 | 0.97 | fallopian tube | 29.5 |
| CRISP1 | 0.0009 | 1.0000 | 0.80 | epididymis | 1536.4 |
| CST9L | 0.7732 | 1.0000 | 0.97 | epididymis | 25.7 |
| DEFB119 | 0.1373 | 0.9360 | 0.90 | epididymis | 1870.7 |
| DEFB123 | 0.4801 | 1.0000 | 0.97 | epididymis | 155 |
| DEFB125 | 0.0000 | 1.0000 | 0.97 | epididymis | 1.5 |
| DEFB126 | 0.0000 | 1.0000 | 0.97 | epididymis | 1251.1 |
| DEFB134 | 0.0408 | 1.0000 | 0.97 | epididymis | 24.5 |
| DNAJC5G | 0.3135 | 0.9777 | 0.95 | epididymis | 45.7 |
| DPPA3 | 0.9545 | 1.0000 | 0.97 | fallopian tube | 1.1 |
| DYDC1 | 0.9576 | 0.9961 | 0.97 | fallopian tube | 5.4 |
| EDDM3A | 0.0002 | 1.0000 | 0.95 | epididymis | 620.9 |
| EDDM3B | 0.0000 | 1.0000 | 0.80 | epididymis | 2502 |
| ELSPBP1 | 0.0000 | 0.9748 | 0.80 | epididymis | 2760 |
| ENSG00000259060 | 0.0000 | 1.0000 | 0.97 | epididymis | 3.7 |
| EPPIN | 0.2014 | 0.9953 | 0.90 | epididymis | 319.7 |
| EVX2 | 0.6333 | 1.0000 | 0.97 | seminal vesicle | 1.5 |
| FAM78B | 0.0000 | 1.0000 | 0.97 | cervix | 1.1 |
| FOXN4 | 0.9266 | 1.0000 | 0.97 | cervix | 1.8 |
| GAPDHS | 0.7938 | 0.8041 | 0.80 | skin | 13.5 |
| GPX5 | 0.0476 | 1.0000 | 0.97 | epididymis | 9 |
| GPX6 | 0.0000 | 1.0000 | 0.97 | epididymis | 3.2 |
| KIAA1210 | 0.3117 | 0.9970 | 0.80 | epididymis | 40.8 |
| MAGEB3 | 0.1931 | 1.0000 | 0.95 | epididymis | 17.3 |
| MEIG1 | 0.7739 | 0.9913 | 0.97 | fallopian tube | 3.5 |
| MSH4 | 0.7949 | 1.0000 | 0.97 | epididymis | 1.8 |
| PAGE4 | 0.4539 | 0.9851 | 0.80 | epididymis | 850.2 |
| PAPPA2 | 0.7570 | 0.8329 | 0.80 | cervix | 6.4 |
| PATE2 | 0.0000 | 1.0000 | 0.97 | epididymis | 439.7 |
| PATE3 | 0.0000 | 1.0000 | 0.97 | epididymis | 154.4 |
| PLSCR2 | 0.9606 | 0.9921 | 0.97 | epididymis | 1.4 |
| POTEG | 0.0000 | 1.0000 | 0.97 | prostate | 1.1 |
| POTEH | 0.0000 | 1.0000 | 0.97 | prostate | 1.6 |
| POTEI | 0.6111 | 1.0000 | 0.97 | prostate | 2.4 |
| POTEM | 0.1450 | 1.0000 | 0.97 | prostate | 18.1 |
| PPP1R42 | 0.6423 | 0.9919 | 0.97 | fallopian tube | 5.3 |
| PRDM9 | 0.8286 | 1.0000 | 0.97 | epididymis | 1.6 |
| RNASE11 | 0.0657 | 1.0000 | 0.97 | epididymis | 643.4 |
| RNASE12 | 0.0004 | 0.9993 | 0.80 | epididymis | 679.2 |
| RNASE13 | 0.0000 | 0.9937 | 0.97 | epididymis | 299.6 |
| SCGB1D1 | 0.0000 | 1.0000 | 0.97 | cervix | 1.2 |
| SEPTIN14 | 0.7932 | 0.8008 | 0.80 | bone marrow | 2.3 |
| SPAG11A | 0.0000 | 1.0000 | 0.90 | epididymis | 4639.4 |
| SPAG11B | 0.0000 | 1.0000 | 0.90 | epididymis | 1228.4 |
| SPATA31C1 | 0.9242 | 1.0000 | 0.97 | epididymis | 1.5 |
| SPATA31C2 | 0.2349 | 1.0000 | 0.97 | epididymis | 22.5 |
| SPATS1 | 0.6476 | 1.0000 | 0.95 | fallopian tube | 7.4 |
| SPINK2 | 0.0580 | 0.9911 | 0.80 | epididymis | 5123.9 |
| SYT10 | 0.0000 | 1.0000 | 0.97 | fallopian tube | 2.1 |
| TDRD12 | 0.8543 | 1.0000 | 0.97 | epididymis | 3.2 |
| TP53TG3B | 0.0000 | 1.0000 | 0.97 | epididymis | 4.3 |
| TTC16 | 0.7143 | 0.9951 | 0.97 | fallopian tube | 6.7 |
| UBQLNL | 0.9426 | 0.9789 | 0.97 | cervix | 1.5 |
| UMODL1 | 1.0000 | 1.0000 | 0.97 | fallopian tube | 1 |
| VRTN | 0.6000 | 1.0000 | 0.97 | epididymis | 1.4 |
| WFDC10A | 0.0000 | 1.0000 | 0.97 | epididymis | 280.5 |
| WFDC6 | 0.1364 | 0.9988 | 0.97 | epididymis | 217.3 |
| ZNF705D | 0.0000 | 1.0000 | 0.97 | epididymis | 3.2 |

The [2,546-gene tissue audit](cta-reproductive-tiers.csv) records all fractions,
HPA gate decisions, thresholds, and outside-core/somatic maxima. The comparison
CSV produced with the figures records both panels' decisions for all 2,532
family-eligible snapshot rows, including the nine archived legacy candidates.
The active source funnel contains 2,523 family-eligible candidates.

## Public access and reproduction

```python
from oncoref import cta

core = cta.cta_gene_names()
extended = cta.cta_extended_gene_names()
cta.cta_extended_gene_ids()
cta.cta_extended_df()
```

Both evidence snapshots ship in the wheel and require no downloads to read.
The regenerator recomputes both from the same candidate identity columns:

```sh
# Sidecars and delta report; add --apply to update both bundled snapshots.
python scripts/regenerate_cta_table.py --tissue-audit docs/audits/cta-reproductive-tiers.csv

# Twelve figures per scope plus a direct funnel comparison; PNG, PDF and CSV.
python -m oncoref.cli plot cta-curation --tissue-scope both --out outputs/cta-curation-both
```

The module command uses this checkout when run from the repository root.
After installing this version, the equivalent `oncoref plot ...` command works
from any directory.

Each panel exports raw and reviewed evidence, source membership/counts,
publication membership/intake counts, stage counts and RNA thresholds. Figure
thresholds use the same constants as curation: Enhanced 0.80, Supported 0.90,
Approved 0.95, Uncertain/missing 0.97. Per-source retention means final panel
membership; HPA-pass rows excluded later are shown separately.

The source importer requires an explicit retrieval date, writes sidecars by
default, preserves historical source biotypes, rejects within-source identity
collisions, and reproduces both snapshots and both provenance tables from the
original Gong workbook without dropping the other landscape sources. Regression coverage checks snapshot regeneration,
missing evidence, versioned identifiers, scope-specific reviews, complete source
intake and plotted counts/thresholds.
