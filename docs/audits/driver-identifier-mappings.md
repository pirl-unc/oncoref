# Driver identifier mapping audit

PR #598 covers both bundled Bailey et al. driver tables: 739 gene-evidence rows
(299 genes) and 579 variant-evidence rows (53 genes). The original symbols,
aliases, Ensembl IDs, transcript IDs, publication locators and source-file
checksums are preserved.

The [machine-readable audit](driver-identifier-mappings.json) records hashes of
all six bundled inputs and verifies that each query returns exactly the record
IDs associated with its expected gene. Reproduce it without downloads or genome
annotations:

```bash
python scripts/audit_driver_identifiers.py --check
```

| Identifier class | Gene table | Variant table |
| --- | ---: | ---: |
| Source symbols | 299 | 53 |
| Current symbols in the bundled Ensembl 115 reference | 299 | 53 |
| Recorded source aliases | 9 | 0 |
| Ensembl IDs | 299 | 53 |
| Versioned Ensembl IDs | 299 | 53 |
| Shipped historical/alternate Ensembl IDs | 35 | 2 |
| Shipped Entrez IDs | 562 | 98 |
| Shipped NCBI synonyms | 1,409 | 337 |

Counts are per identifier class and table; they overlap. The regression suite
also covers padded/lowercase Ensembl IDs, all nine renames, missing symbols,
conflicting aliases, shared aliases and ambiguous case-insensitive queries.

## All renamed source symbols

The variant table has no renamed source symbols. The gene table contains these
nine source/current pairs, each supported by the same recorded Ensembl ID:

| Source symbol | Current symbol | Ensembl gene ID |
| --- | --- | --- |
| CBWD3 | ZNG1C | ENSG00000196873 |
| FAM46D | TENT5D | ENSG00000174016 |
| H3F3A | H3-3A | ENSG00000163041 |
| H3F3C | H3-5 | ENSG00000188375 |
| HIST1H1C | H1-2 | ENSG00000187837 |
| HIST1H1E | H1-4 | ENSG00000168298 |
| RQCD1 | CNOT9 | ENSG00000144580 |
| TCEB1 | ELOC | ENSG00000154582 |
| WHSC1 | NSD2 | ENSG00000109685 |

## CBWD3 requires source context

The pinned NCBI gene-info snapshot lists `CBWD3` as a synonym of both `ZNG1C`
(GeneID 445571, ENSG00000196873) and `ZNG1E` (GeneID 220869,
ENSG00000147996). HGNC:18519 lists `CBWD3` as a previous symbol of `ZNG1C`.
The [source evidence](driver-symbol-source-evidence.json) preserves the two
matching NCBI rows, the relevant HGNC fields, source URLs and response/file
checksums. The NCBI source is the
[pinned June 2026 gene-info snapshot](https://github.com/pirl-unc/pirlygenes/releases/download/ncbi-gene-data-20260601/Homo_sapiens.gene_info.gz);
the HGNC record is available from the
[HGNC API](https://rest.genenames.org/fetch/symbol/ZNG1C).

In the driver table, the recorded ID and the recorded `ZNG1C` alias support that
specific locus. The accessor warns on this fallback and permits queries using
the original source label. These scoped names are not promoted into the global
one-to-one synonym dictionary. This avoids assigning every occurrence of the
shared NCBI alias to a single locus.

## Lookup gaps found during review

The earlier PR validated all source rows, but only tested one renamed symbol.
It also preferred case-insensitive current symbols before exact NCBI aliases.
Consequently, `emb` (a shipped alias of `XPO1`) incorrectly selected the separate
gene `EMB` and returned no driver evidence. The reviewed implementation prefers
exact spelling; `emb` returns `XPO1`, `EMB` resolves to its separate locus, and
`Emb` raises because case folding would combine both identities. Other
ambiguous folds, including `alb` and `met`, also require exact spelling or an
Ensembl ID; exact `ALB` and `MET` remain supported.

Recorded source aliases are now queryable even when absent from the global
synonym table. A collision with another recorded or bundled identity raises.
Unknown symbols return an empty frame consistently, without falling through to
an installed historical genome and selecting an annotation-dependent hit.

This audit covers the pinned driver tables and their shipped mappings. It does
not establish completeness of all human nomenclature or refresh the global
NCBI synonym snapshot.
