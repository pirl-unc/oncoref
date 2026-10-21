# Gene and protein sets for k-mer comparisons

Each gene and each complete exact protein sequence appears in **one group only**.

| Group | Genes | Exact protein sequences |
|---|---:|---:|
| Normal | 18,064 | 126,691 |
| Unassessed | 610 | 1,089 |
| Reproductive | 1,022 | 1,782 |
| Other | 401 | 423 |

## Choose a download

| Group | Gene list | Protein list | Sequences |
|---|---|---|---|
| Normal | [Genes](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/normal-genes.csv.gz) | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/normal-proteins.csv.gz) | [FASTA](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/normal.fasta.gz) |
| Unassessed | [Genes](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/unassessed-genes.csv.gz) | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/unassessed-proteins.csv.gz) | [FASTA](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/unassessed.fasta.gz) |
| Reproductive | [Genes](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/reproductive-genes.csv.gz) | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/reproductive-proteins.csv.gz) | [FASTA](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/reproductive.fasta.gz) |
| Other | [Genes](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/other-genes.csv.gz) | [Proteins](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/other-proteins.csv.gz) | [FASTA](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/other.fasta.gz) |

## How shared proteins are assigned

All complete annotated isoforms are compared by exact amino-acid sequence. If any encoding gene is **Normal**, the shared protein is Normal. Otherwise, any **Unassessed** gene makes it Unassessed; then Reproductive takes precedence over Other. Every original gene, transcript and protein record remains in the [occurrence table](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/protein-occurrences.csv.gz). Verified alternate-contig identities use the existing annotation-aware resolver; matching names alone do not merge loci.

Gene groups follow the approved expression screen. External gene-specific normal-tissue evidence adds MAGEA12 to Normal; its original HPA expression group remains visible. These are sequence-defined protein forms, without PTM-specific distinctions. All isoforms inherit gene-level evidence; isoform-specific expression has not been established.

## Build a k-mer exclusion background

1. Enumerate target k-mers from **Reproductive**.
2. Subtract every k-mer in **Normal**, **Unassessed**, and the [additional background FASTA](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/additional-background.fasta.gz) containing incomplete or unusual translations from Normal/Unassessed genes.
3. Keep record boundaries intact. For incomplete or unusual translations, enumerate only contiguous standard amino-acid windows; never join across an ambiguous residue or a stop.
4. Other is available separately. It contains Disposable/Silent below-threshold expression, and does not supply positive reproductive evidence.

Full-protein sets are disjoint; their k-mers can still overlap. The subtraction step remains necessary. The [machine-readable recipe](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/exclusion-basis.json) lists the exact files. [Complete proteins with any observed normal expression](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/observed-normal.fasta.gz) are an optional stricter background; this can overlap Reproductive because a passing RNA cutoff does not erase a measured trace signal. For incomplete translations, the review table preserves each gene's group; the optional complete-protein file does not establish exhaustive sequence coverage.

## Coverage and provenance

Pinned human Ensembl **115 / GRCh38.p14**, all peptide FASTA records, with transcript completeness checked against start/stop codons and CDS completeness flags in the GTF. The gene categories use HPA v23. No expression cutoffs are transferred to another species.

**34,361 translation records** are incomplete or contain nonstandard/ambiguous residues. They remain [listed with reasons and gene groups](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/translations-for-review.csv.gz). Only records from Normal/Unassessed genes enter the additional background; incompleteness alone does not turn a Reproductive gene into a normal-tissue exclusion. **307 classified genes** have [no complete usable protein](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/genes-without-complete-proteins.csv.gz); that is a coverage gap, not evidence of sequence absence.

- [All gene assignments](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/genes.csv.gz) · [All protein assignments](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/proteins.csv.gz) · [Counts](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/counts.csv)
- [Gene mapping decisions](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/gene-mappings.csv.gz) · [Reference identity](https://github.com/pirl-unc/oncoref/releases/download/v1.13.0/reference.json) · [Checksums](run-manifest.json)
- [CTA curation report](../cta-curation-20261021-r2/index.md) · [Expression groups](../normal-expression-categories-20261021-r2/index.md)

Rebuild with `python scripts/normal_sequence_sets_report.py --out <empty-directory> --categories <expression-report> --gtf <Ensembl-115-GTF.gz> --proteins <Ensembl-115-peptide-FASTA.gz>`.
