# Complete candidate inventory and cancer-restriction accounting

All **19,893 candidate protein records** enter the main assessment under
`oncoref.cta-selection.2026-11-08.1`. Identical selected reference sequences count together;
member loci and original sequences remain traceable.

## Translation is an admission gate

The loci in the 259 formerly omitted records now map to 260 protein
identities: **248 Other:Translation**
and **12 Unsafe** calls. None becomes CTA. Normal exclusion takes
precedence; lack of translation support cannot erase it. Other:Translation requires
sequence review and supplies no direct-target allowance. Comparing the 19,604
previously classified identities with the current inventory gives **452
changed or added/removed identity rows**. [Gene-level membership and default changes](../cta-transcript-selection-20261108/index.md).

The older pre-diffuseness snapshot called 2 of these restored
records Unsafe through DICE immune evidence: P2RX5-TAX1BP3, ENSG00000258461.
They are Other:Translation under the current policy. The export preserves their
historical categories and reasons alongside the current calls; restoring records
does not reinstate superseded thresholds.

[Every restored record and all reasons](restored-records.csv.gz) ·
[Historical omission ledger](formerly-out-of-universe.csv.gz) ·
[Full current classification](../cta-selection-20261028/protein-group-classifications.csv.gz).

The retained translation rule accepts qualifying reference coding annotation or
ncORF support from mass-spectrometry peptides or at least two independent Ribo-seq
studies, followed by exact-record sequence qualification. Annotation is not
experimental proof of production. These two questions remain explicit in the API.

## What can we count as cancer-only with the available data?

The data permit a reproducible **observed-expression screen**, with coverage limits:

| Stored-zero RNA screen | Records |
|---|---:|
| Cancer expression and complete required normal bulk coverage | 19,341 |
| Also zero recorded adult RNA and reproductive/placental bulk RNA | 40 |
| Also have developmental measurements | 22 |
| Also zero recorded developmental RNA | 0 |
| Also no positive adult protein call in the measurement inputs | 0 |

**18** adult-zero records lack developmental measurements;
they remain unresolved for that comparison. [All adult-zero records](adult-zero-rna-records.csv.gz).
A zero here is a stored tissue/lineage summary, not absence from every individual
cell or state. Tiny positive values are also not automatically validated protein
production. RNA and protein claims must remain separate.

## Why an expression-profile match can still be Unsafe

The existing, less literal screen returns **66 preliminary matches**. It uses
replicated tumor expression (at least 10 TPM in two samples, or reviewed tumor
protein), low measured adult RNA, no positive adult protein call in its measurement inputs, and low measured
developmental RNA. Its nominal adult floor is 2 TPM/nTPM, with the separately
specified 0.02 CP10K corneal floor; these scales are not pooled or converted.

These are numerical matches before all curated and rare-cell evidence is applied.
In particular, the profile's protein-detection flag does not incorporate the
separate curated normal-cell knowledge ledger; those findings can establish
normal production even when the profile flag is negative.
They are not 66 established cancer-only proteins:

| Disposition after independent evidence | Records |
|---|---:|
| Normal-context coverage review required | 40 |
| Unresolved normal evidence | 15 |
| Unsafe: normal evidence | 6 |
| Translation review | 5 |

[Every match, identity, classification and review status](preliminary-expression-profile-matches.csv.gz).
The public API exposes `preliminary_tumor_restriction_profile` separately from
`tumor_restriction_status`; an Unsafe record is explicitly `excluded_by_normal_evidence`.
The former undifferentiated `tumor_restricted_lead` evidence flag is no longer exported.

The following reviewed findings demonstrate why a low routine-panel maximum
cannot establish absence from a rare normal population or an activated state:

| Protein | Counterevidence |
|---|---|
| IFNW1 | Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background. [Source](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) |
| MRGPRX4 | Direct normal neuronal protein evidence blocks unrestricted targeting. [Source](https://elifesciences.org/articles/48431) |
| CNGA2 | Established non-permitted normal-cell protein blocks unrestricted direct targeting and shared-sequence clearance. [Source](https://pubmed.ncbi.nlm.nih.gov/14618336/) |
| SCN10A | A rare normal neural population is not represented adequately by the low routine-panel bulk maximum. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC8626086/) |
| GFRAL | Localized normal neuronal protein outweighs a low bulk-panel maximum. [Source](https://pubmed.ncbi.nlm.nih.gov/28846097/) |
| IFNA4 | Established normal immune protein production blocks unrestricted targeting and contributes to the shared-sequence exclusion background. [Source](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.814627/full) |

[Counterevidence with source links and limitations](profile-counterevidence.csv).
None of these numerical screens establishes antigen presentation or therapeutic safety.

[Machine-readable summary](summary.json) · [Reproduction receipt](manifest.json).
