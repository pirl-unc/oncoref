# OncoRef 1.26.0: normal-cell evidence and a complete Other audit

**701 CTA proteins, 18,062 Unsafe, 454 Other and 646 Unassessed.** All 698 previous
CTAs remain included, including PRAME, MAGEA4, NY-ESO-1, CTAG2 and XAGE1A/B.
The 701 selected proteins represent 747 member loci.

[600-protein review and source findings](audits/cta-other-literature-20261030/index.md) ·
[PDF review](audits/cta-other-literature-20261030/other-protein-review.pdf) ·
[Full funnel and figures](audits/cta-source-walkthrough-20261028/index.md) ·
[PDF logic slides](audits/cta-source-walkthrough-20261028/cta-selection-1.26.0.pdf) ·
[Rules and API](cta-selection.md)

| Category | 1.25.0 | 1.26.0 |
|---|---:|---:|
| CTA | 698 | 701 |
| Unsafe | 18,031 | 18,062 |
| Other | 600 | 454 |
| Unassessed | 534 | 646 |

## What changed

All 600 baseline Other proteins, representing 611 genes, have a recorded judgment
and search coverage. They become **3 CTA, 31 Unsafe, 454 Other and 112 review holds**.
Search retrieval is distinct from individual paper adjudication. The report states
review depth and retains the limitations of every accepted finding.

**Normal safety now includes reviewed rare-cell and activation evidence.** Activated
immune proteins are not exempt. Examples include interferons produced by human
plasmacytoid dendritic cells, IL9 in activated T cells, IL25 in airway epithelium,
and KIR3DL3 on rare lymphocytes. GFRAL, SCN10A, MRGPRX4 and RLN3 illustrate
normal neuronal contexts poorly represented by bulk measurements. Exclusion is
an expression-based research decision, not a claim of demonstrated therapy-related injury.

**Unresolved evidence receives a hold.** Sparse olfactory-neuron RNA, uncertain
species/stage attribution and unresolved isoforms cannot silently become negative
normal evidence. OR2B6, SCGB1C2 and OR6C2 also have explicit holds for raw normal
values above 30 TPM/nTPM that lack the required read or recurrence support.

**Cancer presence is now >0 measured TPM/pTPM or reliable tumor protein.** It is
independent of safety and reproductive identity. CDY1/B (1.597 TPM) and OR2T7
(1.184 TPM) now qualify without being described as highly abundant tumor targets.
The source audit retains all 143 registered cohort views, separating incompatible
scales from native TPM. BPY2/B/C qualifies through published human germ-cell
protein evidence. Reliable reproductive protein is an independent route below
the RNA-only floor; safety and concentration still apply.

## Remaining Other categories

| Category | Proteins | Meaning |
|---|---:|---|
| Low | 210 | Below the descriptive normal-expression floor; not proof of absence |
| Normal | 124 | Supported normal expression below the applicable safety ceilings |
| Weak | 67 | Insufficient reproductive abundance or support |
| Tumor-supported | 46 | Strong cancer evidence without established reproductive or tumor-restricted identity |
| Diffuse | 6 | Adult expression is too dispersed for CTA identity |
| Disposable | 1 | Permitted Disposable expression without qualifying CTA identity |

Reproductive and Inducible remain accepted compatibility subtype queries, with
zero current rows. Established normal activation is handled by safety.

## Downstream contract

```python
from oncoref import cta_selection, gene_knowledge, reference_data

roles = cta_selection.sequence_dispositions()
findings = gene_knowledge.findings("GFRAL")
source = reference_data.ensure("olfactory_durante_supplement", "2020-02-18")
```

`direct_target_candidate` is true only for CTA. `shared_sequence` is `exclude` for
Unsafe, `review` for Unassessed, and `allow` for CTA/Other under this named screening
policy. Carry the policy, reference identity and reasons downstream. These roles
are not all-isoform coverage, peptide-presentation evidence or clinical clearance.
The independently versioned legacy forbidden FASTA/k-mer products are unchanged;
consumers must explicitly adopt the current roles rather than mixing policy versions.
Human findings and thresholds never transfer to dog or mouse by orthology.

Sequence realization still requires its own reference check. In particular,
[MROH5 and SERPINA2 (#657)](https://github.com/pirl-unc/oncoref/issues/657) have
primary-assembly peptide records in the Ensembl 112 and 115 peptide FASTAs, but
their transcripts are marked `protein_coding_LoF`; this does not establish a
complete, usable translated target. `identity_status=resolved` means the gene
identity was mapped, not that every downstream translation passes quality
control. Keep these sequence checks explicit; this release does not resolve #657.
