# CTA screening rules

The [current report](audits/cta-screen-20261023/index.md) applies a simple ordered
screen to the pinned human HPA v23 candidate universe. Its categories are
**Retained**, **Variable**, **Excluded** and **Unassessed**. Variable is excluded
from selection, with a concrete reason to reconsider after new evidence.

## Rules

1. Require reproductive expression, usable anatomy/identity/coverage and no
   independent curation or safety exclusion.
2. Positive IHC (including Low) in a Vital, Protected or Manageable compartment
   excludes. Reproductive or Disposable staining does not trigger this veto.
3. Relax RNA limits only when local normal IHC is reliably negative in all routine
   measured cell contexts **and** the same curated gene profile has reliable
   positive reproductive staining. Otherwise use baseline limits.
4. Require every normal tissue to meet the absolute ceiling, reproductive
   selectivity and strict reproductive-higher comparison. The previous numerical
   ladder is tightened by five percentage points; absolute ceilings and the
   fixed 0.1 nTPM allowance remain. The report displays the entire matrix.
5. Apply those **same** inclusion rules to each reviewed individual Vital-tissue
   RNA donor. A failing donor excludes even if the tissue mean passes.
6. Label that exclusion **Variable** when a different donor in the same cohort
   reports zero RNA and the affected tissue has complete reliable negative IHC,
   with no competing protein/safety exclusion or unresolved quality flag.

There is **no separate fixed high-donor RNA threshold**. Zero is the source's
reported estimate, not proof of absence. Missing data are never zero. Cohorts,
source releases, genes, tissues and units are kept separate. Technical duplicate
specimens cannot manufacture two independent donors.

## Evidence scope

The positive reproductive control currently comes from HPA's curated gene
profile. These exports do not establish the same antibody, experimental run or
donor for every comparison. The gate is an interpretable research convention;
it does not calibrate antibody sensitivity. Enhanced validation based partly on
RNA is not an independent RNA replication. Primary staining and curated
interpretation must be distinguished when presumed off-target signal was
disregarded. See [HPA's methods](https://v23.proteinatlas.org/about/assays%2Bannotation)
and [application-specific antibody validation](https://pubmed.ncbi.nlm.nih.gov/27595404/).

Currently **only CTAG2 heart has imported individual-donor RNA**. Its four
reviewed v23 values (13.5, 0, 0, 7.1 nTPM) are a structured extraction of the
existing source-cited donor review. Current image metadata confirm that the
negative heart IHC donors are different individuals. These live IHC records
are supplementary provenance, not replacements for the pinned v23 annotations.
All other genes' donor profiles remain **Not assessed**, including retained genes.

To reconsider Variable: determine the source of donor variation, locate the
expressing cells and relevant isoforms, obtain matched RNA/protein evidence in
RNA-positive and RNA-zero specimens, and validate any proposed subgroup before
an explicit new review. An RNA-zero donor is not automatically an eligible
patient. Variable is not a clinical diagnosis or proof of treatment toxicity.

## Python

```python
from oncoref import cta_screen

decisions = cta_screen.results()
selected = cta_screen.retained_gene_names()
later = cta_screen.variable_gene_names()  # CTAG2; excluded from selected
selected_ids = cta_screen.retained_gene_ids()  # all 740 loci, including unnamed genes
later_ids = cta_screen.variable_gene_ids()    # disjoint from selected_ids
donors = cta_screen.reviewed_donor_rna()

# Evaluate source-shaped candidate/audit frames; all other exclusion gates remain.
decisions, tissues, donor_findings = cta_screen.screen_candidates(
    candidates, original_hpa_v23_tissue_audit, donors=donors
)
```

`Retained` means a point-screen pass on the reviewed evidence. No fitted donor
model, antibody false-negative probability, peptide-HLA absence or clinical
safety probability is implied. The [earlier uncertainty interface](cta-continuous-screen.md)
remains available for independently fitted and validated model outputs.
Symbol accessors omit loci without a current symbol; use the ID accessors to
retain every locus without substituting a paralog.

The published core/extended CTA APIs retain their original policies and already
exclude CTAG2. This new explicit selection also excludes CTAG2. Historical
reports and sequence exports retain their dated provenance; this report does
not generate new k-mer exclusion backgrounds. The tissue taxonomy supports other
species, but human HPA nTPM limits are not transferred to dog/mouse measurements.
