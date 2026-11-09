# OncoRef 1.36.0: cancer-only audit and uniform decisions

The [individual cancer-only audit](audits/cancer-only-20261109/report.html) examines
47 proteins with strong cancer RNA and zero recorded adult reproductive bulk RNA.
It starts from the 66 matches to the existing normal-low tumor-expression screen;
the other 19 have measurable reproductive bulk RNA. No member establishes
cancer-only protein production. All 47 have some developmental RNA signal, and
reproductive-cell coverage remains incomplete.

The reviewed calls are **3 Unsafe, 36 Other:Review and 8 Other:Expression**.
The report links every protein to normal-tissue observations, olfactory assay
results, individual supporting tumor observations, Ensembl, HPA and primary papers.
It includes the 2026 fetal olfactory atlas and keeps human, animal, endogenous and
recombinant experiments distinct. HPA's uncertain shared antibodies remain review
evidence. Missing assays retain missing values.

Two assay details change interpretation: the olfactory study's calibrated RNA-copy
cutoff differs from raw Cq detection, and two probes for one gene do not provide
independent donors. Cross-reactive receptor pairs are retained explicitly.
KRTAP21-1 has HPA MS protein-existence evidence without tissue/peptide attribution
in the XML; it is preserved separately from the retrieved UniProt annotation.

`cancer_only.records()` exposes all 47 individual reviews, while
`cancer_only.audit()` returns their definitions, counts and source links.
`cta_selection.decisions()` uses CTA, Other and Unsafe with three Other reasons:
Review, Expression and Sequence. Every detailed reason and member gene is retained.
The original release labels remain available through `results()`.

This package release publishes code and reviewed evidence. The classification
snapshot and sequence bundle remain **1.35.0**, with **607 CTA, 18,259 Unsafe and
1,027 Other**; expression data remain **5.23.27**. The explicit curation replay has
**600 CTA, 18,259 Unsafe and 1,034 Other**, including the seven previously identified
normal-expression holds. New cancer-only-audit findings refine Other review reasons
and introduce no additional CTA removal or Unsafe assignment. The next scientific
data rebuild must publish updated classifications and sequence roles together.
