# OncoRef 1.31.0: every candidate receives a classification

All **19,863 candidate protein records** now enter the same evidence assessment.
The 259 records previously omitted at entry receive **247 Other:Translation** and
**12 Unsafe** calls. The previous 19,604 classifications are unchanged, including
all **602 CTAs**. Current totals are **602 CTA, 18,255 Unsafe and 1,006 Other**.

Translation support is now a CTA admission requirement. The existing rule is
retained: qualifying reference coding annotation or supported noncanonical ORF
evidence, followed by qualification of the exact selected reference sequence.
Coding annotation does not establish experimental protein production. A failed
translation requirement produces Other:Translation unless independent normal
exclusion evidence establishes Unsafe. Neither result silently permits direct
target use.

The [inventory audit](audits/cta-inventory-20261105/index.md) preserves the historical
omission ledger and lists every restored record, supporting measurements and
current reason. Identical selected reference sequences still count together;
member loci remain traceable. These are candidate protein counts, not an inventory
of every proteoform.

The public profile API now separates `preliminary_tumor_restriction_profile` from
`tumor_restriction_status`. Six numerical matches have independent normal-cell
exclusions and are explicitly `excluded_by_normal_evidence`. The earlier
`tumor_restricted_lead` evidence flag is no longer exported. The audit also counts
stored-zero normal RNA summaries with developmental coverage gaps, without
equating RNA non-detection with protein absence or clinical safety.

The [illustrated report](audits/cta-definition-20261104/report.html),
[PDF](audits/cta-definition-20261104/cta-definition.pdf),
[evidence walkthrough](audits/cta-source-walkthrough-20261028/index.md), and
[verified sequence bundle](sequence-bundle.md) share the current policy and
release receipt. Earlier release assets retain their original policies and counts.
