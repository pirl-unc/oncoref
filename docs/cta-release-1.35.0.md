# OncoRef 1.35.0: MANE-first reference proteins

Reference selection now uses **MANE Select first, Ensembl canonical second**.
The selected CDS must be complete on the primary locus. Incomplete annotation,
conflicting designated sequences and unresolved coding evidence retain a candidate
for review without an approved default. Length and alternate-locus sequences do
not override this choice.

The catalog contains **607 CTA, 18,259 Unsafe and 1,027 Other proteins** across
19,893 protein identities. Eight proteins newly qualify under the representative
rule; GDPD4, SLCO1B3–SLCO1B7, FRG2 and FRG2B move to translation review. CDY1 and
CDY1B now have separate defaults. RBMY1F and RBMY1J retain their shared identity.

The [method report](audits/cta-transcript-selection-20261108/index.md) explains
incomplete cases and how additional evidence can support, hold, resolve or prefer
an exact reference sequence. Its [searchable table](audits/cta-transcript-selection-20261108/table.html)
covers all 20,099 classified loci, including 563 awaiting representative review.
Downloads retain transcript versions, support levels, alternatives, hashes and
every change from the previous catalog. MANE and coding annotation do not establish
tumor-specific protein production.

Normal and tumor expression are reassessed after membership changes. Protein
summaries are rebuilt for 142 cancer cohorts; gene-level results reproduce exactly.
NUTM's stored input contains negative or infinite values in a declared TPM field.
Its historical protein summaries are marked unavailable and rejected by the
reader. Source repair remains in [#604](https://github.com/pirl-unc/oncoref/issues/604)
and [#665](https://github.com/pirl-unc/oncoref/issues/665).

The [sequence bundle](audits/sequence-bundle-20261104/index.md) retains all 245,535
reference peptide records and their source occurrences, with updated default
roles and review requirements. The membrane-coordinate audit transfers coordinates
only to an exact approved sequence. Earlier 602-protein literature reports remain
dated snapshots.

New APIs: `transcript_selection.records()` and `transcript_selection.evidence()`.
`cta_selection.profiles()` distinguishes sequence qualification, representative
approval, candidate sequences, selected transcripts and reviewed findings.
