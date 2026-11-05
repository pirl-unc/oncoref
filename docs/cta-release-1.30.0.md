# OncoRef 1.30.0: CTA and oncofetal candidate definition

This release aligns candidate classification, evidence status and downstream
sequence use under an explicit definition. The [illustrated, searchable report](audits/cta-definition-20261104/report.html)
and [PDF](audits/cta-definition-20261104/cta-definition.pdf) explain the path from
the input universe to CTA, Unsafe and Other, with the complete gate rules,
overlapping evidence and worked decisions. The [release receipt](audits/cta-definition-20261104/release-receipt.json)
binds current counts, policies, source audits and generated outputs.

- Unresolved findings and missing required evidence require review. Unsafe
  requires accepted exclusion evidence. Other alone cannot establish sequence allowance.
- Every historical finding has an explicit migration disposition. MAGEC2,
  LGALS13 and PLCZ1 retain material holds. MAGEA12/EPS8L2 causal explanations remain
  unresolved; MAGEA12 has an independent rare-cell exclusion.
- Current normal-expression breadth uses the tissue registry. The historical
  panel remains reproducible, including its older skin, spleen and thymus tiers.
- Selected primary reference translations must qualify independently of tissue
  hazard. SERPINA2 requires translation review; MROH5 remains Unsafe. Incomplete
  records, alternate loci and unsupported residues retain their original sequences.
- The verified sequence bundle preserves every occurrence and coordinate-bearing
  AA20 segment, explicit coverage gaps and original ORFdata identity. Acquisition
  and cache repair require explicit opt-in. Historical sequence assets keep their policy.
- `ONCOREF_SHARD_CACHE_DIR` supports isolated expression-shard caches with recovery.

The [comparison with 1.29.0](audits/cta-definition-20261104/changes-from-1.29.0.csv)
separates classification changes from consumer-role changes. Earlier reports and
the 1.29.0 tissue-hazard report remain frozen. These are research-policy decisions,
not measured clinical safety probabilities or proof of peptide-HLA recognition.

Reproduction follows the [implementation roadmap](plans/cta-definition-roadmap.md)
and the [sequence build instructions](audits/sequence-bundle-20261104/index.md).
