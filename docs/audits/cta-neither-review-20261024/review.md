# Review of PR #618 and the Neither classification

Reviewed 24 October 2026 by Codex, the implementation agent. This is a recorded code/data review, not an independent human or clinical approval. [PR #618](https://github.com/pirl-unc/oncoref/pull/618) merged as `afcc9019f0cd218a8bed269c1439de61d46a8f10`; its eight CI checks passed, but GitHub contains no submitted PR review. The earlier conversational claim of review referred to implementation checks.

## Findings and fixes

| Finding | Consequence | Resolution |
|---|---|---|
| All added RNA conflicts became Neither regardless of magnitude. | ACRBP, TEX101 and other high normal-expression genes were absent from the forbidden background. | An observed value >30 TPM or nTPM outside mapped Reproductive/Disposable contexts now yields Unsafe / High RNA. |
| Tissue-cluster flags required normal expression at least as high as the reproductive maximum. | A high reproductive peak could hide substantial normal expression, including MAGEA3 and ROPN1. | The new absolute ceiling is independent of enrichment and applies to previous CTA and Neither genes alike. |
| Individual HPA immune samples and the older HPA cell-type summaries were downloaded but did not enter the added screen. | Aggregate summaries could conceal high specimens. | Audit these tables, both native TPM/nTPM sample columns, the original DICE distributions and the reviewed heart donors. |
| Published summary/HTML text contained fixed counts. | Reusing the old renderer after changing rules could display stale totals. | Create a dated report with calculated membership totals, benchmark losses, native-unit maxima and current protein partitions. Preserve the previous report unchanged. |
| Candidate-specific uncertainty was dispersed across exports. | Reconsideration requirements were difficult to find. | A versioned follow-up CSV, readable per-gene record and Python API preserve questions, evidence, limitations and decisions. Notes cannot override rules. |

## Evidence decisions

- The requested >30 rule is a conservative research exclusion, not a validated clinical threshold. One high observed value suffices; low-read or suspected contamination findings remain visible for adjudication.
- TPM and nTPM are evaluated in their native units. CP10K and raw counts are rejected by this rule. Existing corneal review remains separate.
- Reproductive and Disposable expression remain permitted. Unsited/unmapped expression does not establish permission. The report identifies exclusions supported only by unsited observations.
- Both former CTA and Neither candidates are checked. Variable, antibody uncertainty and reproductive-positive IHC cannot waive high normal RNA. Existing Unsafe reasons are retained.
- Every previously non-Unsafe gene must have an explicit audit record. Missing measurements are not zero. Classification verifies agreement between high-RNA flags and the exported native-unit maxima.
- No benchmark allowlist was introduced: MAGEA3 and ROPN1 are explicit benchmark losses. HORMAD1, TDRD1, BRDT, DAZL, FOXR2 and CTAG2 remain under review with recorded next-evidence requirements.
- Shared full sequences and incomplete translations are repartitioned into the forbidden background using existing Unsafe > Neither > CTA precedence. No isoform-specific expression is inferred.

## Validation record

Regression coverage tests the exact 30 boundary, one-value exclusion, permitted contexts, unresolved anatomy, both units, invalid/missing measurements, species/stage rejection, Variable/CTA precedence, complete audit coverage, published Neither maxima and the watched-gene decisions. Published files and implementation inputs are checksum-pinned in the [run manifest](run-manifest.json). The follow-up PR records final local and CI results and merge status.

[Current report](index.md) · [Candidate follow-up](follow-up.md) · [Previous report](../cta-normal-cells-20261023/index.md)
