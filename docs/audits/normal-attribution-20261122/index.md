# Normal Library Peak (#712) and Rare-Cell Attribution (#711) Adjudications

Pinned to `oncoref.target-selection.2026-11-20.2` (**19,893** human protein-equivalence groups; **327 CTA, 69 Cancer Antigen only, 128 Reproductive Antigen only, 1,120 Other, 18,249 Unsafe**). Generated reproducibly by `scripts/audit_normal_attribution_711_712.py` without altering `Unsafe` precedence (`18,249` proteins) or erasing raw normal observations.

## 1. HPA v23 Individual Normal-Library Peak Adjudications (`#712`)

Issue `#712` identified eight former CTAs (`C4orf17`, `TTLL8`, `GOLGA6L2`, `CCDC87`, `RNASE13`, `GAB4`, `PSG7`, `GSTT4`) whose `v1.42.0` `10×` exclusion rested on a single HPA v23 individual biopsy (`bone marrow.97`, `lymph node.175`, `skeletal muscle.339`, `parathyroid gland.393`, `salivary gland.216`, `spleen.378`).

All six HPA libraries pass published cross-tissue carryover and tissue-signature checks (`carryover_P == False`, `carryover_E == False`, `testis_tissue_content == 0.0` in `normal-library-qc.csv`). Under the active within-tissue biopsy `P95` rule:

| Protein | `v1.42` Class (`H₀`) | Active Class (`H₀` TPM) | Next Normal Context (TPM) | Paralog 31-mer % | Adjudication Disposition |
|---|---|---|---|---:|---|
| `TTLL8` | Cancer Antigen only (`1.50`) | **CTA** (`1.35`) | tonsil (`0.80` HPA) | `0.0%` | `restored_to_cta_under_active_p95` |
| `GOLGA6L2` | Reproductive Antigen only (`1.65`) | **CTA** (`1.34`) | cerebral cortex (`0.83` HPA) | `7.6%` | `restored_to_cta_under_active_p95` |
| `GAB4` | Cancer Antigen only (`1.96`) | Cancer Antigen only (`2.85`) | peripheral nerve (`2.85` GTEx P95) | `0.0%` | `superseded_by_independent_gtex_p95` |
| `CCDC87` | Reproductive Antigen only (`2.53`) | Reproductive Antigen only (`2.44`) | thyroid (`2.26` HPA) · pituitary (`1.99` GTEx P95) | `0.0%` | `retained_multi_tissue_cross_atlas_concordant` |
| `RNASE13` | Reproductive Antigen only (`3.48`) | Reproductive Antigen only (`3.19`) | salivary mean (`2.86` HPA) · spleen (`2.48` GTEx P95) | `0.0%` | `retained_multi_tissue_cross_atlas_concordant` |
| `PSG7` | Reproductive Antigen only (`6.79`) | Reproductive Antigen only (`6.12`) | skin (`3.93` HPA) · smooth muscle (`2.33` HPA) | `31.6%` | `retained_paralog_multimapping_unverified` |
| `GSTT4` | Other:Expression (`16.76`) | Other:Expression (`13.50`) | lymph node (`1.00` HPA) | `19.5%` | `retained_paralog_multimapping_unverified` |
| `C4orf17` | Reproductive Antigen only (`2.87`) | Reproductive Antigen only (`2.58`) | skin (`1.09` HPA) · blood (`0.50` Monaco) | `0.2%` | `retained_sole_adult_tissue_panel_biopsy` |

### Genome-Wide HPA Biopsy Statistic Sensitivity

Evaluated across all **19,893** protein groups (`oncoref.hpa_normal_statistic_sensitivity()`), holding `Unsafe` (`18,249`) and unscreened-home gates invariant:

- **Active policy (`GTEx v11 donor P95` + `HPA v23 biopsy P95`):** **327 CTA** (`396 CA`, `455 RA`, `0` delta).
- **Counterfactual — omit each protein's single highest HPA biopsy tissue site:** **334 CTA** (`+7`: `AK4P3`, `C4orf17`, `CCDC87`, `GOLGA6C`, `GSTT4`, `HIGD2B`, `PSG7`).
- **Counterfactual — omit all `hpa_tissue_samples:v23` individual biopsies (keeping HPA tissue means + GTEx P95):** **335 CTA** (`+8`: adds `SPATA31A7`).

Because HPA v23 provides primary adult bone marrow and parathyroid coverage absent from GTEx v11, and no QC defect invalidates the six libraries, all raw observations remain preserved.

---

## 2. Rare-Cell Normal Attribution Adjudications (`#711`)

Issue `#711` prioritized locus-level attribution across the **15** released tenfold-passing candidates held on rare normal-cell UMI evidence (`oncoref.rare_cell_candidate_adjudications()`):

1. **Host-overlap confounder rejected on multi-study signals → promoted to `Unsafe` (`2` proteins):**
   - `RGSL1` (`4` brain studies, `26` donors; upstream `LINC00272` terminates `35.3 kb` upstream of selected `ENST00000294854.13`) and `OOSP4B` (`2` brain studies, `4` donors; upstream `ENSG00000295845` terminates `20.1 kb` upstream of `ENST00000642343.2`).
   - Both are verified present in `antigen_sets.gene_ids("Unsafe")` and assigned `shared_sequence == "exclude"` in `sequence_bundle.target_selection_dispositions()`.
2. **Spliced exon-junction / 3'-polyA BAM attribution required (`4` `Other:Review` proteins):**
   - `CAGE1` (3' exons inside ubiquitous `SSR1`), `RAD21L1` (`ENST00000402452` 3' exons shared with neuronal/retinal `SNPH`), `PPP3R2` (intronless retrogene nested inside neuronal `GRIN3A`), and `TAS2R1` (single-exon CDS nested inside `ENSG00000248525` / `LINC02112`).
3. **Multimapper-excluded unique-locus UMI realignment required (`6` `Other:Review` proteins):**
   - `PSG8` (`54.8%` shared 31-mers), `ZNF280A` (`13.4%` 3'-kb shared 31-mers), `POTEC` (`75.1%`), `VCX3B` (`84.9%`), `PSG3` (`46.4%`), and `TUBB8` (`58.1%`).
4. **Independent second donor study or targeted ISH/IHC required (`3` `Other:Review` proteins):**
   - `DAZL` (host `ENSG00000229271` rejected as confounder, but Sanes eye astrocyte signal appears in `1` study / `4` donors and is unconfirmed in Baylor eye astrocytes), `DCAF8L2` (`1` eye study, `2` donors), and `FBXO47` (`1` eye study, `1` donor).

## Programmatic Access

- `oncoref.normal_library_peak_adjudications(disposition=...)`
- `oncoref.hpa_normal_statistic_sensitivity()`
- `oncoref.rare_cell_attribution_reviews()`
- `oncoref.rare_cell_candidate_adjudications(attribution_status=..., resolution_requirement=...)`
- `oncoref.unresolved_tenfold_rare_cell_entity_ids()`
