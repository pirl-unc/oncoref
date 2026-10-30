# OncoRef 1.25.0: categorical protein decisions

**698 CTA proteins, 600 Other, 18,031 Unsafe and 534 Unassessed.** All 522 previous
CTAs remain included. Unsafe membership is unchanged. CTA represents 741 member
loci after identical selected reference protein sequences are collapsed.

[Full report and funnel](audits/cta-source-walkthrough-20261028/index.md) ·
[PDF slides](audits/cta-source-walkthrough-20261028/cta-selection-1.25.0.pdf) ·
[Every change and individual evidence](audits/cta-neither-20261030/update-1.25.0.md) ·
[Rules and API](cta-selection.md)

| Category | 1.24.0 | 1.25.0 |
|---|---:|---:|
| CTA | 522 | 698 |
| Other (formerly Neither) | 20 | 600 |
| Unsafe | 18,031 | 18,031 |
| Unassessed | 1,290 | 534 |

## What changed

Other now requires the complete normal screen to pass, without requiring cancer
or Reproductive identity. Every Other protein receives a named category:

| Category | Proteins | Definition |
|---|---:|---|
| Other:Diffuse | 8 | Complete HPA profile is dispersed (score >10% with a reference ≥2 nTPM), or the active concentration gate fails; no CTA route remains. |
| Other:Reproductive | 2 | Qualifying Reproductive evidence, but no qualifying cancer expression; the normal screen passes. |
| Other:Inducible | 5 | Cited expression in activated nonreproductive normal cells; a localized indirect atlas signal cannot establish Reproductive identity. |
| Other:Tumor-supported | 66 | Two distinct tumors ≥10 native TPM or reviewed reliable tumor protein, with low measured normal/developmental expression; not established tumor restriction. |
| Other:Disposable | 2 | Positive Disposable RNA/protein; measured RNA elsewhere below the descriptive 2-unit floor and no detected protein elsewhere, including low measured developmental RNA. |
| Other:Low | 254 | All measured normal RNA below the descriptive 2-unit floor and no detected normal protein. Cancer expression need not be low. |
| Other:Normal | 171 | Supported nonreproductive normal RNA reaches the descriptive 2-unit floor or reliable protein is detected; no preceding subtype or CTA definition qualifies. |
| Other:Weak | 92 | Reproductive signals are below admission requirements or lack support; no preceding subtype qualifies and no supported appreciable nonreproductive expression is established. |

The rules apply in table order. Normal describes supported normal expression below
the safety ceilings; Weak describes signals insufficient for an earlier category.
Low means low in measured normal panels, not never expressed. The 66 Tumor-supported
proteins remain reconsideration leads, with rare normal-cell coverage tracked in
[#653](https://github.com/pirl-unc/oncoref/issues/653).

Three admission repairs add 176 CTA proteins: reliable positive reproductive IHC
again counts independently; the bulk maximum uses the full mapped Reproductive
panel; and supported reproductive-cell RNA supplies a presence route. MAGEA11,
FOXR2, CAGE1, PRDM9, SAGE1 and GAGE10 illustrate the corrections. The routes overlap,
so their counts must not be added. MAGE-A4, NY-ESO-1, PRAME, CTAG2 and XAGE1A/B remain CTA.

Generic or immune cells sampled in reproductive tissue may supply indirect
association if the same-source outside panel is complete and every outside cluster
is below 5 nTPM. Published normal induction contradicts indirect-only localization:
IFNA2/14 remain Other:Inducible. A sampling site does not establish the producing cell.
All safety, normal coverage, cancer and measured concentration gates still apply.

## Downstream contract

```python
from oncoref import cta_selection

cta_selection.results(category="CTA")
cta_selection.results(category="Other", subtype="Low")
cta_selection.results()[["label", "classification", "reasons"]]
cta_selection.candidates()  # CTA only
cta_selection.audit().explain("FOXR2")
```

`classification` is CTA, Unsafe, Other:<subtype>, or Unassessed. `category` retains
CTA / Unsafe / Other / null, and `subtype` names each Other category. The deprecated
`category="Neither"` query returns Other and warns. No top-level Tumor-restricted
category is introduced. Missing normal coverage uses `incomplete_evidence`;
unresolved normal-screen holds use `normal_evidence_hold`.

Pin policy **`oncoref.cta-selection.2026-10-30.2`** and the source snapshot. Dog and
mouse retain native nominations with Unassessed status; human decisions do not
transfer through orthology. The separately versioned legacy CTA and exact-protein
background APIs retain their own policies, including k-mer exclusion semantics.

These are research screening assignments, not demonstrated clinical safety or
clearance of every epitope. The report retains source units, missing coverage,
antibody concerns and exact clinical epitope/HLA/product limitations.
