#!/usr/bin/env python3
"""Render the target-set contract and link its measured catalog rebuild."""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyBboxPatch

from oncoref import target_selection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/audits/target-selection-20261110"


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    path = OUT / f"{name}.svg"
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=160)
    plt.close(fig)


def figures():
    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 11, "svg.hashsalt": "oncoref-target-sets-v1"}
    )
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set(xlim=(0, 11), ylim=(0, 7))
    ax.axis("off")
    ax.text(
        5.5, 6.7, "One shared screen, two antigen sets", ha="center", fontsize=19, weight="bold"
    )
    ax.text(
        5.5,
        6.25,
        "All proteins → normal screen + coverage + identity/translation + bounds + sequence",
        ha="center",
    )
    ax.add_patch(
        FancyBboxPatch(
            (0.2, 0.7),
            10.6,
            4.9,
            boxstyle="round,pad=0.02",
            facecolor="#fafafa",
            edgecolor="#8993a0",
        )
    )
    ax.text(5.5, 5.25, "Eligible proteins", ha="center", weight="bold", fontsize=13)
    ax.add_patch(Ellipse((4.1, 3.2), 5.3, 3.3, color="#3e8bba", alpha=0.26))
    ax.add_patch(Ellipse((6.9, 3.2), 5.3, 3.3, color="#d49b42", alpha=0.3))
    ax.text(3.2, 3.7, "Cancer Antigens", ha="center", fontsize=14, weight="bold")
    ax.text(3.05, 2.9, "Qualifying tumor\nevidence", ha="center", linespacing=1.5)
    ax.text(7.8, 3.7, "Reproductive Antigens", ha="center", fontsize=14, weight="bold")
    ax.text(7.95, 2.9, "Reproductive or\ndevelopmental evidence", ha="center", linespacing=1.5)
    ax.text(5.5, 3.2, "CTA", ha="center", fontsize=20, weight="bold", color="#1d3537")
    ax.text(5.5, 2.7, "intersection", ha="center", fontsize=10)
    ax.text(
        5.5,
        1.3,
        "Candidate antigen proteins = the union of both sets",
        ha="center",
        fontsize=12,
        weight="bold",
    )
    ax.text(5.5, 0.95, "Neither association qualifies → Other", ha="center", fontsize=10)
    ax.text(
        5.5,
        5.85,
        "Accepted exclusion → Unsafe     |     Failed shared admission check → Other",
        ha="center",
        color="#8b3039",
    )
    ax.text(
        5.5,
        0.2,
        "Published target support annotates every category. Areas do not represent counts.",
        ha="center",
        fontsize=10,
    )
    save(fig, "sets")

    fig, ax = plt.subplots(figsize=(10, 5.4))
    n = np.linspace(0, 10, 501)
    ax.plot(
        n,
        target_selection.TUMOR_FOLD * np.maximum(n, target_selection.DETECTION_FLOOR),
        color="#a36a20",
        lw=2.5,
        label=f"One tumor sample: ≥{target_selection.TUMOR_FOLD:g} × max(N, Disposable, 1)",
    )
    ax.set(
        xlabel="Highest eligible normal value N (same comparison scale)",
        ylabel="Tumor expression required",
        xlim=(0, 10),
        ylim=(0, 105),
    )
    ax.set_title(
        "Cancer RNA evidence: one tumor at 10 × every non-reproductive normal", loc="left", pad=15
    )
    ax.spines[["right", "top"]].set_visible(False)
    ax.grid(alpha=0.15)
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    fig.text(
        0.13,
        -0.01,
        "Any comparable TPM cohort counts. The shared normal screen applies first.",
        fontsize=10,
    )
    save(fig, "tumor-rule")


def svg(name):
    text = (OUT / f"{name}.svg").read_text()
    # Standalone figures use Matplotlib's repeated figure/font IDs. Namespace
    # them when embedding both SVGs in one HTML document.
    text = re.sub(r' id="([^"]+)"', rf' id="{name}-\1"', text)
    text = text.replace('xlink:href="#', f'xlink:href="#{name}-')
    text = text.replace("url(#", f"url(#{name}-")
    return "<svg" + text.split("<svg", 1)[1]


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    definition = target_selection.definition()
    write_guide(definition)
    (OUT / "definition.json").write_text(json.dumps(definition, indent=2) + "\n")
    figures()
    rows = "".join(
        f"<tr><td>{n:g}</td>"
        f"<td>≥{target_selection.TUMOR_FOLD * max(n, target_selection.DETECTION_FLOOR):g}</td></tr>"
        for n in (0, 0.4, 1, 3, 8)
    )
    rare = definition["rare_cell_rna"]
    fold, cap = target_selection.TUMOR_FOLD, target_selection.DISPOSABLE_CAP
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cancer and Reproductive Antigens — OncoRef</title>
<style>body{{font:17px/1.6 system-ui,sans-serif;color:#23313a;background:#fff;margin:0}}
main{{max-width:1020px;padding:35px 24px;margin:auto}}h1{{font-size:36px;line-height:1.2}}
h2{{margin-top:2em}}.status{{border-left:4px solid #a36a20;background:#faf6ec;padding:15px}}
svg{{width:100%;height:auto}}figure{{margin:25px 0}}table{{width:100%;border-collapse:collapse}}
td,th{{padding:10px;text-align:left;border-bottom:1px solid #d9dfe3}}th{{background:#eff3f5}}
code{{font-size:.9em}}a{{color:#176689}}.note{{color:#53636c}}</style><main>
<h1>Cancer and Reproductive Antigens</h1>
<p><strong>CTA = Cancer Antigens ∩ Reproductive Antigens.</strong><br>
Candidate antigen proteins = Cancer Antigens ∪ Reproductive Antigens.</p>
<p class="status">Rebuilt catalog: <code>{html.escape(definition["policy_id"])}</code>.
{", ".join(f"{target_selection.metadata()['summary']['categories'][c]:,} {c}" for c in target_selection.CATEGORIES)}.
<a href="rebuild/report.html">Filterable protein table, measurements and category changes</a>.
The released 1.36.0 catalog remains a separate snapshot under its earlier policy.</p>
<figure>{svg("sets")}</figure>
<h2>One shared eligibility screen</h2>
<p>Accepted normal-tissue or applicable clinical exclusion → <strong>Unsafe</strong>.
Required coverage, identity, translation or unresolved evidence → <strong>Other:Review</strong>.
Nested, Unsafe-contained or read-through sequence → <strong>Other:Sequence</strong>.
Target maximum or Disposable cap failure → <strong>Other:Expression</strong>.
A recorded normal home that no expression panel resolves (nose, taste bud, hair follicle,
tooth), without reproductive evidence → <strong>Other:Tissue</strong>.
The remaining proteins enter the two independent association checks.</p>
<h2>Rare normal-cell RNA</h2>
<p>Each adult donor is assessed separately in each named Vital or Protected cell population.
A count matrix cannot show where a read aligned: one mis-assigned read per cell from a paralog,
pseudogene or overlapping gene looks like sparse expression. Cells with a single molecule
therefore never change a category.</p>
<table><thead><tr><th>Question</th><th>Requirement</th></tr></thead><tbody>
<tr><td>Is there a signal?</td><td>{html.escape(rare["signal"])}</td></tr>
<tr><td>Is it replicated?</td><td>{html.escape(rare["replication"])} {html.escape(rare["single_study"])}</td></tr>
<tr><td>Could the reads belong elsewhere?</td><td>{html.escape(rare["mapping_risk"])}</td></tr>
<tr><td>Does another assay agree?</td><td>{html.escape(rare["confirmation"])}</td></tr></tbody></table>
<table><thead><tr><th>Outcome</th><th>Rule</th></tr></thead>
<tbody><tr><td>Established</td><td>{html.escape(rare["established"])}</td></tr>
<tr><td>Unresolved</td><td>{html.escape(rare["unresolved"])}</td></tr>
<tr><td>Dismissed</td><td>{html.escape(rare["dismissed"])}</td></tr>
<tr><td>Below the screen</td><td>No replicated or frequent multi-molecule signal. No finding; not evidence of absence or safety.</td></tr></tbody></table>
<p>{html.escape(rare["protected"])} {html.escape(rare["literature"])} 1%, 2% and 5% detection
remain reported comparison metrics. Strong unresolved signals require review.
<a href="../normal-cell-atlas-20261110/rare-cell/">Donor-level evidence</a> ·
<a href="rebuild/literature/">Literature review</a>.</p>
<h2>Cancer evidence</h2>
<p>Let N be the highest qualified nonexempt normal measurement in a complete, compatible
comparison panel and D the highest Disposable value. One tumor sample, from any comparable
TPM cohort, needs T ≥{fold:g} × max(N, D, 1). This is the only cancer route; counts of patients
at the actual tenfold threshold are reported separately from anonymous samples. Sample quality and exact protein identity must pass.
The highest tumor, reproductive or developmental value must also exceed every measured
nonreproductive normal and Disposable value, and Disposable expression must stay at or
below {cap:g} TPM and {cap:g} nTPM, each on its own scale.</p>
<table><thead><tr><th>Non-reproductive maximum max(N, D)</th><th>One tumor sample needs</th></tr></thead>
<tbody>{rows}</tbody></table>
<figure>{svg("tumor-rule")}</figure>
<p>Do not divide tumor TPM by an unrelated normal nTPM or CP10K maximum. Record the
comparison methods, normal contexts, source statistics and sample identities.
Missing required normal data leave the comparison unassessed. All other normal panels
continue to enforce their exclusions.
<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5903355/">Wang et al., 2018</a>.</p>
<p>Reviewed reliable tumor protein is an independent annotation. Qualitative
staining cannot bypass the RNA fold. A high tumor value never cancels a normal exclusion.</p>
<h2>Reproductive and developmental evidence</h2>
<p>Adult reproductive, placental, extraembryonic, embryonic and fetal associations count.
RA requires R ≥{fold:g} × max(N, D, 1), against the same reference as CA. R is an
individual reproductive tissue sample or supported mean within one embryo and lineage
(one population: ≥20 cells, ≥5 positive cells, ≥10% positive; smaller populations:
≥5 cells, ≥50% positive, two independent embryos at the same stage/lineage and fold).
Donors remain distinct from cells.</p>
<p>Single-cell nTPM, UMI counts and qualitative protein staining retain their own evidence
roles; they cannot be pooled into a native-TPM fold. The normal reference includes all
available individual HPA tissue samples plus supported context summaries, not every healthy
donor. <a href="https://v23.proteinatlas.org/about/download">Pinned HPA sources</a>.</p>
<h2>Published target support</h2>
<p>A curated target nomination adds supporting evidence. Preserve the exact protein–paper
link and what it establishes: nomination, biology, tumor expression, protein identity,
immune recognition, or normal-tissue counterevidence. Source-wide roles do not turn
every gene in a paper into a target.</p>
<p>Core papers cover the catalog efficiently; supplementary papers retain independent
experiments, identity limits and counterevidence. Both remain available through
<code>cta_literature.attestations()</code>, including links attached to Other and Unsafe proteins.</p>
<h2>Read the decision precisely</h2>
<p>“Cancer Antigen only” and “Reproductive Antigen only” describe membership in one qualified
set. They do not establish absence everywhere else. Unassessed associations stay visible.
Candidate antigen proteins is the union's name; it does not assert clinical safety.</p>
<p><a href="../../target-selection/">Full definitions and API</a> ·
<a href="definition.json">Machine-readable contract</a></p>
</main></html>"""
    (OUT / "report.html").write_text(page)


def write_guide(definition):
    """Keep the user-facing contract and live counts beside the generated report."""
    from oncoref import target_filters

    summary = target_selection.metadata()["summary"]
    rare = definition["rare_cell_rna"]
    categories = " · ".join(
        f"**{summary['categories'][name]:,} {name}**" for name in target_selection.CATEGORIES
    )
    counts = summary["categories"]
    tiers = summary["cta_tiers"]
    profiles = target_filters.profile_summary()
    profile_rows = "\n".join(
        f"| `{r.profile}` | {r.protected_compartments} | **{r.CTA}** | {r.CA_only} | {r.RA_only} | {r.Union} |"
        for r in profiles.itertuples()
    )
    for index_path in [ROOT / "docs/index.md", ROOT / "README.md"]:
        index_path.write_text(
            re.sub(
                r"\*\*[\d,]+ CTA, [\d,]+ Cancer Antigen only, [\d,]+ Reproductive Antigen only, [\d,]+ Other and\s+[\d,]+ Unsafe\*\*",
                f"**{counts['CTA']:,} CTA, {counts['Cancer Antigen only']:,} Cancer Antigen only, "
                f"{counts['Reproductive Antigen only']:,} Reproductive Antigen only, {counts['Other']:,} Other and\n"
                f"{counts['Unsafe']:,} Unsafe**",
                index_path.read_text(),
                count=1,
            )
        )
    (ROOT / "docs/target-selection.md").write_text(f"""# Cancer and Reproductive Antigens

**CTA = CA ∩ RA. Candidate antigen proteins = CA ∪ RA.** Cancer-only proteins remain
Cancer Antigens; calling them CTAs would erase the reproductive distinction.

Policy `{definition["policy_id"]}`: {categories}.
These are research screening memberships, not demonstrated antigen presentation or safety.

[Browse every protein](audits/target-selection-20261110/rebuild/report.html) ·
[Funnel, Venn plots and counts](audits/target-selection-20261110/rebuild/summary.md) ·
[Host profiles & optional filters](target-filters.md) ·
[Machine-readable rules](audits/target-selection-20261110/definition.json)

## Three sequential questions

OncoRef classifies every human protein group (**19,893** groups from **20,099** Ensembl 115
coding loci) by asking three independent, sequential questions:

1. **Shared somatic safety:** Is the protein safe across shared human vital, protected,
   immune, and disposable organs, with complete 36-tissue coverage, approved translation,
   resolved sequence identity, and bounded Disposable expression? (**18,253** `Unsafe`;
   **779** `Other:Review` / `Other:Sequence`; **15** bounds failures → **846** shared-eligible.)
2. **Symmetric 10× compartment enrichment:** Against the shared somatic baseline
   **H₀ = max(Normal N, Disposable D, 1 TPM)**, does expression reach **≥ 10 × H₀** in
   QC-passing tumors (**CA = {summary["memberships"]["is_cancer_antigen"]}**) and/or
   reproductive-developmental compartments (**RA = {summary["memberships"]["is_reproductive_antigen"]}**)?
   Their intersection is the default broad **CTA ({counts["CTA"]})** set.
3. **Host anatomy & toxicity tolerance (optional profile slice):** Which of the **12**
   conditional reproductive, developmental, and mammary compartments exist in a specific
   patient and are *not* tolerated for toxicity? Selecting a patient profile (`male-patient`,
   `female-patient`, `--lactating`, `--pregnant`, `core`, or custom `--allow-tissue` /
   `--protect-tissue`) moves protected host compartments into the off-target denominator.

![End-to-end curation architecture](audits/target-selection-20261110/rebuild/summary/source-to-tiers-architecture.png)

## One threshold, two compartments

Let **H₀ = max(N, D, 1 TPM)**: N is the highest supported non-reproductive normal
measurement and D the highest Disposable-tissue measurement (`appendix`, `gallbladder`, `tonsil`).
Both sets require **at least 10 × H₀**, independently, after the same eligibility screen.

| Set | Required positive RNA evidence |
|---|---|
| CA | At least one QC-passing tumor sample ≥10 × H₀ (and not an unmeasured home tissue without RA) |
| RA | At least one reproductive tissue sample, or supported embryo-specific lineage mean, ≥10 × H₀ |
| CTA | Both CA and RA (**{counts["CTA"]}** proteins in default `all` mode) |
| CA∪RA | Either set (**{summary["memberships"]["is_target_candidate"]}** candidate proteins; RA-only does not imply tumor expression) |

All comparisons use the declared native-TPM panel. Qualitative protein staining,
single-cell nTPM and UMI counts cannot be divided by TPM or bypass its numerical fold.
Reproductive IHC and older association routes remain visible as supporting annotations.

The normal reference includes **individual HPA v23 tissue samples from all measured
compartments**, plus tissue/subregion means, GTEx v11 donor medians and supported immune
summaries. Reproductive observations are positive candidates; non-reproductive
observations tighten the denominator. Missing required normal coverage never becomes zero.
[Source receipt and limits](audits/target-selection-20261110/rebuild/comparison.json).

Developmental means are computed within one embryo and lineage. One well-sampled
population requires ≥20 cells, ≥5 cells with ≥1 TPM and ≥10% positive cells.
Smaller populations require ≥5 measured cells, ≥50% with ≥1 TPM, and two independent
embryos at the same age and lineage reaching 10 × H₀. Embryos are biological donors;
multiple cells or lineages from one embryo are not independent donors.

![Both enrichment comparisons](audits/target-selection-20261110/rebuild/summary/tumor-normal.png)

## CTA tiers, lactating mammary, and patient toxicity profiles

The default `CTA` set (**{counts["CTA"]}** proteins) intentionally spans **all 12**
reproductive, developmental, and mammary compartments so no sex-specific or surgical
target is discarded upstream. Those 12 compartments divide into four biological groups:

- **Core germline & fetoplacental (`4`):** `embryo`, `placenta`, `testis`, `ovary`
- **Male accessory organs (`3`):** `epididymis`, `prostate`, `seminal vesicle`
  (anatomically absent in XX females; vascularized MHC-I+ organs in intact XY males)
- **Female Müllerian tract (`3`):** `fallopian tube`, `endometrium`, `cervix`
  (anatomically absent in XY males; vascularized MHC-I+ organs in intact XX females)
- **Mammary physiological states (`2`):** `lactating_breast` vs `resting_breast`

### Why `lactating_breast` is distinct from `resting_breast`

Men possess rudimentary mammary ductal epithelium (which expresses constitutive ductal
genes and can undergo gynecomastia), so **resting mammary tissue is protected in both
sexes**. By contrast, secretory milk-protein transcription requires lactogenic hormonal
differentiation: in HPA v23, biopsy `breast.410` expresses coordinated milk proteins at
massive levels (`LALBA` 9,523 TPM, `CSN1S1` 6,058 TPM, `CSN2` 2,530 TPM, `CSN3` 1,582 TPM,
`BTN1A1` 140 TPM), whereas all resting breast biopsies (`breast.373`, `breast.390`,
`breast.405`) are near zero (`LALBA` ≤ 0.5 TPM, `BTN1A1` ≤ 0.2 TPM). Separating
`lactating_breast` (`≥ 3` of `LALBA, CSN1S1, CSN2, CSN3` at `≥ 100 TPM`) from
`resting_breast` keeps `LALBA` and `BTN1A1` available for non-lactating female and male
patients while immediately excluding them when `lactating=True`.

### Three biological CTA tiers (`cta_tier`)

Every CTA row carries an explicit `cta_tier` decomposition (`{tiers["core"]} + {tiers["core_with_extended_exposure"]} + {tiers["extended_only"]} = {counts["CTA"]}`):

- **`core` ({tiers["core"]} CTAs):** Qualifies in `embryo`, `placenta`, `testis`, or `ovary`,
  and every one of the 8 extended adult compartments is ≤ 32 TPM and ≥ 10× below peak tumor
  expression (e.g., `MAGEA4`, `NY-ESO-1`, `PRAME`, `SSX2/B`, `MAGEC2`).
- **`core_with_extended_exposure` ({tiers["core_with_extended_exposure"]} CTAs):** Qualifies
  in a Core compartment, but also exceeds 32 TPM or lacks a 10× tumor margin in an adult
  accessory, tract, or resting-breast sample (e.g., `MAGEA11`, `MAGEB4`, `RNASE11`,
  `KIAA1210`, `ADAM2`, `TCP11X2`).
- **`extended_only` ({tiers["extended_only"]} CTAs):** Qualifies only via extended compartments —
  21 male-organ antigens (`TGM4`, `RLN1`, `SEMG2`, `POTEH`, `SPAG11A/B`, `PATE1`, ...),
  2 lactation antigens (`LALBA`, `BTN1A1`), and 1 shared tract antigen (`SCGB1D4`).

### Selecting which compartments are protected in a patient

Every host profile is defined by **one set P of protected compartments** (tissues present
in the patient where toxicity is *not* accepted):

1. Every compartment in **P** must stay **≤ 32 TPM** (`DISPOSABLE_CAP`),
2. Tumor expression must reach **≥ 10 × max(H₀, max_P TPM)**, and
3. `RA` must qualify (**≥ 10 × H₀**) in at least one **allowed** compartment outside **P**.

![CTA tiers and host profiles](audits/target-selection-20261110/rebuild/summary/cta-tiers-and-profiles.png)

| Profile | Protected compartments (**P**) | CTA | CA only | RA only | CA∪RA |
|---|---|---:|---:|---:|---:|
{profile_rows}

Callers can toggle individual tissues directly (e.g., `allow_tissues=("prostate", "seminal vesicle")`
after radical prostatectomy, `allow_tissues=("breast",)` after bilateral mastectomy, or
`protect_tissues=("testis",)` / `protect_tissues=("ovary",)` for fertility/endocrine preservation).

## Shared eligibility and Other

| Layer | Consequence |
|---|---|
| Accepted normal RNA, normal protein or applicable clinical exclusion | Unsafe |
| Missing coverage, uncertain identity/translation or unresolved normal evidence | Other:Review |
| Unsupported nested protein, ≥95% Unsafe 9-mer containment or read-through locus | Other:Sequence |
| Unmeasured normal home without qualifying reproductive evidence | Other:Tissue |
| Neither enrichment compartment qualifies, or expression bounds fail | Other:Expression |

Unsafe takes precedence. Other:Review then precedes the other Other subtypes.
Every failed check remains on the row, even when another determines its category.
**Other is neither an accepted hazard nor a claim of absent normal expression.**

Disposable expression retains separate ceilings of 32 native TPM and 32 nTPM.
The target maximum must exceed all measured non-reproductive contexts. Tissue
allowances apply to resolved cell compartments: immune or vascular cells inside a
reproductive organ do not automatically inherit that organ's allowance.

## Consequences of tissue loss

Reproductive membership and the consequences of losing the expressing tissue are
separate questions. The [tissue registry](tissue-tiers.md) records **acute loss, lasting
loss, repair and management**. Whole testis and ovary retain Protected assessments
(reflecting Leydig/endocrine and germ-cell loss in bulk measurements); mammary tissue has
a conditional Manageable assessment. Bulk organ measurements do not establish which cell
population expresses a protein or whether only that population would be injured.

```python
from oncoref import antigen_sets, tissue_tiers

loss = tissue_tiers.consequences("testis")
antigen_sets.loss_consequences("LALBA")  # shows breast (lactating_breast): 9523.2 TPM
```

## Rare normal-cell RNA

Single-cell and single-nucleus data are independent normal-expression layers.
They can identify a normal population diluted in bulk tissue. The donor-resolved
UMI rubric is distinct from HPA cell-cluster nTPM and bulk TPM.

| Check | Rule |
|---|---|
| Signal | {rare["signal"]} {rare["protected"]} |
| Replication | {rare["replication"]} |
| Frequent single-study signal | {rare["single_study"]} |
| Orthogonal confirmation | {rare["confirmation"]} |
| Established | {rare["established"]} |
| Unresolved | {rare["unresolved"]} |
| Dismissed | {rare["dismissed"]} |

A reference mapping-risk flag does not show where the observed reads aligned.
DCST1's replicated retinal signal and NKX2-6's frequent single-study cardiac signal
therefore require review, rather than being dismissed as artifacts. Single molecules alone
do not establish expression; reviewed gene-specific evidence can provide independent
confirmation. [Donor records, grades and reasons](audits/normal-cell-atlas-20261110/rare-cell/index.md).

## Protein staining

Reliable, identity-resolved normal protein findings can independently mark a protein
Unsafe under the tissue-tier policy. Antibody specificity, staining reliability,
cell type and source scope remain part of each finding.

## Published normal expression

Reviewed publications are adjudicated through the shared evidence records: accepted
gene-specific normal expression in Vital or Protected tissues excludes; other
observations retain their tissue, scope and rationale.
[Literature rubric and sources](audits/target-selection-20261110/rebuild/literature.md).

## Home tissue the panels do not resolve

Olfactory epithelium, nasal mucosa, taste buds, hair follicle and dental pulp are not
resolved by the comparison panel. A recorded home in one of these tissues yields
Other:Tissue unless reproductive evidence qualifies. The home remains recorded
even when RA qualifies; a candidate label does not establish a measured margin there.

## Tumor quality and replication

Passing TGCT libraries contribute normally. The shared QC screen excludes libraries
with broad olfactory-receptor background, a coordinated late-spermatid signature,
or multiple unrelated tissue signatures.
[Every library measurement](audits/target-selection-20261110/rebuild/library-qc.csv.gz).

Catalog rows expose `samples_at_tenfold`, `independent_patients_at_tenfold`,
`anonymous_samples_at_tenfold` and `cancer_types_at_tenfold`. Unknown donors never
count as known independent patients. `independent_patient_tumor_folds` stores each
qualifying known patient's maximum fold in descending order.

## API, host profiles, and peptide subtraction

```python
from oncoref import antigen_sets, target_filters, target_selection

all_ctas = antigen_sets.proteins("CTA")                              # 290 default CTAs
core_ctas = antigen_sets.proteins("CTA", profile="core")             # 257 Core CTAs
male_ctas = antigen_sets.proteins("CTA", profile="male-patient")     # 260 CTAs (includes LALBA)
female_ctas = antigen_sets.proteins("CTA", profile="female-patient") # 285 CTAs (includes TGM4, LALBA)
lactating_safe = antigen_sets.proteins(
    "CTA", profile="female-patient", lactating=True                  # 283 CTAs (excludes LALBA, BTN1A1)
)
post_rp = antigen_sets.proteins(
    "CTA", profile="male-patient", allow_tissues=("prostate", "seminal vesicle")
)
background = antigen_sets.proteins("Unsafe|Other")
record = target_selection.results("LALBA")
```

All annotated translations can feed exact 8–15mer or optional bounded-distance
subtraction, with occurrences and sequence gaps retained.
[Sequence API](antigen-sequences.md) · [Optional filters & host profiles](target-filters.md) ·
[CLI & curation build](curation.md).

## Reproduction and review status

The [rebuild](audits/target-selection-20261110/rebuild/index.md) publishes inputs,
receipts, every assignment and both sequential and sole-reason attribution.
`release_ready: false` describes outstanding scientific adjudications, including
read-level mapping review and the hair-follicle replay (#687). Open holds stay excluded
from positive sets. Historical catalogs remain dated snapshots.
""")


if __name__ == "__main__":
    build()
