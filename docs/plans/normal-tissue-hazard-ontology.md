# Normal-tissue hazard ontology: 1.29.0

Requested: a centralized, referenced set of Vital / Protected / Manageable /
Disposable calls for bulk tissues and cell subsets, with generated figures and
reports, delivered as the next oncoref PR and release.

Implementation commitments:

- One versioned registry separates anatomical membership, functional consequence,
  regeneration, replacement/support, persistent injury and evidence of recognition.
- Cover every existing adult anatomy root, resolve descendants explicitly, and add
  reviewed cell subsets across organ systems. Unknown labels retain organ concern;
  broad measurements never establish a narrower compartment restriction.
- Retain reproductive/developmental membership independently of hazard. Adult
  human policy is not transferred to mouse/dog by using an anatomical synonym.
- Human clinical, animal ablation, human functional analogy, and mechanistic
  inference remain separate. A tier is not a measured toxicity probability.
- Centralize active bulk and cell decisions, preserve historical releases, and
  report the effect on current CTA classifications before publishing.
- Produce reproducible Markdown/HTML/PDF reports, SVG/PNG figures, complete CSV
  inventories and provenance. Render the PDF for visual review.
- Validate registry coverage, aliases, precedence, cross-layer consistency,
  packaging, regeneration, local release-environment dependencies and CI.
- Open/merge the next PR and publish 1.29.0 after checks. User authorized the release.

Current baseline: main 24cca28, package 1.28.0; branch
`normal-tissue-hazard-ontology`. Repository `.venv` is the release environment;
`/Users/iskander/code/shared-virtual-env` is outside this work.

Key additional lung evidence found during this review:

- Pardo-Saganta 2015, PMC4521991: selective ciliated-cell ablation (~79%) with
  slow replacement and no compensatory progenitor proliferation; not a clinical
  tolerability experiment.
- Dumm 2019, PMC6377627: depletion of infection-surviving, predominantly
  ciliated-like airway cells worsened barrier function; mixed population and
  infection context prevent attribution to healthy mature ciliated cells alone.
- Perl 2011, PMC3056226: persistent club-compartment depletion causes lasting
  peribronchiolar fibrosis.
- Guild 2023, PMC10579698: AT1 loss can be rapidly compensated by AT2 cells;
  repeated depletion yields abnormal repair states.
- Zhao 2000, JCI9786: antigen-specific CD8 recognition of alveolar epithelium
  can initiate lethal inflammatory lung injury.

Implementation complete:

- 107 assessments (39 bulk, 68 cell), 50 references, 123 canonical human anatomy
  concepts and 1,072 native mappings; unreviewed scopes remain explicit.
- Shared RNA/IHC/rare-cell policy; Protected rare-cell review holds, independent
  reproductive membership, conservative fallback and preserved nonhuman policy.
- Complete and cached classification replays are byte-identical: 625 CTA,
  18,351 Unsafe, 628 Other; 104 changes against 1.28.0.
- Deterministic 46-page ontology PDF, searchable HTML, four SVG/PNG figures,
  CSV/JSON exports, and visually reviewed 39-slide current walkthrough.
- Registry/policy/artifact checks, historical-input checks, lint, strict docs,
  wheel/sdist builds, wheel smoke and dependency checks pass. The full local
  suite completed with 2,896 passes and 22 skips; all 12 stale policy/history
  assertions identified by that run were corrected and pass on rerun.
- Release PR: https://github.com/pirl-unc/oncoref/pull/661 . Cross-version CI
  remains the final merge/release gate.
