# Constrained Stochastic Epigenetic Drift

Research on methylation maintenance, functional constraint, selection and identity loss.

The [Dog Aging Project falsification audit](code/results/dog_aging_falsification.md) tests an independent essential-gene proxy in longitudinal PBMC methylation. Protection is not robustly supported or rejected: estimates change sign across context/coverage choices, intervals remain wide, and the promoter comparison fails its support gate. The cohort does not identify intrinsic protection or externally calibrated clone selection; the prior experimental stopping decision remains in force.

The model specification is in [notes/epigenetic_drift_model_spec_v2.md](notes/epigenetic_drift_model_spec_v2.md). The compact Rust implementation and reproducibility commands are in [code/README.md](code/README.md), current results in [code/results/results.md](code/results/results.md), extended biological tests in [code/results/six_test_validation.md](code/results/six_test_validation.md), and the annotated literature review in [notes/reading_list.md](notes/reading_list.md). The [public-data stopping tests](code/results/broad_validation.md) find no incremental recovery prediction from independent local functional scores; known enzyme dynamics remain useful but do not establish our distinctive mechanism. The current [stopping assessment](notes/experiment_brief.md#computational-stopping-assessment-2026-10-04) defers wet experiments and records what would justify reopening the hypothesis.

## Scientific requirement

The equations are revisable; the biological objective is not. A distinctive mechanism must make a restrictive, prospectively specified biological prediction that survives testing in a new condition and cannot be rescued by condition-specific refitting. Better curve fitting, age prediction, or latent-state flexibility alone does not establish biological value. Compare against strong empirical predictors, freeze external functional annotations before outcomes, and distinguish intrinsic restoration from changing clone representation. Negative tests must narrow or retire claims rather than trigger unrestricted complexity. See [the research plan](notes/research_plan.md#biological-prediction-contract-2026-10-03) for the decision rules and current experiment direction.

Revision 3 puts external functional effects in state-dependent cell fitness, with context/exposure-driven molecular transitions and explicit clone birth/death. The [power audit](code/results/functional_power_audit.md) finds no support for the measured functional-recovery prediction, but modest effects were poorly detectable: conditional detection was 37%, 49%, 76% and 98% for 3%, 5%, 10% and 20% oracle error reductions. This does not validate selection. The [novelty assessment](notes/reading_list.md#revision-3-novelty-assessment-2026-10-04) applies two gates: academic originality and a prediction that changes a biological experiment. The broad architecture has prior art; neither a validated transfer prediction nor an actionable target panel has been earned.

## Repository Layout

- `notes/`: research plans, reading notes, project logs, and informal planning.
- `math/`: definitions, theorem targets, proof sketches, examples, and open questions.
- `latex/`: manuscript source, macros, bibliography, and TeX inputs.
- `latex/render/`: rendered PDFs and TeX build products. This directory is ignored by git.
- `code/`: optional computational checks, symbolic experiments, and figure-generation scripts.
- `sources/`: source material used during research.
- `sources/pdfs/`: local PDF references. PDF files are ignored by git.

## Suggested Workflow

1. Start by editing `notes/research_plan.md`.
2. Extract stable definitions, claims, and proof obligations into `math/`.
3. Promote mature material from `math/` into `latex/`.
4. Put references and reading notes under `sources/` and `notes/`.
5. Use `code/` only when computation helps check an example, verify algebra, generate figures, or support reproducibility.

## Theory-First Convention

The central path is:

```text
notes/research_plan.md -> math/ -> latex/
```

Code is optional and should remain subordinate to the mathematical argument unless the project explicitly becomes computational.

## Coding Guidelines

- Prefer Rust for implementation work. Use Python as a backup, with type checking.
- Prefer functional code when possible. Relax this guideline when it causes performance issues.
- Long functions are fine.
- Avoid inheritance. Object-oriented design is generally discouraged.
- Avoid excessive code fragmentation. Longer files are acceptable when they keep related logic together.
- Lint code aggressively. Use Black for Python formatting.
- Shorter code is almost always better code.
