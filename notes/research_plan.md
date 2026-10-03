# Research plan: constrained stochastic methylation drift

Build a falsifiable distinction between intrinsic maintenance/restoration and survival/clone representation at externally defined functionally weighted loci. The design specification is [epigenetic_drift_model_spec_v2.md](epigenetic_drift_model_spec_v2.md); preimplementation corrections are in [model_review.md](model_review.md).

Completed initial work: obtain the three Tier 1 methylation datasets, implement compact Rust model/measurement/simulation/inference components, reproduce empirical age/longitudinal summaries, run held-out approximate baselines and a one-snapshot protection/selection counterexample, and populate [the annotated literature review](reading_list.md). See [the results and limitations](../code/results/results.md).

Remaining scientific gates, in order:

1. Reconcile the 1,131-cell processed Polycomb table and rebuild genomic summaries from raw coverage and versioned Polycomb intervals. Check switch/intercept integration sensitivity against the published solver where available.
2. External ESL, hexamer ranks and platform context are imported and checksum recorded. Extend to lineage-specific regulatory/essentiality scores and genome-wide density, checking correct assemblies and selection-cohort independence.
3. Exploratory ESL F2/F8 runs use matched context and paired endpoint change. Fit covariate-adjusted M0 and held-out F1; test matching sensitivity and additional mask definitions with biological-replicate validation. Stop the functional-constraint claim if this gate fails.
4. Evaluate M1/M2 robustness to pooling and annotation definitions. Fit latent variants and cross-validate residual rank only if the simple distributional checks require them.
5. Extend identifiability simulation to multisite longitudinal/missingness architectures, and obtain independent clone/survival data if required. EPI-Clone and strand-paired maintenance assays are candidate informative modalities.
6. Join paired RNA/physiology without leakage; define an independent lineage-specific threshold; compare weighted, unweighted and matched-mask distances on function, not only age.

No biological conclusion is licensed by numerical simulator correctness. Keep negative outcomes and gates explicit in the append-only falsification log. Rust replaces the specification's suggested Python project scaffold, while preserving the existing repository layout and avoiding extra packages/notebooks/config directories.
