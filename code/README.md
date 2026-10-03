# Rust stochastic methylation model

A single Rust crate keeps related numerical/model code in `src/lib.rs` and data analysis/reporting in `src/main.rs`. Rust 1.88+, `curl`, `unzip`, `pdftotext`, GNU `sha256sum`; no Python runtime. `Cargo.lock` pins the small dependency set. Analyses use seed 20261003; public API simulation seeds are explicit.

From the repository root:

```sh
# First checkout only: populate the Cargo cache if necessary.
cargo fetch --locked --manifest-path code/Cargo.toml
make download
make check
make reproduce
```

Data download needs network access. `make reproduce` verifies existing data checksums and runs offline. Approximately 8 GB of compressed source data and 220 MB of generated results are present; leave adequate temporary space. Source data live under ignored `data/`. Large site tables and kinetic subsets are ignored; compact results/reports are versioned. A frozen derived human feature table is rebuilt under ignored `data/derived/` from checksum-verified ESL, hexamer-ranking and GEO platform sources. CGI/regulatory class is available only for 450K probes; EPIC-only context remains unknown and is excluded from matching. Flanks containing extra CpGs may lack a rank in the 225-motif table; coverage is reported. No raw methylome expansion is needed for the current baseline run. The source matrices are streamed in full.

Commands: `annotate` (ESL/sequence/platform features), `baseline` (Polycomb), `arrays` (both longitudinal matrices), `blood`, `fit` (requires generated subsets), `simulate`, `report`, `reproduce` (all). Each accepts optional raw/results directories after the command. Examples:

```sh
cargo run --release --locked --offline --manifest-path code/Cargo.toml -- simulate
cargo run --release --locked --offline --manifest-path code/Cargo.toml -- reproduce data/raw code/results
```

## Implemented model interfaces

- M0: whole-matrix endpoint means/variance/change and donor/lineage summaries, with detection filtering and biological-replicate bootstrap. External ESL matched-mask contrasts and 500 within-stratum mask permutations run on the full human matrices. These are **unadjusted empirical summaries**, not the complete covariate-adjusted M0 in the design.
- M1/M2: exact two-state marginal probabilities and hierarchical penalized composite-likelihood fitting (`fit_kinetics`), context covariates, constraint effects on away and recovery rates, site-effect shrinkage. `kinetic_prediction` evaluates held-out observations. Rate effects are not mathematically equivalent to independently tuning diffusion restoring strength/noise.
- M3: exact competing-event selection simulation; analytic killed one-locus chain; survivor-normalized Monte Carlo likelihood (`particle_predict` / `fit_particles`). Effective removal does not simulate clone reproduction or lineage sampling.
- M4a: exact irreversible state switching plus altered CpG rates. The published Polycomb aggregate baseline instead uses continuous-switch quadrature and Gaussian Bernoulli-gain approximation, with a shared Gaussian mouse effect integrated by five-point quadrature. Within-mouse likelihood is joint, and training weights animals equally. This is an approximation rather than an exact replication of the MATLAB optimizer.
- M4b: an explicitly chosen OU log-rate state. M5: fixed low-rank loadings driven by OU rate factors. `simulate_grid` exposes time resolution; `residual_factors` supplies a centered low-rank covariance calculation. Rank/covariate adjustment and biological interpretation require independent validation.
- First passage: exact event times in jump models, grid-resolved times for continuous/factor models; death is a competing event. `Site.reference` defines youthful identity distance separately from `Site.target`, the binary away/return and selection target. Weights and site universe must be fixed; thresholds in the numerical checks are illustrative, not validated biological identity boundaries.
- Measurement: binary calls, binomial/beta-binomial bulk counts, noisy aggregate arrays, explicit missing observations. The single-cell particle interface describes **monoallelic** calls/counts with fixed call error; diploid mixtures or allele-specific/hydroxymethylation inference are not implemented. Bulk array noise is not cellular diffusion.
- Matched masks: seeded within-stratum permutations. Caller-supplied strata must encode all prespecified matching covariates. The real-data pipeline matches chromosome, CGI class, coarse regulatory class, early mean beta (bins 0.05), local probe-sequence CpG density (bins 0.01), and reverse-complement-averaged normalized hexamer rank (bins 0.1). Unsupported strata are excluded and counted. Fibroblast mask analyses use controls only. Matching on observed early methylation is exploratory and transductive; significance does not establish held-out mechanistic prediction or functional essentiality. Site-level permutations do not model all spatial dependence.

Particle inference is finite Monte Carlo with common random numbers and a weak optimization penalty; it is not Bayesian posterior sampling. Joint single-cell likelihoods can require very many particles in high dimensions. Increase particles, use independent evaluation seeds, and check stability. `fit_kinetics` balances biological replicates and shrinks site effects; its objective is a penalized composite likelihood, not a calibrated posterior. Synthetic recovery depends on pooling assumptions. The real-data M1 run uses fixed every-65,536th input probe (14 fibroblast, eight blood sites), controls only for fibroblasts, and a training-derived young reference/error scale. It is a small pipeline test, not a genome-wide mechanistic result.

## Scientific status

See [results/results.md](results/results.md), [empirical report](results/empirical_drift_report.html), and [identifiability report](results/identifiability_report.html). Checks cover analytic CTMC/selection/first-passage identities, observation normalization, synthetic constraint recovery, seeded masks, latent/factor variants and low-rank algebra. The one-state/slow-fast comparison holds out whole mice; longitudinal fits hold out donors/people. Bootstrap intervals resample the resulting replicate scores/summaries and **do not refit all training folds**, so they underrepresent full model-fitting uncertainty.

**Not completed:** raw-coverage Polycomb reconstruction, independent essentiality/identity scores, held-out F1 and sensitivity to matching definitions, confounder-adjusted M0, comprehensive multisite longitudinal protection/selection recovery, F5 biological rank selection, F6 RNA/functional identity validation, F7 clock comparisons, full Bayesian uncertainty, diploid observation modeling. Their scientific gates remain open. The code implements mechanisms needed to investigate them; the current run does not validate the complete biological hypothesis.
