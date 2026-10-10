# Rust stochastic methylation model

A single Rust crate keeps numerical/model code in `src/lib.rs` and related data workflows in `src/main.rs`, `src/validation.rs` and `src/broad_validation.rs`. Rust 1.88+, `curl`, `unzip`, `pdftotext`, GNU `sha256sum`; the core model needs no Python runtime. Public-data validation uses a pinned Python format decoder for MATLAB/R/XLSX files. `Cargo.lock` pins the small dependency set. Analyses use seed 20261003; public API simulation seeds are explicit.

From the repository root:

```sh
# First checkout only: populate the Cargo cache if necessary.
cargo fetch --locked --manifest-path code/Cargo.toml
make download
make check
make reproduce
```

Data download needs network access. `make reproduce` verifies existing data checksums and runs offline. Approximately 8 GB of compressed source data and 220 MB of generated results are present; leave adequate temporary space. Source data live under ignored `data/`. Large site tables and kinetic subsets are ignored; compact results/reports are versioned. A frozen derived human feature table is rebuilt under ignored `data/derived/` from checksum-verified ESL, hexamer-ranking and GEO platform sources. CGI/regulatory class is available only for 450K probes; EPIC-only context remains unknown and is excluded from matching. Flanks containing extra CpGs may lack a rank in the 225-motif table; coverage is reported. No raw methylome expansion is needed for the current baseline run. The source matrices are streamed in full.

Commands: `annotate` (ESL/sequence/platform features), `baseline` (Polycomb), `arrays` (both longitudinal matrices), `blood`, `fit` (requires generated subsets), `simulate`, `report`, `reproduce` (core). Extended commands: `forecast`, `kinetic-followup`, `growth-followup`, `matching-followup`, `followup-report`; `make followup` runs the complete extension. Each accepts optional raw/results directories after the command. Examples:

```sh
cargo run --release --locked --offline --manifest-path code/Cargo.toml -- simulate
cargo run --release --locked --offline --manifest-path code/Cargo.toml -- reproduce data/raw code/results
```

## Implemented model interfaces

The current revision-3 API is `PopulationModel` in `src/lib.rs`: independent binary molecular transitions under a fixed piecewise exposure schedule, signed state fitness, explicit birth/death, and clone labels. `population_expectation` evolves expected absolute counts with a positive second-order splitting solver; `population_decomposition` separates intrinsic methylation change from selection covariance. `simulate_population` implements exact branching event times, extinction and explicit resource-limit errors. Exact joint states are limited to 12 sites; birth faithfully inherits state, rather than modeling strand-specific replication errors. External fitness calibration, latent exposure inference and biological identity thresholds are not supplied. The earlier M0–M5 interfaces below remain historical comparators.

`make check` passes strict Clippy and 33 tests, including analytic limits, solver convergence, selection without molecular repair, clone inheritance, branching means and extinction. These are software checks, not biological validation. The primary-source novelty assessment was completed before these revised-model checks.

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

**Not completed:** raw-coverage Polycomb reconstruction, independent essentiality/identity scores, fully adjusted mechanistic F1, confounder-adjusted M0, comprehensive multisite longitudinal protection/selection recovery, F5 biological rank selection, F6 RNA/functional identity validation beyond a small proliferation pilot, F7 clock comparisons, full Bayesian uncertainty, diploid observation modeling. Their scientific gates remain open. The code implements mechanisms needed to investigate them; the current run does not validate the complete biological hypothesis.

## Continued validation

[Extended results](results/validation_followup.md) add empirical donor/chromosome-held-out stability forecasts, restricted exact-CTMC M1/M2 comparisons, full-matrix coarse/fine/density matching sensitivity and an independent measured proliferation endpoint. `make followup` regenerates all outputs after annotations are built; it streams both full matrices four times and retains ignored intermediate runs under `data/derived/matching_runs`. Paired forecast subsets use every 256th probe and are ignored. The held-out empirical log-change ridge model is distinct from mechanistic M2. Blood participants include twin pairs with unavailable family IDs, so person-level cross-validation and bootstrap are provisional. Growth tests have only 11 culture endpoints/eight donors and do not establish identity loss.

The paired fibroblast RNA metadata and processed expression matrix are downloaded as GSE179848. GSE225172 belongs to the mouse study; the fibroblast RNA matrix is now analyzed by the six-test extension; mouse RNA remains unanalyzed. Regression checks cover coefficient recovery and removal of held-out sufficient statistics.

## Six public-data validation tests

[The six-test report](results/six_test_validation.md) is the current biological assessment. One additional Rust module, `src/validation.rs`, keeps the six related data workflows together. It tests future RNA with frozen Reactome ECM promoter weights and 20 external-context-matched masks; observed post-replication counts; independent LARRY clone labels/proteins and representation decomposition; existing DNMT1 perturbations including unseen genotype combinations; full future-array distributions/residual dependence; and exact multisite protection/selection recovery with missingness and survival alternatives.

For a fresh environment, install `uv`, then create the ignored conversion environment with Python 3.10+:

```sh
uv --cache-dir data/tools/uv-cache venv --python /usr/bin/python3 data/tools/format-env
uv --cache-dir data/tools/uv-cache pip install --python data/tools/format-env/bin/python -r code/public-data-format-requirements.txt
make download
make reproduce  # builds the frozen human annotation table
make validation
```

`code/convert_public_data.py` only decodes the observed MATLAB counts and R sparse assays/metadata. Rust performs all matching, calculations, fitting, simulation and reporting. The helper is Black formatted and mypy checked. `make validation` verifies 29 source-file hashes, converts archives directly (no manual MATLAB extraction required), rebuilds ignored functional/trajectory summaries, and runs the complete six-test extension offline. The Rust CLI accepts `validation [raw] [results] [prepare|rna|kinetics|clones|perturbations|distributions|identifiability|report|all]`. New sources add roughly 550 MB compressed; derived caches and the pinned format environment remain ignored. Reactome's upstream `current` URL can change: checksum verification must reject a changed release rather than silently update the tested functional mask. API inventory JSONs are provenance aids, not model inputs or checksum-stable representations.

The RNA test uses only controls, four donors under strict initial-time matching and seven under relaxed matching. Site weights are a curated ECM program proxy; no essentiality or identity-loss threshold is inferred. Mouse/human coordinates are never joined. Clone tests use public mouse objects; human raw reads require controlled access. Future predictive intervals fail calibration despite improved CTMC means. The tested kinetic curves fail to outperform a pooled constant, and an independent-effects perturbation rule fails one genotype combination. Do not summarize these outcomes as full biological validation.

## Functional recovery pilot

`make recovery` tests externally defined BAGEL core-essential versus reference-nonessential gene-body loci in HCT116 after a methylation pulse/washout (GSE51810/51811). It requires the existing `data/derived/human_features.tsv`; run `annotate` first on a fresh checkout. Five new checksum-frozen inputs add about 81 MB compressed. Both the original matching and a later baseline-RNA-detection audit are retained as `recovery_*` and `recovery_detected_*`; tighter matching is exploratory. Only baseline and day5 DNA, baseline RNA and external context determine eligibility/matches. Negative processed RNA intensities are unavailable for log2 transformation, and later RNA missingness never excludes loci from the DNA test. The Rust routine uses titles because GEO characteristic labels are inconsistent. There is one methylation array per condition/time: matched-locus bootstrap intervals describe heterogeneity, not independent biological validation. The proposed experiment and revised mechanistic commitments are in [the experiment brief](../notes/experiment_brief.md).

## Public-data stopping tests

With the pinned format environment already installed:

```sh
make broad-download       # network; only missing frozen inputs, TAR byte ranges
make broad-validation     # checksum verification, conversion, Rust inference offline
make power                # rebuild observed design, then conditional detectability audit
make check
```

`broad_data_manifest.csv` freezes 34 primary inputs (~1 GB compressed); `broad_reference_manifest.csv` freezes the chain and 114 assembly-check sequences. Reference JSON is canonicalized to stable genome/interval/DNA fields, removing API retrieval timestamps. One out-of-bounds hg38 chr19 response is an explicit empty reference check. Changed sequence, chain or source content fails verification; the downloader never updates expected hashes. Archived RRBS members are downloaded in checked 8 MB ranges, avoiding a 6 GB TAR download. Human files use the authors' **five-column BED format: chromosome, zero-based start, end, total count, methylated count**, not standard six-column Bismark coverage. The mouse amplicon tables use fractions; the genome table uses percentages. Species and assemblies are kept separate.

All statistics and fitting live in the single related `broad_validation.rs` module. `convert_public_data.py` handles acquisition, integrity checks, XLSX decoding and intact-chain coordinate translation only. The CLI is `validation data/raw code/results broad`. Coverage/perturbation gates, independent scores and primary splits are documented before outcomes in `notes/research_plan.md`; post-result checks are explicitly exploratory. `broad_validation.md` contains full gate/penalty results and limitations. Primary regional predictions, functional scores, controls and summary tables are versioned; sensitivity forecasts and large CpG count/forecast tables are regenerated under ignored `data/derived/`. No new Rust dependency is required. Black, mypy, formatting, strict Clippy and meaningful numerical/missingness/held-out tests check the implementation.

The `power` subcommand reads the regenerated `data/derived/broad_power_design.csv`, uses 999 global and within-batch score permutations, 500 region bootstraps and 2,000 injections per effect size/generator. Compact outputs are `functional_power_audit.md`, `functional_power_summary.csv` and `functional_power_injection.csv`; full draws remain ignored. Injection sizes describe oracle available recovery-fraction signal, not kinetic-rate effects or guaranteed fitted improvement. Power conditions on this observed single-culture design and treats scores as error-free. Current data do not jointly support the revised external-fitness/longitudinal-clone transfer test; historical biological results must not be relabeled as its validation.

## Dog Aging Project falsification audit

The [report](results/dog_aging_falsification.md) and [frozen protocol](../notes/dog_aging_protocol.md) distinguish an essential-gene stability proxy from actual local methylation necessity and revision-3 clone fitness. Fifteen checksum-frozen sources include 1,640-sample regional methylated/total counts, author sample/batch metadata, canFam4 RefSeq/CpG islands/reference sequence, NCBI orthology and the existing BAGEL lists (~1.5 GB compressed/packed). Mutable orthology and UCSC annotations must not silently update. All downloaded data and observation caches remain ignored.

With the existing pinned format environment:

```sh
make dog-download       # network; acquire missing files, verify all hashes
make dog-validation     # offline verification/metadata decoding, Rust analysis and power audit
make check
data/tools/format-env/bin/black --check code/dog_aging_data.py
data/tools/format-env/bin/mypy --ignore-missing-imports code/dog_aging_data.py
```

`dog_aging_data.py` only acquires/verifies inputs and decodes the R dataframe. Rust in `src/dog_aging.rs` validates every matrix count/ID, verifies one-based CpG coordinates against 2,000 reference starts, joins unique orthologs, matches using training first visits, evaluates withheld dogs, resamples whole dogs/gene pairs and runs conditional endpoint-power simulations. No new package or Rust dependency is introduced. `dog-aging [raw] [results] primary|all` runs the flank-context audit; `dog-aging-audit [raw] [results]` uses regenerated observation caches for diagnostics and conditional power. `all` includes the four prespecified coverage/age sensitivities. Summary files are overwritten on rerun; annotations and matching are deterministic, seeds explicit.

The initial interval-only sequence match had uninformative GC/CpG covariates at singleton CpGs. Its negative estimate is preserved under `results/dog_aging_interval_context/`. To reproduce it separately, run `dog-aging data/raw code/results/dog_aging_interval_context interval`; this also writes its own orthology/dog tables. It uses the same source files and the original RNG sequence, with a separate `interval_observations` cache. The corrected flank results are transparently **post-result**, not an untouched preregistration. Different read gates reverse the point estimate; no preferred sensitivity is promoted to validation or rejection. All promoter comparisons fall far below the 20-pair support gate.

Conditional power uses dog-level wild residuals with the observed control mean, 10,000 calibration draws and 10,000 draws per 0/5/10/20% endpoint effect, a 98.75th-percentile one-sided null cutoff, and seed 20261010. It assumes independent dogs and fixed gene masks; it is not power for clone selection, a kinetic rate or direct locally validated methylation effects. Binomial variance correction cannot remove pooled-CpG read dependence, overdispersion or batch/composition confounding. Numeric reports are regenerated; the narrative assessment and protocol remain authored records.
