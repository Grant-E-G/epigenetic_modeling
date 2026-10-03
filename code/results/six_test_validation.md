# Six public-data validation tests (2026-10-03)

The extension finds useful predictive behavior and several failures of the simplest models. It **does not establish functional constraint or biological identity loss**. The functional RNA comparison is negative in the tested settings; kinetic curves do not outperform the pooled constant comparator; future-array intervals are undercalibrated. Clone structure and non-additive perturbation effects require explicit treatment. These results qualify the earlier stability-only evidence rather than negate its observed transfer.

All metrics and sensitivity settings are retained in six_test_summary.csv and source-specific CSVs. Positive prediction gains mean improvement; negative RNA MSE gains mean worse prediction. The table uses the middle fixed ridge penalty 0.01 for display only.

| Comparison | Units | Mean | 95% resampling interval |
|---|---:|---:|---:|
| RNA_ECM tolerance=32 penalty=0.01 functional minus matched random MSE gain | 7 | 0.000143 | [-0.000247, 0.000664] |
| RNA_ECM tolerance=32 penalty=0.01 functional_ECM MSE gain | 7 | -0.000952 | [-0.002310, 0.000282] |
| RNA_ECM tolerance=8 penalty=0.01 functional minus matched random MSE gain | 4 | 0.000025 | [-0.000044, 0.000095] |
| RNA_ECM tolerance=8 penalty=0.01 functional_ECM MSE gain | 4 | -0.000018 | [-0.000171, 0.000137] |
| RNA_ECM_absolute tolerance=32 penalty=0.01 functional minus matched random MSE gain | 7 | -0.000232 | [-0.000372, -0.000097] |
| RNA_ECM_absolute tolerance=32 penalty=0.01 functional_ECM MSE gain | 7 | -0.000725 | [-0.002188, 0.000335] |
| RNA_ECM_absolute tolerance=8 penalty=0.01 functional minus matched random MSE gain | 4 | 0.000073 | [-0.000064, 0.000308] |
| RNA_ECM_absolute tolerance=8 penalty=0.01 functional_ECM MSE gain | 4 | -0.000540 | [-0.001496, -0.000003] |
| RNA_global tolerance=32 penalty=0.01 functional minus matched random MSE gain | 7 | -0.000068 | [-0.000106, -0.000028] |
| RNA_global tolerance=32 penalty=0.01 functional_ECM MSE gain | 7 | -0.000112 | [-0.000518, 0.000136] |
| RNA_global tolerance=8 penalty=0.01 functional minus matched random MSE gain | 4 | 0.000005 | [-0.000015, 0.000038] |
| RNA_global tolerance=8 penalty=0.01 functional_ECM MSE gain | 4 | -0.000028 | [-0.000099, 0.000010] |
| clone quality=0.95 class=Static log score gain | 2 | 0.056824 | [0.051185, 0.062464] |
| identifiability endpoint_bulk fraction with multiple competitive candidates | 24 | 0.916667 | [0.791667, 1.000000] |
| identifiability longitudinal_cells fraction with multiple competitive candidates | 24 | 0.875000 | [0.708333, 1.000000] |
| identifiability longitudinal_cells_and_survival fraction with multiple competitive candidates | 24 | 0.833333 | [0.666667, 0.958333] |
| identifiability truth=both endpoint_bulk strict highest-score candidate frequency | 8 | 0.500000 | [0.125000, 0.875000] |
| identifiability truth=both longitudinal_cells strict highest-score candidate frequency | 8 | 0.250000 | [0.000000, 0.500000] |
| identifiability truth=both longitudinal_cells_and_survival strict highest-score candidate frequency | 8 | 1.000000 | [1.000000, 1.000000] |
| identifiability truth=protection endpoint_bulk strict highest-score candidate frequency | 8 | 0.250000 | [0.000000, 0.500000] |
| identifiability truth=protection longitudinal_cells strict highest-score candidate frequency | 8 | 0.250000 | [0.000000, 0.625000] |
| identifiability truth=protection longitudinal_cells_and_survival strict highest-score candidate frequency | 8 | 0.375000 | [0.000000, 0.750000] |
| identifiability truth=selection endpoint_bulk strict highest-score candidate frequency | 8 | 0.125000 | [0.000000, 0.375000] |
| identifiability truth=selection longitudinal_cells strict highest-score candidate frequency | 8 | 0.500000 | [0.125000, 0.875000] |
| identifiability truth=selection longitudinal_cells_and_survival strict highest-score candidate frequency | 8 | 0.875000 | [0.625000, 1.000000] |
| kinetics pulse_uniform_CTMC_equivalent hour=16 log score gain vs pooled constant | 1623 | -0.028123 | [-0.051261, -0.004579] |
| kinetics pulse_uniform_CTMC_equivalent hour=4 log score gain vs pooled constant | 1623 | -0.071641 | [-0.096595, -0.046279] |
| perturbation W464A_W465A_W796A independent combination minus calibrated log score | 3 | -0.000005 | [-0.000222, 0.000171] |
| perturbation W465A_W796A independent combination minus calibrated log score | 3 | -0.091342 | [-0.129932, -0.069046] |
| trajectory early_constant RMSE | 9 | 0.111381 | [0.097905, 0.127507] |
| trajectory early_constant fraction outside nominal 95 percent | 9 | 0.122306 | [0.089451, 0.157486] |
| trajectory independent_CTMC RMSE | 9 | 0.065530 | [0.051754, 0.084309] |
| trajectory independent_CTMC fraction outside nominal 95 percent | 9 | 0.211684 | [0.156855, 0.291642] |
| trajectory linear_empirical RMSE | 9 | 0.084149 | [0.063117, 0.109281] |
| trajectory linear_empirical fraction outside nominal 95 percent | 9 | 0.287423 | [0.220506, 0.365619] |

## 1. Independent functional outcomes

A frozen external Reactome extracellular-matrix pathway (R-HSA-1474244) supplies 307 annotated genes and 2,698 supported promoter probes; RNA includes 321 pathway genes. Probe weights balance genes. Twenty control masks sample other gene promoters with replacement in the same chromosome/CGI/probe-density/sequence strata, carrying the same weights. No methylation or RNA outcome selects the mask. This is a functional-program proxy, not essentiality or an identity-loss boundary.

Control fibroblast DNA at time t forecasts RNA at t+1..16 days. Both assays use their independent earlier culture references; no RNA outcome or held-out donor selects methylation weights. Starting expression, earlier beta, culture time, forecast horizon, SURF1 status and oxygen are baseline covariates. The primary baseline-match tolerance is 8 days (33 pairs/four donors); the exploratory 32-day sensitivity retains 63 pairs/seven donors and still enforces RNA baseline earlier than prediction time. HC5/HC6 lack corresponding RNA culture groups. Outcomes are signed mean ECM change, absolute ECM transcript displacement and absolute displacement of an input-order transcript subset. The absolute ECM outcome and wider tolerance were added during the audit, so are exploratory. Most functional-weight gains are negative; matched masks also often fail. This does not support the proposed functional-distance prediction in this cohort/definition. It does not disprove all functional weights. Age/sex/batch/proliferation adjustment and RNA future-increment forecasting remain incomplete.

## 2. Direct post-replication kinetics

The corrected author commit ec46d95b68fad06072702686c248cad55754fc3b provides observed methylated/unmethylated counts at post-pulse 0/1/4/16 hours. Fixed every-32nd-site sampling precedes coverage filtering (at least ten reads at zero and five at each later time). There are 2,959 covered CpGs on chr1 and chr22. Each of 4h and 16h is excluded in turn. We fit f*(1-exp(-k*t)) with either half-hour timing correction or exact one-hour pulse averaging; this marginal is mathematically CTMC-equivalent but cannot separate a mixed plateau from reversible maintenance. A training-pooled constant has better aggregate held-out score than these kinetics. Thus this test gives no predictive advantage for the tested time dependence. Many rates reach boundaries and only one biological cell line is represented. Hour-scale nascent-strand rates are not the year/day-scale drift rates inferred from aging arrays.

## 3. Clone and representation controls

The public EPI-Clone M.1–M.3 object has 28,782 cells. Tests use independent LARRY barcodes and raw antibody UMI counts from the two stained main-experiment recipients. Zero detected antibody molecules are retained only in nonempty stained libraries. Six surface-marker covariates plus antibody-library size and undigested-control quality adjust cell state/measurement. Probe amplification in independent undigested controls must be at least 95%. Cell quality thresholds 0.95 and 0.8 retain 1,903 and 3,285 cells respectively; 370 probes pass assay controls. Methylation-sensitive-digestion amplification is a proxy with dropout, not clean allele-specific bisulfite calls.

Held-out cells within independently labelled clones benefit from a shrunken clone residual, especially at published static CpGs. Published static/dynamic classes use this cohort's proteins and are descriptive, not independent feature validation. Equal-clone and cell-weighted averages use the same qualifying clone universe. A symmetric exact decomposition across shared clones divides mouse3-to-mouse4 bulk differences into within-clone changes and clone-frequency changes; these are **different recipients**, not a longitudinal survival experiment. Clonal representation matters, but selection causing aging stability is not identified. Human raw reads are controlled access and were not used.

## 4. Existing DNMT1 perturbations

GSE145698 contains three WT cultures and three cultures for each of five DNMT1 mutants. Fixed every-32nd-row sampling and WT coverage yield 165,106 autosomal CpGs. WT counts define initial beta; mutant training chromosomes calibrate an effective CTMC loss exposure, then opposite-parity chromosomes are excluded for evaluation. An affine empirical comparator is virtually as predictive: calibrated scaling does not validate a unique mechanism. No elapsed culture exposure is given, so the fitted parameter is a rate-times-exposure product.

An additional genotype-held-out test predicts W465A/W796A from the fitted separate W465A and W796A exposures, and the triple mutant from W464A/W465A plus W796A. The simple independent-effects combination substantially underperforms mutation-specific calibration for W465A/W796A, consistent with domain interactions. This rejects that transfer rule, not every CTMC or the general drift framework. No mutation-specific data calibrate the held-out combination rule.

## 5. Full predictive distribution

Control fibroblast trajectories (fixed every-256th probe) fit their earlier three quarters and predict strictly later arrays. 54,083 covered CpG/culture trajectories produce 50 future array summaries. Independent CTMC means improve RMSE over initial-constant and linear forecasts, but Gaussian predictive intervals miss far more than the nominal 5%. Training residual variance plus a 0.02 array-noise floor omits parameter uncertainty and changing variance; this is a failure of the current predictive distribution, not proof that intrinsic CTMC dynamics are impossible. Covariance analysis uses a fixed 128-probe subset and independent within-donor residual permutations: shared residual structure remains. Batch, culture, oxygen and common proliferation changes can generate it. This is an exploratory rejection of unconditional independence, not a validated biological M5 rank. No threshold is called identity loss.

## 6. What observations identify

An exact eight-state killed chain retains three loci and shared state-dependent killing. Three generating mechanisms, three observation architectures and eight synthetic repeats are compared with independent/protection/selection/both candidates and independent test simulations. There are 400 sampled survivors/timepoint, 40% missing locus calls and 2% call error. Bulk endpoints expose only locus marginals; repeated destructive cell snapshots expose joint patterns; adding independent survival counts supplies 500 at-risk observations/timepoint. Snapshot-only candidates often confuse protection and selection; survival helps recover selection and combined mechanisms. Strict highest-score winner frequencies are sensitive to nearly tied scores; identifiability_winners.csv also lists every candidate within 0.001 mean log score of the best. This tolerance was added during the numerical audit, not preregistered. Multiple candidates remain competitive in 22/24 bulk-endpoint runs, 21/24 repeated-cell runs and 20/24 survival-augmented runs. Survival improves some strict rankings but does not demonstrate robust mechanism identification at this sampling depth. Nested candidates can overfit pure mechanisms even with survival. Baseline killing (0.05), targets, error and missingness are known and no clone reproduction/site-context confounding is simulated. This measures conditional recoverability, not identifiability for unrestricted biology.

## Reproduction and remaining gates

All source files live under ignored data/raw and are checksummed. A pinned Python helper only decodes MATLAB/R objects to TSV; all inference, matching, predictions, simulations and reports are Rust. `make validation` converts formats, rebuilds functional/trajectory summaries, runs all six tests and regenerates this report. Feature assembly is hg19 for human probes; mouse assays are never joined to human coordinates.

Independent functional weights beyond ECM promoters, RNA/identity thresholds, family-aware blood inference, clone reproduction and survival likelihoods for real aging cohorts, full observation uncertainty and intervention-specific molecular recruitment still require work. No broad causal biological-validation claim follows from these tests.
