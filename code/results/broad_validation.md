# Public-data stopping tests (2026-10-04)

The [stopping assessment](../../notes/experiment_brief.md#computational-stopping-assessment-2026-10-04) distinguishes negative functional predictions from generic kinetic plausibility. Protocol frozen in notes/research_plan.md before new recovery outcomes. Inputs verified against broad_data_manifest.csv and broad_reference_manifest.csv. Python handles acquisition, verification, XLSX and reference coordinate formats only; all estimates and predictions are Rust.

## Independent local function versus human recovery

56 screen targets; 3 fail the conservative two-guide/unique-assembly gate. Mean log viability is normalized against negative controls within construct and plate, then averaged equally across available batches. Each target belongs to one supplied batch with eight wells; the three batch IDs distribute targets rather than replicate each target. No biological significance test is inferred from wells. Signed cost is inactive minus active; absolute cost is an explicit sensitivity. VP64 and DLD1 are reported separately.

- WT drug day28: 42 eligible regions; median recovery fraction 0.9789, quartiles 0.9578 / 1.0033.
- WT drug day40: 41 eligible regions; median recovery fraction 0.9931, quartiles 0.9759 / 1.0145.
- WT vehicle day40: 41 eligible regions; median recovery fraction 1.0011, quartiles 0.9814 / 1.0242.
- DNMT3B KO day40: 37 eligible regions; median recovery fraction 0.9098, quartiles 0.8648 / 0.9564.
- DNMT1 KO day40: 32 eligible regions; median recovery fraction 0.8614, quartiles 0.7928 / 0.9138.
- primary: 53 sequence-joined regions, 5 fail early coverage; 39 have sufficient perturbation and training coverage, 39 / 39 have day-28 / day-40 outcomes. Median early methylation loss 0.4142.
- coverage5: 53 sequence-joined regions, 3 fail early coverage; 43 have sufficient perturbation and training coverage, 43 / 43 have day-28 / day-40 outcomes. Median early methylation loss 0.4262.
- coverage20: 53 sequence-joined regions, 9 fail early coverage; 34 have sufficient perturbation and training coverage, 34 / 32 have day-28 / day-40 outcomes. Median early methylation loss 0.4210.
- baseline03: 53 sequence-joined regions, 5 fail early coverage; 41 have sufficient perturbation and training coverage, 41 / 41 have day-28 / day-40 outcomes. Median early methylation loss 0.4142.
- baseline07: 53 sequence-joined regions, 5 fail early coverage; 39 have sufficient perturbation and training coverage, 39 / 39 have day-28 / day-40 outcomes. Median early methylation loss 0.4142.
- loss010: 53 sequence-joined regions, 5 fail early coverage; 39 have sufficient perturbation and training coverage, 39 / 39 have day-28 / day-40 outcomes. Median early methylation loss 0.4142.
- loss020: 53 sequence-joined regions, 5 fail early coverage; 39 have sufficient perturbation and training coverage, 39 / 39 have day-28 / day-40 outcomes. Median early methylation loss 0.4142.

Positive MSE changes mean worse held-out prediction. Primary scores and gate sensitivities were frozen before outcomes. Early recovery-speed, baseline-only context, DLD1 score, VP64 and inactive-construct comparisons below were added as exploratory robustness checks after the first run; none replaces the primary test.

| Gate / score | Day | Penalty | Regions | Context RMSE | Signed cost MSE change | Absolute cost MSE change |
|---|---:|---:|---:|---:|---:|---:|
| primary | 28 | 0.001 | 39 | 0.03563 | 1.57% | 1.67% |
| primary | 28 | 0.01 | 39 | 0.03552 | 1.56% | 1.63% |
| primary | 28 | 0.1 | 39 | 0.03466 | 1.47% | 1.31% |
| primary | 40 | 0.001 | 39 | 0.03550 | 3.89% | 2.97% |
| primary | 40 | 0.01 | 39 | 0.03537 | 3.83% | 2.93% |
| primary | 40 | 0.1 | 39 | 0.03433 | 3.37% | 2.58% |
| exploratory_speed | 6 | 0.001 | 40 | 0.13208 | 4.90% | 2.26% |
| exploratory_speed | 6 | 0.01 | 40 | 0.13201 | 4.88% | 2.24% |
| exploratory_speed | 6 | 0.1 | 40 | 0.13174 | 4.59% | 2.09% |
| exploratory_speed | 13 | 0.001 | 42 | 0.14284 | 3.96% | 3.75% |
| exploratory_speed | 13 | 0.01 | 42 | 0.14270 | 3.93% | 3.73% |
| exploratory_speed | 13 | 0.1 | 42 | 0.14194 | 3.63% | 3.49% |
| exploratory_baseline_context | 28 | 0.001 | 39 | 0.03368 | 2.22% | 2.18% |
| exploratory_baseline_context | 28 | 0.01 | 39 | 0.03365 | 2.19% | 2.14% |
| exploratory_baseline_context | 28 | 0.1 | 39 | 0.03339 | 1.91% | 1.80% |
| exploratory_baseline_context | 40 | 0.001 | 39 | 0.03080 | 3.66% | 2.34% |
| exploratory_baseline_context | 40 | 0.01 | 39 | 0.03077 | 3.61% | 2.31% |
| exploratory_baseline_context | 40 | 0.1 | 39 | 0.03052 | 3.24% | 2.10% |
| exploratory_dld1 | 28 | 0.001 | 39 | 0.03563 | 4.56% | 9.77% |
| exploratory_dld1 | 28 | 0.01 | 39 | 0.03552 | 4.39% | 9.61% |
| exploratory_dld1 | 28 | 0.1 | 39 | 0.03466 | 2.89% | 8.25% |
| exploratory_dld1 | 40 | 0.001 | 39 | 0.03550 | 5.75% | 3.43% |
| exploratory_dld1 | 40 | 0.01 | 39 | 0.03537 | 5.64% | 3.28% |
| exploratory_dld1 | 40 | 0.1 | 39 | 0.03433 | 4.82% | 2.21% |
| exploratory_vp64 | 28 | 0.001 | 39 | 0.03563 | 8.01% | -0.32% |
| exploratory_vp64 | 28 | 0.01 | 39 | 0.03552 | 7.84% | -0.41% |
| exploratory_vp64 | 28 | 0.1 | 39 | 0.03466 | 6.49% | -1.10% |
| exploratory_vp64 | 40 | 0.001 | 39 | 0.03550 | 7.07% | 10.94% |
| exploratory_vp64 | 40 | 0.01 | 39 | 0.03537 | 6.76% | 10.63% |
| exploratory_vp64 | 40 | 0.1 | 39 | 0.03433 | 4.59% | 8.48% |
| exploratory_inactive | 28 | 0.001 | 39 | 0.03563 | 5.97% | 2.61% |
| exploratory_inactive | 28 | 0.01 | 39 | 0.03552 | 5.90% | 2.51% |
| exploratory_inactive | 28 | 0.1 | 39 | 0.03466 | 5.21% | 1.71% |
| exploratory_inactive | 40 | 0.001 | 39 | 0.03550 | 2.50% | 5.37% |
| exploratory_inactive | 40 | 0.01 | 39 | 0.03537 | 2.33% | 5.22% |
| exploratory_inactive | 40 | 0.1 | 39 | 0.03433 | 1.16% | 4.22% |
| coverage5 | 28 | 0.001 | 43 | 0.06391 | 3.21% | 3.96% |
| coverage5 | 28 | 0.01 | 43 | 0.06372 | 3.16% | 3.87% |
| coverage5 | 28 | 0.1 | 43 | 0.06252 | 2.77% | 3.26% |
| coverage5 | 40 | 0.001 | 43 | 0.07644 | -1.03% | -6.46% |
| coverage5 | 40 | 0.01 | 43 | 0.07608 | -1.03% | -6.35% |
| coverage5 | 40 | 0.1 | 43 | 0.07369 | -1.04% | -5.53% |
| coverage20 | 28 | 0.001 | 34 | 0.05872 | 4.53% | -0.77% |
| coverage20 | 28 | 0.01 | 34 | 0.05852 | 4.48% | -0.80% |
| coverage20 | 28 | 0.1 | 34 | 0.05692 | 4.03% | -0.97% |
| coverage20 | 40 | 0.001 | 32 | 0.04728 | 4.95% | 5.82% |
| coverage20 | 40 | 0.01 | 32 | 0.04708 | 4.87% | 5.70% |
| coverage20 | 40 | 0.1 | 32 | 0.04545 | 4.23% | 4.73% |
| baseline03 | 28 | 0.001 | 41 | 0.05099 | 0.74% | -0.49% |
| baseline03 | 28 | 0.01 | 41 | 0.05080 | 0.74% | -0.51% |
| baseline03 | 28 | 0.1 | 41 | 0.04940 | 0.75% | -0.67% |
| baseline03 | 40 | 0.001 | 41 | 0.04438 | 1.75% | 0.89% |
| baseline03 | 40 | 0.01 | 41 | 0.04418 | 1.74% | 0.85% |
| baseline03 | 40 | 0.1 | 41 | 0.04280 | 1.59% | 0.54% |
| baseline07 | 28 | 0.001 | 39 | 0.03563 | 1.57% | 1.67% |
| baseline07 | 28 | 0.01 | 39 | 0.03552 | 1.56% | 1.63% |
| baseline07 | 28 | 0.1 | 39 | 0.03466 | 1.47% | 1.31% |
| baseline07 | 40 | 0.001 | 39 | 0.03550 | 3.89% | 2.97% |
| baseline07 | 40 | 0.01 | 39 | 0.03537 | 3.83% | 2.93% |
| baseline07 | 40 | 0.1 | 39 | 0.03433 | 3.37% | 2.58% |
| loss010 | 28 | 0.001 | 39 | 0.03563 | 1.57% | 1.67% |
| loss010 | 28 | 0.01 | 39 | 0.03552 | 1.56% | 1.63% |
| loss010 | 28 | 0.1 | 39 | 0.03466 | 1.47% | 1.31% |
| loss010 | 40 | 0.001 | 39 | 0.03550 | 3.89% | 2.97% |
| loss010 | 40 | 0.01 | 39 | 0.03537 | 3.83% | 2.93% |
| loss010 | 40 | 0.1 | 39 | 0.03433 | 3.37% | 2.58% |
| loss020 | 28 | 0.001 | 39 | 0.03563 | 1.57% | 1.67% |
| loss020 | 28 | 0.01 | 39 | 0.03552 | 1.56% | 1.63% |
| loss020 | 28 | 0.1 | 39 | 0.03466 | 1.47% | 1.31% |
| loss020 | 40 | 0.001 | 39 | 0.03550 | 3.89% | 2.97% |
| loss020 | 40 | 0.01 | 39 | 0.03537 | 3.83% | 2.93% |
| loss020 | 40 | 0.1 | 39 | 0.03433 | 3.37% | 2.58% |

A coverage-qualified late regional mean requires at least three CpGs and 80% of the frozen early CpG set, each with the specified read coverage. CpGs receive equal weight. Recovery fractions may overshoot; none are clipped. Late vehicle uses the same WT reference and diagnoses untreated culture drift. KO fractions use genotype-specific early reference on the same WT-defined sites. There is one recovery culture per genotype/time: CpGs are not biological replicates.

Prediction CSVs contain whole-region-excluded ridge predictions at three fixed penalties, with training-only feature scaling. Other regions' outcomes at the evaluated day calibrate the regression: this is cross-region generalization, not a simultaneous global future-time holdout or a new-condition transfer test; models 0/1/2 are early recovery/context, plus signed cost, plus absolute cost. Empty prediction tables mean the frozen >=12-region feasibility threshold failed. H3K36me3, expression, division and clone measurements are unavailable in this regional join, so any surviving association would remain mechanistically confounded. Active-versus-inactive editing also lacks a measured per-target editing-efficiency control.

## Replicated mouse enzyme-deletion forecasts

Rates fit only Cre days 4/8/10/13 in two cultures; predictions exclude the third culture and forecast days 17/29 from its measured day-0 baseline. Three rotations share training data and are not three independent experiments. Equal CpG/time/culture squared error is the fitting criterion; read counts are coverage filters, not independent biological samples. No published confidence label or fitted rate selects our subset. The source table itself was filtered by its authors using all time points and excludes one high-variability amplicon; this limits the independence of the retrospective validation. An exploratory two-parameter empirical comparator fits an amplitude and exponential loss rate from day 4 onward, without using the held-out baseline. It was added after the initial amplicon results and does not identify enzyme clearance. WT has TET activity; TTKO lacks TET enzymes. Both undergo de-novo DNMT deletion. Constant CTMC and decaying de-novo activity have two site-specific rates; fixed clearance 0.5/day is the published enzyme-clearance mechanism, 0.25/1.0 are prespecified sensitivities. Grid rates include zero and log10 -4..0 in 0.1 steps, maximum 1/day.

### WT: 336 CpGs in 55 amplicons

| Model | Future RMSE | Equal-amplicon MSE | Upper-rate-bound forecasts |
|---|---:|---:|---:|
| clearance025 | 0.07230 | 0.004665 | 278 / 2016 |
| clearance050 | 0.07543 | 0.004539 | 186 / 2016 |
| clearance100 | 0.08106 | 0.004960 | 398 / 2016 |
| constant_ctmc | 0.13948 | 0.016006 | 840 / 2016 |
| delayed_empirical_loss | 0.07971 | 0.005012 | 0 / 2016 |
| initial_constant | 0.15699 | 0.022729 | 0 / 2016 |
| last_training | 0.11051 | 0.009987 | 0 / 2016 |

### TTKO: 346 CpGs in 52 amplicons

| Model | Future RMSE | Equal-amplicon MSE | Upper-rate-bound forecasts |
|---|---:|---:|---:|
| clearance025 | 0.06695 | 0.004805 | 94 / 2076 |
| clearance050 | 0.05902 | 0.003714 | 54 / 2076 |
| clearance100 | 0.05568 | 0.003308 | 124 / 2076 |
| constant_ctmc | 0.08015 | 0.007213 | 530 / 2076 |
| delayed_empirical_loss | 0.06347 | 0.004654 | 0 / 2076 |
| initial_constant | 0.19574 | 0.043947 | 0 / 2076 |
| last_training | 0.09001 | 0.009452 | 0 / 2076 |

## Exploratory genome-wide enzyme benchmark

Fixed every-256th input-row sampling retains 8507 / 8507 sampled CpGs with >=50 reads at all seven times. The processed mouse TTKO table aggregates cultures, so this is a temporal holdout, not an additional biological replication study. Methylation values are percentages converted to fractions. Published fitted rates, TET labels and confidence classifications are ignored. The processed source already requires coverage across all times and reflects the authors' data filtering; this is not an untouched raw-data evaluation. Training uses days 4/8/10/13; evaluation uses days 17/29. This broader benchmark was added after the amplicon results and is exploratory.

| Model | Future RMSE |
|---|---:|
| clearance025 | 0.05173 |
| clearance050 | 0.04444 |
| clearance100 | 0.04273 |
| constant_ctmc | 0.06648 |
| delayed_empirical_loss | 0.04639 |
| initial_constant | 0.20105 |
| last_training | 0.09597 |

Mock controls are summarized separately by amplicon/culture/time. Clearance fits test an already published intervention mechanism, not functional protection. Amplicons are correlated CpG blocks and were deliberately selected by the original authors; hundreds of CpGs do not constitute hundreds of independent biological replicates. No uncertainty interval assumes otherwise.


## Mock-control diagnostic

Equal amplicon/culture mean changes from day 0, in methylation fractions:

| Genotype | Day | Mean change |
|---|---:|---:|
| TTKO | 10 | 0.0116 |
| TTKO | 13 | 0.0231 |
| TTKO | 17 | 0.0224 |
| TTKO | 29 | -0.0108 |
| TTKO | 4 | 0.0059 |
| TTKO | 8 | 0.0075 |
| WT | 10 | 0.0702 |
| WT | 13 | 0.0394 |
| WT | 17 | 0.0311 |
| WT | 29 | 0.0222 |
| WT | 4 | 0.1012 |
| WT | 8 | 0.1022 |

WT mock methylation rises during the early transduction period. Primary predictions use uncorrected Cre fractions; clearance and empirical transient fits can absorb this background. The near-competitive empirical loss curve and mock drift prevent treating forecast gains as unique evidence for the exact clearance equation or as isolated TET effects. TTKO controls drift much less. No genotype-transfer rule or functional-restoration law is established here.

