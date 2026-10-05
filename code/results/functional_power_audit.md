# Functional-effect detectability audit

Primary audit: signed HCT116 score, day 40, ridge penalty 0.01. Other days/penalties are sensitivities. Positive gain means lower held-out MSE. The linear-weight implementation exactly reproduces the frozen training-only scaling, ridge penalties and whole-region exclusion. The design is regenerated from checksum-verified source inputs; neither injected outcomes nor resampling alter coverage/eligibility or missingness. These estimates condition on one observed recovery culture and noisy observed functional scores; they are not biological-culture power estimates.

## Observed effect and empirical null

| Analysis | Day | Penalty | Regions | Observed gain | Global permutation p | Within-batch permutation p | Bootstrap 95% range |
|---|---:|---:|---:|---:|---:|---:|---:|
| primary | 28 | 0.001 | 39 | -1.57% | 0.296 | 0.411 | -17.45% .. 3.94% |
| primary | 28 | 0.01 | 39 | -1.56% | 0.285 | 0.415 | -16.70% .. 5.44% |
| primary | 28 | 0.1 | 39 | -1.47% | 0.313 | 0.486 | -13.71% .. 5.11% |
| primary | 40 | 0.001 | 39 | -3.89% | 0.542 | 0.798 | -24.56% .. 7.54% |
| primary | 40 | 0.01 | 39 | -3.83% | 0.548 | 0.792 | -22.63% .. 4.98% |
| primary | 40 | 0.1 | 39 | -3.37% | 0.567 | 0.841 | -16.04% .. 6.32% |
| speed | 6 | 0.001 | 40 | -4.90% | 0.665 | 0.723 | -20.37% .. 6.50% |
| speed | 6 | 0.01 | 40 | -4.88% | 0.660 | 0.715 | -22.92% .. 6.30% |
| speed | 6 | 0.1 | 40 | -4.59% | 0.716 | 0.729 | -18.76% .. 4.59% |
| speed | 13 | 0.001 | 42 | -3.96% | 0.576 | 0.534 | -21.12% .. 3.44% |
| speed | 13 | 0.01 | 42 | -3.93% | 0.621 | 0.552 | -19.40% .. 5.27% |
| speed | 13 | 0.1 | 42 | -3.63% | 0.600 | 0.526 | -17.61% .. 4.32% |

## Primary injection results

| Generator | Target oracle reduction | Achieved signal fraction | Effect per score SD | Any MSE gain | Calibrated detection | MC SE |
|---|---:|---:|---:|---:|---:|---:|
| unbounded | 0% | 0.00% | 0.00000 | 16.2% | 4.3% | 0.46% |
| unbounded | 3% | 3.00% | 0.00636 | 42.8% | 34.7% | 1.06% |
| unbounded | 5% | 5.00% | 0.00830 | 57.9% | 50.1% | 1.12% |
| unbounded | 10% | 10.00% | 0.01206 | 80.5% | 74.5% | 0.97% |
| unbounded | 20% | 20.00% | 0.01810 | 98.8% | 98.2% | 0.30% |
| physical_bounds | 0% | 0.00% | 0.00000 | 17.8% | 5.5% | 0.51% |
| physical_bounds | 3% | 3.00% | 0.00636 | 45.4% | 37.4% | 1.08% |
| physical_bounds | 5% | 5.00% | 0.00830 | 58.1% | 49.2% | 1.12% |
| physical_bounds | 10% | 9.99% | 0.01206 | 83.0% | 75.9% | 0.96% |
| physical_bounds | 20% | 19.99% | 0.01810 | 99.0% | 98.4% | 0.28% |

## Fixed-design diagnostics

| Analysis | Day | Regions | Score SD | Residual score SD | Score variance explained by context | Largest region share of residual score energy |
|---|---:|---:|---:|---:|---:|---:|
| primary | 28 | 39 | 0.47264 | 0.46780 | 2.04% | 29.93% |
| primary | 40 | 39 | 0.47264 | 0.46780 | 2.04% | 29.93% |
| speed | 6 | 40 | 0.47714 | 0.47272 | 1.84% | 28.43% |
| speed | 13 | 42 | 0.46920 | 0.46352 | 2.41% | 27.67% |

## Interpretation and limits

An injected beta acts on recovery fraction conditional on the real covariates: y*=context mean + beta*q_res + a signed region residual. It is not a DNMT rate change. Target oracle reductions 3/5/10/20% specify available fixed-design signal relative to the observed context leave-one-out residual variance. They do not assert that ridge achieves those reductions. q_res is the OLS projection residual of the actual observed score; score reliability is treated optimistically as perfect. Context means are ridge fits to this culture. Independent regional wild signs are an assumed error generator, not observed replicate biology. Read-depth-dependent residual magnitudes and actual missingness/filters remain fixed.

Physical-bound sensitivity clips only simulated methylation to [0,1], using each region's observed day3 and induced loss, then converts back to recovery fraction. It does not clip recovery fractions to [0,1]. The achieved signal fraction measures incremental simulated mean variation after context projection relative to simulated noise; clipping can alter the nominal oracle signal. This ratio is an operational signal/noise diagnostic, not a known true biological effect. Calibrated detection requires positive global fitted score coefficient and MSE gain above a separate zero-injection 95th percentile (floored at zero); zero-injection assessment reports realized false positives. MC standard errors describe simulations only.

Global score permutations are not generally exchangeable conditional on context; within-screen-batch permutations address supplied batch grouping but not all biological confounding. Their p-values are empirical diagnostics, not causal significance. Bootstrap intervals resample complete region records and refit, excluding every copy of the held-out gene; they describe regional heterogeneity conditional on this culture. This resampling reduces the unique training-region count, so intervals are not a formally calibrated independent-culture confidence interval. All three penalties and four days are reported without selecting a successful variant.

A low powered audit cannot permanently falsify functional restoration. Even high conditional detectability would only downgrade the tested measured-score/region/regression relationship, not all kinetic effects, functional definitions, tissues or selection. A null protection test does not validate selection.
