# Continued validation (2026-10-03)

These exploratory tests strengthen the external **stability** signal, but do not validate essentiality, allocated maintenance, selection or biological identity loss. All gains compare predictions in excluded donors; positive MSE gain means lower error. Blood participants are twins: family IDs are unavailable in GEO, so person-held-out splits and person-bootstrap intervals are provisional.

| Test | Units | Mean | 95% bootstrap |
|---|---:|---:|---:|
| GSE179847 forecast donor penalty=0.01 MSE gain | 9 | 0.040449 | [0.035722, 0.045022] |
| GSE179847 forecast donor_and_chromosome penalty=0.01 MSE gain | 9 | 0.040943 | [0.036283, 0.045593] |
| GSE179847 full matching coarse ESL-minus-control | 9 | -0.010169 | [-0.011106, -0.009261] |
| GSE179847 full matching fine ESL-minus-control | 9 | -0.003169 | [-0.003373, -0.002945] |
| GSE179847 full matching no_density ESL-minus-control | 9 | -0.005922 | [-0.006311, -0.005518] |
| GSE179847 full matching standard ESL-minus-control | 9 | -0.005734 | [-0.006134, -0.005306] |
| GSE179847 m2_protection log density gain vs M1 | 9 | 0.003209 | [0.001161, 0.005458] |
| GSE179847 m2_signed log density gain vs M1 | 9 | 0.003209 | [0.001260, 0.005380] |
| GSE73115 forecast donor penalty=0.01 MSE gain | 86 | 0.004446 | [0.003829, 0.005099] |
| GSE73115 forecast donor_and_chromosome penalty=0.01 MSE gain | 86 | 0.004728 | [0.004007, 0.005450] |
| GSE73115 full matching coarse ESL-minus-control | 86 | -0.003108 | [-0.003253, -0.002968] |
| GSE73115 full matching fine ESL-minus-control | 86 | -0.001296 | [-0.001352, -0.001242] |
| GSE73115 full matching no_density ESL-minus-control | 86 | -0.001495 | [-0.001548, -0.001442] |
| GSE73115 full matching standard ESL-minus-control | 86 | -0.001611 | [-0.001668, -0.001554] |
| GSE73115 m2_protection log density gain vs M1 | 86 | 0.004952 | [0.004114, 0.005790] |
| GSE73115 m2_signed log density gain vs M1 | 86 | 0.006314 | [0.005267, 0.007287] |
| growth ESL penalty=0.01 MSE gain | 8 | 0.003763 | [0.001269, 0.006166] |
| growth combined penalty=0.01 MSE gain | 8 | -0.016709 | [-0.114721, 0.062785] |
| growth unweighted penalty=0.01 MSE gain | 8 | -0.019694 | [-0.111578, 0.060673] |

## Interpretation and limits

- Forecasts use every 256th input probe, detection filtering, and complete external sequence/context annotations: 1,530 fibroblast probes (85 ESL) and 1,682 blood probes (88 ESL). The response is log(abs(endpoint change)+0.001). The empirical ridge baseline includes starting beta splines, sequence, density, CGI/regulatory class, chromosome, duration/start time and fibroblast condition/oxygen proxies. Its augmented version adds one external ESL indicator. Fixed penalties 0.001/0.01/0.1 are sensitivity settings, not tuned on held-out outcomes. Five chromosome groups and five deterministic donor groups are jointly excluded for the stricter split. Shared CpGs remain in the donor-only test. This is an empirical predictive gate, not a mechanistic M2 likelihood. Sex, array batch, genome-wide density and blood composition remain unadjusted.
- Pooled kinetic tests use the exact binary CTMC marginal, individual earlier beta as the initial probability, and a Gaussian bulk-array residual. They restrict both ESLs and controls to earlier beta <=0.2, and hold out five donor groups. M1 pools gain/return rates; M2 multiplies gain by exp(-alpha*q) and return by exp(beta*q). Both positive-only and signed effects are fit from two starts after three M1 starts. Rates use 100 days/years as their unit. Site/context heterogeneity is not fully modeled here; q coefficients reaching optimization bounds and protection/restoration tradeoffs forbid parameter-level biological interpretation. Small predictive gains are evidence for this restricted pooling model only. In blood, allowing signed effects fits alpha around -0.93 to -0.98 (higher gain at ESLs) and beta around 2.05 to 2.14 (higher return), and scores better than positive-only protection. Thus the observed stability does not specifically support suppressed away rates. Fibroblast beta repeatedly reaches plus/minus 6; restoration is not identified.
- Full-matrix matching varies beta, density and sequence bin widths together (coarse/standard/fine), plus omission of density. These runs use cohort early means and the same external context as the initial analysis; they remain transductive. The sampled, donor-training-only matching exercise is recorded separately in matching_sensitivity.csv: sparse strict strata often have zero support. Unsupported comparisons are NaN, not zero effects.
- Growth is a first independent biological endpoint: log(late/early measured population-doubling time), using the same DNA sample times. Only 11 culture endpoint pairs from eight donors have finite positive doubling measurements at both times. Forecasts compare time/starting-doubling-time covariates against additional unweighted and ESL absolute methylation distances, with donors excluded. This small pilot tests association with proliferation slowing, not identity loss; the ESL distance is not an independently defined functional weight. No threshold, causal mechanism or selection-versus-protection claim follows.
- Bootstrap intervals resample donor-level fixed test scores and omit uncertainty from refitting and spatial probe dependence. Blood twins require family mapping for valid family-held-out inference. No multiple-comparison correction is claimed; all settings/results are retained.

The correct paired fibroblast RNA study is GSE179848 (downloaded and checksummed); the previously downloaded GSE225172 expression file belongs to the mouse study. RNA outcomes are not analyzed in this follow-up.

Run `make followup` after `make reproduce` to regenerate these extended tests. Detailed comparisons are in forecast_comparison.csv, pooled_kinetic_comparison.csv, growth_pairs.csv, growth_prediction.csv, full matching CSVs and followup_summary.csv.
