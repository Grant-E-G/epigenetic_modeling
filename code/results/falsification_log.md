# Falsification log (append-only)

## 2026-10-03: initial run

- Preimplementation review rejects equivalence between independent-site CTMC aggregate noise and Jacobi/Wright–Fisher noise. The microscopic base is a CTMC.
- NCBI DNS access failed inside the sandbox, then approved public-data network access succeeded before coding.
- Published Polycomb data contain 1,131 cells instead of 1,132. Retain observed rows, do not synthesize the missing cell.
- Referenced MATLAB MC helper absent. Baseline uses explicitly labeled continuous-switch/Gaussian approximation; exact optimizer reproduction remains incomplete.
- Blood sample naming includes technical replicate suffixes. Pool duplicate arrays within person/year after detection filtering.
- Initial held-out slow/fast comparisons used animal-balanced marginal fits. A final audit adds a shared Gaussian mouse intercept with joint likelihood quadrature, as required by the design's nested-sampling rule. Final metrics are in replicate_bootstrap.csv; early metrics are superseded.
- F1/F2/F8 real functional constraint and matched covariate-mask tests NOT RUN: no frozen feature table containing sequence-fidelity/context and independent q. No positive or negative constraint result may be claimed.
- F3 one-locus survivor snapshots fail to identify protection versus selection in numerical counterexamples. This is not a universal impossibility proof for richer data.
- F4 approximate Polycomb baseline evaluated with held-out mice. Interpret final uncertainty intervals, not in-sample fits.
- F5/F6/F7 real residual dimension, independent identity/function and clock comparisons NOT RUN. No biological first-passage threshold fitted.
- Continuous OU/factor variants are explicit exploratory choices. Their first-passage observations are grid resolved, with step-size results recorded.

## 2026-10-03: external annotations and matching follow-up

- bioRxiv supplement endpoint returned HTTP 429 and the PMC PDF endpoint returned HTML. A subsequent public bioRxiv DC1 download supplied a valid PDF; its 225 motif ranks were imported. No claim that the source remains inaccessible is warranted.
- ESL supplement supplied all 31,744 unmethylated CpG IDs. Human platform annotations already reside in the downloaded SOFT files. Derived features use hg19 coordinates, probe-sequence density and strand-symmetrized ranks; unknown contexts/ranks are excluded rather than imputed.
- Add exploratory F2/F8 tests specified in the design: within-stratum matched control contrasts and 500 seeded mask permutations, with donor-level effect summaries. This updates the initial NOT RUN status for ESL F2/F8 only; F1 mechanistic prediction, essentiality/identity masks and functional validation remain open.
- Matching strata preserve chromosome, CGI class, coarse promoter/body/other class and bins of early beta (0.05), CpG density (0.01) and sequence rank (0.1). Matching is on cohort early means; do not call it fully independent subject-held-out validation. Check alternative bins, sequence scores and spatial nulls before treating the effect as robust functional evidence.
- Initial follow-up pooled fibroblast treatments and used one strand of sequence ranking. Final follow-up uses control groups and symmetrized rankings; initial matched metrics are superseded.
- Dense held-out random-intercept quadrature gives mean two-minus-one-state score 0.10447 (slow) and 0.01662 (fast), versus 0.10267 and -0.00030 under the fitting quadrature. Individual score discrepancies can be material (up to 0.424 per cell). Fitting-quadrature sensitivity remains unresolved; do not overinterpret the slow-group signal or fast-group sign.

## 2026-10-03: continued real-data validation

- Empirical held-out F1 stability forecasts now run on a fixed denser probe subset with individual earlier beta and external sequence/context covariates. ESL improves the mean donor-held-out log-absolute-change MSE in both cohorts at all three fixed ridge settings; the sign survives joint donor/chromosome exclusion. This passes a limited external-stability predictive gate, not independent functional-constraint validation.
- Full-matrix coarse/standard/fine and no-density matching all retain negative ESL-minus-control contrasts. Finer matching reduces the magnitude and overlap: 14,902 fibroblast and 10,524 blood ESLs retain full fine-stratum support. Matching on cohort early beta remains transductive.
- Restricted low-baseline exact-CTMC M2 has small held-out improvements over pooled M1. Fibroblast restoration estimates hit boundaries. Blood signed-rate M2 scores better by allowing **increased away rates** together with faster return. Suppression of away rates is therefore not specifically supported by these fits. Protection/restoration parameter identification fails under this pooling architecture.
- Independent proliferation endpoint pilot joins actual measured doubling time at the DNA endpoints. Only 11 cultures/eight donors qualify. ESL distance gives a small positive held-out prediction gain across fixed penalties; unweighted/combined distance has no consistent benefit. Gain magnitude varies strongly with penalty, and fixed-score bootstrap excludes training uncertainty. This is suggestive and inadequate for biological identity validation.
- Blood cohort participants are 43 twin pairs rather than 86 independent families. GEO lacks family IDs. Existing person-held-out CV/bootstrap may train on co-twins and underestimate uncertainty. Do not invent a family mapping or treat these intervals as family-valid inference.
- Correct fibroblast paired RNA is GSE179848; its metadata and processed matrix were obtained without access problems and checksummed. GSE225172 is mouse blood RNA. RNA assays have not been analyzed in this follow-up; no transcriptomic validation claim is licensed.
