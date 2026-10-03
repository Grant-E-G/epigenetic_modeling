# Initial results (2026-10-03)

The model components are implemented in Rust and numerical checks pass. The real-data run reproduces aggregate drift summaries, but does **not** validate functional masking or identify maintenance versus selection biologically.

## Data and observations

| Dataset | Analyzed source | Biological units | Result |
|---|---|---|---|
| GSE225171 | Published processed Polycomb summary, 1,131 cells | 11 mice | Mean beta 0.05523 → 0.07019 from 10 to 100 weeks; SD 0.00768 → 0.01971. Matches reported summaries closely. |
| GSE179847 | Full 865,817-probe processed matrix, 479 arrays | Nine donors / 67 culture groups | Control-group donor-averaged endpoint absolute change 0.06889 beta units; 95% replicate bootstrap [0.06047, 0.07917]. |
| GSE73115 | Full 485,577-probe processed matrix, 180 arrays | 86 people / 43 twin pairs | Ten-year mean absolute paired change 0.01767; 95% replicate bootstrap [0.01697, 0.01844]. |

The single-cell raw archives were also downloaded, but genomic Polycomb averages were not independently rebuilt from coverage. The published summary contains one fewer cell than the article reports. Array pairs require detection p ≤ 0.01 at both endpoints; blood technical duplicates are pooled per person/year after this filter. Fibroblast groups preserve treatment/culture identity. Neither aggregate measure is an intrinsic cellular error rate.

## Held-out model comparisons

Final Polycomb comparison includes an animal-level Gaussian random intercept, with five-point quadrature and a continuous-switch/Gaussian weekly-gain approximation. It holds out complete mice and weights training animals equally.

| Comparison | Mean predictive-density gain | 95% replicate bootstrap |
|---|---:|---:|
| Polycomb slow group: two-state minus one-state, joint score per cell | 0.10267 | [0.01728, 0.21456] |
| Polycomb fast group: two-state minus one-state | -0.00030 | [-0.02222, 0.02951] |
| Fibroblast M1 minus constant young reference, 14 fixed probes | 0.08117 | [-0.21724, 0.44143] |
| Blood M1 minus constant young reference, eight fixed probes | 0.00604 | [-0.00475, 0.01831] |

The approximate slow-group fit gives a positive signal; the fast-group result does not reproduce a reliable two-state advantage. Do not interpret either as validation of a global maintenance state. The baseline is numerically approximate, includes broad proliferative classes rather than a full lineage model, and has only 11 mice. The first version omitted a shared mouse effect; those preliminary metrics are superseded. Bootstrap intervals resample fixed held-out scores rather than rerunning fitting and therefore omit training/optimizer uncertainty. M1 gains are inconclusive on the small fixed subsets.

Predictive mean/variance/tail and Wasserstein diagnostics are in `posterior_predictive_checks/polycomb.csv`; despite that specification-required directory name, these are fitted predictive checks, not posterior draws. A 0.089 beta tail cutoff reproduces a published descriptive threshold and is **not** an identity-loss boundary.

## Numerical mechanism checks

At one weighted locus and horizon 5, simulated survivor methylation for pure protection, pure selection and both is respectively 0.1664, 0.1424 and 0.0414. Both protection-only and selection-only candidate fits match each survivor endpoint to about 1e-5, although simulated survival differs greatly (20,000; 3,427; 5,964 survivors out of 20,000). Thus this one-snapshot observation cannot identify the mechanism. This is a counterexample, not comprehensive recovery testing for all multisite longitudinal datasets.

Jump-process tests agree with analytic transition, killed-process survival and first-passage distributions. A synthetic q-dependent kinetic fit recovers endpoint probabilities. Continuous OU simulation at dt 0.1 / 0.025 / 0.00625 gives first-passage fractions 0.6518 / 0.6406 / 0.6390 in 5,000 trajectories; coarse discretization visibly changes crossing estimates. No simulated passage is claimed to be actual biological identity loss.

## External annotation and matched-mask follow-up

ESL and sequence supplements were obtained after initial access errors. The feature table includes 896,166 human probes and all 31,744 published unmethylated ESLs. Motif scores average forward/reverse-complement ranks; CGI/regulatory matching uses external 450K platform context. Unknown ranks/contexts and strata without both groups are excluded, with coverage reported.

Using chromosome, CGI class, coarse regulatory class, early beta (bins 0.05), probe-sequence density (bins 0.01) and sequence score (bins 0.1), ESLs have lower absolute endpoint change than matched controls:

| Dataset | Matched ESL sites | Donor-averaged ESL-minus-control change | 95% donor bootstrap |
|---|---:|---:|---:|
| Fibroblast controls | 21,259 | -0.00573 | [-0.00614, -0.00531] |
| Ten-year blood | 17,450 | -0.00161 | [-0.00167, -0.00155] |

Both observed site-averaged contrasts lie beyond all 500 matched random-mask draws (finite permutation p = 1/501). This supports transfer of a **stability annotation** under this exploratory coarse matching scheme. It does not prove functional essentiality, causal maintenance allocation, or protection versus selection. Matching uses cohort early means rather than fully subject-held-out preprocessing, local density is measured over a probe sequence rather than a full genomic window, and residual batch/composition/spatial effects and matching-bin sensitivity remain unresolved. Source annotations were fixed before inspecting these effects; this follow-up was added after the initial baseline tests.

A denser random-intercept evaluation changes the two-minus-one-state mean scores to 0.10447 (slow) and 0.01662 (fast). Some individual scores shift substantially (maximum 0.424 per cell); the fitted parameters still use five-point quadrature. Treat numerical fitting sensitivity as unresolved.

## Unresolved scientific requirements

Exploratory ESL F2/F8 checks have run with imported external sequence/context/stability annotations. Held-out F1 mechanistic prediction and independently validated essentiality/identity scores have not been established. No outcome-derived functional mask was substituted. F5/F6/F7 remain untested on real residual structure, paired functional/RNA outcomes and clocks. Bulk composition, sex, batch and proliferation adjustments are incomplete. Real-data protection/selection identification, functional threshold validation and comprehensive biological validation therefore remain unfinished.

The literature review is in [notes/reading_list.md](../../notes/reading_list.md), with prior-art distinctions in [novelty_baselines.md](novelty_baselines.md). The combination is candidate novelty; hierarchical fidelity inference, stochastic clocks, latent switching and epigenetic first passage are already established ingredients.

Run `make check` and `make reproduce` from the repository root. The latter recreates core reports from checksum-verified local downloads.

## Continued validation

See [extended results](validation_followup.md) for larger donor/chromosome-held-out stability forecasts, restricted pooled M1/M2 tests, full-matrix matching sensitivity and a measured proliferation pilot. These supersede the earlier statement that all held-out F1/matching sensitivity work was absent, but do not complete the functional-constraint gate. The blood cohort consists of twin pairs: person bootstrap intervals above do not account for family dependence, and person-held-out predictions may train on the co-twin. Family mappings are missing from downloaded GEO metadata.

## Six public-data tests

[The six-test report](six_test_validation.md) adds independent future RNA, direct maintenance-count kinetics, clone-conditioned/representation analyses, existing DNMT1 perturbations, complete forecast distributions and multisite mechanism recovery. Several scientific gates fail: ECM promoter weighting does not improve RNA prediction, kinetic curves do not beat a pooled constant, and the current predictive intervals are too narrow. Clone structure and DNMT1 domain interactions matter. Numerical correctness and the earlier stable-locus signal remain supported; broad biological validity is not established.


## Functional recovery pilot and experiment direction

The subsequent HCT116 recovery pilot uses independently defined gene-fitness labels. Its corrected unfiltered association is positive, but baseline expression detection QC leaves seven matched pairs and tighter expression matching leaves two with inconsistent recovery contrasts. It does not establish preferential functional restoration or expression recovery. [The experiment brief](../../notes/experiment_brief.md) records the full assessment, four candidate anchors, a local-function editing gate, and revised mechanism accounting. `make recovery` reproduces the new Rust analysis without changing the earlier six-test results.
