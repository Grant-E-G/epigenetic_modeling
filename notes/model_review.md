# Preimplementation conceptual review (2026-10-03)

The CTMC foundation and falsification strategy are coherent. The following are material corrections or limits, not evidence against the biological hypothesis.

1. Independent binary CpG switching does **not** produce the specified Jacobi noise term. For N independent alleles, drift is gain*(1-x)-loss*x and instantaneous variance is [gain*(1-x)+loss*x]/N. Wright–Fisher noise requires resampling/genetic drift or another population mechanism. Use the CTMC as the base; do not call the two diffusions equivalent.
2. The protection rate parameterization is not equivalent to the proposed diffusion parameterization: both CTMC rates jointly determine equilibrium and relaxation. Increasing restoration changes equilibrium as well as total rate. Estimate and report both.
3. An allele is binary; a diploid cell can have two different alleles and bisulfite sequencing ordinarily does not distinguish 5mC from 5hmC. Single-cell read counts need an explicit allele/coverage interpretation. Missing coverage is not zero methylation.
4. Survival-weighted cross-sectional distributions usually cannot identify intrinsic maintenance versus selection. Fit survivor observations with survival normalization, and never interpret an apparent protection coefficient as causal. Serial sampling of different cells is not a tracked lineage.
5. Absolute distance and exp(gamma*distance) depend on the number and coverage of sites. Freeze the site universe, youthful reference and weight normalization. Missingness may bias distances. A young-reference quantile is an operational threshold, not established loss of cell identity.
6. The maintenance variable is unspecified in M4b, and M5 needs a bounded generative interpretation. An Ornstein–Uhlenbeck log-rate modifier is an explicit exploratory choice, not an implication of the specification. Latent covariance can also reflect batch, cell cycle, cell composition or read depth.
7. ESLs reflect observed stability, gene essentiality is tissue dependent, and Polycomb binding does not imply that lower methylation error improves fitness. None independently proves a functional constraint allocation.
8. Two time points cannot separate away and return rates without pooling and strong assumptions. Hierarchical pooling should expose sensitivity to those assumptions. Unobserved individual baseline distributions introduce further confounding.

Implementation order: obtain the three priority datasets; reproduce their basic summaries and compare a one-state/slow-fast baseline with mouse-held-out prediction; then implement stochastic model components and simulation checks. Real constraint fitting and identity interpretation remain gated on external annotations and functional outcomes. User-requested Rust replaces the specification's suggested Python project layout.

Initial access failure: sandbox DNS could not resolve ftp.ncbi.nlm.nih.gov. An approved network escalation succeeded before model code was written. Data are in the git-ignored /data directory.
