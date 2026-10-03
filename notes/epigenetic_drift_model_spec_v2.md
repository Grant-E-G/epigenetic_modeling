# Epigenetic Drift as a Constrained Stochastic Process
## Model specification, novelty check, and falsification-first analysis plan

**Purpose:** handoff specification for a coding/research agent.

**Core objective:** build the smallest family of stochastic methylation models that can be *invalidated* on real data. Prediction of chronological age is secondary. The primary questions are whether methylation drift is differentially constrained at functionally important loci, whether apparent constraint reflects better maintenance versus selective removal/clonal filtering, and whether a low-dimensional loss of maintenance state explains coordinated drift.

**Version:** 2026-10-02

---

## 1. Executive decision

Do **not** frame the project as "a novel SDE model of epigenetic drift." That claim is already crowded. Stochastic differential-equation models of age-related methylation maintenance existed by 2021, explicit stochastic CpG-transition models and stochastic clock analyses appeared in 2024, and a two-state slow/fast Polycomb methylation model was published in 2026.

The potentially novel part is narrower and better:

> **Model site-specific functional constraint as a latent/external property that can act through two competing mechanisms - enhanced epigenetic maintenance versus selective loss or clonal under-representation of cells that violate important methylation states - and test these mechanisms by explicit model comparison on longitudinal and single-cell data. Define loss of identity as a first-passage event in a function-weighted epigenetic state space.**

Treat this as **candidate novelty**, not proven novelty. The literature search below is a targeted novelty screen, not a systematic review or patent search.

The project's success criterion is not a good-looking simulation. It is one of the following:

1. the constraint model survives serious falsification and predicts held-out data better than simpler models;
2. the protection-versus-selection distinction becomes empirically identifiable in at least one dataset;
3. the model fails, but the failure exposes a reproducible structure in methylation drift that is not captured by current stochastic-clock or two-state models.

Any of those is useful.

---

## 2. Novelty check

### 2.1 What is already in the literature

| Ingredient | Novelty status | Closest prior work / implication |
|---|---|---|
| SDEs for age-related DNA methylation | **Not novel** | Zagkos et al. (2021) explicitly introduced white-noise/SDE-style stochasticity into an age-related DNA maintenance methylation model and studied site-specific instability. |
| Binary CpG transitions / stochastic clocks | **Not novel** | Tong et al. (2024) model binary methylation switches at individual CpGs and show that a large fraction of clock accuracy can arise from stochastic change. |
| Probabilistic cell-level transition model | **Not novel** | The 2024 Nature Aging "probabilistic inference" paper models methylation transitions at the cellular level and separates acceleration from bias. |
| Slow/fast latent aging state | **Not novel** | Masika et al. (2026) use a two-state kinetic model: cells stochastically switch from a slow Polycomb methylation-gain state to a faster state. |
| Differential site maintenance fidelity | **Not novel** | Site- and region-specific methylation fidelity has been documented for decades. Promoter CpG islands can be more protected from de novo methylation errors than non-promoter islands. |
| Sequence-dependent maintenance fidelity | **Very current, not novel** | Lopez-Moyado et al. (2026 preprint) report that flanking hexanucleotide sequence ranks predict DNMT1/UHRF1 maintenance fidelity and age/replication-associated loss. This must be treated as a baseline covariate, not rediscovered as "functional masking." |
| Stable/protected CpG sets | **Not novel** | A 2026 Nature Communications study identifies tens of thousands of epigenetically stable loci (ESLs) in young healthy blood and shows increasing destabilization with age/disease. |
| Drift plus clonal selection | **Conceptually not novel** | Issa (2014) explicitly argued that methylation drift creates mosaicism on which somatic selection can act, yielding stem-cell extinction and clonal expansion. |
| Aging as an evolving stochastic epigenetic landscape | **Not novel in broad form** | Recent trajectory/landscape approaches model cell state with stochastic dynamics in an age-dependent potential. |

### 2.2 Candidate novelty that survived the search

I did **not** find a paper that clearly combines all of the following as an inferential model of methylation aging:

1. an externally defined or hierarchically inferred **functional constraint score** for CpG sites/regions;
2. two explicitly competing ways that constraint can produce apparent stability:
   - **protection:** lower error rate / stronger restoration at important loci;
   - **selection:** cells that drift at important loci are preferentially lost, fail to expand, or are under-sampled among surviving lineages;
3. a latent **maintenance-fidelity state** producing correlated drift across many sites;
4. an identity-loss boundary defined in a **function-weighted methylation state space**;
5. **first-passage time** to that boundary as the cell-level aging quantity;
6. falsification by model comparison across **single-cell, longitudinal, and bulk** datasets rather than optimization for age prediction.

The strongest candidate scientific contribution is therefore not the SDE. It is the **identifiability and falsification framework for differential epigenetic error tolerance**.

### 2.3 Claims to avoid

Do not claim any of the following unless later citation chaining establishes them:

- "first SDE model of epigenetic aging";
- "first stochastic model of CpG methylation";
- "first model of protected methylation sites";
- "first two-state model of epigenetic aging";
- "first connection between epigenetic drift and selection".

A defensible future claim, if the work succeeds, would look more like:

> We introduce a falsifiable hierarchical stochastic model that separates maintenance-mediated constraint from survival/clonal selection at functionally weighted methylation loci and evaluates identity loss as a first-passage process.

---

## 3. Biological question

The motivating hypothesis is:

> Evolution and cellular regulatory machinery do not tolerate methylation errors equally across the genome. Loci whose methylation state is important for preserving cell identity or viability should exhibit lower effective epigenetic error, stronger restoration after perturbation, stronger selection against deviating cells, or some combination of the three.

The critical distinction is:

$$
\text{observed stability} \neq \text{intrinsic maintenance fidelity}.
$$

A CpG can look stable in old tissue because:

1. it rarely changes;
2. it changes but is rapidly restored;
3. cells that change it die, senesce, differentiate, or fail to expand;
4. population composition changes in a way that hides the drift;
5. measurement architecture makes it appear stable.

The model must keep these alternatives separate for as long as the data permit.

---

## 4. Modeling principles

### 4.1 Separate biological state from measurement

At the microscopic level, a CpG allele in a cell is approximately binary:

$$
M_{c i}(t) \in \{0,1\}.
$$

A bulk-array beta value is **not** a fractional state of a single cell. It is an aggregate/measurement quantity. Therefore use different observation models for:

- single-cell bisulfite data;
- WGBS read counts;
- bulk methylation arrays;
- clonal/crypt-level methylomes.

Do not fit one continuous SDE directly to every beta matrix and interpret its diffusion coefficient as cellular noise.

### 4.2 Prefer a discrete kinetic model as the generative base

For site $i$, use a two-state continuous-time Markov process or a division-indexed transition process:

$$
0 \underset{\lambda_{10,i}}{\overset{\lambda_{01,i}}{\rightleftarrows}} 1.
$$

When an aggregate diffusion approximation is useful, derive it from this process rather than postulating additive Brownian motion.

### 4.3 Keep a bounded phenomenological approximation available

For aggregate methylation $X_i(t) \in [0,1]$, a bounded Jacobi/Wright-Fisher-like diffusion is preferable to unconstrained Brownian motion:

$$
dX_i = \kappa_i(\theta_i-X_i)dt + \sigma_i\sqrt{X_i(1-X_i)}\,dW_i.
$$

Here:

- $\theta_i$: youthful/reference methylation state or local equilibrium;
- $\kappa_i$: restoring/maintenance strength;
- $\sigma_i$: effective stochasticity.

This is a **phenomenological layer**, not automatically the microscopic mechanism.

---

## 5. Functional constraint: replace the binary "mask" with a score

Define a site or region score

$$
q_i \in [0,1]
$$

representing putative functional constraint.

Avoid learning $q_i$ from the same age-associated variance that the model is trying to explain. Start with external annotations and only later allow a hierarchical latent component.

### 5.1 Candidate external components of $q_i$

Use several independent definitions and test robustness rather than committing to one:

- **Epigenetically Stable Loci (ESLs):** loci selected for low variability in young healthy samples, from the 2026 blood instability study.
- **DNMT1/UHRF1 sequence-fidelity score:** flanking hexamer rank from Lopez-Moyado et al. 2026. This is especially important as a mechanistic confounder/baseline.
- **Regulatory element class:** ENCODE cCRE promoter-like, enhancer-like, CTCF-only, etc.
- **Promoter CpG island / imprinting control region / transposable-element context.**
- **Gene essentiality:** map promoter/enhancer-associated CpGs to an independently defined common-essential gene set (for example DepMap common essentials), with tissue relevance treated cautiously.
- **Cell-identity relevance:** lineage-specific transcription-factor and regulatory programs defined from young cells, not from the aged methylation outcome being modeled.
- **Evolutionary conservation:** optional, and preferably used as a separate covariate rather than assumed equivalent to functional importance.

Construct multiple $q_i$ versions:

- `q_esl`
- `q_sequence`
- `q_regulatory`
- `q_essentiality`
- `q_identity`
- a prespecified composite

The model should be considered fragile if the result appears only for one hand-crafted definition.

---

## 6. Competing mechanistic models

Implement these as a nested hierarchy. Do not jump directly to the maximal model.

### M0 - empirical null

No mechanistic SDE. Estimate age-associated mean change and variance change with appropriate covariates and hierarchical structure.

Purpose: establish what needs explanation.

Core outputs:

- $d E[X_i]/dt$
- $d \operatorname{Var}(X_i)/dt$
- age-dependent covariance structure
- dependence on baseline methylation, CpG density, sequence context, cell type, replication status, sex, batch

### M1 - independent site kinetics

Each site has methylation and demethylation rates:

$$
\lambda_{01,i}, \lambda_{10,i}.
$$

Rates may depend on sequence and genomic context but **not functional constraint**.

This is the minimum mechanistic baseline.

### M2 - protection / restoration model

Functional constraint affects transition rates.

For a youthful target state $m_i^*$, parameterize an away-from-target rate and a recovery rate:

$$
\log \lambda_{\mathrm{away},i}
= a_0 + a^T z_i - a_q q_i,
$$

$$
\log \lambda_{\mathrm{return},i}
= b_0 + b^T z_i + b_q q_i.
$$

Here $z_i$ contains known mechanistic/context covariates such as flanking sequence, CpG density and regulatory class.

Protection predicts:

$$
a_q > 0 \quad \text{and/or} \quad b_q > 0.
$$

Equivalent diffusion form:

$$
\kappa_i = \kappa_0\exp(a_q q_i + a^Tz_i),
$$

$$
\sigma_i = \sigma_0\exp(-b_q q_i + b^Tz_i).
$$

### M3 - selection / censoring model

Allow methylation errors to occur but make cell survival/representation depend on function-weighted deviation:

$$
D_c(t)=\sum_i q_i\,d_i(M_{ci}(t),m_i^*),
$$

$$
h_c(t)=h_0(t)\exp(\gamma D_c(t)).
$$

Depending on the dataset, `h` can represent:

- cell death;
- senescence/dropout;
- failure to contribute descendants;
- reduced clone size;
- effective removal from the sampled population.

**Important:** in ordinary cross-sectional methylation data M2 and M3 may be non-identifiable. The code must report this rather than selecting one by force.

### M4 - latent maintenance-state model

Introduce a cell-level maintenance variable $R_c(t)$ that changes many sites simultaneously.

Two useful variants:

#### M4a: two-state switching

$$
R_c(t) \in \{Y,O\}, \qquad Y \xrightarrow{\rho} O.
$$

Rates differ by state:

$$
\lambda_{i}^{(O)} = g_i\lambda_i^{(Y)}.
$$

This directly nests/competes with the 2026 Polycomb slow/fast model.

#### M4b: continuous latent fidelity

$$
dR_c = f(R_c,t)dt + \eta dW_c^R,
$$

with methylation kinetics depending on $R_c$.

Do not assume the continuous version is superior. Let held-out prediction decide.

### M5 - correlated residual drift

After conditioning on $R_c$, test whether substantial residual covariance remains. If yes, use a low-rank factor model rather than a dense covariance matrix:

$$
d\mathbf X = \boldsymbol\mu(\mathbf X,R,q)dt + B\,d\mathbf W.
$$

Keep rank small and cross-validate it.

---

## 7. Identity loss and first-passage formulation

Define a youthful/reference methylation vector $\boldsymbol\theta$, estimated from young cells of the **same lineage** or an independent reference set.

Define function-weighted distance:

$$
D_c(t) = (\mathbf X_c(t)-\boldsymbol\theta)^T Q(\mathbf X_c(t)-\boldsymbol\theta),
$$

where

$$
Q=\operatorname{diag}(q_1,\dots,q_p)
$$

for the first implementation.

Later, if justified, replace diagonal $Q$ with a sparse/network-aware metric.

Define loss of identity as a first-passage event:

$$
\tau_c = \inf\{t : D_c(t)>D_{\mathrm{crit}}\}.
$$

**Do not choose $D_{\mathrm{crit}}$ to maximize age prediction.** Candidate definitions should be externally grounded, for example:

- onset of a transcriptomic identity-loss signature;
- transition to a mesenchymal/stress program;
- loss of a lineage classifier margin;
- a prespecified quantile from young reference cells followed by validation in held-out data.

The first-passage quantity becomes scientifically interesting only if it predicts something other than age.

---

## 8. Observation models by data type

### 8.1 Single-cell methylome

For covered CpGs, use binary or read-count likelihoods directly. Model missingness explicitly enough to avoid treating missing CpGs as unmethylated.

If paired RNA is available, use RNA only for **validation/identity outcome**, not to fit the methylation distance and then claim independent functional relevance.

### 8.2 WGBS / targeted bisulfite counts

For site $i$, sample $s$:

$$
k_{si} \sim \operatorname{BetaBinomial}(n_{si},p_{si},\phi_i)
$$

or a binomial if overdispersion is negligible.

### 8.3 Bulk methylation arrays

Treat beta values as noisy aggregates over cells plus cell-composition effects. At minimum:

$$
y_{si}=x_{si}+\epsilon_{si}
$$

with array/batch error and estimated cellular composition included in the observation/regression layer.

Do not interpret bulk $\sigma_i$ as single-cell stochasticity.

### 8.4 Nested sampling

Cells from the same animal are **not independent biological replicates**. For GSE225171, the 1,132 cells come from only 11 mice. Use mouse-level random effects and leave-one-mouse-out validation. Report uncertainty at the animal level.

---

## 9. Falsification program

The coding agent should attempt to kill the model in the order below.

### Test F1 - does "constraint" explain anything after known local methylation physics?

Regress/fit drift statistics against:

- baseline methylation;
- local CpG density;
- CpG island/shore/shelf/open-sea status;
- flanking sequence maintenance score;
- chromatin/regulatory class;
- replication/proliferation proxy;
- cell type;
- batch/sex where applicable.

Then add $q_i$.

**Invalidate the constraint hypothesis** if $q_i$ provides no stable held-out gain or its coefficient changes sign unpredictably across datasets.

### Test F2 - protected versus permissive sites

Predefine protected/permissive groups using external annotations.

Compare age trajectories of:

- within-person/site variance;
- transition probability;
- restoring behavior after perturbation;
- extreme outlier frequency.

Use matched controls for baseline beta, CpG density, sequence motif, chromosome, and regulatory class.

A simple unmatched comparison is insufficient.

### Test F3 - protection versus survivor/clone selection

Look for signatures that distinguish lower error generation from post-error filtering.

Protection predicts fewer deviations at the point of generation.

Selection predicts that deviations can appear transiently or in small clones but are depleted among persistent/successful lineages.

Possible data signatures:

- clone-size dependence;
- lineage persistence versus methylation deviation;
- truncation/asymmetric tails in cross-sectional distributions;
- stronger deviation in dying/senescent cells than survivors;
- time-lagged disappearance of high-$D$ cells.

If existing datasets do not contain lineage/survival information, label M2 versus M3 **non-identifiable** and stop short of mechanistic claims.

### Test F4 - is a latent maintenance state actually needed?

Fit M1/M2 first. Examine residual covariance and age-dependent tail behavior.

Compare:

- one-state model;
- two-state switching model;
- continuous latent-state model.

Use held-out log likelihood / expected log predictive density, posterior predictive checks, and calibration of distribution tails - not only AIC on the same data.

Invalidate M4 if the latent state does not reproducibly improve prediction of full distributions, especially skewness/bimodality/tails.

### Test F5 - is drift low-dimensional?

After controlling for known covariates, examine the eigenvalue spectrum of residual age-associated covariance.

If one/few factors dominate and replicate across datasets, that supports common maintenance-state variables.

If covariance is high-rank and tissue/site specific, the global $R_t$ idea is wrong or incomplete.

### Test F6 - does function-weighted distance beat unweighted distance?

Compare:

$$
D_q=(X-\theta)^TQ(X-\theta)
$$

against

$$
D_0=\|X-\theta\|^2.
$$

The test outcome should be **cell identity/function**, not chronological age.

Candidate outcomes:

- paired RNA lineage-classifier confidence;
- inappropriate mesenchymal/stress program expression;
- loss of tissue-specific expression programs;
- disease/functional phenotype in external cohorts.

Invalidate the masking concept if weighting by $q$ does not outperform simple distance after proper cross-validation.

### Test F7 - clock comparison

Compare against established aging measures, but do not optimize solely for age correlation.

Ask whether model residuals explain:

- within-age heterogeneity;
- longitudinal change;
- paired expression/function;
- intervention response;
- mortality/disease where available.

A model that merely reconstructs chronological age adds little.

### Test F8 - negative-control masks

Generate many matched random masks preserving:

- number of CpGs;
- baseline methylation distribution;
- CpG density;
- genomic compartment;
- chromosome;
- sequence-fidelity score.

The real functional mask must beat this empirical null distribution.

This is one of the most important tests in the project.

---

## 10. Priority real datasets

### Tier 1 - use first

#### A. GSE225171 - single-cell mouse blood methylome + paired cell typing

Why it matters:

- 1,132 single cells;
- ages 10, 36, 77 and 100 weeks;
- paired transcriptomic information was used for cell type assignment;
- published 2026 Polycomb analysis shows skewed age distributions and motivates a two-state model.

Primary tests:

- reproduce the published one-state versus two-state result;
- then ask whether constraint-weighted site sets improve the model;
- leave-one-mouse-out validation;
- test methylation distance versus RNA-defined identity/stress programs.

Major caveat: only 11 mice. Never treat cells as 1,132 independent animals.

#### B. GSE179847 / GSE179849 - longitudinal human fibroblast lifespan

Why it matters:

- 479 methylation samples;
- roughly 7 time points per donor;
- 6 healthy donors plus 3 mitochondrial-disease donors;
- sampled every ~8 days across replicative lifespan;
- parallel RNA-seq and physiological measurements;
- includes metabolic/endocrine perturbations.

Primary tests:

- estimate within-lineage transition/drift parameters;
- test whether protected sites have lower away-rates or stronger restoration;
- test whether interventions act on $\sigma$, $\kappa$, latent $R$, or equilibrium $\theta$;
- exploit paired RNA as an external identity/function readout.

Caveat: in-vitro replicative aging is not organismal aging.

#### C. GSE73115 - 10-year repeated human blood methylation

Why it matters:

- same individuals measured roughly a decade apart;
- directly useful for within-person drift rather than cross-sectional inference.

Primary tests:

- protected-versus-permissive within-person variance/change;
- matched-mask permutation tests;
- model calibration over a long human interval.

Caveats: whole blood composition and only two main time points.

### Tier 2 - generalization / null-model calibration

#### D. GSE56581 - MESA sorted monocytes

About 1,202 donors across middle-to-older age. Useful because sorted cells reduce composition confounding.

#### E. GSE40279 - Hannum whole blood

Classic broad-age 450K cohort. Useful for external replication, not mechanistic identification.

#### F. GSE87571 - broad lifespan whole blood

421 individuals, ages roughly 14-94, useful for trajectory shape and genomic-context effects.

#### G. GSE130748 - longitudinal leukocyte methylation in older adults

Small, but useful as an additional paired-data stress test.

### Tier 3 - annotation and mechanistic external data

#### H. 2026 Epigenetically Stable Loci source data

Use the published 31,744 normally unmethylated ESLs (plus methylated ESLs if useful) as one independent definition of a protected set.

Important: ESLs were defined empirically by low variability in young blood, so they are a **stability annotation**, not proof of functional essentiality.

#### I. Lopez-Moyado et al. 2026 sequence-fidelity ranking

Use flanking-hexamer preference as a mandatory covariate or mechanistic prior. If "functional protection" disappears after controlling for this, that is a valuable negative result.

#### J. ENCODE Registry of cCREs

Use current cCRE annotations for promoter/enhancer/CTCF context and cell-type-specific regulatory annotation.

---

## 11. Minimum viable implementation for the coding agent

### Phase 0 - repository and reproducibility

Create:

```text
project/
  README.md
  pyproject.toml
  config/
    datasets.yaml
    annotations.yaml
    models.yaml
  src/
    data/
    annotations/
    models/
    inference/
    evaluation/
    plots/
  tests/
  notebooks/
  results/
    preregistered/
    exploratory/
  docs/
```

Requirements:

- Python 3.11+;
- deterministic seeds;
- environment lock file;
- every figure reproducible from a command;
- dataset checksums and versioned annotation sources;
- keep preregistered and exploratory analyses separate.

Suggested libraries: numpy, scipy, pandas/polars, scikit-learn, statsmodels, pymc or numpyro, arviz, xarray/anndata, pyarrow, matplotlib.

Do not make JAX/PyTorch mandatory until model size requires them.

### Phase 1 - reproduce known results before adding novelty

1. Download/process GSE225171.
2. Reproduce age increase and variance/skew in average Polycomb methylation.
3. Reproduce a one-state versus two-state model comparison closely enough to validate the pipeline.
4. Process GSE179847 and reproduce basic age/passage-associated methylation trajectories.
5. Process GSE73115 and reproduce longitudinal drift direction/variance statistics.

**Gate:** no new model work until these sanity checks succeed.

### Phase 2 - build annotations without outcome leakage

For every CpG create a feature table with:

```text
cpg_id
chrom
position
baseline_beta
cpg_density
cgi_context
regulatory_class
sequence_fidelity_score
esl_flag
essential_gene_link
identity_regulatory_score
prc2_flag
polycomb_cgi_flag
```

Freeze this table before fitting the age model.

### Phase 3 - empirical falsification first

Before SDE/CTMC fitting:

1. estimate age-associated mean and variance change;
2. compare protected/permissive sets with matched controls;
3. perform negative-control mask permutations;
4. test residual covariance rank;
5. quantify how much signal is explained by sequence context alone.

If the central functional-constraint signal dies here, document the failure and do not rescue it with a more flexible SDE.

### Phase 4 - fit M1 and M2

Start with hierarchical site kinetics. Use partial pooling by genomic class.

Required comparison:

- M1: no $q_i$;
- M2: $q_i$ modifies away/return rates.

Evaluate with held-out individuals/animals/donors, not random held-out CpGs only.

### Phase 5 - add M4 latent state only if residuals require it

Fit a two-state model first. Compare to a continuous latent process only if justified.

For GSE225171, compare distribution tails and skewness with leave-one-mouse-out prediction.

### Phase 6 - selection model / identifiability

Run simulation-based identifiability experiments before fitting real data:

1. simulate pure protection;
2. simulate pure selection;
3. simulate both;
4. pass through each dataset's observation process;
5. ask whether the inference pipeline can recover the generating mechanism.

If it cannot, mark the corresponding real-data mechanistic distinction as unidentifiable.

This simulation step is mandatory.

### Phase 7 - first-passage / identity analysis

Only after $q_i$ shows independent predictive value:

1. construct $D_q$ and $D_0$;
2. choose identity threshold without optimizing on the test set;
3. compare their prediction of RNA-defined loss of identity/function;
4. estimate first-passage distributions in fitted simulations;
5. assess whether first-passage behavior generalizes across donors/cell types.

---

## 12. Statistical evaluation

Use multiple metrics because different models can match the mean while missing the biology.

Required:

- held-out log predictive density;
- calibration of mean and variance by age;
- tail probability calibration;
- posterior predictive checks;
- distributional distances (for example Wasserstein distance) for cell-level age distributions;
- leave-one-donor/animal-out validation;
- matched permutation nulls for mask effects;
- bootstrap at the **biological replicate** level;
- sensitivity to cell-composition adjustment;
- sensitivity to baseline methylation and sequence-fidelity covariates.

Prefer effect sizes and predictive intervals over gigantic CpG-level P-value counts.

---

## 13. Explicit failure criteria

The project should report a negative result if any of the following hold robustly:

1. functional constraint does not predict drift after sequence/context matching;
2. different plausible definitions of $q_i$ give inconsistent signs;
3. weighted identity distance fails to outperform unweighted distance;
4. latent maintenance state does not improve held-out distributional prediction;
5. all apparent single-cell heterogeneity is explained by cell type, proliferation and animal-level effects;
6. M2 and M3 are non-identifiable under every available data modality;
7. model parameters do not transfer even qualitatively between datasets of comparable cell type;
8. the maximal model improves in-sample fit but loses out-of-sample predictive performance.

Do not weaken these criteria after seeing the data.

---

## 14. Exploratory analyses: where something new may pop out

These are explicitly exploratory and should be labeled as such.

### 14.1 Residual drift modules

Cluster CpGs by **residual drift covariance after known covariates are removed**, not by raw methylation correlation. Determine whether modules map to PRC2, CTCF, enhancer classes, replication timing, lamina association, or sequence-maintenance classes.

### 14.2 Error-budget hypothesis

Ask whether loci near essential/identity genes exhibit a tradeoff:

- low away-rate;
- high recovery-rate;
- low tolerated clone size when perturbed.

A tradeoff between maintenance and selection would be more interesting than a simple protected/unprotected dichotomy.

### 14.3 Asymmetric constraint

A site may tolerate drift in one direction but not the other. Fit separate methylation-gain and methylation-loss penalties.

### 14.4 Tissue-specific masks

A CpG can be vital for identity in one lineage and irrelevant in another. Test

$$
q_{i,c}\neq q_i.
$$

A strong cell-type interaction would argue against a universal genomic mask and for lineage-specific regulatory constraint.

### 14.5 Intervention parameter fingerprints

In the fibroblast perturbation dataset, ask whether different interventions primarily change:

- transition bias;
- stochasticity;
- restoration;
- latent state-switching probability.

Distinct parameter fingerprints could be more mechanistically informative than a clock-age change.

### 14.6 Extreme-value / first-passage signatures

If tissue dysfunction depends on a minority of cells crossing an identity boundary, tissue aging may depend on the distribution of earliest or cumulative passage events rather than mean methylation age. Explore this only after the cell-level model is validated.

---

## 15. What the coding agent should produce

The first complete run should output:

1. `novelty_baselines.md` - exact prior models reproduced and how the current model differs;
2. `data_manifest.csv` - datasets, accessions, sample counts, biological replicate counts, cell types, ages, data type;
3. `annotation_manifest.csv` - every annotation source/version and leakage status;
4. `empirical_drift_report.html` - M0 results and matched-mask tests;
5. `model_comparison.csv` - M1/M2/M4 predictive metrics;
6. `identifiability_report.html` - simulated protection versus selection recovery;
7. `posterior_predictive_checks/` - distributions, tails, residual covariance;
8. `falsification_log.md` - every failed prediction and model change, append-only;
9. `exploratory_findings.md` - discoveries not part of the preregistered tests;
10. a single `make reproduce` or equivalent command that regenerates core results.

The agent should not silently expand model complexity. Every added mechanism must correspond to a failed posterior predictive check or a preregistered biological hypothesis.

---

## 16. Recommended starting order

If computation is cheap and agent time is effectively free, parallelize data ingestion, but keep the inferential order disciplined:

1. **GSE225171:** reproduce the 2026 two-state Polycomb result and test animal-level robustness.
2. **GSE179847:** estimate actual longitudinal within-donor dynamics and test external functional masks.
3. **GSE73115:** test long-interval within-person drift.
4. Build the external constraint feature table.
5. Run matched empirical falsification before fitting the full stochastic model.
6. Fit M1 versus M2.
7. Add M4 only if residual distributional structure demands it.
8. Test M2 versus M3 identifiability in simulation before interpreting real data.
9. Only then construct the first-passage identity model.

This order makes it hard for a flexible model to talk us into believing the hypothesis.

---

## 17. Bottom line

The original idea - "an SDE with certain vital methylation sites masked" - is too close to existing stochastic methylation work and too biologically ambiguous if implemented literally.

The stronger project is:

> **A falsification-first model of differential epigenetic error tolerance. Functional constraint is allowed to manifest as enhanced maintenance, restoration, selection against deviant cells, or combinations thereof. A latent maintenance state accounts for correlated drift only if the data demand it. Loss of cell identity is tested as a function-weighted first-passage event rather than assumed.**

The most interesting result may be that the model is wrong. For example, if externally defined essential/identity loci are not unusually protected once flanking sequence and chromatin context are controlled, that would argue that apparent epigenetic "masking" is largely biophysical rather than allocated according to functional importance. Conversely, if functional constraint retains predictive power after these controls - especially in longitudinal and single-cell data - that is a much more consequential observation than another accurate aging clock.

---

## 18. Key references and data links

1. Zagkos L. et al. **A mathematical model which examines age-related stochastic fluctuations in DNA maintenance methylation.** Experimental Gerontology (2021). PubMed: https://pubmed.ncbi.nlm.nih.gov/34774717/
2. Tong H. et al. **Quantifying the stochastic component of epigenetic aging.** Nature Aging 4, 886-901 (2024). https://www.nature.com/articles/s43587-024-00600-8
3. **Probabilistic inference of epigenetic age acceleration from cellular dynamics.** Nature Aging (2024). https://www.nature.com/articles/s43587-024-00700-5
4. Schumacher B., Meyer D. **Aging clocks based on accumulating stochastic variation.** Nature Aging (2024). https://www.nature.com/articles/s43587-024-00619-x
5. Masika H. et al. **Cell-to-cell variability and gain of methylation at polycomb CpG islands as a hallmark of aging.** Nature Communications 17, 7318 (2026). https://doi.org/10.1038/s41467-026-74118-5
6. Lopez-Moyado I.F. et al. **Flanking DNA sequences determine DNA methylation maintenance in proliferation, cancer and aging.** bioRxiv (2026 preprint). https://doi.org/10.64898/2026.04.09.717557
7. Ushijima T. et al. **Fidelity of the Methylation Pattern and Its Variation in the Genome.** Genome Research 13, 868-874 (2003). https://genome.cshlp.org/content/13/5/868
8. Issa J.-P. **Aging and epigenetic drift: a vicious cycle.** Journal of Clinical Investigation 124, 24-29 (2014). https://www.jci.org/articles/view/69735
9. **Blood-based epigenetic instability linked to human aging and disease.** Nature Communications (2026). https://www.nature.com/articles/s41467-026-69430-z
10. **DNA methylation rates scale with maximum lifespan across mammals.** Nature Aging (2024 issue / online 2023). https://www.nature.com/articles/s43587-023-00535-6
11. **Meta-analysis of DNA methylation aging signatures in 17 human tissues.** Nature Aging (2026). https://doi.org/10.1038/s43587-026-01164-5
12. **Stochastic modeling of epigenetic memory.** npj Systems Biology and Applications (2026). https://www.nature.com/articles/s41540-026-00664-9
13. **Expanded ENCODE Registry of candidate cis-regulatory elements.** Nature (2026). https://www.nature.com/articles/s41586-025-09909-9
14. GSE225171 - single-cell mouse blood methylome/RNA source dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225171
15. GSE179847 - longitudinal primary human fibroblast methylation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179847
16. GSE179849 - fibroblast methylation + RNA super-series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179849
17. GSE73115 - ten-year repeated human blood methylation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73115
18. GSE56581 - MESA sorted monocytes: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56581
19. GSE40279 - Hannum blood dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279
20. GSE87571 - lifespan whole-blood methylation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571
21. GSE130748 - longitudinal leukocyte methylation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130748
22. ENCODE SCREEN / cCRE downloads: https://screen.wenglab.org/

