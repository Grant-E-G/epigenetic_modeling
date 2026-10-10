# Dog Aging Project: attempted falsification of functional protection

2026-10-09. **No robust support for the distinctive protection claim, and no defensible rejection of the full model.** The essential-gene proxy comparison is unstable across sequence matching and coverage choices, its intervals admit both harm and substantial protection, and the promoter arm fails its support gate. The acquired data cannot adjudicate independently calibrated clone selection. Keep the previous decision to defer wet experiments; no target panel is earned.

This is an executed data analysis, not an access proposal. Fifteen source files are downloaded/checksum-frozen, all processed counts and identifiers are validated, and executable Rust matching, evaluation, resampling and conditional-power analyses are supplied. The [protocol](../../notes/dog_aging_protocol.md) was written before matrix outcomes. A necessary sequence-context correction occurred **after** the initial contrast; both versions are retained. This is retrospective analysis, not externally registered prospective validation.

## What can be falsified here

The historical protection claim predicts that independently vital regions resist methylation disruption beyond molecular context. Its operational proxy here is lower within-dog methylation change in regions of core-essential genes, compared with reference-nonessential genes, after matching. The practical hypothesis fixed before outcomes is at least a 10% reduction in noise-adjusted mean squared annual change. The 10% threshold is an explicit decision threshold, not an effect derived from the equations. Also report compatibility with 5% and 20% protection.

Gene-knockout essentiality is **not** local methylation necessity. Human cell-culture essentiality transferred through verified orthology does not establish that a particular dog PBMC gene-body mark is functionally important. This distinction was fixed before outcomes and restricts either a positive or negative interpretation. A null proxy test cannot be renamed a definitive test of vital methylation marks; a positive proxy test cannot validate them either.

Revision 3 instead puts independently measured local consequences in signed state fitness, explicitly excluding importance from molecular transitions. Its distinctive candidate prediction concerns future clone expansion/loss conditional on molecular state transitions. The acquired dog resources contain no jointly usable local-effect/edit-efficacy calibration, longitudinal clone-linked methylation, clone abundance/absolute counts or division history. **That claim fails the data-compatibility gate.** It is untested here, not confirmed and not rejected. Failure of direct protection supplies no evidence for selection.

## Acquired data and independent biological units

The paper is [Mariner et al., Science, DOI 10.1126/science.aeb2986](https://doi.org/10.1126/science.aeb2986), published 2026-10-08. Primary data are [GSE306794](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306794), its [processed-count subseries GSE306792](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306792/suppl/) and [sample metadata GSE306793](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306793/suppl/). The [authors' repository](https://github.com/smacklab/DAP_DNAm/tree/8fbf9084f633a2f14717b6b18f2dc321d99f14dc) supplies preparation-batch metadata and a count-model workflow. We did not obtain the inaccessible Science supplementary PDF and do not claim to have reproduced its entire analysis.

| Audit item | Observed |
|---|---:|
| Samples / unique dogs | 1,640 / 894 |
| Dogs with 1 / 2 / 3 / 4 visits | 356 / 339 / 190 / 9 |
| Dogs with repeated ages | 538 |
| Adult first-to-next pairs, first age >=2, interval 0.5–2.5 years | 490 |
| Training / withheld evaluation dogs | 390 / 100 |
| Unsuffixed CpG regions | 184,152 |
| Excluded labeled promoter/other aggregate rows | 20,338 |
| Validated methylated/total count pairs | 335,363,600 |
| Reference CpG concordance under one-based / zero-based starts | 2,000/2,000 / 0/2,000 |

Matrix entries are integer methylated counts and total coverage, not percentages or array betas. All count pairs satisfy 0 <= methylated <= total, row/column identifiers match, regional identifiers are unique, and all matrix columns match the sample metadata and author dog IDs. Zero coverage is missing. One-based inclusive region coordinates are converted to zero-based half-open for canFam4 reference joins. Sex chromosomes and multiply assigned genes are excluded. The matrix concatenates region and labeled aggregate tables; the latter are excluded to avoid reusing the same underlying CpGs as extra observations.

The GEO metadata repeats the same BAM filename in all rows. We therefore join by the 1,640 unique LID/PID identifiers, never by filename. The public author metadata includes preparation dates, sex, age, size/breed-status category, predicted height, reproductive status and heterozygosity. It does not provide individual breed labels, a genomic relationship matrix or usable longitudinal blood-cell fractions in the acquired files. The authors' code references a GRM that is not in the published repository; dog ID splitting does not guarantee independence of relatives or a breed-held-out test. Broader Terra resources may add covariates, but we did not access an account or sign a data-use agreement.

## Frozen external annotation and matching

Use the previously frozen [BAGEL CEGv2/NEGv1 lists](https://github.com/hart-lab/bagel/tree/53388adbb4fb0931e5c9dda135502be19e4555f0), [NCBI human–dog orthology](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz), and exact-assembly [canFam4 RefSeq/CpG islands](https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/database/) plus [reference sequence](https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.2bit). Mutable reference inputs are frozen by hashes in [dog_aging_sources.csv](dog_aging_sources.csv).

Of 684 essential and 927 control list entries, 666 essential and 546 control genes have unique human-to-dog/dog-to-human relations in the acquired orthology table. Fifteen essential/377 control entries lack ortholog relations and three entries in each list lack a human gene ID. No symbol-based substitution is used. The lists define the proxy independently of dog age effects, LINE1 hits, clock coefficients or disease outcomes.

Prespecified promoter arm: wholly inside a +/-1-kb RefSeq TSS window and training baseline beta <=0.3. Gene-body arm: wholly inside a unique gene span, outside every promoter window, and training baseline beta >=0.7. Exact matching is on chromosome, arm and CpG-island overlap; calipers are baseline beta 0.05, mean log coverage 0.5, length ratio <=2, GC fraction 0.05 and CpG fraction 0.02. At least 20 training first visits must support each regional baseline. A deterministic greedy closest-distance match permits one region pair per distinct gene pair, without gene reuse.

The initial implementation used sequence only inside measured intervals. Single-CpG regions then have GC=1 and CpG density=0, because the G falls outside the one-base interval. That is inadequate local context, even though reference coordinates are correct. We preserved [the original estimate and matches](dog_aging_interval_context/dog_aging_summary.csv), then used a fixed 250-bp flank on each side, retaining every other rule. The corrected results below are a transparently post-result context audit. This correction reverses the point estimate; neither version is selected as evidence for a preferred mechanism.

For the corrected primary gene-body matching, standardized mean differences are 0.133 for baseline beta, 0.083 for log coverage, -0.068 for interval length, 0.036 for GC, 0.054 for CpG density and 0.069 for the number of covered training dogs. Exact chromosome/island balance is enforced. These covariates are reasonably balanced on average, but unmeasured transcription/chromatin, genetic ancestry and cell composition remain. Essential genes may be expressed while reference-nonessential controls are inactive in PBMCs; sequence balance cannot remove that possibility. [Full balance](dog_aging_balance.csv).

## Endpoint, uncertainty and results

For each dog/region pair, calculate

`adjusted_change = [(p_next-p_first)^2 - p_first(1-p_first)/(n_first-1) - p_next(1-p_next)/(n_next-1)] / elapsed_years^2`.

The subtraction is unbiased under independent binomial calls; retain negative adjusted observations. Pooled regional counts may include multiple CpGs/read and biological/technical dependence. Thus this subtraction is an idealized measurement correction, not an established noise model for the source assay. Raw squared change is a prespecified secondary endpoint and also does not isolate intrinsic drift.

Average region-pair outcomes equally within dogs and dogs equally across the available paired records. Report protection fraction `1 - essential_mean/control_mean`, so positive values indicate less essential-associated change and negative values indicate more. Missing regions are not assigned zero and no endpoint is calculated from a nonpositive control mean.

Dog resampling keeps complete observed dog summaries together. Gene resampling and crossed dog/gene resampling assess sensitivity to genomic support; 10,000 draws each, seed 20261009. The table uses **crossed** 95% descriptive intervals, which capture more uncertainty than treating genes or dogs alone. The planned two-arm multiplicity-aware interval is 97.5% two-sided; it is also retained in [dog_aging_summary.csv](dog_aging_summary.csv). These are conditional resampling intervals, not causal or fully adjusted biological confidence guarantees.

| Flank-context gene-body analysis | Matched pairs | Observed dogs | Estimated protection | Crossed 95% interval |
|---|---:|---:|---:|---:|
| >=20 reads, first age >=2 | 75 | 96 | +32.6% | -62.0% to +59.3% |
| >=10 reads, first age >=2 | 72 | 100 | -2.1% | -81.1% to +40.7% |
| >=30 reads, first age >=2 | 67 | 87 | -12.3% | -237.9% to +52.8% |
| >=20 reads, first age >=1 | 75 | 106 | +26.7% | -75.3% to +56.5% |
| >=20 reads, first age >=4 | 68 | 64 | +13.9% | -111.3% to +57.4% |

The primary multiplicity-aware crossed interval is **-81.5% to +63.2%**. The raw-squared primary estimate is +26.5%, crossed 95% interval -34.4% to +51.4% (97.5%: -46.2% to +55.5%). The original interval-context gene-body estimate was **-41.2%**, crossed 95% interval -197.8% to +34.5%. Large negative lower bounds are possible because the endpoint ratio is unbounded below and adjusted control denominators can approach zero under resampling. All primary gene-body bootstrap draws have positive denominators; some insufficient-support promoter draws do not.

**No primary or prespecified sensitivity interval excludes 10% protection, 20% protection or no protection.** None earns a robust positive protection claim. Changing the measurement gate changes the matched locus set as well as precision; these are not independent replications. The original negative point estimate and corrected positive estimate cannot be averaged into an independent meta-analysis.

The corrected promoter arm has only **three** matched gene pairs (two at >=10 reads), far below the prespecified 20-pair gate. Its positive point estimates cannot compensate for failed genomic support. Keep the failed arm, do not nominate a replacement mask or infer general promoter protection from a handful of genes.

## Why the large cohort still gives a weak test

The primary gene-body comparison has 942 jointly covered dog/pair observations among 100 eligible evaluation dogs x 75 matched pairs: **12.56%** of the potential records. Only 65 of the 75 matched pairs are observed at all in evaluation dogs, with a median of **eight dogs per observed pair**; the median observed dog has seven pairs. Thus 96 observed dogs does not mean a complete 96 x 75 biological design. Coverage availability could depend on genotype, library or region properties. No missing-at-random assumption has been validated. [Measurement audit](dog_aging_measurement_audit.csv).

All known first/next preparation batches differ, and three evaluation dogs have missing preparation-date information for at least one visit. Paired essential/control observations share the dog's batch transition, but regional batch effects need not cancel after squaring. Elapsed time and technical changes are not separately identified. The data do not measure clone representation or within-cell repair.

The conditional power simulation centers observed dog-level paired differences, preserves their residual magnitudes, randomizes signs, then adds a fixed 0/5/10/20% effect relative to the observed mean control endpoint. It calibrates a one-sided 98.75th-percentile null statistic with 10,000 draws and evaluates each effect with 10,000 new draws (seed 20261010). This is **planning under a conditional independent-dog endpoint generator**, not power for the generative model or externally validated local-effect annotations. It also conditions on the observed gene mask and missingness rather than simulating new genome-wide sampling. Its test statistic is not identical to the percentile-bootstrap ratio decision.

| Injected gene-body protection | Conditional detection probability | Monte Carlo SE |
|---|---:|---:|
| 0% | 1.19% | 0.11 percentage points |
| 5% | 2.35% | 0.15 percentage points |
| 10% | 4.33% | 0.20 percentage points |
| 20% | 10.36% | 0.30 percentage points |

The raw-squared conditional generator detects 10%/20% effects with 5.95%/15.47% probability. These low detection rates make a null or unstable estimate weak evidence against modest protection. They cannot be compared numerically with the earlier recovery power audit, whose effect scale was oracle forecast-MSE improvement. [All conditional simulations](dog_aging_power.csv).

Sex/size summaries are retained as [descriptive subgroup diagnostics](dog_aging_subgroups.csv); none chooses the primary population or rescues the claim. The acquired covariates do not justify a breed/kinship-controlled transfer claim.

## Falsifiability assessment and decision

1. **Operational essential-gene protection >=10%:** measured, but unresolved. The intervals do not rule it out; the power audit is weak. The same holds for 5% and 20%. A sign-changing point estimate is evidence of fragility, not statistical proof that the biological effect is zero.
2. **Vital local methylation marks are intrinsically protected:** the assay lacks independently validated local necessity and cell-conditioned state observations. This dataset does not establish or reject that mechanism.
3. **External local fitness predicts clone selection (revision 3):** incompatible measurements. Untested, with no support transferred from the stability comparison.
4. **A robust practical signal from this vital-gene proxy:** not earned. Do not present the +32.6% corrected estimate as validation; do not present the -41.2% original estimate as definitive rejection.

The broad phrase “vital regions are protected” is insufficiently restrictive without a fixed local importance definition, reference state, effect size/direction and observable condition. If the model allows an arbitrarily small effect or unspecified weights, a finite noisy bulk dataset cannot falsify that unrestricted family. Revision 3 also nests zero functional selection. This is a limitation of the model's current empirical specification, not a reason to infer a hidden protective effect or fit extra parameters.

The observable population derivative combines molecular transitions with selection covariance; differing cell composition can generate either sign of bulk change. Without external local-effect calibration and lineage/count observations, assigning the residual to repair or selection would evade the scientific contract. No condition-specific fitness landscape, alternate importance mask or learned identity threshold is fitted here.

**Stopping decision:** retain these results as a failed attempt to obtain decisive evidence. No distinctive biological mechanism is validated, and the available data do not justify whole-model invalidation. Continue to treat protection/selection as unsupported hypotheses and defer wet experiments. A future test must independently fix the local consequence and quantitative prediction before new outcomes; increasing nominal sample size alone does not fix the present functional-definition and observation gaps.

## Reproduction and software verification

See [code instructions](../README.md#dog-aging-project-falsification-audit). `make dog-download` acquires only missing checksum-frozen sources. `make dog-validation` verifies sources, decodes author R metadata, executes all Rust settings and regenerates the power/balance/missingness/subgroup outputs offline. No new dependencies are introduced; matrices/reference sequence and full dog/pair caches remain ignored under `data/`. The initial interval-context analysis is separately reproducible and retained.

All 15 input hashes pass. Rust formatting/strict Clippy and 37 tests pass, including exact-binomial correction, overlapping-interval lookup and an explicit guard against choosing matches from future counts or treating missing coverage as zero; Python Black/mypy pass. Those checks verify implementation behavior, not the biological model. No outreach, data-use agreement, spending or wet experiment has occurred.

A complete second offline `make dog-validation` reproduces the summary, power, balance and measurement-audit files byte-for-byte (SHA-256 comparison). No additional statistical settings were selected after that reproduction.
