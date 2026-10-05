# Context-dependent methylation transitions and functional selection

**Current specification: revision 3, 2026-10-04.** The historical filename is retained to avoid breaking repository links. This supersedes the original protection-centered specification; its numerical APIs and earlier results remain historical comparators. The [review-directed power audit](../code/results/functional_power_audit.md), [scientific contract](research_plan.md#biological-prediction-contract-2026-10-03) and [novelty assessment](reading_list.md#revision-3-novelty-assessment-2026-10-04) govern interpretation.

## Decision and scope

Keep stochastic bounded states and intervention-dependent maintenance. Place independently measured functional consequences in cell fitness, not in intrinsic restoration or noise. Do not infer selection from the failure of protection: it is a separate, unvalidated mechanism. Constant rates remain a nested control within conditions where a time-dependent effect is unsupported; they are not universally biologically invalid.

The power audit does **not** establish that modest protection effects are absent. The primary physical-bound simulation detects available oracle MSE reductions of 3/5/10/20% with approximately 37/49/76/98% probability. Observed prediction deterioration is ordinary under score permutations. These are conditional regression detectability results, not power for a kinetic-rate effect or independent cultures. The tested functional predictor remains unsupported; permanent biological falsification would overstate the evidence.

The smallest primary model is a finite-state, time-inhomogeneous transition process coupled to birth/death and inherited clone labels. Its expected-population equation is a standard mutation-selection system. A bounded diffusion is optional and subordinate to the discrete observation architecture. No extra arbitrary latent layers are added to rescue previous fits.

## 1. Separate state, kinetics, fitness and identity

For selected alleles or explicitly coarse regulatory states, let

\[
M_{ci}(t)\in\{0,1\},\qquad s_c(t)\in\{0,1\}^{L}.
\]

Bulk beta is a sampled population methylation fraction, not a continuous state of an individual cell. A region-level functional edit cannot automatically be assigned to every CpG in that region. Diploid states, strand-specific maintenance and multiple CpGs require an explicit mapping if those mechanisms are tested.

Local molecular covariates \(Z_i\) include sequence, CpG context, chromatin, transcription and enzyme recruitment. Measured or separately calibrated maintenance exposure \(R(t)\) and division history supply nonstationarity. Initial implementation conditions on a frozen piecewise maintenance schedule; it does not infer a latent maintenance trajectory from fitness annotations.

A kinetic equilibrium \(b_i^{\mathrm{kin}}\), a fitness reference \(b_i^{\mathrm{fit}}\) and an independent identity reference \(b_i^{\mathrm{id}}\) are different quantities. None is automatically youthful, optimal or causally necessary.

## 2. Molecular transitions contain no functional score

Conditional on context and maintenance exposure,

\[
0\xrightarrow{a_i(t;Z_i,R)}1,\qquad
1\xrightarrow{d_i(t;Z_i,R)}0.
\]

For the minimal implementation, independently frozen context-specific rates are multiplied by epoch-specific gain/loss modifiers:

\[
a_i(t)=a_i^0g_{e(t)},\qquad d_i(t)=d_i^0\ell_{e(t)}.
\]

A globally shared exposure schedule creates common time dependence; uncertain common exposure can create dependence between loci. Given a known schedule, local jumps are conditionally independent. This is an explicit starting restriction, not a claim of unconditional CpG independence.

The joint generator \(Q(t)\) flips one coordinate at a time with the corresponding molecular rate and has row sums zero. Functional scores and fitness references are absent from the transition parameter type. No \(q_i\)-dependent error suppression, restoring coefficient or diffusion term is part of the primary model. Older protection/recovery coefficients remain available only to reproduce and compare historical hypotheses.

Gain/loss are effective transition hazards; they must not be mislabeled as uniquely identified DNMT1/DNMT3/TET activities. If methylation error is division-linked, its per-time hazard depends on division exposure. A functional effect on proliferation can therefore affect cumulative drift **indirectly**, even when per-division maintenance fidelity is independent of function. The minimal calendar-time implementation uses effective hazards and faithful inheritance at birth; it does not separately simulate replication-strand errors. Interpret claims conditional on division history, and require a division/inheritance extension before making those mechanistic predictions.

## 3. Independently specified state fitness determines representation

Let \(q_i\geq0\) denote magnitude of a validated functional effect and let \(h_i(s_i;b_i^{\mathrm{fit}})\) specify its **signed, state-dependent** fitness consequence. Define

\[
C(s)=\sum_i q_i h_i(s_i;b_i^{\mathrm{fit}}),\qquad
r_s(t)=r_0(t)-\gamma C(s),\quad\gamma\geq0.
\]

An equivalent implementation uses signed per-coordinate costs and a binary fitness reference. Positive costs penalize departures; negative costs favor them. The blanket assumption that every methylation deviation reduces proliferation is not justified, particularly in cancer. Tissue function, cell identity and competitive cell fitness can disagree.

Use independent editing/viability/growth assays to define costs, with measured edit success, direction and uncertainty. Do not infer those costs from methylation stability or clone expansion being predicted. A TET active/dead viability contrast without confirmed editing efficacy is only a proxy. Region effects, epistasis and state dependence must be declared before testing; an unrestricted fitted fitness landscape would violate the scientific contract.

Specify nonnegative state birth/death rates \(B_s,D_s\), with \(r_s=B_s-D_s\). The helper decomposition \(B_s=\nu+\max(r_s,0)\), \(D_s=\nu+\max(-r_s,0)\) is a mathematical choice with fixed turnover \(\nu\), not identification of proliferation versus killing. Equal net growth can have different extinction probabilities. Rate and turnover units are inverse time; arbitrary normalization of q rescales gamma.

Cells change state via Q, reproduce at B and disappear at D. Offspring retain clone labels and, in the minimal implementation, parental state; subsequent molecular transitions can change either lineage's state. Birth can expand an unfavorable or favorable state: functional selection is not represented solely by retrospective deletion of survivors.

## 4. Population dynamics and observations

For an exogenous schedule and no interactions, expected state counts satisfy

\[
\dot n_s=\sum_u n_uQ_{us}+r_sn_s.
\]

For their normalized composition \(p_s=n_s/\sum_u n_u\),

\[
\dot p_s=\sum_u p_uQ_{us}+p_s(r_s-\bar r),
\qquad \bar r=\sum_s p_sr_s.
\]

This describes normalized **expected counts**, not generally the expectation of a normalized finite stochastic population. Finite branching populations can become extinct; their observed conditional composition differs. Keep the birth/death simulator separate from the deterministic expectation solver.

For any state function f,

\[
\frac{d}{dt}\mathbb E_p[f]=\mathbb E_p[Qf]+\operatorname{Cov}_p(f,r).
\]

Thus bulk methylation can change through within-lineage state transitions or state-dependent reproduction/loss. The covariance identity is established mathematics, not project novelty. Common additive net growth cancels from composition but changes absolute counts. Common additive crowding and a common multiplicative brake on proliferation are different ecological assumptions and can produce different selection behavior.

Observations remain separate: binary allele calls, bisulfite counts with biological overdispersion, arrays with measurement error, and independently assigned clone labels with sampling probabilities. DNA or cell yield and clone counts constrain growth; methylation means alone generally cannot. Different transplant recipients are not longitudinal clone survival observations. Cell-state composition, genetics, treatment toxicity and collection bias can all imitate selection.

## 5. Derived bounded approximation

For N independent alleles without selection,

\[
\mu(x,t)=a(t)(1-x)-d(t)x,\qquad
v(x,t)=\frac{a(t)(1-x)+d(t)x}{N}.
\]

These transition-derived moments do not imply an arbitrary Jacobi noise term proportional to \(x(1-x)\). Functional scores do not enter either moment. The bulk diffusion approximation requires sampling/population assumptions; it is not a single-cell SDE. Selection adds representation and demographic stochasticity and cannot be absorbed into a fitted restoring coefficient with a biological interpretation.

## 6. Restrictive predictions and novelty gate

The broad architecture is **not novel**. Methylation-pattern transitions plus state-dependent selection appear in Stromberg (2012); reversible epigenetic switching with different growth rates and a Price decomposition appear in Qian (2014). Nonstationary methylation kinetics also predate this project. See the primary-source novelty assessment before testing any new empirical claim.

A narrower candidate contribution would transfer independently measured local methylation consequences into a held-out prediction of **clone expansion/loss conditional on molecular transitions**, rather than fit clone fitness to its own outcomes. For approximately stable states and a common environment,

\[
\log\frac{n_c(t)}{n_d(t)}-\log\frac{n_c(0)}{n_d(0)}
=\int_0^t(r_c-r_d)\,du.
\]

With transitions, predict the complete branching/state process, not a static final-state cost substituted into that identity. This mathematical relation is familiar; novelty, if any, is the external biological calibration and successful condition transfer. It is not established by combining familiar ingredients.

Compare context-only transitions, unrestricted empirical clone effects and external-cost selection with frozen exposure/cost rules. Exclude whole cultures/conditions, clones and relevant genomic regions. Keep absolute counts and within-clone methylation so pure repair and representation changes can disagree observably. Replication, division and editing efficiency must be measured or controlled. No dataset currently used jointly supplies those measurements and valid independent local functional weights; current clone data cannot test that transfer claim.

## 7. Identity loss is downstream

An independently validated identity distance may be

\[
D_{\mathrm{id}}(s)=\sum_i w_i h_{\mathrm{id},i}(s_i;b_i^{\mathrm{id}}),\qquad
\tau=\inf\{t:D_{\mathrm{id}}(s_c(t))>B\}.
\]

Its weights and boundary cannot be calibrated from survival or chronological age alone. Identity loss need not lower cell fitness; malignant expansion is an obvious possible mismatch. Death is a competing event, and descendants can carry distinct first-passage histories. The legacy single-cell first-passage simulator remains a numerical tool; a tissue-level validated identity threshold is not claimed.

## 8. Implementation and evaluation order

The revision-3 public Rust API separates context-only sites, exposure epochs, state fitness and birth/death. A positivity-preserving expected-count solver and exact event simulation for piecewise-constant exposures provide complementary checks. Exact state enumeration is restricted to small selected state systems; it is not a whole-genome inference engine. The existing historical APIs and data analyses are retained for reproduction, with protection labeled a comparator.

Complete power interpretation and primary-source novelty assessment before revised-model tests. Analytic reduction, covariance decomposition, q/transition separation, exposure-boundary handling, growth and extinction checks verify the software. They are not biological confirmation. No unrestricted new fit or full biological-suite rerun is justified until a dataset can evaluate the narrower external-fitness transfer claim.

Completed numerical verification: `make check` passes formatting, strict Clippy and 33 tests. Six new population checks cover exact time-dependent CTMC limits, the Price decomposition, fitness/transition separation, second-order solver convergence, stochastic branching means/clone inheritance, extinction and resource limits. A deliberately synthetic zero-transition example changes bulk methylation from 0.5 to 0.818 over five time units solely through a 0.3 net-growth difference. This demonstrates the bookkeeping distinction, not an estimated biological effect or novel experimental prediction.
