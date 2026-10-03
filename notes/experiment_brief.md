# From functional constraint to a tractable experiment

## Decision and collaborator-facing question

The equations are revisable. The objective is a discriminating biological prediction: functionally consequential methylation states should exhibit a predictable maintenance/restoration response that cannot be explained entirely by local sequence, transcription/chromatin context, or differential cell expansion. A better methylation fit is not the endpoint.

**Recommended first laboratory question:** does localized demethylation of an expressed, externally essential gene's candidate gene-body region change its transcriptional output, and does subsequent restoration accompany reversal of that change? Start with four regions and targeted assays. This is an exploratory functional feasibility test, not a request for a mouse aging study, genome-wide single-cell profiling, or validation of a universal identity threshold.

The current public-data evidence does not justify a definitive experiment asserting that essential loci recover faster. The initial candidate assay should establish whether the marks themselves have a measurable function. If it fails after successful editing, these candidates cannot support that interpretation regardless of their recovery curves.

## What the new public-data pilot actually shows

We selected [GSE51810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51810) and [GSE51811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51811), an HCT116 pulse/washout study with methylation, expression, and DNMT knockout comparisons. The authors already reported DNMT3B-dependent rebound; this supplies a positive control rather than a novelty claim. The system is cancer-cell culture, not normal colon or aging. [Yang et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4224113/)

Functional membership comes from 684 core-essential and 927 reference-nonessential genes in the authors' [BAGEL repository](https://github.com/hart-lab/bagel/tree/53388adbb4fb0931e5c9dda135502be19e4555f0), independent of these methylation outcomes. This is gene knockout fitness, not evidence that any particular methylation mark is important. An unexpressed nonessential gene is an inadequate control for an expressed essential gene.

Eligibility: a single annotated gene, BODY-only probe, WT baseline beta >=0.7 and induced baseline-minus-day5 loss >=0.15. Match chromosome and CGI class exactly, and use baseline beta, induced loss, baseline expression, sequence rank and probe-flank CpG density with fixed calipers. No later DNA/RNA outcome chooses the match. Use one probe pair per distinct gene pair, without control-gene reuse. The primary endpoint was specified before analyzing recovery: day42 regained methylation divided by induced loss; values outside [0,1] remain observable.

| Analysis | Matched gene pairs | Mean essential-minus-control day42 recovery fraction |
|---|---:|---:|
| Original matching, without RNA detection filtering | 122 | +0.312 |
| Same matching after baseline RNA detection QC | 7 | +0.475 |
| Detection QC plus expression difference <=0.5 log2 units | 2 | +0.079 |
| Detection QC plus all original calipers halved | 1 | +0.519 |

RNA detection QC requires p<=0.01 in both WT baseline arrays, with only detected probes contributing to gene expression. This and tighter matching were added after inspecting the initial result and are exploratory measurement/confounding audits. There is severe lack of overlap: after detection QC, 453 essential genes but only ten reference-nonessential genes have eligible loci. The narrow two-pair contrast is inconsistent across pairs: CPSF4/DTX2 is -0.105, PTPN23/DRD3 is +0.264. We retain both, including the negative pair.

Only one detection-QC pair has an absolute RNA change >=0.25 log2 units at day5 in both genes, so there is no supported functional-expression recovery comparison. There is one methylation array per genotype/time. Resampling loci measures heterogeneity and cannot establish biological reproducibility. Genomic sequence rank and probe-flank density do not replace chromatin accessibility, histone marks, replication timing, or genome-wide density. Cell growth changes and clone composition remain uncontrolled. Methylation detection p-values and cross-reactive/SNP probe exclusions are unavailable in the current matrix analysis. Some GEO genotype/treatment characteristic labels are inconsistent; sample titles were audited and used instead.

An initial coding error excluded genes based on later RNA values. That run's 48-pair result is superseded. The corrected routine checks baseline RNA eligibility only; negative later intensities become unavailable log2 measurements, never DNA eligibility exclusions. The corrected original result and all sensitivity results above are retained. See [the data audit](../code/results/recovery_data_audit.md), [the detection-QC audit](../code/results/recovery_detected_data_audit.md), and [the sensitivity table](../code/results/recovery_detected_matching_sensitivity.csv). Reproduce with `make recovery` after building the frozen human feature table.

## Candidate regions for a small functional pilot

These are the two pairs surviving detection QC and the 0.5-log2 expression caliper. They are selected by baseline/context matching, not favorable recovery or RNA outcomes. Array CpGs are anchors for a locally assayed region; editing one CpG is not assumed to perturb a functional element.

| Pair | External category | Gene | CpG anchor | Frozen hg19 coordinate |
|---|---|---|---|---|
| A | Core essential | CPSF4 | cg01611679 | chr7:99051465 |
| A | Reference nonessential | DTX2 | cg08801436 | chr7:76111937 |
| B | Core essential | PTPN23 | cg19986876 | chr3:47423585 |
| B | Reference nonessential | DRD3 | cg26358350 | chr3:113872166 |

The categories concern the published reference screens, not verified HCT116 essentiality of each local mark. Reconfirm expression and regional methylation in the actual cell stock. Before reagent ordering, align each anchor/amplicon uniquely, inspect transcript and alternative-promoter annotations, check cell-line variants/copy number and guide compatibility, and record any coordinate conversion. Do not design guides from an unverified hg19-to-hg38 conversion. No guide sequences or primers have yet been validated. If a region overlaps an alternative promoter in the stock's relevant transcript annotation, analyze it as such rather than retaining a gene-body interpretation.

## Staged laboratory design

**Stage A: establish that the local marks affect gene output.** Use a lab already experienced with targeted methylation editing. Deliver a transient/inducible dCas9-TET catalytic editor separately to each of the four regions, paired with a catalytically inactive editor using the same guides. Add active and inactive editors with nontargeting guides. Targeted editing methods have primary experimental support, but these four regions have not been demonstrated editable. [Liu et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27662091/), [Morita et al., 2016](https://www.nature.com/articles/nbt.3658)

Start with three independent delivery/culture batches, blocking all conditions within batch. Four targets x two editor arms x three batches is 24 cultures; two nontargeting arms add six, for 30 cultures. An established positive-control editing locus in the collaborator's own system adds six cultures. At two destructive sampling times this is 60 candidate/control DNA-RNA sample pairs, or 72 with the positive control. Technical wells/reads are not independent replicates. This is a feasibility sample size, not a powered confirmatory study.

Collect a pretreatment stock reference and two matched active/dead-editor endpoints: one during verified editing and one approximately 48 hours later. Define the actual editing window from editor-expression measurements before the biological test. Measure regional 5mC across multiple CpGs, the intended transcript and canonical splice junctions, neighboring/alternative transcripts where annotation warrants them, and viable cell counts/divisions. Standard bisulfite alone cannot distinguish residual 5mC from TET-generated 5hmC: use a paired assay such as oxidative bisulfite with appropriate conversion controls. Measure editor persistence; removal of inducer alone does not establish removal of catalytic activity.

Primary functional hypothesis for this panel: localized demethylation changes the intended gene's transcriptional output beyond dead-editor and nontargeting controls, with a stronger consequence at the essential candidates. Do not assume promoter-like activation for gene-body demethylation. Specify transcript assays and the predicted direction for each region after transcript annotation review, before experiments. Where no defensible sign is available, register a two-sided functional screen; it does not confirm a directional mechanistic prediction.

Provisional assay gates, to be finalized using control-assay variance before candidate outcomes: at least 0.15 absolute regional 5mC reduction, a reproducible transcript effect of at least 20%, and concordant direction across independent batches. The 20% threshold is an operational screen, not a known biological identity boundary. Confirm any hit with a second independently targeting guide set. A targeted remethylation rescue can strengthen methylation-specific interpretation, provided a matched inactive writer control and regional/promoter-spreading measurements exclude simple occupancy effects. TET-associated changes without such confirmation are evidence for an editing intervention, not clean proof of methylation causality. Keep failed edits separate from successful edits with no functional effect.

**Stage B: test restoration only for functionally supported regions.** Use new cultures after Stage A, freeze region importance and endpoint definitions, stop the editor and verify its clearance, then sample at recovery baseline and days 2/4/7. Compare recovery at comparable induced methylation loss with context-matched regions. RNA return toward its independent pre-edit reference is an additional outcome, not a way to choose methylation weights. A shared early-time recovery model must predict withheld later measurements without condition-specific refitting; compare carry-forward, linear and context-based empirical predictors.

Track absolute viable counts, division histories, and lineage abundance. Bulk rebound alone remains ambiguous. Before claiming protection over selection, obtain clone-conditioned methylation or a sufficiently informative lineage-linked observation design and simulate its recoverability using measured variation. Barcodes plus bulk methylation do not by themselves establish within-clone restoration. A short window reduces opportunities for selection but does not eliminate it. If those measurements are unavailable, keep the conclusion at local function/recovery, not a selection-versus-protection mechanism.

A DNMT3B perturbation with matched viability/division controls and rescue is a subsequent mechanism test if Stage B supports restoration. Constitutive knockout comparisons can contain adaptation and baseline differences; they cannot substitute for a timed, reversible intervention. No DNMT1-only recovery mechanism is assumed.

## Revised mathematical commitments

Use a cell-state transition model plus cell/clone abundance, keeping observations separate from biology. For an effective within-state generator Q and net growth rate r, a normalized state distribution satisfies

`dp_s/dt = sum_u p_u Q_us + p_s (r_s - mean_p(r))`.

Consequently, change in mean methylation contains both a transition contribution and `Cov_p(methylation, r)`. This is why a bulk restoring curve cannot identify intrinsic restoration. If division changes inherited methylation states, that inheritance must be included in the transition accounting, not silently counted again as growth. A molecular refinement can use unmethylated/hemimethylated/fully methylated dyads and a measured division inheritance kernel; binary array beta cannot identify those states on its own.

A minimal descriptive recovery curve is `x_i(t)=b_i-(b_i-x_i(0))*exp(-k_i*t)`, with b independently defined and k predicted from context plus a frozen functional annotation. It does not imply that DNMT1 causes recovery, that there is a fixed equilibrium in every condition, or that a cell senses gene essentiality. Separate intrinsic loss and restoration, active demethylation, and net expansion. Add rates/states only when a discriminating observable requires them. Days and cell divisions are different clocks; the proposed assays measure both. These are revised theory commitments, not newly validated or fully fitted model components.

## Go/no-go interpretation

Successful local editing without a functional consequence rejects these regions as the needed functional link; retain them only as methylation kinetics controls. A transcript effect without reproducible restoration supports local epigenetic function but not preferential repair. Bulk recovery explained by changes in lineage abundance supports representation, not protection. Reproducible within-lineage recovery that transfers to withheld conditions and predicts functional return would support the distinctive hypothesis. A stronger fit after adding arbitrary latent states would not.

The immediate justified collaborator request is Stage A, with four anchor regions and these controls, not a claim that the model is validated. Guide/assay design and a technical editing screen remain necessary before execution. No collaborator has been contacted and no experiment has been run.
