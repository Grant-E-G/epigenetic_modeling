# From functional constraint to a tractable experiment

## Current decision: public data and broad assays first

The user rejects the value of a $15k consumables pilot at two candidate regions and directs us to assume C5R has no mammalian workflow. The immediate priority is the broader public-data route and cost comparison below. The earlier local-editing budgets are retained as reference, not the current spending recommendation. No mammalian workflow should be built at C5R for this pilot.

## Biological question

The equations are revisable. The objective is a discriminating biological prediction: functionally consequential methylation states should exhibit a predictable maintenance/restoration response that cannot be explained entirely by local sequence, transcription/chromatin context, or differential cell expansion. A better methylation fit is not the endpoint.

**Previously considered local-editing question:** does localized demethylation of an expressed, externally essential gene's candidate gene-body region change its transcriptional output, and does subsequent restoration accompany reversal of that change? The full panel has four regions; the costed smaller first step has one candidate pair and retains the same controls. This is an exploratory functional feasibility test, not a request for a mouse aging study, genome-wide single-cell profiling, or validation of a universal identity threshold.

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

The earlier collaborator proposal was Stage A at four anchor regions. The current broader-data-first direction below supersedes it as the initial experiment. Guide/assay design and a technical editing screen remain necessary before execution. No collaborator has been contacted and no experiment has been run.

## Cost and C5R feasibility assessment (2026-10-03)

The user wants a defensible cost and workload estimate before approaching personal acquaintances at C5R. The earlier $2k–5k / half scientist-day / one operator-day caps were arbitrary negotiation targets, not forecasts, and are withdrawn. The earlier $25k–55k service allowances are superseded by the counted budget below. No inquiry, purchase or experiment has occurred.

**Historical local-editing budget, superseded as the first spending recommendation:** reserve **$15,000 cash plus about 90 lab staff hours** for a controlled two-region Stage A pilot. At an assumed $100/hour staff value, this is about **$24,000 total**, so $25,000 is a reasonable initial funding envelope if existing equipment, usable editors and assay qualification suffice. For the full four-region panel, reserve **$20,000 cash plus about 110 lab hours**, approximately $31,000 at the same rate. These are conditional planning estimates, not C5R quotes or hard ceilings. Budget 4–6 elapsed weeks, potentially 6–8 with procurement/qualification delays; this is not a one-day favor.

### Scope and actual workload

Start with CPSF4/DTX2, the first existing candidate pair, subject to baseline assay qualification. Its public recovery contrast was negative; choosing it does not select the favorable outcome. One pair can establish local functional feasibility, not a general essentiality effect. The full panel adds PTPN23/DRD3. Neither initial scope includes Stage B recovery/clone tracking, a DNMT3B intervention, remethylation rescue, or independent-guide confirmation of a hit. Those require subsequent funding. The initial panel can reject these candidate regions as a functional link, but cannot validate the distinctive stochastic restoration mechanism.

| Workload | Two regions | Four regions |
|---|---:|---:|
| Candidate regions × active/dead editor × 3 independent batches | 12 culture arms | 24 culture arms |
| Active/dead nontargeting controls × 3 batches | 6 arms | 6 arms |
| Active/dead qualified positive control × 3 batches | 6 arms | 6 arms |
| Total arms × 2 destructive endpoints | **48 separate culture harvests** | **72 separate culture harvests** |
| Pretreatment references + qualification harvest allowance | 3 + 12 | 3 + 12 |
| Planned DNA/RNA paired extractions, before repeats | **63** | **87** |
| Two separate modification conversions per DNA sample | **126** | **174** |
| Purchased conversion capacity | 96 per chemistry = 192 total | 96 per chemistry = 192 total |
| Target + flanking amplicons at each candidate and positive-control region | 6 per converted sample | 10 per converted sample |
| Singleplex locus PCRs + one indexing PCR per converted sample | 756 + 126 = **882** | 1,740 + 174 = **1,914** |
| Approximate RNA reaction purchase, including qualification/QC | **1,000** | **1,500** |

Culture arms are conditions, not reusable harvested wells. Use matched wells for each endpoint; allocate adequate cell yield for paired DNA/RNA and both conversions. A 24-well format is a budgeting assumption, not a validated plate layout. The three batches are independent culture/delivery starts, not technical replicates. Qualification can consume its allocation and end the project without candidate execution if delivery or the assay fails.

RNA costing uses six sample-relevant assays in duplicate: target abundance, a transcript-junction assay, a neighboring transcript, editor expression and two reference genes. Nontargeting and baseline samples need the wider candidate panel; remaining capacity covers these extra reactions, dilution standards, no-template/no-RT controls and limited repeats. The two-region estimate starts at 63 × 6 × 2 = 756 reactions; the four-region estimate starts at 87 × 6 × 2 = 1,044. Assay identity and reference stability still require validation. No RNA-seq is included.

The DNA budget deliberately measures the full small panel in every sample, including flanks. It assumes singleplex locus amplification, followed by sample/chemistry-specific indexing and pooling, rather than free, already validated multiplex PCR. Aim initially for roughly 2,000 usable reads per sample/chemistry/amplicon: about 1.5M reads for the smaller panel and 3.5M for the larger before QC overhead. Sequence a compatible prepared pool with low-diversity compensation; a 5M-class run may suffice for the smaller panel, while the full panel may require more capacity or a rerun. Read count is not biological replication or independent starting-molecule count.

### Cash bill of materials

USD; full purchased kit/bottle capacity is charged, not just the aliquots consumed. Existing culture, PCR/qPCR, quantification and library-QC equipment are assumed. New editor constructs and substantial new method development are addressed separately below. Ranges marked allowance are estimates tied to the listed workload, not supplier quotes.

| Item / quantity | Two regions | Four regions | Basis |
|---|---:|---:|---|
| One HCT116 stock | $577 | $577 | [ATCC CCL-247](https://www.atcc.org/products/ccl-247), displayed list price; subtract if a suitable stock is already available |
| Culture medium/serum, plates/flasks, cell-stock identity/mycoplasma QC | $800–1,400 | $1,000–1,700 | Allowance for stock expansion, qualification and three experimental batches; typically 1–2 L medium and shared serum stock |
| One 0.75 mL delivery-reagent bottle | $556 | $556 | [Thermo Fisher L3000008](https://www.thermofisher.com/order/catalog/product/kr/en/L3000008), displayed USD reference price; local checkout may differ |
| Guide cloning, plasmid preparation and sequence QC | $500–1,000 | $800–1,500 | Allowance: two guides per target/positive locus plus one nontargeting guide = 7 / 11 guide constructs; includes DNA prep of existing editor pair |
| RNA/DNA primers and sample-index oligos | $700–1,200 | $1,000–1,600 | Allowance for 6 / 10 DNA primer pairs, RNA assay panel, controls and 126 / 174 indexed preparations; custom sequence quote required |
| Two paired DNA/RNA extraction kits, 50 samples each | $960 | $960 | 2 × $480, [Zymo D7003](https://www.zymoresearch.com/products/quick-dna-rna-miniprep-plus-kit); 100 extraction-pair capacity |
| Total-modification conversion module, 96 reactions | $696 | $696 | [NEB E8020L](https://www.neb.com/en-us/products/e8020-nebnext-enzymatic-methyl-seq-v2-conversion-module), displayed USD list price |
| 5hmC-specific conversion module, 96 reactions | $467 | $467 | [NEB E3365L](https://www.neb.com/en-gb/products/e3365nebnext-enzymatic-methyl-seq-5hmc-conversion-module), displayed USD list price |
| Uracil-compatible locus/index PCR, cleanup beads and reaction consumables | $1,200–2,000 | $2,200–3,400 | Allowance for ~1,000 / ~2,200 PCR reactions including QC, plus conversion/PCR cleanup; additional library chemistry if needed is not established |
| One-step RNA assay reagent | $1,158 | $1,827 | [NEB E3005](https://www.neb.com/en/luna/~/link.aspx?_id=ECCE5D577424448D983A70BA9BDA8613&_z=z): 1,000 reactions $1,158; additional 500 $669 |
| DNA/RNA quantification and prepared-library QC | $300–600 | $400–800 | Consumables allowance using existing instruments; no external per-sample core fee assumed |
| Prepared-pool sequencing | $1,000–2,000 | $1,000–2,000 | Allowance requiring external eligibility/run quote; see price anchor below; not a per-sample library service |
| Shipping and tax allowance | $500–900 | $600–1,100 | Multiple reagent orders and sequencing transfer |
| **Subtotal** | **$9,414–13,514** | **$12,083–17,183** | Before contingency |
| **With 25% contingency** | **$11,768–16,893** | **$15,104–21,479** | Round to **$12k–17k / $15k–22k** |

Prices observed on 2026-10-03. Manufacturer pages sometimes localize currency; only explicitly displayed dollar values above were used, not pounds/euros relabeled as USD. Account, country and eligibility differences still need checkout confirmation. The precise subtotal is arithmetic, not precision about execution cost. The contingency is a planning reserve, not a probability statement, and does not fund unlimited assay development.

[UT Austin GSAF's July 2026 prices](https://utexas.atlassian.net/wiki/spaces/GSAF/pages/1416560685/GSAF%2BPrice%2BList) list internal MiSeq i100 5M/600-cycle sequencing at $1,026.84 and 25M/600-cycle at $1,584.54. These anchor run scale, not an external commercial entitlement. The $1k–2k allowance depends on an eligible sequencing provider accepting our properly prepared pool; externally priced sequencing can exceed it. Library preparation is already allocated to our lab work above and is not included in these run prices.

### Material assumptions that can change the estimate

**Current modification chemistry needs qualification.** [Tecan's phase-out notice](https://www.tecan.com/phaseouts) lists the TrueMethyl oxBS Module 0414-32 as shipped until December 2024, with support ending March 2025. The prior proposal cannot assume a newly purchasable turnkey TrueMethyl kit. The costed alternative pairs total-modification EM-seq with 5hmC-specific E5hmC-seq and estimates 5mC from their difference. Validate conversion controls, subtraction uncertainty and amplification bias; the chemistries have different template conversion patterns. They are separate measurements, not simultaneous marks on the same molecule.

Conversion-only modules do not include sequencing-library preparation. Published [targeted EM-seq methods](https://www.nature.com/articles/s41467-025-57920-5) establish a precedent for conversion followed by locus PCR and indexed sequencing. They do not establish that our paired E5hmC targeted panel works without adaptation. The low-cost route is conditional on qualification of both chemistries, compatible amplification/indexing and controls before the full batch. If they require the complete manufacturer library kits or outsourcing, obtain a replacement budget; do not present this reagent estimate as an executable supplier protocol. An experienced lab's existing oxBS method is another acceptable route, with its own quote.

**Editors are not guaranteed to be cheaply procurable for C5R.** [Active Addgene #84475](https://www.addgene.org/84475/) and [inactive #84479](https://www.addgene.org/84479/) list $94 bacterial stabs but academic/nonprofit-only availability. Those prices are not included as industrial procurement. The main estimate assumes the executing lab already has a usable, permitted editor/control pair. If not, carry an additional **$2k–5k provisional construct-production reserve plus 15–30 staff hours and 1–3 elapsed weeks**, pending a real supplier/material-access quote. This reserve does not establish license availability or price; it is not covered by the routine contingency. Large-editor delivery may fail qualification even if the constructs are obtained. Stable-cell-line generation, viral production, sorting and extensive off-target profiling are outside this initial transient-delivery budget.

Existing mammalian culture and measurement equipment are prerequisites. Building those capabilities is a separate facility project, not an incidental project surcharge. Guide/primer designs, exact endpoint timing and the positive-control locus are not yet experimentally qualified. A 12-harvest qualification allowance is bounded screening capacity, not assurance of success. A failed qualification stops the pilot; repeated redesign needs a new allocation.

### Staff burden, distinct from cash

Task-based estimates for experienced staff using existing equipment; no automation time reduction is credited before workflow compatibility is demonstrated. Hours are hands-on/review time across the project, not continuous attendance during incubations.

| Work | Two regions, hours | Four regions, hours |
|---|---:|---:|
| Coordinate/annotation checks, guide/primer design and manifest | 6–10 | 8–16 |
| Guide cloning, existing-editor DNA preparation and QC | 10–18 | 12–24 |
| Cell stock establishment and stock QC | 6–10 | 6–10 |
| Delivery, primer and paired-assay qualification | 10–18 | 12–24 |
| Three culture/delivery batches, handling and harvests | 14–20 | 16–24 |
| Paired DNA/RNA extraction | 6–10 | 8–12 |
| Two conversions, locus PCRs, indexing, cleanup and pool QC | 16–24 | 18–30 |
| RNA plate setup and assay QC | 5–8 | 6–10 |
| Analysis/QC report | 6–10 | 6–10 |
| Procurement/coordination and deviations | 4–8 | 6–10 |
| **Total project work** | **83–136** | **98–180** |
| **Lab share if we provide design, analysis and coordination** | **67–108** | **78–144** |

The lab-share subtraction removes only the first, ninth and tenth rows. Our proposed contributions are not finished or accepted by a lab yet; experimental review/debugging remains in qualification and cannot be delegated to a document. Round to **70–110 / 80–145 lab hours**, approximately 2–3 / 2–4 staff workweeks spread across 4–6 calendar weeks. A rough staffing split is 15–30 specialist hours plus 50–80 operator hours for the smaller pilot; those labels overlap and are not an independently measured allocation.

At an assumed $100/hour staff value, cash including contingency plus the lab share gives **$19k–28k** for two regions and **$23k–36k** for four. At $150/hour, those ranges become approximately **$22k–33k** and **$27k–43k**. [Rice's published staff rates](https://research.rice.edu/bml/facilities/fee-service) are $110/hour external nonprofit and $120/hour for-profit, illustrating the scale; that biomaterials service is neither an epigenetics service nor a C5R quote. Contract-lab overhead, profit, instrument charges and formal deliverables can increase the paid quote beyond this collaboration budget. Covering reagents alone still asks them to contribute several thousand dollars of skilled time.

### C5R fit and a responsible eventual approach

[C5R's SciUniverse page](https://c5r.net/sciuniverse/) describes equipment and human-operator direction, PCR and cell-free protein tasks, and welcomes research collaborations. Its reported benchmark costs exclude lab and labor. The public examples reviewed do not establish mammalian culture, local epigenetic editing or paired 5mC/5hmC capability. Existing instruments and a permitted editor/assay workflow determine whether the above collaboration estimate applies.

The earlier proposed discussion concerned a two-region functional pilot, approximately **$15k consumables/services and 70–110 lab hours**, with design/analysis/report preparation supplied by us; whether they have the required workflow or an experienced partner; and whether a controlled computational-hypothesis-to-experiment case is worth that allocation. If editors need building, disclose the additional reserve and work rather than calling it a small favor. Publicity value is speculative, and the experiment must permit a negative result. A positive Stage A result supports local function; it does not yet establish preferential restoration or validate the model.

Before execution, the remaining provider-specific inputs are the editor source, paired-assay implementation, prepared-pool sequencing eligibility/price, available cell stock/equipment and staff estimate. They can refine an already counted scope rather than design the project from scratch. We should fund a defined qualification stop point and require a separate decision for repeats, hit confirmation or Stage B. No message has been sent and no capacity or cost commitment made.

## Broader, cheaper tests: revised priority (2026-10-03)

**Recommendation:** spend no experimental money yet. Use independent local-intervention evidence to define functional importance, then test its ability to predict recovery in existing human data across all supported regions. We have downloaded the public inputs below into ignored `data/raw/broad_validation/`. This turn audits availability and cost, not a newly executed test of functional restoration. No model-validation result is claimed from downloading files.

### Public data that can change the decision

| Priority | Resource | What is actually available | Role and limitation |
|---|---|---|---|
| 1, functional definition | [Tejedor et al. 2023](https://link.springer.com/article/10.1186/s13148-023-01546-1), [author deposit](https://zenodo.org/records/7761423) | Downloaded guide/region table and six well-level screen sheets: 56 promoter targets in HCT116/DLD1, active TET1, inactive TET1 and VP64 | Independent local-intervention fitness effects instead of gene-knockout essentiality; no recovery measurements. Do not classify null effects as unimportant marks without edit-success evidence. |
| 1, recovery outcome | [GSE158406 / Masalmeh et al.](https://www.nature.com/articles/s41467-020-20716-w) | Downloaded 27 processed HCT116 RRBS count files: WT, DNMT3B KO and DNMT1 KO, baselines, pulse/washout days 3/6/13/16/22/28/40 where available, late controls | Human recovery/enzyme-dependence and an H3K36me3 recruitment comparator. Coverage determines how many screened regions can join; seven listed days are experiment days, not seven days of continuous drug. Sparse trajectory replication and bulk selection remain limitations. |
| 2, mechanistic benchmark | [GSE129470 / Ginno et al.](https://www.nature.com/articles/s41467-020-16354-x) | Downloaded two 405-CpG amplicon time courses with replicate coverage and a 2,177,608-row genome table, including 860,406 `lowNoise` rate estimates | Large-scale enzyme-deletion dynamics, maintenance/de novo/TET assumptions and context dependence. Mouse ESCs, not human promoter fitness or a restoration pulse/washout. Author-fitted rates are not independent validation targets. |
| 3, cross-context replication | [GSE216989 / breast-cancer recovery study](https://www.nature.com/articles/s41594-023-01181-7) | Inventory verified: eight EPIC sample pairs of IDAT files and processed promoter-contact files; small series matrices contain metadata, not methylation values | Independent tissue context with early/late vehicle and drug/recovery samples, two replicates; processed beta matrices must be obtained from supplemental data or derived from IDATs. Full analysis inputs have not been downloaded. |

The priority-1 join is a **new promoter-level functional hypothesis**, distinct from the earlier gene-body essentiality test. It does not rehabilitate the failed/confounded gene-body annotation. Functional loss is defined using the independent active-versus-inactive intervention effect, corrected against each construct's controls, before examining recovery. Keep signed fitness effects; a cell-fitness penalty following loss of a methylated state is a clearer constraint measure than an absolute change of either sign. Cross-construct standardized z-scores have different scales and must not simply be subtracted as methylation-specific effects. Plate, batch and raw-signal columns are present; determine their biological replication meaning from methods before using them for uncertainty.

Freeze the test before joining outcome values:

1. Verify guide-region assembly by sequence/annotation checks, promoter overlap and local RRBS coverage. Use the edited region, not arbitrary CpGs elsewhere in the same gene. Screen regions and RRBS data cannot be joined blindly by gene symbol. GEO human metadata includes a contradictory boilerplate `mm10` alignment sentence, but its explicit genome-build field and publication identify hg38; document this discrepancy. The screen-region build still needs independent confirmation.
2. Define the functional score from HCT116 editing effects and controls, retaining uncertainty and edit-failure ambiguity. DLD1 is a sensitivity context, not another HCT116 replicate. Catalogue methylation-specific validation evidence for each candidate; a TET phenotype alone is not proof of mark causality.
3. Select recovery loci using only pretreatment/postpulse values and fixed coverage gates. Require overlap across early and late outcomes; report support, excluded regions and independent genes. If overlap or effect precision is poor, stop rather than expanding the definition after looking at recovery.
4. Compare context-only prediction against context plus frozen functional effect. Control baseline modified-C level, induced loss, CpG/CGI context, transcription and H3K36me3 where independent data are available. Train on early recovery and withhold days 28/40; hold out whole regions/genes, not neighboring CpGs. Keep a model with enzyme clearance distinct from a constant-rate benchmark.
5. Test DNMT3B dependence without re-estimating a separate functional relationship for every genotype. Constitutive KOs alter baseline context and cannot cleanly establish acute dependence. Bulk rebound, even if predicted, cannot distinguish within-lineage restoration from differential expansion. Cross-locus different slopes reject only narrow uniform-outgrowth explanations, not general selection.

An incremental held-out prediction from independently measured local consequences would be materially more useful than another age fit. Failure to add predictive value over chromatin/context weakens this version of functional constraint. A positive result motivates a much narrower mechanism experiment; it does not establish cell sensing of importance, an aging identity threshold, or protection over selection.

The new files total roughly 1 GB compressed, rather than tens of GB of raw sequencing. Individual RRBS files were unavailable at direct URLs; byte-range extraction from the public 6.26 GB archive retrieved only the requested processed counts. One interrupted connection was retried in smaller chunks. Source URLs, archive byte ranges, sizes and SHA-256 hashes are frozen in [the acquisition manifest](../code/results/broad_data_manifest.csv). Gzip/ZIP integrity is checked. These sources are not yet wired into `make validation` or analyzed by new Rust code; existing numerical tests are not evidence for these new hypotheses. Data purchase cost is zero; existing local compute avoids a new cloud bill. Annotation/QC and analysis still require work.

### Cost optimization for actual wet work

The economical design removes per-target editors, guide cloning and many tiny custom assay workflows. Use a global methylation perturbation with broad profiling at a lab that **already** cultures mammalian cells. The tradeoff is causal specificity: drug effects and differential growth require controls, and bulk methylation remains 5mC + 5hmC. This route screens recovery predictions across many loci; it does not individually perturb each site or establish methylation causality.

[Johns Hopkins' public EPIC 2.0 pricing](https://igc.jhmi.edu/service/methylation/) lists **$354.45/sample for 14–92 samples**, with nominal content of approximately **935,000 sites/sample**. External/commercial eligibility, included processing and surcharge changes need confirmation. It also lists a smaller approximately 270k-site screening array at $255/sample for 22 samples: $5,610 total, exceeding the $4,962.30 cost of 14 EPIC assays. For this small batch, the broader array is the better listed-price option. Passing and biologically informative probes will be fewer than nominal content; sites are not independent experimental replicates.

| Budget / USD | Archived material, 14 DNA assays | New broad pilot, 14 DNA + 12 RNA assays |
|---|---:|---:|
| EPIC, 14 × $354.45 | $4,962.30 | $4,962.30 |
| Paired extraction, up to 50 samples | Included in allowance below if needed | $480, existing Zymo price anchor |
| New cell stock if needed | $0 | $577, existing ATCC anchor |
| Medium/serum, drug, culture plastics and stock QC | $0 | $800–1,400 allowance |
| Gene-level 3-prime RNA profiling, 12 × $100–200 | Existing paired RNA required for functional-return analysis | $1,200–2,400 allowance, provider quote needed |
| Sample QC, preparation consumables and shipping | $500–1,200 allowance | $300–700 allowance |
| Cash subtotal | $5,462–6,162 | $8,319–10,519 |
| Cash with 20% contingency | **$6,555–7,395 (~$6k–8k)** | **$9,983–12,623 (~$10k–13k)** |
| Additional lab staff estimate | **8–16 hours** | **30–50 hours** |
| Including staff at assumed $100/hour | **~$7k–9k** | **~$13k–18k** |

RNA price is a planning allowance, not a commercial entitlement: [UT Austin lists internal TagSeq preparation/sequencing at $48.59 for 3–5M reads](https://utexas.atlassian.net/wiki/spaces/GSAF/pages/1416560685/GSAF%2BPrice%2BList). Confirm an external gene-level assay quote, adequate RNA quality and count-depth requirements. Gene-level 3-prime counts do not establish splice/cryptic-transcript mechanisms. The archived-material route assumes correctly collected, accessible DNA from a documented perturbation/recovery experiment; samples have not been secured and donor lab retrieval effort is not free. Without paired expression/fitness evidence it is a recovery-association test only.

The 14-DNA design uses **12 biological harvests plus two assay/QC duplicates**, not 14 biological replicates: three independent culture batches × four groups (vehicle and drug at pulse end; vehicle and recovered drug at one later endpoint). Measure division/viability. Two endpoints identify rebound relative to matched controls, not a recovery rate or held-out trajectory. C5R needs an experienced partner to supply the cultures/samples. New culture setup at C5R, targeted 5hmC discrimination, multiple cell lines, KOs, lineage tracing and paid computation/reporting are excluded. Existing-lab execution is a prerequisite, not a discount assumed for an inexperienced provider.

For a genuine early-to-late prediction test, use **24 biological harvests**: drug/vehicle × four endpoint times × three culture batches, including pulse end and three recovery times. At the same allowances, 24 DNA/RNA assays cost about **$17k–21k cash** and **45–75 lab hours**, approximately **$21k–29k including staff**. This is similar to the earlier two-region project's total cost while observing hundreds of thousands of loci, with new biological replication. It still requires an existing mammalian lab. Quality/coverage gates and functional annotations are frozen before outcomes, with later recovery withheld.

### What fits a lab without mammalian culture

A purified-enzyme/substrate panel is physically more compatible with PCR/protein/plate-reader infrastructure, but tests sequence-dependent chemistry rather than cellular functional protection. It is not an inexpensive substitute for the thesis test.

A counted example is 32 sequence contexts × three reaction times × three technical replicates = 288 reactions, plus 64 controls/qualification reactions. [Active Motif DNMT1](https://www.activemotif.com/catalog/details/31404/recombinant-dnmt1-protein) is listed at $505/20 micrograms. At an unqualified 0.1–0.5 micrograms/reaction, 352 reactions require 2–9 vials: **$1,010–4,545**. At the displayed synthesis scale, thirty-two individually methylated strands incur **$2,592** in 5mC modification charges at [IDT's displayed $81/100 nmol modification price](https://www.idtdna.com/site/Catalog/Modifications/Product/1101), before base synthesis and complementary strands; lower synthesis scales require a quote.

Allow $3,200–3,800 total for the synthetic substrate panel, $700–1,000 for a 400-reaction activity readout (unquoted reserve; [MTase-Glo](https://www.promega.com/products/epigenetics/methylation-analysis/mtase_glo-methyltransferase-assay/) currently requires a price inquiry), and $600–1,000 for buffers/consumables/shipping. Including 25% contingency gives roughly **$7k–13k cash plus 20–40 staff hours**, strongly dependent on qualified enzyme dose and oligo pricing. Luminescence reports total methyl transfer, not individual CpG restoration; assay backgrounds and activity linearity require qualification. Technical repeats are not biological replication.

Cheaper bacterial M.SssI can prepare methylated DNA, but its activity is not a test of human DNMT1 restoration. A standard cell-free sfGFP expression setup does not reproduce mammalian chromatin maintenance. We should not spend this money unless public data reveal a distinctive enzyme/substrate prediction that existing biochemical studies have not answered.

**Spending decision:** public-data analysis first; archived broad profiling only if it fills a specific data gap; a replicated broad live-cell recovery study at an existing specialist lab if the prediction survives. Retire the two-site pilot as the initial value proposition. No outreach, quote request or purchase has been sent.
