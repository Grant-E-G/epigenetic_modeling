# Research Log

Use this file for dated progress notes.

## 2026-10-09: Dog Aging Project attempted falsification

User requests an adversarial test of distinctive vital-region protection, explicitly permitting model invalidation. Read the repository scientific contract and separated the historical direct-protection claim from revision-3 external-fitness clone transfer. Freeze [the protocol](dog_aging_protocol.md) before matrix outcomes; use the existing external BAGEL lists and unique NCBI human/dog orthologs, exact canFam4 annotations/reference sequence, training-only matching, and held-out dogs. Acquire/hash 15 sources and validate all 335,363,600 matrix counts. Cohort: 894 dogs/1,640 samples, 538 repeated; metadata-eligible adult longitudinal split: 390 training/100 evaluation dogs.

Source-format audit finds 184,152 unsuffixed CpG regions plus 20,338 labeled aggregates; exclude aggregate rows to avoid duplicating underlying CpGs. All 2,000 audited starts are CpGs under one-based coordinates, none under zero-based. Exclude six functional-list entries lacking gene IDs rather than guessing symbols. Gene-level essentiality does not prove local methylation importance.

Initial interval-only context match gives apparent increased essential-gene-body disruption, but singleton regions have uninformative sequence covariates. Preserve that result and protocol clarification; the post-result +/-250-bp flank correction reverses the estimate. Corrected adult >=20-read result: 75 matched gene pairs, 96 observed evaluation dogs, 942 dog/pair records, estimated protection 32.6%; crossed 95% interval -62.0% to +59.3%, simultaneous 97.5% interval -81.5% to +63.2%. At >=10 and >=30 reads the point estimates are -2.1% and -12.3%, respectively, also inconclusive. Promoter support is only 2–3 matched gene pairs. Do not select a favorable result.

Only 12.56% of potential primary gene-body dog/pair records are jointly covered; median observed gene pair has eight dogs. All known first/next preparation batches differ; no relatedness matrix, direct local-effect calibration or longitudinal lineage/count observations are provided by the acquired resources. Conditional independent-dog endpoint-power simulation detects 10%/20% protection only 4.33%/10.36% of the time. This cannot exclude modest protection or adjudicate intrinsic repair/clone selection.

Decision: no robust support, no defensible whole-model invalidation, no new experiment/target panel. The narrowly specified practical proxy prediction remains unresolved. Broad vital-region protection needs independently fixed local importance, state/reference and an observable minimum effect; allowing unknown scores/arbitrarily small effects cannot earn empirical support. Prior stopping decision remains. [Full report](../code/results/dog_aging_falsification.md). Software checks: 37 Rust tests pass; strict Clippy, formatting, Python Black/mypy and all source hashes pass.

## YYYY-MM-DD

- 
