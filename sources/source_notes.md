# Source and data provenance

The annotated literature review is in [notes/reading_list.md](../notes/reading_list.md). It distinguishes experimental mechanism evidence, association, mathematical precedent and preprints.

Data downloaded on 2026-10-03 into the ignored `data/raw/` directory. Exact bytes are recorded in [SHA256SUMS](../code/results/SHA256SUMS). NCBI download access initially failed in the network sandbox and succeeded with an approved escalation before implementation.

| Source | Version / use | Access limitations |
|---|---|---|
| [GSE225171](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225171) | Complete 3.0 GB RAW tar plus separate 805 MB coverage archive and GEO SOFT metadata downloaded. | Raw coverage retained; the current real-data run uses the published processed Polycomb summary rather than rebuilding genomic averages. |
| [GSE179847](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179847) | Full processed beta/detection matrix and SOFT metadata. 479 samples, nine donors. | No physiology/RNA joins in this run. |
| [GSE73115](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73115) | Full processed beta/detection matrix and SOFT metadata. 180 arrays, 86 people. | Whole-blood composition unadjusted. Technical duplicate arrays pooled per person/year. |
| [Published Polycomb repository](https://github.com/shmuel-ruppo/cellular-aging) | Git tree 805eb382cc65f91b9953661550af325c48582ecd; processed cell summary and annotation downloaded. | 1,131 rows versus 1,132 reported cells. Referenced MATLAB MC helper absent from the inspected tree. |
| [scAgingPaper](https://github.com/EpigenomeClock/scAgingPaper) | Git tree 07fc2b743e6cdd04a0392a4ac0ed97202c543eed; inspected mouse-ID derivation. | Sample prefixes verified against preprocessing convention, not newly reconciled with all raw libraries. |
| [GSE225172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225172) | Paired expression matrix downloaded for future validation. | Not analyzed or used to fit q or a biological threshold. |

All dataset discovery was read-only. No external accounts, private datasets or communications were used. Publisher citations in the review were checked through publisher pages, indexed primary abstracts or PubMed/PMC. Full-text access was not successful for every paper; mathematical baseline reproduction is limited to the explicitly described Polycomb approximation.

External annotation follow-up: the ESL supplementary workbook was downloaded from the publisher and the valid hexamer-ranking supplement from a bioRxiv DC1 URL after initial 429/browser-challenge responses. Their checksums are included in `SHA256SUMS`. `epidrift annotate` validates 31,744 ESL IDs and 225 ranks, then merges hg19 GEO platform sequences/context into ignored `data/derived/human_features.tsv`. Motif scores average forward and reverse-complement ranks. Annotation coverage is recorded in `code/results/annotation_coverage.csv`. These source definitions were fixed before inspecting the matched ESL effect; the exploratory matching choices were added after the initial baseline run.
