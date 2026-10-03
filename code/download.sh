#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
fetch() {
  local name="$1" url="$2"
  if [[ ! -s "data/raw/$name" ]]; then
    curl -fL --retry 3 --connect-timeout 20 -o "data/raw/$name.partial" "$url"
    mv "data/raw/$name.partial" "data/raw/$name"
  fi
}
base=https://ftp.ncbi.nlm.nih.gov/geo/series
fetch GSE225171_RAW.tar "$base/GSE225nnn/GSE225171/suppl/GSE225171_RAW.tar"
fetch GSE225171_extra_cov.tar.gz "$base/GSE225nnn/GSE225171/suppl/GSE225171_sample5498-5505_processed_cov.tar.gz"
fetch GSE179847_matrix.csv.gz "$base/GSE179nnn/GSE179847/suppl/GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz"
fetch GSE73115_processed.txt.gz "$base/GSE73nnn/GSE73115/suppl/GSE73115_processed.txt.gz"
for accession in GSE225171 GSE179847 GSE179848 GSE73115; do
  group="${accession:0:${#accession}-3}nnn"
  fetch "${accession}_family.soft.gz" "$base/$group/$accession/soft/${accession}_family.soft.gz"
done
fetch polycomb_cells.txt https://raw.githubusercontent.com/shmuel-ruppo/cellular-aging/main/Figures/poly_CpG_islands_meth_cell.txt
fetch annotations.csv https://raw.githubusercontent.com/shmuel-ruppo/cellular-aging/main/data/annotations.csv
fetch esl_annotation.xlsx 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69430-z/MediaObjects/41467_2026_69430_MOESM4_ESM.xlsx'
fetch sequence_fidelity_supplement.pdf 'https://www.biorxiv.org/content/biorxiv/early/2026/04/11/2026.04.09.717557/DC1/embed/media-1.pdf?download=true'
fetch GSE225172_metadata.txt.gz "$base/GSE225nnn/GSE225172/suppl/GSE225172_sc_MT_Babraham_Blood_aging.txt.gz"
fetch GSE179848_expression.csv.gz "$base/GSE179nnn/GSE179848/suppl/GSE179848_processed_cell_lifespan_RNAseq_data.csv.gz"
sha256sum -c code/results/SHA256SUMS
