#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${1:-}" == "--broad" ]]; then
  data/tools/format-env/bin/python code/convert_public_data.py broad-download
  exit
fi
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
for accession in GSE225171 GSE179847 GSE179848 GSE73115 GSE145698 GSE282971; do
  group="${accession:0:${#accession}-3}nnn"
  fetch "${accession}_family.soft.gz" "$base/$group/$accession/soft/${accession}_family.soft.gz"
done
fetch polycomb_cells.txt https://raw.githubusercontent.com/shmuel-ruppo/cellular-aging/main/Figures/poly_CpG_islands_meth_cell.txt
fetch annotations.csv https://raw.githubusercontent.com/shmuel-ruppo/cellular-aging/main/data/annotations.csv
fetch esl_annotation.xlsx 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69430-z/MediaObjects/41467_2026_69430_MOESM4_ESM.xlsx'
fetch sequence_fidelity_supplement.pdf 'https://www.biorxiv.org/content/biorxiv/early/2026/04/11/2026.04.09.717557/DC1/embed/media-1.pdf?download=true'
fetch GSE225172_metadata.txt.gz "$base/GSE225nnn/GSE225172/suppl/GSE225172_sc_MT_Babraham_Blood_aging.txt.gz"
fetch GSE179848_expression.csv.gz "$base/GSE179nnn/GSE179848/suppl/GSE179848_processed_cell_lifespan_RNAseq_data.csv.gz"
fetch ReactomePathways.gmt.zip https://reactome.org/download/current/ReactomePathways.gmt.zip
fetch GSE145698_counts.txt.gz "$base/GSE145nnn/GSE145698/suppl/GSE145698_counts_matrix.txt.gz"
maintenance=https://raw.githubusercontent.com/Read-Lab-UCI/DNA-methylation-maintenance-kinetics/ec46d95b68fad06072702686c248cad55754fc3b
fetch maintenance_chr1.tar.gz "$maintenance/data/ReadData/AllDat_chr1.tar.gz"
fetch maintenance_chr22.tar.gz "$maintenance/data/ReadData/AllDat_chr22.tar.gz"
fetch maintenance_README.md "$maintenance/README.md"
fetch maintenance_MLEInference.m "$maintenance/MLEInference.m"
clone=https://raw.githubusercontent.com/veltenlab/EPI-clone/0e612a5dffe060316dfc2a25ebce7fcedcf5c952
fetch EPIClone_M123.rds https://ndownloader.figshare.com/files/42479346
fetch EPIClone_Figure1.Rmd "$clone/figures/Figure1/Figure1.Rmd"
fetch EPIClone_cpg_selection.csv "$clone/figures/Figure1/cpg_selection.csv"
fetch EPIClone_processing.Rmd "$clone/processing_vignette.Rmd"
fetch EPIClone_panel.tsv "$clone/infos/panel_info_dropout_pwm.tsv"
fetch EPIClone_mouse3_umi.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8653nnn/GSM8653743/suppl/GSM8653743_LARRY_mouse3-umi-counts.tsv.gz
fetch EPIClone_mouse4_umi.tsv.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8653nnn/GSM8653744/suppl/GSM8653744_LARRY_mouse4-umi-counts.tsv.gz
fetch GSE51810_series_matrix.txt.gz "$base/GSE51nnn/GSE51810/matrix/GSE51810_series_matrix.txt.gz"
fetch GSE51811_series_matrix.txt.gz "$base/GSE51nnn/GSE51811/matrix/GSE51811_series_matrix.txt.gz"
fetch GSE51811_family.soft.gz "$base/GSE51nnn/GSE51811/soft/GSE51811_family.soft.gz"
fitness=https://raw.githubusercontent.com/hart-lab/bagel/53388adbb4fb0931e5c9dda135502be19e4555f0
fetch BAGEL_CEGv2.txt "$fitness/CEGv2.txt"
fetch BAGEL_NEGv1.txt "$fitness/NEGv1.txt"
sha256sum -c code/results/SHA256SUMS
