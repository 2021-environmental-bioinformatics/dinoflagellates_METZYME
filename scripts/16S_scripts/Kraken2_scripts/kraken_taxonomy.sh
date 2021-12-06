#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=kraken_taxon               # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=10gb                      # Job memory request
#SBATCH --time=02:00:00                  # Time limit hrs:min:sec
#SBATCH --output=kraken_taxon.log             # Standard output/error
#export OMP_NUM_THREADS=1

module load anaconda/5.1
source activate kraken2

#We want to use KrakenTools to generate a kraken report from the subsampled output (since the original report is valid for all samples)
#First, use make_ktaxonomy.py, then make_kreport.py

#Set directories
script_dir=/vortexfs1/home/kaabbott/.conda/pkgs/krakentools-1.2-pyh5e36f6f_0/python-scripts
database_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva
data_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/kraken_taxonomy

cd $data_dir

for taxonomy in WHOI*_no_plastids_6000_subsampled
do
	python ${script_dir}/make_kreport.py -i $taxonomy -t ${database_dir}/16S_SILVA132_taxonomy_file -o ${taxonomy}_report
done

#rename
mv WHOI040_16S_TCGACGAG_40m_no_plastids_6000_subsampled_report 40m_16S_no_plastids_subsample_6000.kreport
mv WHOI020_16S_TCGACGAG_70m_no_plastids_6000_subsampled_report 70m_16S_no_plastids_subsample_6000.kreport
mv WHOI010_16S_AGAGTCAC_380m_no_plastids_6000_subsampled_report 380m_16S_no_plastids_subsample_6000.kreport
