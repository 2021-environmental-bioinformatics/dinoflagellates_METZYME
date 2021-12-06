#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=fasta36_plastids                  # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=180gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=fasta36_plastids%j.log              # Standard output/error

export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate fasta36

fasta36 -U WHOI010_16S_AGAGTCAC_380m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/pr2_version_4.14.0_SSU_taxo_long.fasta.gz > WHOI010_16S_AGAGTCAC_380m_plastids_4.14.0
