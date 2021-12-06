#!/bin/bash

# Katy - editing to investigate source of job failure

#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=fasta36                           # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu                     # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=50gb                                   # Job memory request
#SBATCH --time=20:00:00                              # Time limit hrs:min:sec
#SBATCH --output=fasta36_16S%j.log                   # Standard output/error

export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate fasta36

cd /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/merged_16S

#40m, SSU
fasta36 -U /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI010_16S_AGAGTCAC_380m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/SSURef_111_tax_silva.fasta > testfile_SSU_taxonomy


