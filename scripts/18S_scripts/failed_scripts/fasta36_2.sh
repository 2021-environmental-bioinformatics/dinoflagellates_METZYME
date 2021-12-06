#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=fasta36_18S2                       # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=180gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=fasta36_18S%j.log                   # Standard output/error

export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate fasta36

cd /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/18S_reads/trimmed_18S

fasta36 -U  WHOI040-V9_S124_L001_70m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/Protist_Ribosomal_reference_db_18S/pr2_version_4.11.1_taxo_long.fasta > WHOI040-V9_S124_L001_70m_18S_taxonomy
