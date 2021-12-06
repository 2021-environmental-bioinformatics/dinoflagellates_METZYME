#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=fasta36_2                           # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=elks@mit.edu                     # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=180gb                                   # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=fasta36_2_16S%j.log                   # Standard output/error

export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate fasta36

cd /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/merged_fasta

#70m head, SSU
fasta36 -U head_WHOI020_16S_TCGACGAG_70m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/SSURef_111_tax_silva.fasta > head_WHOI020_16S_TCGACGAG_70m_SSU_taxonomy

#70 tail, SSU 
fasta36 -U tail_WHOI020_16S_TCGACGAG_70m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/SSURef_111_tax_silva.fasta > tail_WHOI020_16S_TCGACGAG_70m_SSU_taxonomy
