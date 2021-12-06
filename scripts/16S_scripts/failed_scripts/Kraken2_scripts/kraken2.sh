#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=kraken2               # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=kraken2.log             # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate kraken2

#40m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db --output WHOI040_16S_TCGACGAG_40m_taxonomy --report 40m.kreport WHOI040_16S_TCGACGAG_40m_merged.fasta

#70m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db --output WHOI020_16S_TCGACGAG_70m_taxonomy --report 70m.kreport WHOI020_16S_TCGACGAG_70m_merged.fasta

#380m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db --output WHOI010_16S_AGAGTCAC_380m_taxonomy --report 380m.kreport WHOI010_16S_AGAGTCAC_380m_merged.fasta
