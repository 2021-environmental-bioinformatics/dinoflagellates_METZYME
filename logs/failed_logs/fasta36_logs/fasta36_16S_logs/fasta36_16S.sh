#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=fasta36                           # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=elks@mit.edu                     # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=120gb                                   # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=fasta36_16S%j.log                   # Standard output/error

export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate fasta36

cd /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/merged_16S

#40m, SSU
fasta36 -U  WHOI040_16S_TCGACGAG_40m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/SSURef_111_tax_silva.fasta > WHOI040_16S_TCGACGAG_40m_SSU_taxonomy

#40m, LSU
fasta36 -U WHOI040_16S_TCGACGAG_40m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/LSURef_111_tax_silva.fasta > WHOI040_16S_TCGACGAG_40m_LSU_taxonomy

#70m, SSU
fasta36 -U WHOI020_16S_TCGACGAG_70m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/SSURef_111_tax_silva.fasta > WHOI020_16S_TCGACGAG_70m_SSU_taxonomy

#70m, LSU
fasta36 -U WHOI020_16S_TCGACGAG_70m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/LSURef_111_tax_silva.fasta > WHOI020_16S_TCGACGAG_70m_LSU_taxonomy

#380m, SSU
fasta36 -U WHOI010_16S_AGAGTCAC_380m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/SSURef_111_tax_silva.fasta > WHOI010_16S_AGAGTCAC_380m_SSU_taxonomy

#380m, LSU
fasta36 -U WHOI010_16S_AGAGTCAC_380m_merged.fasta /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/LSURef_111_tax_silva.fasta > WHOI010_16S_AGAGTCAC_380m_LSU_taxonomy

mv *taxonomy* /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/merged_16S/taxonomy
