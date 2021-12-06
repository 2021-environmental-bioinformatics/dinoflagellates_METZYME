#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=kraken2               # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=kraken2_paired.log             # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate kraken2

#40m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva --output WHOI040_16S_TCGACGAG_40m_taxonomy --report 40m.kreport --paired /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI040_16S_TCGACGAG_R1_val_1.fq /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI040_16S_TCGACGAG_R2_val_2.fq

#70m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva --output WHOI020_16S_TCGACGAG_70m_taxonomy --report 70m.kreport --paired /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI020_16S_TCGACGAG_R1_val_1.fq /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI020_16S_TCGACGAG_R2_val_2.fq

#380m, SSU
kraken2 --db /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/silva111/16S_SILVA132_k2db/kraken_silva --output WHOI010_16S_AGAGTCAC_380m_taxonomy --report 380m.kreport --paired /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI010_16S_AGAGTCAC_R1_val_1.fq /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/16S_reads/trimmed_16S/WHOI010_16S_AGAGTCAC_R2_val_2.fq

#Flags
#The --db flag specifies the database to check reads against to id taxonomy
#The --output flag names the output file
#The --report flag should generate a kreport, though more processing is necessary to generate this.
#The —-paired flag indicated the input are two paired forward and reverse reads
