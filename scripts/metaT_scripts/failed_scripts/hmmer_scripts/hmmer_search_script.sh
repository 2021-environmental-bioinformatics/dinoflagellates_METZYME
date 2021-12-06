#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=hmmsearch             # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=hmmsearch%j.log         # Standard output/error
#export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate hmmer

#40m
hmmsearch --domtblout=30B8Z_S11_001_40m_hmmersearch_out -E 0.01 /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/Pfam-A.hmm 30B8Z_S11_001_40m_mRNA_annotated.faa

#70m
hmmsearch --domtblout=30B90_S12_001_70m_hmmersearch_out -E 0.01 /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/Pfam-A.hmm 30B90_S12_001_70m_mRNA_annotated.faa

#380m
hmmsearch --domtblout=30B91_S28_001_380m_hmmersearch_out -E 0.01 /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/Pfam-A.hmm 30B91_S28_001_380m_mRNA_annotated.faa

#coassembly
hmmsearch --domtblout=dino_metzyme_annotated_coassembly_hmmersearch_out -E 0.01 /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/Pfam-A.hmm dino_metzyme_annotated_coassembly.faa

