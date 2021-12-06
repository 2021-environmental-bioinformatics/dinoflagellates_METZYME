#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=FragGeneScan_annot    # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=12               # Number of CPU cores per task
#SBATCH --mem=80gb                       # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=FragGeneScan_annot.log  # Standard output/error
#export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate FragGeneScan

#Coassembly
FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/metaT_assembly/transcripts.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/dino_metzyme_annotated_coassembly -w 0 -t illumina_5 -p 8

#mRNA file for each depth. This will be necessary for the diamond step later on in the pipeline. 

FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B8Z_S11_001_40m_mRNA.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B8Z_S11_001_40m_mRNA_annotated -w 0 -t illumina_5 -p 8

FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B90_S12_001_70m_mRNA.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B90_S12_001_70m_mRNA_annotated -w 0 -t illumina_5 -p 8

FragGeneScan -s /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B91_S28_001_380m_mRNA.fasta -o /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/30B91_S28_001_380m_mRNA_annotated -w 0 -t illumina_5 -p 8

#Flags
#The -s flag is the assembly path. Please note that the output of RNASPAdes that you want to use as the input for FragGeneScan is the transcripts.fasta file in the main assembly folder. 
#The -o flag is the output directory 
#The -w flag is the read type; 0 if the sequence file has short sequence reads and 1 if the sequence file has complete genomic sequences
#The -t flag is the train file that FragGeneScan should use based on the sequencer of the reads
#The -p flag is the number of threads to use. default is 1. 

