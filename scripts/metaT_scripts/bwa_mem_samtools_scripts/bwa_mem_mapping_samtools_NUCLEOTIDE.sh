#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=bwa_mem_samtools      # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=bwa_mem_samtools.log    # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate bwa_samtools 

#did this before running — bwa index dino_metzyme_annotated_coassembly.ffn

#bwa code to map to coassembly  

#40m
bwa mem -t 8 dino_metzyme_annotated_coassembly.ffn 30B8Z_S11_001_40m_mRNA.fasta  -o 30B8Z_S11_001_40m_mRNA_mapped.sam
#70m
bwa mem -t 8 dino_metzyme_annotated_coassembly.ffn 30B90_S12_001_70m_mRNA.fasta -o 30B90_S12_001_70m_mRNA_mapped.sam
#380m
bwa mem -t 8 dino_metzyme_annotated_coassembly.ffn 30B91_S28_001_380m_mRNA.fasta -o 30B91_S28_001_380m_mRNA_mapped.sam

#Flags 
#The -t flag specifies the number of threads to use
#The -o flag specifies the sam output name

#use samtools to convert sam to bam
#40m
samtools view -S -b 30B8Z_S11_001_40m_mRNA_mapped.sam > 30B8Z_S11_001_40m_mRNA_mapped.bam
#70m
samtools view -S -b 30B90_S12_001_70m_mRNA_mapped.sam > 30B90_S12_001_70m_mRNA_mapped.bam
#380m
samtools view -S -b 30B91_S28_001_380m_mRNA_mapped.sam > 30B91_S28_001_380m_mRNA_mapped.bam

#Flags
#The -S flag indicates you have a sam input
#The -b flag indicates you want to have a bam output

#use samtools to sort bam files
#40m
samtools sort 30B8Z_S11_001_40m_mRNA_mapped.bam -o 30B8Z_S11_001_40m_mRNA_mapped_sorted.bam
#70m
samtools sort 30B90_S12_001_70m_mRNA_mapped.bam -o 30B90_S12_001_70m_mRNA_mapped_sorted.bam
#380m
samtools sort 30B91_S28_001_380m_mRNA_mapped.bam -o 30B91_S28_001_380m_mRNA_mapped_sorted.bam

#Flags
#The -o flag specifies output file
