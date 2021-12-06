#Station 9 40m
trim_galore --paired --fastqc WHOI040_16S_TCGACGAG_R1.fastq.gz WHOI040_16S_TCGACGAG_R2.fastq.gz
#Station 9 70m
trim_galore --paired --fastqc WHOI020_16S_TCGACGAG_R1.fastq.gz WHOI020_16S_TCGACGAG_R2.fastq.gz
#Station 9 380m
trim_galore --paired --fastqc WHOI010_16S_AGAGTCAC_R1.fastq.gz  WHOI010_16S_AGAGTCAC_R2.fastq.gz

#The --paired flag indicates the input is paired-end reads
#The --fastqc flag tells trimgalore! to run a quality check on the reads post-trimming
