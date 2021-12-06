#Station 9 40m
trim_galore --paired --fastqc WHOI021-V9_S105_L001_R1_001.fastq WHOI021-V9_S105_L001_R2_001.fastq
#Station 9 70m
trim_galore --paired --fastqc WHOI040-V9_S124_L001_R1_001.fastq WHOI040-V9_S124_L001_R2_001.fastq
#Station 9 380m
trim_galore --paired fastqc WHOI010-V9_S94_L001_R1_001.fastq WHOI010-V9_S94_L001_R2_001.fastq

#The --paired flag indicates the input is paired-end reads
#The --fastqc flag tells trimgalore! to run a quality check on the reads post-trimming
