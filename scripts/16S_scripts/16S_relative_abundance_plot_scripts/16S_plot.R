library(tidyverse)
library(reshape2) 
library(RColorBrewer)
library(ggplot2)
library(dplyr)

#PLOTTING
##With NEW designations
station_9_16S_tax_table <- read.delim("/Users/sabrinaelkassas/Desktop/16S/all_16S_count_table_normalized.txt", header = TRUE, sep = "\t")

#This selects the colors for the graphs
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


#ggplot code - 16S
ggplot(station_9_16S_tax_table, aes(x = fct_relevel(Sample, "9_40m", "9_70m", "9_380"), y = Percent, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x.bottom = element_text(angle = 45)) + scale_fill_manual(values=col_vector)

#PLOTTING
##With OLD designation of deltaproteobacteria
station_9_16S_delta_tax_table <- read.delim("/Users/sabrinaelkassas/Desktop/16S/all_16S_count_table_deltaproteobacteria_grouped_normalized.txt", header = TRUE, sep = "\t")

#This selects the colors for the graphs
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#ggplot code - 16S
ggplot(station_9_16S_delta_tax_table, aes(x = fct_relevel(Sample, "9_40m", "9_70m", "9_380"), y = Percent, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x.bottom = element_text(angle = 45)) + scale_fill_manual(values=col_vector)
