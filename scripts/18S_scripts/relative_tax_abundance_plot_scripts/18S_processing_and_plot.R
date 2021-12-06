library(tidyverse)
library(reshape2) 
library(RColorBrewer)
library(ggplot2)
library(dplyr)

#Read in each 18S .csv file, but will have to make formatting changes, hence the variable name "bad"
bad40m <- read.csv("/Users/sabrinaelkassas/Desktop/WHOI021_V9_S105_L001_40m_18S_taxonomy_2500_subsampled.csv", header = FALSE)
bad70m <- read.csv("/Users/sabrinaelkassas/Desktop/WHOI040_V9_S124_L001_70m_18S_taxonomy_2500_subsampled.csv", header = FALSE)
bad380m <- read.csv("/Users/sabrinaelkassas/Desktop/WHOI010_V9_S94_L001_380m_18S_taxonomy_2500_subsampled.csv", header = FALSE)

#This then reformats the columns correctly for the tax_table code
good_40m <- bad40m %>% filter(!is.na(V1)) %>% 
  select(x = V1)
good_70m <- bad70m %>% filter(!is.na(V1)) %>% 
  select(x = V1)
good_380m <- bad380m %>% filter(!is.na(V1)) %>% 
  select(x = V1)

write.csv(good_40m, "/Users/sabrinaelkassas/Desktop/WHOI021_V9_S105_L001_40m_18S_taxonomy_2500_subsampled.NEW.csv")
write.csv(good_70m, "/Users/sabrinaelkassas/Desktop/WHOI040_V9_S124_L001_70m_18S_taxonomy_2500_subsampled.NEW.csv")
write.csv(good_380m, "/Users/sabrinaelkassas/Desktop/WHOI010_V9_S94_L001_380m_18S_taxonomy_2500_subsampled.NEW.csv")

#Generating the tax_table - Did this on poseidon

#for loop to import rRNA read taxonomy assignments
#Credits: Sarah Hu
##12/01/21 18S rRNA for dino_METZYME Project 

#install packages
library(tidyverse)

#Set up files
#select correct files for HPC
taxa_raw <- list.files(path = "/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/18S_reads/trimmed_18S/mothur_18S_taxonomy/subsampled_18S", pattern = "subsampled.NEW.csv", full.names = FALSE)

#path to files for HPC
path_data <- "/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/18S_reads/trimmed_18S/mothur_18S_taxonomy/subsampled_18S/"

#For loop

for(a in taxa_raw){
  #import files, use paste to string together path and file names
  imported_tax <- read.csv(paste(path_data, a, sep = ""))
  #Extract sample name from "a", and split at ".wang"
  sample_names <- unlist(strsplit(a, ".wang"))
  #Modify imported data
  output_tmp <- imported_tax %>%
    #Adding in the sample name
    mutate(SAMPLE = sample_names[1]) %>%
    #filter out unknowns
    filter(!(grepl("unknown_unclassified", x))) %>%
    select(useful = x) %>%
    separate(useful, into = c("ACCESSION_NUMBER", "taxonomy"), sep = "\t") %>%
    # use regex to modify taxonomy column
    mutate(new_tax = str_replace_all(taxonomy, pattern = "\\(\\d+\\)", replacement = "")) %>%
    # parse taxonomy lineage name by semicolon
    separate(new_tax, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
    # add artifical count column
    add_column(COUNT = 1) %>%
    # add sample information
    group_by(sample_names[1], Phylum, Class, Order, Family, Genus, Species) %>%
    summarise(SUM = sum(COUNT))
  cat("Processing...", sample_names[1], "/n/n")
  # if else statement to facilitate row bind
  if (!exists("tax_table")){
    tax_table <- output_tmp
  } else {
    tax_table <- bind_rows(tax_table, output_tmp)
  }
  rm(output_tmp)
}
#rm(output_tmp)
#run the rm(tax_table) every time you run the for-loop so it doesn't add files twice.
#rm(tax_table)

write.csv(tax_table, file  = "output-tax-table_18S.csv")


#PLOTTING
station_9_tax_table <- read.delim("/Users/sabrinaelkassas/Desktop/output-tax-table_18S_normalized.txt", header = TRUE, sep = "\t")

This selects the colors for the graphs
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


#ggplot code - 18S
ggplot(station_9_tax_table, aes(x = fct_relevel(Sample, "9_40m", "9_70m", "9_380"), y = Percent, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) + theme(axis.text.x.bottom = element_text(angle = 45)) + scale_fill_manual(values=col_vector)
