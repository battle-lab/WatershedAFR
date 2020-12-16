#!/usr/bin/env Rscript

## R script to split the entire GTEx matrices into files for each tissue
## then process RPKM and read count matrices so that columns match covariate files.
## Prepare matrices for input to PEER.

## Load required packages
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

##------------- MAIN

cov_dir='/scratch/groups/abattle4/victor/WatershedAFR/raw_data/GTEx/GTEx_Analysis_v8_eQTL_covariates'
outdir='/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep'

# open covariate files for each tissue
cov_files = list.files(path=cov_dir, pattern="*covariates.txt", full.names=TRUE)
cov_list = lapply(cov_files, fread)

# bind them together
combined_cov = rbindlist(cov_list, fill = TRUE)

# get PC1 - PC5 and sex for each subject
combined_cov = combined_cov %>% arrange(ID) %>% filter(str_detect(ID, 'PC|sex')) %>%
  fill(colnames(combined_cov), .direction = "downup") %>% distinct(ID, .keep_all=TRUE)

# transpose
combined_cov = t(combined_cov)
rownames(combined_cov)[1] = "SUBJID"

# write
filename = paste0(outdir,'/gtex_v8_eQTL_covariates.txt')
write.table(combined_cov, filename, col.names = FALSE, row.names = TRUE, quote = FALSE, sep = '\t')
