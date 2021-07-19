#!/usr/bin/env Rscript

## Rscript to process TPM and read count matrices so that columns match covariate files.
## Prepare matrices for input to PEER.

## Load required packages
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(optparse)

##------------- FUNCTIONS



##------------- MAIN

## Read command line arguments
option_list = list(make_option(c('--OUT'), type = 'character', default = NULL, help = 'output file path'),
                   make_option(c('--ANNOT'), type = 'character', default = NULL, help = 'path to gene-level feature annotation file'),
                   make_option(c('--ZSCORES'), type = 'character', default = 70, help = 'path to zscore file'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

zscore.file = opt$ZSCORES
out.file = opt$OUT 
annot.file = opt$ANNOT


# annot.file = '/work-zfs/abattle4/bstrober/random_projects/feature_generation_for_victor_and_jessica/river_input_gene_level.txt'
# zscore.file = ${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt
# out.file = ${datadir}/data_prep/RIVER/river_input_v8_african_all_07-19-2021.txt


# read in gene-level annotation file
annot_df <- read.table(annot.file, header=TRUE) %>% rename(GENE=GeneName)

# read in zscores
z_scores <- read.table(zscore.file, header=TRUE, check.names = FALSE)

# reformat zscore table to wide format with a row for every (Subject ID, Gene) pair
zscore_table <- z_scores %>% 
  gather(key="SubjectID",value="Zscore",-GENE) %>% 
  rename(GeneName=GENE)
zscore_table <- as.data.table(zscore_table)

# create a lookup
setkey(zscore_table,GeneName,SubjectID)

# For every row of the annotation file, retrieve the zscore for the (subject, gene) pair
zscore_list <- vector("numeric",nrow(annot_df))
for (i in seq(nrow(annot_df))) {
  subject <- as.character(annot_df[i,1])
  gene <- as.character(annot_df[i,2])
  zscore_list[i] <- zscore_table[.(gene,subject)]$Zscore
}

# add zscores to annot_df
annot_df$Zscore <- zscore_list

## move N2pair column to be the last column
annot_df <- annot_df %>% select(-N2pair,N2pair) 

# drop rows with NA in zscores (NEED TO ACCOUNT FOR N2)
annot_df <- annot_df[complete.cases(zscore_list), ]

# save version that keeps all features
write.table(x = annot_df, row.names = FALSE, quote = FALSE, file = out.file)

# keep features that are common between european and african data
# afr_df <- afr_df %>% select(SubjectID, GeneName, GC, CpG, bStatistic, priPhCons, mamPhCons, verPhCons, priPhyloP, mamPhyloP, verPhyloP, GerpN, GerpS, PHRED, SIFTcat_deleterious, SIFTcat_tolerated, SIFTval, PolyPhenCat_benign, PolyPhenCat_possibly_damaging, PolyPhenCat_probably_damaging, PolyPhenCat_unknown, PolyPhenVal, dist_to_tss, Zscore,  N2pair)

# save version with common features
# write.table(x = afr_df, row.names = FALSE, quote = FALSE, file = afr_data_common)

