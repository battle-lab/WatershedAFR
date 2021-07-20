#!/usr/bin/env Rscript

rm(list = ls())

## R script to process VCFtools output containing the allele counts for each V8 individual at the V8 cis-eQTL calls

library(dplyr)
library(reshape2)


dir = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep'
prefix = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs.012'
outpath = paste0(dir,'/','/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt')
#--------------- MAIN

# Load the list of individuals, the site position info, and the allele counts
inds = read.table(paste(prefix, '.indv', sep = ''), header = F, stringsAsFactors = F)[, 1]
pos = read.table(paste(prefix, '.pos', sep = ''), header = F)
genos = t(as.data.frame(fread(prefix, na.strings = c('-1'))))[-1, ]

# Combine the position and genotype info into a single data frame
# Add the individual IDs as headers
out = cbind(pos, genos)
colnames(out) = c('Chrom', 'Pos', inds)

# Remove duplicated positions
out = out %>% mutate(ID = paste(Chrom, '_', Pos, sep = ''))
ids = as.data.frame(table(out$ID))
ids.to.keep = ids %>% filter(Freq == 1)
out = out %>% filter(ID %in% ids.to.keep$Var1) %>% select(-ID)

# Write out the combined data frame
write.table(out, outpath, sep = '\t', col.names = T, row.names = F, quote = F)
