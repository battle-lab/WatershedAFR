#!/usr/bin/env Rscript

## Makes file with following columns
## Col1: Gene
## Col2: Indiv
## Col3: Comma separated list of variants
## 
## Col1 is restricted to genes that have at least one multi-tissue outlier individual
## Col2 is all individuals that have at least one rare variant within 10kb TSS
## of the gene in that row
## Col3 is the rare variant of the form chrom;position;major_allele;variant_allele

library(dplyr)
library(tidyr)
library(data.table)

# Read command line arguments
outdir = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants'
rv_sites_file = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants/rv_bed_EUR/all_rv_sites.bed'

rv_sites = read.table(rv_sites_file, header = FALSE, check.names = FALSE)

df = rv_sites[,c("V10","V11","V1","V2","V4","V5")]
colnames(df) = c("Gene","Ind","Chrom","Pos","Ref","Alt")

df = arrange(df, Gene)

df = unite(df, RV, sep = ";", Chrom:Alt)

filename = paste0(outdir,'/','gene-EUR-rv.txt')
write.table(df, filename, sep = '\t', row.names = FALSE, quote = FALSE)