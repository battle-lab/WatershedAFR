#!/usr/bin/env Rscript

## Makes file with following columns
## Col1: Gene
## Col2: Indiv
## Col3: Chrom
## Col4: Start
## Col5: End
## Col6: Ref
## Col7: Alt

## Col1 is restricted to genes that have at least one multi-tissue outlier individual
## Col2 is all individuals that have at least one rare variant within 10kb +/- window
## around the gene in that row
## Col3 and Col4 follow bed coordinates of "0-start, half-open"

library(dplyr)
library(tidyr)
library(data.table)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--rv_sites'), type = 'character', default = NULL, help = 'bed file with rare variants overlapping 10kb +/- window around the genes'),
                   make_option(c('--popname'), type = 'character', default = NULL, help = 'population of individuals for output file name'),
                   make_option(c('--outdir'), type = 'character', default = NULL, help = 'directory to save list of rare variants per gene-individual pair'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

rv_sites_file = opt$rv_sites
popname = opt$popname
outdir = opt$outdir

# Open file
rv_sites = read.table(rv_sites_file, header = FALSE, check.names = FALSE)

# get relevant columns
df = rv_sites[,c("V10","V11","V1","V2","V3","V4","V5")]
colnames(df) = c("Gene","Ind","Chrom","Start","End","Ref","Alt")

# sort by gene
df = arrange(df, Gene)

# make each rare variant of the form chrom;position;major_allele;variant_allele
df = unite(df, RV, sep = ";", Chrom:Alt)

# save
filename = paste0(outdir,'/gene-', popname, '-rv.txt')
write.table(df, filename, sep = '\t', row.names = FALSE, quote = FALSE)
print(paste0("Saved to ",filename))
