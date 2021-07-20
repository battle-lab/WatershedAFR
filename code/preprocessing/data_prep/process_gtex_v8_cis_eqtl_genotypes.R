#!/usr/bin/env Rscript

rm(list = ls())

## R script to process VCFtools output containing the allele counts for each V8 individual at the V8 cis-eQTL calls

library(dplyr)
# library(reshape2)


# dir = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep'
# prefix = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs.012'
# outpath = paste0(dir,'/','/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt')

#--------------- MAIN

## Read command line arguments
option_list = list(make_option(c('--INDIV'), type = 'character', default = NULL, help = 'path to list of individuals'),
                   make_option(c('--POS'), type = 'character', default = NULL, help = 'path to site positions'),
                   make_option(c('--COUNT'), type = 'character', default = NULL, help = 'path to allele counts'),
                   make_option(c('--POS'), type = 'character', default = NULL, help = 'prefix of path to all three files with .indiv , .pos, .count extensions'),
                   make_option(c('--OUT'), type = 'character', default = NULL, help = 'output file path'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

outpath = opt$OUT
indiv_file = opt$INDIV
pos_file = opt$POS
count_file = opt$COUNT

if (! is.null(opt$PREFIX)){
  prefix = opt$PREFIX
  if (is.null(indiv_file)){ indiv_file = paste0(prefix,'.indiv')} 
  if (is.null(pos_file )){ pos_file = paste0(prefix,'.pos')}
  if (is.null(count_file)){ count_file = paste0(prefix,'.count')}
  if (is.null(outpath)) {outpath = paste0(prefix,'_processed.txt')}
} 
# else {
#   indiv_file = opt$INDIV
#   pos_file = opt$POS
#   count_file = opt$COUNT
# }
if (is.null(outpath) ){
  stop("No output file path provided. Without prefix, this cannot be inferred.")
}



# Load the list of individuals, the site position info, and the allele counts
inds = read.table(indiv_file, header = F, stringsAsFactors = F)[, 1]
pos = read.table(pos_file, header = F)
genos = t(as.data.frame(fread(count_file, na.strings = c('-1'))))[-1, ]

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
