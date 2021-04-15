#!/usr/bin/env Rscript

library(data.table)
library(reshape2)
library(plyr)
library(optparse)
library(doMC)

###
### Setup parallel processing
###
doMC::registerDoMC(cores=2)

############ FUNCTIONS

###
### Functions for calculating number of tissues/sample and meta-analysis z-score
###
meta.n = function(values) {
  length(values) - sum(is.na(values))
}

meta.median = function(values) {
  median(values, na.rm=T)
}

meta.analysis = function(x) {
  samples=colnames(x)[3:ncol(x)]
  y = t(x[,3:ncol(x)]) # individuals (rows) x tissues (columns)
  n = apply(y, 1, meta.n)
  m1 = apply(y, 1, meta.median)
  data.frame(sample = samples, n.tissues = n, median.z = m1)
}

############ MAIN

##-- Read command line arguments
option_list = list(make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the normalized expression data'),
                   make_option(c('--OUT'), type = 'character', default = NULL, help = 'path to the outlier file output'),
                   make_option(c('--POP'), type = 'character', default = NULL, help = 'path to list of individuals belonging to specific population'),
                   make_option(c('--N.PHEN'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier'),
                   make_option(c('--ZTHRESH'), type = 'numeric', default = 3, help = 'threshold for abs(MedZ) outliers'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


expression_file = opt$Z.SCORES
outfile = opt$OUT
pop_file = opt$POP
nphen = opt$N.PHEN
zthresh = opt$ZTHRESH

###
### Loading data
###
## Load file with normalized expression data
data = fread(expression_file, fill = TRUE, header = TRUE)

# load population list and keep individuals belonging to specified population
if (!is.null(pop_file)){
  pop_list = as.character(read.table(pop_file, header = FALSE)$V1)
  
  cols_to_keep = c(c("Tissue","Gene"),pop_list)
  data = data[,..cols_to_keep]
}

# sort the data by gene
setkey(data, Gene)

## Get sample list
individs = colnames(data)[-c(1,2)]

## Calculate meta-analysis test statistics
results = ddply(data, .(Gene), meta.analysis, .parallel = F)

## Remove samples with < n tissues
results = results[results$n.tissues >= nphen,]

# Find outliers
results$Y = ifelse(abs(results$median.z) > zthresh, 'outlier', 'control')

# Rename columns and save
colnames(results) = c("Gene","Ind","Df","MedZ","Y")
write.table(results, outfile, sep = '\t', col.names = T, row.names = F, quote = F)
