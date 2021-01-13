#!/usr/bin/env Rscript

## Script to call outliers given a file with Z-scores.
## The two first columns have the gene and tissue or phenotype, and the gene columne is named "Gene".
## Subsequent columns have Z-score data for each individual and the individual ID is the column name.

## Load required packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(optparse)
# library(reshape2)
library(foreach)
library(doMC)
# library(robustbase)

## Register parallel backend
registerDoMC(cores = 2)

#--------------- FUNCTIONS


#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data in bed format with the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers <- function(data, nphen, zthresh) {
    cat('Melting...')
    data.melted = melt(data, id.vars = names(data)[1:2], variable.name = 'Ind', value.name = 'Z')
    colnames(data.melted)[3:4] = c('Ind', 'Z')
    rm(data)
    
    
    cat('Computing median Z scores across tissues...')
    test.stats = data.melted %>% group_by(Ind, Gene) %>%
        summarise(MedZ = median(Z, na.rm = T),
                  Df = sum(!is.na(Z)))
    
    # pick outliers
    cat('Restricting to individuals with at least 5 tissues...')
    test.stats = test.stats %>% filter(Df >= nphen) %>%
      group_by(Gene) %>%
      mutate(N = n())
    cat('Picking outliers...')
    outliers = test.stats %>% arrange(desc(abs(MedZ))) %>%
      mutate(Y = ifelse(abs(MedZ) > zthresh, 'outlier', 'control')) %>%
      ungroup() %>% select(Ind, Gene, N, Df, MedZ, Y) %>% as_tibble()
    
    return(outliers)
}

#### Function to write the outlier information to file.
## Input: data frame with data to write and output filename
## Output: None, writes directly to file.
write.outliers <- function(outliers, filename) {
    write.table(outliers, filename, sep = '\t', col.names = T, row.names = F, quote = F)
}

#--------------- MAIN

##-- Read command line arguments
option_list = list(make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'),
	make_option(c('--OUT'), type = 'character', default = NULL, help = 'path to the outlier file output'),
	make_option(c('--GLOBAL'), type = 'character', default = NA, help = 'global outlier file'),
	make_option(c('--N.PHEN'), type = 'numeric', default = 5, help = 'number of observed phenotypes required to test for outlier'),
	make_option(c('--ZTHRESH'), type = 'numeric', default = 3, help = 'threshold for abs(MedZ) outliers'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

nphen = opt$N.PHEN
zthresh = opt$ZTHRESH
outFile = opt$OUT
globalFile = opt$GLOBAL
# popFile = opt$POP

##-- Analysis

## Read in the normalized data
cat('Reading normalized expression data')
data = as.data.frame(fread(opt$Z.SCORES, fill = TRUE))
# data = as.data.frame(fread('/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep/norm_expr_test.txt.gz', fill = TRUE))


##-- Call outliers using median Z-score
# Filtering global outliers
if (!is.na(globalFile)) {
    goutliers = fread(globalFile, header=F)
    rinds = which(colnames(data) %in% goutliers$V1)
    data = data[,-1*rinds]
}

# Filtering for individuals in a specified population

## Median Z-score
outliers.medz = call.outliers(data, nphen, zthresh)
write.outliers(outliers.medz, outFile)
