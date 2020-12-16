#!/usr/bin/env Rscript

## Rscript to process TPM and read count matrices so that columns match covariate files.
## Prepare matrices for input to PEER.

## Load required packages
library(data.table)
library(stringr)

##------------- FUNCTIONS

## For each tissue, read in the TPM and read files from the given directory.
## Subset and reorder columns of TPM file to match the given covariates.
## Filter for genes with > 20% individuals with TPM > 0.1 and read count > 6.
## Log2(tpm + 2) transform the data, then z-transform.
## Finally, output the transformed matrix to a new file.
ztrans.tissue = function(tissue, dir, covs, read.filt = 6, tpm.filt = 0.1) {
    tpm = fread(paste0(dir, '/', tissue, '.tpm.txt'))
    reads = fread(paste0(dir, '/', tissue, '.reads.txt'))
    ## sanity checks
    stopifnot(sum(colnames(tpm) != colnames(reads)) == 0)
    stopifnot(sum(tpm$Gene != reads$Gene) == 0)

    covariates.subset = covs$SUBJID[covs$SUBJID %in% colnames(tpm)]
    genes = tpm$Gene
    tpm = tpm[, covariates.subset, with = F]
    reads = reads[, covariates.subset, with = F]
    ind.filt = round(0.2*ncol(tpm))
    zero.filt = round(0.05*ncol(tpm))
    indices.keep = rowSums(tpm > tpm.filt & reads > read.filt) >= ind.filt
    zero.keep = rowSums(tpm > 0) >= zero.filt
    indices.keep = indices.keep[which(indices.keep %in% zero.keep)]
    tpm = tpm[indices.keep, ]
    tpm.out = scale(t(log2(tpm + 2))) 
    colnames(tpm.out) = genes[indices.keep]
    write.table(tpm.out, paste0(dir, '/', tissue, '.log2.ztrans.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
    return()
}

##------------- MAIN

dir = '/scratch/groups/abattle4/victor/WatershedAFR/data/data_prep'
peer.dir = paste0(dir,'/PEER')
cov.file = paste0(dir,'/gtex_v8_eQTL_covariates.txt')
map.file = paste0(dir, '/gtex_v8_samples_tissues.txt')

## Read in list of tissues and keep those with more than 50 samples
tissue.list = read.table(map.file, header = F, stringsAsFactors = F)[,2]
tissues = names(table(tissue.list)[table(tissue.list) > 50]) 
tissues = tissues[tissues != 'Cells_Leukemia_cell_line_CML'] # exclude K-652 samples

## Read in covariates (PC's 1-5 and sex)
covariates = read.table(cov.file, header = T)

sapply(tissues, ztrans.tissue, dir = peer.dir, covs = covariates)

print(paste("log2(tpm + 2) transformed data saved to", peer.dir))
print(paste(length(tissues),"out of",length(unique(tissue.list)), "tissues had more than 50"))
print("The remaining tissues were not log2-transformed")