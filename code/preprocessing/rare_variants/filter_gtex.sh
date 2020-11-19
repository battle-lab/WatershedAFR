#!/bin/bash

# file paths
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
datadir=${rootdir}/data
rawdir=${rootdir}/raw_data

rvdir=${datadir}/rare_variants

gtex_raw=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
regions=${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed

# GTEx variants (keep SNPs only)
echo "**** Filtering GTEx variants"

gtex=${rvdir}/gtex.vcf.gz

bcftools view --regions-file $regions --types snps $gtex_raw | bcftools sort | bcftools norm --rm-dup snps \
        --output $gtex --output-type z


