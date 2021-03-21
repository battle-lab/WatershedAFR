#!/bin/bash

## read in arguments
while getopts d:g:r: flag
do
    case "${flag}" in
        d) rvdir=${OPTARG};;
        g) gtex_raw=${OPTARG};;
        r) regions=${OPTARG};;
    esac
done


# GTEx variants (keep SNPs only)
echo "**** Filtering GTEx variants"

gtex=${rvdir}/gtex.vcf.gz

bcftools view --regions-file $regions --types snps $gtex_raw | bcftools sort | bcftools norm --rm-dup snps \
        --output $gtex --output-type z


