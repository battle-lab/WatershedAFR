#!/bin/bash

## read in arguments
while getopts d:f:r: flag
do
    case "${flag}" in
        d) rvdir=${OPTARG};;
        f) gnomad_raw=${OPTARG};;
        r) regions=${OPTARG};;
    esac
done


#gnomAD variants (keep SNPs only)
echo "**** Filtering gnomAD variants"

gnomad=${rvdir}/gnomad.padded10kb_PCandlinc_only.vcf.gz


bcftools view --regions-file $regions --types snps -Oz -o $gnomad $gnomad_raw 
bcftools index --tbi $gnomad


