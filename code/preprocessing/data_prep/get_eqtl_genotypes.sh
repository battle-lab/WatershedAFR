#!/bin/bash

rawdir='/scratch/groups/abattle4/victor/WatershedAFR/raw_data'
datadir='/scratch/groups/abattle4/victor/WatershedAFR/data'

### Check if vcftools and bedtools are available. If not, exit the script
if ! command -v vcftools &> /dev/null
then
    echo "vcftools could not be found"
    exit
fi

if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found"
    exit
fi

### Generate list of top eQTLs for each gene in each tissue
eQTL_bed=${datadir}/data_prep/gtex_2017-06-05_v8_cis_eQTLs.bed

zcat ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL/*.v8.egenes.txt.gz | \
        cut -f14-17 | grep -v variant_pos | sed 's/^X/23/g' | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' | sort -k 1,1 -k2,2n | \
        sed 's/^23/X/g' | uniq > $eQTL_bed

wait

### Extract these sites from the GTEx v8 VCF using bcftools
GTEX_WGSv8=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
eQTL_vcf=${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf.gz

bcftools view --regions-file $eQTL_bed $GTEX_WGSv8 | bcftools sort --output-file $eQTL_vcf --output-type z
bcftools index --tbi $eQTL_vcf

### Convert the cis-eQTL genotypes in VCF format to the number of alternate alleles using VCFTools
vcftools --gzvcf $eQTL_vcf --out ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs --012 --maf 0.01

Rscript process_gtex_v8_cis_eqtl_genotypes.R
