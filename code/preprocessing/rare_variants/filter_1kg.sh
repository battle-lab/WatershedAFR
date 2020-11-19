#!/bin/bash

# file paths
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
datadir=${rootdir}/data
rawdir=${rootdir}/raw_data

rvdir=${datadir}/rare_variants

gtex_raw=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
regions=${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed

# 1KG variants
echo "**** Filtering 1KG variants"

_1kg=${rvdir}/1KG.padded10kb_PCandlinc_only.vcf.gz
outdir_1kg=${rvdir}/1KG/genes_padded10kb_PCandlinc_only

for i in {1..22}
do
  infile=${rawdir}/1KG/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
  outfile=`basename $infile | sed 's/.vcf.gz//'`
  outfile=${outdir_1kg}/${outfile}.genes_padded10kb_PCandlinc_only.vcf.gz
  bcftools view --drop-genotypes --regions-file $regions --output-file $outfile --output-type z $infile
done

# concatenate chromosome VCFs from 1KG into one VCF
ls -1 ${outdir_1kg}/*.vcf.gz | sort -V > ${outdir_1kg}/chrom_list.txt

bcftools concat --file-list ${outdir_1kg}/chrom_list.txt | bcftools sort | \
bcftools norm --rm-dup snps --output $_1kg --output-type z


