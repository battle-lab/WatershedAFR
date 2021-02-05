#!/bin/bash

### file paths
# rootdir=/scratch/groups/abattle4/victor/WatershedAFR
# datadir=${rootdir}/data
# rawdir=${rootdir}/raw_data
# 
# rvdir=${datadir}/rare_variants
# 
# gtex_raw=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
# regions=${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed
# 
# pop_list=${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt
# pop="EUR"
# 

## read in arguments
while getopts d:g:r:l:p:o: flag
do
    case "${flag}" in
        d) rvdir=${OPTARG};;
        g) gtex_raw=${OPTARG};;
        r) regions=${OPTARG};;
        l) pop_list=${OPTARG};;
        p) pop=${OPTARG};;
    esac
done
echo $rvdir
echo $gtex_raw
echo $regions
echo $pop_list;
echo $pop

## Check if bcftools and bedtools are available. If not, exit the script
if ! command -v bcftools &> /dev/null
then
    echo "bcftools could not be found"
    exit
fi

if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found"
    exit
fi

# Make output directory
mkdir -p $rvdir

## Filter for rare variants within 10kb +/- window around gene body of genes that are protein coding and lincRNA coding
# Save the resulting VCFs to `gtex.vcf.gz` and `1KG.padded10kb_PCandlinc_only.vcf.gz` in `${datadir}/rare_variants`

echo "*** Filter for rare variants within 10kb +/- window around gene body of genes that are protein coding and lincRNA coding"

# GTEx variants (keep SNPs only)
gtex=${rvdir}/gtex.vcf.gz

if [ -f "$gtex" ]; then
        echo "**** Filtered GTEx VCF already exists"
else
        bash code/preprocessing/rare_variants/filter_gtex.sh
        -d $rvdir \
        -g $gtex_raw \
        -r $regions
fi

if [ -f "$gtex.tbi" ]; then
        echo "**** Filtered GTEx VCF already indexed"
else
        echo "**** Indexing GTEx VCF"
        bcftools index --tbi $gtex
fi


# 1KG variants
_1kg=${rvdir}/1KG.padded10kb_PCandlinc_only.vcf.gz

if [ -f "$_1kg" ]; then
        echo "**** Filtered 1KG VCF already exists"
else
        bash code/preprocessing/rare_variants/filter_1kg.sh \

fi


if [ -f "$_1kg.tbi" ]; then
        echo "**** Filtered 1KG VCF already indexed"
else
        echo "**** Indexing 1KG VCF"
        bcftools index --tbi $_1kg
fi

echo "*** Done"


## Subset GTEx VCF by population
# also recomputes allele frequencies within that population
# save the resulting VCFs as `gtex_${pop}.vcf.gz` in `{datadir}/rare_variants`
echo "*** Subset GTEx VCF by population"

gtex_pop=${rvdir}/gtex_${pop}.vcf.gz

# bcftools +fill-tags recomputes AF after we remove samples
bcftools view --samples-file $pop_list $gtex | bcftools +fill-tags --output $gtex_pop --output-type z -- -t AF

echo "*** Done"

## Filter GTEx VCF for rare variants (MAF < 0.01)
# Save the resulting VCFs as `gtex_${pop}_rare.vcf.gz` in `${datadir}/rare_variants`
echo "*** Filter GTEx VCF for rare variants (MAF < 0.01)"

gtex_pop_rare=${rvdir}/gtex_${pop}_rare.vcf.gz

bcftools view --include 'AF<0.01 & AF>0' --output-file $gtex_pop_rare --output-type z $gtex_pop

echo "*** Done"

## Confirm rarity of GTEx variants in 1KG
# Quality control to check that rare variants in GTEx are also rare in 1KG population.
# Save the resulting VCFs to `gtex_${pop}_rare.QC.vcf.gz` in ${datadir}/rare_variants
echo "*** Confirm rarity of GTEx variants in 1KG"

gtex_rare_bed=${rvdir}/gtex_${pop}_rare.bed
_1kg_common=${rvdir}/1KG.${pop}_common.vcf.gz
_1kg_common_bed=${rvdir}/1KG.${pop}_common.bed
gtex_pop_rareQC=${rvdir}/gtex_${pop}_rare.QC.vcf.gz

# Filter 1KG for variants in GTEx that have AF >= 0.01
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $gtex_pop_rare > $gtex_rare_bed
bcftools view --regions-file $gtex_rare_bed \
--include 'INFO/EUR_AF >= 0.01' --output-file $_1kg_common --output-type z $_1kg

bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $_1kg_common > $_1kg_common_bed

# Remove common variants from GTEx
bedtools intersect -v -a $gtex_pop_rare -b $_1kg_common_bed -header | \
bcftools convert --output $gtex_pop_rareQC --output-type z

echo "*** Done"

## Get list of rare variants per each gene-individual pair
# Save the resulting gene-individual-variant files to `gene-${pop}-rv.txt` in `${datadir}/rare_variants`
echo "*** Get list of rare variants per each gene-individual pair"

# Make list of samples with each rare variant. Position is 0-based like the start position in bed file format
echo "**** List samples that have each rare variant"

indiv_at_rv=${rvdir}/gtex_${pop}_rare.QC.indiv.txt

bcftools query -f'[%CHROM\t%POS0\t%END\t%REF\t%ALT\t%SAMPLE\n]' --include 'GT="alt"' $gtex_pop_rareQC \
> $indiv_at_rv


# Use bedtools to intersect list of rare variants with protein and lincRNA coding genes padded by 10kb around gene body
echo "**** Subset rare variants that lie within 10kb around protein coding and lincRNA coding genes"
rv_sites_raw=${rvdir}/gene-${pop}-rv.raw.txt
bedtools intersect -wa -wb -a $indiv_at_rv -b $regions > $rv_sites_raw

#Get list of rare variants per each gene-individual pair
echo "**** Building file with list of rare variants per each gene-individual pair"

Rscript code/preprocessing/rare_variants/gene_indiv_rare_variants.R \
--rv_sites=$rv_sites_raw \
--popname=$pop \
--outdir=$rvdir

echo "*** Done"
date +"%r"