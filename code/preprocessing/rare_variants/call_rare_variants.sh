#!/bin/bash

### file paths
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
datadir=${rootdir}/data
rawdir=${rootdir}/raw_data

rvdir=${datadir}/rare_variants

gtex_raw=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
regions=${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed

pop_list=${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt
pop="EUR"

outliers=${datadir}/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt

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

## Filter for rare variants within 10kb +/- window around gene body of genes that are protein coding and lincRNA coding
# Save the resulting VCFs to `gtex.vcf.gz` and `1KG.padded10kb_PCandlinc_only.vcf.gz` in `${datadir}/rare_variants`

echo "*** Filter for rare variants within 10kb +/- window around gene body of genes that are protein coding and lincRNA coding"

# GTEx variants (keep SNPs only)
echo "**** Filtering GTEx variants"

gtex=${rvdir}/gtex.vcf.gz

bcftools view --regions-file $regions --types snps $gtex_raw | bcftools sort | bcftools norm --rm-dup snps \
--output $gtex --output-type z
bcftools index --tbi $gtex

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
bcftools index --tbi $_1kg

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

echo -e "CHROM\tPOS\tREF\tALT\tSAMPLE" > $indiv_at_rv
bcftools query -f'[%CHROM\t%POS0\t%REF\t%ALT\t%SAMPLE\n]' --include 'GT="alt"' $gtex_pop_rareQC \
>> $indiv_at_rv

#Create two bed files for each individual:
#- 10kb window before and after the gene in the gene-individual pair
#- all rare variants belonging to the individual
echo "**** Running bed_prep.R to create bed files for each sample"

bed_outdir=${rvdir}/rv_bed_${pop}

Rscript code/preprocessing/rare_variants/bed_prep.R \
--outliers=$outliers \
--population=$pop_list \
--regions=$regions \
--indiv_at_rv=$indiv_at_rv
--outdir=$bed_outdir

# Find overlap between gene-individual pair's region and rare variants
echo "**** Matching rare variants to their gene-individual pairs"

rv_sites=${bed_outdir}/all_rv_sites.bed

while IFS='' read -r indiv || [ -n "${indiv}" ]; do
    regions_gene_indiv=${bed_outdir}/${indiv}.gene-indiv.bed
    rv=${bed_outdir}/${indiv}.rv.bed

    bedtools intersect -wa -wb -a $rv -b $regions_gene_indiv > ${bed_outdir}/${indiv}.rv_sites.bed

done < ${bed_outdir}/indiv_list

cat $(ls ${bed_outdir}/*.rv_sites.bed) > $rv_sites

#Get list of rare variants per each gene-individual pair
echo "Running gene_indiv_rare_variants.R"

Rscript code/preprocessing/rare_variants/gene_indiv_rare_variants.R \
--rv_sites=$rv_sites \
--popname=${pop} \
--outdir=$rvdir

echo "*** Done"
