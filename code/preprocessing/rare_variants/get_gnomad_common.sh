#!/bin/bash

## read in arguments
while getopts g:r:d:p:o: flag
do
    case "${flag}" in
        g) gnomad=${OPTARG};;
        r) gtex_rare_bed=${OPTARG};;
        d) gnomad_dir=${OPTARG};;
        p) pop=${OPTARG};;
        o) gnomad_common=${OPTARG};;
    esac
done

# environment variables
# gnomad=${datadir}/rare_variants_gnomad/gnomad.padded10kb_PCandlinc_only.vcf.gz
# gtex_rare_bed=${rvdir}/gtex_${pop}_rare.bed
# gnomad_dir=${rvdir}/gnomad_by_chr
# pop=EUR
# gnomad_common=${rvdir}/gnomad.${pop}_common.vcf.gz


# split gnomAD by chromosome
echo "split gnomAD by chromosome"

mkdir -p $gnomad_dir
for i in {1..22};
do
        bcftools view -r chr$i --output-file ${gnomad_dir}/gnomad.padded10kb_PCandlinc_only_chr${i}.vcf.gz --output-type z $gnomad
        bcftools index --tbi ${gnomad_dir}/gnomad.padded10kb_PCandlinc_only_chr${i}.vcf.gz
        echo $i
done

# get common gnomAD variants
echo "Filtering for common variants in gnomAD"
get_common() {
        gnomad=$1
        chr=$(basename ${gnomad##*_} .vcf.gz)

        if [ $pop == "AFR" ]; then
                bcftools view --regions-file $gtex_rare_bed --include 'INFO/AF_afr >= 0.01' \
                        -Oz -o ${gnomad_dir}/gnomad.${pop}_common_$chr.vcf.gz $gnomad
                
        elif [ $pop == "EUR" ]; then
                bcftools view --regions-file $gtex_rare_bed --include 'INFO/AF_nfe >= 0.01' \
                        -Oz -o ${gnomad_dir}/gnomad.${pop}_common_$chr.vcf.gz $gnomad
        
        else
                echo "Population is neither AFR nor EUR. Exiting"
                exit 1
        fi 
        
        
        echo "Done: $1"

}

export -f get_common
export pop
export gtex_rare_bed
export gnomad_dir

parallel --jobs 12 get_common ::: ${gnomad_dir}/gnomad.padded10kb_PCandlinc_only_chr*.vcf.gz

# combine back into one file
by_chr_list=${gnomad_dir}/by_chr_list.txt

ls ${gnomad_dir}/gnomad.${pop}_common_chr*.vcf.gz | sort -V > $by_chr_list
bcftools concat --file-list $by_chr_list --output $gnomad_common --output-type z

echo "Done"
