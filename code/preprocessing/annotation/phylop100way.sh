#!/bin/bash

# Obtain conservation scores for rare variants from PhyloP 100way scores in bedgraph format

## read in arguments
while getopts b:r: flag
do
    case "${flag}" in
        b) phylop=${OPTARG};; #/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/hg38.phyloP100way.bedGraph
        r) rarevariantsbed=${OPTARG};; #/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.bed
    esac
done

# defining directories
annodir=$(dirname $phylop)
outdir=$(dirname $rarevariantsbed)

# split by chromosome
for i in {1..22}
do
  chrom=chr$i
  bed=${annodir}/hg38.phyloP100way.${chrom}.bed
  bedsorted=${annodir}/hg38.phyloP100way.${chrom}.sorted.bed
 echo "Extracting $chrom"
 grep -w $chrom $phylop > $bed

 # sort on each split bed file
 echo "Sorting $chrom"
 sort --parallel=4 -T $annodir -k1,1 -k2,2n $bed > $bedsorted

done

# combine by chromosomes
sortedlist=${annodir}/sorted.bed.files.list
combinedbed=${annodir}/hg38.phyloP100way.sorted.bed
echo "Combining sorted bed files by chromosome"
ls ${annodir}/hg38.phyloP100way.*.sorted.bed | sort -V > $sortedlist
cat $(grep -v '^#' $sortedlist) > $combinedbed

# compress with bigzip
bgzippedcombined=${combinedbed}.gz
echo "bigzipping"
bgzip -c $combinedbed > $bgzippedcombined
# index with tabix
echo "indexing"
tabix -p bed $bgzippedcombined

# query phylop scores for the rare variants
echo "querying phylop 100 way scores"
conda activate watershed
python phylop100way_query.py --phylop $bgzippedcombined --rarevariant $rarevariantsbed