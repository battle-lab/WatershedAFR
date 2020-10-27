#!/bin/bash

# manipulate annotation files in various ways to make them easier to use downstream
# takes in GTEx annotation file from the command line
# Not gzipped, unlike previous

set -o nounset -o errexit -o pipefail

if [ $# -ne 2 ]; then
	echo "Usage: process.reference.files.sh <gtex annotation, output directory>"
	exit
fi

gtex=$1
outdir=$2

rootname=`basename $gtex | sed 's/.gtf//'`
gtexprefix=${outdir}/${rootname}

# gene bed file
cat $gtex | awk 'BEGIN{OFS="\t"}{
if (substr($1,1,1)=="#") {next};
if ($3!="gene") {next};

if (substr($1,1,1)!="c") {
  chr="chr"$1
} else {
  chr=$1
}

strand=$7
if (strand=="+") {
  start=$4-1
  end=$4
} else if (strand=="-") {
  start=$5-1
  end=$5
} else {
  print "no strand information: exiting"
  exit 1
}

name=substr($10,2,length($10)-3);
print chr,start,end,name,0,strand
}' > ${gtexprefix}.bed

# gene bed file with 10kb added before the TSS
cat ${gtexprefix}.bed | awk 'BEGIN{OFS="\t"}{
strand=$6
if (strand=="+") {
  low = $2-10000;
  if (low < 0) {low = 0};
  print $1, low, $3, $4, $5, strand
} else if (strand=="-") {
  print $1, $2, $3+10000, $4, $5, strand
} else {
  print "no strand information: exiting"
  exit 1
}

}' > ${gtexprefix}_padded10kb.bed

# genetypes
cat $gtex | awk '{
if($3=="gene" && $1 ~ /^[chr0-9]+$/){
  print substr($10,2,length($10)-3)"\t"substr($14,2,length($14)-3)
}}' > ${gtexprefix}_genetypes_autosomal.txt

# gene bed file with PC and lincRNA only from padded 10kb
grep -E 'lincRNA|protein_coding' ${gtexprefix}_genetypes_autosomal.txt > ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt
## sort files by gene
sort -k 1 ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt > tmp.txt && mv tmp.txt ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt
sort -k 4 ${gtexprefix}_padded10kb.bed > tmp.txt && mv tmp.txt ${gtexprefix}_padded10kb.bed
join --nocheck-order -j1 1 -j2 4 ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt ${gtexprefix}_padded10kb.bed | \
awk '{print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $6 "\t" $7}' | sort -V > ${gtexprefix}_padded10kb_PCandlinc_only.bed
