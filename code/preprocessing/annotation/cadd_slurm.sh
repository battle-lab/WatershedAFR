##!/bin/bash
chrom=$SLURM_ARRAY_TASK_ID
corenum=$SLURM_CPUS_PER_TASK 	
# set -o nounset -o errexit -o pipefail  <--- Note to self that this is what broke everything
source ~/.bashrc
module load R
# chrom=$1
# corenum=$2

conda init 
# . /software/apps/anaconda/5.2/python/3.6/etc/profile.d/conda.sh

datadir=$1 #/scratch/groups/abattle4/victor/WatershedAFR/data/
envloc=$2 #~/.conda/envs/cadd
cadd_loc=$3 #/scratch/groups/abattle4/jessica/RareVar_AFR/cadd/CADD-scripts
code_loc=$4 #/scratch/groups/abattle4/jessica/RareVar_AFR/WatershedAFR/code/preprocessing/annotation/

# /software/apps/anaconda/5.2/python/3.6/bin/conda activate ${envloc}
# conda activate ${envloc} <- this one worked most recently

vcf_loc=${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf
chr_vcf_loc=$datadir/rare_variants_gnomad/gene-AFR-rv.CADD.chr${chrom}.vcf
outprefix=$datadir/annotation/gene-AFR-rv.CADD.chr${chrom}

# awk -v XX="$chrom" '$1==XX' ${vcf_loc} > ${chr_vcf_loc}
# bash $cadd_loc/CADD.sh -a -g GRCh38 -c ${corenum} \
#   -o ${outprefix}_tmp.tsv.gz \
#   ${chr_vcf_loc} 2>&1 | tee ${outprefix}.log
# 


Rscript --vanilla $code_loc/cadd_annot_munge.R --cadd ${outprefix}_tmp.tsv.gz --out ${outprefix}.tsv.gz

