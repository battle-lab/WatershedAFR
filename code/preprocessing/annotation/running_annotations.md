# Annotating rare variants
```bash
datadir=/scratch/groups/abattle4/victor/WatershedAFR/data/
envloc=/home-2/jbonnie1@jhu.edu/.conda/envs/cadd
cadd_loc=/scratch/groups/abattle4/jessica/RareVar_AFR/cadd/CADD-scripts

```



## CADD (Combined Annotation Dependent Depletion)


Convert rare variants to .vcf format
```bash
Rscript rare_variants_to_vcf.R --RV ${datadir}/rare_variants_gnomad/gene-AFR-rv.txt
ml anaconda
conda activate $envloc
bash $cadd_loc/CADD.sh -a -g GRCh38 -o $datadir/annotation/gene-AFR-rv.CADD.tsv.gz ${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf

```


## VEP (Variant Annotation Predictor)

## LoF (Loss of Function) from loftee

## UCSC Conservation scores (Gerp, PhyloP, Phastcons)