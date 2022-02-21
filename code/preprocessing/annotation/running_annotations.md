# Annotating rare variants
```bash
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
rawdir=${rootdir}/raw_data
datadir=${rootdir}/data
envloc=/home-2/jbonnie1@jhu.edu/.conda/envs/cadd
cadd_loc=/scratch/groups/abattle4/jessica/RareVar_AFR/cadd/CADD-scripts
```



## CADD (Combined Annotation Dependent Depletion)


Convert rare variants to .vcf format. Getting slurm to cooperate with the environmental chicanery used by CADD-scripts was complicated, but scripts are included here.
```bash
Rscript rare_variants_to_vcf.R --RV ${datadir}/rare_variants_gnomad/gene-AFR-rv.txt

sbatch ./cadd.slurm
# bash $cadd_loc/CADD.sh -a -g GRCh38 -o $datadir/annotation/gene-AFR-rv.CADD.tsv.gz ${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf

```


## VEP (Variant Annotation Predictor)
If you are planning on getting LoF annotations, skip this section. Running LOFTEE will get both VEP and LoF annotations. Follow instructions to install and setup VEP [here](https://github.com/battle-lab/battle-lab-guide/blob/master/marcc_guide/software/VEP-singularity-docker.md)

### Run VEP from singularity container
```bash
vcf_input=${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf
# output=${datadir}/rare_variants_gnomad/gene-AFR-rv.vep.vcf
output=${datadir}/rare_variants_gnomad/gene-AFR-rv.vep.txt

# singularity exec ensembl-vep.simg vep -i $vcf_input --format vcf --output_file $output --vcf --cache
singularity exec ensembl-vep.simg vep -i $vcf_input --format vcf --output_file $output --cache --show_ref_allele --regulatory
```
## LoF (Loss of Function) from LOFTEE
Running the LOFTEE plugin from VEP will return both VEP and LoF annotations. Follow instructions to install and setup VEP with LOFTEE [here](https://github.com/battle-lab/battle-lab-guide/blob/master/marcc_guide/software/VEP-singularity-docker.md)

### Run VEP with LOFTEE
```bash
vcf_input=${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf
output=${datadir}/annotation/gene-AFR-rv.vep.loftee.vcf

singularity exec ensembl-vep.simg vep -i $vcf_input --format vcf --output_file $output --vcf --cache --regulatory --plugin LoF,loftee_path:$HOME/.vep/Plugins/loftee/ --dir_plugins $HOME/.vep/Plugins/loftee/

# Tabix
bgzip $output
tabix -p vcf $output.gz

# Parsing annotations from output VCF
vep_loftee_file=${datadir}/annotation/gene-AFR-rv.vep.loftee.vcf.gz
python parse_vep_loftee.py --anno $vep_loftee_file
```
## UCSC Conservation scores (PhyloP 100way)
PhyloP 100way scores can be downloaded from [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/). The scores are in bigwig format (.bw). They can be converted to the bedGraph format by using [bigWigToBedGraph](http://hgdownload.soe.ucsc.edu/admin/exe/)

### Downloading PhyloP as bigwig
```bash
wget -c -P ${rawdir} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
```

### Convert hg38.phyloP100way.bw to bedGraph format
```bash
bigwig=${rawdir}/hg38.phyloP100way.bw
annodir=${datadir}/annotation
bedgraph=${annodir}/hg38.phyloP100way.bedGraph
bigWigToBedGraph $bigwig $bedgraph
```

### Convert rare variants to bed format
```bash
Rscript rare_variants_to_bed.R --RV ${datadir}/rare_variants_gnomad/gene-AFR-rv.txt
```

### Use tabix to query scores
Phylop scores for the rare variants are saved to ${datadir}/rare_variants_gnomad
The following script requires `bgzip` and `tabix`, and makes use of python package `pysam`
```bash
bash phylop100way.sh -b $bedgraph -r ${datadir}/rare_variants_gnomad/gene-AFR-rv.bed
```

## Gencode annotations
distTSS - absolute distance to transcription start site
distTES - absolute distance to transcription end site
```bash

```

# Parsing annotations to be SNV level annotations


## VEP and LoF from LOFTEE
```bash
vep_loftee_file=${datadir}/annotation/gene-AFR-rv.vep.loftee.vcf.gz
python parse_vep_loftee.py --anno $vep_loftee_file
```
## UCSC Conservation scores



# Parsing annotations to be GENE level annotations  

## CADD (Also combines individual chromosome annotation files while we are at it)
```bash
Rscript --vanilla collapse_cadd.R

```


# Standardize All Annotations

## This script expects only two non-numeric non-annotation columns

```bash
annotation_file=xxx
outloc=test.tsv
Rscript --vanilla standardize_annotation.R --ANNOT ${annotation_file} --OUTFILE ${outloc}

```



## Convert from z-scores to p-values (this should not be here, but for now it is.)

```bash
zscore_file=${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt
outfile=${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved_PVAL.txt
zcol=MedZ
Rscript ../outlier_calling/convert_z_to_p.R --ZFILE ${zscore_file} --OUTFILE ${outfile} --ZCOL $zcol

```

