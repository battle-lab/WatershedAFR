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
If you are planning on getting LoF annotations, skip this section. Running loftee will get both VEP and LoF annotations. Follow instructions to install and setup VEP [here](https://github.com/battle-lab/battle-lab-guide/blob/master/marcc_guide/software/VEP-singularity-docker.md)

### Run VEP from singularity container
```bash
vcf_input=${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf
output=${datadir}/rare_variants_gnomad/gene-AFR-rv.vep.vcf

singularity exec ensembl-vep.simg vep -i $vcf_input --format vcf --output_file $output --vcf --cache
```
## LoF (Loss of Function) from loftee
Running the loftee plugin from VEP will return both VEP and LoF annotations. Follow instructions to install and setup VEP with loftee [here](https://github.com/battle-lab/battle-lab-guide/blob/master/marcc_guide/software/VEP-singularity-docker.md)

### Run VEP with loftee
```bash
vcf_input=${datadir}/rare_variants_gnomad/gene-AFR-rv.CADD.vcf
output=${datadir}/rare_variants_gnomad/gene-AFR-rv.vep.loftee.vcf

singularity exec ensembl-vep.simg vep -i $vcf_input --format vcf --output_file $output --vcf --cache --plugin LoF,loftee_path:$HOME/.vep/Plugins/loftee/ --dir_plugins $HOME/.vep/Plugins/loftee/
```
## UCSC Conservation scores (PhyloP 100way)
PhyloP 100way scores can be downloaded from [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/). The scores are in bigwig format (.bw). They can be converted to the bedGraph format by using [bigWigToBedGraph](http://hgdownload.soe.ucsc.edu/admin/exe/)

### Downloading PhyloP as bigwig
```bash
wget -c -P ${rawdir} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
```

### Convert to bedGraph and use tabix for querying
```bash
bigwig=${rawdir}/hg38.phyloP100way.bw
annodir=${datadir}/annotation
bedgraph=${annodir}/hg38.phyloP100way.bedGraph

# convert to bedgraph
bigWigToBedGraph $bigwig $bedgraph

# split by chromosome
for i in {1..22}
do
  chrom=chr$i
  bed=${annodir}/hg38.phyloP100way.${chrom}.bed
  bedsorted=${annodir}/hg38.phyloP100way.${chrom}.sorted.bed
  bgzipped=${annodir}/hg38.phyloP100way.${chrom}.sorted.bed.gz
  echo "Extracting $chrom"
  grep -w $chrom $bedgraph > $bed
  
  # sort on each split bed file
  echo "Sorting $chrom"
  sort --parallel=4 -T $annodir -k1,1 -k2,2n $bed > $bedsorted
  
  # compress with bgzip
  echo "bgzipping $chrom"
  bgzip -c $bedsorted > $bgzipped
  
  # index with tabix
  echo "indexing $chrom"
  tabix -p bed $bgzipped
  
done

```
