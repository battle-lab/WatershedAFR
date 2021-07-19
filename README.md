# WatershedAFR
Watershed implementation for African American samples from GTEx


Parts of this code were originally written by:
* Emily Tsang (https://github.com/ektsang)
* Joe Davis (https://github.com/joed3)
* Yungil Kim (https://github.com/ipw012)
* Nicole Ferraro (https://github.com/nmferraro5)
* Ben Strober (https://github.com/BennyStrobes)

## Setup

### Dependencies

Install python and dependencies from `environment.yml`
```
conda env create -f environment.yml
```


### Directory structure
```
.
├── data
│   ├── data_prep
│   ├── figures
│   ├── models
│   ├── outlier_calling
│   └── rare_variants
│   └── RIVER
├── raw_data
│   ├── 1KG
│   └── GTEx
└── WatershedAFR
    ├── code
    │   ├── analysis
    │   ├── figures
    │   ├── preprocessing
    │   │   ├── data_prep
    │   │   ├── outlier_calling
    │   │   └── rare_variants
    │   └── Watershed
    ├── data_sources.md
    ├── environment.yml
    ├── LICENSE
    ├── pipelines
    └── README.md
```

No absolute paths are hard coded into scripts. `data` and `raw_data` directories are assigned local variable names.

```bash
# Defining root, data, and raw data directories
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
datadir=${rootdir}/data
rawdir=${rootdir}/raw_data

## Code location for executable Watershed Repo
watershed_dir=/scratch/groups/abattle4/jessica/RareVar_AFR/Watershed
```

```
# Making directories
mkdir ${rootdir}

mkdir ${datadir}
mkdir ${datadir}/data_prep
mkdir ${datadir}/figures
mkdir ${datadir}/models
mkdir ${datadir}/outlier_calling
mkdir ${datadir}/rare_variants
mkdir ${datadir}/rare_variants/1KG
mkdir ${datadir}/rare_variants/1KG/genes_padded10kb_PCandlinc_only
mkdir ${datadir}/rare_variants/rv_bed_EUR
mkdir ${datadir}/rare_variants/rv_bed_AFR


mkdir ${rawdir}
mkdir ${rawdir}/1KG
mkdir ${rawdir}/GTEx

# Clone the repo
git clone https://github.com/battle-lab/WatershedAFR.git ${rootdir}/WatershedAFR

```

### Raw data
See [data_sources.md](https://github.com/battle-lab/WatershedAFR/blob/master/data_sources.md)

## Preprocessing pipeline
Goal: Preprocess raw data to be used as input for training Watershed models.

#### Make lists of individuals
Make list of all individuals from GTEx v8. Make list of all individuals with reported African American ancestry and European ancestry. Requires `bcftools`
* `${datadir}/data_prep/gtex_v8_individuals_all.txt` - All 948 individuals from GTEx v8
* `${datadir}/data_prep/gtex_v8_wgs_individuals.txt` - All 838 individuals from GTEx v8 with WGS data
* `${datadir}/data_prep/gtex_v8_individuals_AFR.txt` - 121 African American individuals from GTEx v8
* `${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt` - 103 African American individuals from GTEx v8 with WGS data
* `${datadir}/data_prep/gtex_v8_individuals_EUR.txt` - 804 European individuals from GTEx v8
* `${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt` - 714 European individuals from GTEx v8 with WGS 

```bash
# All individuals
## list of all individuals from GTEx v8 (some do not have WGS data)
cat ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F "\t" 'NR>1{print $1}' \
> ${datadir}/data_prep/gtex_v8_individuals_all.txt


## list of all individuals from GTEx v8 with WGS data
bcftools query -l ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
> ${datadir}/data_prep/gtex_v8_wgs_individuals_all.txt

# African American individuals
## list of all African American individuals from GTEx v8 (some do not have WGS data)
cat ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F '\t' '{if ($5 == 2) print $1;}' \
> ${datadir}/data_prep/gtex_v8_individuals_AFR.txt

## list of all African American individuals from GTEx v8 with WGS data
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' ${datadir}/data_prep/gtex_v8_wgs_individuals_all.txt \
${datadir}/data_prep/gtex_v8_individuals_AFR.txt > ${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt

# European individuals
# list of all European individuals from GTEx v8 (some do not have WGS data)
cat ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F '\t' '{if ($5 == 3) print $1;}' \
> ${datadir}/data_prep/gtex_v8_individuals_EUR.txt

# list of all European individuals from GTEx v8 with WGS data
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' ${datadir}/data_prep/gtex_v8_wgs_individuals_all.txt \
${datadir}/data_prep/gtex_v8_individuals_EUR.txt > ${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt
```

### Expression data correction and normalization

#### Generate tpm and read count matrices from GTEx V8

Generate file mapping sample identifiers to tissues. Restrict to samples that pass RNA-seq QC (marked as RNASEQ in SMAFRZE column).
```bash
cat ${rawdir}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | tail -n+2 | cut -f1,7,17 | \
  sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | \
  awk '$3=="RNASEQ" {print $1"\t"$2}' | sort -k 1 > ${datadir}/data_prep/gtex_v8_samples_tissues.txt
```

Splits the combined TPM file and read counts file by tissue.
```bash
# split TPM
GTEX_tpm=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
OUT=${datadir}/data_prep/PEER
SAMPLE=${datadir}/data_prep/gtex_v8_samples_tissues.txt
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX_tpm --out $OUT --sample $SAMPLE --end .tpm.txt

# split read counts
GTEX_reads=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX_reads --out $OUT --sample $SAMPLE --end .reads.txt

```

Population specific TPM and read counts
```bash
# Subset `gtex_v8_samples_tissues.txt` by African and European populations
Rscript code/preprocessing/data_prep/population_sample_tissue_map.R \
  --MAP ${datadir}/data_prep/gtex_v8_samples_tissues.txt \
  --POP.LIST ${datadir}/data_prep/gtex_v8_individuals_AFR.txt \
  --OUT ${datadir}/data_prep/gtex_v8_samples_tissues_AFR.txt
  
Rscript code/preprocessing/data_prep/population_sample_tissue_map.R \
  --MAP ${datadir}/data_prep/gtex_v8_samples_tissues.txt \
  --POP.LIST ${datadir}/data_prep/gtex_v8_individuals_EUR.txt \
  --OUT ${datadir}/data_prep/gtex_v8_samples_tissues_EUR.txt

# split TMP
GTEX_tpm=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

OUT_AFR=${datadir}/data_prep/PEER_AFR
SAMPLE_AFR=${datadir}/data_prep/gtex_v8_samples_tissues_AFR.txt
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX_tpm --out $OUT_AFR --sample $SAMPLE_AFR --end .tpm.txt

OUT_EUR=${datadir}/data_prep/PEER_EUR
SAMPLE_EUR=${datadir}/data_prep/gtex_v8_samples_tissues_EUR.txt
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX_tpm --out $OUT_EUR --sample $SAMPLE_EUR --end .tpm.txt

# split read counts
GTEX_reads=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX_reads --out $OUT_AFR --sample $SAMPLE_AFR --end .reads.txt
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX_reads --out $OUT_EUR --sample $SAMPLE_EUR --end .reads.txt

```

#### Generate PEER corrected data (includes non-EAs) with covariates removed

Build covariate matrix with PC's 1 - 5 and sex from the eQTL covariates by tissue.
```bash
Rscript code/preprocessing/data_prep/combine_covariates_across_tissues.R \
  --COV ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL_covariates \
  --OUT ${datadir}/data_prep
```

For each tissue, filter for genes with > 20% individuals with TPM > 0.1 and read count > 6, Log2(tpm + 2) transfrom the data, and then z-transform.
```bash
# All populations
Rscript code/preprocessing/data_prep/preprocess_expr.R \
  --COV ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  --MAP ${datadir}/data_prep/gtex_v8_samples_tissues.txt \
  --PEER ${datadir}/data_prep/PEER

# African
Rscript code/preprocessing/data_prep/preprocess_expr.R \
  --COV ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  --MAP ${datadir}/data_prep/gtex_v8_samples_tissues_AFR.txt \
  --MIN.SAMPLE 11 \
  --PEER ${datadir}/data_prep/PEER_AFR
  
# European
Rscript code/preprocessing/data_prep/preprocess_expr.R \
  --COV ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  --MAP ${datadir}/data_prep/gtex_v8_samples_tissues_EUR.txt \
  --PEER ${datadir}/data_prep/PEER_EUR
```

Rename eQTL files from "cervical_c-1" to "cervical_c1" for consistency
```bash
rename c-1 c1 ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL/*c-1*
```

Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles
```bash
bash code/preprocessing/data_prep/get_eqtl_genotypes.sh
```

Run PEER correction and compute residuals
```bash
## All populations
# Compute PEER factors (use PEER dockerimage if unable to install peer natively https://hub.docker.com/r/bryancquach/peer)
bash code/preprocessing/data_prep/calculate_PEER_factors.sh \
  -p ${datadir}/data_prep/PEER \
  -e ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL \
  -c ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  -g ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
  -t ~/scratch \
  -d ~/code/PEER/peer-1.3.simg
  
# Compute residuals (does not need PEER package)
bash code/preprocessing/data_prep/calculate_PEER_residuals.sh \
  -p ${datadir}/data_prep/PEER \
  -e ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL \
  -c ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  -g ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
  -t ~/scratch


## African
# Compute PEER factors (use PEER dockerimage if unable to install peer natively https://hub.docker.com/r/bryancquach/peer)
bash code/preprocessing/data_prep/calculate_PEER_factors.sh \
  -p ${datadir}/data_prep/PEER_AFR \
  -e ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL \
  -c ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  -g ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
  -t ~/scratch \
  -d ~/code/PEER/peer-1.3.simg
  
# Compute residuals (does not need PEER package)
bash code/preprocessing/data_prep/calculate_PEER_residuals.sh \
  -p ${datadir}/data_prep/PEER_AFR \
  -e ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL \
  -c ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  -g ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
  -t ~/scratch


## European
# Compute PEER factors (use PEER dockerimage if unable to install peer natively https://hub.docker.com/r/bryancquach/peer)
bash code/preprocessing/data_prep/calculate_PEER_factors.sh \
  -p ${datadir}/data_prep/PEER_EUR \
  -e ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL \
  -c ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  -g ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
  -t ~/scratch \
  -d ~/code/PEER/peer-1.3.simg
  
# Compute residuals (does not need PEER package)
bash code/preprocessing/data_prep/calculate_PEER_residuals.sh \
  -p ${datadir}/data_prep/PEER_EUR \
  -e ${rawdir}/GTEx/GTEx_Analysis_v8_eQTL \
  -c ${datadir}/data_prep/gtex_v8_eQTL_covariates.txt \
  -g ${datadir}/data_prep/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
  -t ~/scratch
```

#### Combine PEER-corrected data into a single flat file

Generate files with data on what tissues are available per individual after PEER correction
```bash
# All populations
bash code/preprocessing/data_prep/get_tissue_by_individual.sh \
  -p ${datadir}/data_prep/PEER \
  -o ${datadir}/data_prep/gtex_2017-06-05_tissue_by_ind.txt \
  -t ${datadir}/data_prep/gtex_2017-06-05_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_2017-06-05_individuals_all_normalized_samples.txt

# African
bash code/preprocessing/data_prep/get_tissue_by_individual.sh \
  -p ${datadir}/data_prep/PEER_AFR \
  -o ${datadir}/data_prep/gtex_AFR_tissue_by_ind.txt \
  -t ${datadir}/data_prep/gtex_AFR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_AFR_individuals_all_normalized_samples.txt
  
# European
bash code/preprocessing/data_prep/get_tissue_by_individual.sh \
  -p ${datadir}/data_prep/PEER_EUR \
  -o ${datadir}/data_prep/gtex_EUR_tissue_by_ind.txt \
  -t ${datadir}/data_prep/gtex_EUR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_EUR_individuals_all_normalized_samples.txt

```

Combine the PEER-corrected data
```bash
# All populations
# output located in ${datadir}/data_prep/gtex_2017-06-05_normalized_expression.txt.gz
python code/preprocessing/data_prep/gather_filter_normalized_expression.py \
  -p ${datadir}/data_prep/PEER \
  -t ${datadir}/data_prep/gtex_2017-06-05_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_2017-06-05_individuals_all_normalized_samples.txt \
  -o ${datadir}/data_prep/gtex_2017-06-05_normalized_expression.txt

gzip ${datadir}/data_prep/gtex_2017-06-05_normalized_expression.txt


# African
# output located in ${datadir}/data_prep/gtex_AFR_normalized_expression.txt.gz
python code/preprocessing/data_prep/gather_filter_normalized_expression.py \
  -p ${datadir}/data_prep/PEER_AFR \
  -t ${datadir}/data_prep/gtex_AFR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_AFR_individuals_all_normalized_samples.txt \
  -o ${datadir}/data_prep/gtex_AFR_normalized_expression.txt

gzip ${datadir}/data_prep/gtex_AFR_normalized_expression.txt


# European
# output located in ${datadir}/data_prep/gtex_EUR_normalized_expression.txt.gz
python code/preprocessing/data_prep/gather_filter_normalized_expression.py \
  -p ${datadir}/data_prep/PEER_EUR \
  -t ${datadir}/data_prep/gtex_EUR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_EUR_individuals_all_normalized_samples.txt \
  -o ${datadir}/data_prep/gtex_EUR_normalized_expression.txt

gzip ${datadir}/data_prep/gtex_EUR_normalized_expression.txt
```

### Outlier calling

#### Call outliers on all individuals

Saved to `${datadir}/outlier_calling/test/gtexV8.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt`
```bash
Rscript code/preprocessing/outlier_calling/call_outliers.R \
  --Z.SCORES=${datadir}/data_prep/gtex_2017-06-05_normalized_expression.txt.gz \
  --OUT=${datadir}/outlier_calling/test/gtexV8.outlier.controls.v8ciseQTLs.txt \
  --N.PHEN=5 --ZTHRESH=3


# Remove global outliers
Rscript code/preprocessing/outlier_calling/identify_global_outliers.R \
  --OUTLIERS=${datadir}/outlier_calling/test/gtexV8.outlier.controls.v8ciseQTLs.txt \
  --METHOD=proportion
```

#### Call outliers on African individuals  
Saved to `${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt` 

```bash
Rscript code/preprocessing/outlier_calling/call_outliers.R \
  --Z.SCORES=${datadir}/data_prep/gtex_AFR_normalized_expression.txt.gz \
  --OUT=${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs.txt \
  --POP=${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt \
  --N.PHEN=5 --ZTHRESH=3

Rscript code/preprocessing/outlier_calling/identify_global_outliers.R \
  --OUTLIERS=${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs.txt \
  --METHOD=proportion
```

#### Call outliers on European individuals

Saved to `${datadir}/outlier_calling/EUR/gtexV8.EUR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt`
```bash
Rscript code/preprocessing/outlier_calling/call_outliers.R \
  --Z.SCORES=${datadir}/data_prep/gtex_EUR_normalized_expression.txt.gz \
  --OUT=${datadir}/outlier_calling/EUR/gtexV8.EUR.outlier.controls.v8ciseQTLs.txt \
  --POP=${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt \
  --N.PHEN=5 --ZTHRESH=3


# Remove global outliers
Rscript code/preprocessing/outlier_calling/identify_global_outliers.R \
  --OUTLIERS=${datadir}/outlier_calling/EUR/gtexV8.EUR.outlier.controls.v8ciseQTLs.txt \
  --METHOD=proportion

```

GTEx v8 outliers from Watershed paper are located in 
`/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/outlier_calls/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt`
Save as `${datadir}/outlier_calling/WATERSHED.gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt`

### Prepare gene annotations from gencode

Manipulate annotation files in various ways to make them easier to use downstream.
Outputs
* `${datadir}/data_prep/gencode.v26.GRCh38.genes.bed` - Gene bed file
* `${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb.bed` - Gene bed file with 10kb added on either side
* `${datadir}/data_prep/gencode.v26.GRCh38.genes_genetypes_autosomal.txt` - Genetypes
* `${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed` - Gene bed file with protein coding and lincRNA coding genes padded by 10kb
```
bash ${rootdir}/WatershedAFR/code/preprocessing/data_prep/process_reference_files.sh \
${rawdir}/GTEx/gencode.v26.GRCh38.genes.gtf \
${datadir}/data_prep
```

### Finding rare variants
Requires `bcftools` and `bedtools`

#### European Individuals

Filter with gnomAD. Rare variants file will be saved to `{datadir}/rare_variants_gnomad/gene-EUR-rv.txt`  
arguments to `find_rare_variants_gnomad.sh`:  
-g: raw gtex vcf  
-r: bed file with protein coding and lnc rna coding regions (generated during outlier calling earlier)  
-f gnomad raw file  
-l list of individuals in population  
-p: prefix/suffix also used to specify population. Select from "EUR" or "AFR"  

```bash
bash code/preprocessing/rare_variants/find_rare_variants_gnomad.sh \
-d ${datadir}/rare_variants_gnomad \
-g ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
-r ${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-f /work-zfs/abattle4/lab_data/gnomAD_v2.1/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz \
-l ${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt \
-p "EUR" 

```

#### African Individuals
Filter with gnomAD. Rare variants file will be saved to `{datadir}/rare_variants_gnomad/gene-AFR-rv.txt`
```bash
bash code/preprocessing/rare_variants/find_rare_variants_gnomad.sh \
-d ${datadir}/rare_variants_gnomad \
-g ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
-r ${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-f /work-zfs/abattle4/lab_data/gnomAD_v2.1/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz \
-l ${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt \
-p "AFR" 

```


## Training Watershed models  
Prep files for RIVER and watershed.

```bash
Rscript code/preprocessing/data_prep/dataprep_watershed.R \
--ZSCORES ${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt \
--ANNOT /work-zfs/abattle4/bstrober/random_projects/feature_generation_for_victor_and_jessica/river_input_gene_level.txt \
--OUT ${datadir}/data_prep/RIVER/river_input_v8_african_all_07-19-2021.txt
```

look at RIVER folder in RIVER repo for RIVER.Rmd

## Analysis and Figures
look at RIVER folder in RIVER repo for enrichment_AFR.Rmd


