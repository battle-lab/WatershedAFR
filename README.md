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

```
# Defining root, data, and raw data directories
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
datadir=${rootdir}/data
rawdir=${rootdir}/raw_data
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

### Expression data correction and normalization


#### Generate tpm and read count matrices from GTEx V8

Generate file mapping sample identifiers to tissues. Restrict to samples that pass RNA-seq QC (marked as RNASEQ in SMAFRZE column).
```
cat ${rawdir}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | tail -n+2 | cut -f1,7,17 | \
  sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | \
  awk '$3=="RNASEQ" {print $1"\t"$2}' | sort -k 1 > ${datadir}/data_prep/gtex_v8_samples_tissues.txt
```

Split combined TPM and read counts by tissue.
```
# split TPM
OUT=${datadir}/data_prep/PEER
GTEX=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
SAMPLE=${datadir}/data_prep/gtex_v8_samples_tissues.txt
END='.tpm.txt'
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX --out $OUT --sample $SAMPLE --end $END


# split read counts
GTEX=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
END='.reads.txt'
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX --out $OUT --sample $SAMPLE --end $END

```

#### Generate PEER corrected data (includes non-EAs) with covariates removed
to-do

#### Make lists of individuals
Make list of all individuals from GTEx v8. Make list of all individuals with reported African American ancestry and European ancestry.
* `${RAREVARDIR}/reference/gtex_v8_individuals_all_normalized_samples.txt` - All 948 individuals from GTEx v8
* `${RAREVARDIR}/reference/gtex_v8_wgs_individuals.txt` - All 838 individuals from GTEx v8 with WGS data
* `${RAREVARDIR}/reference/gtex_v8_individuals_AFA.txt` - 121 African American individuals from GTEx v8
* `${RAREVARDIR}/reference/gtex_v8_wgs_individuals_AFA.txt` - 103 African American individuals from GTEx v8 with WGS data
* `${RAREVARDIR}/reference/gtex_v8_individuals_EUR.txt` - 804 European individuals from GTEx v8
* `${RAREVARDIR}/reference/gtex_v8_wgs_individuals_EUR.txt` - 714 European individuals from GTEx v8 with WGS 
```{bash gtex_individuals, eval=FALSE, cache=TRUE}
# All individuals
## list of all individuals from GTEx v8 (some do not have WGS data)
cat ${RAREVARDIR}/data/GTEx/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F "\t" 'NR>1{print $1}' \
> ${RAREVARDIR}/reference/gtex_v8_individuals_all_normalized_samples.txt

## list of all individuals from GTEx v8 with WGS data
vcf-query -l ${RAREVARDIR}/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
> ${RAREVARDIR}/reference/gtex_v8_wgs_individuals.txt

# African American individuals
## list of all African American individuals from GTEx v8 (some do not have WGS data)
cat ${RAREVARDIR}/data/GTEx/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F '\t' '{if ($5 == 2) print $1;}' \
> ${RAREVARDIR}/reference/gtex_v8_individuals_AFA.txt

## list of all African American individuals from GTEx v8 with WGS data
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' ${RAREVARDIR}/reference/gtex_v8_wgs_individuals.txt ${RAREVARDIR}/reference/gtex_v8_individuals_AFA.txt > ${RAREVARDIR}/reference/gtex_v8_wgs_individuals_AFA.txt


# European individuals
# list of all European individuals from GTEx v8 (some do not have WGS data)
cat ${RAREVARDIR}/data/GTEx/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt | awk -F '\t' '{if ($5 == 3) print $1;}' \
> ${RAREVARDIR}/reference/gtex_v8_individuals_EUR.txt

# list of all African American individuals from GTEx v8 with WGS data
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' ${RAREVARDIR}/reference/gtex_v8_wgs_individuals.txt ${RAREVARDIR}/reference/gtex_v8_individuals_EUR.txt > ${RAREVARDIR}/reference/gtex_v8_wgs_individuals_EUR.txt
```

#### Make flat files from normalized data
```{bash, eval=FALSE, cache=TRUE}
python preprocessing/data_prep/gather_filter_normalized_expression_v8.py
```

#### Prepare gene annotations from gencode
This list will be used to filter out genes for those that are protein coding and lincRNA coding
**Resulting `reference/gencode.v26.genes.v8.patched_contigs_genetypes_autosomal.txt` may be empty for some reason**
```{bash gencode_genes, eval=FALSE, cache=TRUE}
bash preprocessing/process.reference.files_v8.sh ${RAREVARDIR}/data/GTEx/gencode.v26.GRCh38.genes.gtf.gz ${RAREVARDIR}/download/gencode/gencode.v26.annotation.gtf
```


### Outlier calling
Using outlier calling method from watershed paper

https://github.com/nmferraro5/correlation_outliers

### Finding rare variants

#### Get allele frequency from 1KG for AFR

#### Get allele frequency from gnomAD for AFR
to-do

#### Confirm rarity of GTEx variants in 1KG and gnomAD
python scripts to do this


## Training Watershed models
look at RIVER folder in RIVER repo for RIVER.Rmd

## Analysis and Figures
look at RIVER folder in RIVER repo for enrichment_AFR.Rmd


