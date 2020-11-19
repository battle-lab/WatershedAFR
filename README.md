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

```bash
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

<!---
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


For each tissue, filter for genes with > 20% individuals with TPM > 0.1 and read count > 6, Log2(tpm + 2) transfrom the data, and then z-transform.
```
Rscript code/preprocessing/data_prep/preprocess_expr.R
```

Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles
```


bash code/preprocessing/data_prep/get_genotypes.sh
```
--->


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

<!---
#### Make flat files from normalized data

```{bash, eval=FALSE, cache=TRUE}
python preprocessing/data_prep/gather_filter_normalized_expression_v8.py
```
--->

#### Prepare gene annotations from gencode

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


### Outlier calling
<!---
Using outlier calling method from watershed paper

https://github.com/nmferraro5/correlation_outliers
--->

GTEx v8 outliers are located in 
`/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/outlier_calls/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt`
Save to `${datadir}/outlier_calling`

### Finding rare variants
Requires `bcftools` and `bedtools`

#### European Individuals
Rare variants file will be saved to `{datadir}/rare_variants/gene-EUR-rv.txt`

```bash
bash code/preprocessing/rare_variants/find_rare_variants.sh \
-d ${datadir}/rare_variants \
-g ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
-r ${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-l ${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt \
-p "EUR" \
-o ${datadir}/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt

```

#### African Individuals
Rare variants file will be saved to `{datadir}/rare_variants/gene-AFR-rv.txt`

```bash
bash code/preprocessing/rare_variants/find_rare_variants.sh \
-d ${datadir}/rare_variants \
-g ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
-r ${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-l ${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt \
-p "AFR" \
-o ${datadir}/outlier_calling/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt

```


## Training Watershed models
look at RIVER folder in RIVER repo for RIVER.Rmd

## Analysis and Figures
look at RIVER folder in RIVER repo for enrichment_AFR.Rmd


