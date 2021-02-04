# Population Specific Finding Rare Variants



```bash
# Defining root, data, and raw data directories
rootdir=/scratch/groups/abattle4/victor/WatershedAFR
datadir=${rootdir}/data
rawdir=${rootdir}/raw_data
```



### Expression data correction and normalization


#### Generate tpm and read count matrices from GTEx V8

Generate file mapping sample identifiers to tissues. Restrict to samples that pass RNA-seq QC (marked as RNASEQ in SMAFRZE column).
```bash
cat ${rawdir}/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | tail -n+2 | cut -f1,7,17 | \
  sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | \
  awk '$3=="RNASEQ" {print $1"\t"$2}' | sort -k 1 > ${datadir}/data_prep/gtex_v8_samples_tissues.txt
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
OUT_AFR=${datadir}/data_prep/PEER_AFR
OUT_EUR=${datadir}/data_prep/PEER_EUR
GTEX=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
SAMPLE_AFR=${datadir}/data_prep/gtex_v8_samples_tissues_AFR.txt
SAMPLE_EUR=${datadir}/data_prep/gtex_v8_samples_tissues_EUR.txt
END='.tpm.txt'
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX --out $OUT_AFR --sample $SAMPLE_AFR --end $END
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX --out $OUT_EUR --sample $SAMPLE_EUR --end $END

# split read counts
GTEX=${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
END='.reads.txt'
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX --out $OUT_AFR --sample $SAMPLE_AFR --end $END
python code/preprocessing/data_prep/split_expr_by_tissues.py --gtex $GTEX --out $OUT_EUR --sample $SAMPLE_EUR --end $END

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
#African
bash code/preprocessing/data_prep/get_tissue_by_individual.sh \
  -p ${datadir}/data_prep/PEER_AFR \
  -o ${datadir}/data_prep/gtex_AFR_tissue_by_ind.txt \
  -t ${datadir}/data_prep/gtex_AFR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_AFR_individuals_all_normalized_samples.txt
  
#European
bash code/preprocessing/data_prep/get_tissue_by_individual.sh \
  -p ${datadir}/data_prep/PEER_EUR \
  -o ${datadir}/data_prep/gtex_EUR_tissue_by_ind.txt \
  -t ${datadir}/data_prep/gtex_EUR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_EUR_individuals_all_normalized_samples.txt

```

Combine the PEER-corrected data

Combined files will be located in ${datadir}/data_prep/gtex_AFR_normalized_expression.txt.gz and ${datadir}/data_prep/gtex_EUR_normalized_expression.txt.gz
```bash
#African
python code/preprocessing/data_prep/gather_filter_normalized_expression.py \
  -p ${datadir}/data_prep/PEER_AFR \
  -t ${datadir}/data_prep/gtex_AFR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_AFR_individuals_all_normalized_samples.txt \
  -o ${datadir}/data_prep/gtex_AFR_normalized_expression.txt

gzip ${datadir}/data_prep/gtex_AFR_normalized_expression.txt

#European
python code/preprocessing/data_prep/gather_filter_normalized_expression.py \
  -p ${datadir}/data_prep/PEER_EUR \
  -t ${datadir}/data_prep/gtex_EUR_tissues_all_normalized_samples.txt \
  -i ${datadir}/data_prep/gtex_EUR_individuals_all_normalized_samples.txt \
  -o ${datadir}/data_prep/gtex_EUR_normalized_expression.txt

gzip ${datadir}/data_prep/gtex_EUR_normalized_expression.txt
```

### Outlier calling

Call outliers on African individuals

Saved to `${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt`
```bash
Rscript code/preprocessing/outlier_calling/call_outliers.R \
  --Z.SCORES=${datadir}/data_prep/gtex_AFR_normalized_expression.txt.gz \
  --OUT=${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs.txt \
  --POP=${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt \
  --N.PHEN=5 --ZTHRESH=3


# Remove global outliers
Rscript code/preprocessing/outlier_calling/identify_global_outliers.R \
  --OUTLIERS=${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs.txt \
  --METHOD=proportion
```

Call outliers on European individuals

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

### Finding rare variants
Requires `bcftools` and `bedtools`

#### African Individuals
Rare variants file will be saved to `{datadir}/rare_variants/gene-AFR-rv.txt`

```bash
bash code/preprocessing/rare_variants/find_rare_variants_changed_outliers_only.sh \
-d ${datadir}/rare_variants_pop_norm \
-g ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
-r ${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-l ${datadir}/data_prep/gtex_v8_wgs_individuals_AFR.txt \
-p "AFR" \
-o ${datadir}/outlier_calling/AFR/gtexV8.AFR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt

```

#### European Individuals
Rare variants file will be saved to `{datadir}/rare_variants/gene-EUR-rv.txt`

```bash
bash code/preprocessing/rare_variants/find_rare_variants_changed_outliers_only.sh \
-d ${datadir}/rare_variants_pop_norm \
-g ${rawdir}/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
-r ${datadir}/data_prep/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed \
-l ${datadir}/data_prep/gtex_v8_wgs_individuals_EUR.txt \
-p "EUR" \
-o ${datadir}/outlier_calling/EUR/gtexV8.EUR.outlier.controls.v8ciseQTLs_globalOutliersRemoved.txt

```