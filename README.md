# WatershedAFR
Watershed implementation for African American samples from GTEx


## Install the dependencies
* have YAML for python dependencies


## Directory structure
```
├── data
│   ├── data_prep
│   ├── figures
│   ├── models
│   ├── outlier_calling
│   └── rare_variants
├── raw_data
└── WatershedAFR
├── code
│   ├── analysis
│   ├── figures
│   ├── preprocessing
│   │   ├── data_prep
│   │   ├── outlier_calling
│   │   └── rare_variants
│   └── watershed
├── data_sources.md
├── LICENSE
├── pipelines
└── README.md
```

```
├── _config.yml
├── _drafts
│   ├── begin-with-the-crazy-ideas.textile
│   └── on-simplicity-in-technology.markdown
├── _includes
│   ├── footer.html
│   └── header.html
├── _layouts
│   ├── default.html
│   └── post.html
├── _posts
│   ├── 2007-10-29-why-every-programmer-should-play-nethack.textile
│   └── 2009-04-26-barcamp-boston-4-roundup.textile
├── _data
│   └── members.yml
├── _site
└── index.html
```

datadir=~/data/
mkdir ${datadir}
mkdir ${datadir}/data_prep/
...
mkdir ${datadir}figures/
Create raw data directory to contain (either files or links to) GTEx, 1k Genomes, gnomAD
Sub repo of Watershed (for now until we make it into an R package)



pipelines/ 
	AA_full.Rmd State of the art references ../code/data_prep.R
	Experiments in thresholding
	Enrichment.Rmd → update to call the new function
	

code/
	preprocessing/  (remember to create output sub directories as needed)
		data_prep
		outlier_calling
		rare_variants
	analysis/ (model stuff)
		watershed/ → symbolic link to our branch of watershed repo
		Function to compute enrichment of rare variants
	Watershed/ -> subrepo of the Watershed repo
		
	figures/
README detailing data of download and source
raw_data/ ←- all data sources needed to run everything in the repo; can be softlinks for large files e.g. gnomAD
	Expression TPM from gtex portal
	Genotype vcf for gtex
?PEER?
1KG
	gnomAD→ softlink only if we can’t use API and have to download the vcfs 
GENCODE gene info
	
data/ 	
data_prep/←  ${datadir}/data_prep
outlier_calling/
	rare_variants/
	models/ ← save trained models
figures/ → reminder to self to output the table and the figure from the figure scripts
	E.g. frequency table used for boxplot
