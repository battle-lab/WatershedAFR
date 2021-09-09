"""Parses vcf output from VEP/LOFTEE to a table where
rows are rare variants and columns are the annotations
"""
import os
import argparse

import pandas as pd
from pysam import VariantFile

def get_info_fields(vcf):
    """Get info fields in vcf

    Parameters
    ----------
    vcf : pysam.libcbcf.VariantFile
        VCF file read in by pysam

    Returns
    -------
    list
    """
    info_fields=  vcf.header.info['CSQ'].record['Description'] \
        .split('Format: ')[-1].strip('"').split('|')
    return(info_fields)


def get_info_dictionary(rec, info_fields):
    """Parses info into a dictionary with info fields as the keys

    Parameters
    ----------
    rec : pysam.libcbcf.VariantRecord
        Obatined from iterating over pysam.VariantFile.fetch()
    info_fields : list
        List of info fields in the vcf
    
    Returns
    -------
    dict
        Info field parsed into a dictionary
    """

    info_data = rec.info['CSQ']
    info_dict = [dict(zip(info_fields, x.split('|'))) \
                 for x in info_data if len(x.split('|')) == len(info_fields)]
    return(info_dict)

def vcf_to_dataframe(vcf_filepath, write_path=None):
    """Parse info from vcf into a dataframe with genomic coordinates and 
    info fields as the columns

    Parameters
    ----------
    vcf_filepath : str
        Path to the VEP/LOFTEE output vcf
    write_path : None or str
        If not None, the function will write the returned dataframe to write_path

    Returns
    -------
    pandas.DataFrame
    """

    # Read vcf
    vcf = VariantFile(vcf_filepath)
    info_fields = get_info_fields(vcf)

    # Parse info from vcf into dataframe 
    df_list = []
    for i, rec in enumerate(vcf.fetch()):

        info_dict = get_info_dictionary(rec, info_fields)

        df = pd.DataFrame(info_dict)
        df.insert(0, 'Alt', rec.alts[0])
        df.insert(0, 'Pos', rec.pos)
        df.insert(0, 'Chrom', 'chr' + rec.chrom)
        df_list.append(df)

    info_df = pd.concat(df_list, ignore_index=True)

    # Write to file
    if write_path != None:
        info_df.to_csv(write_path, sep='\t')
        print(f'VCF info dataframe saved to\n{write_path}')

    return(info_df)


def main():
    filename = '/scratch/groups/abattle4/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.vep.loftee.vcf'

    info_df_filepath = os.path.splitext(filename)[0] + '.info.tsv'

    if not os.path.isfile(info_df_filepath):
        info_df = vcf_to_dataframe(filename, info_df_filepath)
    else:
        info_df = pd.read_csv(info_df_filepath, sep='\t', index_col=0)

    # Retrieve VEP annotations and LoF as binary values
    #coord_df = info_df[['Chrom','Pos','Alt','Allele']]
    #vep_df = pd.get_dummies(info_df['Consequence'])
    #lof_df = pd.get_dummies(info_df['LoF'])

    #anno_df = coord_df.join(vep_df).join(lof_df)

    # rename and select for annotation columns to be consistent with the Watershed paper table S3

    # collapse transcripts so we get SNV level annotations



if __name__ == '__main__':
    main()

