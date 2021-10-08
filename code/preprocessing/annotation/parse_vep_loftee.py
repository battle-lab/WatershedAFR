"""Parses vcf output from VEP/LOFTEE to a table where
rows are rare variants and columns are the annotations
"""
import os
import sys
import argparse

import numpy as np
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
    info_fields = (
        vcf.header.info["CSQ"]
        .record["Description"]
        .split("Format: ")[-1]
        .strip('"')
        .split("|")
    )

    return info_fields


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

    info_data = rec.info["CSQ"]
    info_dict = [
        dict(zip(info_fields, x.split("|")))
        for x in info_data
        if len(x.split("|")) == len(info_fields)
    ]

    return info_dict


def vcf_to_dataframe(vcf_filepath):
    """Parse info from vcf into dataframes with genomic coordinates and
    info fields as the columns. Makes a dataframe for each chromosome
    and returns list of paths to the saved dataframe.

    Parameters
    ----------
    vcf_filepath : str
        Path to the VEP/LOFTEE output vcf

    Returns
    -------
    list
    """

    # Read vcf
    vcf = VariantFile(vcf_filepath)
    info_fields = get_info_fields(vcf)

    # Parse info from vcf into dataframe
    df_list = []
    cur_chrom = None
    write_path_list = []

    for i, rec in enumerate(vcf.fetch()):

        # Keep track of current chrom
        chrom = "chr" + rec.chrom

        if cur_chrom is None:  # Starting chromosome
            cur_chrom = chrom

        if chrom != "chr22":
            continue

        cur_chrom = chrom

        if cur_chrom != chrom:  # Save and print the chromosome dataframe
            info_df = pd.concat(df_list, ignore_index=True)

            write_path = (
                vcf_filepath.split(".vcf.gz")[0] + "." + cur_chrom + ".info.tsv"
            )
            info_df.to_csv(write_path, sep="\t")
            print(f"VCF info dataframe saved to\n{write_path}")

            # Reset list of dataframes and update current chromosome
            del info_df  # Delete to save on memory
            df_list = []
            cur_chrom = chrom

        # Parse info from vcf
        info_dict = get_info_dictionary(rec, info_fields)

        df = pd.DataFrame(info_dict)
        df.insert(0, "Alt", rec.alts[0])
        df.insert(0, "Pos", rec.pos)
        df.insert(0, "Chrom", chrom)
        df_list.append(df)

    # Save and print the final chromosome
    info_df = pd.concat(df_list, ignore_index=True)
    write_path = vcf_filepath.split(".vcf.gz")[0] + "." + cur_chrom + ".info.tsv"
    info_df.to_csv(write_path, sep="\t")
    print(f"VCF info dataframe saved to\n{write_path}")


def get_annotations_from_df(info_df):
    """Extracts annotations from the `Consequence` and `LoF` columns

    Parameters
    ----------
    info_df : pandas.Dataframe
        Dataframe created by vcf_to_dataframe

    Returns
    -------
    pandas.Dataframe
    """
    # Genomic coordinates
    coord_df = info_df[["Chrom", "Pos", "Alt", "Allele"]]
    # VEP annotations
    vep_df = pd.get_dummies(info_df["Consequence"])
    # Obtain splice_region_variant from combined consequences
    splice_col_list = [col for col in vep_df.columns if "splice_region_variant" in col]
    splice_region_variant = np.amax(vep_df[splice_col_list].values, axis=1)
    vep_df["splice_region_variant"] = splice_region_variant
    # LOFTEE annotations
    lof_df = pd.get_dummies(info_df["LoF"])

    return coord_df.join(vep_df).join(lof_df)


def test_get_annotations_from_df(info_df, anno_df):
    """Check that output of `get_annoations_from_df` is correctly representing the
    annotations from info_df (Excluding `splice_region_variant` since it is
    combined with other consequences)
    """

    # Get list of annoations from info_df
    vep_list = list(info_df["Consequence"].unique())
    lof_list = [x for x in info_df["LoF"].unique() if isinstance(x, str)]
    anno_list = vep_list + lof_list

    test_result = {}

    for anno in anno_list:
        anno_df_num = len(anno_df[anno_df[anno] == 1])

        if anno in lof_list:
            info_df_num = len(info_df[info_df["LoF"] == anno])
        else:
            info_df_num = len(info_df[info_df["Consequence"] == anno])

        test_result[anno] = anno_df_num == info_df_num

    for anno in anno_list:
        if not test_result[anno]:
            print(test_result)

            return test_result[anno]

    return True


def dataframe_to_tidy(info_df):
    """Extracts annotations from the `Consqeuence` and `LoF` columns.
    Returns a tidy dataframe with the columns
    Chrom | Pos | Alt | Allele | annotations...

    Parameters
    ----------
    info_df : pandas.Dataframe
        Dataframe created by vcf_to_dataframe

    Returns
    -------
    pandas.Dataframe
    """
    anno_list = [
        "Chrom",
        "Pos",
        "Alt",
        "Allele",
        "3_prime_UTR_variant",
        "5_prime_UTR_variant",
        "TF_binding_site_variant",
        "downstream_gene_variant",
        "intergenic_variant",
        "intron_variant",
        "missense_variant",
        "non_coding_transcript_exon_variant",
        "regulatory_region_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant",
        "stop_gained",
        "synonymous_variant",
        "upstream_gene_variant",
        "LoF_HC",
        "LoF_LC",
    ]
    anno_df = get_annotations_from_df(info_df)

    # check that annotations were properly extracted
    if test_get_annotations_from_df(info_df, anno_df):
        # Rename and select for annotation columns to be consistent with the
        # Watershed paper table S3
        anno_df.rename(columns={"HC": "LoF_HC", "LC": "LoF_LC"}, inplace=True)
        anno_df = anno_df[anno_list]

        # Collapse transcripts by taking maximum so we get SNV level
        # annotations
        anno_snv_df = anno_df.groupby(
            ["Chrom", "Pos", "Alt", "Allele"], as_index=False
        ).max()

        # Set LoF_LC to 0 if LoF_HC is 1
        lof_hc = anno_snv_df["LoF_HC"].values
        lof_lc = anno_snv_df["LoF_LC"].values
        anno_snv_df["LoF_LC"] = [0 if hc == 1 else lc for hc, lc in zip(lof_hc, lof_lc)]

        return anno_snv_df

    return None


def main():
    #    filename = "/scratch/groups/abattle4/victor/WatershedAFR/data/annotation/gene-AFR-rv.vep.loftee.vcf.gz"

    usage = "Parses VCF output from VEP/LOFTEE to a table where rows are rare variants and columns are the annotations"
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument(
        "--anno",
        required=True,
        help="Annotations from VEP and LOFTEE as vcf.gz format. Needs to be tbi indexed.",
    )
    args = parser.parse_args()

    filename = args.anno
    if not filename.endswith(".vcf.gz"):
        print("Annotation file needs to be bgzipped and tabix indexed")
        sys.exit(1)

    # Info field from VCF parsed to dataframe
    df_path_list = [
        filename.split(".vcf.gz")[0] + ".chr" + str(x) + ".info.tsv"
        for x in range(1, 23)
    ]
    if not os.path.isfile(df_path_list[0]):
        vcf_to_dataframe(filename)

    # Combine tidy dataframes for each chromosome into one file
    anno_snv_all_df = None
    for df_path in df_path_list:
        # Open info field from vcf parsed to dataframe
        info_df = pd.read_csv(df_path, sep="\t", index_col=0)

        # Retreive VEP and LoF annotations in tidy format
        anno_snv_df = dataframe_to_tidy(info_df)

        # Concatenate
        if anno_snv_all_df is None:
            anno_snv_all_df = anno_snv_df
        else:
            anno_snv_all_df = pd.concat([anno_snv_all_df, anno_snv_df])

    write_path = filename.split(".vcf.gz")[0] + ".snv.tsv"
    anno_snv_all_df.to_csv(write_path, sep="\t")
    print(f"SNV level VEP and LOFTEE annotations are saved to\n{write_path}")


if __name__ == "__main__":
    main()
