"""Given a CSV file containing summary statistics for a trait on GWAS Catalog, generates a normalized version of the data."""

import re
from os import listdir

import pandas as pd

from dbsnp_api import get_mafs_for_refsnps

# Input file directory.
INPUT_DIR = "../data/raw/"
OUTPUT_DIR = "../data/clean/"

# Files with this suffix pattern are assumed to be trait data files,
# assumed these file name prefixes are the main/parent trait name.
# E.g. `schizophrenia_gwas_catalog_2022.csv` -> schizophrenia
TRAIT_FILE_SUFFIX = "_gwas_catalog_2022.csv"
METADATA_FILE = "gwas_trait_metadata.csv"

# From https://www.ebi.ac.uk/gwas/docs/methods/curation
# "? for unknown risk allele"
UNKNOWN_ALLELE = "?"
UNKNOWN_GENE = "UNKNOWN"
UNKNOWN_AF = -1.0
CHILD_TRAIT_DELIMITER = ";"
ALLELE_REGEX_PATTERN = r"<b>(.*)</b>"

# Caches mean allele frequencies (MAF) for alleles for each given SNP ID.
mafs_ = {}


def main():
    trait_files = get_input_trait_files_()
    for input_file in trait_files:
        process_file_(INPUT_DIR + input_file)

    print("Done")


def get_input_trait_files_():
    files = listdir(INPUT_DIR)
    trait_files = [f for f in files if TRAIT_FILE_SUFFIX in f]
    return trait_files


def process_file_(input_file_path):
    main_trait = get_trait_from_file_path_(input_file_path)
    if main_trait in done:
        return
    traits_to_retain = get_child_traits_(main_trait)
    traits_to_retain.append(main_trait)
    print(f"Cleaning file {input_file_path} for trait {main_trait}.")

    clean_df = clean_file_(input_file_path, traits_to_retain)
    clean_df = append_maf_data_(clean_df)
    write_to_file_(main_trait, clean_df)


def get_child_traits_(main_trait):
    metadata_df = pd.read_csv(INPUT_DIR + METADATA_FILE)
    trait_row = metadata_df.loc[metadata_df["Trait"] == main_trait]
    traits = (
        trait_row["Child traits"].astype(str).tolist()[0].split(CHILD_TRAIT_DELIMITER)
    )
    return list(map(normalize_trait_name_, traits))


def get_trait_from_file_path_(file_path):
    file_name = file_path[len(INPUT_DIR) :]
    trait = file_name.split(TRAIT_FILE_SUFFIX)[0].replace("_", " ")
    return trait


def clean_file_(input_file_path, traits_to_retain):
    df = pd.read_csv(input_file_path)

    # Parse P-values as numbers.
    df["P-value_norm"] = df["P-value"].apply(pval_to_num_)

    # Filter to only records with the canonical trait.
    # It's important this step occurs before the gene normalization
    # as that step intentionally creates duplicate variant entries.
    df = filter_by_traits_(df, traits_to_retain)

    # Sanity-check that all duplicated variants have same mapped-gene value.
    sanity_check_duplicate_variants_(df)
    # Remove duplicate variants by taking only the one with lowest p-value.
    df = df.sort_values("P-value_norm").drop_duplicates(
        "Variant and risk allele", keep="first"
    )

    # Normalize 'Mapped gene' column by creating an extra row for each individual gene
    # (some have multiple).
    df = normalize_mapped_genes_(df)

    return df


def append_maf_data_(input_df):
    """Looks up MAF info for all variants and appends column 'maf' with this data to a new DataFrame."""
    all_variants = input_df["Variant and risk allele"].tolist()
    rs_variants_and_alleles = [var for var in all_variants if "rs" in var]
    rs_variants = [parse_variant_(var) for var in rs_variants_and_alleles]
    global mafs_
    mafs_ = {}
    mafs_ = get_mafs_for_refsnps(rs_variants)
    out_df = input_df.copy()
    out_df["maf"] = out_df["Variant and risk allele"].map(
        try_get_maf_for_variant_and_allele_
    )
    return out_df


def try_get_maf_for_variant_and_allele_(variant_and_allele):
    """Finds AF for input if it is known, otherwise returns UNKNOWN_AF."""
    variant, allele = parse_variant_and_allele_(variant_and_allele)
    if variant not in mafs_:
        return UNKNOWN_AF

    var_mafs = mafs_[variant]
    if allele == UNKNOWN_ALLELE or allele not in var_mafs:
        return UNKNOWN_AF

    return var_mafs[allele]


def parse_variant_and_allele_(variant_and_allele):
    """Given a string like 'rs1001780-<b>G</b>', returns ('rs1001780', 'G')."""
    parts = variant_and_allele.split("-")
    if len(parts) < 2:
        return variant_and_allele, UNKNOWN_ALLELE
    variant = parts[0]
    allele_matches = re.findall(ALLELE_REGEX_PATTERN, parts[1], flags=0)
    if len(allele_matches) == 0:
        return variant, UNKNOWN_ALLELE
    return variant, allele_matches[0]


def parse_variant_(variant_and_allele):
    """Given a string like 'rs1001780-<b>G</b>', returns 'rs1001780'."""
    return parse_variant_and_allele_(variant_and_allele)[0]


def pval_to_num_(pval):
    """P-values are reported as strings; convert to numbers for easier processing."""
    parts = pval.split(" x 10-")
    return float(parts[0]) * pow(10, -float(parts[1]))


def normalize_trait_name_(trait):
    """Ensure consistency in formatting/spelling."""
    return trait.strip().lower()


def filter_by_traits_(df, traits_to_retain):
    """Filters the DataFrame such that only rows with trait value in the given list are kept."""
    # First normalize trait column (un-capitalize)
    df["Trait(s)"] = df["Trait(s)"].map(normalize_trait_name_)
    original_total_rows = len(df)
    df["is_retained_trait"] = df["Trait(s)"].map(lambda t: t in traits_to_retain)
    filtered_df = df[df["is_retained_trait"] == True]
    filtered_rows = len(df) - len(filtered_df)
    print(
        f"Trait filtering removed {filtered_rows} rows ({original_total_rows} originally)."
    )
    return filtered_df


def sanity_check_duplicate_variants_(df):
    duplicate_variants = (
        df.groupby("Variant and risk allele")
        .filter(lambda x: len(x) > 1)["Variant and risk allele"]
        .unique()
    )
    all_good = True
    for variant in duplicate_variants:
        all_mapped_genes = df[df["Variant and risk allele"] == variant][
            "Mapped gene"
        ].unique()
        if len(all_mapped_genes) > 1:
            print(f"Found variant, {variant}, with differing mapped gene values.")
            all_good = False

    if all_good:
        print("No repeated variants with differing mapped gene values.")


def normalize_mapped_genes_(df):
    df["gene_norm"] = df["Mapped gene"].apply(lambda val: val.split(", "))
    df = df.explode("gene_norm")
    df["gene_norm"] = df["gene_norm"].apply(replace_unknown_gene_)
    return df


def replace_unknown_gene_(gene):
    """I'm pretty sure '- indicates unknown gene association."""
    return UNKNOWN_GENE if gene == "'-" else gene


def write_to_file_(trait, df):
    # First keep only the relevant, normalized columns for brevity.
    out_df = df[
        [
            "Variant and risk allele",
            "P-value_norm",
            "Trait(s)",
            "gene_norm",
            "Location",
            "maf",
        ]
    ]
    column_remapping = {
        "Variant and risk allele": "variant_and_allele",
        "P-value_norm": "p_value",
        "Trait(s)": "trait",
        "gene_norm": "gene",
        "Location": "location",
    }
    out_df = out_df.rename(columns=column_remapping)

    output_file = f'{OUTPUT_DIR}{trait.replace(" ", "_")}_cleaned.csv'
    out_df.to_csv(output_file)
    print(f"Wrote {output_file}.")

    # Cache dbSNP data to avoid lengthy API calls on re-runs.
    maf_output_file = f'{OUTPUT_DIR}{trait.replace(" ", "_")}_variant_mafs.csv'
    out_df[["variant_and_allele", "maf"]].to_csv(maf_output_file)
    print(f"Wrote {maf_output_file}.")


main()
