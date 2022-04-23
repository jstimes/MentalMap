"""Generates normalized versions of summary statistics for traits from GWAS Catalog."""

from os import listdir
from typing import List, Set

import pandas as pd

# Input file directory.
INPUT_DIR = "../data/gwas/raw/"
OUTPUT_DIR = "../data/gwas/clean/"

# Files with this suffix pattern are assumed to be trait data files,
# assumed these file name prefixes are the main/parent trait name.
# E.g. `schizophrenia_gwas_catalog_2022.csv` -> schizophrenia
TRAIT_FILE_SUFFIX = "_gwas_catalog_2022.csv"
METADATA_FILE_PATH = "../data/gwas/gwas_trait_metadata.csv"

UNKNOWN_GENE = "UNKNOWN"
CHILD_TRAIT_DELIMITER = ";"
# Some variants report this value in 'Variant and risk allele', and have no location information.
BAD_VARIANT_VALUE = "undefined"

# The 'Location ' column is mostly fine. Some associations don't report an rs ID in 'Variant
# and risk allele' and instead report in the format 'chrX:location-<b>?</b>',
# and then their 'Location' value is 'Mapping not available'. Just copy over the
# location info in same format as other variants. Note some rs ID variants also
# use 'Mapping not available' though, so we leave those as is.
NON_RSVAR_FORMAT_LOCATION_VALUE = "Mapping not available"


def main():
    trait_files = get_input_trait_files_()
    for input_file in trait_files:
        process_file_(INPUT_DIR + input_file)

    print("Done")


def get_input_trait_files_() -> List[str]:
    files = listdir(INPUT_DIR)
    trait_files = [f for f in files if TRAIT_FILE_SUFFIX in f]
    return trait_files


def process_file_(input_file_path: str) -> None:
    main_trait = get_trait_from_file_path_(input_file_path)
    traits_to_retain = get_child_traits_(main_trait)
    traits_to_retain.add(main_trait)
    print(f"Cleaning file {input_file_path} for trait {main_trait}.")

    clean_df = clean_file_(input_file_path, traits_to_retain)
    write_to_file_(main_trait, clean_df)


def get_child_traits_(main_trait: str) -> Set[str]:
    metadata_df = pd.read_csv(METADATA_FILE_PATH)
    trait_row = metadata_df.loc[metadata_df["Trait"] == main_trait]
    traits = (
        trait_row["Child traits"].astype(str).tolist()[0].split(CHILD_TRAIT_DELIMITER)
    )
    return set(map(normalize_trait_name_, traits))


def get_trait_from_file_path_(file_path: str) -> str:
    file_name = file_path[len(INPUT_DIR) :]
    trait = file_name.split(TRAIT_FILE_SUFFIX)[0].replace("_", " ")
    return trait


def clean_file_(input_file_path: str, traits_to_retain: Set[str]) -> pd.DataFrame:
    df = pd.read_csv(input_file_path)

    # Parse P-values as numbers.
    df["P-value_norm"] = df["P-value"].apply(pval_to_num_)

    # Filter to only records with the canonical trait.
    # It's important this step occurs before the gene normalization
    # as that step intentionally creates duplicate variant entries.
    df = filter_by_traits_(df, traits_to_retain)

    # Sanity-check that all duplicated variants have same mapped-gene value.
    sanity_check_duplicate_variants_(df)

    # Remove rows that don't really have the minimum necessary data.
    df = df[df["Variant and risk allele"] != BAD_VARIANT_VALUE]

    # Remove duplicate variants by taking only the one with lowest p-value.
    df = df.sort_values("P-value_norm").drop_duplicates(
        "Variant and risk allele", keep="first"
    )

    # Normalize 'Mapped gene' column by creating an extra row for each individual gene
    # (some have multiple).
    df = normalize_mapped_genes_(df)

    # Populate location field for variants where 'Variant and risk allele' actually reports location.
    df["Location"] = df.apply(try_fix_variant_location_, axis=1)

    return df


def pval_to_num_(pval: str) -> float:
    """P-values are reported as strings; convert to numbers for easier processing."""
    parts = pval.split(" x 10-")
    return float(parts[0]) * pow(10, -float(parts[1]))


def normalize_trait_name_(trait: str) -> str:
    """Ensure consistency in formatting/spelling."""
    return trait.strip().lower()


def filter_by_traits_(df: pd.DataFrame, traits_to_retain: Set[str]) -> pd.DataFrame:
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


def sanity_check_duplicate_variants_(df: pd.DataFrame) -> None:
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


def normalize_mapped_genes_(df: pd.DataFrame) -> pd.DataFrame:
    df["gene_norm"] = df["Mapped gene"].apply(lambda val: val.split(", "))
    df = df.explode("gene_norm")
    df["gene_norm"] = df["gene_norm"].apply(replace_unknown_gene_)
    return df


def replace_unknown_gene_(gene: str) -> str:
    """I'm pretty sure '- indicates unknown gene association."""
    return UNKNOWN_GENE if gene == "'-" else gene


def try_fix_variant_location_(variant_row: pd.Series) -> str:
    """Tries setting location value if not specified but can be derived from variant column.

    For example, given a variant of form 'chr6:55564517-<b>?</b>' returns '6:55564517'.
    Some are also in form 'chr7_140700006_I-<b>?</b>'.
    """
    if variant_row["Location"] != NON_RSVAR_FORMAT_LOCATION_VALUE:
        return variant_row["Location"]

    variant = variant_row["Variant and risk allele"]

    if "rs" in variant or "chr" not in variant:
        return NON_RSVAR_FORMAT_LOCATION_VALUE

    if "chr" in variant and ":" in variant:
        parts = variant.split(":")
        chr_num = parts[0][3]
        location = parts[1].split("-")[0]
        return f"{chr_num}:{location}"

    parts = variant.split("_")
    chr_num = parts[0][3]
    location = parts[1]
    return f"{chr_num}:{location}"


def write_to_file_(trait: str, df: pd.DataFrame) -> None:
    # First keep only the relevant, normalized columns for brevity.
    out_df = df[
        [
            "Variant and risk allele",
            "P-value_norm",
            "Trait(s)",
            "gene_norm",
            "Location",
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

    output_file = f'{OUTPUT_DIR}{trait.replace(" ", "_")}.csv'
    out_df.to_csv(output_file, index=False)
    print(f"Wrote {output_file}.")


main()
