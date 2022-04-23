"""Loads cleaned GWAS data and combines it with extra data sources (allele frequency, tissue associations."""

import re
from os import listdir
from typing import List, Set, Tuple

import pandas as pd

from dbsnp_api import get_afs_for_refsnps, print_progress_bar

CLEAN_GWAS_DIR = "../data/gwas/clean/"
AF_DIR = "../data/af/"
GTEX_DIR = "../data/gtex/clean/"
OUTPUT_DIR = "../data/joined/"

# From https://www.ebi.ac.uk/gwas/docs/methods/curation
# "? for unknown risk allele"
UNKNOWN_ALLELE = "?"
UNKNOWN_AF = -1.0
ALLELE_REGEX_PATTERN = r"<b>(.*)</b>"

# Separates distinct tissues in the to-be-added 'tissues' column.
TISSUE_DELIM = "&"

# Caches allele frequencies (AF) for alleles for each given SNP ID.
afs_ = {}
# Caches mappings from variants to sets of significantly associated tissues.
variant_to_tissues_ = {}


def main() -> None:
    gwas_files = get_input_trait_files_()
    for gwas_file in gwas_files:
        process_file_(CLEAN_GWAS_DIR + gwas_file)

    print("Done")


def get_input_trait_files_() -> List[str]:
    files = listdir(CLEAN_GWAS_DIR)
    return [file for file in files if ".csv" in file]


def get_tissue_variant_files_() -> List[str]:
    files = listdir(GTEX_DIR)
    return [file for file in files if ".txt" in file]


def process_file_(file_path: str) -> None:
    trait = get_trait_from_file_path_(file_path)
    print(f"Processing {trait} from {file_path}")
    df = pd.read_csv(file_path)

    df = append_af_data_(df, trait)
    df = append_tissue_data_(df)

    output_file = f'{OUTPUT_DIR}{trait.replace(" ", "_")}.csv'
    df.to_csv(output_file, index=False)
    print(f"Wrote {output_file}")


def get_trait_from_file_path_(file_path: str) -> str:
    file_name = file_path[len(CLEAN_GWAS_DIR) :]
    trait = file_name.split(".")[0].replace("_", " ")
    return trait


def get_tissue_from_file_(file: str) -> str:
    return file.split(".")[0]


def append_af_data_(input_df: pd.DataFrame, trait: str) -> pd.DataFrame:
    """Looks up AF info for all variants and appends column 'af' with this data to a new DataFrame."""
    out_df = input_df.copy()
    trait_af_file = f'{trait.replace(" ", "_")}.csv'
    if trait_af_file in listdir(AF_DIR):
        print(f"AF data cached, loading from file.")
        af_df = pd.read_csv(f"{AF_DIR}{trait_af_file}")

        out_df = out_df.merge(af_df, on="variant_and_allele")
        return out_df

    print(f"AF data not cached, loading from dbSNP.")
    all_variants = out_df["variant_and_allele"].tolist()
    rs_variants_and_alleles = [var for var in all_variants if "rs" in var]
    rs_variants = [parse_variant_(var) for var in rs_variants_and_alleles]
    global afs_
    afs_ = {}

    afs_ = get_afs_for_refsnps(rs_variants)
    out_df["af"] = out_df["variant_and_allele"].map(try_get_af_for_variant_and_allele_)

    # Cache dbSNP data to avoid lengthy API calls on re-runs.
    maf_output_file = f'{AF_DIR}{trait.replace(" ", "_")}.csv'
    out_df[["variant_and_allele", "af"]].to_csv(maf_output_file)
    print(f"Wrote {maf_output_file}.")

    return out_df


def try_get_af_for_variant_and_allele_(variant_and_allele: str) -> float:
    """Finds AF for input if it is known, otherwise returns UNKNOWN_AF."""
    variant, allele = parse_variant_and_allele_(variant_and_allele)
    if variant not in afs_:
        return UNKNOWN_AF

    var_afs = afs_[variant]
    if allele == UNKNOWN_ALLELE or allele not in var_afs:
        return UNKNOWN_AF

    return var_afs[allele]


def parse_variant_and_allele_(variant_and_allele: str) -> Tuple[str, str]:
    """Given a string like 'rs1001780-<b>G</b>', returns ('rs1001780', 'G')."""
    parts = variant_and_allele.split("-")
    if len(parts) < 2:
        return variant_and_allele, UNKNOWN_ALLELE
    variant = parts[0]
    allele_matches = re.findall(ALLELE_REGEX_PATTERN, parts[1], flags=0)
    if len(allele_matches) == 0:
        return variant, UNKNOWN_ALLELE
    return variant, allele_matches[0]


def parse_variant_(variant_and_allele: str) -> str:
    """Given a string like 'rs1001780-<b>G</b>', returns 'rs1001780'."""
    return parse_variant_and_allele_(variant_and_allele)[0]


def append_tissue_data_(df: pd.DataFrame) -> pd.DataFrame:
    """Determines tissue associations of all variants in the data frame.

    Opens processed tissue files which contain variants significantly associated with gene expression in those tissues,
    finds matching variants in the data frame, and adds the matches to a new data frame.
    """
    variant_positions = set(df["location"])
    tissue_files = get_tissue_variant_files_()
    global variant_to_tissues_
    variant_to_tissues_ = {}
    print(f"Finding tissue associations...")
    for idx, tissue_file in enumerate(tissue_files):
        print_progress_bar(idx, len(tissue_files))
        tissue = get_tissue_from_file_(tissue_file)
        significant_tissue_variants = find_tissue_associations_(
            variant_positions, f"{GTEX_DIR}{tissue_file}"
        )
        for variant in significant_tissue_variants:
            if variant not in variant_to_tissues_:
                variant_to_tissues_[variant] = set()
            variant_to_tissues_[variant].add(tissue)

    print_progress_bar(len(tissue_files), len(tissue_files))
    print("Writing associations to DF...")
    out_df = df.copy()
    out_df["tissues"] = out_df["location"].map(get_associated_tissues_)
    return out_df


def gtex_location_to_gwas_location_(gtex_location: str) -> str:
    """Converts variant locations in GWAS catalog format to GTEx format.

    i.e. given 'chr1_64764_C_T_b38', returns '4:79296443'.
    """
    parts = gtex_location.split("_")
    chr = parts[0][3]
    return f"{chr}:{parts[1]}"


def find_tissue_associations_(variant_pos_set: Set, tissue_filepath: str) -> Set[str]:
    """Identifies set of variants with significant tissue associations."""
    tissue_associations = set()
    with open(tissue_filepath) as tissue_file:
        for variant_pos_with_allele in tissue_file.readlines():
            gwas_variant_pos = gtex_location_to_gwas_location_(variant_pos_with_allele)
            if gwas_variant_pos in variant_pos_set:
                tissue_associations.add(gwas_variant_pos)
    return tissue_associations


def get_associated_tissues_(variant_location: str) -> str:
    """Returns concatenated string of all tissues which this variant is associated with."""
    if variant_location not in variant_to_tissues_:
        return ""

    tissues_set = variant_to_tissues_[variant_location]
    return TISSUE_DELIM.join(tissues_set)


main()
