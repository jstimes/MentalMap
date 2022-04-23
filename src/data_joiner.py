"""Loads cleaned GWAS data and combines it with extra data sources (allele frequency, tissue associations."""

import re
from os import listdir
from typing import Any, Dict, List, Mapping, Sequence, Set, Text, Tuple

import pandas as pd

from dbsnp_api import get_mafs_for_refsnps

CLEAN_GWAS_DIR = "../data/gwas/clean/"
AF_DIR = "../data/af/"
OUTPUT_DIR = "../data/joined/"

# From https://www.ebi.ac.uk/gwas/docs/methods/curation
# "? for unknown risk allele"
UNKNOWN_ALLELE = "?"
UNKNOWN_AF = -1.0
ALLELE_REGEX_PATTERN = r"<b>(.*)</b>"

# Caches allele frequencies (AF) for alleles for each given SNP ID.
afs_ = {}


def main():
    gwas_files = get_input_trait_files_()
    for gwas_file in gwas_files:
        process_file_(CLEAN_GWAS_DIR + gwas_file)

    print("Done")


def get_input_trait_files_() -> List[Text]:
    files = listdir(CLEAN_GWAS_DIR)
    return [file for file in files if ".csv" in file]


def process_file_(file_path: Text):
    trait = get_trait_from_file_path_(file_path)
    print(f"Processing {trait} from {file_path}")
    df = pd.read_csv(file_path)
    df = append_maf_data_(df, trait)
    output_file = f'{OUTPUT_DIR}{trait.replace(" ", "_")}.csv'
    df.to_csv(output_file, index=False)
    print(f"Wrote {output_file}")


def get_trait_from_file_path_(file_path: Text) -> Text:
    file_name = file_path[len(CLEAN_GWAS_DIR) :]
    trait = file_name.split(".")[0].replace("_", " ")
    return trait


def append_maf_data_(input_df: pd.DataFrame, trait: Text) -> pd.DataFrame:
    """Looks up MAF info for all variants and appends column 'af' with this data to a new DataFrame."""
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

    afs_ = get_mafs_for_refsnps(rs_variants)
    out_df["af"] = out_df["variant_and_allele"].map(try_get_maf_for_variant_and_allele_)

    # Cache dbSNP data to avoid lengthy API calls on re-runs.
    maf_output_file = f'{AF_DIR}{trait.replace(" ", "_")}.csv'
    out_df[["variant_and_allele", "af"]].to_csv(maf_output_file)
    print(f"Wrote {maf_output_file}.")

    return out_df


def try_get_maf_for_variant_and_allele_(variant_and_allele: Text) -> float:
    """Finds AF for input if it is known, otherwise returns UNKNOWN_AF."""
    variant, allele = parse_variant_and_allele_(variant_and_allele)
    if variant not in afs_:
        return UNKNOWN_AF

    var_afs = afs_[variant]
    if allele == UNKNOWN_ALLELE or allele not in var_afs:
        return UNKNOWN_AF

    return var_afs[allele]


def parse_variant_and_allele_(variant_and_allele: Text) -> Tuple[Text, Text]:
    """Given a string like 'rs1001780-<b>G</b>', returns ('rs1001780', 'G')."""
    parts = variant_and_allele.split("-")
    if len(parts) < 2:
        return variant_and_allele, UNKNOWN_ALLELE
    variant = parts[0]
    allele_matches = re.findall(ALLELE_REGEX_PATTERN, parts[1], flags=0)
    if len(allele_matches) == 0:
        return variant, UNKNOWN_ALLELE
    return variant, allele_matches[0]


def parse_variant_(variant_and_allele: Text) -> Text:
    """Given a string like 'rs1001780-<b>G</b>', returns 'rs1001780'."""
    return parse_variant_and_allele_(variant_and_allele)[0]


main()
