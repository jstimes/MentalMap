"""Given a CSV file containing summary statistics for a trait on GWAS Catalog, generates a normalized version of the data."""

import pandas as pd
import sys


# Command line arguments.
# TODO: clean this up a little.
# First argument is expected to be input file
INPUT_FILE = sys.argv[1]
# Second argument is the canonical trait
TRAIT = sys.argv[2]

UNKNOWN_GENE = "UNKNOWN"
OUTPUT_DIR = "clean/"


def main():
    if INPUT_FILE is None or TRAIT is None:
        print("Need to specify input file and trait.")
        return

    print(f"Cleaning file {INPUT_FILE} for trait {TRAIT}.")

    clean_df = clean_file_()
    write_to_file_(clean_df)


def clean_file_():
    df = pd.read_csv(INPUT_FILE)

    # Parse P-values as numbers.
    df['P-value_norm'] = df['P-value'].apply(pval_to_num_)

    # Filter to only records with the canonical trait.
    # It's important this step occurs before the gene normalization
    # as that step intentionally creates duplicate variant entries.
    df = filter_by_trait_(df)

    # Sanity-check that all duplicated variants have same mapped-gene value.
    sanity_check_duplicate_variants_(df)
    # Remove duplicate variants by taking only the one with lowest p-value.
    df = df.sort_values('P-value_norm').drop_duplicates(
        'Variant and risk allele', keep='first')

    # Normalize 'Mapped gene' column by creating an extra row for each individual gene
    # (some have multiple).
    df = normalize_mapped_genes_(df)

    return df


def pval_to_num_(pval):
  parts = pval.split(" x 10-")
  return float(parts[0]) * pow(10, -float(parts[1]))


def filter_by_trait_(df):
    original_total_rows = len(df)
    filtered_df = df[df['Trait(s)'] == TRAIT]
    filtered_rows = len(df) - len(filtered_df)
    print(f"Trait filtering removed {filtered_rows} rows ({original_total_rows} originally).")
    return filtered_df


def sanity_check_duplicate_variants_(df):
    duplicate_variants = df.groupby('Variant and risk allele').filter(lambda x: len(x) > 1)[
        'Variant and risk allele'].unique()
    all_good = True
    for variant in duplicate_variants:
        all_mapped_genes = df[df['Variant and risk allele'] == variant]['Mapped gene'].unique()
        if len(all_mapped_genes) > 1:
            print(f"Found variant, {variant}, with differing mapped gene values.")
            all_good = False

    if all_good:
        print("No repeated variants with differing mapped gene values.")


def normalize_mapped_genes_(df):
    df['gene_norm'] = df['Mapped gene'].apply(lambda val: val.split(", "))
    df = df.explode('gene_norm')
    df['gene_norm'] = df['gene_norm'].apply(replace_unknown_gene_)
    return df


def replace_unknown_gene_(gene):
    """I'm pretty sure '- indicates unknown gene association."""
    # TODO: should double-check this.
    return UNKNOWN_GENE if gene == "'-" else gene


def write_to_file_(df):
    # First keep only the relevant, normalized columns for brevity.
    out_df = df[['Variant and risk allele', 'P-value_norm', 'Trait(s)', 'gene_norm']]
    column_remapping = {
        'Variant and risk allele': 'variant_and_allele',
        'P-value_norm': 'p_value',
        'Trait(s)': 'trait',
        'gene_norm': 'gene',
    }
    out_df = out_df.rename(columns=column_remapping)

    output_file = INPUT_FILE.split(".")[0] + "_cleaned.csv"
    output_file = output_file.replace("raw/", OUTPUT_DIR)
    out_df.to_csv(output_file)
    print(f"Wrote {output_file}.")


main()
