"""Given a CSV file containing summary statistics for a trait on GWAS Catalog, generates a normalized version of the data."""

import pandas as pd
from os import listdir


# Input file directory.
INPUT_DIR = '../data/raw/'
OUTPUT_DIR = '../data/clean/'

# Files with this suffix pattern are assumed to be trait data files,
# assumed these file name prefixes are the main/parent trait name.
# E.g. `schizophrenia_gwas_catalog_2022.csv` -> schizophrenia
TRAIT_FILE_SUFFIX = '_gwas_catalog_2022.csv'

UNKNOWN_GENE = 'UNKNOWN'


def main():
    trait_files = get_input_trait_files_()
    for input_file in trait_files:
        process_file_(INPUT_DIR + input_file)

    print('Done')


def get_input_trait_files_():
    files = listdir(INPUT_DIR)
    trait_files = [f for f in files if TRAIT_FILE_SUFFIX in f]
    return trait_files


def process_file_(input_file_path):
    trait = get_trait_from_file_path_(input_file_path)
    print(f'Cleaning file {input_file_path} for trait {trait}.')

    clean_df = clean_file_(input_file_path, trait)
    write_to_file_(trait, clean_df)


def get_trait_from_file_path_(file_path):
    file_name = file_path[len(INPUT_DIR):]
    trait = file_name.split(TRAIT_FILE_SUFFIX)[0].replace("_", " ")
    return trait


def clean_file_(input_file_path, trait):
    df = pd.read_csv(input_file_path)

    # Parse P-values as numbers.
    df['P-value_norm'] = df['P-value'].apply(pval_to_num_)

    # Filter to only records with the canonical trait.
    # It's important this step occurs before the gene normalization
    # as that step intentionally creates duplicate variant entries.
    df = filter_by_trait_(df, trait)

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
  parts = pval.split(' x 10-')
  return float(parts[0]) * pow(10, -float(parts[1]))


def filter_by_trait_(df, trait):
    # First normalize trait column (un-capitalize)
    df['Trait(s)'] = df['Trait(s)'].map(lambda t: t.lower())
    original_total_rows = len(df)
    filtered_df = df[df['Trait(s)'] == trait]
    filtered_rows = len(df) - len(filtered_df)
    print(f'Trait filtering removed {filtered_rows} rows ({original_total_rows} originally).')
    return filtered_df


def sanity_check_duplicate_variants_(df):
    duplicate_variants = df.groupby('Variant and risk allele').filter(lambda x: len(x) > 1)[
        'Variant and risk allele'].unique()
    all_good = True
    for variant in duplicate_variants:
        all_mapped_genes = df[df['Variant and risk allele'] == variant]['Mapped gene'].unique()
        if len(all_mapped_genes) > 1:
            print(f'Found variant, {variant}, with differing mapped gene values.')
            all_good = False

    if all_good:
        print('No repeated variants with differing mapped gene values.')


def normalize_mapped_genes_(df):
    df['gene_norm'] = df['Mapped gene'].apply(lambda val: val.split(", "))
    df = df.explode('gene_norm')
    df['gene_norm'] = df['gene_norm'].apply(replace_unknown_gene_)
    return df


def replace_unknown_gene_(gene):
    """I'm pretty sure '- indicates unknown gene association."""
    # TODO: should double-check this.
    return UNKNOWN_GENE if gene == "'-" else gene


def write_to_file_(trait, df):
    # First keep only the relevant, normalized columns for brevity.
    out_df = df[['Variant and risk allele', 'P-value_norm', 'Trait(s)', 'gene_norm']]
    column_remapping = {
        'Variant and risk allele': 'variant_and_allele',
        'P-value_norm': 'p_value',
        'Trait(s)': 'trait',
        'gene_norm': 'gene',
    }
    out_df = out_df.rename(columns=column_remapping)

    output_file = f'{OUTPUT_DIR}{trait.replace(" ", "_")}_cleaned.csv'
    out_df.to_csv(output_file)
    print(f'Wrote {output_file}.')


main()
