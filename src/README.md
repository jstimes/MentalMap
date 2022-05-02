These scripts clean data, load more data, and join data together. They are expected to be run from this directory with the data directory as configured in this project.

Order of processing:

1.    Run `gwas_catalog_data_cleaner.py` to convert the [data/gwas/raw/](https://github.com/jstimes/MentalMap/tree/main/data/gwas/raw) data into [data/gwas/clean/](https://github.com/jstimes/MentalMap/tree/main/data/gwas/clean) format.
2.    Run `extract_significant_variants.sh` (after downloading the GTEx Portal data as described [here](https://github.com/jstimes/MentalMap/blob/main/data/gtex/README.md))
3.    Run `data_joiner.py` to generated the [data/joined/](https://github.com/jstimes/MentalMap/tree/main/data/joined) dataset using the clean GWAS data, dbSNP data (to be cached to [data/af](https://github.com/jstimes/MentalMap/tree/main/data/af) when run), and GTEx Portal data.

`extract_tissues.sh` is an optional script to generate a [list](https://github.com/jstimes/MentalMap/blob/main/data/gtex/tissue_metadata.txt) of all tissues to help with analysis.
