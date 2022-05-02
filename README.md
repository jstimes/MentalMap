# MentalMap
Clustering-based meta-analysis of GWAS summary statistics for mental health conditions.

Data comes from [GWAS Catalog](https://www.ebi.ac.uk/gwas/home), [GTEx Portal](https://gtexportal.org/home/datasets), and [dbSNP](https://www.ncbi.nlm.nih.gov/snp/).

There are three main components of this project:

1.   Exploratory analysis of raw data from GWAS catalog and auxilliary sources.
     1.    GWAS summary statisics were retrieved from the links specified in the [gwas_metadata.csv file](https://github.com/jstimes/MentalMap/blob/main/data/gwas/gwas_trait_metadata.csv) and the data for each condition is stored in [data/gwas/raw](https://github.com/jstimes/MentalMap/tree/main/data/gwas/raw).
     2.    Data for a single trait was explored in the notebook [GWAS_catalog_EDA_and_normalization](https://github.com/jstimes/MentalMap/blob/main/notebooks/GWAS_catalog_EDA_and_normalization.ipynb) to understand the format, properties of the data, and certain limitations/inconsistencies.
     3.    Additional data from dbSNP (allele frequencies) and GTEx Portal (significant variant-gene associations in tissues) were retrieved and joined with the GWAS summary statistic data. See the GTEx [README](https://github.com/jstimes/MentalMap/blob/main/data/gtex/README.md) for more details.
     4.    The exploratory cleaning and data-joining was converted to scripts in the [src folder](https://github.com/jstimes/MentalMap/tree/main/src) and used to generate the [joined data](https://github.com/jstimes/MentalMap/tree/main/data/joined).
2.   Analysis of cleaned/normalized data in the notebook [GWAS_post_cleaning_analysis](https://github.com/jstimes/MentalMap/blob/main/notebooks/GWAS_post_cleaning_analysis.ipynb).
3.   Encoding and clustering exploration in the notebook [GWAS_encode_and_cluster](https://github.com/jstimes/MentalMap/blob/main/notebooks/GWAS_encode_and_cluster.ipynb).     
