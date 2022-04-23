GTEx data is quite large, and already hosted by... GTEx, so will be excluded from this repo.

Instructions for setup:

1.    On the [datasets](https://gtexportal.org/home/datasets) page, find the 'Single-Tissue cis-QTL Data' section.
1.    Download the eGene and significant variant-gene associations data: the file is named "GTEx_Analysis_v8_eQTL.tar"
1.    Extract the tarball contents into the folder 'data/gtex/raw'
1.    Create a folder 'data/gtex/clean'
1.    Execute the script 'src/extract_significant_variants.sh' which will unzip the significant variant files, parse out
      the variants only, and write them to text files. This is to reduce overall file size and speed up processing for
      later use.
