*    `af`: allele frequency data from dbSNP. In CSV format where one column defines variant, other is the allele frequency.
*    `gtex`: eQTL data from GTEx. See inner readme for more details.
*    `gwas`: summary statistics from GWAS catalog. See metadata file for exact retrieval details. `gwas/raw/` stores data that came directly from GWAS catalog; `gwas/clean/` stores normalized, filtered summary statistics.
*    `joined`: this is the cleaned summary statistic data joined with allele frequency info from dbSNP and tissue association data from GTEx.

