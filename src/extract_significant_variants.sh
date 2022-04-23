#!/bin/bash

gtex_dir=../data/gtex/raw/
out_dir=../data/gtex/clean/
for filepath in "$gtex_dir"/*signif_variant_gene_pairs.txt.gz
do
  file="${filepath#*//}"
  tissue="${file%%.*}"
  out_file="$out_dir$tissue.significant_variant_locations.txt"
  gzip -cd "$filepath" | awk '{print $1 }' > "$out_file"
  echo "Wrote $out_file"
done

