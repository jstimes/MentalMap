#!/bin/bash

gtex_dir=../data/gtex/clean/
out_file=../data/gtex/tissue_metadata.txt
for filepath in "$gtex_dir"/*.significant_variant_locations.txt
do
  file="${filepath#*//}"
  tissue="${file%%.*}"
  echo $tissue >> $out_file
done

