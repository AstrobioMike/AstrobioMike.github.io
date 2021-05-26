#!/usr/bin/env bash

## this is not meant for widespread use, it is a helper script in a process of generating single-copy gene HMM sets for GToTree (https://github.com/AstrobioMike/GToTree/wiki/What-is-GToTree%3F)
## it's usage is shown in the code presented here: https://github.com/AstrobioMike/GToTree/wiki/SCG-sets#code

rm -rf all_genes.faa # in case file exists already

for assembly in $(cat $1)
do

  gunzip ${assembly}_protein.faa.gz

  bit-rename-fasta-headers -i ${assembly}_protein.faa -w $assembly -o genes.tmp # this just renames the headers and appends a number to keep them unique (https://github.com/AstrobioMike/bioinf_tools)

  cat genes.tmp >> all_genes.faa

  rm genes.tmp ${assembly}_protein.faa

done < $1
