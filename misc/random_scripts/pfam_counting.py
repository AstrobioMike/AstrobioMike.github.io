#!/usr/bin/env python3

## this is not meant for widespread use, it is a helper script in a process of generating single-copy gene HMM sets for GToTree (https://github.com/AstrobioMike/GToTree/wiki/What-is-GToTree%3F)
## it's usage is shown in the code presented here: https://github.com/AstrobioMike/GToTree/wiki/SCG-sets#code


import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="This script is for counting pfam hits to genomes' gene files.")

required = parser.add_argument_group('required arguments')

required.add_argument("-p", "--input_pfam_list", help="Single-column file with all pfam accessions", action="store", dest="input_pfam_accs", required=True)
required.add_argument("-g", "--input_genome_ids", help="Single-column file with all genome IDs", action="store", dest="input_genome_ids", required=True)
required.add_argument("-H", "--input_hmm_results_tab", help="HMM results tab", action="store", dest="hmm_results", required=True)

parser.add_argument("-o", "--output_counts_tsv", help='Output file with pfam-hit counts per genome. (default: "pfam_counts.tsv")', action="store", dest="output_file", default="pfam_counts.tsv")

args = parser.parse_args()

acc_list = []

with open(args.input_pfam_accs) as accs:
    for line in accs:
        acc = line.strip()
        acc_list.append(acc)

genome_list = []

with open(args.input_genome_ids) as genomes:
    for line in genomes:
        genome = line.strip()
        genome_list.append(genome)

array = np.zeros( (len(genome_list),len(acc_list)), dtype=int )
colnames = acc_list
rownames = genome_list
df = pd.DataFrame(array, index=rownames, columns=colnames)

num = 0

with open(args.hmm_results) as hits:
    for line in hits:
        num += 1
        if num % 1000 == 0:
            print(num)

        line = line.strip()

        if line.startswith("#"):
            continue

        line_list = line.split()
        curr_acc_list = line_list[0].split('_')[:2]
        curr_acc = "_".join(curr_acc_list)
        curr_pfam = line_list[3]

        df.loc[curr_acc,curr_pfam] = df.loc[curr_acc,curr_pfam] + 1

df.to_csv(args.output_file, index=True, header=True, sep="\t")
