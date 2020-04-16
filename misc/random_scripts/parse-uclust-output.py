#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='This script parses the uclust-formatted output table from `vsearch --derep_fulllength ... --uc output.uc` into a new table with 3 columns: 1) representative sequence; 2) # of seqs represented; 3) all original seq IDs that are represented by this representative sequence, delimited by ";".')

required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--uclust_tsv", help='input uclust-formatted table', action="store", dest="uclust_tsv", required=True)

parser.add_argument("-o", "--output_tsv", help='Output tsv as described in main help (default: "derep-map.tsv")', action="store", dest="output_tsv", default="derep-map.tsv")

args = parser.parse_args()

# reading in uclust table
uclust = pd.read_csv(args.uclust_tsv, sep="\t", names=["record_type", "cluster_num", "length_or_cluster_size", "identity_with_rep", "strand_with_rep", "query", "rep"], usecols=[0,1,2,3,4,8,9])

rep_seqs_dict = {}

# adding rep seq ids to initial rep_seqs_dict as key, and values will be list of all ids represented by that rep seq, starting now with just the one that was chosen as the rep seq
for index, row in uclust.iterrows():
    if row['record_type'] == "S":
        rep_seqs_dict[row['query']] = [row['query']]

# now adding the ids represented by each rep sequence
for index, row in uclust.iterrows():
    if row['record_type'] == "H":
        rep_seqs_dict[row['rep']].append(row['query'])

# writing out as table
  # header
out_tab = open(args.output_tsv, "w")
out_tab.write("representative_seq_id" + "\t" + "number_seqs_represented" + "\t" + "seq_ids_represented" + "\n")

for key in rep_seqs_dict:
    out_tab.write(str(key) + "\t" + str(len(rep_seqs_dict[key])) + "\t" + ";".join(rep_seqs_dict[key]) + "\n")

out_tab.close()
