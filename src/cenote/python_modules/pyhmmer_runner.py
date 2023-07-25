#!/usr/bin/env python

import subprocess
import os
import sys
import pyhmmer
from pyhmmer import hmmscan as hmmscan
import pandas as pd
import multiprocessing.pool
import math
import re
import time

input_dir = sys.argv[1]

out_dir = sys.argv[2]

which_DB = sys.argv[3]


if not os.path.isdir(out_dir):
    os.makedirs(out_dir)




print("pyhmmscan pool part")

starttime = time.time()

def hmmscanner(seqs):
    scanout = list(hmmscan(pyhmmer.easel.SequenceFile(seqs, digital=True), pyhmmer.plan7.HMMFile(which_DB)))
    return scanout



splitAA_list = []
for splitAA in os.listdir(input_dir):
    if splitAA.endswith('.faa'):
        f = os.path.join(input_dir, splitAA)

        if os.path.isfile(f):
            splitAA_list.append(f)

if not splitAA_list:
    print("no files found for pyhmmer in " + str(input_dir))
    exit

hmmscan_list = []
with multiprocessing.pool.ThreadPool() as pool:
    for alignments in pool.map(hmmscanner, splitAA_list):
        for model in alignments:
            quer1 = model.query_name.decode()
            for hit in model:
                target_name = hit.name.decode()
                target_acc = hit.accession
                full_seq_evalue = hit.evalue
                seq_pvalue = hit.pvalue        
                hmmscan_list.append([quer1, target_name, full_seq_evalue, seq_pvalue])

hmmscan_pools_df = pd.DataFrame(hmmscan_list, columns=["ORFquery", "target", "evalue", "pvalue"])\
    .sort_values('evalue').drop_duplicates('ORFquery').query("evalue <= 1e-8")


hmmscan_output_file = os.path.join(out_dir, "pyhmmer_orig_split_AAs.tsv")

hmmscan_pools_df.to_csv(hmmscan_output_file,
                        sep = "\t", index = False)

hmmscan_pools_df["pos"] = hmmscan_pools_df["ORFquery"].str.rfind("_")

hmmscan_pools_df["contig"] = hmmscan_pools_df.apply(lambda x: x['ORFquery'][0:x['pos']],axis=1)

hmmscan_contig_sum = hmmscan_pools_df.groupby("contig").size().reset_index(name='count')

contig_sum_file = os.path.join(out_dir, "contig_hallmark_count.tsv")

hmmscan_contig_sum.to_csv(contig_sum_file,
                        sep = "\t", index = False)

endtime = time.time()

time_taken = endtime - starttime

print("pyhmmscan part took: " + str(time_taken) + " seconds")