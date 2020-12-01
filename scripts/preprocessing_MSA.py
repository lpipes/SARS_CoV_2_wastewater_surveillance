#!/usr/bin/python3
from Bio import SeqIO
import pandas as pd, numpy as np, time, sys

strain_col = 'SN'
nuc_col = 'read'
wuhanHu1 = "NC_045512v2"
q_name = "qName"
t_start = "tStart"
t_end = "tEnd"

output_reads_sequence_col = 9

t_first = time.time()

reference_msa_filepath = sys.argv[1]
output_filepath = sys.argv[2]
output_reads_filepath = sys.argv[3]

t0 = time.time()
output_table = pd.read_csv(output_filepath, sep="\t", header=0)
print("Read in output", time.time() - t0)

t0 = time.time()
output_reads_table = pd.read_csv(output_reads_filepath, sep="\t", header=None, skiprows=3)
print("Read in output reads", time.time() - t0)

t0 = time.time()
strain_nucs = {}
with open(reference_msa_filepath) as fasta_file:  # Will close handle cleanly
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        strain_nucs[seq_record.id] = str(seq_record.seq)
print("Read in reference msa", time.time() - t0)

wuhanHu1_sequence = strain_nucs[wuhanHu1]

###############Preprocess strains so that we know what areas have same values or not 
seq = strain_nucs[list(strain_nucs.keys())[0]]
same_indices = [str(i) for i in range(len(seq)) if seq[i] == wuhanHu1_sequence[i]]
preprocessing_path = "strain_similarities_MSA"
with open(preprocessing_path, "w+") as fo:
	for id in strain_nucs.keys():
		seq = strain_nucs[id]
		sim = [str(i) for i in range(len(seq)) if seq[i] == wuhanHu1_sequence[i]]
		for s in same_indices:
			if s not in sim:
				same_indices.remove(s)
		fo.write(id + ",".join(sim))

preprocessing_path = "wuhan_similarities_MSA"
with open(preprocessing_path, "w+") as fo:
	fo.write(",".join(same_indices))