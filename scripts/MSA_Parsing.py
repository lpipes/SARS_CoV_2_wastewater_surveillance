#!/usr/bin/python3
from Bio import SeqIO
import pandas as pd, numpy as np, time, sys, re

strain_col = 'SN'
nuc_col = 'read'
wuhanHu1 = "NC_045512v2"
q_name = "qName"
t_start = "tStart"
t_end = "tEnd"
q_start = "qStart"
q_end = "qEnd"
sam_start = 3
output_reads_sequence_col = 9
min_num_pos = 2; freq_tresh = 0.03

t_first = time.time()

reference_msa_filepath = sys.argv[1]
output_filepath = sys.argv[2]
output_reads_filepath = sys.argv[3]
similarities_filepath = sys.argv[3]

t0 = time.time()
output_table = pd.read_csv(output_filepath, sep="\t", header=0)
print("Read in output", time.time() - t0)

t0 = time.time()
output_reads_table = pd.read_csv(output_reads_filepath, sep="\t", header=None, skiprows=3)
print("Read in output reads", time.time() - t0)
# print("output reads table: ", output_reads_table)

t0 = time.time()
with open(reference_msa_filepath) as fasta_file:  # Will close handle cleanly
    strain_ids = []
    strain_nucs = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        strain_ids.append(seq_record.id)
        strain_nucs.append(str(seq_record.seq))
# print(type(strain_nucs[0]))
ref_msa_table = pd.DataFrame(np.column_stack([strain_ids, strain_nucs]), columns=['SN', 'read'])
print("Read in reference msa", time.time() - t0)

t0 = time.time()
with open(similarities_filepath, "r") as similarities_file:
	similarities = similarities_file.read().split(",")
print("Read in similarities between read", time.time() - t0)

# print(ref_msa_table)
# print(output_table)

## Finding allele freequencies for reads using alignment and mismatches 
## prepare pandas dataframe and write in data 
## get wuhanHu1 to compare 
wuhanHu1_sequence = ref_msa_table.loc[ref_msa_table[strain_col] == wuhanHu1][nuc_col]
wuhanHu1_sequence = wuhanHu1_sequence.iloc[0]
read_seq = wuhanHu1_sequence
# print("strains: ", ref_msa_table[strain_col])
# print("wuhanHu1: ", ref_msa_table.loc[ref_msa_table[strain_col] == wuhanHu1])
# print("wuhanHu1: ", wuhanHu1_sequence)

def convert_bp_to_array(bp):
	return_array = np.array([0, 0, 0, 0, 1]) #Final element of array indicates total number of reads overlapping this position
	if bp == 'A':
		return_array[0] += 1
	if bp == 'C':
		return_array[1] += 1
	if bp == 'G':
		return_array[2] += 1
	if bp == 'T':
		return_array[3] += 1
	return return_array 

pd.set_option('display.max_columns', None)
pd.set_option("display.max_colwidth", 10000000)
pd.set_option('max_seq_item', None)

def translate_seq(m):
	return '0' if m == 'A' else ('1' if m == 'C' else 
		('2' if m == 'G' else ('3' if m == 'T' else '4')))

## Need to preprocess the strains to figure out where everything is the same 
def find_matches(output_row):
	o_tPositions = "tPositions"
	alignments = output_row[o_tPositions] 
	alignments = list(map(int, alignments.split(","))) #convert string to int array
	# print("type of alignment: ", type(alignments))
	# total_matches = 0
	# total_mismatches = 0
	# read_seq = pd.Series.to_string(read_seq)
	# read_seq = read_seq.strip().split(" ")[-1] #remove any spaces in beginning
	# print("type of read_seq: ", type(read_seq), len(read_seq))
	# print("read_seq: ", read_seq)
	allele_freq = np.zeros((len(read_seq), 5))
	# print("alignment: ", alignments)

	# i = 0
	# last_match = output_row[t_start]
	# last_match = output_row[q_start] ########################################### t_start too big??? using alignment start
	start = output_reads_row[sam_start]
	# print(type(start), "start", len(start))
	start = start.iloc[0]
	# print(type(start), "start", start)
	end = start+output_row[q_end]-output_row[q_start]+1
	for m in range(start, end): #m in alignment: ## will alignment work in case of gaps? #
		# check if reference alignment with read_seq at location m 
		# m = alignments[i]
		# val = m - last_match
		# print(m, ", ", val)
		# print(len(read_seq))
		# print(len(alignments))
		# print("read seq m: ", read_seq[i])
		# if val in alignment:
		if m in alignments and m not in similarities:
			allele_freq[m] = convert_bp_to_array(read_seq[i])
			#np.append(allele_freq, convert_bp_to_array(read_seq[i]))
		# last_match = m
		# i = val
		else: #ignore mismatches for now
			allele_freq[m] = np.array([0, 0, 0, 0, 1])
			#np.append(allele_freq, np.array([0, 0, 0, 0, 1]))
		# if read_seq[m] == wuhanHu1_sequence[m]:
		# 	total_matches += 1
		# else:
		# 	total_mismatches += 1
	return allele_freq # total_matches, total_mismatches, allele_freq

# all_allele_freq = {}
total_freq_nuc = np.zeros((len(wuhanHu1_sequence), 5))
t0 = time.time()
for i in range(len(output_table)):
	output_row = output_table.iloc[i]
	output_reads_row = output_reads_table.loc[output_reads_table[0] == output_row[q_name]]
	# print("a: ", a)
	# read_seq =  a[output_reads_sequence_col]
	allele_freq = find_matches(output_row) #alignment, mismatches, 
	# all_allele_freq[read_name] = allele_freq
	# print("len wuhan: ", len(wuhanHu1_sequence), type(wuhanHu1_sequence))
	# print(len(total_freq_nuc))
	# total_freq_nuc[output_row[t_start], output_row[t_end]] += allele_freq
	np.add(total_freq_nuc, allele_freq)

print("Calculate frequencies for each read (# reads = ", len(output_table), ")", time.time() - t0)

## Actually calculate frequencies using sum
print("total_freq_nuc: ", total_freq_nuc)
print("total_freq_nuc size: ", len(total_freq_nuc))
print("total_freq_nuc size: ", total_freq_nuc.shape)

freq_nuc = np.zeros((len(wuhanHu1_sequence), 5))
for row in range((len(wuhanHu1_sequence))):
	# print("total freq_nuc row ", total_freq_nuc[row].shape)
	# print("total freq_nuc ", total_freq_nuc[row].shape)
	a = np.divide(total_freq_nuc[row], total_freq_nuc[row, 4])
	print(a.shape)
	freq_nuc[row] = a

## Removing uncessary sequeces using allele frequencies
acceptable = ref_msa_table.keys().copy()
eliminated = set()
t0 = time.time()
for index, row in ref_msa_table.iterrows(): 
	strain, strain_nuc = row[strain_col], row[nuc_col]
	strain_nuc_list = []
	strain_nuc = strain_nuc.upper()
	for m in strain_nuc:
		strain_nuc_list.append(translate_seq(m))
	strain_nuc = np.array(strain_nuc_list, dtype=np.int8)
	print(strain_nuc, type(strain_nuc), len(strain_nuc))
	if strain != wuhanHu1:
		freq_nuc_i = freq_nuc[range(freq_nuc.shape[0]), strain_nuc]
		if sum(freq_nuc_i < freq_tresh) >= min_num_pos:
			eliminated.add(strain)
			acceptable.remove(strain)
print("Filter strains (# strains = ", len(ref_msa_table), ")", time.time() - t0)

print("Eliminated", eliminated)

f = open("eliminated_new_test.txt", "a")
f.write(eliminated)
f.close()
