#!/usr/bin/python3
import pandas as pd, numpy as np, time, sys, re
import feather as ft

ref_msa_path = 'test_input/reference_virus_MSA.feather'
output_mismatch_path = "test_input/output_mismatches.feather"
case = "test"
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
similarities_filepath = sys.argv[4]

t0 = time.time()
output_table = pd.read_csv(output_filepath, sep="\t", header=0)
print("Read in output", time.time() - t0)

t0 = time.time()
output_reads_table = pd.read_csv(output_reads_filepath, sep="\t", header=None, skiprows=3)
print("Read in output reads", time.time() - t0)


t0 = time.time()
output_mismatch_table = pd.read_csv(output_reads_filepath, sep="\t", header=None, skiprows=3)
print("Read in output mismatch reads", time.time() - t0)
ft.write_dataframe(output_mismatch_table, output_mismatch_path)

######## parse using new lines and >
#############> name of sequence then new line and sequence itself then new line 
############### save the name of the names of the reference sequence 
############## could be limit of length of sequence reading in (default buffer )
################ can manipulate matrix however i want
reference_sequences = []
reference_strains = []
t0 = time.time()
with open(reference_msa_filepath, "r") as fasta_file: 
     msa = fasta_file.read().split(">")
     for strains in msa:
            print(strains) 
            if strains != "":
                s0 = strains.split("\n")
                reference_strains.append(s0[0])
                reference_sequences.append(s0[1]) # use np.fromiter(s0[1], (np.unicode,1))) for char array
ref_msa_table = pd.DataFrame(np.column_stack([reference_strains, reference_sequences]), columns=[strain_col, nuc_col])
ft.write_dataframe(ref_msa_table, ref_msa_path)
print("Read in reference msa", time.time() - t0)

t0 = time.time()
with open(similarities_filepath, "r") as similarities_file:
        similarities = similarities_file.read().split(",")
print("Read in similarities between read", time.time() - t0)


## Finding allele freequencies for reads using alignment and mismatches 
## prepare pandas dataframe and write in data 
## get wuhanHu1 to compare 
wuhanHu1_sequence = ref_msa_table.loc[ref_msa_table[strain_col] == wuhanHu1][nuc_col]
wuhanHu1_sequence = wuhanHu1_sequence.iloc[0]
print(wuhanHu1_sequence, "wuhan", type(wuhanHu1_sequence), len(wuhanHu1_sequence))

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

total_freq_nuc = np.repeat(np.array([0, 0, 0, 0, 1])[None, :], len(wuhanHu1_sequence), axis=0)
## Need to preprocess the strains to figure out where everything is the same 
def find_matches(output_row, read_seq):
        read_seq = read_seq.iloc[0]
        o_tPositions = "tPositions"
        alignments = output_row[o_tPositions] 
        alignments = list(map(int, alignments.split(","))) 
        allele_freq = np.repeat(np.array([0, 0, 0, 0, 1])[None, :], len(wuhanHu1_sequence), axis=0)#np.zeros((len(wuhanHu1_sequence), 5))
        
        print(allele_freq)
        a = output_row["qStarts"]
        print(a, type(a))
        
        starts = output_row["qStarts"].split(",")
        block_sizes = output_row["blockSizes"].split(",")
        for j in range(len(starts)-1):
            start = int(starts[j]) 
            end = int(start) + int(block_sizes[j])
            k = 0
            for i in range(start, end):
                print(i, k)
                tPosition = alignments[k]
                #print(output_row[0], i, tPosition, np.size(allele_freq), len(read_seq), type(read_seq), len(alignments), read_seq[i], wuhanHu1_sequence[tPosition])
                # if tPosition not in similarities:
                a = convert_bp_to_array(read_seq[i])
                print("a ", a, "tPosition", tPosition, " allele_freq[tPosition]", allele_freq[tPosition])
                allele_freq[tPosition] = a
                print(output_row["qName"])
                print("allele_freq[tPosition]", allele_freq[tPosition])
                k += 1
                # else: #ignore mismatches for now
                #       allele_freq[m] = np.array([0, 0, 0, 0, 1])
        return allele_freq # total_matches, total_mismatches, allele_freq

# all_allele_freq = {}
t0 = time.time()
for i in range(len(output_table)):
        output_row = output_table.iloc[i] 
        ############## pass in output stuff to the 
        ############ pass in read to find_matches and go through those alleletes, use the posistions to figure out where to put the values in teh total sequence (care of gaps, will need to move everythig over )
        output_reads_row = output_reads_table.loc[output_reads_table[0] == output_row[q_name]]
        read_seq =  output_reads_row[output_reads_sequence_col]
        allele_freq = find_matches(output_row, read_seq) 
        print("allele freq", allele_freq)
        np.add(total_freq_nuc, allele_freq)

print("Calculate frequencies for each read (# reads = ", len(output_table), ")", time.time() - t0)

## Actually calculate frequencies using sum
print("total_freq_nuc: ", total_freq_nuc)
print("total_freq_nuc size: ", len(total_freq_nuc))
print("total_freq_nuc size: ", total_freq_nuc.shape)

print(total_freq_nuc)
freq_nuc = np.zeros((len(wuhanHu1_sequence), 5))
for row in range((len(wuhanHu1_sequence))):
        a = np.divide(total_freq_nuc[row], total_freq_nuc[row, 4])
        print("a shape: ", a.shape, a, "     total_freq_nuc[row]: ", total_freq_nuc[row])
        freq_nuc[row] = a
print("freq nuc: ", freq_nuc)

# output_n_df = ///////////////////////////////////////
# output_d_df = ///////////////////////////////////////
# ft.write_dataframe(output_df, "output_table.feather")

## Removing uncessary sequeces using allele frequencies
acceptable = set(ref_msa_table[strain_col]) #ref_msa_table.keys().copy()
print("acceptable: ", acceptable)
eliminated = set()
t0 = time.time()
for index, row in ref_msa_table.iterrows(): 
        strain, strain_nuc = row[strain_col], row[nuc_col]
        strain_nuc_list = []
        strain_nuc = strain_nuc.upper()
        for m in strain_nuc:
                strain_nuc_list.append(translate_seq(m))
        strain_nuc = np.array(strain_nuc_list, dtype=np.int8)
        #print(strain_nuc, type(strain_nuc), len(strain_nuc))
        if strain != wuhanHu1:
                freq_nuc_i = freq_nuc[range(freq_nuc.shape[0]), strain_nuc]
                print("freq_nuc_i: ", freq_nuc_i, " freq_tresh: ", freq_tresh, " min_num_pos:", min_num_pos)
                print(sum(freq_nuc_i < freq_tresh))
                if sum(freq_nuc_i < freq_tresh) >= min_num_pos:
                        eliminated.add(strain)
                        acceptable.remove(strain)
print("Filter strains (# strains = ", len(ref_msa_table), ")", time.time() - t0)

print("Eliminated", eliminated)

f = open("test_input/eliminated.txt", "a")
for e in eliminated:
        f.write(e + "\n")
f.close()


