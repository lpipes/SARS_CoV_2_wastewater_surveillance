from Bio import SeqIO

import pandas as pd 
import numpy as np 
import csv 
import subprocess

NO_GAPS = True
fname = 'filtered_virus_MSA_short.fasta' #file with strains
strain_ids, strain_nucs = read_lines(fname)

num_strains_list = [1] #, 2, 5, 10, 20, 50, len(strain_ids)]
num_reads_list = [10]# , 100, 1000, 10000]

len_read = 100
seed = 2

for n_strain in num_strains_list:
    for n_reads in num_reads_list:
        create_tests(num_strains = n_strain \
             len_read = len_read, \
             num_reads = n_reads, seed = seed)
                                 

def read_lines(fname):
    with open(fname) as fasta_file:  # Will close handle cleanly
        strain_ids = []
        strain_nucs = []
        strain_nucs_og = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            strain_ids.append(seq_record.id)
            if NO_GAPS:
                strain_nucs.append(str(seq_record.seq).replace('-', ''))
            else:
                strain_nucs.append(str(seq_record.seq))
        
    return strain_ids, strain_nucs


def write(value, fname='rfile', append=True):
    permissions = 'w+'
    if append: 
        permissions = 'a+'
        
    with open(fname, permissions) as file:
        file.write(value + "\n")
        

def create_tests(num_strains, proportion = None, len_read = 300, num_reads = 1000, seed = 5):
    names = []

    np.random.seed(seed)
    indices = np.random.randint(size = num_strains, low = 0, high = len(strain_ids)) ## Create vector of random indices 

    case_num = f"t{num_strains}_{num_reads}_{len_read}_{seed}"
    rname = f"test_input/reads_case{case_num}.fasta"
    sname = f"test_input/strains_case{case_num}.txt"
    
    total_prop = sum(proportion)

    for i in range(num_strains): 
        s = strain_nucs[indices[i]]
        n = strain_ids[indices[i]]

        write(f"{i}:{n}", fname=sname, append=True)
        
        n_reads = num_reads if proportion == None else int(num_reads * proportion[i]) 
        
        for j in range(n_reads): 
             # random start position of the read in the strain
            idx = np.random.randint(low = 0, high = (len(s) - len_read))
            r = s[idx:(idx + len_read - 1)] # extract the read from the strain

            name = f">r{j}s{i}p{idx}"
            if name not in names:
                # write the name of read in the format: ">r2s1p27015", 
                ## where 2 is the index of the read, 1 is the index of the strain, 
                ## 27015 is the start position of the read in the strain
                write(name, fname=rname, append=True)
                 # write the read to the file
                write(r, fname=rname, append=True)
                names.append(name)

