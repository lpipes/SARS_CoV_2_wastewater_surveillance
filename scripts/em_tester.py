from Bio import SeqIO

import pandas as pd 
import numpy as np 
import csv 
import subprocess

fname = 'filtered_virus_MSA_short.fasta' #file with strains
strain_ids, strain_nucs_fixed, strain_nucs_og = read_lines(fname)

a = 0.75/9
b = 0.1/9

proportions_list = [[[1]], 
             [[0.5, 0.5], [0.75, 0.25], [0.9, 0.1]],
             [[0.33, 0.33, 0.33], [0.5, 0.25, 0.25], [0.9, 0.05, 0.05]],
             [[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1], 
                  [0.25, a, a, a, a, a, a, a, a, a], 
                  [0.9, b, b, b, b, b, b, b, b, b]]]

c = [0.01, 0.1, 0.75, 0.9]
for i in range(len(c)):
    props = [c[i]]
    val = (1-c[i])/99
    props.extend([val for j in range(99)])
    
# Get sequences ignoring gaps
num_strains_list = [1] #[1, 2, 3, 5, 10, 100]
len_read_list = 300
num_reads_list = [100] #[100, 1000, 10000, 20000]
seeds_list = [2,3]


for i in range(len(seeds_list)):
    for j in range(len(num_strains_list)):
        for k in range(len(num_reads_list)):
            if not (num_strains_list[j] == 100 and num_reads_list[k] == 100):
                a = proportions_list[j]
                for l in range(len(a)):
                    print(i, j, k)
                    print(seeds_list[i], num_strains_list[j], num_reads_list[k])
                    create_tests(num_strains = num_strains_list[j], proportion = a[l], \
                                 len_read = len_read, \
                                 num_reads = num_reads_list[k], seed = seeds_list[i])
                                 

def read_lines(fname):
    with open(fname) as fasta_file:  # Will close handle cleanly
        strain_ids = []
        strain_nucs_fixed = []
        strain_nucs_og = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            strain_ids.append(seq_record.id)
            strain_nucs_fixed.append(str(seq_record.seq).replace('-', ''))
            strain_nucs_og.append(str(seq_record.seq))
        
    return strain_ids, strain_nucs_fixed, strain_nucs_og
    

def write(value, fname='rfile', append=True):
    permissions = 'w+'
    if append: 
        permissions = 'a+'
        
    with open(fname, permissions) as file:
        file.write(value + "\n")
        

def create_tests(num_strains, proportion, len_read = 300, num_reads = 1000, seed = 5):
    names = []

    np.random.seed(seed)
    indices = np.random.randint(size = num_strains, low = 0, high = len(strain_ids)) ## Create vector of random indices 

    case_num = f"t{num_strains}_{num_reads}_{seed}"
    rname = f"test_input/reads_case{case_num}.fasta"
    sname = f"test_input/strains_case{case_num}.txt"

    for i in range(num_strains): 
        s = strain_nucs_fixed[indices[i]]
        n = strain_ids[indices[i]]

        print(n)
        write(f"{i}:{n}", fname=sname, append=True)
        
        print(num_strains, num_reads, proportion[i], proportion)
        n_reads = int(num_reads * proportion[i])
        
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
                with open(f"{case_num}.out", 'w') as f:
                    subprocess.call(['Rscript', 'EM_new.R', rname])

