#!/usr/bin/Rscript
args <- commandArgs(TRUE)

fname = 'filtered_virus_MSA.fasta' #file with strains
tmp = readLines(fname)
strains = tmp[seq(2,length(tmp),2)] #viral strains

names(strains) = tmp[seq(1,length(tmp),2)]
names(strains) = gsub(pattern = ">", replacement = "", names(strains))

num_strains = 1 #args[1] #3 #number of the strains the reads should be generated from
num_reads = 25000 # args[2] # 1000 #number of reads
len_read = 300 #length of the read

case_num = "t1-25000"
fname = paste0("test_input/reads_case", case_num, ".fasta") #"~/readsNCase1.fasta"#name of the file with reads
r_file = file(fname, "w")
fname = paste0("test_input/strains_case", case_num, ".txt")
s_file = file(fname, "w")

indices = runif(num_strains, min = 0, max = length(names(strains))) ## Create vector of random indices 
print(indices)

for (i in 1:num_strains){
  s = strains[indices[i]]
  n = names(strains)[i]
  print(n)
  write(paste0(i, ": ", n), file=s_file, append=T)
  for (j in 1:num_reads){
    idx = sample(1:(nchar(s)-len_read), 1) #random start position of the read in the strain
    r = substr(s,idx,idx+len_read-1) #extract the read from the strain
    write(paste0(">r", j, "s", i, "p", idx), file=r_file, append=T) #write the name of read in the format:
#>r2s1p27015, where 2 is the index of the read, 1 is the index of the strain, 27015 is the start position of the read in the strain
    write(r, file=r_file, append=T)#write the read to the file
    }

  }

  close(r_file) 
