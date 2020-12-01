library(Matrix);library(plyr); library(data.table); library(dplyr); library(stringr); 
library(turboEM); library(hash); 
setwd('~/')


args <- commandArgs(trailingOnly = TRUE)

case_num = "sample" #<- as.numeric(strsplit(grep('--case*', args, value = TRUE), split = '=')[[1]][[2]])
output_file = "output.txt" #<- toString(strsplit(grep('--output_file*', args, value = TRUE), split = '=')[[1]][[2]])
msa_file = "reference_virus_MSA.fasta" #<- toString(strsplit(grep('--msa_file*', args, value = TRUE), split = '=')[[1]][[2]])

syspath = "/Users/selina/Desktop/sars-Cov-2/Viral-Strains/" #/space/s1/selina/

cat("case ", case_num)
cat("output file ", file)
cat("msa file ", msa_file)

#--------------------------load viral strains---------------------------------
#fname = '10_strains.fasta'
ptm <- proc.time() 
fname = paste(syspath, msa_file, sep="") 
tmp = readLines(fname)
strains = tmp[seq(2,length(tmp),2)] #viral strains
names(strains) = tmp[seq(1,length(tmp),2)] ## set names
names(strains) = gsub(pattern = ">", replacement = "", names(strains)) ## fix names


cat("\nLoad in sequences:", proc.time() - ptm)
ptm = proc.time()
alignments_file = "sample.sam" #paste('alignments_case_', case_num, '.sam', sep="")
fname=paste(syspath, alignments_file, sep="") 

csv_alignments = read.csv(fname, sep='\t', skip=0, header=F)

df_align = setNames(csv_alignments, c('Type','SN','LN'))
df_align_pg = subset(df_align, Type=='@PG')
SN = sapply(df_align$SN[df_align$Type=='@SQ'], function(s) strsplit(toString(s),":")[[1]][2], USE.NAMES=F) # sequence names
names(SN) = 1:length(SN)
SID = setNames(1:length(SN), SN)
nSkip = sum(grepl("@", df_align$Type))


cat("\nLoad in alignments:", proc.time() - ptm)
ptm = proc.time()
df_align = setNames(read.csv(fname, sep='\t', skip=nSkip#length(SN)+length(df_align_pg) 
			     +2, header=F), c('read','SN','LN'))
df_align = setNames(df_align[,c(1,3,4,6,10)], c('read','SN','pos','cigar','seq'))
df_align$SID = unname(SID[df_align$SN])
df_align$SN = NULL

cat("\nSet names of alignments:", proc.time() - ptm)
ptm = proc.time()
#-------------------------- create matrix Q (Pr(read_i | read_i from strain_j))----
#calculate the alignment length as a sum of numbers in CIGAR field
#and the number of differences as the sum of numbers before character "X" in CIGAR field
isEmpty <- function(a) {
  return(length(a)==0)
}

n = unname(sapply(df_align$cigar, function(a) sum(as.numeric(str_extract_all(a, "[0-9]+")[[1]]))))

d = unname(sapply(df_align$cigar, function(a) {
  x = str_extract_all(a, "[^0-9]+")[[1]]
  y = as.numeric(str_extract_all(a, "[0-9]+")[[1]])
  if (!isEmpty(y) && !isEmpty(x)) {
    v = setNames(as.numeric(str_extract_all(a, "[0-9]+")[[1]]), str_extract_all(a, "[^0-9]+")[[1]])
    sum(v[names(v)=='X'])
  }
}))
nd_mtx = cbind(n=n,d=d)
dfN = cbind(df_align[,c(1,5)], n=c(nd_mtx[,1]))
dfD = cbind(df_align[,c(1,5)], d=c(nd_mtx[,2]))
print("dfD")
dfN = reshape(dfN, direction = "wide", idvar = "read", timevar = "SID")
dfD = reshape(dfD, direction = "wide", idvar = "read", timevar = "SID")
colnames(dfN)[2:ncol(dfN)] = substring(colnames(dfN)[2:ncol(dfN)], 3)
colnames(dfD)[2:ncol(dfD)] = substring(colnames(dfD)[2:ncol(dfD)], 3)
sum(is.na(dfN))
sum(is.na(dfD))

#if there are NA values, then fill them
dfN[is.na(dfN)] = 300
dfD[is.na(dfD)] = 300

n_mtx = data.matrix(dfN[,2:ncol(dfN)])
d_mtx = data.matrix(dfD[,2:ncol(dfD)])
# print('created matrices')
k = ncol(n_mtx) #number of viral strains

cat("\nSet up for creating Q matrix", proc.time() - ptm)
ptm = proc.time()
#calculate the likelihood P(X=x_i|x_i \in y_j) (matrix Q)
error = 0.01   # error
Q = ((error ** d_mtx) * (1 - error) ** (n_mtx - d_mtx)) ** (1/n_mtx)

cat("\nCreate Q-matrix", proc.time() - ptm)
ptm = proc.time()
#-------------------------------------TurboEM---------------------
#objective function
logL = function(p, y) return(sum(log(colSums(t(y)*p))))

#EM
EM = function(p, y){
  w = t(t(y) * p)    # numerator of Eq.12 in Malone's paper
  E_step = w/rowSums(w) # denominator of Eq.12 in Malone's paper (normalization of rows)
  #return(colSums(E_step)/nrow(y))# M-step to compute p, Eq.33
  return(colMeans(E_step))
}

p0 = runif(k)
p0 = p0 / sum(p0)
res = turboem(par=p0, fixptfn=EM, objfn=logL, method="squarem", y=Q)
options(digits=13)
res
p = pars(res) #[c(2),]

par(mar=c(13,4,4,1))

names.arg = paste0(tolower(names(SID)[as.numeric(colnames(n_mtx))]), ", ", 1:k)

em_output_df = cbind(p, names.arg)
em_output_df = as.data.frame(em_output_df)
colnames(em_output_df) = c("p", "names.arg")

print("EM Output")
print(head(em_output_df))
write.csv(em_output_df, file = paste(syspath, "em_output_", case_num, ".csv", sep=""))

d <- data.table(em_output_df, key="p")
d <- tail(d, 50)

print("k:")
print(k)
print(length(names(SID)))
pdf(paste(syspath, "rplot_case_", case_num, ".pdf", sep="")) 
barplot(as.numeric(factor(d[["p"]])), beside = T, ylab = "p vec values", names.arg=d[["names.arg"]],
        cex.names=0.8, las=2, main='Viral strains estimation')
dev.off()

#grid()
cat("\nRun EM algorithm", proc.time() - ptm)
cat("\n")

