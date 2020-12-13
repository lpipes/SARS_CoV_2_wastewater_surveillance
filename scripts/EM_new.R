library(Matrix);library(plyr); library(data.table); library(dplyr); library(stringr); 
library(turboEM); library(hash); library(feather); library(glue); library(arrow); 
setwd('~/')

#create matrix of alignments' lengths and matrix of numbers of differences
# dfN = cbind(df_align[,c(1,5)], n=c(nd_mtx[,1]))
# dfD = cbind(df_align[,c(1,5)], d=c(nd_mtx[,2]))

# ##A ‘wide’ longitudinal dataset will have one record for each individual with some time-constant variables that occupy single columns and some time-varying variables that occupy a column for each time point. 
# ##for ex: reads have multiple duplicates, so then reshape will place duplicates into the same row and exand the columns; for SID, only the first one will be taken
# ## row = reads are unique as id var  and col = SID changes per read 
ptm = proc.time()
#syspath = "/space/s1/selina/SARS_CoV_2_wastewater_surveillance/scripts/test_input/"#"test_input/test_input/"#/space/s1/selina/SARS_CoV_2_wastewater_surveillance/scripts/test_input" #Viral-Strains/" #~/Desktop/sars-Cov-2/SARS_CoV_2_wastewater_surveillance/
syspath = "~/Desktop/sars-Cov-2/SARS_CoV_2_wastewater_surveillance/scripts/test_input/"

dfD_path = glue("{syspath}dfD.feather")
dfN_path = glue("{syspath}dfN.feather")
SID_path = glue("{syspath}SID.feather")

print(dfD_path)
dfD = arrow::read_feather(dfD_path)# compression='uncompressed')
print(dfN_path)
dfN = arrow::read_feather(dfN_path)# compression='uncompressed')
print(SID_path)
SID = arrow::read_feather(SID_path)#, compression='uncompressed')
cat("\nRead in files", proc.time() - ptm)
ptm = proc.time()
# dfN = reshape(dfN, direction = "wide", idvar = "read", timevar = "SID")
# dfD = reshape(dfD, direction = "wide", idvar = "read", timevar = "SID")

##replace column names to remove the d.
colnames(dfN)[2:ncol(dfN)] = substring(colnames(dfN)[2:ncol(dfN)], 3)
colnames(dfD)[2:ncol(dfD)] = substring(colnames(dfD)[2:ncol(dfD)], 3)

## get the total number of na values 
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

