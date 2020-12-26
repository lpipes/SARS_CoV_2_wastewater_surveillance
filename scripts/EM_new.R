library(Matrix);library(plyr); library(data.table); library(dplyr); library(stringr); 
library(turboEM); library(hash); library(feather); library(glue); library(arrow); 
setwd('~/')

case_num = "t4_5000_3"
syspath = "/space/s1/selina/SARS_CoV_2_wastewater_surveillance/scripts/test_input/"#"test_input/test_input/"#/space/s1/selina/SARS_CoV_2_wastewater_surveillance/scripts/test_input" #Viral-Strains/" #~/Desktop/sars-Cov-2/SARS_CoV_2_wastewater_surveillance/

# ##A ‘wide’ longitudinal dataset will have one record for each individual with some time-constant variables that occupy single columns and some time-varying variables that occupy a column for each time point. 
# ##for ex: reads have multiple duplicates, so then reshape will place duplicates into the same row and exand the columns; for SID, only the first one will be taken
# ## row = reads are unique as id var  and col = SID changes per read 
# dfN = reshape(dfN, direction = "wide", idvar = "read", timevar = "SID")
# dfD = reshape(dfD, direction = "wide", idvar = "read", timevar = "SID")
ptm = proc.time()
data <- read.table("/space/s1/selina/SARS_CoV_2_wastewater_surveillance/scripts/times/mismatches_output_.txt",fill = TRUE, header=T)
cat("data: ", dim(data))
d_mtx <- data[,3:ncol(data)]
reads<-data[,1]

k = ncol(d_mtx) #number of viral strains
n_mtx <- matrix(rep(data[,2],k),nrow=nrow(d_mtx),ncol=k)
colnames(n_mtx) <- colnames(d_mtx)
cat("Loading matrices in from file", proc.time() - ptm)
ptm = proc.time()
d_mtx <- d_mtx[ , colSums(is.na(d_mtx)) != nrow(d_mtx)]
n_mtx <- n_mtx[ , colSums(is.na(n_mtx)) != nrow(n_mtx)]

k = ncol(d_mtx) #number of viral strains
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

names.arg = colnames(n_mtx) #paste0(tolower(names(SID)[as.numeric(colnames(n_mtx))]), ", ", 1:k)
if (length(names.arg) != length(p)){
	names.arg = colnames(d_mtx)
}
em_output_df = cbind(p, names.arg)
em_output_df = as.data.frame(em_output_df)
colnames(em_output_df) = c("p", "names.arg")

write.csv(em_output_df, file = paste(syspath, "em_output_", case_num, ".csv", sep=""))

d <- data.table(em_output_df, key="p")
d <- tail(d, 50)

pdf(paste(syspath, "rplot_case_", case_num, ".pdf", sep="")) 
barplot(as.numeric(factor(d[["p"]])), beside = T, ylab = "p vec values", names.arg=d[["names.arg"]],
        cex.names=0.8, las=2, main='Viral strains estimation')
dev.off()

#grid()
cat("\nRun EM algorithm", proc.time() - ptm)
cat("\n")

