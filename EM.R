#!/usr/bin/env Rscript
# Commands [mismatch matrix](.txt) [error rate](default = 0.005)
library(Matrix)
library(plyr)
library(data.table)
library(dplyr)
library(stringr)
library(hash)
library(feather)
library(glue)
library(arrow)
library(turboEM)
library(ggplot2)
library(stringr)
# library(parallel)
# library(doParallel)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  print("Not enough arguments!")
  quit()
}
if (!file.exists(args[1])) {
  print("No such file exists")
  quit()
}

# Read in data
ptm <- proc.time()
data <- read.table(args[1], header = T, sep = "\t")
cat("Loading matrices in from file", (proc.time() - ptm)[3])

if (nrow(data) < 20) {
  print("Not enough reads!")
  quit()
}

args[1] <- gsub(".txt", "", args[1], fixed = T)

ptm <- proc.time()
strains <- colnames(data)[3:ncol(data)]
names(strains) <- strains
k <- ncol(data) - 2 # number of viral strains

if (k==1){
  par(mar = c(13, 4, 4, 1))
  names.arg <- strains
  em_output_df <- cbind(1, names.arg)
  em_output_df <- data.frame(em_output_df)
  colnames(em_output_df) <- c("p", "names.arg")
  
  write.csv(em_output_df, file = paste("./em_output_", args[1], ".csv", sep = ""))
  
  # d <- data.table(em_output_df, key="p")
  d <- data.frame(1, names.arg)
  colnames(d)[1] <- "p"
  
  pdf(file = paste("./proportion_plot_", args[1], ".pdf", sep = ""), width = 15, height = 10)
  ggplot(data = d, mapping = aes(x = p, y = reorder(names.arg, p))) +
    geom_bar(stat = "identity", fill = "lightblue") +
    labs(x = "Proportion", y = NULL) +
    geom_text(aes(label = p, vjust = 0.45, hjust = 0)) +
    theme(
      axis.text.y = element_text(size = 11),
      plot.title = element_text(
        size = 15,
        hjust = 0.4
      )
    ) +
    labs(title = "Viral Strains Estimation") +
    xlim(c(0, 1))
  # barplot(as.numeric(d$p), beside = T, xlab = "proportion", names.arg=d[["names.arg"]],
  #         cex.names=0.8, las=2, main='Viral strains estimation',horiz = T,
  #         xlim = c(0,max(p)+0.1),col = "lightblue")
  dev.off()
} else{

# Mismatch matrix
d_mtx <- data[, 3:ncol(data)]
rownames(d_mtx) <- data[, 1]
n_mtx <- matrix(rep(data[, 2], k), nrow = nrow(d_mtx), ncol = k)
rownames(n_mtx) <- data[, 1]
colnames(n_mtx) <- colnames(d_mtx)
cat("\nSet up for creating Q matrix", (proc.time() - ptm)[3])

# Calculating q_ij

if (length(args) == 1) {
  error <- 0.005
} else if (is.na(as.numeric(args[2])) == T) {
  print("Wrong format for error rate!")
  quit()
} else {
  error <- as.numeric(args[2])
  if (error < 0.001 | error > 1) {
    print("\nError rate should be smaller than 0.01 or larger than 1!")
    quit()
  }
}
ptm <- proc.time()
Q <- ((error)^d_mtx) * ((1 - error)^(n_mtx - d_mtx))
Q <- as.matrix(Q)
cat("\nCreate Q-matrix", (proc.time() - ptm)[3])
ptm <- proc.time()

#-------------------------------------SQUAREM---------------------
# objective function
logL <- function(p, y) {
  return(sum(log(colSums(t(y+1e-20) * p))))
}
# objective function for SQUAREM
neg.logL <- function(p, y) {
  return((-logL(p, y)))
}

# Update Function
EM <- function(p, y) {
  w <- t(t(y+1e-20) * p) # numerator of Eq.12 in Malone's paper
  E_step <- w / rowSums(w) # denominator of Eq.12 in Malone's paper (normalization of rows)
  # return(colSums(E_step)/nrow(y))# M-step to compute p, Eq.33
  return(colMeans(E_step))
}

round.output <- function(p) {
  p[p < 0] <- 0
  p <- round(p, digits = 3)
  p <- as.matrix(p)
  p <- as.vector(p)
  p
}

# Combining unidentifiable strains
ptm <- proc.time()
loglikelihoods <- colSums(log(Q))
p.one <- rep(0, k)
p.one[which.max(loglikelihoods)] <- 1
uni.like <- unique(loglikelihoods)
if (length(uni.like) != length(loglikelihoods)) {
  r_file <- file(paste0("Unidentifiable_Strains_", args[1], ".txt"), "w")
  unident <- 0
  for (i in uni.like) {
    if (sum(loglikelihoods == i) > 1) {
      tmp <- which(loglikelihoods == i)
      # The first two column of d_mtx aren't mismatch number
      d_tmp <- unique(as.matrix(d_mtx[, tmp]), MARGIN = 2)
      if (ncol(d_tmp) < length(tmp)) {
        for (j in 1:ncol(d_tmp)) {
          unident.names <- colnames(d_mtx[, tmp])[apply(d_mtx[, tmp], 2, identical, y = d_tmp[, j])]
          if (length(unident.names) > 1) {
            unident <- unident + 1
            write(paste0("Group ", unident, ":"), file = r_file, append = T)
            write(unident.names, file = r_file, append = T)
            Q.index <- which(colnames(Q) %in% unident.names)
            colnames(Q)[Q.index[1]] <- paste0("Group ", unident)
            Q <- as.matrix(Q[, -Q.index[-1]])
          }
        }
      }
    }
  }
  close(r_file)
  if (ncol(Q) == 1) {
    colnames(Q) <- "Group 1"
  }
}
cat("\nTime cost for finding unidentifiable strains: ", (proc.time() - ptm)[3])

# EM Algorithm
p0 <- colSums(Q)
p0 <- p0 / sum(p0)
options(digits = 13)

# cl.cores = detectCores(logical = F)
# cl <- makeCluster(cl.cores)
# registerDoParallel(cl,cores = cl.cores)

res <- turboem(
  par = p0, fixptfn = EM, objfn = neg.logL, method = "squarem", y = Q, parallel = F,
  control.run = list(
    convtype = "objfn", tol = 1.0e-7,
    stoptype = "maxtime", maxtime = 10000
  )
)

p <- pars(res)
p <- round.output(p)

# stopImplicitCluster()
cat("\nSQUAREM algorithm ", (proc.time() - ptm)[3])

# One Strain
loglikelihoods <- colSums(log(Q))
cat("\nAssuming only one strain: ", colnames(n_mtx)[loglikelihoods == max(loglikelihoods)], sep = "\n")
p.one <- rep(0, ncol(Q))
p.one[which.max(loglikelihoods)] <- 1
cat("\nlog(Likelihood of multiple strains/Likelihood of one strain): ", logL(p, Q) - logL(p.one, Q), "\n")

par(mar = c(13, 4, 4, 1))
names.arg <- colnames(Q)
em_output_df <- cbind(p, names.arg)
em_output_df <- data.frame(em_output_df)
colnames(em_output_df) <- c("p", "names.arg")

write.csv(em_output_df, file = paste("./em_output_", args[1], ".csv", sep = ""))

# d <- data.table(em_output_df, key="p")
d <- data.frame(p, names.arg)
d <- d[order(d$p, decreasing = T), ]
d <- d[d$p > 0.01, ]
d <- head(d, 30)

pdf(file = paste("./proportion_plot_", args[1], ".pdf", sep = ""), width = 15, height = 10)
ggplot(data = d, mapping = aes(x = p, y = reorder(names.arg, p))) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(x = "Proportion", y = NULL) +
  geom_text(aes(label = p, vjust = 0.45, hjust = 0)) +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(
      size = 15,
      hjust = 0.4
    )
  ) +
  labs(title = "Viral Strains Estimation") +
  xlim(c(0, max(p) * 1.05))
# barplot(as.numeric(d$p), beside = T, xlab = "proportion", names.arg=d[["names.arg"]],
#         cex.names=0.8, las=2, main='Viral strains estimation',horiz = T,
#         xlim = c(0,max(p)+0.1),col = "lightblue")
dev.off()
}