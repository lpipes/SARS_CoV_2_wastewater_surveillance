#!/usr/bin/env Rscript
# Commands [mismatch matrix](.txt) [error rate](default = 0.005)
# library(turboEM)
library(Matrix)
library(plyr)
library(data.table)
library(dplyr)
library(stringr)
library(hash)
library(feather)
library(glue)
library(arrow)
library(daarem)
library(ggplot2)
library(stringr)


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
data <- read.table(args[1], header = T)
cat("Loading matrices in from file", (proc.time() - ptm)[3])

if (nrow(data)<20){
  print("Not enough reads!")
  quit()
}

ptm <- proc.time()
strains <- colnames(data[, 3:ncol(data)])
names(strains) <- strains
k <- ncol(data) - 2 # number of viral strains

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
}
ptm <- proc.time()
Q <- ((error)^d_mtx) * ((1 - error)^(n_mtx - d_mtx))
cat("\nCreate Q-matrix", (proc.time() - ptm)[3])
ptm <- proc.time()

#-------------------------------------DAAREM---------------------
# objective function
logL <- function(p, y) {
  return(sum(log(colSums(t(y) * p))))
}
# #objective function for SQUAREM
# neg.logL <- function(p, y) return((-logL(p,y)))

# Update Function
EM <- function(p, y) {
  w <- t(t(y) * p) # numerator of Eq.12 in Malone's paper
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

# EM Algorithm
p0 <- runif(k)
p0 <- p0 / sum(p0)
res <- daarem(p0,
  fixptfn = EM, objfn = logL, y = Q,
  control = list(tol = 1e-07, convtype = "objfn", maxiter = 200000)
)
p <- res$par
p <- round.output(p)
cat("\nDAAREM algorithm ", (proc.time() - ptm)[3])

# One Strain
loglikelihoods <- apply(diag(k), 1, logL, y = Q)
cat("\nAssuming only one strain: ", colnames(n_mtx)[which.max(loglikelihoods)])
p.one <- rep(0, k)
p.one[which.max(loglikelihoods)] <- 1
cat("\nlog(Likelihood of multiple strains/Likelihood of one strain): ", logL(p, Q) - logL(p.one, Q),"\n")

par(mar = c(13, 4, 4, 1))
names.arg <- colnames(n_mtx) # paste0(tolower(names(SID)[as.numeric(colnames(n_mtx))]), ", ", 1:k)
if (length(names.arg) != length(p)) {
  names.arg <- colnames(d_mtx)
}
em_output_df <- cbind(p, names.arg)
em_output_df <- as.data.frame(em_output_df)
colnames(em_output_df) <- c("p", "names.arg")

write.csv(em_output_df, file = paste("./em_output_", args[1], ".csv", sep = ""))
reads_caset10_30000_423_.005.txt

# d <- data.table(em_output_df, key="p")
d <- data.frame(p, names.arg)
d <- d[order(d$p, decreasing = T), ]
d <- head(d, 50)
d <- d[d$p != 0, ]

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
