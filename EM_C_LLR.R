#!/usr/bin/env Rscript
# -h for help
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(hash))
suppressMessages(library(feather))
suppressMessages(library(glue))
suppressMessages(library(arrow))
suppressMessages(library(turboEM))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(CovidEM))
suppressMessages(library(getopt))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
suppressMessages(library(snow))

spec <- matrix(
  c("mismatch",  "i", 2, "character", "Name of the file containing the mismatch matrix. (This should be a TXT file.)",
    "error_rate", "e", 1, "double",  "Assumed error rate in the EM algorithm (default = 0.005).",
    "filter",  "f", 2, "double",  "The allele frequency cut-off we used to remove ‘unlikely’ strains.",
    "llr","l",0,"logical","Use this to perform the LLR procedure for each strain.",
    "num_show","n",1,"integer","Maximal number of strains to show in the plot (default = 10).",
    "help", "h", 0, "logical",  "Show help.",
    "site_LLR","s",0,"logical", "Perform the LLR procedure for each site.",
    "variant","v",1,"character" , "Name of the file containing variant sites.",
    "read_count","b",1,"character", "Name of the file containing allele counts from reads.",
    "ref_file","r",1,"character","MSA FASTA of SARS-CoV-2 reference strains.",
    "deletion_file","d",1,"character","Deletion proportion for each site.",
    "core_num","c",1,"integer","Number of cores to use"),
  byrow=TRUE, ncol=5)

opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$mismatch) || is.null(opt$filter)){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

if (!is.null(opt$site_LLR)){
  if (!file.exists(opt$variant)){
    print("Please provide a file containing the variant sites")
    quit()
  }
  if (!file.exists(opt$read_count)){
    print("Please provide a file containing the allele counts from reads")
    quit()
  }
  if (!file.exists(opt$ref_file)){
    print("Please provide a MSA FASTA of SARS-CoV-2 reference strains")
    quit()
  }
  if (!file.exists(opt$deletion_file)){
    print("Please provide a file containing deletion proportion for each site")
    quit()
  }
}

if (!file.exists(opt$mismatch)) {
  print("No such file exists")
  quit()
}


# Read in data
ptm <- proc.time()
data <- read.table(opt$mismatch, header = T, sep = "\t",quote = "",check.names = F)
cat("Loading matrices in from file", (proc.time() - ptm)[3])

if (nrow(data) < 20) {
  print("Not enough reads!")
  quit()
}
if (ncol(data)<3){
  r_file <- file(paste0("output_", opt$mismatch, ".txt"), "w")
  write("There is no strain left after filtering", file = r_file, append = T)
  close(r_file)
  quit()
}

opt$mismatch <- gsub(".txt", "", opt$mismatch, fixed = T)

ptm <- proc.time()
strains <- colnames(data)[3:ncol(data)]
names(strains) <- strains
k <- ncol(data) - 2 # number of viral strains

if (k==1){
  # par(mar = c(13, 4, 4, 1))
  names.arg <- strains
  em_output_df <- cbind(1, names.arg)
  em_output_df <- data.frame(em_output_df)
  colnames(em_output_df) <- c("p", "names.arg")
  
  write.csv(em_output_df, file = paste("./em_output_", opt$mismatch, ".csv", sep = ""))
  
  # d <- data.table(em_output_df, key="p")
  d <- data.frame(1, names.arg)
  colnames(d)[1] <- "p"
  
  pdf(file = paste("./proportion_plot_", opt$mismatch, ".pdf", sep = ""), width = 15, height = 10)
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
    xlim(c(0, 1))
  # barplot(as.numeric(d$p), beside = T, xlab = "proportion", names.arg=d[["names.arg"]],
  #         cex.names=0.8, las=2, main='Viral strains estimation',horiz = T,
  #         xlim = c(0,max(p)+0.1),col = "lightblue")
  dev.off()
} else{

# Mismatch matrix
d_mtx <- data[, 3:ncol(data)]
rownames(d_mtx) <- data[, 1]
strain_name_copy <- data[,1]
n_mtx <- matrix(rep(data[, 2], k), nrow = nrow(d_mtx), ncol = k)
rownames(n_mtx) <- data[, 1]
colnames(n_mtx) <- colnames(d_mtx)
cat("\nSet up for creating Q matrix", (proc.time() - ptm)[3])

# Calculating q_ij

if (is.null(opt$error_rate)) {
  error <- 0.005
} else {
  error <- opt$error_rate
  if (error < 0.001 | error > 1) {
    print("\nError rate should be smaller than 0.001 or larger than 1!")
    quit()
  }
}
ptm <- proc.time()
d_mtx <- as.matrix(d_mtx)
Q <- calcQ(as.matrix(n_mtx),as.matrix(d_mtx),error)
colnames(Q) <- colnames(d_mtx)
cat("\nCreate Q-matrix", (proc.time() - ptm)[3])
ptm <- proc.time()

#-------------------------------------SQUAREM---------------------
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
  r_file <- file(paste0("Unidentifiable_Strains_", opt$mismatch, ".txt"), "w")
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
Q[Q==0] <- .Machine$double.xmin
p0 <- colSums(Q)
# p0 <- runif(ncol(Q))
p0 <- p0 / sum(p0)
p1 <- runif(ncol(Q))
p1 <- p1 / sum(p1)
options(digits = 13)

# cl.cores = detectCores(logical = F)

ptm <- proc.time()
res <- turboem(
  par = p0, fixptfn = EM, objfn = neg_logL, method = "squarem", y = Q, parallel = F,
  control.run = list(
    convtype = "parameter", tol = 1.0e-7,
    stoptype = "maxtime", maxtime = 14400
  )
)
res2 <- turboem(
  par = p1, fixptfn = EM, objfn = neg_logL, method = "squarem", y = Q, parallel = F,
  control.run = list(
    convtype = "parameter", tol = 1.0e-7,
    stoptype = "maxtime", maxtime = 14400
  )
)
cat(paste0("\nDifference between two optimazations from different starting points ",
             signif(sum(abs(pars(res)-pars(res2)))),2))

p <- pars(res)
log.likelihood.with <- -neg_logL(p,Q)
p <- round.output(p)

# stopImplicitCluster()
cat("\nSQUAREM algorithm ", (proc.time() - ptm)[3],"\n")

# LLR
if (!is.null(opt$llr)){
  if(is.null(opt$core_num)){
    cl.cores <- 1
  } else {
    cl.cores <- opt$core_num
  }
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl,cores = cl.cores)
p_tmp <- pars(res)
# LLR <- rep(NA,length(p_tmp))
ptm <- proc.time()
# too_large_flag <- logical(length(p_tmp))
# for (ii in 1:length(p_tmp)){
#   if (p[ii]>=opt$filter){
#     p0 <- p_tmp[-ii]
#     p0 <- p0/sum(p0)
#     res.rm <- turboem(
#       par = p0, fixptfn = EM, objfn = neg_logL, method = "squarem", y = Q[,-ii], parallel = F,
#       control.run = list(
#         convtype = "objfn", tol = 1.0e-7,
#         stoptype = "maxtime", maxtime = 10000
#       )
#     )
#     log.likelihood.without <- -neg_logL(pars(res.rm),Q[,-ii])
#     if (is.na(log.likelihood.without)){
#       LLR[ii] <- 100
#       too_large_flag[ii] <- T
#     }else {
#       LLR[ii] <- 2*(log.likelihood.with-log.likelihood.without)
#     }
#   }
# }

tmp <- foreach (ii=1:length(p_tmp), .combine = cbind) %dopar% {
  require(turboEM)
  require(CovidEM)
  if (p[ii]>=opt$filter){
    p0 <- p_tmp[-ii]
    p0 <- p0/sum(p0)
    res.rm <- turboem(
      par = p0, fixptfn = EM, objfn = neg_logL, method = "squarem", y = Q[,-ii], parallel = F,
      control.run = list(
        convtype = "objfn", tol = 1.0e-7,
        stoptype = "maxtime", maxtime = 10000
      )
    )
    log.likelihood.without <- -neg_logL(pars(res.rm),Q[,-ii])
    if (is.na(log.likelihood.without)){
      # LLR[ii] <- 100
      # too_large_flag[ii] <- T
      return(c(100,T))
    }else {
      # LLR[ii] <- 2*(log.likelihood.with-log.likelihood.without)
      return(c(2*(log.likelihood.with-log.likelihood.without),F))
    }
  }
}
stopCluster(cl)
cat("\nLLR ", (proc.time() - ptm)[3],"\n")
LLR <- tmp[1,]
too_large_flag <- tmp[2,]
} else {
  LLR <- rep(NA,length(p))
  too_large_flag <- rep(F,length(p))
}


names.arg <- colnames(Q)
em_output_df <- cbind(names.arg,p,LLR,too_large_flag)
em_output_df <- data.frame(em_output_df,check.names = F)
colnames(em_output_df) <- c("names.arg","p", "LLR","flag")
em_output_df$p <- as.numeric(em_output_df$p)
em_output_df <- em_output_df[order(em_output_df$p,decreasing = T),]
em_output_df <- em_output_df[em_output_df$p>=opt$filter,]
em_output_df$p <- em_output_df$p/sum(em_output_df$p)
p_forsite <- em_output_df
em_output_df$p <- round(em_output_df$p,digits = 3)
em_output_df$LLR <- as.character(em_output_df$LLR)
em_output_df$flag <- as.logical(em_output_df$flag)
em_output_df$LLR[em_output_df$flag] <- ">>100"

if (!is.null(opt$llr)){
  write.csv(em_output_df[,1:3], file = paste("./em_output_", opt$mismatch, ".csv", sep = ""),row.names = F)
} else {
  write.csv(em_output_df[,1:2], file = paste("./em_output_", opt$mismatch, ".csv", sep = ""),row.names = F)
}


# d <- data.table(em_output_df, key="p")
d <- data.frame(p, LLR, names.arg,check.names = F)
d <- d[order(d$p, decreasing = T), ]
d <- d[d$p >= opt$filter, ]
d$p <- d$p/sum(d$p)
d$p <- round(d$p,digits = 3)
if (is.null(opt$num_show)){
  d <- head(d, 30)
} else{
  d <- head(d, opt$num_show)
}
d <- cbind(d,flag = factor((d$LLR>2)+(d$LLR>4),levels = c(0,1,2)))
levels(d$flag) <- c("LLR<=2","2<LLR<=4","LLR>4")

pdf(file = paste("./proportion_plot_", opt$mismatch, ".pdf", sep = ""), width = 15, height = 10)
if (!is.null(opt$llr)){
  prop_plot <- ggplot(data = d, mapping = aes(x = p, y = reorder(names.arg, p), fill = flag)) +
    geom_bar(stat = "identity") +
    labs(x = "Proportion", y = NULL) +
    geom_text(aes(label = p, vjust = 0.45, hjust = 0)) +
    theme(
      axis.text.y = element_text(size = 11),
      plot.title = element_text(
        size = 15,
        hjust = 0.4
      )
    ) +
    theme(legend.position = "bottom",legend.title=element_blank())+
    xlim(c(0, max(d$p) * 1.05))+ scale_fill_manual(values=c("LLR>4" = "orange","2<LLR<=4" = "lightblue","LLR<=2" = "grey"))
} else {
  prop_plot <- ggplot(data = d, mapping = aes(x = p, y = reorder(names.arg, p))) +
    geom_bar(stat = "identity",fill = "lightblue") +
    labs(x = "Proportion", y = NULL) +
    geom_text(aes(label = p, vjust = 0.45, hjust = 0)) +
    theme(
      axis.text.y = element_text(size = 11),
      plot.title = element_text(
        size = 15,
        hjust = 0.4
      )
    ) +
    theme(legend.position = "none",legend.title=element_blank())+
    xlim(c(0, max(d$p) * 1.05))
}

# barplot(as.numeric(d$p), beside = T, xlab = "proportion", names.arg=d[["names.arg"]],
#         cex.names=0.8, las=2, main='Viral strains estimation',horiz = T,
#         xlim = c(0,max(p)+0.1),col = "lightblue")
print(prop_plot)
dev.off()


##### site LLR
if (!is.null(opt$site_LLR)){
  # p <- read.csv(paste("./em_output_", opt$mismatch, ".csv", sep = ""))
  p <- p_forsite
  ref <- readLines(opt$ref_file)
  ref_name <- gsub(">","",ref[seq(1,length(ref),2)])
  ref <- ref[seq(2,length(ref),2)]
  names(ref) <- ref_name
  if (file.exists(paste0("Unidentifiable_Strains_", opt$mismatch, ".txt"))){
    un_group <- readLines(paste0("Unidentifiable_Strains_", opt$mismatch, ".txt"))
    if (length(un_group)>=2){
      group_tmp <- which(str_detect(un_group,"Group"))
      un_group_list <- vector("list",length(group_tmp))
      for (u in 1:length(group_tmp)){
        if (u!=length(group_tmp)){
          un_group_list[[u]] <- un_group[(group_tmp[u]+1):(group_tmp[u+1]-1)]
        } else {
          un_group_list[[u]] <- un_group[(group_tmp[u]+1):length(un_group)]
        }
      }
    }
  }
  
  
  ref_allele <- matrix(0,nrow = nchar(ref[1]),ncol = 4,
                       dimnames = list(1:nchar(ref[1]),
                                       c("A","C","G","T")))
  is_group <- str_detect(p$names.arg,"Group")
  update_ref_allele <- function(x,t_p){
    tmp <- c("A"=0,"C"=0,"G"=0,"T"=0)
    tmp[x] <- t_p
    tmp
  }
  for (t in p$names.arg[!is_group]){
    t_seq <- str_split(ref[t],pattern = "")[[1]]
    t_p <- p[p$names.arg==t,2]
    t_tmp <- t(apply(rbind(t_seq),2,update_ref_allele,t_p=t_p))
    ref_allele <- ref_allele+t_tmp
  }
  # save(ref_allele,file = "ref_allele.RData")
  
  un_p <- rep(0,nchar(ref[1]))
  if(sum(is_group)!=0){
    for (u in p$names.arg[is_group]){
      p_u <- p[p$names.arg==u,2]
      u <- substr(u,str_locate(u,"\\d")[1,1],nchar(u))
      u <- as.numeric(u)
      un_allele <- matrix(0,nrow = nchar(ref[1]),ncol = 4,
                          dimnames = list(1:nchar(ref[1]),
                                          c("A","C","G","T")))
      for (t in un_group_list[[u]]){
        t_seq <- str_split(ref[t],pattern = "")[[1]]
        t_tmp <- t(apply(rbind(t_seq),2,update_ref_allele,t_p=1))
        un_allele <- un_allele+t_tmp
      }
      for (l in 1:length(t_seq)){
        ref_a <- colnames(ref_allele)[which.max(ref_allele[l,])]
        if(sum(un_allele[l,]!=0)==1){
          ref_allele[l,un_allele[l,]!=0] <- ref_allele[l,un_allele[l,]!=0]+p_u
        } else {
          un_p[l] <- un_p[l]+p_u
        }
      }
    }
  }
  
  
  
  read_allele <- read.table(opt$read_count,header = T,
                            sep = "\t",row.names = 1)
  site_LLR <- rep(NA,length(t_seq))
  variant_site <- read.table(opt$variant)[-1,1]+1
  for (l in variant_site){
    if (sum(read_allele[l,])!=0){
      ref_a <- colnames(ref_allele)[which.max(ref_allele[l,])]
      nr <- read_allele[l,ref_a]
      nm <- sum(read_allele[l,])-nr
      f_hat <- min(ref_allele[l,ref_a],1)
      f_mle <- (nr*(1-opt$error_rate/3)-nm*opt$error_rate/3)/(nr+nm)/(1-4*opt$error_rate/3)
      f_mle <- max(min(f_mle,1),0)
      if (nm==0){
        p_hat <- f_hat*(1-opt$error_rate)+(1-f_hat)*opt$error_rate/3
        p_mle <- f_mle*(1-opt$error_rate)+(1-f_mle)*opt$error_rate/3
        p_hat2 <- (f_hat+un_p[l])*(1-opt$error_rate)+(1-f_hat-un_p[l])*opt$error_rate/3
        site_LLR[l] <- min(nr * log(p_hat)-nr * log(p_mle),
                           nr * log(p_hat2)-nr * log(p_mle))
      } else {
        # p_tmp <- sum(p[is_group,2])
        p_hat <- f_hat*(1-opt$error_rate)+(1-f_hat)*opt$error_rate/3
        p_mle <- f_mle*(1-opt$error_rate)+(1-f_mle)*opt$error_rate/3
        p_hat2 <- (f_hat+un_p[l])*(1-opt$error_rate)+(1-f_hat-un_p[l])*opt$error_rate/3
        site_LLR[l] <- min(nr * log(p_hat)+nm * log(1-p_hat)
                           -nr * log(p_mle)-nm * log(1-p_mle),
                           nr * log(p_hat2)+nm * log(1-p_hat2)
                           -nr * log(p_mle)-nm * log(1-p_mle))
      }
    }
  }
  # save(site_LLR,file = paste0("./site_LLR_",opt$mismatch,".RData"))
  deletion <- read.table(opt$deletion_file,header = T,sep = "\t")
  site_LLR[deletion[,1]+1] <- NA
  site_frame <- data.frame(index = 1:length(site_LLR),
                           LLR = abs(site_LLR))
  site_frame <- filter(site_frame,!is.na(LLR))
  site_frame <- site_frame[order(site_frame$LLR,decreasing = T),]
  write.csv(site_frame,file = paste0(opt$mismatch,"_siteLLR.csv"))
}

}