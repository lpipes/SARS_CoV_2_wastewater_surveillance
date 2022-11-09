library(ape)
args = commandArgs(trailingOnly=TRUE)
if ( length(args)==0 ){
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
tree<-read.tree(args[1])
leaves_to_prune<-read.table(args[2],header=F)
tree<-drop.tip(tree,leaves_to_prune$V1,trim.internal=T)
tree<-multi2di(tree)
tree<-collapse.singles(tree)
write.tree(tree,file="global_out.tree")
