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
.write.tree3 <- function(phy, digits = 10, tree.prefix = "", check_tips = TRUE)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
    if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste0(":%.", digits, "g")
    n <- length(phy$tip.label)

    ## terminal branches:
    terms <- phy$edge[, 2] <= n
    TERMS <- phy$tip.label[phy$edge[terms, 2]]
    if (brl) TERMS <- paste0(TERMS, sprintf(f.d, phy$edge.length[terms]))

    ## internal branches, including root edge:
    INTS <- rep(")", phy$Nnode)
    if (nodelab) INTS <- paste0(INTS, phy$node.label)
    if (brl) {
        tmp <- phy$edge.length[!terms][order(phy$edge[!terms, 2])]
        tmp <- c("", sprintf(f.d, tmp))
        if (!is.null(phy$root.edge)) tmp[1L] <- sprintf(f.d, phy$root.edge)
        INTS <- paste0(INTS, tmp)
    }
    INTS[1] <- paste0(INTS[1], ";")

###    ## borrowed from phangorn:
###    parent <- phy$edge[, 1]
###    children <- phy$edge[, 2]
###    kids <- vector("list", n + phy$Nnode)
###    for (i in 1:length(parent))
###        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])
###    Nkids <- lengths(kids, FALSE)
###    root <- parent[! parent %in% children][1]
###
    o <- postorder(phy)
    ANC <- phy$edge[o, 1L]
    DESC <- phy$edge[o, 2L]
    NEWICK <- character(n + phy$Nnode)
    NEWICK[1:n] <- TERMS
    root <- n + 1L
    from <- to <- 1L
    repeat {
        thenode <- ANC[from]
        if (thenode == root) {
            to <- length(ANC)
        } else {
            while (ANC[to + 1L] == thenode) to <- to + 1L
        }
        tmp <- paste(NEWICK[DESC[from:to]], collapse = ",")
        tmp <- paste0("(", tmp, INTS[thenode - n])
        NEWICK[thenode] <- tmp
        if (thenode == root) break
        from <- to + 1L
    }
    NEWICK[root]
}
.write.tree3(tree,tree.prefix="global_out.tree")
