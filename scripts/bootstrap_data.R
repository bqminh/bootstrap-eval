library(distory)
library(reshape2)
library(dplyr)
library(ggplot2)
library(parallel)
datadir = "/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulate_data/datasets/"
dirs = list.dirs(path = datadir)[-c(1)]


one_tree = function(dir){
    
    fbp = read.tree(file.path(dir, "alignment.fasta.fbptree"))
    tbe = read.tree(file.path(dir, "alignment.fasta.tbetree"))
    true = read.tree(file.path(dir, "phylogram.phy"))

    # classify ML tree branches into true.false
    
    # FALSE edges (present in FBP tree but not the true tree)
    # and warning: this function is super slow with 500 taxa
    d = distinct.edges(fbp, true)

    # the node labels only apply to internal edges. So to get the index of the false edges in teh node.labl, 
    # we get the child node from edge.label, then subtract the number of taxa in the tree
    ntax = length(fbp$tip.label)
    false_nodes = fbp$edge[,2][d] - ntax
        
    fbps = fbp$node.label
    tbes = tbe$node.label
    
    r = data.frame(FBP = fbps, TBE = tbes, recovered = TRUE)
    
    # set the false nodes to false
    r$recovered[false_nodes] = FALSE

    # trim off the first row, because ape adds it on for the root node
    r = r[-1,]
    
    # add extra info
    r$taxa = ntax
    aln = read.dna(file.path(dir, "alignment.fasta"), format = "fasta", as.character = T)
    seqlen = length(aln[1,])
    r$sites = seqlen
    
    return(r)
    
}

branches = do.call("rbind", mclapply(dirs, one_tree, mc.cores = 4))
branches$id = 1:nrow(branches)

branches$TBE = as.numeric(as.character(branches$TBE))
branches$FBP = as.numeric(as.character(branches$FBP))

#branches$diff = branches$TBE - branches$FBP

b = melt(branches, id.vars = c("recovered", "taxa", "id", "sites"))
b$bootstrap = round(as.numeric(b$value)*100)

write.csv(b, "/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulate_data/b.csv")
