library(NELSI)
library(TreeSim)
library(phangorn)
library(phytools)
library(plyr)

random.paired.path.dist = function(N){
  
  a = rtree(N, rooted = FALSE, tip.label = NULL, br = NULL)
  b = rtree(N, rooted = FALSE, tip.label = NULL, br = NULL)

  return(path.dist(a, b))
}

random.oneway.path.dist = function(N, t){
  
  a = rtree(N, rooted = FALSE, tip.label = NULL, br = NULL)
  return(path.dist(a, t))
  
}

normalised.oneway.path.dist = function(t1, t2){
  
  observed = path.dist(t1, t2)
  
  N = length(t1$tip.label)
  
  mean.pd = mean(replicate(n = 1000, expr = random.oneway.path.dist(N)))
  
  norm.pd = observed / mean.pd
  
  return(norm.pd)    
}

normalised.path.dist = function(t1, t2){
  
    observed = path.dist(t1, t2)
    
    N = length(t1$tip.label)
    
    mean.pd = mean(replicate(n = 1000, expr = random.paired.path.dist(N)))
  
    norm.pd = observed / mean.pd

    return(norm.pd)    
}



analyse_data <- function(id = 1, output_folder, iqtree_path, n_taxa = 10, n_sites = 100, birth_rate = 0.5, death_rate = 0.0, mrate_sd = 0.1, tree_age = 1, mol_rate = 0.1){

  print(paste("analysing", id, ": n_taxa", n_taxa, "; n_sites", n_sites))

  analysis_folder = file.path(output_folder, id)

  # run if not already done
  logfile = file.path(analysis_folder, "alignment.fasta.log")
  if(file.exists(logfile)){
    done = grep("Analysis results written to: ", readLines(logfile))
    if(length(done)==0){
      dir.create(analysis_folder)
      dataset = get_data(n_taxa, n_sites, birth_rate, death_rate, mrate_sd, tree_age, mol_rate)
      write_data(analysis_folder, dataset)
      run_iqtree(iqtree_path, analysis_folder)
    }
  }else{
    dir.create(analysis_folder)
    dataset = get_data(n_taxa, n_sites, birth_rate, death_rate, mrate_sd, tree_age, mol_rate)
    write_data(analysis_folder, dataset)
    run_iqtree(iqtree_path, analysis_folder)			
  }

}
    

# simulate a data set that includes the chonogram, phylogram, sequence data, Q matrix, tip states, and number of transitions
get_data <- function(n_taxa = 10, n_sites = 100, birth_rate = 0.5, death_rate = 0.0, mrate_sd = 0.1, tree_age = 1, mol_rate = 0.1){

  # Simulate chronogram
  tr_sim <- sim.bd.taxa.age(n = n_taxa, numbsim = 1, lambda = birth_rate, mu = death_rate, frac = 1, age = tree_age, mrca = FALSE)[[1]]
  tr_sim$edge.length <- tr_sim$edge.length * (tree_age / max(branching.times(tr_sim)))

  # Simulate phylogram
  phylo_sim <- tr_sim
  phylo_sim$edge.length <- tr_sim$edge.length * rlnorm(length(tr_sim$edge.length), meanlog = log(mol_rate), sdlog = mrate_sd)

  # Simulate sequence data
  if(n_sites > 0){
    dna_sim <- as.DNAbin(simSeq(phylo_sim, l = n_sites))
  }else{
    dna_sim <- matrix('n', ncol = 1, nrow = n_taxa)
    rownames(dna_sim) <- tr_sim$tip.label
    dna_sim <- as.DNAbin(dna_sim)
  }
    
  return(list(chronogram = tr_sim, phylogram = phylo_sim, seq_data = dna_sim))
}

write_data <- function(analysis_folder, dataset){

  write.dna(dataset$seq_data, file = file.path(analysis_folder, "alignment.fasta"), format = "fasta")
  write.tree(dataset$chronogram, file = file.path(analysis_folder, "chronogram.phy"))
  write.tree(unroot(dataset$phylogram), file = file.path(analysis_folder, "phylogram.phy"))

}

run_iqtree <- function(iqtree_path, analysis_folder){

  alignment = file.path(analysis_folder, "alignment.fasta")
    phylogram = file.path(analysis_folder, "phylogram.phy")
      command = paste(iqtree_path, "-s", alignment, "-m JC", "-redo -quiet -b 100")
        system(command)

}


