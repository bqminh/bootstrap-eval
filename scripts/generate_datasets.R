#NOTE: change the path to functions.R script of the git repo

source('/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulate_data/functions.R')

library(lhs)
library(parallel)

# parameters
# id: integer number of simulation
# output_folder: base folder for putting new datasets 
# iqtree_path: path to iq_tree executable 
# n_taxa: number of taxa in tree: 20-1000 
# n_sites: number of sites in aln: 20-1000 
# birth_rate: birth rate for BD tree: 0.1
# death_rate: death rate for BD tree: 0.01-0.99
# mrate_sd: standard deviation of molecular rate: 0-1
# tree_age: tree age: always fixed to 1 
# mol_rate: molecular rate in subst. per million years: 0.01-1

# first we generate data on all the points we need
num_samples = 3200

# generate our simulation parameters using LHS
raw_lhs = as.data.frame(improvedLHS(num_samples,3))
samples = raw_lhs
names(samples) = c("death_rate", "mrate_sd", "mol_rate")
# now we fix them up to be the right ranges, because they all start as 0-1 continuous
samples$mol_rate = (samples$mol_rate*0.99)+0.01
samples$death_rate = (samples$death_rate*0.98)+0.01
samples$id = 1:num_samples

# reps - fewer reps for larger trees because they have more branches of
samples$n_taxa = c(rep(20, 2500), rep(100,500), rep(500,100))
samples$n_sites  = 100

# sanity checks
head(samples)
summary(samples)

#NOTE: change the path to the output samples.csv file

write.csv(samples, file = "/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulate_data/samples.csv", row.names = FALSE)

samples = read.csv("/home/rob/bootstrap/samples.csv", header = TRUE)

output_folder = "/home/rob/bootstrap/datasets"
iqtree_path = "/home/rob/bootstrap/iqtree"

analyse_one_sample <- function(rownum, samples){

  row = samples[rownum,]

  analyse_data(id = row$id, 
                output_folder = output_folder, 
                iqtree_path = iqtree_path, 
                n_taxa = row$n_taxa, 
                n_sites = row$n_sites, 
                birth_rate = 0.1, 
                death_rate = row$death_rate, 
                mrate_sd = row$mrate_sd, 
                tree_age = 1, 
                mol_rate = row$mol_rate)

}

# do the harder ones first, to optimise CPU packing
samples <- samples[order(-samples$n_taxa, -samples$n_sites),]

# generate the datasets
mclapply(1:num_samples, analyse_one_sample, samples = samples, mc.cores = 30)
