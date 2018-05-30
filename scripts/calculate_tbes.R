library(parallel)

datadir = "/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulate_data/datasets/"
booster = "/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulate_data/booster_macos64"

dirs = list.dirs(path = datadir)[-c(1)]

calculate_tbe = function(dir){
    
    command = paste("",
                    booster,
                    " -i ",
                    file.path(dir, "alignment.fasta.treefile"),
                    " -b ", 
                    file.path(dir, "alignment.fasta.boottrees"),
                    " -o ",
                    file.path(dir, "alignment.fasta.tbetree"),
                    sep = ""
                    )
    
    command2 = paste("",
                    booster,
                    " -i ",
                    file.path(dir, "alignment.fasta.treefile"),
                    " -b ", 
                    file.path(dir, "alignment.fasta.boottrees"),
                    " -o ",
                    file.path(dir, "alignment.fasta.fbptree"),
                    " -a fbp ",
                    sep = ""
    )
    
    # don't redo ones we've already done
    if(file.exists(file.path(dir, "alignment.fasta.tbetree"))==FALSE){
        system(command)
    }
    if(file.exists(file.path(dir, "alignment.fasta.fbptree"))==FALSE){
        system(command2)
    }
}

mclapply(dirs, calculate_tbe, mc.cores = 4)
