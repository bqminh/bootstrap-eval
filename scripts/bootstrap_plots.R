library(distory)
library(reshape2)
library(dplyr)
library(ggplot2)
library(parallel)

s = read.csv("/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/simulated.csv")
p = read.csv("/Users/roblanfear/Dropbox/Projects_Current/bootstrap_reply/pandit.csv")

s$category = paste(s$taxa, "taxon birth death trees with 100 simulated sites")
p$category = "PANDIT database simulations (4-403 taxa, 24-6891 sites)"

# do the analysis for the simulated data
b = rbind(s, p)

# now for each bootstrap proportion, let's get the percentage of branches with that proportion that were actually true.
true = subset(b, recovered == TRUE)
c = b %>% group_by(category, variable, bootstrap) %>% summarise(n = n())
d = true %>% group_by(category, variable, bootstrap) %>% summarise(n = n())
names(d) = c("category", "variable", "bootstrap", "true")
e = merge(c, d, by = c("category", "variable", "bootstrap"), all = T)
e$true[which(is.na(e$true))] = 0
e$proportion = e$true*100 / e$n

e$category <- factor(e$category, levels = c("20 taxon birth death trees with 100 simulated sites", "100 taxon birth death trees with 100 simulated sites", "500 taxon birth death trees with 100 simulated sites", "PANDIT database simulations (4-403 taxa, 24-6891 sites)"))


ggplot(e, aes(x = bootstrap, y = proportion)) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    annotate("rect", xmin=70, xmax=72, ymin=0, ymax=100, alpha = 0.2, fill = 'red') +
    annotate("rect", xmin=90, xmax=92, ymin=0, ymax=100, alpha = 0.2, fill = 'blue') +
    geom_hline(yintercept = 95, linetype = 'dotted', alpha = 1) +
    geom_line(aes(colour = variable), size = 1.5, alpha = 0.75) + 
    facet_wrap(~category) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    labs(colour="bootstrap") +
    xlab("Bootstrap percentage for FBP and TBE") +
    ylab("Proportion of branches that are correct")

    