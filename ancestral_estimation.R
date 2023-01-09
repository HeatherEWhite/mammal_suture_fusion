# Title: Ancestral state estimation

# Name: Heather White

# Date created: 02/12/21

# Last modified: 25/10/22

# License: MIT license


# Estimation of the ancestral state using ML for total suture closure scores

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ape)
library(nlme)
library(ggrepel)
library(broom)
library(phytools)


#########################################################################################################

# STEP 1: Load the data

# Import trees in Nexus format - branch lengths needed
tree <- "Data/my_mammal_tree.nexus"  
# Read the tree for analysis
tree_species <- read.nexus(tree)
Phylogeny_species_list <- read.csv("Data/tree_taxa_names.csv")

# Plot tree
plotTree(tree_species,ftype="i")


# Load total suture closure score data - adults only
TSCS_adults <- read.csv("Data/total_suture_closure_scores_adults.csv")
# Select the closure scores only
TSCS_adults <- TSCS_adults$closure_score_.
# Rename the column heading
names(TSCS_adults)[names(TSCS_adults) == "closure_score_."] <- "total_closure_score"
# Associate phylogeny names with the scores
names(TSCS_adults) <- Phylogeny_species_list$Taxa_names


#########################################################################################################

# STEP 2: Perform ancestral state estimation

# Ancestral estimation of total suture closure score (%) for adults only
# TSCS calculated by combining all suture scores for one specimen 


# Calculate ancestral states by using ML
anc_adult_suture_closure_score <- anc.ML(tree_species, TSCS_adults, CI = F)

# Plot the phylogeny with node labels to check where each is
plot(tree_species, cex=0.5, show.node.label = T)
nodelabels(cex=0.5) # add ancestral node information to the tree

# Combine the ancestral states with the original MSC_species results
combined <- data.frame(adult_suture_closure_score = c(TSCS_adults, anc_adult_suture_closure_score$ace))

# Name the ancestral clade nodes
rownames(combined)[rownames(combined) == "23"] <- "Ancestral mammal"
rownames(combined)[rownames(combined) == "24"] <- "Ancestral therian mammal"
rownames(combined)[rownames(combined) == "25"] <- "Ancestral placental mammal"
rownames(combined)[rownames(combined) == "27"] <- "Ancestral Laurasiatheria"
rownames(combined)[rownames(combined) == "31"] <- "Ancestral Euarchontoglires"
rownames(combined)[rownames(combined) == "36"] <- "Ancestral Afrotheria"
rownames(combined)[rownames(combined) == "37"] <- "Ancestral Xenarthra"
rownames(combined)[rownames(combined) == "39"] <- "Ancestral Marsupialia"

# Save only the ancestral states I am interested in (the ones named above)
combined <- combined %>%
  slice(c(1:25,27,31,36,37,39))
combined

# Save data if necessary
#write.csv(combined, file = "Analysis/Suture_closure/reconstructed_ancestral_states_suture_closure_adults.csv")


#########################################################################################################

# STEP 3: Map the traits onto the phylogeny

# Map the traits onto the tree and plot
obj <- contMap(tree_species, TSCS_adults, plot=FALSE)
plot(obj, legend = 0.7 * max(nodeHeights(tree_species)), fsize=c(0.7,0.9))


### END ###
