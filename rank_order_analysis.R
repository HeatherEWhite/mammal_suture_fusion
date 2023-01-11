# Title: Rank order analysis

# Name: Heather White

# Date created: 02/12/21

# Last modified: 25/10/22

# License: MIT license


# Rank order analysis - between species and compared with the Krogman order (Krogman, 1930)

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


#########################################################################################################

# STEP 1: Load and format the data


# Order the sutures based on the percentage of fusion (suture closure scores) calculated (as described in manuscript)
# Suture closure scores calculated in "suture_closure_scores.R" script
# Load this rank order data
rank <- read.csv("Data/rank_order_sutures.csv")

# Load the adult specimen info
adults_data <- read.csv("Data/suture_closure_status_adults.csv")

# Reorder and rename headings so they match the species names
rownames(rank) <- rank$X
rank <- rank %>% select(2:23)
colnames(rank) <- adults_data$Species


#########################################################################################################

# STEP 2: Species pairwise comparison in rank order closure - Kendall's tau


# Iteratively create two species pairs and pairwise Kendalls test on rank order of suture closure and save output data from comparison
PairwiseComparisons=list()
i=1
n=2
for (i in 1:22){ 
  for (h in n:22){ 
    if (i==length(adults_data$Species)){
      break
    }else
      SpeciesA <- adults_data$Species[i]
    SpeciesB <- adults_data$Species[h]
    toMatch <- c(SpeciesA, SpeciesB)
    filename <- paste(SpeciesA, SpeciesB, sep=" vs. ")
    cor_pairwise <-cor.test(rank[[SpeciesA]], rank[[SpeciesB]], method = "kendall", conf.level = .95)
    PairwiseComparisons[[filename]]<-cor_pairwise
  }
  n=n+1
}

# Create matrix and fill with values from pairwise comparisons
Pvalues = matrix (nrow = length(names(PairwiseComparisons)), 
                  ncol = length(c(paste("Kendall", names(PairwiseComparisons[[i]]$p.value[1]), sep = "_"), names(PairwiseComparisons[[i]]$estimate[1]))),
                  dimnames = list(c(names(PairwiseComparisons)), c(paste("Kendall", names(PairwiseComparisons[[i]]$p.value[2]), sep = "_"), names(PairwiseComparisons[[i]]$estimate[1])))
)

# Fills the matrix with all pairwise comparisons from PairwiseComparisons object using $aov.table
for (i in 1:length(names(PairwiseComparisons))){
  Pvalues[i,] <- rbind(as.numeric(paste(c(PairwiseComparisons[[i]]$p.value[1], PairwiseComparisons[[i]]$estimate[1]))))
}

# Save if necessary
#write.csv(Pvalues, "Pairwise_comparison_results_rank_order_suture_closure.csv")

########

# Perform Bonferroni correction to account for pairwise comparisons
# Create the empty matrix
CorrectedPvalues <- matrix(nrow = length(names(PairwiseComparisons)),
                           ncol = length(c("corrected Kendall p.value")),
                           dimnames = list(c(names(PairwiseComparisons)), c("corrected_Kendall_p.value"))
)

# Fill the matrix with the corrected p-values
CorrectedPvalues[,1] <- p.adjust(Pvalues[,1], method = "bonferroni", n = length(PairwiseComparisons))

# Save if necessary
#write.csv(CorrectedPvalues, "Pairwise_comparisons_corrected_pvalues_rank_order_suture_closure.csv")


#########################################################################################################

# STEP 3: Average rank order closure vs Krogman pattern - Kendall's tau


# Load the data
# Average across sutures in each Krogman region (adult specimens only) 
# Then averaged across all species, to calculate average closure and average rank order
average_adults <- read.csv("Data/average_adult_score_Krogman_rank.csv")
# Same as above but seperated for marsupials and placentals
marsupial_placental_pattern <- read.csv("Data/average_total_suture_closure_adults_marsupials_vs_placentals_Krogman.csv")


# Average total suture closure score for each Krogman region calculated across the adult only species
# This average score is ranked from 1-6 for the Krogman pattern 
# and compared using Kendall's tau to the Krogman pattern
average_adult_vs_Krogman_rank <-cor.test(average_adults$Adult_rank, average_adults$Krogman.rank, 
                                         method = "kendall", exact = F)
average_adult_vs_Krogman_rank


# Average total suture closure score is calculated across the adult only species for each Krogman region
# Divided into the marsupials as one group and the placentals as another
# This average score is ranked from 1-6 for the Krogman pattern 
# and compared using Kendall's tau to the Krogman pattern
# For both marsupials and placentals seperately
marsupial_adult_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Marsupials_adults_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)
marsupial_adult_vs_Krogman_rank
placental_adult_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Placentals_adults_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)
placental_adult_vs_Krogman_rank

# A significant p-value (p < 0.05) indicates a significant correlation with the Krogman pattern


#########################################################################################################

# STEP 4: Plot the order of suture fusion across all 31 sutures


# Load the data - order of suture fusion
# A .csv containing the order of fusion for each suture and for each species - with all specimens combined for each species
# Data below is reformatted from the data loaded in STEP 1
rank_box <- read.csv("Data/rank_order_sutures_boxplot.csv")

# Create a variable containing the suture names
suture_names <- rank$X

# Divide out the marsupials and placentals to make seperate box plots for each
marsupials_rank <- rank_box %>%
  filter(Species %in% c('Bettongia penicillata', 'Monodelphis domestica',
                        'Phascolarctos cincereus', 'Setonix brachyurus', 'Sminthopsis macroura', 
                        'Trichosaurus vulpecula'))
placentals_rank <- rank_box[!(rank_box$Species %in% marsupials_rank$Species),]

# Plot the box plot - ALL
# Suture name on x axis and closure rank (1-31) on the y-axis
# All specimens here
suture_rank_box <- ggplot(rank_box, aes(x=Suture_number, y=closure_rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ 
  xlab("Suture")+ 
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Suture closure rank for all specimens combined")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box

# Plot the box plot - MARSUPIALS
# Suture name on x axis and closure rank (1-31) on the y-axis
# Marsupials specimens only
suture_rank_box_marsupials <- ggplot(marsupials_rank, aes(x=Suture_number, y=closure_rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ 
  xlab("Suture")+ 
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Suture closure rank for all specimens combined - Marsupials only - monotreme excluded")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_marsupials

# Plot the box plot - PLACENTALS
# Suture name on x axis and closure rank (1-31) on the y-axis
# Marsupials specimens only
suture_rank_box_placentals <- ggplot(placentals_rank, aes(x=Suture_number, y=closure_rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ 
  xlab("Suture")+ 
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Suture closure rank for all specimens combined - Placentals only")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_placentals


#########################################################################################################

# STEP 5: Correlation between developmental origin and suture fusion


# Developmental origin = mesoderm, neural crest, boundary

# Plot box plot - coloured by developmental origin
# Suture name on x axis and closure rank (1-31) on the y-axis
# All specimens plotted here
suture_rank_box_dev <- ggplot(rank_box, aes(x=Suture_number, y=closure_rank, fill = Dev_origin))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ 
  xlab("Suture")+ 
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  scale_fill_brewer(palette = "Set2")+
  ggtitle("Suture closure rank for all specimens combined - coloured by developmental origin")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_dev

# Suture closure scores as calculated from "suture_closure_scores.R" code
all_suture_closure <- read.csv("Data/suture_closure_scores_all_specimens_dev_origin.csv")
adult_suture_closure <- read.csv("Data/suture_closure_scores_adults_dev_origin.csv")

# Correlation between closure score (%) and dev origin - all specimens combined
cor_dev_all <- lm(SCS ~ Dev_origin, data = all_suture_closure)
anova(cor_dev_all)

# Correlation between closure score (%) and dev origin - all adults combined
cor_dev_adults <- lm(SCS ~ Dev_origin, data = adult_suture_closure)
anova(cor_dev_adults)


### END ###
