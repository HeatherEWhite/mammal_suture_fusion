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

# Suture rank order closure - pairwise comparisons across species

# Load the data
rank <- read.csv("Analysis/Suture_closure/rank_order_sutures.csv")

# Reorder and rename headings so they match the species names
rank <- rank %>% select(2:23)
rownames(rank) <- suture_names
colnames(rank) <- adults_data$Species

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
# Fills the matrix with all pairwise comparisons from PairwiseComparisons object using $aov.table - pairwise Procrustes ANOVA of ontogenetic trajectories
Pvalues = matrix (nrow = length(names(PairwiseComparisons)), 
                  ncol = length(c(paste("Kendall", names(PairwiseComparisons[[i]]$p.value[1]), sep = "_"), names(PairwiseComparisons[[i]]$estimate[1]))),
                  dimnames = list(c(names(PairwiseComparisons)), c(paste("Kendall", names(PairwiseComparisons[[i]]$p.value[2]), sep = "_"), names(PairwiseComparisons[[i]]$estimate[1])))
)

for (i in 1:length(names(PairwiseComparisons))){
  Pvalues[i,] <- rbind(as.numeric(paste(c(PairwiseComparisons[[i]]$p.value[1], PairwiseComparisons[[i]]$estimate[1]))))
}

write.csv(Pvalues, "Analysis/Suture_closure/Pairwise_comparison_results_rank_order_suture_closure.csv")

# Bonferroni corrected p-values to account for pairwise comparisons
# Create the empty matrix
CorrectedPvalues <- matrix(nrow = length(names(PairwiseComparisons)),
                           ncol = length(c("corrected Kendall p.value")),
                           dimnames = list(c(names(PairwiseComparisons)), c("corrected_Kendall_p.value"))
)

# Fill the matrix with the corrected p-values
CorrectedPvalues[,1] <- p.adjust(Pvalues[,1], method = "bonferroni", n = length(PairwiseComparisons))

write.csv(CorrectedPvalues, "Analysis/Suture_closure/Pairwise_comparisons_corrected_pvalues_rank_order_suture_closure.csv")


#########################################################################################################

# Suture rank order closure - pairwise comparisons between each species and Krogman pattern

# Load the data
rank_krog <- read.csv("Analysis/Suture_closure/Krogman/rank_order_Krogman_region.csv")
rank_krog_embryos <- read.csv("Analysis/Suture_closure/Krogman/rank_order_embryos_Krogman_region.csv")
rank_krog_infants <- read.csv("Analysis/Suture_closure/Krogman/rank_order_infants_Krogman_region.csv")
rank_krog_SA <- read.csv("Analysis/Suture_closure/Krogman/rank_order_subadults_Krogman_region.csv")
rank_krog_adults <- read.csv("Analysis/Suture_closure/Krogman/rank_order_adults_Krogman_region.csv")
primates_krog_adults <- read.csv("Analysis/Suture_closure/Krogman/primate_clousre_score_rank_adults.csv")

average_adults <- read.csv("Analysis/Suture_closure/Krogman/average_adult_score_Krogman_rank.csv")

# Reorder and rename headings so they match the species names
krog_names <- c("Vault", 'Cranial base', "Circum-meatal", "Palate", "Facial", "Cranio-facial")
rank_krog_embryos <- rank_krog_embryos %>% select(2:18)
rank_krog_infants <- rank_krog_infants %>% select(2:24)
rank_krog_SA <- rank_krog_SA %>% select(2:22)
rank_krog_adults <- rank_krog_adults %>% select(2:24)
primates_krog <- primates_krog_adults %>% select(5:6)
rownames(rank_krog_embryos) <- krog_names
rownames(rank_krog_infants) <- krog_names
rownames(rank_krog_SA) <- krog_names
rownames(rank_krog_adults) <- krog_names
rownames(primates_krog) <- krog_names

# Compare the primates order with the Krogman order
primate_cor <- cor.test(primates_krog$Primate_rank, primates_krog$Krogman_rank, method = "kendall", exact =F)
primate_cor


species_embryos <- c("Bettongia penicillata", "Bradypus tridactylus", "Sapajus apella",            
                     "Cyclopes didactylus", "Dasyprocta leporina", "Dasypus novemcinctus",
                     "Epomops franqueti", "Felis catus", "Macroscelides proboscideus",
                     "Phataginus tricuspis", "Monodelphis domestica",
                     "Phacochoerus aethiopicus", "Phascolarctos cincereus",       
                     "Setonix brachyurus", "Sminthopsis macroura",        
                     "Trichosaurus vulpecula", "Krogman pattern")
species_infants <- c("Bettongia penicillata", "Bradypus tridactylus", "Sapajus apella",            
                     "Cyclopes didactylus", "Dasyprocta leporina", "Dasypus novemcinctus",
                     "Epomops franqueti", "Felis catus", "Macroscelides proboscideus",
                     "Phataginus tricuspis", "Microcebus murinus", "Monodelphis domestica",
                     "Mus musculus", "Ornithorhynchus anatinus", "Phacochoerus aethiopicus", 
                     "Phascolarctos cincereus", "Rattus rattus", "Setifer setosus",       
                     "Setonix brachyurus", "Sminthopsis macroura", "Talpa europaea",          
                     "Trichosaurus vulpecula", "Krogman pattern")
species_SA <- c("Bettongia penicillata", "Bradypus tridactylus", "Sapajus apella",            
                "Cyclopes didactylus", "Dasyprocta leporina",
                "Felis catus", "Macroscelides proboscideus",
                "Phataginus tricuspis", "Microcebus murinus", "Monodelphis domestica",
                "Mus musculus", "Ornithorhynchus anatinus", "Phacochoerus aethiopicus", 
                "Phascolarctos cincereus", "Rattus rattus", "Setifer setosus",       
                "Setonix brachyurus", "Sminthopsis macroura", "Talpa europaea",          
                "Trichosaurus vulpecula", "Krogman pattern")
species_adults <- c("Bettongia penicillata", "Bradypus tridactylus", "Sapajus apella",            
                    "Cyclopes didactylus", "Dasyprocta leporina", "Dasypus novemcinctus",
                    "Epomops franqueti", "Felis catus", "Macroscelides proboscideus",
                    "Phataginus tricuspis", "Microcebus murinus", "Monodelphis domestica",
                    "Mus musculus", "Ornithorhynchus anatinus", "Phacochoerus aethiopicus", 
                    "Phascolarctos cincereus", "Rattus rattus", "Setifer setosus",       
                    "Setonix brachyurus", "Sminthopsis macroura", "Talpa europaea",          
                    "Trichosaurus vulpecula", "Krogman pattern")
colnames(rank_krog_embryos) <- species_embryos
colnames(rank_krog_infants) <- species_infants
colnames(rank_krog_SA) <- species_SA
colnames(rank_krog_adults) <- species_adults
species_embryos <- as.data.frame(species_embryos)
species_infants <- as.data.frame(species_infants)
species_SA <- as.data.frame(species_SA)
species_adults <- as.data.frame(species_adults)


# Iteratively create two species pairs and pairwise Kendalls test on rank order of suture closure and save output data from comparison
PairwiseComparisons_krog_adults=list()
i=1
n=2
for (i in 1:17){ 
  for (h in n:17){ 
    if (i==length(species_adults$species_adults)){
      break
    }else
      SpeciesA <- species_adults$species_adults[i]
    SpeciesB <- species_adults$species_adults[h]
    toMatch <- c(SpeciesA, SpeciesB)
    filename <- paste(SpeciesA, SpeciesB, sep=" vs. ")
    cor_pairwise <-cor.test(rank_krog_adults[[SpeciesA]], rank_krog_adults[[SpeciesB]], method = "kendall", exact = F)
    PairwiseComparisons_krog_adults[[filename]]<-cor_pairwise
  }
  n=n+1
}

# Create matrix and fill with values from pairwise comparisons
# Fills the matrix with all pairwise comparisons from PairwiseComparisons object using $aov.table - pairwise Procrustes ANOVA of ontogenetic trajectories
Pvalues_krog_adults = matrix (nrow = length(names(PairwiseComparisons_krog_adults)), 
                              ncol = length(c(paste("Kendall", names(PairwiseComparisons_krog_adults[[i]]$p.value[1]), sep = "_"), names(PairwiseComparisons_krog_adults[[i]]$estimate[1]))),
                              dimnames = list(c(names(PairwiseComparisons_krog_adults)), c(paste("Kendall", names(PairwiseComparisons_krog_adults[[i]]$p.value[2]), sep = "_"), names(PairwiseComparisons_krog_adults[[i]]$estimate[1])))
)

for (i in 1:length(names(PairwiseComparisons_krog_adults))){
  Pvalues_krog_adults[i,] <- rbind(as.numeric(paste(c(PairwiseComparisons_krog_adults[[i]]$p.value[1], PairwiseComparisons_krog_adults[[i]]$estimate[1]))))
}

write.csv(Pvalues_krog_adults, "Analysis/Suture_closure/Krogman/Pairwise_comparison_results_rank_order_adults_Krogman_pattern.csv")

# Bonferroni corrected p-values to account for pairwise comparisons
# Create the empty matrix
CorrectedPvalues_krog_adults <- matrix(nrow = length(names(PairwiseComparisons_krog_adults)),
                                       ncol = length(c("corrected Kendall p.value")),
                                       dimnames = list(c(names(PairwiseComparisons_krog_adults)), c("corrected_Kendall_p.value"))
)

# Fill the matrix with the corrected p-values
CorrectedPvalues_krog_adults[,1] <- p.adjust(Pvalues_krog_adults[,1], method = "bonferroni", n = length(PairwiseComparisons_krog_adults))

write.csv(CorrectedPvalues_krog_adults, "Analysis/Suture_closure/Krogman/Pairwise_comparisons_corrected_pvalues_rank_order_adults_Krogman_pattern.csv")


############

# Average total suture closure score for each Krogman region calculated across the adult only species
# This average score ranked from 1-6 for the Krogman pattern and compared using Kendall's tau to the Krogman pattern
average_adult_vs_Krogman_rank <-cor.test(average_adults$Adult_rank, average_adults$Krogman.rank, method = "kendall", exact = F)

sink("Analysis/Suture_closure/Krogman/Kendalls_tau_average_adult_total_closure_score_ranked_vs_Krogman_rank.txt")
print(average_adult_vs_Krogman_rank)
sink() 


###########

# Average total suture closure score for each Krogman region calculated across the adult only species
# Divided into the marsupials and monotreme as one group and the placentals as another
# This average score ranked from 1-6 for the Krogman pattern and compared using Kendall's tau to the Krogman pattern
# For both marsupials and placentals

marsupial_placental_pattern <- read.csv("Analysis/Suture_closure/Krogman/average_total_suture_closure_adults_marsupials_vs_placentals_Krogman.csv")

marsupial_adult_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Marsupials_adults_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)
placental_adult_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Placentals_adults_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)

sink("Analysis/Suture_closure/Krogman/Kendalls_tau_marsupial_placental_vs_Krogman_rank_adult_average.txt")
print(marsupial_adult_vs_Krogman_rank)
print(placental_adult_vs_Krogman_rank)
sink() 

marsupial_embryo_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Marsupials_embryos_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)
placental_embryo_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Placentals_embryos_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)

marsupial_infant_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Marsupials_infant_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)
placental_infant_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Placentals_infant_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)

marsupial_SA_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Marsupials_SA_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)
placental_SA_vs_Krogman_rank <- cor.test(marsupial_placental_pattern$Placentals_SA_order, marsupial_placental_pattern$Krogman_order, method = "kendall", exact = F)



#########################################################################################################

# Box plot of rank suture closure - like Figure 2 in Goswami et al 2013

# Read in the data
# The csv with the rank order (1-31) for every suture for each species, with all specimens combined for each species
rank_box <- read.csv("Analysis/Suture_closure/rank_order_sutures_boxplot.csv")
# The csv with adult only data, combining the sutures for each Krogman region to get a rank order (1-6)
rank_adults_Krogman_box <- read.csv("Analysis/Suture_closure/Krogman/rank_order_adults_Krogman_region_boxplot.csv")
# The csv with the raw suture closure scores (1-5) for the adults only
OG_suture_score_adults <- read.csv("Data/Suture_closure/suture_closure_status_adults_boxplot.csv")

krog_names <- c("Vault", 'Cranial base', "Circum-meatal", "Palate", "Facial", "Cranio-facial")

# Divide out the marsupials and placentals to make seperate box plots for each
# Monotreme included in marsupials
marsupials_monotreme_rank <- rank_box %>%
  filter(Species %in% c('Bettongia penicillata', 'Monodelphis domestica', 'Ornithorhynchus anatinus',
                        'Phascolarctos cincereus', 'Setonix brachyurus', 'Sminthopsis macroura', 
                        'Trichosaurus vulpecula'))
marsupials_rank <- rank_box %>%
  filter(Species %in% c('Bettongia penicillata', 'Monodelphis domestica',
                        'Phascolarctos cincereus', 'Setonix brachyurus', 'Sminthopsis macroura', 
                        'Trichosaurus vulpecula'))
placentals_rank <- rank_box[!(rank_box$Species %in% marsupials_rank$Species),]

# Box plot for the suture on x axis and closure rank (1-31) on the y-axis
# All specimens and species
suture_rank_box <- ggplot(rank_box, aes(x=Suture_number, y=closure_rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Suture closure rank for all specimens combined")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box

# Box plot for the suture on x axis and closure rank (1-31) on the y-axis
# Marsupials only - monotreme excluded
suture_rank_box_marsupials <- ggplot(marsupials_rank, aes(x=Suture_number, y=closure_rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Suture closure rank for all specimens combined - Marsupials only - monotreme excluded")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_marsupials

# Box plot for the suture on x axis and closure rank (1-31) on the y-axis
# Placentals only
suture_rank_box_placentals <- ggplot(placentals_rank, aes(x=Suture_number, y=closure_rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Suture closure rank for all specimens combined - Placentals only")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_placentals

###########

# Coloured by developmental origin

# Box plot for the suture on x axis and closure rank (1-31) on the y-axis
# All specimens and species
suture_rank_box_dev <- ggplot(rank_box, aes(x=Suture_number, y=closure_rank, fill = Dev_origin))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  scale_fill_brewer(palette = "Set2")+
  ggtitle("Suture closure rank for all specimens combined - coloured by developmental origin")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_dev

# Box plot for the suture on x axis and closure rank (1-31) on the y-axis
# Marsupials only - monotreme excluded
suture_rank_box_marsupials_dev <- ggplot(marsupials_rank, aes(x=Suture_number, y=closure_rank, fill = Dev_origin))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  scale_fill_brewer(palette = "Set2")+
  ggtitle("Suture closure rank for all specimens combined - Marsupials only - monotreme excluded - coloured by developmental origin")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_marsupials_dev

# Box plot for the suture on x axis and closure rank (1-31) on the y-axis
# Placentals only
suture_rank_box_placentals_dev <- ggplot(placentals_rank, aes(x=Suture_number, y=closure_rank, fill = Dev_origin))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = suture_names)+
  scale_fill_brewer(palette = "Set2")+
  ggtitle("Suture closure rank for all specimens combined - Placentals only - coloured by developmental origin")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
suture_rank_box_placentals_dev

###########

# Adults only

# Box plot for the Krogman region on x axis and closure rank (1-6) on the y-axis
# Adults only
krog_rank_adults_box <- ggplot(rank_adults_Krogman_box, aes(x=Region_number, y=Rank))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Krogman region")+ # changing the x-label, can add the measurement details here
  ylab("Closure rank")+
  scale_x_discrete(labels = krog_names)+
  ggtitle("Suture closure rank for adults only")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
krog_rank_adults_box

# Box plot for adults only for all sutures on x-axis and raw suture closure score (1-5) on the y-axis
# Adults only
raw_suture_closure_adults_box <- ggplot(OG_suture_score_adults, aes(x=Suture_number, y=Closure_score))+ 
  geom_boxplot()+
  theme_classic(base_size = 14)+ # change the background of the graph, theme_classic = white background, theme_bw = white with grid lines, saying base size here changes the size of everything
  xlab("Suture")+ # changing the x-label, can add the measurement details here
  ylab("Raw suture closure score")+
  scale_x_discrete(labels = suture_names)+
  ggtitle("Raw suture closure score (1-5) for adults only")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(angle = 90, size = 10, vjust = 0, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
raw_suture_closure_adults_box

#######

# Relationship between developmental origin and suture fusion

all_suture_closure <- read.csv("Analysis/Suture_closure/average_suture_closure_per_suture_all_specimens.csv")
adult_suture_closure <- read.csv("Analysis/Suture_closure/average_suture_closure_per_suture_adults.csv")

# Correlation between closure score (%) and dev origin - all specimens combined
cor_dev_all <- lm(V1 ~ Dev_origin, data = all_suture_closure)
anova(cor_dev_all)

# Correlation between closure score (%) and dev origin - all adults combined
cor_dev_adults <- lm(V1 ~ Dev_origin, data = adult_suture_closure)
anova(cor_dev_adults)

sink("Analysis/Suture_closure/ANOVA_closure_score_vs_dev_origin.txt")
print("All specimens")
print(anova(cor_dev_all))
print("Adults only")
print(anova(cor_dev_adults))
sink

# Dividing the dataset into placentals
placental_sutures <- all_data %>% 
  filter(Clade == "Placental") %>%
  select(c(13:43))
placental_names <- all_data %>%
  filter(Clade == "Placental") %>%
  select(,1)
rownames(placental_sutures) <- placental_names$Specimen

# Dividing the dataset into marsupials
marsupial_sutures <- all_data %>% 
  filter(Clade == "Marsupial") %>%
  select(c(13:43))
marsupial_names <- all_data %>%
  filter(Clade == "Marsupial") %>%
  select(,1)
rownames(marsupial_sutures) <- marsupial_names$Specimen

# Dividing the adult dataset into placentals
placental_sutures_adults <- adults_data %>% 
  filter(Clade == "Placental") %>%
  select(c(13:43))
placental_names_adults <- adults_data %>%
  filter(Clade == "Placental") %>%
  select(,2)
rownames(placental_sutures_adults) <- placental_names_adults$Species

# Dividing the dataset into marsupials
marsupial_sutures_adults <- adults_data %>% 
  filter(Clade == "Marsupial") %>%
  select(c(13:43))
marsupial_names_adults <- adults_data %>%
  filter(Clade == "Marsupial") %>%
  select(,2)
rownames(marsupial_sutures_adults) <- marsupial_names_adults$Species

# Calculate the score for placentals
# Average suture closure
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_placentals_suture <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == placental_sutures[,i])))/98)*100
  score4[i,] <- (sum(length(which(4 == placental_sutures[,i])))/98)*100
  score3[i,] <- (sum(length(which(3 == placental_sutures[,i])))/98)*100
  score2[i,] <- (sum(length(which(2 == placental_sutures[,i])))/98)*100
  MSC_placentals_suture[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_placentals_suture
rownames(MSC_placentals_suture) <- suture_names
colnames(MSC_placentals_suture) <- "MSC"
MSC_placentals_suture <- as_tibble(MSC_placentals_suture) %>% 
  mutate(Dev_origin = all_suture_closure$Dev_origin)
MSC_placentals_suture <- as.data.frame(MSC_placentals_suture)
rownames(MSC_placentals_suture) <- suture_names

# Calculate the score for marsupials
# Average suture closure
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_marsupials_suture <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == marsupial_sutures[,i])))/63)*100
  score4[i,] <- (sum(length(which(4 == marsupial_sutures[,i])))/63)*100
  score3[i,] <- (sum(length(which(3 == marsupial_sutures[,i])))/63)*100
  score2[i,] <- (sum(length(which(2 == marsupial_sutures[,i])))/63)*100
  MSC_marsupials_suture[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_marsupials_suture
rownames(MSC_marsupials_suture) <- suture_names
colnames(MSC_marsupials_suture) <- "MSC"
MSC_marsupials_suture <- as_tibble(MSC_marsupials_suture) %>% 
  mutate(Dev_origin = all_suture_closure$Dev_origin)
MSC_marsupials_suture <- as.data.frame(MSC_marsupials_suture)
rownames(MSC_marsupials_suture) <- suture_names

# Calculate the score for placentals - adults only
# Average suture closure
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_placentals_suture_adults <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == placental_sutures_adults[,i])))/15)*100
  score4[i,] <- (sum(length(which(4 == placental_sutures_adults[,i])))/15)*100
  score3[i,] <- (sum(length(which(3 == placental_sutures_adults[,i])))/15)*100
  score2[i,] <- (sum(length(which(2 == placental_sutures_adults[,i])))/15)*100
  MSC_placentals_suture_adults[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_placentals_suture_adults
rownames(MSC_placentals_suture_adults) <- suture_names
colnames(MSC_placentals_suture_adults) <- "MSC"
MSC_placentals_suture_adults <- as_tibble(MSC_placentals_suture_adults) %>% 
  mutate(Dev_origin = all_suture_closure$Dev_origin)
MSC_placentals_suture_adults <- as.data.frame(MSC_placentals_suture_adults)
rownames(MSC_placentals_suture_adults) <- suture_names

# Calculate the score for marsupials - adults only
# Average suture closure
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_marsupials_suture_adults <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == marsupial_sutures_adults[,i])))/6)*100
  score4[i,] <- (sum(length(which(4 == marsupial_sutures_adults[,i])))/6)*100
  score3[i,] <- (sum(length(which(3 == marsupial_sutures_adults[,i])))/6)*100
  score2[i,] <- (sum(length(which(2 == marsupial_sutures_adults[,i])))/6)*100
  MSC_marsupials_suture_adults[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_marsupials_suture_adults
rownames(MSC_marsupials_suture_adults) <- suture_names
colnames(MSC_marsupials_suture_adults) <- "MSC"
MSC_marsupials_suture_adults <- as_tibble(MSC_marsupials_suture_adults) %>% 
  mutate(Dev_origin = all_suture_closure$Dev_origin)
MSC_marsupials_suture_adults <- as.data.frame(MSC_marsupials_suture_adults)
rownames(MSC_marsupials_suture_adults) <- suture_names

# Correlation between closure score (%) and dev origin - all specimens combined
cor_dev_all_placentals <- lm(MSC ~ Dev_origin, data = MSC_placentals_suture)
anova(cor_dev_all_placentals)
cor_dev_all_marsupials <- lm(MSC ~ Dev_origin, data = MSC_marsupials_suture)
anova(cor_dev_all_marsupials)

# Correlation between closure score (%) and dev origin - adults only
cor_dev_all_placentals_adults <- lm(MSC ~ Dev_origin, data = MSC_placentals_suture_adults)
anova(cor_dev_all_placentals_adults)
cor_dev_all_marsupials_adults <- lm(MSC ~ Dev_origin, data = MSC_marsupials_suture_adults)
anova(cor_dev_all_marsupials_adults)


sink("Analysis/Suture_closure/ANOVA_developmental_origin_vs_closure_score.txt")
print("MSC of all specimens placentals vs developmental origin")
print(anova(cor_dev_all_placentals))
print("MSC of all specimens marsupials vs developmental origin")
print(anova(cor_dev_all_marsupials))
print("MSC of all ADULT specimens placentals vs developmental origin")
print(anova(cor_dev_all_placentals_adults))
print("MSC of all ADULT specimens marsupials vs developmental origin")
print(anova(cor_dev_all_marsupials_adults))
sink() 




#########################################################################################################

# Maximum suture closure - like fig 3 in Goswami et al 2013

MSC_species <- as.data.frame(MSC_species) %>% mutate(logCS = adults_data$CS_logged, Species = adults_data$Species)
MSC_names <- c("MSC", "logCS", "Species")
colnames(MSC_species) <- MSC_names

# Scatter plot CS vs MSC
CS_vs_MSC <- ggplot(MSC_species, aes(x = logCS, y = MSC, label = Species))+
  geom_point(size = 2)+
  geom_text_repel(aes(fontface="italic"), size = 3)+
  theme_classic(base_size = 12)+
  xlab("Size (logged centroid size)")+
  ylab("Maximum suture closure")+
  ggtitle("Log CS vs MSC for adults only")
CS_vs_MSC


### END ###
