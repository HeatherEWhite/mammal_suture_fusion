# Title: Fusion vs skull size

# Name: Heather White

# Date created: 02/12/21

# Last modified: 25/10/22

# License: MIT license


# Correlation between suture fusion and skull size

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

# Total suture closure scores for adults and all specimens
# Summing all sutures per specimen
TSCS_adults <- read.csv("Data/total_suture_closure_scores_adults.csv")
TSCS_all <- read.csv("Data/total_suture_closure_scores_all.csv")

# Load the specimen info data
adults_data <- read.csv("Data/Specimen_info_adults.csv")
all_data <- read.csv("Data/Specimen_info.csv")


# Format the TSCS with the specimen info (skull size - CS logged)
# Adults
combined_data_adults <- as.data.frame(TSCS_adults) %>%
  mutate(logCS = adults_data$CS_logged)
names(combined_data_adults)[names(combined_data_adults) == "closure_score_."] <- "total_closure_score"
names(combined_data_adults)[names(combined_data_adults) == "X"] <- "Species"
# All specimens
combined_data <- as.data.frame(TSCS_all) %>%
  mutate(logCS = all_data$CS_logged, Species = all_data$Species)
names(combined_data)[names(combined_data) == "closure_score_."] <- "total_closure_score"
combined_data_tbl <- as_tibble(combined_data)

# Create a colour palette for species
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                       "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                       "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                       "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:22, 1, as.matrix(1:22), col = mypalette_species, xlab = "Species",
      ylab = "", yaxt = "n")

# Order the species to match the colours above
combined_data_tbl$Species <- factor(combined_data_tbl$Species,
                                    levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                               "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                               "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                               "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                               "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                               "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                               "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                               "Dasypus novemcinctus"),
                                    ordered = TRUE)

#########################################################################################################

# STEP 2: Plot the data - skull size vs total suture closure score 


# Scatter plot CS vs TSCS
CS_vs_MSC <- ggplot(combined_data_adults, aes(x = logCS, y = total_closure_score, label = Species))+
  geom_point(size = 2)+
  geom_text_repel(aes(fontface="italic"), size = 3)+
  theme_classic(base_size = 12)+
  xlab("Size (logged centroid size)")+
  ylab("Total suture closure score")+
  ggtitle("Log CS vs MSC for adults only")
CS_vs_MSC


#########################################################################################################

# STEP 3: Spearman's rank correlation between total suture closure scores and skull size - ADULTS ONLY


# Plot the results - adults only
ggscatter(combined_data_adults, x = "logCS", y = "total_closure_score", 
          add = "reg.line", conf.int = TRUE, 
          add.params = list(linetype=c("dashed")),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "Total suture closure score",
          title = "Correlation between suture fusion and skull size - adults only")


# Spearmans rank correlation for adults only
cor_adults <-cor.test(combined_data_adults$total_closure_score, combined_data_adults$logCS, method = "spearman", conf.level = .95, exact = F)
cor_adults



#########################################################################################################

# STEP 4: Spearman's rank correlation between total suture closure scores and skull size - FULL DEVELOPMENTAL DATASET

# Plot the results - all specimens
ggscatter(combined_data, x = "logCS", y = "total_closure_score", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "Total suture closure score",
          title = "Correlation between suture fusion and skull size - all specimens")


# Spearmans rank correlation for all specimens
cor_all <-cor.test(combined_data$total_closure_score, combined_data$logCS, method = "spearman", conf.level = .95)
cor_all


#########################################################################################################

# STEP 5: Spearman's rank correlation between total suture closure scores and skull size - for each species seperately


# Plot the results as regressions for each species seperately and colour by species
ggscatter(combined_data_tbl, x = "logCS", y = "total_closure_score", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(linetype=c("dashed")),
          color = "Species",
          palette = mypalette_species,
          legend = "right",
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Logged Centroid Size", 
          ylab = "Total suture closure score",
          title = "Correlation between suture fusion and skull size - all specimens")

# Spearmans rank correlation for all specimens, done by species
Spearman_species <- combined_data %>% 
  group_by(Species) %>% 
  do(tidy(cor.test(.$total_closure_score, .$logCS, method = "spearman", conf.level = .95, exact = F)))
as.data.frame(Spearman_species)

# Save data if necessary
#write.csv(Spearman_species, file = "Analysis/Suture_closure/Spearmans_rank_suture_closure_score_vs_CS_within_species.csv")


#########################################################################################################

# STEP 6: Phylogenetic correction for adult correlation


# Read in phylogeny and taxa names
my_tree <- read.nexus("Data/my_mammal_tree.nexus")
tree_names <- read.csv("Data/tree_taxa_names.csv")

# Add the phylogeny taxa names to the adult only dataframe
pglsdata <- combined_data_adults %>% mutate(Taxa_names = tree_names$Taxa_names)


# Run the phylogenetic generalised least squares (pgls) model
pglsModel <- gls(total_closure_score ~ logCS, correlation = corBrownian(phy = my_tree),
                 data = pglsdata, method = "ML")
summary(pglsModel)


# Plot the phylogenetic generalised least squares regression
plot(total_closure_score ~ logCS, data = pglsdata)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])
title("Phylogenetic generalisaed least squares for total suture closure vs logCS of adults only")


### END ###
