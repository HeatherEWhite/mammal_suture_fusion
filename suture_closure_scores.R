# Title: Calculation of suture closure scores

# Name: Heather White

# Date created: 02/12/21

# Last modified: 25/10/22

# License: MIT license


# Calculating suture closure scores and total suture closure scores 
# following protocol in Goswami et al. 2013

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


#######################################################################################################

# STEP 1: Load and format the data

# Read in the raw suture closure values as scored using CT scans - full dataset
all_data <- read.csv("Data/raw_closure_status.csv")
# Read in the raw suture closure values as scored using CT scans - adults only
adults_data <- read.csv("Data/suture_closure_status_adults.csv")

# Select the marsupials and placentals
marsupials <- adults_data %>% filter(Clade == "Marsupial")
placentals <- adults_data %>% filter(Clade == "Placental")


# Input the closure information into a dataframe
# Adults
adults_closure <- adults_data %>%
  select(c(13:43))
rownames(adults_closure) <- adults_data$Species
# Marsupials
marsupials_closure <- marsupials %>%
  select(c(13:43))
rownames(marsupials_closure) <- marsupials$Species
# Placentals
placentals_closure <- placentals %>%
  select(c(13:43))
rownames(placentals_closure) <- placentals$Species
# All
all_closure <- all_data %>%
  select(c(13:43))
rownames(all_closure) <- all_data$Specimen


# Save the suture names as a variable
suture_names <- rownames(t(adults_closure))


#######################################################################################################

# STEP 2 - Calculate total suture closure scores (adults only)

# Total suture closure scores combine all sutures for one specimen 
# to get an overall level of closure across the skull

# This calculation uses equation 1 in White et al. 2023


# Calculate total suture fusion combining all sutures for each species seperately
score5 <- array(dim = c(22,1))
score4 <- array(dim = c(22,1))
score3 <- array(dim = c(22,1))
score2 <- array(dim = c(22,1))
MSC_species <- array(dim = c(22,1))
for(i in 1:22)
{
  score5[i,] <- (sum(length(which(5 == adults_closure[i,])))/31)*100
  score4[i,] <- (sum(length(which(4 == adults_closure[i,])))/31)*100
  score3[i,] <- (sum(length(which(3 == adults_closure[i,])))/31)*100
  score2[i,] <- (sum(length(which(2 == adults_closure[i,])))/31)*100
  MSC_species[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_species
rownames(MSC_species) <- adults_data$Species

# Save the data if necessary
#write.csv(MSC_species, "total_suture_closure_adult_species.csv")

#######################################################################################################

# STEP 3 - Calculate suture closure scores (adults only)

# Suture closure score combines all specimens for one suture
# to get an level of closure across for one suture

# This calculation uses equation 1 in White et al. 2023

# Calculate suture closure score combining all specimens for each suture seperately - full adult dataset
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_sutures <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == adults_closure[,i])))/22)*100
  score4[i,] <- (sum(length(which(4 == adults_closure[,i])))/22)*100
  score3[i,] <- (sum(length(which(3 == adults_closure[,i])))/22)*100
  score2[i,] <- (sum(length(which(2 == adults_closure[,i])))/22)*100
  MSC_sutures[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_sutures
rownames(MSC_sutures) <- suture_names

# Save the data if necessary
#write.csv(MSC_sutures, "Analysis/Suture_closure/average_suture_closure_per_suture_adults.csv")

########

# Calculate suture closure score combining all adult marsupial specimens for each suture seperately
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_sutures_marsupials <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == marsupials_closure[,i])))/6)*100
  score4[i,] <- (sum(length(which(4 == marsupials_closure[,i])))/6)*100
  score3[i,] <- (sum(length(which(3 == marsupials_closure[,i])))/6)*100
  score2[i,] <- (sum(length(which(2 == marsupials_closure[,i])))/6)*100
  MSC_sutures_marsupials[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_sutures_marsupials
rownames(MSC_sutures_marsupials) <- suture_names

# Save the data if necessary
#write.csv(MSC_sutures_marsupials, "Analysis/Suture_closure/average_suture_closure_per_suture_adults_marsupials.csv")

########

# Calculate maximum suture fusion combining all adult placental specimens for each suture seperately
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_sutures_placentals <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == placentals_closure[,i])))/15)*100
  score4[i,] <- (sum(length(which(4 == placentals_closure[,i])))/15)*100
  score3[i,] <- (sum(length(which(3 == placentals_closure[,i])))/15)*100
  score2[i,] <- (sum(length(which(2 == placentals_closure[,i])))/15)*100
  MSC_sutures_placentals[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_sutures_placentals
rownames(MSC_sutures_placentals) <- suture_names

# Save the data if necessary
#write.csv(MSC_sutures_placentals, "Analysis/Suture_closure/average_suture_closure_per_suture_adults_placentals.csv")


#######################################################################################################

# STEP 4: Calculate total suture closure scores (full developmental dataset)

# Total suture closure score combines all sutures for one specimen 
# to get an overall level of closure across the skull

# This calculation uses equation 1 in White et al. 2023


# Calculate total suture closure scores by combining all sutures for each specimen seperately
score5 <- array(dim = c(165,1))
score4 <- array(dim = c(165,1))
score3 <- array(dim = c(165,1))
score2 <- array(dim = c(165,1))
MSC_all_specimens <- array(dim = c(165,1))
for(i in 1:165)
{
  score5[i,] <- (sum(length(which(5 == all_closure[i,])))/31)*100
  score4[i,] <- (sum(length(which(4 == all_closure[i,])))/31)*100
  score3[i,] <- (sum(length(which(3 == all_closure[i,])))/31)*100
  score2[i,] <- (sum(length(which(2 == all_closure[i,])))/31)*100
  MSC_all_specimens[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_all_specimens
rownames(MSC_all_specimens) <- all_data$Specimen

# Save the data if necessary
#write.csv(MSC_all_specimens, "Analysis/Suture_closure/overall_suture_closure_all_specimens.csv")


#######################################################################################################

# STEP 5 - Calculate suture closure scores (full developmental dataset)

# Suture closure score combines all specimens for one suture
# to get an level of closure across for one suture

# This calculation uses equation 1 in White et al. 2023

# Calculate suture closure scores by combining all specimens for each suture seperately
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_all_specimens_suture <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == all_closure[,i])))/165)*100
  score4[i,] <- (sum(length(which(4 == all_closure[,i])))/165)*100
  score3[i,] <- (sum(length(which(3 == all_closure[,i])))/165)*100
  score2[i,] <- (sum(length(which(2 == all_closure[,i])))/165)*100
  MSC_all_specimens_suture[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_all_specimens_suture
rownames(MSC_all_specimens_suture) <- suture_names

# Save the data if necessary
#write.csv(MSC_all_specimens_suture, "Analysis/Suture_closure/average_suture_closure_per_suture_all_specimens.csv")


#######################################################################################################

# STEP 6: Subset the dataset into developmental categories and individual species


# Subset out the different age categories
# Fetal
fetal <- all_data %>%
   filter(Discrete_age == "E")
fetal_closure <- fetal %>%
  select(c(12:42))
rownames(fetal_closure) <- fetal$Specimen
# Infant
infant <- all_data %>%
  filter(Discrete_age == "I")
infant_closure <- infant %>%
  select(c(12:42))
rownames(infant_closure) <- infant$Specimen
# Subadult
SA <- all_data %>%
  filter(Discrete_age == "SA")
SA_closure <- SA %>%
  select(c(12:42))
rownames(SA_closure) <- SA$Specimen

#######

# Sort the data by species, to get a summary of each suture closure for every species
all_closure_df <- all_closure %>% 
  mutate(Species = all_data$Species)
Bettongia <- all_closure_df %>%
  filter(Species == "Bettongia penicillata")
Bradypus <- all_closure_df %>%
  filter(Species == "Bradypus tridactylus")
Sapajus <- all_closure_df %>%
  filter(Species == "Sapajus apella")
Cyclopes <- all_closure_df %>%
  filter(Species == "Cyclopes didactylus")
Dasprocta <- all_closure_df %>%
  filter(Species == "Dasyprocta leporina")
Dasypus <- all_closure_df %>%
  filter(Species == "Dasypus novemcinctus")
Epomops <- all_closure_df %>%
  filter(Species == "Epomops franqueti")
Felis <- all_closure_df %>%
  filter(Species == "Felis catus")
Macroscelides <- all_closure_df %>%
  filter(Species == "Macroscelides proboscideus")
Phataginus <- all_closure_df %>%
  filter(Species == "Phataginus tricuspis")
Microcebus <- all_closure_df %>%
  filter(Species == "Microcebus murinus")
Monodelphis <- all_closure_df %>%
  filter(Species == "Monodelphis domestica")
Mus <- all_closure_df %>%
  filter(Species == "Mus musculus")
Ornithorhynchus <- all_closure_df %>%
  filter(Species == "Ornithorhynchus anatinus")
Phacochoerus <- all_closure_df %>%
  filter(Species == "Phacochoerus aethiopicus")
Phascolarctos <- all_closure_df %>%
  filter(Species == "Phascolarctos cincereus")
Rattus <- all_closure_df %>%
  filter(Species == "Rattus rattus")
Setifer <- all_closure_df %>%
  filter(Species == "Setifer setosus")
Setonix <- all_closure_df %>%
  filter(Species == "Setonix brachyurus")
Sminthopsis <- all_closure_df %>%
  filter(Species == "Sminthopsis macroura")
Talpa <- all_closure_df %>%
  filter(Species == "Talpa europaea")
Trichosaurus <- all_closure_df %>%
  filter(Species == "Trichosaurus vulpecula")


#######################################################################################################

# STEP 7: Calculate total suture closure scores for each specimen at each age category


# Fetal
# Calculate total suture closure scores by combining all sutures for each fetal specimen
score5 <- array(dim = c(28,1))
score4 <- array(dim = c(28,1))
score3 <- array(dim = c(28,1))
score2 <- array(dim = c(28,1))
MSC_fetal <- array(dim = c(28,1))
for(i in 1:28)
{
  score5[i,] <- (sum(length(which(5 == fetal_closure[i,])))/31)*100
  score4[i,] <- (sum(length(which(4 == fetal_closure[i,])))/31)*100
  score3[i,] <- (sum(length(which(3 == fetal_closure[i,])))/31)*100
  score2[i,] <- (sum(length(which(2 == fetal_closure[i,])))/31)*100
  MSC_fetal[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_fetal
rownames(MSC_fetal) <- fetal$Specimen

# Add the species names into the dataframe
MSC_fetal <- as.data.frame(MSC_fetal) %>%
  mutate(Species = fetal$Species)

# Calculate the average total suture closure score for each species - for the fetal specimens only
mean_fetal <- MSC_fetal %>% 
  group_by(Species) %>% 
  summarize(Mean_closure = mean(V1))

# Save the data if necessary
#write.csv(mean_embryo, "Analysis/Suture_closure/mean_closure_species_embryos.csv")


#######

# Infant
# Calculate total suture closure scores by combining all sutures for each infant specimen
score5 <- array(dim = c(95,1))
score4 <- array(dim = c(95,1))
score3 <- array(dim = c(95,1))
score2 <- array(dim = c(95,1))
MSC_infant <- array(dim = c(95,1))
for(i in 1:95)
{
  score5[i,] <- (sum(length(which(5 == infant_closure[i,])))/31)*100
  score4[i,] <- (sum(length(which(4 == infant_closure[i,])))/31)*100
  score3[i,] <- (sum(length(which(3 == infant_closure[i,])))/31)*100
  score2[i,] <- (sum(length(which(2 == infant_closure[i,])))/31)*100
  MSC_infant[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_infant
rownames(MSC_infant) <- infant$Specimen

# Add the species names into the dataframe
MSC_infant <- as.data.frame(MSC_infant) %>%
  mutate(Species = infant$Species)

# Calculate the average total suture closure score for each species - for infants only
mean_infant <- MSC_infant %>% 
  group_by(Species) %>% 
  summarize(Mean_closure = mean(V1))

# Save the data if necessary
#write.csv(mean_infant, "Analysis/Suture_closure/mean_closure_species_infants.csv")


#######

# Subadults
# Calculate total suture closure scores by combining all sutures for each subadult specimen
score5 <- array(dim = c(20,1))
score4 <- array(dim = c(20,1))
score3 <- array(dim = c(20,1))
score2 <- array(dim = c(20,1))
MSC_SA <- array(dim = c(20,1))
for(i in 1:20)
{
  score5[i,] <- (sum(length(which(5 == SA_closure[i,])))/31)*100
  score4[i,] <- (sum(length(which(4 == SA_closure[i,])))/31)*100
  score3[i,] <- (sum(length(which(3 == SA_closure[i,])))/31)*100
  score2[i,] <- (sum(length(which(2 == SA_closure[i,])))/31)*100
  MSC_SA[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_SA
rownames(MSC_SA) <- SA$Specimen

# Add the species names into the dataframe
MSC_SA <- as.data.frame(MSC_SA) %>%
  mutate(Species = SA$Species)

# Calculate the average total suture closure score for each species - for subadult specimens only
mean_SA <- MSC_SA %>% 
  group_by(Species) %>% 
  summarize(Mean_closure = mean(V1))

# Save the data if necessary
#write.csv(mean_SA, "Analysis/Suture_closure/mean_closure_species_subadults.csv")


#######################################################################################################

# STEP 8: Calculate suture closure scores for each suture seperately for each species 

# Calculate for each species for each suture 
# Need to change the name of the variable and species dataset each time
# Also need to change the number each score is divided by - number of specimens for each species
# The example below is for Trichosurus
score5 <- array(dim = c(31,1))
score4 <- array(dim = c(31,1))
score3 <- array(dim = c(31,1))
score2 <- array(dim = c(31,1))
MSC_Trichosaurus <- array(dim = c(31,1))
for(i in 1:31)
{
  score5[i,] <- (sum(length(which(5 == Trichosaurus[,i])))/9)*100
  score4[i,] <- (sum(length(which(4 == Trichosaurus[,i])))/9)*100
  score3[i,] <- (sum(length(which(3 == Trichosaurus[,i])))/9)*100
  score2[i,] <- (sum(length(which(2 == Trichosaurus[,i])))/9)*100
  MSC_Trichosaurus[i,] <- score5[i,] + (0.75 * score4[i,]) + (0.5 * score3[i,]) + (0.25 * score2[i,])
} 
MSC_Trichosaurus
rownames(MSC_Trichosaurus) <- suture_names

# Once run for each species then combine into a dataset
MSC_sutures_species <- cbind(MSC_Bettongia, MSC_Bradypus, MSC_Sapajus, MSC_Cyclopes, MSC_Dasprocta, MSC_Dasypus,
                             MSC_Epomops, MSC_Felis, MSC_Macroscelides, MSC_Phataginus, MSC_Microcebus,
                             MSC_Monodelphis, MSC_Mus, MSC_Ornithorhynchus, MSC_Phacochoerus, 
                             MSC_Phascolarctos, MSC_Rattus, MSC_Setifer, MSC_Setonix,
                             MSC_Sminthopsis, MSC_Talpa, MSC_Trichosaurus)
colnames(MSC_sutures_species) <- adults_data$Species

# Save the data if necessary
#write.csv(MSC_sutures_species, "Analysis/Suture_closure/suture_closure_score_species_average_per_suture.csv")



#######################################################################################################

# STEP 9: Plot the total suture closure scores against age


# Setup the data for the plots

# Load the data
# Total suture closure scores for each age category
# This .csv combines the results from step 7 and adults only (step 2)
suture_closure_av_age <- read.csv("Data/suture_closure_average_by_age.csv")

# Calculated total suture closure scores for each specimen
# Combining all sutures for each Krogman region (Krogman, 1930)
# Suture closure score for every specimen for each of the six regions
suture_closure_av_krog <- read.csv("Data/suture_closure_score_per_specimen_per_Krogman_region.csv")


# Create a colour palette for species
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
  "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
  "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
  "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
  "gold1", # Monotremata
  "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:22, 1, as.matrix(1:22), col = mypalette_species, xlab = "Species",
      ylab = "", yaxt = "n")

# Convert the dataset dataframe to a tibble so I can reorder the species
MSC_all_specimens_tbl <- as.data.frame(MSC_all_specimens) %>%
  mutate(Percent_adult = all_data$CS_percent_adult, Species = all_data$Species, logCS = all_data$CS_logged)
names(MSC_all_specimens_tbl)[names(MSC_all_specimens_tbl) == "V1"] <- "total_closure_score"
MSC_all_specimens_tbl <- as_tibble(MSC_all_specimens_tbl)
suture_closure_av_age_tbl <- as_tibble(suture_closure_av_age)
suture_closure_av_krog_tbl <- as_tibble(suture_closure_av_krog)

# Order the species to match the colours above
MSC_all_specimens_tbl$Species <- factor(MSC_all_specimens_tbl$Species,
                                   levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                              "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                              "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                              "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                              "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                              "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                              "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                              "Dasypus novemcinctus"),
                                   ordered = TRUE)

suture_closure_av_age_tbl$Species <- factor(suture_closure_av_age_tbl$Species,
                                        levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                                   "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                                   "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                                   "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                                   "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                                   "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                                   "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                                   "Dasypus novemcinctus"),
                                        ordered = TRUE)

suture_closure_av_krog_tbl$Species <- factor(suture_closure_av_krog_tbl$Species,
                                        levels = c("Macroscelides proboscideus", "Setifer setosus", "Sapajus apella",            
                                                   "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                                   "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                                   "Phataginus tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                                   "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                                   "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                                   "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                                   "Dasypus novemcinctus"),
                                        ordered = TRUE)


#######

# Plot the scatter plots


# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
suture_closure_ontogeny <- ggplot(MSC_all_specimens_tbl, aes(x = Percent_adult, y = total_closure_score, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age")
suture_closure_ontogeny

# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
# For each Krogman region 
# 1) Cranial vault suture closure
suture_closure_vault <- ggplot(suture_closure_av_krog_tbl, aes(x = CS_percent_adult, y = Vault, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age - for the cranial vault sutures only")
suture_closure_vault

# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
# For each Krogman region 
# 2) Cranial base suture closure
suture_closure_base <- ggplot(suture_closure_av_krog_tbl, aes(x = CS_percent_adult, y = Cranial_base, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age - for the cranial base sutures only")
suture_closure_base

# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
# For each Krogman region 
# 3) Circum-meatal suture closure
suture_closure_circum <- ggplot(suture_closure_av_krog_tbl, aes(x = CS_percent_adult, y = Circum.meatal, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age - for the circum-meatal sutures only")
suture_closure_circum

# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
# For each Krogman region 
# 4) Palate suture closure
suture_closure_palate <- ggplot(suture_closure_av_krog_tbl, aes(x = CS_percent_adult, y = Palate, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age - for the palate sutures only")
suture_closure_palate

# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
# For each Krogman region 
# 5) Facial suture closure
suture_closure_facial <- ggplot(suture_closure_av_krog_tbl, aes(x = CS_percent_adult, y = Facial, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age - for the facial sutures only")
suture_closure_facial

# Scatter graph for continuous age (% of adult) vs total suture closure score - coloured by species
# For each Krogman region 
# 6) Craniofacial suture closure
suture_closure_craniofacial <- ggplot(suture_closure_av_krog_tbl, aes(x = CS_percent_adult, y = Cranio.facial, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age - for the craniofacial sutures only")
suture_closure_craniofacial

# Scatter graph for discrete age categories vs total suture closure score - coloured by species
suture_closure_age_cat <- ggplot(suture_closure_av_age_tbl, aes(x = Age, y = total_closure_score, colour = Species))+
  geom_point(size = 2)+
  geom_line()+
  theme_classic(base_size = 12)+
  scale_colour_manual(values = mypalette_species)+
  xlab("Age (% of adult CS)")+
  ylab("Total suture closure score %")+
  ggtitle("Suture closure against age")
suture_closure_age_cat

### END ###
