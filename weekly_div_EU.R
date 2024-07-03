# Genetic Distance over Time: Weekly Data
# All EU sequences
# German vs. EU sequences

library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(plotly)
library(adegenet)
library(usedist)
library(lubridate)
library(MMWRweek)
library(Rmisc)

# Read metadata

meta_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv") #already in alphabetical order
meta_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv")

meta_rsvA$Collection_Date <- as.Date(meta_rsvA$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
meta_rsvB$Collection_Date <- as.Date(meta_rsvB$Collection_Date, format = "%Y-%m-%d")

# Metadata: Add Year/Week for window
meta_rsvA$Collection_YearWeek <- paste(meta_rsvA$MMWRyear, 
                                       ifelse(meta_rsvA$MMWRweek < 10, paste0("0", meta_rsvA$MMWRweek), meta_rsvA$MMWRweek), 
                                       sep = "/")

meta_rsvB$Collection_YearWeek <- paste(meta_rsvB$MMWRyear, 
                                       ifelse(meta_rsvB$MMWRweek < 10, paste0("0", meta_rsvB$MMWRweek), meta_rsvB$MMWRweek), 
                                       sep = "/")
RefSeq_rsvA <- "NC_038235.1"
RefSeq_rsvB <- "NC_001781.1"

##########################################################################################

# All EU sequences

## German sequences
meta_rsvA_GER <- subset(meta_rsvA, Country == "Germany" | Accession == RefSeq_rsvA) #same as NCBI_rsvA_wgs_germany_2015.csv
meta_rsvB_GER <- subset(meta_rsvB, Country == "Germany"| Accession == RefSeq_rsvB) #same as NCBI_rsvB_wgs_germany_2015.csv

## Non-German EU sequences
meta_rsvA_EU <- subset(meta_rsvA, Country != "Germany"| Accession == RefSeq_rsvA)
meta_rsvB_EU <- subset(meta_rsvB, Country != "Germany"| Accession == RefSeq_rsvB)

# Read distance matrices 

# Pairwise
pairdist_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/pairdist_rsvA_EU.csv")
pairdist_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/pairdist_rsvB_EU.csv")

# SNP
snpdist_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_EU.csv")
snpdist_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_EU.csv")

# Hamming 

hamdist_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_protein.csv")
hamdist_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_protein.csv")

# Sort row and col 

#Pairwise
pairdist_rsvA <- arrange(pairdist_rsvA, pairdist_rsvA[,1])
pairdist_rsvA <- pairdist_rsvA[,-1]
pairdist_rsvA$Accession <- meta_rsvA$Accession #shorter name
rownames(pairdist_rsvA) <- pairdist_rsvA$Accession
pairdist_rsvA$Accession <- NULL
pairdist_rsvA <- pairdist_rsvA %>% select(order(colnames(.)))

pairdist_rsvB <- arrange(pairdist_rsvB, pairdist_rsvB[,1])
pairdist_rsvB <- pairdist_rsvB[,-1]
pairdist_rsvB$Accession <- meta_rsvB$Accession #shorter name
rownames(pairdist_rsvB) <- pairdist_rsvB$Accession
pairdist_rsvB$Accession <- NULL
pairdist_rsvB <- pairdist_rsvB %>% select(order(colnames(.)))

#SNP
snpdist_rsvA <- arrange(snpdist_rsvA, snpdist_rsvA[,1])
snpdist_rsvA <- snpdist_rsvA[,-1]
snpdist_rsvA$Accession <- meta_rsvA$Accession #shorter name
rownames(snpdist_rsvA) <- snpdist_rsvA$Accession
snpdist_rsvA$Accession <- NULL
snpdist_rsvA <- snpdist_rsvA %>% select(order(colnames(.)))

snpdist_rsvB <- arrange(snpdist_rsvB, snpdist_rsvB[,1])
snpdist_rsvB <- snpdist_rsvB[,-1]
snpdist_rsvB$Accession <- meta_rsvB$Accession #shorter name
rownames(snpdist_rsvB) <- snpdist_rsvB$Accession
snpdist_rsvB$Accession <- NULL
snpdist_rsvB <- snpdist_rsvB %>% select(order(colnames(.)))

#Hamming
hamdist_rsvA <- arrange(hamdist_rsvA, hamdist_rsvA[,1])
hamdist_rsvA <- hamdist_rsvA[,-1]
hamdist_rsvA$Accession <- meta_rsvA$Accession #shorter name
rownames(hamdist_rsvA) <- hamdist_rsvA$Accession
hamdist_rsvA$Accession <- NULL
hamdist_rsvA <- hamdist_rsvA %>% select(order(colnames(.)))

hamdist_rsvB <- arrange(hamdist_rsvB, hamdist_rsvB[,1])
hamdist_rsvB <- hamdist_rsvB[,-1]
hamdist_rsvB$Accession <- meta_rsvB$Accession #shorter name
rownames(hamdist_rsvB) <- hamdist_rsvB$Accession
hamdist_rsvB$Accession <- NULL
hamdist_rsvB <- hamdist_rsvB %>% select(order(colnames(.)))

# Distgroups

#Pairwise    
group_pairdist_rsvA <- dist_groups(pairdist_rsvA, meta_rsvA$Collection_YearWeek)
within_pairdist_rsvA <- subset(group_pairdist_rsvA, Group1 == Group2)

group_pairdist_rsvB <- dist_groups(pairdist_rsvB, meta_rsvB$Collection_YearWeek)
within_pairdist_rsvB <- subset(group_pairdist_rsvB, Group1 == Group2)

colnames(within_pairdist_rsvA)[1] <- "Accession"
within_pairdist_rsvA$Year <- meta_rsvA$MMWRyear[match(within_pairdist_rsvA$Accession, meta_rsvA$Accession)] #SVERWEIS/VLOOKUP
within_pairdist_rsvA$Week <- meta_rsvA$MMWRweek[match(within_pairdist_rsvA$Accession, meta_rsvA$Accession)]
within_pairdist_rsvA$Date <- decimal_date(MMWRweek2Date(within_pairdist_rsvA$Year, within_pairdist_rsvA$Week))

colnames(within_pairdist_rsvB)[1] <- "Accession"
within_pairdist_rsvB$Year <- meta_rsvB$MMWRyear[match(within_pairdist_rsvB$Accession, meta_rsvB$Accession)]
within_pairdist_rsvB$Week <- meta_rsvB$MMWRweek[match(within_pairdist_rsvB$Accession, meta_rsvB$Accession)]
within_pairdist_rsvB$Date <- decimal_date(MMWRweek2Date(within_pairdist_rsvB$Year, within_pairdist_rsvB$Week))

stats_pairdist_rsvA <- summarySE(data = within_pairdist_rsvA, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval
stats_pairdist_rsvB <- summarySE(data = within_pairdist_rsvB, measurevar = "Distance", groupvars = "Date")

#SNP
group_snpdist_rsvA <- dist_groups(snpdist_rsvA, meta_rsvA$Collection_YearWeek)
within_snpdist_rsvA <- subset(group_snpdist_rsvA, Group1 == Group2)

group_snpdist_rsvB <- dist_groups(snpdist_rsvB, meta_rsvB$Collection_YearWeek)
within_snpdist_rsvB <- subset(group_snpdist_rsvB, Group1 == Group2)

colnames(within_snpdist_rsvA)[1] <- "Accession"
within_snpdist_rsvA$Year <- meta_rsvA$MMWRyear[match(within_snpdist_rsvA$Accession, meta_rsvA$Accession)]
within_snpdist_rsvA$Week <- meta_rsvA$MMWRweek[match(within_snpdist_rsvA$Accession, meta_rsvA$Accession)]
within_snpdist_rsvA$Date <- decimal_date(MMWRweek2Date(within_snpdist_rsvA$Year, within_snpdist_rsvA$Week))

colnames(within_snpdist_rsvB)[1] <- "Accession"
within_snpdist_rsvB$Year <- meta_rsvB$MMWRyear[match(within_snpdist_rsvB$Accession, meta_rsvB$Accession)]
within_snpdist_rsvB$Week <- meta_rsvB$MMWRweek[match(within_snpdist_rsvB$Accession, meta_rsvB$Accession)]
within_snpdist_rsvB$Date <- decimal_date(MMWRweek2Date(within_snpdist_rsvB$Year, within_snpdist_rsvB$Week))

stats_snpdist_rsvA <- summarySE(data = within_snpdist_rsvA, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval
stats_snpdist_rsvB <- summarySE(data = within_snpdist_rsvB, measurevar = "Distance", groupvars = "Date")

#Ham
group_hamdist_rsvA <- dist_groups(hamdist_rsvA, meta_rsvA$Collection_YearWeek)
within_hamdist_rsvA <- subset(group_hamdist_rsvA, Group1 == Group2)

group_hamdist_rsvB <- dist_groups(hamdist_rsvB, meta_rsvB$Collection_YearWeek)
within_hamdist_rsvB <- subset(group_hamdist_rsvB, Group1 == Group2)

colnames(within_hamdist_rsvA)[1] <- "Accession"
within_hamdist_rsvA$Year <- meta_rsvA$MMWRyear[match(within_hamdist_rsvA$Accession, meta_rsvA$Accession)]
within_hamdist_rsvA$Week <- meta_rsvA$MMWRweek[match(within_hamdist_rsvA$Accession, meta_rsvA$Accession)]
within_hamdist_rsvA$Date <- decimal_date(MMWRweek2Date(within_hamdist_rsvA$Year, within_hamdist_rsvA$Week))

colnames(within_hamdist_rsvB)[1] <- "Accession"
within_hamdist_rsvB$Year <- meta_rsvB$MMWRyear[match(within_hamdist_rsvB$Accession, meta_rsvB$Accession)]
within_hamdist_rsvB$Week <- meta_rsvB$MMWRweek[match(within_hamdist_rsvB$Accession, meta_rsvB$Accession)]
within_hamdist_rsvB$Date <- decimal_date(MMWRweek2Date(within_hamdist_rsvB$Year, within_hamdist_rsvB$Week))

stats_hamdist_rsvA <- summarySE(data = within_hamdist_rsvA, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval
stats_hamdist_rsvB <- summarySE(data = within_hamdist_rsvB, measurevar = "Distance", groupvars = "Date")

# Plot

#Pairwise
ggplot(stats_pairdist_rsvA, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/Plots/weekly_rsvA_pair_EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)

ggplot(stats_pairdist_rsvB, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/weekly_rsvB_pair_EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)

#SNP
ggplot(stats_snpdist_rsvA, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/weekly_rsvA_snp_EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)

ggplot(stats_snpdist_rsvB, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/weekly_rsvB_snp_EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)

#Hamming
ggplot(stats_hamdist_rsvA, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/weekly_rsvA_ham_EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)

ggplot(stats_hamdist_rsvB, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#ggsave(filename = "~/Yale_Projects/Genetic_Diversity_RSV/weekly_rsvB_ham_EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)
