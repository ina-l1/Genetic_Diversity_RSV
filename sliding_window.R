# Genetic Distance over Time: Sliding Window
# EU sequences

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
meta_rsvA <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/rsvA_ref_metadata_EU.csv") #already in alphabetical order
meta_rsvB <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/rsvB_ref_metadata_EU.csv")

# Read distance matrices 

# Pairwise
pairdist_rsvA <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/pairdist_rsvA_EU.csv")
pairdist_rsvB <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/pairdist_rsvB_EU.csv")

# SNP
snpdist_rsvA <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/snpdist_rsvA_EU.csv")
snpdist_rsvB <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/snpdist_rsvB_EU.csv")

# Metadata: Add Year/Week for window

meta_rsvA$Collection_YearWeek <- paste(meta_rsvA$MMWRyear, 
                                       ifelse(meta_rsvA$MMWRweek < 10, paste0("0", meta_rsvA$MMWRweek), meta_rsvA$MMWRweek), 
                                       sep = "/")

meta_rsvB$Collection_YearWeek <- paste(meta_rsvB$MMWRyear, 
                                       ifelse(meta_rsvB$MMWRweek < 10, paste0("0", meta_rsvB$MMWRweek), meta_rsvB$MMWRweek), 
                                       sep = "/")

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


# Distgroups

#Pairwise    
group_pairdist_rsvA <- dist_groups(pairdist_rsvA, meta_rsvA$Collection_YearWeek)
within_pairdist_rsvA <- subset(group_pairdist_rsvA, Group1 == Group2)

group_pairdist_rsvB <- dist_groups(pairdist_rsvB, meta_rsvB$Collection_YearWeek)
within_pairdist_rsvB <- subset(group_pairdist_rsvB, Group1 == Group2)

colnames(within_pairdist_rsvA)[1] <- "Accession"
within_pairdist_rsvA$Year <- meta_rsvA$MMWRyear[match(within_pairdist_rsvA$Accession, meta_rsvA$Accession)]
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

# Plot

#Pairwise
ggplot(stats_pairdist_rsvA, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

ggplot(stats_pairdist_rsvB, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

#SNP
ggplot(stats_snpdist_rsvA, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()

ggplot(stats_snpdist_rsvB, aes(x = Date, y = Distance)) +
  geom_errorbar(aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point()