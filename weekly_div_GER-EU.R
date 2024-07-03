# Genetic Distance over Time: Weekly Data
# German vs. EU sequences

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
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

meta_rsvA$Collection_YearWeek <- paste(meta_rsvA$MMWRyear, 
                                       ifelse(meta_rsvA$MMWRweek < 10, paste0("0", meta_rsvA$MMWRweek), meta_rsvA$MMWRweek), 
                                       sep = "/")

meta_rsvB$Collection_YearWeek <- paste(meta_rsvB$MMWRyear, 
                                       ifelse(meta_rsvB$MMWRweek < 10, paste0("0", meta_rsvB$MMWRweek), meta_rsvB$MMWRweek), 
                                       sep = "/")
RefSeq_rsvA <- "NC_038235.1"
RefSeq_rsvB <- "NC_001781.1"

## German sequences
meta_rsvA_GER <- subset(meta_rsvA, Country == "Germany" | Accession == RefSeq_rsvA) #same as NCBI_rsvA_wgs_germany_2015.csv
meta_rsvB_GER <- subset(meta_rsvB, Country == "Germany"| Accession == RefSeq_rsvB) #same as NCBI_rsvB_wgs_germany_2015.csv

## Non-German EU sequences
meta_rsvA_EU <- subset(meta_rsvA, Country != "Germany"| Accession == RefSeq_rsvA)
meta_rsvB_EU <- subset(meta_rsvB, Country != "Germany"| Accession == RefSeq_rsvB)

# Distance matrix 
## SNP
### RSV A
snpdist_rsvA_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvA.csv")
snpdist_rsvA_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_EU_noGer.csv")

snpdist_rsvA_GER <- arrange(snpdist_rsvA_GER, snpdist_rsvA_GER[,1])
snpdist_rsvA_GER <- snpdist_rsvA_GER[,-1]
rownames(snpdist_rsvA_GER) <- meta_rsvA_GER$Accession
snpdist_rsvA_GER <- snpdist_rsvA_GER %>% select(order(colnames(.)))
colnames(snpdist_rsvA_GER) <- meta_rsvA_GER$Accession

snpdist_rsvA_EU <- arrange(snpdist_rsvA_EU, snpdist_rsvA_EU[,1])
snpdist_rsvA_EU <- snpdist_rsvA_EU[,-1]
rownames(snpdist_rsvA_EU) <- meta_rsvA_EU$Accession
snpdist_rsvA_EU <- snpdist_rsvA_EU %>% select(order(colnames(.)))
colnames(snpdist_rsvA_EU) <- meta_rsvA_EU$Accession

### RSV B
snpdist_rsvB_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvB.csv")
snpdist_rsvB_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_EU_noGer.csv")

snpdist_rsvB_GER <- arrange(snpdist_rsvB_GER, snpdist_rsvB_GER[,1])
snpdist_rsvB_GER <- snpdist_rsvB_GER[,-1]
rownames(snpdist_rsvB_GER) <- meta_rsvB_GER$Accession
snpdist_rsvB_GER <- snpdist_rsvB_GER %>% select(order(colnames(.)))
colnames(snpdist_rsvB_GER) <- meta_rsvB_GER$Accession

snpdist_rsvB_EU <- arrange(snpdist_rsvB_EU, snpdist_rsvB_EU[,1])
snpdist_rsvB_EU <- snpdist_rsvB_EU[,-1]
rownames(snpdist_rsvB_EU) <- meta_rsvB_EU$Accession
snpdist_rsvB_EU <- snpdist_rsvB_EU %>% select(order(colnames(.)))
colnames(snpdist_rsvB_EU) <- meta_rsvB_EU$Accession

# Distgroups 
## SNP
### RSV A
group_snpdist_rsvA_GER <- dist_groups(snpdist_rsvA_GER, meta_rsvA_GER$Collection_YearWeek)
within_snpdist_rsvA_GER <- subset(group_snpdist_rsvA_GER, Group1 == Group2)
between_snpdist_rsvA_GER <- subset(group_snpdist_rsvA_GER, Group1 != Group2)
#between_snpdist_rsvA_GER$Between_Season <- paste(between_snpdist_rsvA_GER$Group1, "&", between_snpdist_rsvA_GER$Group2, sep = "")

colnames(within_snpdist_rsvA_GER)[1] <- "Accession"
within_snpdist_rsvA_GER$Year <- meta_rsvA_GER$MMWRyear[match(within_snpdist_rsvA_GER$Accession, meta_rsvA_GER$Accession)]
within_snpdist_rsvA_GER$Week <- meta_rsvA_GER$MMWRweek[match(within_snpdist_rsvA_GER$Accession, meta_rsvA_GER$Accession)]
within_snpdist_rsvA_GER$Date <- decimal_date(MMWRweek2Date(within_snpdist_rsvA_GER$Year, within_snpdist_rsvA_GER$Week))

stats_snpdist_rsvA_GER <- summarySE(data = within_snpdist_rsvA_GER, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval

group_snpdist_rsvA_EU <- dist_groups(snpdist_rsvA_EU, meta_rsvA_EU$Collection_YearWeek)
within_snpdist_rsvA_EU <- subset(group_snpdist_rsvA_EU, Group1 == Group2)
between_snpdist_rsvA_EU <- subset(group_snpdist_rsvA_EU, Group1 != Group2)

colnames(within_snpdist_rsvA_EU)[1] <- "Accession"
within_snpdist_rsvA_EU$Year <- meta_rsvA_EU$MMWRyear[match(within_snpdist_rsvA_EU$Accession, meta_rsvA_EU$Accession)]
within_snpdist_rsvA_EU$Week <- meta_rsvA_EU$MMWRweek[match(within_snpdist_rsvA_EU$Accession, meta_rsvA_EU$Accession)]
within_snpdist_rsvA_EU$Date <- decimal_date(MMWRweek2Date(within_snpdist_rsvA_EU$Year, within_snpdist_rsvA_EU$Week))

stats_snpdist_rsvA_EU <- summarySE(data = within_snpdist_rsvA_EU, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval

### RSV B
group_snpdist_rsvB_GER <- dist_groups(snpdist_rsvB_GER, meta_rsvB_GER$Collection_YearWeek)
within_snpdist_rsvB_GER <- subset(group_snpdist_rsvB_GER, Group1 == Group2)
between_snpdist_rsvB_GER <- subset(group_snpdist_rsvB_GER, Group1 != Group2)

colnames(within_snpdist_rsvB_GER)[1] <- "Accession"
within_snpdist_rsvB_GER$Year <- meta_rsvB_GER$MMWRyear[match(within_snpdist_rsvB_GER$Accession, meta_rsvB_GER$Accession)]
within_snpdist_rsvB_GER$Week <- meta_rsvB_GER$MMWRweek[match(within_snpdist_rsvB_GER$Accession, meta_rsvB_GER$Accession)]
within_snpdist_rsvB_GER$Date <- decimal_date(MMWRweek2Date(within_snpdist_rsvB_GER$Year, within_snpdist_rsvB_GER$Week))

stats_snpdist_rsvB_GER <- summarySE(data = within_snpdist_rsvB_GER, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval

group_snpdist_rsvB_EU <- dist_groups(snpdist_rsvB_EU, meta_rsvB_EU$Collection_YearWeek)
within_snpdist_rsvB_EU <- subset(group_snpdist_rsvB_EU, Group1 == Group2)
between_snpdist_rsvB_EU <- subset(group_snpdist_rsvB_EU, Group1 != Group2)

colnames(within_snpdist_rsvB_EU)[1] <- "Accession"
within_snpdist_rsvB_EU$Year <- meta_rsvB_EU$MMWRyear[match(within_snpdist_rsvB_EU$Accession, meta_rsvB_EU$Accession)]
within_snpdist_rsvB_EU$Week <- meta_rsvB_EU$MMWRweek[match(within_snpdist_rsvB_EU$Accession, meta_rsvB_EU$Accession)]
within_snpdist_rsvB_EU$Date <- decimal_date(MMWRweek2Date(within_snpdist_rsvB_EU$Year, within_snpdist_rsvB_EU$Week))

stats_snpdist_rsvB_EU <- summarySE(data = within_snpdist_rsvB_EU, measurevar = "Distance", groupvars = "Date") #mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval

# Plots
## SNP
### RSV A
plot_snpdist_rsvA <- ggplot(stats_snpdist_rsvA_GER, aes(x = Date, y = Distance)) +
  geom_errorbar(color = "grey", aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point(color = "red") +
  geom_errorbar(color = "grey", data = stats_snpdist_rsvA_EU, aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point(data = stats_snpdist_rsvA_EU, color = "blue")
plot_snpdist_rsvA
#ggsave(filename = "~/RSV/git/Plots/weekly_rsvA_snp_GER-EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)

### RSV B
plot_snpdist_rsvB <- ggplot(stats_snpdist_rsvB_GER, aes(x = Date, y = Distance)) +
  geom_errorbar(color = "grey", aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point(color = "red") +
  geom_errorbar(color = "grey", data = stats_snpdist_rsvB_EU, aes(ymin = Distance-se, ymax = Distance+se)) +
  geom_point(data = stats_snpdist_rsvB_EU, color = "blue")
plot_snpdist_rsvB
#ggsave(filename = "~/RSV/git/Plots/weekly_rsvB_snp_GER-EU.png", width = 50, height = 20, units = "cm", limitsize = FALSE)