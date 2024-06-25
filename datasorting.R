# Project: Change in genetic diversity of RSV isolates over the course of one pandemic season in Germany
# Data Sorting 

library(dplyr)
library(tidyr)
library(MMWRweek)

##################################################

# Reference sequences

RefSeq_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/RefSeq/RefSeq_rsvA.csv")
RefSeq_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/RefSeq/RefSeq_rsvB.csv")

##################################################

# German sequences 

# Sequence information for NCBI WGS sequences starting from 2014/2015 season 
rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/NCBI_rsvA_wgs_europe_2015.csv")
rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/NCBI_rsvB_wgs_europe_2015.csv")

rsvA_nextclade <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/nextclade_rsvA_EU.csv") #Nextclade output
rsvB_nextclade <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/nextclade_rsvB_EU.csv")

# Reformat Dates
rsvA$Collection_Date <- as.Date(rsvA$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
rsvB$Collection_Date <- as.Date(rsvB$Collection_Date, format = "%Y-%m-%d")

rsvA$Collection_Year <- as.numeric(format(rsvA$Collection_Date,"%Y"))
rsvA$Collection_Month <- as.numeric(format(rsvA$Collection_Date,"%m"))

rsvB$Collection_Year <- as.numeric(format(rsvB$Collection_Date, "%Y"))
rsvB$Collection_Month <- as.numeric(format(rsvB$Collection_Date, "%m"))

# RSV epidemic season in Germany from Week 40 to Week 39 of next year (October)
rsvAMMWR <- MMWRweek(rsvA$Collection_Date)
rsvBMMWR <- MMWRweek(rsvB$Collection_Date)
rsvA$MMWRyear <- rsvAMMWR$MMWRyear
rsvA$MMWRweek <- rsvAMMWR$MMWRweek
rsvA$MMWRday <- rsvAMMWR$MMWRday
rsvB$MMWRyear <- rsvBMMWR$MMWRyear
rsvB$MMWRweek <- rsvBMMWR$MMWRweek
rsvB$MMWRday <- rsvBMMWR$MMWRday

rsvA$Collection_Season <- ifelse(rsvA$MMWRweek >= 1 & rsvA$MMWRweek <= 39, paste(rsvA$MMWRyear-1,"/",rsvA$MMWRyear, sep = ""), paste(rsvA$MMWRyear,"/",rsvA$MMWRyear+1, sep = ""))
rsvB$Collection_Season <- ifelse(rsvB$MMWRweek >= 1 & rsvB$MMWRweek <= 39, paste(rsvB$MMWRyear-1,"/",rsvB$MMWRyear, sep = ""), paste(rsvB$MMWRyear,"/",rsvB$MMWRyear+1, sep = ""))

# A and B in one table
rsvA$Type <- "A"
rsvB$Type <- "B"
rsvAB <- rbind(rsvA, rsvB)

# rsvA + RefSeqA, rsvB + RefSeqB
rsvA_ref <- full_join(rsvA, RefSeq_rsvA, by = join_by(Accession, Organism_Name, Organization, Org_location, Isolate,
                                               Molecule_type, Length, Nuc_Completeness, Country, Collection_Date, GenBank_Title))
rsvB_ref <- full_join(rsvB, RefSeq_rsvB, by = join_by(Accession, Organism_Name, Organization, Org_location, Isolate,
                                               Molecule_type, Length, Nuc_Completeness, Country, Collection_Date, GenBank_Title))

# Sort dataframes
rsvA <- rsvA[order(rsvA$Accession),]
rsvB <- rsvB[order(rsvB$Accession),]
rsvAB <- rsvAB[order(rsvAB$Accession),]
rsvA_ref <- rsvA_ref[order(rsvA_ref$Accession),]
rsvB_ref <- rsvB_ref[order(rsvB_ref$Accession),]

# Add nextclade clade information to metadata
rsvA_nextclade <- arrange(rsvA_nextclade, seqName)
rsvB_nextclade <- arrange(rsvB_nextclade, seqName)
rsvA_clades <- rsvA_nextclade[c("clade", "G_clade")]
rsvB_clades <- rsvB_nextclade[c("clade", "G_clade")]
rsvA_ref <- cbind(rsvA_ref, rsvA_clades)
rsvB_ref <- cbind(rsvB_ref, rsvB_clades)

# Export metadata
#write.csv(rsvA_ref, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata.csv", row.names = FALSE)
#write.csv(rsvB_ref, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata.csv", row.names = FALSE)

#write.csv(rsvAB, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvAB_metadata_EU.csv")

##########################################################

# European sequences

rsvA_meta_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/NCBI_rsvA_wgs_europe_2015.csv")
rsvB_meta_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/NCBI_rsvB_wgs_europe_2015.csv")

rsvA_meta_EU$Collection_Date <- as.Date(rsvA_meta_EU$Collection_Date, format = "%Y-%m-%d") #WARNING: Format changes if opened in Excel before
rsvB_meta_EU$Collection_Date <- as.Date(rsvB_meta_EU$Collection_Date, format = "%Y-%m-%d")

#rsvA_meta_EU <- subset(rsvA_meta_EU, is.na(Collection_Date) == FALSE)
#rsvB_meta_EU <- subset(rsvB_meta_EU, is.na(Collection_Date) == FALSE)

#write.csv(rsvA_meta_EU, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/NCBI_rsvA_wgs_europe_2015.csv", row.names = FALSE)
#write.csv(rsvB_meta_EU, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/NCBI_rsvB_wgs_europe_2015.csv", row.names = FALSE)

rsvA_nextclade_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/nextclade_rsvA_EU.csv")
rsvB_nextclade_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/nextclade_rsvB_EU.csv")

# Generating metadata

MMWR_rsvA_meta_EU <- MMWRweek(rsvA_meta_EU$Collection_Date)
MMWR_rsvB_meta_EU <- MMWRweek(rsvB_meta_EU$Collection_Date)

rsvA_meta_EU$MMWRyear <- MMWR_rsvA_meta_EU$MMWRyear
rsvA_meta_EU$MMWRweek <- MMWR_rsvA_meta_EU$MMWRweek
rsvB_meta_EU$MMWRyear <- MMWR_rsvB_meta_EU$MMWRyear
rsvB_meta_EU$MMWRweek <- MMWR_rsvB_meta_EU$MMWRweek

rsvA_meta_EU$Collection_Season <- ifelse(rsvA_meta_EU$MMWRweek >= 1 & rsvA_meta_EU$MMWRweek <= 39, 
                                         paste(rsvA_meta_EU$MMWRyear-1,"/",rsvA_meta_EU$MMWRyear, sep = ""), 
                                         paste(rsvA_meta_EU$MMWRyear,"/",rsvA_meta_EU$MMWRyear+1, sep = ""))
rsvB_meta_EU$Collection_Season <- ifelse(rsvB_meta_EU$MMWRweek >= 1 & rsvB_meta_EU$MMWRweek <= 39, 
                                         paste(rsvB_meta_EU$MMWRyear-1,"/",rsvB_meta_EU$MMWRyear, sep = ""), 
                                         paste(rsvB_meta_EU$MMWRyear,"/",rsvB_meta_EU$MMWRyear+1, sep = ""))

rsvA_meta_EU$Type <- "A"
RefSeq_rsvA$Type <- "A"
rsvB_meta_EU$Type <- "B"
RefSeq_rsvB$Type <- "B"

rsvA_ref_meta_EU <- full_join(rsvA_meta_EU, RefSeq_rsvA, by = join_by(Accession, Organism_Name, Organization, Org_location, Isolate,
                                                                      Molecule_type, Length, Nuc_Completeness, Country, Collection_Date, GenBank_Title, Type))
rsvB_ref_meta_EU <- full_join(rsvB_meta_EU, RefSeq_rsvB, by = join_by(Accession, Organism_Name, Organization, Org_location, Isolate,
                                                                      Molecule_type, Length, Nuc_Completeness, Country, Collection_Date, GenBank_Title, Type))

rsvA_ref_meta_EU <- arrange(rsvA_ref_meta_EU, Accession)
rsvB_ref_meta_EU <- arrange(rsvB_ref_meta_EU, Accession)
rsvA_nextclade_EU <- arrange(rsvA_nextclade_EU, seqName)
rsvB_nextclade_EU <- arrange(rsvB_nextclade_EU, seqName)

rsvA_ref_meta_EU$clade <- rsvA_nextclade_EU$clade
rsvA_ref_meta_EU$G_clade <- rsvA_nextclade_EU$G_clade
rsvB_ref_meta_EU$clade <- rsvB_nextclade_EU$clade
rsvB_ref_meta_EU$G_clade <- rsvB_nextclade_EU$G_clade

#write.csv(rsvA_ref_meta_EU, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv", row.names = FALSE)
#write.csv(rsvB_ref_meta_EU, file = "~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv", row.names = FALSE)

#Number of sequences for each country and season 

rsvA_seqNum_EU <- as.data.frame(table(rsvA_meta_EU$Country))
rsvB_seqNum_EU <- as.data.frame(table(rsvB_meta_EU$Country))

rsvA_seqNum_season_EU <- as.data.frame.matrix(table(rsvA_meta_EU$Country, rsvA_meta_EU$Collection_Season))
rsvB_seqNum_season_EU <- as.data.frame.matrix(table(rsvB_meta_EU$Country, rsvB_meta_EU$Collection_Season))