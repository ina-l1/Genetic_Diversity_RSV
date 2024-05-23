#Project: Change in genetic diversity of RSV isolates over the course of one pandemic season in Germany
#Data Sorting 

library(dplyr)
library(MMWRweek)

# German sequences 

'NCBI_wgs_df = read.csv("~/RSV/Genetic Diversity Project/NCBI_rsvA_wgs_germany_2015.csv") #NCBI RSV-A data

NCBI_wgs_df$Collection_Year = as.numeric(format(as.Date(NCBI_wgs_df$Collection_Date, format = "%m/%d/%Y"),"%Y"))
NCBI_wgs_df$Collection_Month = as.numeric(format(as.Date(NCBI_wgs_df$Collection_Date, format = "%m/%d/%Y"),"%m"))
NCBI_wgs_df$Collection_Season = ifelse(NCBI_wgs_df$Collection_Month >= 7 & NCBI_wgs_df$Collection_Month <= 12, paste(NCBI_wgs_df$Collection_Year,"/",NCBI_wgs_df$Collection_Year+1, sep = ""), paste(NCBI_wgs_df$Collection_Year-1,"/",NCBI_wgs_df$Collection_Year, sep = ""))

numberofisolates = as.data.frame.matrix(table(NCBI_wgs_df$Collection_Season))
colnames(numberofisolates) = c("Season", "number_of_isolates")'

# European sequences

rsvA_meta_EU <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/NCBI_rsvA_wgs_europe_2015.csv")
rsvB_meta_EU <- read.csv("~/RSV/git/RSV Genetic Diversity/Europe/NCBI_rsvB_wgs_europe_2015.csv")

rsvA_meta_EU$Collection_Date <- as.Date(rsvA_meta_EU$Collection_Date, format = "%m/%d/%Y")
rsvB_meta_EU$Collection_Date <- as.Date(rsvB_meta_EU$Collection_Date, format = "%m/%d/%Y")

rsvA_meta_EU <- subset(rsvA_meta_EU, is.na(Collection_Date) == FALSE)
rsvB_meta_EU <- subset(rsvB_meta_EU, is.na(Collection_Date) == FALSE)

# --> Number of sequences: RSV-A 731, RSV-B 699

write.csv(rsvA_meta_EU, file = "~/RSV/git/RSV Genetic Diversity/Europe/NCBI_rsvA_wgs_europe_2015.csv", row.names = FALSE)
write.csv(rsvB_meta_EU, file = "~/RSV/git/RSV Genetic Diversity/Europe/NCBI_rsvB_wgs_europe_2015.csv", row.names = FALSE)

#Collection Season

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

#Number of sequences for each country and season 

rsvA_seqNum_EU <- as.data.frame(table(rsvA_meta_EU$Country))
rsvB_seqNum_EU <- as.data.frame(table(rsvB_meta_EU$Country))

rsvA_seqNum_season_EU <- as.data.frame.matrix(table(rsvA_meta_EU$Country, rsvA_meta_EU$Collection_Season))
rsvB_seqNum_season_EU <- as.data.frame.matrix(table(rsvB_meta_EU$Country, rsvB_meta_EU$Collection_Season))