#Project: Change in genetic diversity of RSV isolates over the course of one pandemic season in Germany

library(dplyr)

NCBI_wgs_df = read.csv("~/RSV/Genetic Diversity Project/NCBI_rsvA_wgs_germany_2015.csv") #NCBI RSV-A data

NCBI_wgs_df$Collection_Year = as.numeric(format(as.Date(NCBI_wgs_df$Collection_Date, format = "%m/%d/%Y"),"%Y"))
NCBI_wgs_df$Collection_Month = as.numeric(format(as.Date(NCBI_wgs_df$Collection_Date, format = "%m/%d/%Y"),"%m"))
NCBI_wgs_df$Collection_Season = ifelse(NCBI_wgs_df$Collection_Month >= 7 & NCBI_wgs_df$Collection_Month <= 12, paste(NCBI_wgs_df$Collection_Year,"/",NCBI_wgs_df$Collection_Year+1, sep = ""), paste(NCBI_wgs_df$Collection_Year-1,"/",NCBI_wgs_df$Collection_Year, sep = ""))

numberofisolates = as.data.frame.matrix(table(NCBI_wgs_df$Collection_Season))
colnames(numberofisolates) = c("Season", "number_of_isolates")