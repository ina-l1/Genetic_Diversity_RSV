# Genetic distance to reference sequence

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(MMWRweek)

#################################################################################

# RSV-A and RSV-B: Replace rsvA/rsvB
# Diversity Measures: Replace snpdist/evodist/hamdist

rsvAB_choose <- "rsvB" ##

dist_rsvA_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvA.csv") ##
dist_rsvA_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_EU_noGer.csv") ##
dist_rsvB_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvB.csv") ##
dist_rsvB_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_EU_noGer.csv") ##

#################################################################################

# Read metadata

meta_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv") #already in alphabetical order
meta_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv")

meta_rsvA$Collection_Date <- as.Date(meta_rsvA$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
meta_rsvB$Collection_Date <- as.Date(meta_rsvB$Collection_Date, format = "%Y-%m-%d")

RefSeq_rsvA <- "NC_038235.1"
RefSeq_rsvB <- "NC_001781.1"

## German sequences
meta_rsvA_GER <- subset(meta_rsvA, Country == "Germany" | Accession == RefSeq_rsvA) #same as NCBI_rsvA_wgs_germany_2015.csv
meta_rsvB_GER <- subset(meta_rsvB, Country == "Germany" | Accession == RefSeq_rsvB) #same as NCBI_rsvB_wgs_germany_2015.csv

## Non-German EU sequences
meta_rsvA_EU <- subset(meta_rsvA, Country != "Germany" | Accession == RefSeq_rsvA)
meta_rsvB_EU <- subset(meta_rsvB, Country != "Germany" | Accession == RefSeq_rsvB)

#################################################################################

# Timetable (df) with weeks and years

# Start: 2014 W40
# End: 2023 W39
# 2015 and 2020 have 53 instead of 53 weeks

week <- rep(c(1:53), times = 10)
year <- rep(c(2014:2023), each = 53)
dates_df <- data.frame(year, week)

vec <- which(dates_df$week == 53 & (dates_df$year != 2015 & dates_df$year!=2020)) #years with no 53rd week
vec <- append(vec, which((between(dates_df$week, 1, 39) & dates_df$year == 2014)|(between(dates_df$week, 40, 53) & dates_df$year == 2023)))

dates_df <- dates_df[-vec,]

dates_df$date <- decimal_date(MMWRweek2Date(dates_df$year, dates_df$week))

dates_df$index <- 1:nrow(dates_df)
rownames(dates_df) <- dates_df$index

# Determining window size and start/end date

## SET WINDOW SIZE (weeks) ##
sliding_window_size <- 8 

sliding_window <- data.frame(matrix(ncol = 3))
colnames(sliding_window) <- c("index","start_date","end_date")
sliding_window <- sliding_window[-1,]

for(i in 1:(nrow(dates_df)-(sliding_window_size - 1))) { 
  sliding_window[i,"index"] <- i
  sliding_window[i,"start_date"] <- dates_df$date[i]
  sliding_window[i,"end_date"] <- dates_df$date[i+(sliding_window_size - 1)]
}

#################################################################################

# Distance matrix 

dist_matrix_func <- function(distance_matrix) {
  distance_matrix <- arrange(distance_matrix, distance_matrix[,1])
  distance_matrix <- distance_matrix[,-1]
  distance_matrix <- distance_matrix %>% select(order(colnames(.)))
}

dist_rsvA_GER <- dist_matrix_func(distance_matrix = dist_rsvA_GER)
dist_rsvA_EU <- dist_matrix_func(distance_matrix = dist_rsvA_EU)
dist_rsvB_GER <- dist_matrix_func(distance_matrix = dist_rsvB_GER)
dist_rsvB_EU <- dist_matrix_func(distance_matrix = dist_rsvB_EU)

rownames(dist_rsvA_GER) <- meta_rsvA_GER$Accession
colnames(dist_rsvA_GER) <- meta_rsvA_GER$Accession

rownames(dist_rsvA_EU) <- meta_rsvA_EU$Accession
colnames(dist_rsvA_EU) <- meta_rsvA_EU$Accession

rownames(dist_rsvB_GER) <- meta_rsvB_GER$Accession
colnames(dist_rsvB_GER) <- meta_rsvB_GER$Accession

rownames(dist_rsvB_EU) <- meta_rsvB_EU$Accession
colnames(dist_rsvB_EU) <- meta_rsvB_EU$Accession

if(rsvAB_choose == "rsvA") {
  dist_rsvAB_GER <- dist_rsvA_GER
  dist_rsvAB_EU <- dist_rsvA_EU
  RefSeq <- RefSeq_rsvA
  
} else {
  dist_rsvAB_GER <- dist_rsvB_GER
  dist_rsvAB_EU <- dist_rsvB_EU
  RefSeq <- RefSeq_rsvB
}

#################################################################################

# Distance to Ref

refdist_mean_func <- function(dist_matrix, EUGER) {
  if(EUGER == "GER" & rsvAB_choose == "rsvA"){
    meta <- meta_rsvA_GER
    dist_matrix <- dist_rsvAB_GER
  } else if(EUGER == "EU" & rsvAB_choose == "rsvA") {
    meta <- meta_rsvA_EU
    dist_matrix <- dist_rsvAB_EU
  } else if(EUGER == "GER" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_GER
    dist_matrix <- dist_rsvAB_GER
  } else if(EUGER == "EU" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_EU
    dist_matrix <- dist_rsvAB_EU
  } else {
    print("invalid arguments: Choose GER or EU.")
  }
  
  refdist_matrix <- select(dist_matrix, contains(RefSeq))
  refdist_matrix <- add_column(refdist_matrix, decimal_date(meta$Collection_Date))
  refdist_matrix <- refdist_matrix[-(which(is.na(refdist_matrix[,2]))),]
  colnames(refdist_matrix) <- c("Refdist", "Collection_Date")
  refdist_matrix$Window <- NA
  
  for(refdist_index in 1:nrow(refdist_matrix)) {
    for(window_index in 1:nrow(sliding_window)) {
      if(refdist_matrix[refdist_index, "Collection_Date"] >= sliding_window[window_index, "start_date"] &
         refdist_matrix[refdist_index, "Collection_Date"] <= sliding_window[window_index, "end_date"]) {
        refdist_matrix[refdist_index, "Window"] <- sliding_window[window_index, "index"]
      }
    }
  }
  
  refdist_win_df <- data.frame("Window" = sliding_window$index)
  refdist_df <- full_join(refdist_win_df, refdist_matrix, by = "Window")
  
  refdist_mean_df <- data.frame("Window" =  sliding_window$index, "start_date" = sliding_window$start_date, "mean" = NA, "std" = NA)
  for(win_index in 1:nrow(refdist_mean_df)) {
    refdist_mean_df$mean[win_index] <- mean(refdist_df$Refdist[which(refdist_df$Window == win_index)])
    refdist_mean_df$std[win_index] <- sd(refdist_df$Refdist[which(refdist_df$Window == win_index)])
  }
  
  return(refdist_mean_df)
}

## GER

refdist_GER <- refdist_mean_func(dist_rsvAB_GER, "GER")

'ggplot(refdist_GER, aes(x = start_date, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin = mean-std, ymax = mean+std))'

## EU

refdist_EU <- refdist_mean_func(dist_rsvAB_EU, "EU")

'ggplot(refdist_EU, aes(x = start_date, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin = mean-std, ymax = mean+std))'

## EU vs. GER

sliding_window_plot <- ggplot(refdist_GER, aes(x = start_date, y = mean)) +
  geom_errorbar(data = refdist_EU, aes(ymin = mean-std, ymax = mean+std), colour = "lightgrey") +
  geom_errorbar(aes(ymin = mean-std, ymax = mean+std), colour = "black", alpha = 0.45) +
  geom_point(data = refdist_EU, aes(x = start_date, y = mean), colour = "blue") +
  geom_point(colour = "red") +
  scale_x_continuous(breaks = c(2014:2023)) +
  xlab("Time") +
  ylab("Average pairwise distance to RefSeq") +
  theme_minimal()
sliding_window_plot
