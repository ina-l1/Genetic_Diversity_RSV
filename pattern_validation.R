# Validating temporal genetic diversity pattern of German sequences
# Countries: Netherlands, Spain, UK

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(MMWRweek)

#################################################################################

rsvAB_choose <- "rsvB" ##

# Read metadata

meta_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv") #already in alphabetical order
meta_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv")

meta_rsvA$Collection_Date <- as.Date(meta_rsvA$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
meta_rsvB$Collection_Date <- as.Date(meta_rsvB$Collection_Date, format = "%Y-%m-%d")

RefSeq_rsvA <- "NC_038235.1"
RefSeq_rsvB <- "NC_001781.1"

## Germany (GER)
meta_rsvA_GER <- subset(meta_rsvA, Country == "Germany" | Accession == RefSeq_rsvA) #same as NCBI_rsvA_wgs_germany_2015.csv
meta_rsvB_GER <- subset(meta_rsvB, Country == "Germany" | Accession == RefSeq_rsvB) #same as NCBI_rsvB_wgs_germany_2015.csv

## Netherlands (NLD)
meta_rsvA_NLD <- subset(meta_rsvA, Country == "Netherlands" | Accession == RefSeq_rsvA)
meta_rsvB_NLD <- subset(meta_rsvB, Country == "Netherlands" | Accession == RefSeq_rsvB)

## Spain (ESP)
meta_rsvA_ESP <- subset(meta_rsvA, Country == "Spain" | Accession == RefSeq_rsvA)
meta_rsvB_ESP <- subset(meta_rsvB, Country == "Spain" | Accession == RefSeq_rsvB)

## UK (GBR)
meta_rsvA_GBR <- subset(meta_rsvA, Country == "United Kingdom" | Accession == RefSeq_rsvA)
meta_rsvB_GBR <- subset(meta_rsvB, Country == "United Kingdom" | Accession == RefSeq_rsvB)

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

dist_rsvA_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/evodist_rsvA.csv") ##
dist_rsvA_NLD <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_netherlands.csv") ##
dist_rsvA_ESP <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_spain.csv") ##
dist_rsvA_GBR <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_UK.csv") ##

dist_rsvB_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/evodist_rsvB.csv") ##
dist_rsvB_NLD <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_netherlands.csv") ##
dist_rsvB_ESP <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_spain.csv") ##
dist_rsvB_GBR <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_UK.csv") ##

dist_matrix_func <- function(distance_matrix) {
  distance_matrix <- arrange(distance_matrix, distance_matrix[,1])
  distance_matrix <- distance_matrix[,-1]
  distance_matrix <- distance_matrix %>% select(order(colnames(.)))
}

dist_rsvA_GER <- dist_matrix_func(distance_matrix = dist_rsvA_GER)
dist_rsvB_GER <- dist_matrix_func(distance_matrix = dist_rsvB_GER)

dist_rsvA_NLD <- dist_matrix_func(distance_matrix = dist_rsvA_NLD)
dist_rsvB_NLD <- dist_matrix_func(distance_matrix = dist_rsvB_NLD)

dist_rsvA_ESP <- dist_matrix_func(distance_matrix = dist_rsvA_ESP)
dist_rsvB_ESP <- dist_matrix_func(distance_matrix = dist_rsvB_ESP)

dist_rsvA_GBR <- dist_matrix_func(distance_matrix = dist_rsvA_GBR)
dist_rsvB_GBR <- dist_matrix_func(distance_matrix = dist_rsvB_GBR)

rownames(dist_rsvA_GER) <- meta_rsvA_GER$Accession
colnames(dist_rsvA_GER) <- meta_rsvA_GER$Accession

rownames(dist_rsvA_NLD) <- meta_rsvA_NLD$Accession
colnames(dist_rsvA_NLD) <- meta_rsvA_NLD$Accession

rownames(dist_rsvA_ESP) <- meta_rsvA_ESP$Accession
colnames(dist_rsvA_ESP) <- meta_rsvA_ESP$Accession

rownames(dist_rsvA_GBR) <- meta_rsvA_GBR$Accession
colnames(dist_rsvA_GBR) <- meta_rsvA_GBR$Accession

rownames(dist_rsvB_GER) <- meta_rsvB_GER$Accession
colnames(dist_rsvB_GER) <- meta_rsvB_GER$Accession

rownames(dist_rsvB_NLD) <- meta_rsvB_NLD$Accession
colnames(dist_rsvB_NLD) <- meta_rsvB_NLD$Accession

rownames(dist_rsvB_ESP) <- meta_rsvB_ESP$Accession
colnames(dist_rsvB_ESP) <- meta_rsvB_ESP$Accession

rownames(dist_rsvB_GBR) <- meta_rsvB_GBR$Accession
colnames(dist_rsvB_GBR) <- meta_rsvB_GBR$Accession

#################################################################################

if(rsvAB_choose == "rsvA") {
  dist_rsvAB_GER <- dist_rsvA_GER
  dist_rsvAB_GBR <- dist_rsvA_GBR
  RefSeq <- RefSeq_rsvA
  
} else {
  dist_rsvAB_GER <- dist_rsvB_GER
  dist_rsvAB_GBR <- dist_rsvB_GBR
  RefSeq <- RefSeq_rsvB
}

# Sliding window function

sliding_window_func <- function(dist_mean_func_df, GBRGER) {
  
  if(GBRGER == "GER" & rsvAB_choose == "rsvA"){
    meta <- meta_rsvA_GER
    dist_matrix <- dist_rsvAB_GER
  } else if(GBRGER == "GBR" & rsvAB_choose == "rsvA") {
    meta <- meta_rsvA_GBR
    dist_matrix <- dist_rsvAB_GBR
  } else if(GBRGER == "GER" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_GER
    dist_matrix <- dist_rsvAB_GER
  } else if(GBRGER == "GBR" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_GBR
    dist_matrix <- dist_rsvAB_GBR
  } else {
    print("invalid arguments: Choose GER or GBR.")
  }
  
  for(window_index in 1:nrow(sliding_window)) {
    dist_sum <- 0
    dist_vec <- c()
    acc_list <- c()
    
    for(acc_index0 in 1:nrow(meta)) {
      if(meta$Accession[acc_index0] != RefSeq & 
         decimal_date(meta$Collection_Date[acc_index0]) >= sliding_window[window_index, "start_date"] & 
         decimal_date(meta$Collection_Date[acc_index0]) <= sliding_window[window_index, "end_date"]) {
        acc_list <- append(acc_list, meta$Accession[acc_index0])
      }
    }
    print(window_index)
    print(acc_list)
    
    acc_length <- length(acc_list)
    if(acc_length <= 1) {
      dist_mean_func_df[window_index, "dist_mean"] <- NA
    } else {
      for(acc_index1 in 1:acc_length) {
        for(acc_index2 in acc_index1:acc_length) {
          dist_sum <- dist_matrix[acc_list[acc_index1], 
                                  acc_list[acc_index2]] + dist_sum
          if(acc_index1 != acc_index2)  {
            dist_vec <- append(dist_vec, dist_matrix[acc_list[acc_index1], acc_list[acc_index2]])
          }
        }
      }
      print(dist_vec)
      
      number_pairwise_dist <- (acc_length*(acc_length-1))/2
      mean <- dist_sum/number_pairwise_dist
      dist_mean_func_df[window_index, "dist_mean"] <- mean
      
      dist_mean_func_df[window_index, "dist_std"] <- sd(unlist(dist_vec))
      dist_mean_func_df[window_index, "num"] <- acc_length
    }
  }
  return(dist_mean_func_df)
}

#################################################################################

# Distances
dist_mean_df <- data.frame("window" =  1:nrow(sliding_window), "num" = 0, "dist_mean" = 0, "dist_std" = 0, "start_date" = sliding_window$start_date)

## GER
dist_mean_GER <- sliding_window_func(dist_mean_func_df = dist_mean_df, GBRGER = "GER")

## GBR
dist_mean_GBR <- sliding_window_func(dist_mean_func_df = dist_mean_df, GBRGER = "GBR")

#################################################################################

## Evo
sliding_window_plot <- ggplot(dist_mean_GER, aes(x = start_date, y = dist_mean)) +
  geom_errorbar(data = dist_mean_GBR, aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "lightgrey") +
  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "black", alpha = 0.45) +
  geom_point(colour = "red") +
  geom_point(data = dist_mean_GBR, aes(x = start_date, y = dist_mean), colour = "darkgreen") +
  geom_line(aes(y = num/20000), colour = "red") +
  geom_line(aes(y = dist_mean_GBR$num/20000), colour = "darkgreen") +
  geom_line(aes(y = (num + dist_mean_GBR$num)/20000), colour = "black", alpha = 0.4) +
  scale_x_continuous(breaks = c(2014:2023)) +
  theme_minimal() +
  scale_y_continuous("distance (mean)", sec.axis = sec_axis(~.*20000, name="sequence count"), limits = c(0, 0.02)) # IF NECESSARY ADJUST Y-AXIS
sliding_window_plot

## SNP
sliding_window_plot <- ggplot(dist_mean_GER, aes(x = start_date, y = dist_mean)) +
  geom_errorbar(data = dist_mean_GBR, aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "lightgrey") +
  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "black", alpha = 0.45) +
  geom_point(colour = "red") +
  geom_point(data = dist_mean_GBR, aes(x = start_date, y = dist_mean), colour = "blue") +
  geom_line(aes(y = num/1.33), colour = "red") +
  geom_line(aes(y = dist_mean_GBR$num/1.33), colour = "blue") +
  geom_line(aes(y = (num + dist_mean_GBR$num)/1.33), colour = "black", alpha = 0.4) +
  scale_x_continuous(breaks = c(2014:2023)) +
  theme_minimal() +
  scale_y_continuous("distance (mean)", sec.axis = sec_axis(~.*1.33, name="sequence count"), limits = c(0, 300)) # IF NECESSARY ADJUST Y-AXIS
sliding_window_plot
