# Sliding Window
# GER vs. EU sequences

# RSV-A and RSV-B
# SNP, pairwise and Hamming distance

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

##'*SET WINDOW SIZE (weeks)*
sliding_window_size <- 8 

sliding_window <- data.frame(matrix(ncol = 3))
colnames(sliding_window) <- c("index","start_date","end_date")
sliding_window <- sliding_window[-1,]

for(i in 1:(nrow(dates_df)-(sliding_window_size - 1))) { 
  sliding_window[i,"index"] <- i
  sliding_window[i,"start_date"] <- dates_df$date[i]
  sliding_window[i,"end_date"] <- dates_df$date[i+(sliding_window_size - 1)]
}

# Distance matrix 

# RSV-A and RSV-B: Replace rsvA/rsvB
# Diversity Measures: Replace snpdist/pairdist/hamdist

dist_rsvB_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvB.csv")
dist_rsvB_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_EU_noGer.csv")

dist_rsvB_GER <- arrange(dist_rsvB_GER, dist_rsvB_GER[,1])
dist_rsvB_GER <- dist_rsvB_GER[,-1]
rownames(dist_rsvB_GER) <- meta_rsvB_GER$Accession
dist_rsvB_GER <- dist_rsvB_GER %>% select(order(colnames(.)))
colnames(dist_rsvB_GER) <- meta_rsvB_GER$Accession

dist_rsvB_EU <- arrange(dist_rsvB_EU, dist_rsvB_EU[,1])
dist_rsvB_EU <- dist_rsvB_EU[,-1]
rownames(dist_rsvB_EU) <- meta_rsvB_EU$Accession
dist_rsvB_EU <- dist_rsvB_EU %>% select(order(colnames(.)))
colnames(dist_rsvB_EU) <- meta_rsvB_EU$Accession

# Distances

## GER
dist_mean_GER <- data.frame("window" =  1:nrow(sliding_window), "num" = 0, "dist_mean" = 0, "dist_std" = 0, "start_date" = sliding_window$start_date)
violin_rsvB_GER_df <- data.frame(matrix(ncol = 3, nrow = 0))

for(window_index in 1:nrow(sliding_window)) {
  dist_sum <- 0
  dist_vec <- c()
  acc_list <- c()
  
  for(acc_index0 in 1:nrow(meta_rsvB_GER)) {
    if(meta_rsvB_GER$Accession[acc_index0] != RefSeq_rsvB & 
       decimal_date(meta_rsvB_GER$Collection_Date[acc_index0]) >= sliding_window[window_index, "start_date"] & 
       decimal_date(meta_rsvB_GER$Collection_Date[acc_index0]) <= sliding_window[window_index, "end_date"]) {
      acc_list <- append(acc_list, meta_rsvB_GER$Accession[acc_index0])
    }
  }
  print(window_index)
  print(acc_list)
  
  acc_length <- length(acc_list)
  if(acc_length <= 1) {
    dist_mean_GER[window_index, "dist_mean"] <- NA
  } else {
    for(acc_index1 in 1:acc_length) {
      for(acc_index2 in acc_index1:acc_length) {
        dist_sum <- dist_rsvB_GER[acc_list[acc_index1], 
                                      acc_list[acc_index2]] + dist_sum
        if(acc_index1 != acc_index2)  {
          dist_vec <- append(dist_vec, dist_rsvB_GER[acc_list[acc_index1], acc_list[acc_index2]])
        }
      }
    }
    print(dist_vec)
    temp_df <- data.frame(matrix(ncol = 3, nrow = length(dist_vec)))
    temp_df[,1] <- window_index
    temp_df[,2] <- sliding_window[window_index, "start_date"]
    temp_df[,3] <- dist_vec
    
    violin_rsvB_GER_df <- rbind(violin_rsvB_GER_df, temp_df)
    
    number_pairwise_dist <- (acc_length*(acc_length-1))/2
    mean <- dist_sum/number_pairwise_dist
    dist_mean_GER[window_index, "dist_mean"] <- mean
    
    dist_mean_GER[window_index, "dist_std"] <- sd(unlist(dist_vec))
    dist_mean_GER[window_index, "num"] <- acc_length
  }
}

colnames(violin_rsvB_GER_df) <- c("window", "start_date", "dist")

'ggplot(dist_mean_GER, aes(x = start_date, y = dist_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std))

ggplot(violin_rsvB_GER_df, aes(x = as.factor(start_date), y = dist, colour = "red")) +
  geom_violin()'

## EU
dist_mean_EU <- data.frame("window" =  1:nrow(sliding_window), "num" = 0, "dist_mean" = 0, "dist_std" = 0, "start_date" = sliding_window$start_date)
violin_rsvB_EU_df <- data.frame(matrix(ncol = 3, nrow = 0))

for(window_index in 1:nrow(sliding_window)) {
  dist_sum <- 0
  dist_vec <- c()
  acc_list <- c()
  for(acc_index0 in 1:nrow(meta_rsvB_EU)) {
    if(meta_rsvB_EU$Accession[acc_index0] != RefSeq_rsvB & 
       decimal_date(meta_rsvB_EU$Collection_Date[acc_index0]) >= sliding_window[window_index, "start_date"] & 
       decimal_date(meta_rsvB_EU$Collection_Date[acc_index0]) <= sliding_window[window_index, "end_date"]) {
      acc_list <- append(acc_list, meta_rsvB_EU$Accession[acc_index0])
    }
  }
  print(window_index)
  print(acc_list)
  
  acc_length <- length(acc_list)
  if(acc_length <= 1) {
    dist_mean_EU[window_index, "dist_mean"] <- NA
  } else {
    for(acc_index1 in 1:acc_length) {
      for(acc_index2 in acc_index1:acc_length) {
        dist_sum <- dist_rsvB_EU[acc_list[acc_index1], 
                                    acc_list[acc_index2]] + dist_sum
        if(acc_index1 != acc_index2)  { #dist_rsvB_EU[acc_list[acc_index1], acc_list[acc_index2]] != 0
          dist_vec <- append(dist_vec, dist_rsvB_EU[acc_list[acc_index1], acc_list[acc_index2]])
        }
      }
    }
    number_pairwise_dist <- (acc_length*(acc_length-1))/2
    mean <- dist_sum/number_pairwise_dist
    dist_mean_EU[window_index, "dist_mean"] <- mean
    
    print(dist_vec)
    
    temp_df <- data.frame(matrix(ncol = 3, nrow = length(dist_vec)))
    temp_df[,1] <- window_index
    temp_df[,2] <- sliding_window[window_index, "start_date"]
    temp_df[,3] <- dist_vec
    
    violin_rsvB_EU_df <- rbind(violin_rsvB_EU_df, temp_df)
    
    dist_mean_EU[window_index, "dist_std"] <- sd(unlist(dist_vec))
    #sqrt(sum(unlist(lapply(dist_vec, function(x) (x - mean)**2)))/(number_pairwise_dist))
    dist_mean_EU[window_index, "num"] <- acc_length
  }
}

colnames(violin_rsvB_EU_df) <- c("window", "start_date", "dist")
  
'ggplot(dist_mean_EU, aes(x = start_date, y = dist_mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std))

ggplot(violin_rsvB_EU_df, aes(x = as.factor(window), y = dist)) +
  geom_violin()'

# Plot EU + GER

sliding_window_plot <- ggplot(dist_mean_GER, aes(x = start_date, y = dist_mean)) +
  geom_errorbar(data = dist_mean_EU, aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "grey") +
  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "black", alpha = 0.45) +
  geom_point(colour = "red") +
  geom_point(data = dist_mean_EU, aes(x = start_date, y = dist_mean), colour = "blue")
sliding_window_plot

# Violin Plot

violin_rsvB_GER_df$'EU/GER' <- "GER"
violin_rsvB_EU_df$'EU/GER' <- "EU"

violin_GER_EU_df <- rbind(violin_rsvB_GER_df, violin_rsvB_EU_df)
violin_GER_EU_df <- arrange(violin_GER_EU_df, window)

violin_rsvB_df <- data.frame("window" = sliding_window$index)
violin_rsvB_df <- full_join(violin_rsvB_df, violin_GER_EU_df, by = "window")

ggplot(violin_rsvB_df, aes(x = factor(window), y = dist, fill = `EU/GER`, colour = `EU/GER`)) +
  geom_violin() +
  scale_colour_manual(values = c("GER" = "red", "EU" = "blue")) +
  scale_fill_manual(values = c("GER" = "red", "EU" = "blue")) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))

#########################################
'# Assign window to collection date
sliding_window_meta_rsvB <- subset(meta_rsvB, Accession != RefSeq_rsvB)
sliding_window_meta_rsvB <- sliding_window_meta_rsvB[, c("Accession", "Country", "Collection_Date")]
sliding_window_meta_rsvB$Collection_Date_dec <- decimal_date(sliding_window_meta_rsvB$Collection_Date)

sliding_windows_index <- list()
for(j in 1:(nrow(sliding_window_meta_rsvB))) {
  sliding_windows_index[j] <- list(which(sliding_window_meta_rsvB$Collection_Date_dec[j] >= sliding_window$start_date & 
                                           sliding_window_meta_rsvB$Collection_Date_dec[j] <= sliding_window$end_date))
}

for(n in 1:length(sliding_windows_index)) {
  for(m in 1:length(sliding_windows_index[[n]])) {
    sliding_windows_index[[n]][m] <- ifelse(as.integer(sliding_windows_index[[n]][m]) < 100, ifelse(as.integer(sliding_windows_index[[n]][m]) < 10, paste0(00, sliding_windows_index[[n]][m]), paste0(0, sliding_windows_index[[n]][m])), as.character(sliding_windows_index[[n]][m]))
  }
}

for(l in 1:(nrow(sliding_window_meta_rsvB))) {
  window_temp <- ""
  for(k in 1:length(sliding_windows_index[[l]])) {
    window_temp <- paste(window_temp, sliding_windows_index[[l]][k], sep = "W")
  }
  sliding_window_meta_rsvB$Window[l] <- window_temp
}

RefA_df <- data.frame("Accession" = RefSeq_rsvB, "Country" = NA, "Collection_Date" = NA, 
                      "Collection_Date_dec" = NA, "Window" = NA,
                      stringsAsFactors = FALSE)

sliding_window_meta_rsvB_GER <- sliding_window_meta_rsvB %>% subset(Country == "Germany") %>% rbind(RefA_df) %>% arrange(Accession)
sliding_window_meta_rsvB_EU <- sliding_window_meta_rsvB %>% subset(Country != "Germany") %>% rbind(RefA_df) %>% arrange(Accession)'