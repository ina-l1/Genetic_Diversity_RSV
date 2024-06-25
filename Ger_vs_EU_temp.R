# German sequences vs. sequences from other EU countries
# Temporary: Code should be merged with other files 

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
####################
snpdist_rsvA_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvA.csv")
snpdist_rsvA_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_EU_noGer.csv")

snpdist_rsvA_GER <- arrange(snpdist_rsvA_GER, snpdist_rsvA_GER[,1])
snpdist_rsvA_GER <- snpdist_rsvA_GER[,-1]
#snpdist_rsvA_GER$Accession <- meta_rsvA_GER$Accession #shorter name
rownames(snpdist_rsvA_GER) <- meta_rsvA_GER$Accession
snpdist_rsvA_GER <- snpdist_rsvA_GER %>% select(order(colnames(.)))
colnames(snpdist_rsvA_GER) <- meta_rsvA_GER$Accession
######################

snpdist_rsvA_EU <- arrange(snpdist_rsvA_EU, snpdist_rsvA_EU[,1])
snpdist_rsvA_EU <- snpdist_rsvA_EU[,-1]
#snpdist_rsvA_EU$Accession <- meta_rsvA_EU$Accession #shorter name
rownames(snpdist_rsvA_EU) <- meta_rsvA_EU$Accession
snpdist_rsvA_EU <- snpdist_rsvA_EU %>% select(order(colnames(.)))
colnames(snpdist_rsvA_EU) <- meta_rsvA_EU$Accession


### RSV B
snpdist_rsvB_GER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/snpdist_rsvB.csv")
snpdist_rsvB_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvB_EU_noGer.csv")

snpdist_rsvB_GER <- arrange(snpdist_rsvB_GER, snpdist_rsvB_GER[,1])
snpdist_rsvB_GER <- snpdist_rsvB_GER[,-1]
snpdist_rsvB_GER$Accession <- meta_rsvB_GER$Accession #shorter name
rownames(snpdist_rsvB_GER) <- snpdist_rsvB_GER$Accession
snpdist_rsvB_GER$Accession <- NULL
snpdist_rsvB_GER <- snpdist_rsvB_GER %>% select(order(colnames(.)))

snpdist_rsvB_EU <- arrange(snpdist_rsvB_EU, snpdist_rsvB_EU[,1])
snpdist_rsvB_EU <- snpdist_rsvB_EU[,-1]
snpdist_rsvB_EU$Accession <- meta_rsvB_EU$Accession #shorter name
rownames(snpdist_rsvB_EU) <- snpdist_rsvB_EU$Accession
snpdist_rsvB_EU$Accession <- NULL
snpdist_rsvB_EU <- snpdist_rsvB_EU %>% select(order(colnames(.)))

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

################################################################
# Sliding Window
# Timeframe: 3 weeks 

# Define windows 
# Start: 2014 W40
# End: 2023 W39
# 2015 and 2020 have 53 instead of 53 weeks

snpdist_rsvA_EUGER <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/snpdist_rsvA_EU.csv")

snpdist_rsvA_EUGER <- arrange(snpdist_rsvA_EUGER, snpdist_rsvA_EUGER[,1])
snpdist_rsvA_EUGER <- snpdist_rsvA_EUGER[,-1]
#snpdist_rsvA_EUGER$Accession <- meta_rsvA_GER$Accession #shorter name
rownames(snpdist_rsvA_EUGER) <- meta_rsvA$Accession
snpdist_rsvA_EUGER <- snpdist_rsvA_EUGER %>% select(order(colnames(.)))
colnames(snpdist_rsvA_EUGER) <- meta_rsvA$Accession

# Creating df with weeks for each year

week <- rep(c(1:53), times = 10)
year <- rep(c(2014:2023), each = 53)
dates_df <- data.frame(year, week)

vec <- which(dates_df$week == 53 & (dates_df$year != 2015 & dates_df$year!=2020)) #years with no 53rd week
vec <- append(vec, which((between(dates_df$week, 1, 39) & dates_df$year == 2014)|(between(dates_df$week, 40, 53) & dates_df$year == 2023)))

dates_df <- dates_df[-vec,]

dates_df$date <- decimal_date(MMWRweek2Date(dates_df$year, dates_df$week))

dates_df$index <- 1:nrow(dates_df)
rownames(dates_df) <- dates_df$index

# Defining window with start and end date
sliding_window <- data.frame(matrix(ncol = 3))
colnames(sliding_window) <- c("index","start_date","end_date")
sliding_window <- sliding_window[-1,]

for(i in 1:(nrow(dates_df)-2)) {
  sliding_window[i,"index"] <- i
  sliding_window[i,"start_date"] <- dates_df$date[i]
  sliding_window[i,"end_date"] <- dates_df$date[i+2]
}

# Assign window to collection date
sliding_window_meta_rsvA <- subset(meta_rsvA, Accession != RefSeq_rsvA)
sliding_window_meta_rsvA <- sliding_window_meta_rsvA[, c("Accession", "Country", "Collection_Date")]
sliding_window_meta_rsvA$Collection_Date_dec <- decimal_date(sliding_window_meta_rsvA$Collection_Date)

sliding_windows_index <- list()
for(j in 1:(nrow(sliding_window_meta_rsvA))) {
  sliding_windows_index[j] <- list(which(sliding_window_meta_rsvA$Collection_Date_dec[j] >= sliding_window$start_date & 
                       sliding_window_meta_rsvA$Collection_Date_dec[j] <= sliding_window$end_date))
}

for(n in 1:length(sliding_windows_index)) {
  for(m in 1:length(sliding_windows_index[[n]])) {
  sliding_windows_index[[n]][m] <- ifelse(as.integer(sliding_windows_index[[n]][m]) < 100, ifelse(as.integer(sliding_windows_index[[n]][m]) < 10, paste0(00, sliding_windows_index[[n]][m]), paste0(0, sliding_windows_index[[n]][m])), as.character(sliding_windows_index[[n]][m]))
  }
}

for(l in 1:(nrow(sliding_window_meta_rsvA))) {
  window_temp <- ""
  for(k in 1:length(sliding_windows_index[[l]])) {
    window_temp <- paste(window_temp, sliding_windows_index[[l]][k], sep = "W")
  }
  sliding_window_meta_rsvA$Window[l] <- window_temp
}

RefA_df <- data.frame("Accession" = RefSeq_rsvA, "Country" = NA, "Collection_Date" = NA, 
                      "Collection_Date_dec" = NA, "Window" = NA,
                      stringsAsFactors = FALSE)

sliding_window_meta_rsvA_GER <- sliding_window_meta_rsvA %>% subset(Country == "Germany") %>% rbind(RefA_df) %>% arrange(Accession)
sliding_window_meta_rsvA_EU <- sliding_window_meta_rsvA %>% subset(Country != "Germany") %>% rbind(RefA_df) %>% arrange(Accession)

# SNP distance
snp_window <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Sliding Window/sliding_window_rsvA_assign.csv")

#meta_rsvA$Accession[which(snp_window$W172, TRUE)]
#snpdist_rsvA_EU["LR699315.1", "LR699315.1"]
#snpdist_rsvA_EU["LR699315.1", "MZ515623.1"]

#test <- paste0("W", 172)
#snpdist_rsvA_EU[meta_rsvA$Accession[which(snp_window[, 174], TRUE)][1], meta_rsvA$Accession[which(snp_window[, 174], TRUE)][4]]

snp_mean <- data.frame("window" =  1:nrow(sliding_window), "snp_mean" = 0, "snp_std" = 0, "start_date" = sliding_window$start_date)

for(window_index in 1:(nrow(sliding_window)-2)) {
  snp_sum <- 0
  snp_vec <- c()
  acc_list <- c()
  #to improve performance maybe sort and CUT 
  for(acc_index0 in 1:nrow(meta_rsvA)) {
    if(meta_rsvA$Accession[acc_index0] != RefSeq_rsvA & decimal_date(meta_rsvA$Collection_Date[acc_index0]) >= sliding_window[window_index, "start_date"] & decimal_date(meta_rsvA$Collection_Date[acc_index0]) <= sliding_window[window_index, "end_date"]) {
      acc_list <- append(acc_list, meta_rsvA$Accession[acc_index0])
      }  #loop metadaten: check if collection date in window, if yes append to acc_list, acc_list <- meta_rsvA$Accession[which(snp_window[, window_index+2], TRUE)]
  }
  print(window_index)
  print(acc_list)
  
  acc_length <- length(acc_list)
  if(acc_length <= 1) {
    snp_mean[window_index, "snp_mean"] <- NA
  } else {
    for(acc_index1 in 1:acc_length) {
      for(acc_index2 in acc_index1:acc_length) {
        snp_sum <- snpdist_rsvA_EUGER[acc_list[acc_index1], 
                                      acc_list[acc_index2]] + snp_sum
        snp_vec <- append(snp_vec, snpdist_rsvA_EUGER[acc_list[acc_index1], acc_list[acc_index2]])
      }
    }
    number_pairwise_dist <- (acc_length*(acc_length-1))/2
    snp_mean[window_index, "snp_mean"] <- snp_sum/number_pairwise_dist
    snp_mean[window_index, "snp_std"] <- lapply(snp_vec, function(x) sqrt(sum(abs(x - snp_mean[window_index, "snp_mean"])**2))/number_pairwise_dist)
  }
}


ggplot(snp_mean, aes(x = start_date, y = snp_mean)) +
  geom_point()

#PLAYGROUND