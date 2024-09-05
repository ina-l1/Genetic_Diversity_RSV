# Sliding Window

# Validating temporal genetic diversity pattern of German sequences
# Countries: Netherlands, Spain, UK
# All included EU countries: Austria, Finland, France, Germany, Italy, Netherlands, Russia, Spain, United Kingdom 

library(tidyr)
library(lubridate)
library(MMWRweek)
library(ggplot2)
library(patchwork) #combines plots in one figure

#################################################################################

# Sequence data: Read metadata

meta_rsvA <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvA_ref_metadata_EU.csv") #already in alphabetical order
meta_rsvB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvB_ref_metadata_EU.csv")

meta_rsvA$Collection_Date <- as_date(meta_rsvA$Collection_Date) 
meta_rsvB$Collection_Date <- as_date(meta_rsvB$Collection_Date)

RefSeq_rsvA <- "NC_038235.1"
RefSeq_rsvB <- "NC_001781.1"

## Germany (DEU)
meta_rsvA_DEU <- subset(meta_rsvA, Country == "Germany" | Accession == RefSeq_rsvA) #same as NCBI_rsvA_wgs_germany_2015.csv
meta_rsvB_DEU <- subset(meta_rsvB, Country == "Germany" | Accession == RefSeq_rsvB) #same as NCBI_rsvB_wgs_germany_2015.csv

## Netherlands (NLD)
meta_rsvA_NLD <- subset(meta_rsvA, Country == "Netherlands" | Accession == RefSeq_rsvA)
meta_rsvB_NLD <- subset(meta_rsvB, Country == "Netherlands" | Accession == RefSeq_rsvB)

## Spain (ESP)
meta_rsvA_ESP <- subset(meta_rsvA, Country == "Spain" | Accession == RefSeq_rsvA)
meta_rsvB_ESP <- subset(meta_rsvB, Country == "Spain" | Accession == RefSeq_rsvB)

## UK (GBR)
meta_rsvA_GBR <- subset(meta_rsvA, Country == "United Kingdom" | Accession == RefSeq_rsvA)
meta_rsvB_GBR <- subset(meta_rsvB, Country == "United Kingdom" | Accession == RefSeq_rsvB)

## Non-German EU sequences
meta_rsvA_EU <- subset(meta_rsvA, Country != "Germany" | Accession == RefSeq_rsvA)
meta_rsvB_EU <- subset(meta_rsvB, Country != "Germany" | Accession == RefSeq_rsvB)

# Case count: Read metadata

case_flunet <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/VIW_FNT.csv")
case_flunet$ISO_WEEKSTARTDATE <- as_date(case_flunet$ISO_WEEKSTARTDATE)
case_flunet$MMWR_WEEKSTARTDATE <- as_date(case_flunet$MMWR_WEEKSTARTDATE)

start_date <- MMWRweek2Date(2014, 40)
end_date <- MMWRweek2Date(2023, 39) #Most recent sequence data available from 2022/2023 season

case_DEU <- case_flunet %>% subset(COUNTRY_CODE == "DEU") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_NLD <- case_flunet %>% subset(COUNTRY_CODE == "NLD") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_ESP <- case_flunet %>% subset(COUNTRY_CODE == "ESP") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_GBR <- case_flunet %>% subset(COUNTRY_CODE == "X09") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date)

case_EU <- case_flunet %>% subset(COUNTRY_CODE == "AUT" | COUNTRY_CODE == "FIN" | COUNTRY_CODE == "FRA" | COUNTRY_CODE == "ITA" |
                                    COUNTRY_CODE == "NLD" | COUNTRY_CODE == "RUS" | COUNTRY_CODE == "ESP" | COUNTRY_CODE == "X09") %>% 
  subset(select = c(COUNTRY_CODE, COUNTRY_AREA_TERRITORY, ISO_WEEKSTARTDATE, ISO_YEAR, ISO_WEEK, MMWR_WEEKSTARTDATE, MMWR_YEAR, MMWR_WEEK, RSV, MMWRYW, ISOYW)) %>%
  subset(MMWR_WEEKSTARTDATE >= start_date & MMWR_WEEKSTARTDATE <= end_date) #GenBank sequences only available for these countries

# Distance matrices

dist_rsvA_DEU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/evodist_rsvA.csv") ##
dist_rsvA_NLD <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_netherlands.csv") ##
dist_rsvA_ESP <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_spain.csv") ##
dist_rsvA_GBR <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_UK.csv") ##
dist_rsvA_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvA_EU_noGer.csv") ##

dist_rsvB_DEU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/evodist_rsvB.csv") ##
dist_rsvB_NLD <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_netherlands.csv") ##
dist_rsvB_ESP <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_spain.csv") ##
dist_rsvB_GBR <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_UK.csv") ##
dist_rsvB_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/evodist_rsvB_EU_noGer.csv") ##

dist_matrix_func <- function(distance_matrix) {
  distance_matrix <- arrange(distance_matrix, distance_matrix[,1])
  distance_matrix <- distance_matrix[,-1]
  distance_matrix <- distance_matrix %>% select(order(colnames(.)))
}

dist_rsvA_DEU <- dist_matrix_func(distance_matrix = dist_rsvA_DEU)
dist_rsvB_DEU <- dist_matrix_func(distance_matrix = dist_rsvB_DEU)

dist_rsvA_NLD <- dist_matrix_func(distance_matrix = dist_rsvA_NLD)
dist_rsvB_NLD <- dist_matrix_func(distance_matrix = dist_rsvB_NLD)

dist_rsvA_ESP <- dist_matrix_func(distance_matrix = dist_rsvA_ESP)
dist_rsvB_ESP <- dist_matrix_func(distance_matrix = dist_rsvB_ESP)

dist_rsvA_GBR <- dist_matrix_func(distance_matrix = dist_rsvA_GBR)
dist_rsvB_GBR <- dist_matrix_func(distance_matrix = dist_rsvB_GBR)

dist_rsvA_EU <- dist_matrix_func(distance_matrix = dist_rsvA_EU)
dist_rsvB_EU <- dist_matrix_func(distance_matrix = dist_rsvB_EU)

rownames(dist_rsvA_DEU) <- meta_rsvA_DEU$Accession
colnames(dist_rsvA_DEU) <- meta_rsvA_DEU$Accession

rownames(dist_rsvA_NLD) <- meta_rsvA_NLD$Accession
colnames(dist_rsvA_NLD) <- meta_rsvA_NLD$Accession

rownames(dist_rsvA_ESP) <- meta_rsvA_ESP$Accession
colnames(dist_rsvA_ESP) <- meta_rsvA_ESP$Accession

rownames(dist_rsvA_GBR) <- meta_rsvA_GBR$Accession
colnames(dist_rsvA_GBR) <- meta_rsvA_GBR$Accession

rownames(dist_rsvA_EU) <- meta_rsvA_EU$Accession
colnames(dist_rsvA_EU) <- meta_rsvA_EU$Accession

rownames(dist_rsvB_DEU) <- meta_rsvB_DEU$Accession
colnames(dist_rsvB_DEU) <- meta_rsvB_DEU$Accession

rownames(dist_rsvB_NLD) <- meta_rsvB_NLD$Accession
colnames(dist_rsvB_NLD) <- meta_rsvB_NLD$Accession

rownames(dist_rsvB_ESP) <- meta_rsvB_ESP$Accession
colnames(dist_rsvB_ESP) <- meta_rsvB_ESP$Accession

rownames(dist_rsvB_GBR) <- meta_rsvB_GBR$Accession
colnames(dist_rsvB_GBR) <- meta_rsvB_GBR$Accession

rownames(dist_rsvB_EU) <- meta_rsvB_EU$Accession
colnames(dist_rsvB_EU) <- meta_rsvB_EU$Accession

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

# Sliding window function

sliding_window_func <- function(country) {
  
  if(country == "DEU" & rsvAB_choose == "rsvA"){
    meta <- meta_rsvA_DEU
    dist_matrix <- dist_rsvAB_DEU
  } else if(country == "ESP" & rsvAB_choose == "rsvA") {
    meta <- meta_rsvA_ESP
    dist_matrix <- dist_rsvAB_ESP
  } else if(country == "NLD" & rsvAB_choose == "rsvA") {
    meta <- meta_rsvA_NLD
    dist_matrix <- dist_rsvAB_NLD
  } else if(country == "GBR" & rsvAB_choose == "rsvA") {
    meta <- meta_rsvA_GBR
    dist_matrix <- dist_rsvAB_GBR
  } else if(country == "EU" & rsvAB_choose == "rsvA") {
    meta <- meta_rsvA_EU
    dist_matrix <- dist_rsvAB_EU
  } else if(country == "DEU" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_DEU
    dist_matrix <- dist_rsvAB_DEU
  } else if(country == "ESP" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_ESP
    dist_matrix <- dist_rsvAB_ESP
  } else if(country == "NLD" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_NLD
    dist_matrix <- dist_rsvAB_NLD
  } else if(country == "GBR" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_GBR
    dist_matrix <- dist_rsvAB_GBR
  } else if(country == "EU" & rsvAB_choose == "rsvB") {
    meta <- meta_rsvB_EU
    dist_matrix <- dist_rsvAB_EU
  }
  
  dist_mean_df <- data.frame("window" =  1:nrow(sliding_window), "num" = 0, "dist_mean" = 0, "dist_std" = 0, "start_date" = sliding_window$start_date)
  
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
      dist_mean_df[window_index, "dist_mean"] <- NA
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
      dist_mean_df[window_index, "dist_mean"] <- mean
      
      dist_mean_df[window_index, "dist_std"] <- sd(unlist(dist_vec))
      dist_mean_df[window_index, "num"] <- acc_length
    }
  }
  
  return(dist_mean_df)
}

#################################################################################

# Case count function

case_count_DEU <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_DEU, sum)
case_count_DEU$start_date <- decimal_date(case_count_DEU$MMWR_WEEKSTARTDATE)

case_count_ESP <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_ESP, sum)
case_count_ESP$start_date <- decimal_date(case_count_ESP$MMWR_WEEKSTARTDATE)

case_count_NLD <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_NLD, sum)
case_count_NLD$start_date <- decimal_date(case_count_NLD$MMWR_WEEKSTARTDATE)

case_count_GBR <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_GBR, sum)
case_count_GBR$start_date <- decimal_date(case_count_GBR$MMWR_WEEKSTARTDATE)

case_count_EU <- aggregate(RSV ~ MMWR_WEEKSTARTDATE, data = case_EU, sum)
case_count_EU$start_date <- decimal_date(case_count_EU$MMWR_WEEKSTARTDATE)

case_count_func <- function(country) {
  
  if(country == "DEU") {
    case_df <- case_count_DEU
  } else if(country == "ESP") {
    case_df <- case_count_ESP
  } else if(country == "NLD") {
    case_df <- case_count_NLD
  } else if(country == "GBR") {
    case_df <- case_count_GBR
  } else if(country == "EU") {
    case_df <- case_count_EU
  }
  
  case_count_df <- data.frame("window" =  1:nrow(sliding_window), "start_date" = sliding_window$start_date, "sum" = 0)
  
  for(window_index in 1:nrow(sliding_window)) {
    case_count_sum <- 0
    print(window_index)
    for(case_index in 1:nrow(case_df)) {
      if(case_df$start_date[case_index] >= sliding_window[window_index, "start_date"] & 
         case_df$start_date[case_index] <= sliding_window[window_index, "end_date"]) {
        case_count_sum <- case_count_sum + case_df$RSV[case_index]
      }
    }
    case_count_df[window_index, "sum"] <- case_count_sum
  }
  return(case_count_df)
}

#################################################################################

# Choose RSV-A or RSV-B ("rsvA" or "rsvB")

rsvAB_choose <- "rsvA" ##

if(rsvAB_choose == "rsvA") {
  dist_rsvAB_DEU <- dist_rsvA_DEU
  dist_rsvAB_ESP <- dist_rsvA_ESP
  dist_rsvAB_NLD <- dist_rsvA_NLD
  dist_rsvAB_GBR <- dist_rsvA_GBR
  dist_rsvAB_EU <- dist_rsvA_EU
  RefSeq <- RefSeq_rsvA
  
} else if(rsvAB_choose == "rsvB") {
  dist_rsvAB_DEU <- dist_rsvB_DEU
  dist_rsvAB_ESP <- dist_rsvB_ESP
  dist_rsvAB_NLD <- dist_rsvB_NLD
  dist_rsvAB_GBR <- dist_rsvB_GBR
  dist_rsvAB_EU <- dist_rsvB_EU
  RefSeq <- RefSeq_rsvB
} else {
  print("Choose rsvA or rsvB!")
}

#################################################################################

# Distances

#dist_mean_DEU <- sliding_window_func(country = "DEU")
dist_mean_ESP <- sliding_window_func(country = "ESP")
dist_mean_NLD <- sliding_window_func(country = "NLD")
#dist_mean_GBR <- sliding_window_func(country = "GBR")
#dist_mean_EU <- sliding_window_func(country = "EU")

#################################################################################

# Case count (FluNet)

#case_win_DEU <- case_count_func("DEU")
case_win_ESP <- case_count_func("ESP")
case_win_NLD <- case_count_func("NLD")
#case_win_GBR <- case_count_func("GBR")
#case_win_EU <- case_count_func("EU")

#################################################################################

# Plots

dist_mean_1 <- dist_mean_DEU ##
dist_mean_2 <- dist_mean_NLD ##

case_win_1 <- case_win_DEU ##
case_win_2 <- case_win_NLD ##

## Evo (changed 20000 to 30000!)
sliding_window_plot <- ggplot(dist_mean_1, aes(x = start_date, y = dist_mean)) +
  geom_errorbar(data = dist_mean_2, aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "lightgrey") +
#  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "black", alpha = 0.45) +
  geom_point(data = dist_mean_2, aes(x = start_date, y = dist_mean), colour = "blue") + #darkgreen
#  geom_point(colour = "red") +
#  geom_line(aes(y = num/30000), colour = "red", linewidth = 0.8) +
#  geom_line(aes(y = dist_mean_2$num/30000), colour = "darkgreen", linewidth = 0.8) +
#  geom_line(aes(y = (num + dist_mean_2$num)/30000), colour = "black", alpha = 0.4) +
#  geom_line(data = case_win_1, aes(x = start_date, y = sum/200000), colour = "red", alpha = 0.45) +
#  geom_line(data = case_win_2, aes(x = start_date, y = sum/200000), colour = "darkgreen", alpha = 0.45) +
  scale_x_continuous(breaks = c(2014:2023)) +
#  scale_y_continuous("distance (mean)", sec.axis = sec_axis(~.*30000, name="sequence count"), limits = c(0, 0.02)) + # IF NECESSARY ADJUST Y-AXIS
  theme_minimal() +
  xlab("time") +
  ylab("evolutionary distance (mean)")
sliding_window_plot

## SNP
sliding_window_plot <- ggplot(dist_mean_1, aes(x = start_date, y = dist_mean)) +
  geom_errorbar(data = dist_mean_2, aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "lightgrey") +
  geom_errorbar(aes(ymin = dist_mean-dist_std, ymax = dist_mean+dist_std), colour = "black", alpha = 0.45) +
  geom_point(colour = "red") +
  geom_point(data = dist_mean_2, aes(x = start_date, y = dist_mean), colour = "darkgreen") +
#  geom_line(aes(y = num/1.33), colour = "red") +
#  geom_line(aes(y = dist_mean_2$num/1.33), colour = "darkgreen") +
#  geom_line(aes(y = (num + dist_mean_2$num)/1.33), colour = "black", alpha = 0.4) +
  scale_x_continuous(breaks = c(2014:2023)) +
#  scale_y_continuous("distance (mean)", sec.axis = sec_axis(~.*1.33, name="sequence count"), limits = c(0, 300)) + # IF NECESSARY ADJUST Y-AXIS
  theme_minimal() +
  xlab("time") +
  ylab("SNP distance (mean)")
#sliding_window_plot

## Sequence count
sequence_count_plot <- ggplot(data = dist_mean_1) +
  geom_line(aes(x = start_date, y = num), colour = "red") +
  geom_line(aes(x = start_date, y = dist_mean_2$num), colour = "blue") +
  scale_x_continuous(breaks = c(2014:2023)) + 
  theme_minimal() +
  xlab("time") +
  ylab("sequence count")
sequence_count_plot

## Case count
case_count_plot <- ggplot() + 
  geom_line(data = case_win_1, aes(x = start_date, y = sum*40), colour = "red") +
  geom_line(data = case_win_2, aes(x = start_date, y = sum), colour = "blue") +
  scale_x_continuous(breaks = c(2014:2023)) +
  theme_minimal() +
  scale_y_continuous("case count (EU)", sec.axis = sec_axis(~./40, name="case count (DEU)"), limits = c(0, 45000)) +
  xlab("time") +
  ylab("case count")
case_count_plot

# Combined plot

combined_plot <- sliding_window_plot / sequence_count_plot / case_count_plot + plot_layout(heights = c(2, 1, 1))
combined_plot

ggsave("~/Yale_Projects/Genetic_Diversity_RSV/Plots/sliding_window_rsvB_evo_EU_combined.png", width = 50, height = 25, units = "cm", limitsize = FALSE)
