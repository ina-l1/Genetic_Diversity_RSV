# Data visualization: Time series

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# Read metadata
rsvAB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Germany/rsvAB_metadata_GER.csv", row.names = 1)
rsvAB$Collection_Date <- as.Date(rsvAB$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
rsvAB$Date_round <- floor_date(rsvAB$Collection_Date, "week")

rsvAB_EU <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvAB_metadata_EU.csv", row.names = 1)
rsvAB_EU$Collection_Date <- as.Date(rsvAB_EU$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y
rsvAB_EU$Date_round <- floor_date(rsvAB_EU$Collection_Date, "week")

rsvAB_EU$'EU/GER' <- ifelse(rsvAB_EU$Country == "Germany", "GER", "EU")

#rsvAB_EU <- subset(rsvAB_EU, Country != "Germany")

# Number of sequences per season
tab_rsvAB_season <- as.data.frame.matrix(table(rsvAB$Collection_Season, rsvAB$Type))
tab_rsvAB_EU_season <- as.data.frame.matrix(table(rsvAB_EU$Collection_Season, rsvAB_EU$Type))

# Time series
timeseries_plot <- ggplot(rsvAB, aes(x = Collection_Date, color = Type, fill = Type))+
  geom_histogram(binwidth = 7)+
  scale_x_date(
    date_breaks = "3 months", 
    date_minor_breaks = "1 week", 
    date_labels = "%W"
  ) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  facet_grid(Type ~ .)
timeseries_plot

timeseries_stack_plot <- ggplot() +
  geom_bar(data = rsvAB, aes(x = Date_round, fill = Type)) +
  scale_fill_manual(values = c("red", "blue")) +
  xlab(label = "collection year") +
  ylab(label = "sequence count") +
  theme_minimal()
timeseries_stack_plot

#EU + GER
ggplot() +
  geom_bar(data = rsvAB_EU, aes(x = Date_round)) +
  geom_bar(data = rsvAB, aes(x = Date_round, fill = Type)) +
  scale_fill_manual(values = c("red", "blue")) +
  xlab(label = "collection year") +
  ylab(label = "sequence count") +
  theme_minimal()
