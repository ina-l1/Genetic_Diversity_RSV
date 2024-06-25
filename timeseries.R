# Data visualization: Time series

library(dplyr)
library(tidyr)
library(ggplot2)

# Read metadata
rsvAB <- read.csv("~/Yale_Projects/Genetic_Diversity_RSV/Europe/rsvAB_metadata_EU.csv", row.names = 1)
rsvAB$Collection_Date <- as.Date(rsvAB$Collection_Date, format = "%Y-%m-%d") #%m/%d/%Y

# Number of sequences per season
tab_rsvAB_season <- as.data.frame.matrix(table(rsvAB$Collection_Season, rsvAB$Type))

# Time series
timeseries_plot <- ggplot(rsvAB, aes(x = Collection_Date, color = Type, fill = Type))+
  geom_histogram(binwidth = 7)+
  scale_x_date(
    date_breaks = "3 months", 
    date_minor_breaks = "1 week", 
    date_labels = "%W"
  )+
  facet_grid(Type ~ .)
timeseries_plot

timeseries_stack_plot <- ggplot(rsvAB, aes(x = Collection_Date, color = Type, fill = Type))+
  geom_histogram(binwidth = 7, position = "stack", alpha = 0.5)+
  scale_x_date(
    date_breaks = "3 months",
    date_minor_breaks = "1 week",
    date_labels = "%W"
  )
timeseries_stack_plot

timeseries_bar_plot <- ggplot(rsvAB, aes(x = Collection_Date, color = Type, fill = Type))+
  geom_bar(position = "stack")+
  scale_x_date(
    date_breaks = "3 months",
    date_minor_breaks = "1 week",
    date_labels = "%W"
  )
timeseries_bar_plot