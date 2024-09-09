# SkyGrid reconstruction
# BEAST v1.10.4

# DEU vs. EU 

library(tidyr)
library(ggplot2)

#################################################################################

# Base directory

base_dir <- "~/Yale_Projects/Genetic_Diversity_RSV/"

#################################################################################

# Read data

path_skygrid_rsvA_DEU_data <- file.path(base_dir, "BEAST_output", "skygrid", "skygrid_rsvA_GER_data.csv")
skygrid_rsvA_DEU_data <- read.csv(path_skygrid_rsvA_DEU_data, header = TRUE)

path_skygrid_rsvB_DEU_data <- file.path(base_dir, "BEAST_output", "skygrid", "skygrid_rsvB_GER_data.csv")
skygrid_rsvB_DEU_data <- read.csv(path_skygrid_rsvB_DEU_data, header = TRUE)

## Combine data
skygrid_rsvA_DEU_data$type <- "RSV-A"
skygrid_rsvB_DEU_data$type <- "RSV-B"

skygrid_DEU_data <- rbind(skygrid_rsvA_DEU_data, skygrid_rsvB_DEU_data)

# Plot SkyGrid

# DEU

## RSV-A
ggplot(skygrid_rsvA_DEU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  xlim(2000, 2022.93)

## RSV-B
ggplot(skygrid_rsvB_DEU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  xlim(2000, 2022.93)

## RSV-A + RSV-B
ggplot(skygrid_DEU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  ylab("log(mean)") +
  scale_x_continuous(breaks = seq(2000, 2022.93, 1)) +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  facet_grid(type ~ .)