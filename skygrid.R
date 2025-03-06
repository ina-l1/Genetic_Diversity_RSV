# SkyGrid reconstruction
# BEAST v1.10.4

# DEU vs. EU 

library(tidyr)
library(ggplot2)
library(MMWRweek)
library(lubridate)
library(patchwork)
library(scales)

#################################################################################

# Base directory

base_dir <- "~/Yale_Projects/Genetic_Diversity_RSV/"

#################################################################################

# Read data

path_skygrid_rsvA_DEU_data <- file.path(base_dir, "BEAST_output", "skygrid", "skygrid_rsvA_GER_data.csv")
skygrid_rsvA_DEU_data <- read.csv(path_skygrid_rsvA_DEU_data, header = TRUE)

path_skygrid_rsvB_DEU_data <- file.path(base_dir, "BEAST_output", "skygrid", "skygrid_rsvB_GER_data.csv")
skygrid_rsvB_DEU_data <- read.csv(path_skygrid_rsvB_DEU_data, header = TRUE)

path_skygrid_rsvA_EU_data <- file.path(base_dir, "BEAST_output", "skygrid", "skygrid_rsvA_noGER_data.csv")
skygrid_rsvA_EU_data <- read.csv(path_skygrid_rsvA_EU_data, header = TRUE)

path_skygrid_rsvB_EU_data <- file.path(base_dir, "BEAST_output", "skygrid", "skygrid_rsvB_noGER_data.csv")
skygrid_rsvB_EU_data <- read.csv(path_skygrid_rsvB_EU_data, header = TRUE)

# Reference lines for each season

season_start_df <- data.frame("year" = 2000:2023, "week" = 40, "season_start" = NA)
season_start_df$season_start <- MMWRweek2Date(season_start_df$year, season_start_df$week)
season_start_df$season_start_dec <- decimal_date(season_start_df$season_start)

## Combine data
skygrid_rsvA_DEU_data$type <- "RSV-A"
skygrid_rsvB_DEU_data$type <- "RSV-B"

skygrid_DEU_data <- rbind(skygrid_rsvA_DEU_data, skygrid_rsvB_DEU_data)

skygrid_rsvA_EU_data$type <- "RSV-A"
skygrid_rsvB_EU_data$type <- "RSV-B"

skygrid_EU_data <- rbind(skygrid_rsvA_EU_data, skygrid_rsvB_EU_data)

# Plot SkyGrid

## RSV-A
ggplot(skygrid_rsvA_DEU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  xlim(2000, 2022.93)

ggplot(skygrid_rsvA_EU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  xlim(2011, 2022.86)

## RSV-B
ggplot(skygrid_rsvB_DEU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  xlim(2000, 2022.93)

ggplot(skygrid_rsvB_EU_data, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70") +
  geom_line(aes(y = mean)) +
  scale_y_continuous(trans='log2') +
  xlim(2011, 2022.86)

## DEU + EU 

### RSV-A
skygrid_rsvA_DEU_EU <- ggplot(skygrid_rsvA_DEU_data, aes(x = time)) +
  #theme_minimal() +
  geom_ribbon(data = skygrid_rsvA_EU_data, aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +
  geom_line(data = skygrid_rsvA_EU_data, aes(y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.5) +
  geom_line(aes(y = mean)) +
  geom_vline(xintercept = season_start_df$season_start_dec[1:24], color = "gray", linetype = "dashed") +
  scale_y_continuous(trans='log2', labels = label_scientific(digits = 2)) +
  scale_x_continuous(breaks = seq(2010,2023, by = 1), limits = c(2011, 2022.93)) +
  ylab("Effective population size, log(mean)") +
  xlab("Collection date") +
  theme(axis.line.x = element_line(color = "black"),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.title.x = element_blank(),
        text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_blank())
skygrid_rsvA_DEU_EU

### RSV-B
skygrid_rsvB_DEU_EU <- ggplot(skygrid_rsvB_DEU_data, aes(x = time)) +
  #theme_minimal() +
  geom_ribbon(data = skygrid_rsvB_EU_data, aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +
  geom_line(data = skygrid_rsvB_EU_data, aes(y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.5) +
  geom_line(aes(y = mean)) +
  geom_vline(xintercept = season_start_df$season_start_dec[1:24], color = "gray", linetype = "dashed") +
  scale_y_continuous(trans='log2', labels = label_scientific(digits = 2)) +
  scale_x_continuous(breaks = seq(2000,2023, by = 1), limits = c(2006, 2022.93)) +
  ylab("Effective population size, log(mean)") +
  xlab("Collection date") +
  theme(axis.line.x = element_line(color = "black"),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.title.x = element_blank(),
        text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_blank())
skygrid_rsvB_DEU_EU
  

## RSV-A + RSV-B

ggplot(skygrid_rsvA_DEU_data, aes(x = time)) +
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) +
  geom_ribbon(data = skygrid_rsvB_DEU_data, aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +
  geom_line(data = skygrid_rsvB_DEU_data, aes(y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.5) +
  geom_line(aes(y = mean)) +
  geom_vline(xintercept = season_start_df$season_start_dec[1:24], color = "gray", linetype = "dashed") +
  geom_vline(xintercept = decimal_date(as_date("2020-9-1"))) +
  geom_vline(xintercept = decimal_date(as_date("2020-3-16"))) +
  scale_y_continuous(trans='log2') +
  scale_x_continuous(breaks = seq(2000,2023, by = 1), limits = c(2006, 2022.93)) +
  ylab("Effective population size, log(mean)") +
  xlab("Collection date") +
  theme(axis.line.x = element_line(color = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_line(),
        #axis.title.x = element_blank(),
        text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.background = element_blank())


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


combined_skygrid <- skygrid_rsvA_DEU_EU / skygrid_rsvB_DEU_EU + 
  plot_layout(widths = c(1, 1.3)) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 30))
combined_skygrid_path <- file.path(base_dir, "BEAST_output", "phylogenetic_tree", "skygrid_combined.png")
#ggsave(filename = combined_skygrid_path, width = 35, height = 35, units = "cm", limitsize = FALSE)
