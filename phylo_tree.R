#Add important metadata (Season, Clade), Data Sorting and Tree Visualization 

library(dplyr)
library(MMWRweek)
library(ggplot2)
library(ape)
library(treeio)
library(ggtree)
library(lubridate)
library(patchwork)
library(here)

##################################################

# Base directory

base_dir <- here()

##################################################

# Metadata

path_RefSeq_rsvA <- file.path(base_dir, "RefSeq", "RefSeq_rsvA.csv")
RefSeq_rsvA <- read.csv(path_RefSeq_rsvA)

path_RefSeq_rsvB <- file.path(base_dir, "RefSeq", "RefSeq_rsvB.csv")
RefSeq_rsvB <- read.csv(path_RefSeq_rsvB)

path_rsvA_ref <- file.path(base_dir, "Europe", "rsvA_ref_metadata_EU.csv")
rsvA_ref <- read.csv(path_rsvA_ref)

path_rsvB_ref <- file.path(base_dir, "Europe", "rsvB_ref_metadata_EU.csv")
rsvB_ref <- read.csv(path_rsvB_ref)

rsvA_ref$Collection_Date <- as_date(rsvA_ref$Collection_Date) 
rsvB_ref$Collection_Date <- as_date(rsvB_ref$Collection_Date) 

rsvA_ref <- subset(rsvA_ref, Accession != "NC_038235.1") ###
rsvB_ref <- subset(rsvB_ref, Accession != "NC_001781.1") ###
rsvA_ref$'EU/GER' <- ifelse(rsvA_ref$Country == "Germany", "GER", "EU") ###
rsvB_ref$'EU/GER' <- ifelse(rsvB_ref$Country == "Germany", "GER", "EU") ###

# Reference lines for each season

season_start_df <- data.frame("year" = 2000:2023, "week" = 40, "season_start" = NA)
season_start_df$season_start <- MMWRweek2Date(season_start_df$year, season_start_df$week)
season_start_df$season_start_dec <- decimal_date(season_start_df$season_start)

##################################################

# Tree visualization: RSV A

path_beast_tree <- file.path(base_dir, "BEAST", "tree_files", "tree_rsvA_EU.trees")
beast_tree <- read.beast(path_beast_tree)

## Create tip labels

beast_tree_df <- as_tibble(beast_tree)
beast_tree_label <- sort(beast_tree_df$label)
beast_label_df <- data.frame(treelabel = beast_tree_label,
                             accession = rsvA_ref$Accession,
                             type = rsvA_ref$Type,
                             collection_date = rsvA_ref$Collection_Date,
                             season = rsvA_ref$Collection_Season,
                             country = rsvA_ref$Country,
                             EUGER = rsvA_ref$`EU/GER`)
beast_label_df$tip.label <- paste(beast_label_df$accession, beast_label_df$type, beast_label_df$season, beast_label_df$country, sep = "_")
beast_tree <- rename_taxa(beast_tree, data = beast_label_df, key = 1, value = tip.label)
beast_label_df <- relocate(beast_label_df, tip.label)

## Plot tree

beast_rsvA_plot <- ggtree(beast_tree, mrsd = "2022-12-07", aes(color = country), size = 1) %<+% beast_label_df + 
  geom_tippoint(aes(color = country), size=3) + #color = season or color = EUGER
  geom_tiplab(
    geom = "text",
    size = 5,
    color = "black"
  ) +
  geom_vline(xintercept = season_start_df$season_start_dec[12:24], color = "gray", linetype = "dashed") +
  #scale_color_manual(values = c("GER" = "red", "EU" =  "blue")) + #enable if EUGER
  guides(color = guide_legend(title = "Country")) + #title = "Season" or "EU/GER"
  theme_tree2() +
  theme(axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(2011,2023, by = 1)) +
  geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)>0.5, x=branch), vjust=-1, size = 4, color = "black")

path_treeA_plot_new <- file.path(base_dir, "BEAST_output", "phylogenetic_tree", "tree_rsvA_EU_country.pdf")
#ggsave(filename = path_treeA_plot_new, width = 100, height = 320, units = "cm", limitsize = FALSE)

##################################################

# Tree visualization: RSV B

path_beast_tree <- file.path(base_dir, "BEAST", "tree_files", "tree_rsvB_EU.trees")
beast_tree <- read.beast(path_beast_tree)

## Create tip labels
beast_tree_df <- as_tibble(beast_tree)
beast_tree_label <- sort(beast_tree_df$label)
beast_label_df <- data.frame(treelabel = beast_tree_label,
                             accession = rsvB_ref$Accession,
                             type = rsvB_ref$Type,
                             collection_date = rsvB_ref$Collection_Date,
                             season = rsvB_ref$Collection_Season,
                             country = rsvB_ref$Country,
                             EUGER = rsvB_ref$`EU/GER`)
beast_label_df$tip.label <- paste(beast_label_df$accession, beast_label_df$type, beast_label_df$season, beast_label_df$country, sep = "_")
beast_tree <- rename_taxa(beast_tree, data = beast_label_df, key = 1, value = tip.label)
beast_label_df <- relocate(beast_label_df, tip.label)

## Plot tree

beast_rsvB_plot <- ggtree(beast_tree, mrsd = "2022-12-07", aes(color = country), size = 2) %<+% beast_label_df + 
  geom_tippoint(aes(color = country), size=4) + #color = season or color = EUGER
  geom_tiplab(
    geom = "text",
    size = 4,
    color = "black"
  ) +
  geom_vline(xintercept = season_start_df$season_start_dec[7:24], color = "gray", linetype = "dashed") +
  #scale_color_manual(values = c("GER" = "red", "EU" =  "blue")) + #enable if EUGER
  guides(color = guide_legend(title = "Country")) + #title = "Season" or "EU/GER"
  theme_tree2() +
  theme(axis.text.x = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(2006,2023, by = 1)) +
  geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)>0.5, x=branch), vjust=-1, size = 4, color = "black")

path_treeB_plot_new <- file.path(base_dir, "BEAST_output", "phylogenetic_tree", "tree_rsvB_EU_country.pdf")
#ggsave(filename = path_treeB_plot_new, width = 200, height = 300, units = "cm", limitsize = FALSE)
